"""Construction of the activity-weighted ionising-radiation emission source term.

Given an inventory ``(i; A_i)`` of nuclides and their activities (Bq) and the
ICRP-107 emission data, this module builds the collective source term
``(E_g; Y_g)_j`` for each particle class ``j``.

Weighting
---------
Every nuclide contributes per *average decay of the mixture*: each emission
yield is scaled by the activity fraction ``af_i = A_i / sum(A)``.  Discrete
lines that share an energy are summed.  The per-class PAR yield is then exactly
the sum of that class's energy-distribution weights, so the SDEF is internally
consistent.

Particle classes (MCNP distribution numbers)
--------------------------------------------
=====  ===========================================  =======================
Dist.  Class                                        Source
=====  ===========================================  =======================
21     Neutrons (N)                                 RAD ``N``
22     Discrete-energy photons (DEP)                RAD ``X G PG DG``
23     Discrete-energy electrons (DEE)              RAD ``IE`` + Auger
24     Beta electrons, negative (BEN)               ``betaspectra`` (B- branch)
25     Beta positrons (BEP)                         ``betaspectra`` (B+ branch)
26     Discrete-energy alpha particles (DEA)        RAD ``A``
27     Discrete-energy recoil atoms (DER)           RAD ``AR``
28     Discrete-energy heavy fission fragment       ``fission`` (RAD ``FF``)
29     Discrete-energy light fission fragment       ``fission`` (RAD ``FF``)
=====  ===========================================  =======================

Beta branches
-------------
The two beta classes are **not** read straight from ICRP-07.BET.  That file
holds one combined spectrum per nuclide, which is the sum of both branches for
a nuclide that emits beta-minus and beta-plus together; using it directly gives
both branches the wrong yield and the positron branch the wrong shape.  See
``betaspectra.py`` for the separation and its verification.  Each branch arrives
here already carrying the authoritative ICRP-07.RAD yield and mean energy.

Fission fragments
-----------------
ICRP-07.RAD carries one ``FF`` line per spontaneously fissioning nuclide, giving
the number of fragments per decay and the mean kinetic energy of *one* fragment.
Binary fission produces one light and one heavy fragment, so each of the two
classes takes **half** the fragment yield, at its own energy.  See
``fission.py`` for the partition and for the Cf-252 yield anomaly it detects.

Auger electrons
---------------
The detailed ICRP-07.ACK Auger / Coster-Kronig spectrum is the fine breakdown
of the *same* electrons the RAD file aggregates as ``AE`` lines (the per-decay
sums are identical).  Where ACK data exist for a nuclide they therefore
*replace* its RAD ``AE`` lines (never add to them) to avoid double-counting; the
RAD ``AE`` lines are used only as a fallback.  This restores the low-energy
electrons the 2009 code discarded -- essential for Auger-emitter dosimetry.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import betaspectra
from . import databases as db
from . import fission
from .betaspectra import BETA_E_MIN_MEV, BetaBranch
from .elements import Nuclide
from .inventory import Inventory

# Lowest emission energy retained, in MeV (= 10 eV).  Electrons and photons are
# kept down to this energy; the historical 1 keV floor is gone.
INCLUDE_THRESHOLD_MEV = 1.0e-5

# RAD Z-designators grouped by emission class.
PHOTON_TYPES = frozenset({"X", "G", "PG", "DG"})   # AQ deliberately excluded
DEE_LINE_TYPES = frozenset({"IE"})                  # IC electrons (Auger handled separately)


@dataclass
class DiscreteDist:
    """A discrete-line distribution: parallel energy/yield lists, energy-descending."""

    energy: list[float] = field(default_factory=list)
    yield_: list[float] = field(default_factory=list)
    zaid: int | None = None        # heavy-ion ZAID for recoil/FF classes

    @property
    def total(self) -> float:
        return float(sum(self.yield_))

    def __bool__(self) -> bool:
        return bool(self.yield_)


@dataclass
class BetaDist:
    """A continuous beta distribution sampled on a log energy grid (ascending).

    ``expected_total`` and ``expected_mean`` are the activity-weighted ICRP-107
    values the tabulation is required to reproduce.  The tabulated pdf is scaled
    so that :attr:`total` equals ``expected_total`` exactly; :attr:`mean` is then
    a free result and is checked against ``expected_mean`` by ``qa.py``.
    """

    energy: list[float] = field(default_factory=list)
    pdf: list[float] = field(default_factory=list)
    expected_total: float = 0.0
    expected_mean: float = 0.0

    @property
    def total(self) -> float:
        """Trapezoidal integral of the pdf = beta particles per average decay."""
        if len(self.energy) < 2:
            return 0.0
        return float(np.trapezoid(self.pdf, self.energy))

    @property
    def mean(self) -> float:
        """Mean energy of the tabulated spectrum (MeV)."""
        t = self.total
        if t <= 0.0:
            return 0.0
        return float(np.trapezoid(np.asarray(self.pdf) * np.asarray(self.energy),
                                  self.energy) / t)

    def __bool__(self) -> bool:
        return any(p > 0.0 for p in self.pdf)


@dataclass
class BetaRecord:
    """How one nuclide's beta branch was reconstructed -- for the QA report."""

    nuclide: str
    particle: str                 # "BEN" or "BEP"
    activity_fraction: float
    yield_: float                 # per decay of THIS nuclide (ICRP-07.RAD)
    mean_energy: float            # MeV (ICRP-07.RAD)
    origin: str
    mean_dev: float               # relative error of the reconstructed mean energy


@dataclass
class FissionRecord:
    """One nuclide's fission-fragment partition -- for the QA report."""

    split: fission.FissionSplit
    activity_fraction: float


@dataclass
class SourceTerm:
    neutron: DiscreteDist
    photon: DiscreteDist
    dee: DiscreteDist
    ben: BetaDist
    bep: BetaDist
    alpha: DiscreteDist
    recoil: DiscreteDist
    hff: DiscreteDist
    lff: DiscreteDist
    total_activity: float
    dominant_alpha: Nuclide | None = None
    dominant_fissioner: Nuclide | None = None
    auger_from_ack: list[str] = field(default_factory=list)   # nuclides using ACK detail
    no_beta_data: list[str] = field(default_factory=list)     # B emitters lacking BET data
    beta_records: list[BetaRecord] = field(default_factory=list)
    fission_records: list[FissionRecord] = field(default_factory=list)


def _merge_discrete(accum: dict[float, float]) -> DiscreteDist:
    """Turn an ``{energy: yield}`` map into an energy-descending DiscreteDist."""
    dist = DiscreteDist()
    for e in sorted(accum, reverse=True):
        y = accum[e]
        if y > 0.0:
            dist.energy.append(e)
            dist.yield_.append(y)
    return dist


def build_source_term(
    inv: Inventory,
    rad: dict[str, db.RadIsotope],
    bet: dict[str, db.BetIsotope],
    ack: dict[str, db.AckIsotope],
    corr_bep: dict[str, db.BetIsotope] | None = None,
) -> SourceTerm:
    """Assemble the collective source term from the inventory and ICRP-107 data.

    ``corr_bep`` holds the evaluated positron spectra of ICRP-07-CORR.BEP; pass
    ``None`` only in tests -- without it every dual beta emitter falls back to an
    allowed-shape positron spectrum.
    """
    corr_bep = corr_bep or {}
    total_activity = inv.total_activity
    if total_activity <= 0.0:
        raise ValueError("Total inventory activity is zero; nothing to build.")

    # Discrete-line accumulators keyed by exact energy.
    neutron: dict[float, float] = {}
    photon: dict[float, float] = {}
    dee: dict[float, float] = {}
    alpha: dict[float, float] = {}
    recoil: dict[float, float] = {}
    lff: dict[float, float] = {}
    hff: dict[float, float] = {}

    auger_from_ack: list[str] = []
    no_beta_data: list[str] = []
    beta_records: list[BetaRecord] = []
    fission_records: list[FissionRecord] = []

    # Reconstructed beta branches awaiting activity-weighted mixing.
    ben_contrib: list[tuple[BetaBranch, float]] = []
    bep_contrib: list[tuple[BetaBranch, float]] = []
    best_alpha_yield = 0.0
    dominant_alpha: Nuclide | None = None
    best_fission_yield = 0.0
    dominant_fissioner: Nuclide | None = None
    dominant_split: fission.FissionSplit | None = None

    def add(accum: dict[float, float], energy: float, value: float) -> None:
        if energy >= INCLUDE_THRESHOLD_MEV and value > 0.0:
            accum[energy] = accum.get(energy, 0.0) + value

    for item in inv.items:
        riso = rad.get(item.key)
        if riso is None:
            continue
        af = item.activity / total_activity

        rad_alpha = 0.0
        has_ack = item.key in ack

        for e in riso.entries:
            w = e.yield_ * af
            t = e.ptype
            if t == "N":
                add(neutron, e.energy, w)
            elif t in PHOTON_TYPES:
                add(photon, e.energy, w)
            elif t == "IE":
                add(dee, e.energy, w)
            elif t == "AE":
                if not has_ack:                      # fallback only
                    add(dee, e.energy, w)
            elif t == "A":
                add(alpha, e.energy, w)
                rad_alpha += e.yield_
            elif t == "AR":
                add(recoil, e.energy, w)
            # B- / B+ are handled by betaspectra.decompose below, FF by
            # fission.split -- neither can be read line-by-line from RAD.

        # Detailed Auger / Coster-Kronig electrons replace the RAD AE lines.
        if has_ack:
            aiso = ack[item.key]
            for energy, y in zip(aiso.energy, aiso.yield_):
                add(dee, energy, y * af)
            auger_from_ack.append(riso.name)

        # Dominant alpha contributor (activity-weighted) fixes the recoil ZAID.
        alpha_contrib = rad_alpha * af
        if alpha_contrib > best_alpha_yield and item.nuclide is not None:
            best_alpha_yield = alpha_contrib
            dominant_alpha = item.nuclide

        # --- Beta branches, correctly separated --------------------------
        branches = betaspectra.decompose(
            item.nuclide, riso, bet.get(item.key), corr_bep.get(item.key)
        )
        if (branches.ben is not None or branches.bep is not None) and item.key not in bet:
            no_beta_data.append(riso.name)
        for branch, sink in ((branches.ben, ben_contrib), (branches.bep, bep_contrib)):
            if branch is None:
                continue
            sink.append((branch, af))
            beta_records.append(BetaRecord(
                nuclide=riso.name,
                particle=branch.particle,
                activity_fraction=af,
                yield_=branch.yield_,
                mean_energy=branch.mean_energy,
                origin=branch.origin,
                mean_dev=branch.mean_dev,
            ))

        # --- Fission fragments, partitioned light / heavy -----------------
        fsplit = fission.split(item.nuclide, riso)
        if fsplit is not None:
            add(lff, fsplit.light.energy, fsplit.light.yield_ * af)
            add(hff, fsplit.heavy.energy, fsplit.heavy.yield_ * af)
            fission_records.append(FissionRecord(split=fsplit, activity_fraction=af))
            # Dominant fission contributor (activity-weighted) fixes the
            # fragment ZAIDs, exactly as the dominant alpha fixes the recoil.
            contrib = fsplit.fissions * af
            if contrib > best_fission_yield and item.nuclide is not None:
                best_fission_yield = contrib
                dominant_fissioner = item.nuclide
                dominant_split = fsplit

    recoil_dist = _merge_discrete(recoil)
    if dominant_alpha is not None:
        recoil_dist.zaid = 1000 * (dominant_alpha.z - 2) + (dominant_alpha.a - 4)

    lff_dist = _merge_discrete(lff)
    hff_dist = _merge_discrete(hff)
    if dominant_split is not None:
        lff_dist.zaid = dominant_split.light.zaid
        hff_dist.zaid = dominant_split.heavy.zaid

    return SourceTerm(
        neutron=_merge_discrete(neutron),
        photon=_merge_discrete(photon),
        dee=_merge_discrete(dee),
        ben=_mix_beta(ben_contrib, inv.num_beta_points),
        bep=_mix_beta(bep_contrib, inv.num_beta_points),
        alpha=_merge_discrete(alpha),
        recoil=recoil_dist,
        hff=hff_dist,
        lff=lff_dist,
        total_activity=total_activity,
        dominant_alpha=dominant_alpha,
        dominant_fissioner=dominant_fissioner,
        auger_from_ack=auger_from_ack,
        no_beta_data=sorted(set(no_beta_data)),
        beta_records=beta_records,
        fission_records=fission_records,
    )


def _mix_beta(contrib: list[tuple[BetaBranch, float]], num_points: int) -> BetaDist:
    """Combine per-nuclide beta branches onto one log energy grid.

    Each branch is linearly interpolated onto the shared grid, zero outside its
    own energy range, and weighted by the nuclide's activity fraction.  The grid
    spans ``BETA_E_MIN_MEV`` to the largest beta energy among the contributors.

    Resampling onto a fixed grid is not exactly area-preserving, so the mixed
    spectrum is finally rescaled to the activity-weighted ICRP-107 yield.  That
    makes the deck exactly self-consistent -- the PAR yield on the ``SP7`` card
    equals the integral of the ``SI24``/``SP24`` tabulation beneath it -- and
    costs nothing physically, since the absolute scale of an MCNP ``SP`` card is
    normalised internally by MCNP anyway.  ``qa.py`` reports the size of the
    correction and the resulting mean-energy error.
    """
    dist = BetaDist()
    contrib = [(b, af) for b, af in contrib if b.energy and b.yield_ > 0.0]
    if not contrib:
        return dist

    dist.expected_total = sum(af * b.yield_ for b, af in contrib)
    if dist.expected_total > 0.0:
        dist.expected_mean = (
            sum(af * b.yield_ * b.mean_energy for b, af in contrib) / dist.expected_total
        )

    e_max = max(b.energy[-1] for b, _ in contrib)
    if e_max <= BETA_E_MIN_MEV:
        return dist

    n = max(int(num_points), 2)
    grid = np.logspace(np.log10(BETA_E_MIN_MEV), np.log10(e_max), n)
    spectrum = np.zeros_like(grid)
    for branch, af in contrib:
        xp = np.asarray(branch.energy)
        fp = np.asarray(branch.pdf)
        spectrum += af * np.interp(grid, xp, fp, left=0.0, right=0.0)

    area = float(np.trapezoid(spectrum, grid))
    if area > 0.0:
        spectrum *= dist.expected_total / area

    dist.energy = grid.tolist()
    dist.pdf = spectrum.tolist()
    return dist
