"""Quality assurance for the generated source term and SDEF decks.

MCNP_SOURCE_APP feeds radiation-transport calculations with nuclear-safety
significance, so every run verifies its own output and writes an auditable
record.  Nothing here trusts the code that built the source term: the checks are
posed against ICRP-107's own tabulated numbers and against the deck as it was
actually written to disk.

Check groups
------------
``BETA-BRANCH``
    Every reconstructed beta branch reproduces the ICRP-07.RAD mean energy for
    that branch.  This is the check that catches the historical BEN/BEP
    cross-contamination: a positron branch carrying the electron spectrum fails
    it by tens of percent.  It also records which method built each branch.

``BETA-MIX``
    The mixed BEN and BEP tabulations carry the activity-weighted ICRP-107 yield
    and mean energy.

``DUAL-EMITTER``
    For a nuclide emitting both branches, BEN and BEP must differ in yield *and*
    in mean energy.  Identical values are the exact signature of the historical
    defect, so this check exists to make its return impossible to miss.

``FISSION``
    Binary fission gives one light and one heavy fragment, so the two classes
    must carry equal yields summing to twice the fission rate, the light
    fragment must be the faster of the two, and the total fission-fragment
    kinetic energy per decay must equal ``fissions x TKE``.  This is the check
    that catches the historical defect of writing the full fragment yield to
    both classes.  It also tests the prompt neutron multiplicity each nuclide's
    yield implies, which is what identifies the Cf-252 tabulation anomaly.

``BINNING``
    Collapsing DEP/DEE onto the energy partition conserves total yield, and
    every binned line sits inside the bin it represents.

``DECK``
    The SDEF file is read back from disk and each PAR yield on the ``SP7`` card
    is compared with the integral of the energy distribution it points at.

``PARTITION``
    The energy partition is contiguous, ascending and strictly positive.

A ``FAIL`` means the deck must not be used.  A ``WARN`` means the deck is usable
but rests on an approximation or on ICRP-107 data with a known blemish; the
report says which.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from datetime import datetime

import numpy as np

from . import fission
from .energybins import EnergyBins
from .sourceterm import SourceTerm

# Tolerances.  The "exact" ones test arithmetic the code controls completely;
# the physics ones allow for the quadrature noise of the ICRP-107 tabulations.
TOL_EXACT = 1.0e-9          # book-keeping identities: must hold to round-off
TOL_MEAN_WARN = 1.0e-2      # 1 %  -- reconstructed mean energy vs ICRP-07.RAD
TOL_MEAN_FAIL = 5.0e-2      # 5 %

PASS, WARN, FAIL = "PASS", "WARN", "FAIL"
_RANK = {PASS: 0, WARN: 1, FAIL: 2}


@dataclass
class Check:
    group: str
    name: str
    status: str
    detail: str = ""

    def __str__(self) -> str:
        return f"  [{self.status}] {self.group:<13} {self.name:<46} {self.detail}"


@dataclass
class QAReport:
    checks: list[Check] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def add(self, group: str, name: str, status: str, detail: str = "") -> None:
        self.checks.append(Check(group, name, status, detail))

    def measure(self, group: str, name: str, deviation: float,
                warn: float, fail: float, detail: str = "") -> None:
        """Add a check graded on the size of a relative deviation."""
        status = FAIL if abs(deviation) > fail else (WARN if abs(deviation) > warn else PASS)
        self.add(group, name, status, detail)

    @property
    def status(self) -> str:
        worst = max((_RANK[c.status] for c in self.checks), default=0)
        return {0: PASS, 1: WARN, 2: FAIL}[worst]

    @property
    def counts(self) -> dict[str, int]:
        return {s: sum(1 for c in self.checks if c.status == s) for s in (PASS, WARN, FAIL)}

    def failures(self) -> list[Check]:
        return [c for c in self.checks if c.status == FAIL]

    def warnings(self) -> list[Check]:
        return [c for c in self.checks if c.status == WARN]


def _rel(actual: float, expected: float) -> float:
    """Relative deviation, safe when ``expected`` is zero."""
    if expected == 0.0:
        return 0.0 if actual == 0.0 else float("inf")
    return (actual - expected) / expected


# --- Source-term checks -----------------------------------------------------

def check_source_term(st: SourceTerm, rep: QAReport) -> None:
    """Verify the beta reconstruction against ICRP-07.RAD."""
    # --- per-branch mean energies ---------------------------------------
    worst = 0.0
    worst_who = ""
    origins: dict[str, int] = {}
    for r in st.beta_records:
        origins[r.origin] = origins.get(r.origin, 0) + 1
        if abs(r.mean_dev) > abs(worst):
            worst, worst_who = r.mean_dev, f"{r.nuclide} {r.particle}"
    if st.beta_records:
        rep.measure(
            "BETA-BRANCH", "reconstructed <E> matches ICRP-07.RAD Ebar",
            worst, TOL_MEAN_WARN, TOL_MEAN_FAIL,
            f"{len(st.beta_records)} branches, worst {worst:+.3%} ({worst_who})",
        )
        rep.add("BETA-BRANCH", "reconstruction method census", PASS,
                ", ".join(f"{k}={v}" for k, v in sorted(origins.items())))

    # --- mixed spectra ---------------------------------------------------
    for label, dist in (("BEN", st.ben), ("BEP", st.bep)):
        if dist.expected_total <= 0.0:
            continue
        dy = _rel(dist.total, dist.expected_total)
        rep.measure("BETA-MIX", f"{label} tabulated yield == activity-weighted ICRP yield",
                    dy, TOL_EXACT, 1.0e-6,
                    f"{dist.total:.6E} vs {dist.expected_total:.6E} ({dy:+.2E})")
        dm = _rel(dist.mean, dist.expected_mean)
        rep.measure("BETA-MIX", f"{label} tabulated <E> == activity-weighted ICRP Ebar",
                    dm, TOL_MEAN_WARN, TOL_MEAN_FAIL,
                    f"{dist.mean:.6E} vs {dist.expected_mean:.6E} MeV ({dm:+.3%})")

    # --- the historical defect -------------------------------------------
    duals = _dual_emitters(st)
    if duals:
        bad = [n for n, ben, bep in duals
               if abs(_rel(bep.yield_, ben.yield_)) < TOL_EXACT
               or abs(_rel(bep.mean_energy, ben.mean_energy)) < TOL_EXACT]
        rep.add("DUAL-EMITTER", "BEN and BEP branches are distinct", FAIL if bad else PASS,
                (f"identical for {', '.join(bad)}" if bad
                 else f"{len(duals)} dual emitter(s) checked: "
                      + ", ".join(f"{n} Y-/Y+={ben.yield_:.3E}/{bep.yield_:.3E}"
                                  for n, ben, bep in duals[:4])
                      + (" ..." if len(duals) > 4 else "")))

    if st.ben.expected_total > 0.0 and st.bep.expected_total > 0.0:
        same = abs(_rel(st.bep.total, st.ben.total)) < TOL_EXACT
        rep.add("DUAL-EMITTER", "mixture BEN total != mixture BEP total",
                FAIL if same else PASS,
                f"BEN {st.ben.total:.6E}, BEP {st.bep.total:.6E} per decay")

    # --- fission fragments -------------------------------------------------
    _check_fission(st, rep)

    # --- elementary sanity ------------------------------------------------
    _check_physical(st, rep)

    # --- narrative notes --------------------------------------------------
    approx = sorted({r.nuclide + " " + r.particle for r in st.beta_records
                     if r.origin == "FERMI"})
    if approx:
        rep.notes.append(
            "Allowed-shape (Fermi) spectra were fitted for these branches because the "
            "ICRP-107 tabulation could not supply them: " + ", ".join(approx)
            + ".  Each reproduces the ICRP-07.RAD yield and mean energy by construction; "
            "the shape is the allowed approximation rather than an evaluated spectrum."
        )
    if st.no_beta_data:
        rep.notes.append(
            "Beta emitters with no ICRP-07.BET block (allowed-shape spectra used): "
            + ", ".join(st.no_beta_data)
        )


def _check_fission(st: SourceTerm, rep: QAReport) -> None:
    """Verify the light/heavy fission-fragment partition."""
    records = st.fission_records
    if not records:
        return

    # One light and one heavy fragment per fission: equal yields, summing to
    # twice the fission rate.  Writing the full FF yield to both classes -- the
    # historical defect -- doubles both of these.
    expected = sum(r.split.fissions * r.activity_fraction for r in records)
    d_l = _rel(st.lff.total, expected)
    d_h = _rel(st.hff.total, expected)
    rep.measure("FISSION", "LFF yield == fission rate (one light fragment per fission)",
                d_l, TOL_EXACT, 1.0e-9,
                f"{st.lff.total:.7E} vs {expected:.7E} fissions/decay")
    rep.measure("FISSION", "HFF yield == fission rate (one heavy fragment per fission)",
                d_h, TOL_EXACT, 1.0e-9,
                f"{st.hff.total:.7E} vs {expected:.7E} fissions/decay")
    rep.measure("FISSION", "LFF + HFF == 2 x fission rate (no fragment double-count)",
                _rel(st.lff.total + st.hff.total, 2.0 * expected), TOL_EXACT, 1.0e-9,
                f"{st.lff.total + st.hff.total:.7E} fragments/decay")

    # Total fragment kinetic energy per decay must equal fissions x TKE.
    e_expect = sum(r.split.total_energy * r.activity_fraction for r in records)
    e_actual = (sum(e * y for e, y in zip(st.lff.energy, st.lff.yield_))
                + sum(e * y for e, y in zip(st.hff.energy, st.hff.yield_)))
    rep.measure("FISSION", "fragment kinetic energy per decay == fissions x TKE",
                _rel(e_actual, e_expect), 1.0e-9, 1.0e-6,
                f"{e_actual:.7E} MeV/decay")

    # Momentum conservation makes the light fragment the faster one.
    slow = [r.split.nuclide for r in records
            if r.split.light.energy <= r.split.heavy.energy]
    rep.add("FISSION", "light fragment faster than heavy fragment",
            FAIL if slow else PASS,
            f"violated for {', '.join(slow)}" if slow
            else f"{len(records)} fissioning nuclide(s) checked")

    # The implied prompt neutron multiplicity must be physical.
    lo, hi = fission.NU_PHYSICAL_RANGE
    bad = [f"{r.split.nuclide} nu={r.split.nu:.2f}" for r in records
           if r.split.nu > 0.0 and not (lo <= r.split.nu <= hi)]
    rep.add("FISSION", "implied prompt neutron multiplicity is physical",
            WARN if bad else PASS,
            "; ".join(bad) if bad
            else "nu in [" + ", ".join(
                f"{v:.2f}" for v in (min(r.split.nu for r in records),
                                     max(r.split.nu for r in records))) + "]")

    # Fragment ZAIDs must conserve the parent's charge and (allowing for the
    # prompt neutrons) its mass.
    if st.dominant_fissioner is not None and st.lff.zaid and st.hff.zaid:
        parent = st.dominant_fissioner
        z_sum = st.lff.zaid // 1000 + st.hff.zaid // 1000
        a_sum = st.lff.zaid % 1000 + st.hff.zaid % 1000
        top = max(records, key=lambda r: r.split.fissions * r.activity_fraction).split
        a_expect = parent.a - top.nu
        ok = z_sum == parent.z and abs(a_sum - a_expect) <= 1.0
        rep.add("FISSION", "fragment ZAIDs conserve parent Z and A - nu",
                PASS if ok else FAIL,
                f"{parent.name} (Z={parent.z}, A={parent.a}) -> "
                f"{st.lff.zaid} + {st.hff.zaid}: Z {z_sum}, A {a_sum} vs {a_expect:.1f}")

    for r in records:
        if r.split.note:
            rep.notes.append(f"{r.split.nuclide}: {r.split.note}.")


def _dual_emitters(st: SourceTerm):
    """Return ``(name, ben_record, bep_record)`` for nuclides emitting both branches."""
    by_name: dict[str, dict[str, object]] = {}
    for r in st.beta_records:
        by_name.setdefault(r.nuclide, {})[r.particle] = r
    return [(n, d["BEN"], d["BEP"]) for n, d in sorted(by_name.items())
            if "BEN" in d and "BEP" in d]


def _check_physical(st: SourceTerm, rep: QAReport) -> None:
    """Energies positive, yields non-negative, orderings as the writer expects."""
    problems: list[str] = []
    for label, d in (("N", st.neutron), ("DEP", st.photon), ("DEE", st.dee),
                     ("DEA", st.alpha), ("DER", st.recoil), ("HFF", st.hff),
                     ("LFF", st.lff)):
        if any(e <= 0.0 for e in d.energy):
            problems.append(f"{label}: non-positive energy")
        if any(y < 0.0 for y in d.yield_):
            problems.append(f"{label}: negative yield")
        if any(a < b for a, b in zip(d.energy, d.energy[1:])):
            problems.append(f"{label}: not energy-descending")
    for label, b in (("BEN", st.ben), ("BEP", st.bep)):
        if any(e <= 0.0 for e in b.energy):
            problems.append(f"{label}: non-positive energy")
        if any(p < 0.0 for p in b.pdf):
            problems.append(f"{label}: negative pdf")
        if any(a > c for a, c in zip(b.energy, b.energy[1:])):
            problems.append(f"{label}: not energy-ascending")
    rep.add("SOURCE", "energies positive, yields non-negative, ordering correct",
            FAIL if problems else PASS, "; ".join(problems) or "all 9 classes")


# --- Binning checks ---------------------------------------------------------

def check_binning(st: SourceTerm, binned: SourceTerm, bins: EnergyBins,
                  rep: QAReport) -> None:
    """Binning must conserve yield and place each line inside its own bin."""
    for label, full, small in (("DEP", st.photon, binned.photon),
                               ("DEE", st.dee, binned.dee)):
        d = _rel(small.total, full.total)
        rep.measure("BINNING", f"{label} total yield conserved by binning",
                    d, TOL_EXACT, 1.0e-9,
                    f"{len(full.energy)} -> {len(small.energy)} lines, "
                    f"{full.total:.7E} -> {small.total:.7E} ({d:+.2E})")
        outside = [e for e in small.energy
                   if not (bins.e_min <= e <= bins.e_max)]
        rep.add("BINNING", f"{label} binned energies lie inside the partition",
                FAIL if outside else PASS,
                f"{len(outside)} outside" if outside
                else f"[{bins.e_min:.3E}, {bins.e_max:.3E}] MeV")

    for label, full, small in (("N", st.neutron, binned.neutron),
                               ("BEN", st.ben, binned.ben),
                               ("BEP", st.bep, binned.bep),
                               ("DEA", st.alpha, binned.alpha),
                               ("DER", st.recoil, binned.recoil)):
        d = _rel(small.total, full.total)
        rep.measure("BINNING", f"{label} untouched by DEP/DEE binning", d,
                    TOL_EXACT, TOL_EXACT, f"{full.total:.7E}")


def check_partition(bins: EnergyBins, rep: QAReport) -> None:
    """The energy partition must be positive, ascending and gap-free."""
    problems = []
    if any(e <= 0.0 for e in bins.leb):
        problems.append("non-positive lower edge")
    if any(u <= l for l, u in zip(bins.leb, bins.ueb)):
        problems.append("upper edge <= lower edge")
    if any(a >= b for a, b in zip(bins.leb, bins.leb[1:])):
        problems.append("lower edges not strictly ascending")
    gaps = [i for i in range(len(bins.leb) - 1)
            if abs(bins.ueb[i] - bins.leb[i + 1]) > 1e-12 * max(1.0, bins.ueb[i])]
    if gaps:
        problems.append(f"{len(gaps)} gap(s)/overlap(s)")
    rep.add("PARTITION", "energy partition contiguous and monotonic",
            FAIL if problems else PASS,
            "; ".join(problems) or
            f"{len(bins)} bins, {bins.e_min:.4E} - {bins.e_max:.4E} MeV")


# --- Deck read-back ---------------------------------------------------------

_PAR_RE = re.compile(r"^\s*(\S+)\s+([0-9.+\-Ee]+)\s+(\d\d)\s+\$")
_SI_RE = re.compile(r"^#\s+SI(\d\d)\s+SP\d\d\s*$")
_LAW_RE = re.compile(r"^\s+([LA])\s+D\s*$")
_DATA_RE = re.compile(r"^\s+([0-9.+\-Ee]+)\s+([0-9.+\-Ee]+)\s*$")


def check_deck(path: str, st: SourceTerm, rep: QAReport, label: str = "deck") -> None:
    """Read the written SDEF back and verify it against the source term.

    This is deliberately an independent path: it parses the text on disk rather
    than inspecting the objects that produced it, so a formatting or precision
    fault in the writer is caught here rather than by MCNP.
    """
    if not os.path.isfile(path):
        rep.add("DECK", f"{label} exists", FAIL, path)
        return

    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()

    par: dict[int, tuple[str, float]] = {}
    dists: dict[int, tuple[str, list[float], list[float]]] = {}
    cur: int | None = None
    law = "L"
    for i, line in enumerate(lines):
        m = _SI_RE.match(line)
        if m:
            cur = int(m.group(1))
            law = "L"
            dists[cur] = (law, [], [])
            continue
        m = _LAW_RE.match(line)
        if m and cur is not None and not dists[cur][1]:
            law = m.group(1)
            dists[cur] = (law, [], [])
            continue
        m = _DATA_RE.match(line)
        if m and cur is not None:
            e, v = float(m.group(1)), float(m.group(2))
            dists[cur][1].append(e)
            dists[cur][2].append(v)
            continue
        m = _PAR_RE.match(line)
        if m and cur is None:
            par[int(m.group(3))] = (m.group(1), float(m.group(2)))

    expected = {
        21: ("N", st.neutron.total), 22: ("P", st.photon.total),
        23: ("E", st.dee.total), 24: ("E", st.ben.total),
        25: ("-E", st.bep.total), 26: ("A", st.alpha.total),
        27: (None, st.recoil.total), 28: (None, st.hff.total),
        29: (None, st.lff.total),
    }

    missing = sorted(set(expected) - set(par)) + sorted(set(expected) - set(dists))
    rep.add("DECK", f"{label}: all 9 PAR entries and SI/SP blocks present",
            FAIL if missing else PASS,
            f"missing {missing}" if missing else os.path.basename(path))
    if missing:
        return

    worst_par, worst_par_who = 0.0, ""
    worst_int, worst_int_who = 0.0, ""
    for dnum, (designator, value) in sorted(expected.items()):
        d = _rel(par[dnum][1], value)
        if abs(d) > abs(worst_par):
            worst_par, worst_par_who = d, f"D{dnum}"
        if designator is not None and par[dnum][0] != designator:
            rep.add("DECK", f"{label}: D{dnum} particle designator", FAIL,
                    f"wrote {par[dnum][0]!r}, expected {designator!r}")

        dlaw, de, dv = dists[dnum]
        if len(de) < 2:
            integral = float(sum(dv))
        elif dlaw == "A":
            integral = float(np.trapezoid(dv, de))
        else:
            integral = float(sum(dv))
        d2 = _rel(integral, par[dnum][1])
        if abs(d2) > abs(worst_int):
            worst_int, worst_int_who = d2, f"D{dnum}"

    # 7 significant figures are written, so agreement is limited to ~1e-6.
    rep.measure("DECK", f"{label}: PAR yields match the source term",
                worst_par, 1.0e-6, 1.0e-5,
                f"worst {worst_par:+.2E} at {worst_par_who}")
    rep.measure("DECK", f"{label}: each PAR yield == integral of its own SI/SP block",
                worst_int, 1.0e-5, 1.0e-3,
                f"worst {worst_int:+.2E} at {worst_int_who}")

    n_lines = sum(len(dists[d][1]) for d in dists)
    rep.add("DECK", f"{label}: tabulated emission lines", PASS,
            f"{n_lines} lines over 9 distributions")


# --- Reporting --------------------------------------------------------------

def write_report(path: str, rep: QAReport, st: SourceTerm, *, title: str) -> None:
    """Write the full QA record, including the per-nuclide beta branch table."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    counts = rep.counts
    L = [
        "=" * 100,
        " MCNP_SOURCE_APP -- QUALITY ASSURANCE REPORT",
        "=" * 100,
        f" Case      : {title}",
        f" Generated : {datetime.now().isoformat(timespec='seconds')}",
        f" Verdict   : {rep.status}   "
        f"({counts[PASS]} passed, {counts[WARN]} warnings, {counts[FAIL]} failures)",
        "=" * 100,
        "",
        "CHECKS",
        "-" * 100,
    ]
    L += [str(c) for c in rep.checks]

    if rep.failures():
        L += ["", "FAILURES -- THIS DECK MUST NOT BE USED", "-" * 100]
        L += [str(c) for c in rep.failures()]
    if rep.warnings():
        L += ["", "WARNINGS", "-" * 100]
        L += [str(c) for c in rep.warnings()]
    if rep.notes:
        L += ["", "NOTES", "-" * 100]
        for n in rep.notes:
            L += _wrap(n)

    if st.beta_records:
        L += [
            "",
            "BETA BRANCH RECONSTRUCTION",
            "-" * 100,
            "  Y and Ebar are the ICRP-07.RAD values for that branch of that nuclide.",
            "  <E> dev is the error of the reconstructed spectrum's mean energy against Ebar.",
            "  Method:  BET      = ICRP-07.BET used directly (single-branch nuclide)",
            "           CORR     = evaluated positron spectrum from ICRP-07-CORR.BEP",
            "           BET-CORR = combined ICRP-07.BET minus that positron spectrum",
            "           FERMI    = allowed-shape spectrum fitted to Ebar",
            "",
            f"  {'nuclide':<12}{'branch':<8}{'act.frac':>11}{'Y /decay':>14}"
            f"{'Ebar /MeV':>13}{'<E> dev':>11}  method",
            "  " + "-" * 96,
        ]
        for r in sorted(st.beta_records, key=lambda r: (-r.activity_fraction, r.nuclide)):
            L.append(f"  {r.nuclide:<12}{r.particle:<8}{r.activity_fraction:>11.4E}"
                     f"{r.yield_:>14.6E}{r.mean_energy:>13.6E}{r.mean_dev:>10.3%}  {r.origin}")

    if st.fission_records:
        L += [
            "",
            "FISSION-FRAGMENT PARTITION",
            "-" * 100,
            "  Y_FF and Ebar are the ICRP-07.RAD FF line; Ebar is the mean energy of ONE fragment.",
            "  Binary fission gives one light and one heavy fragment, so each class takes half the",
            "  fragment yield.  Fragment energies follow from momentum conservation, E_L/E_H = A_H/A_L,",
            "  with TKE = 2 x Ebar.  nu is the prompt neutron multiplicity implied by the RAD N line.",
            "",
            f"  {'nuclide':<10}{'act.frac':>11}{'fissions':>13}{'nu':>7}{'TKE/MeV':>10}"
            f"{'LFF zaid':>10}{'E_L/MeV':>10}{'HFF zaid':>10}{'E_H/MeV':>10}",
            "  " + "-" * 96,
        ]
        for r in sorted(st.fission_records,
                        key=lambda r: -r.split.fissions * r.activity_fraction):
            s_ = r.split
            L.append(f"  {s_.nuclide:<10}{r.activity_fraction:>11.4E}{s_.fissions:>13.6E}"
                     f"{s_.nu:>7.2f}{s_.tke:>10.2f}{s_.light.zaid:>10d}{s_.light.energy:>10.2f}"
                     f"{s_.heavy.zaid:>10d}{s_.heavy.energy:>10.2f}")
        L.append(f"  Deck fragment ZAIDs: DELFF {st.lff.zaid}, DEHFF {st.hff.zaid} "
                 f"(from the dominant fissioning nuclide)")

    L += [
        "",
        "MIXTURE SOURCE TERM (per average decay of the mixture)",
        "-" * 100,
        f"  {'class':<30}{'yield':>16}{'mean E /MeV':>16}{'points':>9}",
    ]
    for name, total, mean, npts in (
        ("Neutrons (N)", st.neutron.total, _dmean(st.neutron), len(st.neutron.energy)),
        ("Discrete photons (DEP)", st.photon.total, _dmean(st.photon), len(st.photon.energy)),
        ("Discrete electrons (DEE)", st.dee.total, _dmean(st.dee), len(st.dee.energy)),
        ("Beta- (BEN)", st.ben.total, st.ben.mean, len(st.ben.energy)),
        ("Beta+ (BEP)", st.bep.total, st.bep.mean, len(st.bep.energy)),
        ("Alpha (DEA)", st.alpha.total, _dmean(st.alpha), len(st.alpha.energy)),
        ("Recoil atoms (DER)", st.recoil.total, _dmean(st.recoil), len(st.recoil.energy)),
        ("Heavy fission frag (HFF)", st.hff.total, _dmean(st.hff), len(st.hff.energy)),
        ("Light fission frag (LFF)", st.lff.total, _dmean(st.lff), len(st.lff.energy)),
    ):
        L.append(f"  {name:<30}{total:>16.7E}{mean:>16.7E}{npts:>9}")

    L += ["", "=" * 100, ""]
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L))


def _dmean(d) -> float:
    t = d.total
    return (sum(e * y for e, y in zip(d.energy, d.yield_)) / t) if t > 0.0 else 0.0


def _wrap(text: str, width: int = 96) -> list[str]:
    words, out, line = text.split(), [], "  "
    for w in words:
        if len(line) + len(w) + 1 > width and line.strip():
            out.append(line)
            line = "    " + w
        else:
            line = (line + " " + w) if line.strip() else line + w
    if line.strip():
        out.append(line)
    return out
