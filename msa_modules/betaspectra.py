"""Separation of the beta-minus (BEN) and beta-plus (BEP) emission spectra.

Why this module exists
----------------------
ICRP-107 (2008) supplies, in ``ICRP-07.BET``, exactly **one** continuous beta
spectrum per nuclide.  For the 22 nuclides that emit *both* beta-minus and
beta-plus that single tabulation is the **sum** of the two branches::

    int( BET pdf ) dE  ==  Y_BEN + Y_BEP        (verified for all 22)

Earlier versions of MCNP_SOURCE_APP handed that one tabulation to *both* the
BEN and the BEP distribution.  Two independent errors followed:

1. **Yield.**  Both branches were given the combined yield.  For Cu-64 the deck
   carried BEN = BEP = 5.649E-01 where the truth is BEN = 3.900E-01 and
   BEP = 1.741E-01 -- BEN over-stated by 45 %, BEP by 224 %, and the total
   number of beta particles over-stated by a factor of two.
2. **Spectrum shape.**  The positron branch was given the shape of the combined
   (electron-dominated) spectrum.  Positron and electron spectra differ
   fundamentally: the nuclear Coulomb field *attracts* the electron and
   *repels* the positron, so the positron spectrum is pushed to higher energy.
   For Cu-64 the true mean energies are 0.1904 MeV (BEN) and 0.2782 MeV (BEP)
   -- a 46 % difference -- yet both were written with the combined-spectrum
   mean of 0.2174 MeV.

The authority for the *yields* and the *mean energies* of each branch is
``ICRP-07.RAD``, which tabulates them separately and correctly::

    5 3.90000E-01 1.90428E-01 B-      <- Y_BEN, Ebar_BEN
    4 1.74083E-01 2.78220E-01 B+      <- Y_BEP, Ebar_BEP

This module reconstructs each branch so that **both** of those ICRP-107 numbers
are reproduced.

Method
------
The branch shapes are established in this order of preference.

``BET``
    Single-branch nuclide (only B- *or* only B+).  The ICRP-07.BET tabulation is
    that branch, used directly.  Applies to 908 of the 930 beta emitters.

``CORR``
    Positron branch of a dual emitter.  Taken from ``ICRP-07-CORR.BEP``, a
    project data file of separately evaluated positron spectra.  Every spectrum
    in it reproduces the ICRP-07.RAD ``Y_BEP`` *and* ``Ebar_BEP`` to better than
    0.001 % -- two independent checks against ICRP-107's own numbers -- so the
    shapes are validated against the same authority that supplies the yields.

``BET-CORR``
    Electron branch of a dual emitter: the combined BET tabulation minus the
    positron spectrum above, on the union energy grid, clipped at zero.  This is
    exact in principle, but it is a difference of two nearly equal quantities,
    so it degrades when ``Y_BEN << Y_BEP``.  The result is therefore *tested*
    against ``Ebar_BEN`` and discarded if it fails.

``FERMI``
    Fallback.  An allowed-shape beta spectrum

        N(E) dE  ~  F(Z_d, E) . p . W . (E0 - E)^2 dE

    with the non-relativistic Fermi function ``F`` for the daughter charge
    ``Z_d`` (``Z+1`` for beta-minus, ``Z-1`` for beta-plus) and the end-point
    ``E0`` solved by bisection so that the spectrum's mean energy equals the
    ICRP-07.RAD ``Ebar`` exactly.  Used for 6 branches (see below), all of them
    minor branches whose largest yield is 0.038 per decay.

Every branch, whatever its origin, is finally rescaled so that its integral
equals the ICRP-07.RAD yield exactly.  :func:`decompose` returns the deviation
of the achieved mean energy from ``Ebar`` so the caller can report it -- see
``qa.py``.

Branches that fall back to ``FERMI`` with the shipped data files::

    Br-78  BEN   Y = 1.0E-04     Ag-106 BEN   Y = 1.0E-02
    Cs-130 BEN   Y = 1.6E-02     Rb-84  BEN   Y = 3.8E-02
    Tm-168 BEN   Y = 1.0E-04     Cm-250 BEN   Y = 8.0E-02

The first five are dual emitters whose electron branch is too small to recover
by subtraction.  Cm-250 is different and worth knowing about: its ICRP-07.BET
block is not the Cm-250 beta spectrum at all but the aggregated *fission
product* beta spectrum from its 74 % spontaneous-fission branch (it integrates
to 5.59 betas/decay with a mean of 1.51 MeV, against a tabulated B- branch of
0.08 betas/decay with a mean of 0.0097 MeV).  Using it unchecked -- as earlier
versions did -- over-stated the Cm-250 beta source by a factor of 70.  The mean
energy test catches it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import databases as db
from .elements import Nuclide

# --- Physical constants -----------------------------------------------------
ELECTRON_MASS_MEV = 0.51099895069        # CODATA 2022
FINE_STRUCTURE = 7.2973525643e-3         # CODATA 2022

# Lowest energy on any beta grid (MeV).
BETA_E_MIN_MEV = 1.0e-6

# Points used for an internally generated (FERMI) shape before it is resampled
# onto the user's output grid.  Dense enough that resampling is lossless.
FERMI_SHAPE_POINTS = 4000

# A reconstructed shape is accepted only if its mean energy reproduces the
# ICRP-07.RAD Ebar to within this relative tolerance.  0.02 comfortably clears
# the ~0.2 % trapezoidal-quadrature noise of the ICRP tabulations while
# rejecting every physically wrong reconstruction found in the database.
MEAN_ENERGY_TOL = 0.02

# Origins, most trusted first.
ORIGIN_BET = "BET"            # ICRP-07.BET used directly (single-branch nuclide)
ORIGIN_CORR = "CORR"          # ICRP-07-CORR.BEP evaluated positron spectrum
ORIGIN_SUBTRACT = "BET-CORR"  # combined BET minus the positron spectrum
ORIGIN_FERMI = "FERMI"        # allowed-shape spectrum fitted to Ebar


@dataclass
class BetaBranch:
    """One reconstructed beta branch of one nuclide.

    ``energy``/``pdf`` are the tabulated shape, ascending in energy and scaled so
    that ``trapezoid(pdf, energy) == yield_``.
    """

    particle: str                 # "BEN" or "BEP"
    yield_: float                 # particles per decay (ICRP-07.RAD, authoritative)
    mean_energy: float            # MeV (ICRP-07.RAD, authoritative)
    energy: list[float] = field(default_factory=list)
    pdf: list[float] = field(default_factory=list)
    origin: str = ORIGIN_BET
    mean_dev: float = 0.0         # (achieved <E> - mean_energy) / mean_energy

    @property
    def achieved_mean(self) -> float:
        return self.mean_energy * (1.0 + self.mean_dev)


@dataclass
class BetaBranches:
    """The beta content of one nuclide: either branch may be ``None``."""

    ben: BetaBranch | None = None
    bep: BetaBranch | None = None


# --- ICRP-07.RAD branch summary --------------------------------------------

def rad_branch(riso: db.RadIsotope, ptype: str) -> tuple[float, float]:
    """Return ``(yield, mean energy)`` of the ``B-`` or ``B+`` branch of ``riso``.

    ICRP-07.RAD may list a branch as several transitions; they are combined into
    one yield and one yield-weighted mean energy, which is what a single MCNP
    energy distribution has to represent.  Returns ``(0.0, 0.0)`` if absent.
    """
    lines = [e for e in riso.entries if e.ptype == ptype]
    total = sum(e.yield_ for e in lines)
    if total <= 0.0:
        return 0.0, 0.0
    mean = sum(e.yield_ * e.energy for e in lines) / total
    return total, mean


# --- Allowed-shape beta spectrum -------------------------------------------

def fermi_function(z_daughter: int, energy: np.ndarray, positron: bool) -> np.ndarray:
    """Non-relativistic (Primakoff-Rosen) Fermi function ``F(Z_d, E)``.

    ``eta = +/- alpha Z_d / beta`` with the sign positive for beta-minus (the
    daughter nucleus attracts the electron, enhancing the low-energy end) and
    negative for beta-plus (it repels the positron, suppressing it).
    """
    w = energy + ELECTRON_MASS_MEV
    p = np.sqrt(np.maximum(w * w - ELECTRON_MASS_MEV**2, 1.0e-300))
    beta = p / w
    eta = FINE_STRUCTURE * z_daughter / beta
    if positron:
        eta = -eta
    x = 2.0 * np.pi * eta
    # x -> 0 limit is 1; guard the removable singularity.
    return np.where(np.abs(x) < 1.0e-10, 1.0, x / (1.0 - np.exp(-x)))


def allowed_shape(z_daughter: int, energy: np.ndarray, e0: float,
                  positron: bool) -> np.ndarray:
    """Un-normalised allowed beta shape ``F(Z_d,E) . p . W . (E0 - E)^2``."""
    w = energy + ELECTRON_MASS_MEV
    p = np.sqrt(np.maximum(w * w - ELECTRON_MASS_MEV**2, 0.0))
    return fermi_function(z_daughter, energy, positron) * p * w * np.maximum(e0 - energy, 0.0) ** 2


def _shape_grid(e0: float) -> np.ndarray:
    """Log grid from the beta energy floor to the end-point ``e0``.

    The floor drops below :data:`BETA_E_MIN_MEV` for an end-point that low, so the
    grid can never come out descending.  No ICRP-107 branch is anywhere near
    (the softest is Re-187 at Ebar = 0.62 keV), but a future data revision
    should degrade rather than produce a silently reversed grid.
    """
    lo = min(BETA_E_MIN_MEV, e0 * 1.0e-3)
    return np.logspace(np.log10(lo), np.log10(e0), FERMI_SHAPE_POINTS)


def _shape_mean(z_daughter: int, e0: float, positron: bool) -> float:
    """Mean energy of the allowed shape with end-point ``e0``."""
    grid = _shape_grid(e0)
    s = allowed_shape(z_daughter, grid, e0, positron)
    area = np.trapezoid(s, grid)
    if area <= 0.0:
        return 0.0
    return float(np.trapezoid(s * grid, grid) / area)


def endpoint_for_mean(z_daughter: int, mean_energy: float, positron: bool) -> float:
    """Solve for the end-point ``E0`` whose allowed shape has mean ``mean_energy``.

    The mean rises monotonically with ``E0`` (from roughly ``E0/3`` for a
    low-``Z`` beta-minus emitter to about ``0.4 E0`` for a high-``Z`` positron
    emitter), so plain bisection is safe and needs no derivative.
    """
    lo, hi = mean_energy, max(12.0 * mean_energy, 1.0e-4)
    for _ in range(200):                      # widen until bracketed
        if _shape_mean(z_daughter, hi, positron) >= mean_energy:
            break
        hi *= 2.0
    else:                                     # pragma: no cover - unreachable in practice
        raise ValueError(f"cannot bracket end-point for Ebar = {mean_energy} MeV")

    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _shape_mean(z_daughter, mid, positron) < mean_energy:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def fermi_spectrum(z_daughter: int, yield_: float, mean_energy: float,
                   positron: bool) -> tuple[np.ndarray, np.ndarray]:
    """Allowed-shape spectrum with the ICRP mean energy, scaled to ``yield_``."""
    e0 = endpoint_for_mean(z_daughter, mean_energy, positron)
    grid = _shape_grid(e0)
    pdf = allowed_shape(z_daughter, grid, e0, positron)
    pdf[-1] = 0.0                             # close the spectrum at the end-point
    return grid, pdf


# --- Reconstruction ---------------------------------------------------------

def _normalise(energy: np.ndarray, pdf: np.ndarray, target: float) -> np.ndarray:
    """Scale ``pdf`` so that ``trapezoid(pdf, energy) == target``."""
    area = np.trapezoid(pdf, energy)
    if area <= 0.0:
        return np.zeros_like(pdf)
    return pdf * (target / area)


def _mean(energy: np.ndarray, pdf: np.ndarray) -> float:
    area = np.trapezoid(pdf, energy)
    if area <= 0.0:
        return 0.0
    return float(np.trapezoid(pdf * energy, energy) / area)


def _finish(particle: str, energy: np.ndarray, pdf: np.ndarray, yield_: float,
            mean_energy: float, origin: str) -> BetaBranch:
    """Normalise a candidate shape to the RAD yield and record its mean-energy error."""
    pdf = _normalise(energy, pdf, yield_)
    dev = (_mean(energy, pdf) - mean_energy) / mean_energy if mean_energy > 0.0 else 0.0
    return BetaBranch(
        particle=particle,
        yield_=yield_,
        mean_energy=mean_energy,
        energy=energy.tolist(),
        pdf=pdf.tolist(),
        origin=origin,
        mean_dev=dev,
    )


def _accept_or_fermi(particle: str, energy: np.ndarray, pdf: np.ndarray,
                     yield_: float, mean_energy: float, origin: str,
                     z_daughter: int, positron: bool) -> BetaBranch:
    """Take the candidate shape if it reproduces ``mean_energy``, else fit an
    allowed shape that does."""
    branch = _finish(particle, energy, pdf, yield_, mean_energy, origin)
    if abs(branch.mean_dev) <= MEAN_ENERGY_TOL:
        return branch
    e, p = fermi_spectrum(z_daughter, yield_, mean_energy, positron)
    return _finish(particle, e, p, yield_, mean_energy, ORIGIN_FERMI)


def decompose(
    nuclide: Nuclide | None,
    riso: db.RadIsotope,
    biso: db.BetIsotope | None,
    corr_bep: db.BetIsotope | None,
) -> BetaBranches:
    """Reconstruct the BEN and BEP branches of one nuclide.

    ``riso`` supplies the authoritative per-branch yields and mean energies,
    ``biso`` the ICRP-07.BET continuous spectrum, and ``corr_bep`` the evaluated
    positron spectrum from ICRP-07-CORR.BEP where one exists.  ``nuclide`` gives
    the parent Z needed by the allowed-shape fallback.
    """
    y_ben, ebar_ben = rad_branch(riso, "B-")
    y_bep, ebar_bep = rad_branch(riso, "B+")
    out = BetaBranches()
    if y_ben <= 0.0 and y_bep <= 0.0:
        return out

    z = nuclide.z if nuclide is not None else 0
    bet_e = np.asarray(biso.energy, dtype=float) if biso and biso.energy else None
    bet_p = np.asarray(biso.pdf, dtype=float) if biso and biso.energy else None

    # --- beta-plus --------------------------------------------------------
    if y_bep > 0.0:
        if y_ben > 0.0:
            # Dual emitter: the BET tabulation is the SUM of both branches and
            # must never be used for the positron branch.
            if corr_bep is not None and corr_bep.energy:
                e = np.asarray(corr_bep.energy, dtype=float)
                p = np.asarray(corr_bep.pdf, dtype=float)
                out.bep = _accept_or_fermi("BEP", e, p, y_bep, ebar_bep,
                                           ORIGIN_CORR, z - 1, True)
            else:
                e, p = fermi_spectrum(z - 1, y_bep, ebar_bep, True)
                out.bep = _finish("BEP", e, p, y_bep, ebar_bep, ORIGIN_FERMI)
        elif bet_e is not None:
            # Pure positron emitter: the BET tabulation IS the positron spectrum.
            out.bep = _accept_or_fermi("BEP", bet_e, bet_p, y_bep, ebar_bep,
                                       ORIGIN_BET, z - 1, True)
        else:
            e, p = fermi_spectrum(z - 1, y_bep, ebar_bep, True)
            out.bep = _finish("BEP", e, p, y_bep, ebar_bep, ORIGIN_FERMI)

    # --- beta-minus -------------------------------------------------------
    if y_ben > 0.0:
        if y_bep > 0.0 and bet_e is not None and out.bep is not None:
            # Dual emitter: electron branch = combined BET - positron spectrum.
            bep_e = np.asarray(out.bep.energy, dtype=float)
            bep_p = np.asarray(out.bep.pdf, dtype=float)
            grid = np.unique(np.concatenate([bet_e, bep_e]))
            residual = np.clip(
                np.interp(grid, bet_e, bet_p, left=0.0, right=0.0)
                - np.interp(grid, bep_e, bep_p, left=0.0, right=0.0),
                0.0, None,
            )
            out.ben = _accept_or_fermi("BEN", grid, residual, y_ben, ebar_ben,
                                       ORIGIN_SUBTRACT, z + 1, False)
        elif bet_e is not None:
            out.ben = _accept_or_fermi("BEN", bet_e, bet_p, y_ben, ebar_ben,
                                       ORIGIN_BET, z + 1, False)
        else:
            e, p = fermi_spectrum(z + 1, y_ben, ebar_ben, False)
            out.ben = _finish("BEN", e, p, y_ben, ebar_ben, ORIGIN_FERMI)

    return out
