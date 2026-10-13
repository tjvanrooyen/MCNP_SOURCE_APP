"""Partition of the spontaneous-fission source into light and heavy fragments.

What ICRP-07.RAD gives
----------------------
One ``FF`` line per spontaneously fissioning nuclide::

    10 1.48000E+00 9.06223E+01 FF        <- Cm-250: yield, mean energy /MeV

* the **energy** is the mean kinetic energy of **one fragment**;
* the **yield** is the number of **fragments** per decay, i.e. twice the
  spontaneous-fission branch.

Both readings are established from the data itself, not assumed:

* Twice the tabulated energy reproduces the measured total fission kinetic
  energy (TKE) for every nuclide with a published value -- Cf-252 188.1 vs
  185.9 MeV, Cm-250 181.2 vs 182.0, U-238 169.7 vs 169.0, Pu-240 176.4 vs
  177.9 -- all within 1.2 %.
* Dividing the yield by two and forming ``Y_n / fissions`` gives a prompt
  neutron multiplicity in the physical range 2.1 - 4.0 for 34 of the 35 nuclides,
  matching the published values (Cf-254 3.89, Fm-256 4.01, Cm-250 3.31,
  Pu-244 2.30).

What was wrong before
---------------------
Both the DEHFF and the DELFF distribution were given the **full** fragment
yield, at the **same** energy.  A deck therefore emitted two fragments for every
fragment ICRP-107 lists -- for Cm-250, 2.96 per decay against a physical 1.48 --
and put both at the average fragment energy instead of at the light and heavy
fragment energies.  The fission-fragment energy deposition was consequently
double the truth.

The partition used here
-----------------------
Binary fission gives exactly one light and one heavy fragment, so each class
takes **half** the fragment yield.

Masses follow the well-established behaviour of actinide spontaneous fission:
the heavy fragment peak sits at an almost constant ``A ~ 140`` while the light
peak moves with the fissioning nucleus.  The fragments share ``A - nu`` between
them, ``nu`` being the prompt neutron multiplicity implied by the RAD neutron
line.  Charges follow the unchanged-charge-distribution rule ``Z_f / A_f = Z / A``
and are forced to sum to the parent ``Z``.

Energies follow from momentum conservation -- the fragments recoil against each
other, so ``E_L / E_H = A_H / A_L`` -- with ``TKE = 2 x`` the tabulated mean::

    E_L = TKE . A_H / (A_L + A_H)        E_H = TKE . A_L / (A_L + A_H)

For Cm-250 this gives 102.9 / 78.4 MeV against measured 103 / 79, and for Cf-252
106.1 / 82.0 MeV against measured 105.7 / 80.2 -- agreement to about 2 %.

The Cf-252 anomaly
------------------
Cf-252 is the one nuclide whose ``FF`` yield is **not** twice its fission
branch: ICRP-107 gives 3.0920E-02, which is the fission branch itself.  Read as
fragments, it implies a neutron multiplicity of 7.52 -- roughly twice any
multiplicity ever measured.  Read as fissions it gives 3.760, against the
published Cf-252 value of 3.7676, one of the best-determined constants in
nuclear physics because Cf-252 is the standard neutron-source calibration.

:func:`fission_rate` therefore tests the implied multiplicity and falls back to
the fissions-per-decay reading when the fragments-per-decay reading is
unphysical.  The correction is recorded in :class:`FissionSplit.note` and
reported by ``qa.py``; it doubles the Cf-252 fission-fragment source, which
would otherwise be half of the truth.
"""

from __future__ import annotations

from dataclasses import dataclass

from . import databases as db
from .elements import Nuclide

# Binary fission: one light plus one heavy fragment.
FRAGMENTS_PER_FISSION = 2

# Prompt neutron multiplicity of actinide spontaneous fission never leaves this
# range; used to decide how the RAD ``FF`` yield is to be read.  The measured
# span across ICRP-107's 35 fissioning nuclides is 2.1 to 4.0.
NU_PHYSICAL_RANGE = (1.0, 5.0)

# Heavy-fragment peak mass for actinide spontaneous fission (nearly independent
# of the fissioning nucleus).
HEAVY_FRAGMENT_MASS = 140.0

# Guard rails on the heavy-fragment mass fraction, so the "heavy" fragment can
# never come out lighter than the "light" one for an unusually light fissioner.
HEAVY_FRACTION_MIN = 0.50
HEAVY_FRACTION_MAX = 0.65

# Fallback heavy mass fraction when the parent nuclide cannot be identified.
FALLBACK_HEAVY_FRACTION = 0.60


@dataclass
class Fragment:
    """One fission-fragment class of one fissioning nuclide."""

    mass: float               # fragment mass number A
    charge: int               # fragment atomic number Z
    energy: float             # mean kinetic energy /MeV
    yield_: float             # fragments per decay of this nuclide

    @property
    def zaid(self) -> int:
        return 1000 * self.charge + int(round(self.mass))


@dataclass
class FissionSplit:
    """The light/heavy fragment partition of one nuclide's fission source."""

    nuclide: str
    fissions: float           # fissions per decay
    nu: float                 # implied prompt neutron multiplicity
    mean_energy: float        # ICRP-07.RAD mean energy per fragment /MeV
    tke: float                # total kinetic energy per fission /MeV
    light: Fragment
    heavy: Fragment
    note: str = ""            # non-empty when the yield reading was corrected

    @property
    def total_energy(self) -> float:
        """Fission-fragment kinetic energy released per decay /MeV."""
        return self.light.energy * self.light.yield_ + self.heavy.energy * self.heavy.yield_


def fission_yield(riso: db.RadIsotope) -> tuple[float, float, float]:
    """Return ``(FF yield, mean fragment energy /MeV, neutron yield)`` from RAD."""
    ff = [e for e in riso.entries if e.ptype == "FF"]
    y_ff = sum(e.yield_ for e in ff)
    e_ff = (sum(e.yield_ * e.energy for e in ff) / y_ff) if y_ff > 0.0 else 0.0
    y_n = sum(e.yield_ for e in riso.entries if e.ptype == "N")
    return y_ff, e_ff, y_n


def fission_rate(y_ff: float, y_n: float) -> tuple[float, float, str]:
    """Fissions per decay implied by the RAD ``FF`` yield.

    Returns ``(fissions, nu, note)``.  Normally the yield counts fragments, so
    there are ``y_ff / 2`` fissions.  That reading is checked against the prompt
    neutron multiplicity it implies; if the result is unphysical the yield is
    taken to count fissions instead (the Cf-252 case -- see the module
    docstring) and ``note`` says so.
    """
    if y_ff <= 0.0:
        return 0.0, 0.0, ""

    fissions = y_ff / FRAGMENTS_PER_FISSION
    if y_n <= 0.0:                       # no neutron line: nothing to test against
        return fissions, 0.0, ""

    lo, hi = NU_PHYSICAL_RANGE
    nu = y_n / fissions
    if lo <= nu <= hi:
        return fissions, nu, ""

    nu_alt = y_n / y_ff
    if lo <= nu_alt <= hi:
        return y_ff, nu_alt, (
            f"RAD FF yield read as fissions/decay, not fragments/decay: the "
            f"fragments reading implies nu = {nu:.2f}, outside the physical "
            f"range {lo:g}-{hi:g}; the fissions reading gives nu = {nu_alt:.2f}"
        )
    return fissions, nu, (
        f"implied neutron multiplicity nu = {nu:.2f} is outside the physical "
        f"range {lo:g}-{hi:g} under either reading of the RAD FF yield"
    )


def split(nuclide: Nuclide | None, riso: db.RadIsotope) -> FissionSplit | None:
    """Partition a nuclide's fission source into light and heavy fragments.

    Returns ``None`` if the nuclide does not spontaneously fission.
    """
    y_ff, e_ff, y_n = fission_yield(riso)
    if y_ff <= 0.0 or e_ff <= 0.0:
        return None

    fissions, nu, note = fission_rate(y_ff, y_n)
    if fissions <= 0.0:
        return None

    tke = FRAGMENTS_PER_FISSION * e_ff

    if nuclide is not None:
        a_f = max(nuclide.a - nu, 1.0)             # mass shared by the fragments
        a_h = min(max(HEAVY_FRAGMENT_MASS, HEAVY_FRACTION_MIN * a_f),
                  HEAVY_FRACTION_MAX * a_f)
        a_l = a_f - a_h
        z_h = int(round(nuclide.z * a_h / a_f))    # unchanged charge distribution
        z_h = min(max(z_h, 1), nuclide.z - 1)
        z_l = nuclide.z - z_h
    else:                                          # unparseable name: shape only
        a_f = 1.0
        a_h = FALLBACK_HEAVY_FRACTION
        a_l = a_f - a_h
        z_h = z_l = 0

    # Momentum conservation: the lighter fragment carries the larger share.
    return FissionSplit(
        nuclide=riso.name,
        fissions=fissions,
        nu=nu,
        mean_energy=e_ff,
        tke=tke,
        light=Fragment(mass=a_l, charge=z_l, energy=tke * a_h / a_f, yield_=fissions),
        heavy=Fragment(mass=a_h, charge=z_h, energy=tke * a_l / a_f, yield_=fissions),
        note=note,
    )
