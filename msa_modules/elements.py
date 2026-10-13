"""Periodic table and radionuclide-name handling.

Two jobs:

1.  Map an element symbol to its atomic number Z (built in, so no external
    ``Zlist.txt`` dependency is required).
2.  Parse and normalise radionuclide names.  The user may write a nuclide with
    or without a hyphen and in any case -- ``Na-24``, ``Na24``, ``na24`` and
    ``na-24`` must all be accepted -- so every name is reduced to a canonical
    key for matching against the ICRP-107 database names.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# --- Atomic numbers Z=1..118 ------------------------------------------------
SYMBOL_TO_Z: dict[str, int] = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
    "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
    "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
    "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
    "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57,
    "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
    "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
    "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
    "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85,
    "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92,
    "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99,
    "Fm": 100, "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105,
    "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111,
    "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117,
    "Og": 118,
}

# Element symbol (1-2 letters), mass number, optional metastable tag (m, n,
# m1, m2, ...).  Matched case-insensitively; the hyphen is optional.
_NUCLIDE_RE = re.compile(
    r"^\s*(?P<sym>[A-Za-z]{1,2})-?(?P<a>\d{1,3})(?P<meta>[mn]\d*)?\s*$"
)


@dataclass(frozen=True)
class Nuclide:
    """A parsed radionuclide identity."""

    symbol: str          # canonical element symbol, e.g. "Pa"
    a: int               # mass number A
    meta: str            # metastable tag in lower case ("", "m", "m2", "n")
    z: int               # atomic number Z

    @property
    def name(self) -> str:
        """Canonical hyphenated name, e.g. ``Pa-234m``."""
        return f"{self.symbol}-{self.a}{self.meta}"

    @property
    def zaid(self) -> int:
        """Plain ZAID ``1000*Z + A`` (metastable state not encoded)."""
        return 1000 * self.z + self.a


def canonical_key(name: str) -> str:
    """Reduce any spelling of a nuclide name to a match key.

    ``Na-24``, ``Na24``, ``na24`` and ``NA-24`` all map to ``"na24"``.
    The key is simply the name lower-cased with hyphens and surrounding
    whitespace removed, so it is robust without needing a full parse.
    """
    return name.strip().lower().replace("-", "")


def parse_nuclide(name: str) -> Nuclide | None:
    """Parse a nuclide name into a :class:`Nuclide`, or ``None`` if invalid.

    The element symbol is normalised to standard capitalisation (first letter
    upper, second lower) and validated against the periodic table.
    """
    m = _NUCLIDE_RE.match(name)
    if not m:
        return None
    raw_sym = m.group("sym")
    sym = raw_sym[0].upper() + raw_sym[1:].lower()
    z = SYMBOL_TO_Z.get(sym)
    if z is None:
        return None
    meta = (m.group("meta") or "").lower()
    return Nuclide(symbol=sym, a=int(m.group("a")), meta=meta, z=z)
