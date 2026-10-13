"""Parsers for the three ICRP-107 (2008) data files.

* ``ICRP-07.RAD`` -- discrete radiations.  Per-isotope block::

      <name> <half-life> <n_entries>
      <type_num> <yield> <energy/MeV> <Z-designator>     (x n_entries)

* ``ICRP-07.BET`` -- continuous beta(+/-) spectra.  Per-isotope block::

      <name> <n_entries>
      <energy/MeV> <pdf>                                 (x n_entries)

  ``pdf`` is the beta probability density; its integral over energy equals the
  number of beta particles emitted per nuclear decay.

* ``ICRP-07.ACK`` -- detailed Auger / Coster-Kronig electron lines.  Block::

      <name> <half-life> <n_entries>
      <yield> <energy/eV> <shell> <shell> <shell>        (x n_entries)

  Note the (yield, energy) column order and the energy unit of eV.  The trailing
  shell labels carry no transport information and are ignored.

* ``ICRP-07-CORR.BEP`` -- corrected beta-plus (positron) spectra for the
  nuclides that emit *both* beta-minus and beta-plus.  Same block layout as the
  BET file, so it is read by :func:`load_bet`.  See ``betaspectra.py`` for why
  it is needed.

Each block is delimited purely by its declared entry count, so blocks are
skipped by counting -- no fragile header/data line discrimination is needed.
Lines whose first token is ``c`` are comments and are ignored, which lets a
data file carry its own provenance and verification header.  Only the isotopes
requested by the caller are retained in memory.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .elements import canonical_key

# Energy unit conversions.
EV_TO_MEV = 1.0e-6


@dataclass
class RadEntry:
    """One discrete radiation line from the RAD file."""

    type_num: int        # particle designator number (1..11, redundant with z)
    yield_: float        # emission probability per decay
    energy: float        # line energy in MeV
    ptype: str           # Z designator: X, G, PG, DG, B-, B+, AE, IE, A, AR, N, FF, ...


@dataclass
class RadIsotope:
    name: str
    half_life: str
    entries: list[RadEntry] = field(default_factory=list)


@dataclass
class BetIsotope:
    name: str
    energy: list[float] = field(default_factory=list)   # MeV, monotonically increasing
    pdf: list[float] = field(default_factory=list)       # beta probability density


@dataclass
class AckIsotope:
    name: str
    half_life: str
    energy: list[float] = field(default_factory=list)   # MeV (converted from eV)
    yield_: list[float] = field(default_factory=list)    # emission probability per decay


def _iter_blocks(path: str):
    """Yield ``(header_tokens, [data_line_tokens, ...])`` for every block.

    Blank lines and ``c``-prefixed comment lines are skipped; a lone ``END``
    token terminates the file.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    i = 0
    n = len(lines)
    while i < n:
        header = lines[i].split()
        if not header or header[0] == "c":
            i += 1
            continue
        if header[0] == "END":
            break
        try:
            count = int(header[-1])
        except ValueError as exc:                       # pragma: no cover - corrupt file
            raise ValueError(
                f"{path}: malformed block header at line {i + 1!r}: {lines[i]!r}"
            ) from exc
        data = [lines[i + 1 + k].split() for k in range(count)]
        yield header, data
        i += 1 + count


def load_rad(path: str, wanted: set[str] | None) -> dict[str, RadIsotope]:
    """Load requested isotopes from the RAD file, keyed by canonical name key.

    ``wanted`` is a set of canonical name keys; pass ``None`` to load every
    isotope in the file (used by the on-demand ALL mode).
    """
    out: dict[str, RadIsotope] = {}
    for header, data in _iter_blocks(path):
        key = canonical_key(header[0])
        if wanted is not None and key not in wanted:
            continue
        iso = RadIsotope(name=header[0], half_life=header[1])
        for tok in data:
            iso.entries.append(
                RadEntry(
                    type_num=int(tok[0]),
                    yield_=float(tok[1]),
                    energy=float(tok[2]),
                    ptype=tok[3],
                )
            )
        out[key] = iso
    return out


def load_bet(path: str, wanted: set[str] | None) -> dict[str, BetIsotope]:
    """Load requested isotopes from the BET file, keyed by canonical name key.

    Pass ``wanted=None`` to load every isotope.
    """
    out: dict[str, BetIsotope] = {}
    for header, data in _iter_blocks(path):
        key = canonical_key(header[0])
        if wanted is not None and key not in wanted:
            continue
        iso = BetIsotope(name=header[0])
        for tok in data:
            iso.energy.append(float(tok[0]))
            iso.pdf.append(float(tok[1]))
        out[key] = iso
    return out


def load_ack(path: str, wanted: set[str] | None) -> dict[str, AckIsotope]:
    """Load requested isotopes from the ACK file, keyed by canonical name key.

    Columns are ``(yield, energy_eV, ...)``; energy is converted to MeV and the
    pair is stored as (energy, yield) ready for transport use.  Pass
    ``wanted=None`` to load every isotope.
    """
    out: dict[str, AckIsotope] = {}
    for header, data in _iter_blocks(path):
        key = canonical_key(header[0])
        if wanted is not None and key not in wanted:
            continue
        iso = AckIsotope(name=header[0], half_life=header[1])
        for tok in data:
            iso.yield_.append(float(tok[0]))
            iso.energy.append(float(tok[1]) * EV_TO_MEV)
        out[key] = iso
    return out
