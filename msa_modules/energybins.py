"""Energy binning of the discrete photon (DEP) and electron (DEE) source terms.

A spectroscopically accurate source term can contain tens of thousands of
discrete photon and electron lines -- too many for a Monte-Carlo code to sample
equitably.  This module collapses those two classes onto a fixed energy
partition (``data_EnergyBins/EnergyBins_PE.dat``): every emission falling inside
a bin is merged into a single representative line placed at the **yield-weighted
mean energy** of the bin, carrying the **summed yield**.  Empty bins are dropped.

Only DEP and DEE are binned; all other classes (neutrons, betas, alphas,
recoils, fission fragments) remain spectroscopically exact.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, replace

from .sourceterm import DiscreteDist, SourceTerm


@dataclass(frozen=True)
class EnergyBins:
    """An energy partition: contiguous ``[leb, ueb)`` bins covering a range."""

    leb: list[float]      # lower edges, ascending
    ueb: list[float]      # upper edges, ascending (ueb[i] == leb[i+1])

    @property
    def e_min(self) -> float:
        return self.leb[0]

    @property
    def e_max(self) -> float:
        return self.ueb[-1]

    def __len__(self) -> int:
        return len(self.leb)


def read_energy_bins(path: str) -> EnergyBins:
    """Parse ``EnergyBins_PE.dat`` (``Group  LEB(MeV)  UEB(MeV)``, descending).

    The file lists bins from high to low energy; they are reordered ascending by
    lower edge for fast lookup.  The partition is validated for ``ueb > leb`` and
    contiguity (no gaps or overlaps).
    """
    rows: list[tuple[float, float]] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            tok = line.split()
            if len(tok) < 3 or not tok[0].lstrip("-").isdigit():
                continue                                  # header / blank
            leb, ueb = float(tok[1]), float(tok[2])
            rows.append((leb, ueb))

    if not rows:
        raise ValueError(f"No energy-bin rows parsed from {path!r}")

    rows.sort(key=lambda r: r[0])                          # ascending by lower edge
    leb = [r[0] for r in rows]
    ueb = [r[1] for r in rows]

    for i, (lo, hi) in enumerate(zip(leb, ueb)):
        if not hi > lo:
            raise ValueError(f"{path}: bin {i} has UEB <= LEB ({hi} <= {lo})")
    for i in range(len(leb) - 1):
        if abs(ueb[i] - leb[i + 1]) > 1e-12 * max(1.0, ueb[i]):
            raise ValueError(f"{path}: non-contiguous bins near {ueb[i]} MeV")

    return EnergyBins(leb=leb, ueb=ueb)


def bin_discrete(dist: DiscreteDist, bins: EnergyBins) -> DiscreteDist:
    """Collapse a discrete distribution onto the energy partition.

    Each non-empty bin yields one line at the yield-weighted mean energy with the
    summed yield.  Lines outside the partition are clamped to the nearest edge
    bin (this should not occur for DEP/DEE, which already respect the energy
    floor; gammas above 15 MeV are exceedingly rare).  Output is energy-descending
    to match the spectroscopic sections, and the per-class total is preserved.
    """
    sum_y: dict[int, float] = {}
    sum_ey: dict[int, float] = {}
    n = len(bins)

    for e, y in zip(dist.energy, dist.yield_):
        idx = bisect.bisect_right(bins.leb, e) - 1
        if idx < 0:
            idx = 0
        elif idx >= n:
            idx = n - 1
        sum_y[idx] = sum_y.get(idx, 0.0) + y
        sum_ey[idx] = sum_ey.get(idx, 0.0) + e * y

    out = DiscreteDist(zaid=dist.zaid)
    for idx in sorted(sum_y, reverse=True):               # high energy first
        y = sum_y[idx]
        if y <= 0.0:
            continue
        out.energy.append(sum_ey[idx] / y)
        out.yield_.append(y)
    return out


def binned_source_term(st: SourceTerm, bins: EnergyBins) -> SourceTerm:
    """Return a copy of ``st`` with only the DEP and DEE classes binned."""
    return replace(
        st,
        photon=bin_discrete(st.photon, bins),
        dee=bin_discrete(st.dee, bins),
    )
