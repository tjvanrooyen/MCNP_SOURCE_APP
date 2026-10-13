"""Reader for the user activity inventory file ``Input.inp``.

File layout::

    Num_Beta_Eng_Points  500
    @Start
    U-238      6.17400E+07
    U-234      6.17400E+07
    ...
    @End

* The ``Num_Beta_Eng_Points`` line sets how many energy points are used when
  the net beta- and beta+ spectra are sampled.
* The ``@Start`` / ``@End`` block is the ``(nuclide, activity[Bq])`` inventory.

Nuclide names may be written with or without a hyphen and in any case.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .elements import Nuclide, canonical_key, parse_nuclide

DEFAULT_BETA_POINTS = 500


@dataclass
class InventoryItem:
    raw_name: str            # exactly as written by the user
    activity: float          # Bq
    nuclide: Nuclide | None  # parsed identity (None if the name is unparseable)

    @property
    def key(self) -> str:
        return canonical_key(self.raw_name)


@dataclass
class Inventory:
    num_beta_points: int = DEFAULT_BETA_POINTS
    items: list[InventoryItem] = field(default_factory=list)
    unparseable: list[str] = field(default_factory=list)

    @property
    def total_activity(self) -> float:
        return sum(it.activity for it in self.items)

    @property
    def keys(self) -> set[str]:
        return {it.key for it in self.items}


def read_inventory(path: str) -> Inventory:
    """Parse ``Input.inp`` into an :class:`Inventory`.

    Lines outside the ``@Start``/``@End`` block (other than the beta-points
    setting) and blank lines are ignored.  Duplicate nuclides are accumulated
    so that the same nuclide listed twice simply sums its activity.
    """
    inv = Inventory()
    in_block = False
    merged: dict[str, InventoryItem] = {}

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue

            low = line.lower()
            if low.startswith("num_beta_eng_points"):
                parts = line.split()
                if len(parts) >= 2:
                    inv.num_beta_points = int(float(parts[1]))
                continue
            if low.startswith("@start"):
                in_block = True
                continue
            if low.startswith("@end"):
                in_block = False
                continue
            if not in_block:
                continue

            parts = line.split()
            if len(parts) < 2:
                continue
            name = parts[0]
            try:
                activity = float(parts[1])
            except ValueError:
                continue

            nuc = parse_nuclide(name)
            if nuc is None:
                inv.unparseable.append(name)

            key = canonical_key(name)
            if key in merged:
                merged[key].activity += activity
            else:
                merged[key] = InventoryItem(raw_name=name, activity=activity, nuclide=nuc)

    inv.items = list(merged.values())
    return inv


def check_availability(inv: Inventory, rad_keys: set[str]) -> list[InventoryItem]:
    """Return inventory items whose nuclide is absent from the RAD database."""
    return [it for it in inv.items if it.key not in rad_keys]
