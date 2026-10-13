"""Writer for the MCNP6.3 SDEF (source definition) construct.

The layout follows ``Input_Output/Output_sample.out``:

* a ``mode`` card and commented charged-particle cut-off cards;
* an ``SDEF`` card whose particle type is drawn from distribution ``D7`` and
  whose energy is the per-particle distribution ``D10`` (``ERG = FPAR = D10``);
* the ``SI7 / SP7 / DS10`` block giving the per-class emission probabilities --
  the recoil-atom and fission-fragment entries carry a nuclide-derived heavy-ion
  ZAID, or ``XX`` when that class is absent from the inventory;
* the nine energy distributions (``SI21/SP21`` ... ``SI29/SP29``).

Discrete classes use an ``L`` (line) interpolation law; the continuous beta
spectra use an ``A`` (point-wise pdf) law.  Every distribution referenced by
``DS10`` is always written -- empty classes get a single inert placeholder line
-- so the deck is self-consistent for MCNP.
"""

from __future__ import annotations

from .sourceterm import (
    BETA_E_MIN_MEV,
    INCLUDE_THRESHOLD_MEV,
    BetaDist,
    DiscreteDist,
    SourceTerm,
)

# Commented cut-off energies (MeV).  The electron floor matches MCNP6.3's
# reliable low-energy electron transport limit (~15 eV).
CUT_E_MEV = 1.5e-5
CUT_A_MEV = 1.0e-3
CUT_HI_MEV = 1.0e-3


def _num(x: float) -> str:
    return f"{x:.7E}"


def _data_line(energy: float, value: float) -> str:
    return f"        {_num(energy)}       {_num(value)}"


def _par_line(designator: str, yield_: float, dist: int, desc: str) -> str:
    return f"{designator:>9}         {_num(yield_)}       {dist:<2}     $ {desc}"


def _zaid(value: int | None) -> str:
    """Heavy-ion ZAID for the PAR card, or ``XX`` when the class is empty."""
    return str(value) if value is not None else "XX"


def _placeholder() -> str:
    """An inert line so an empty distribution still parses."""
    return _data_line(INCLUDE_THRESHOLD_MEV, 0.0)


def write_sdef(path: str, st: SourceTerm, *, title: str = "") -> None:
    """Write the complete SDEF deck for source term ``st`` to ``path``."""
    L: list[str] = []

    # Heavy-ion classes carry a nuclide-derived ZAID; "XX" marks a class with no
    # emission in this inventory, so the placeholder is visible rather than
    # looking like a real ion.
    recoil_id = _zaid(st.recoil.zaid)
    hff_id = _zaid(st.hff.zaid)
    lff_id = _zaid(st.lff.zaid)

    # --- Header ----------------------------------------------------------
    L += [
        "c",
        "c Generated SDEF for incorporation into an MCNP6 calculation model",
    ]
    if title:
        L.append(f"c {title}")
    L += [
        "c",
        "mode   N  P  E  A  #",
        "c",
        "c Cutoff Energies for Charged Particle Transport:",
        f"c cut:E      J    {_num(CUT_E_MEV)}  J  J     $ Cutoff energy for electron       transport",
        f"c cut:A      J    {_num(CUT_A_MEV)}  J  J     $ Cutoff energy for alpha-particle transport",
        f"c cut:#      J    {_num(CUT_HI_MEV)}  J  J     $ Cutoff energy for heavy ion      transport",
        "c",
        "SDEF   POS = 0  0  0",
        "       PAR = D7",
        "       ERG = FPAR  = D10",
        "c",
    ]

    # --- PAR / energy-distribution selector block ------------------------
    L += [
        "#      SI7       SP7                 DS10",
        "       L         D                   S",
        _par_line("N", st.neutron.total, 21, "Neutrons"),
        _par_line("P", st.photon.total, 22, "Discrete-Eng Photons (DEP)"),
        _par_line("E", st.dee.total, 23, "Discrete-Eng Electrons (DEE)"),
        _par_line("E", st.ben.total, 24, "Beta-Electrons, negative (BEN)"),
        _par_line("-E", st.bep.total, 25, "Beta-positrons, positive (BEP)"),
        _par_line("A", st.alpha.total, 26, "Discrete-Eng Alpha-particles (DEA)"),
        _par_line(recoil_id, st.recoil.total, 27, "Discrete-Eng Recoil atoms (DER)"),
        _par_line(hff_id, st.hff.total, 28, "Discrete-Eng Heavy Fission Fragment (DEHFF)"),
        _par_line(lff_id, st.lff.total, 29, "Discrete-Eng Light Fission Fragment (DELFF)"),
        "c",
    ]

    # --- Energy distributions 21..29 -------------------------------------
    _discrete_section(L, "Neutron emissions:", 21, st.neutron)
    _discrete_section(L, "Discrete-Energy Photon (DEP) emissions:", 22, st.photon)
    _discrete_section(L, "Discrete-Energy Electron (DEE) emissions:", 23, st.dee)
    _beta_section(L, "Beta-Electron, negative (BEN) emissions:", 24, st.ben)
    _beta_section(L, "Beta-positron (BEP) emissions:", 25, st.bep)
    _discrete_section(L, "Discrete-Energy Alpha-particle (DEA) emissions:", 26, st.alpha)
    _discrete_section(L, "Discrete-energy Recoil Atom (DER) from Alpha-particle emissions:", 27, st.recoil)
    _discrete_section(L, "Discrete Energy Heavy Fission Fragment (DEHFF) emission:", 28, st.hff)
    _discrete_section(L, "Discrete Energy Light Fission Fragment (DELFF) emission:", 29, st.lff)

    L += ["c", "c END OF SDEF", "c", ""]

    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L))


def _discrete_section(L: list[str], heading: str, dist: int, d: DiscreteDist) -> None:
    L += ["c", f"c {heading}", f"#       SI{dist}                SP{dist}", "        L                   D"]
    if d.energy:
        for e, y in zip(d.energy, d.yield_):
            L.append(_data_line(e, y))
    else:
        L.append(_placeholder())
    L.append("c")


def _beta_section(L: list[str], heading: str, dist: int, b: BetaDist) -> None:
    L += ["c", f"c {heading}", f"#       SI{dist}                SP{dist}", "        A                   D"]
    if b.energy:
        for e, y in zip(b.energy, b.pdf):
            L.append(_data_line(e, y))
    else:
        L.append(_data_line(BETA_E_MIN_MEV, 0.0))
    L.append("c")
