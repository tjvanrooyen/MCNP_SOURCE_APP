"""MCNP_SOURCE_APP -- user-facing driver.

Reads a radionuclide activity inventory (``Input/Input.inp``) and the ICRP-107
emission databases, builds the collective ionising-radiation emission source
term, and writes:

* ``Output_Spectroscopic_Case/SDEF_MCNP6_spectroscopic.out`` -- the full,
  spectroscopically exact SDEF deck (can be >20 000 lines);
* ``Output_Binned/SDEF_MCNP6_binned.out``                    -- the same deck but
  with the discrete photon (DEP) and electron (DEE) classes collapsed onto the
  fixed energy partition in ``data_EnergyBins/EnergyBins_PE.dat``;
* ``GraphBank/SDEF_IREST_ALL.svg``                           -- an overview graph
  (binned DEP/DEE plus the spectroscopic remaining classes).

This is the modern Python port of the 2009 FORTRAN ``mcnpsource`` utility.  Key
corrections over the legacy code: electron and photon emissions are retained all
the way down to 1.0E-5 MeV (10 eV), and the detailed ICRP-07.ACK Auger /
Coster-Kronig spectrum is used in place of the coarse RAD ``AE`` lines where
available -- both essential for Auger-emitter dosimetry.

-------------------------------------------------------------------------------
ON-DEMAND "ALL" MODE  (per-nuclide spectroscopic SDEF for the WHOLE database)
-------------------------------------------------------------------------------
The legacy ``ALL`` mode writes one SDEF file per radionuclide for every nuclide
in the ICRP-107 RAD/BET/ACK databases (~1250 files) into
``Output_Spectroscopic_ALL/``.  This is expensive and rarely needs repeating, so
it is **off by default**.  Enable it explicitly:

    python MCNP_SOURCE_APP.py --all              # case run + regenerate ALL
    python MCNP_SOURCE_APP.py --all-only         # regenerate ALL only, skip case

The generated ALL files supersede the older ones in ``Output_Spectroscopic_ALL/
ALL/`` because they now include the low-energy photons and electrons.
-------------------------------------------------------------------------------

Usage
-----
    python MCNP_SOURCE_APP.py [--input PATH] [--emissions DIR] [--bins FILE]
                              [--case-out PATH] [--binned-out PATH]
                              [--all-dir DIR] [--graphdir DIR]
                              [--all] [--all-only] [--no-graph] [--force]
                              [--graph-energy-max MEV]
"""

from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from msa_modules import databases as db                                 # noqa: E402
from msa_modules import energybins, plotting, sdef, sourceterm          # noqa: E402
from msa_modules.elements import canonical_key, parse_nuclide           # noqa: E402
from msa_modules.inventory import (                                     # noqa: E402
    Inventory,
    InventoryItem,
    check_availability,
    read_inventory,
)

# Default project layout (relative to this script).
DEF_INPUT = "Input/Input.inp"
DEF_EMISSIONS = "data_Emissions_RadioNuclides"
DEF_BINS = "data_EnergyBins/EnergyBins_PE.dat"
DEF_CASE_OUT = "Output_Spectroscopic_Case/SDEF_MCNP6_spectroscopic.out"
DEF_BINNED_OUT = "Output_Binned/SDEF_MCNP6_binned.out"
DEF_ALL_DIR = "Output_Spectroscopic_ALL"
DEF_GRAPHDIR = "GraphBank"

# Beta sampling points used for per-nuclide ALL output.
ALL_BETA_POINTS = 500


def _resolve(path: str) -> str:
    return path if os.path.isabs(path) else os.path.join(HERE, path)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="MCNP_SOURCE_APP",
        description="Build MCNP6.3 SDEF source terms from an ICRP-107 radionuclide inventory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", default=DEF_INPUT, help="activity inventory file")
    p.add_argument("--emissions", default=DEF_EMISSIONS, help="folder with the ICRP-07.* data files")
    p.add_argument("--bins", default=DEF_BINS, help="energy-partition file for DEP/DEE binning")
    p.add_argument("--case-out", default=DEF_CASE_OUT, help="spectroscopic SDEF output")
    p.add_argument("--binned-out", default=DEF_BINNED_OUT, help="binned SDEF output")
    p.add_argument("--all-dir", default=DEF_ALL_DIR, help="folder for per-nuclide ALL output")
    p.add_argument("--graphdir", default=DEF_GRAPHDIR, help="folder for output graphs")
    p.add_argument("--all", action="store_true", help="ALSO (re)generate per-nuclide ALL output")
    p.add_argument("--all-only", action="store_true", help="ONLY (re)generate ALL output, skip the case run")
    p.add_argument("--no-graph", action="store_true", help="skip graph generation")
    p.add_argument("--force", action="store_true", help="drop unsupported nuclides and continue")
    p.add_argument("--graph-energy-max", type=float, default=plotting.ENERGY_MAX_MEV,
                   help="upper energy (MeV) of the graph x-axis; raise to view fission fragments")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    emissions_dir = _resolve(args.emissions)
    rad_path = os.path.join(emissions_dir, "ICRP-07.RAD")
    bet_path = os.path.join(emissions_dir, "ICRP-07.BET")
    ack_path = os.path.join(emissions_dir, "ICRP-07.ACK")

    print("=" * 70)
    print(" MCNP_SOURCE_APP  --  ICRP-107 -> MCNP6.3 SDEF source-term builder")
    print("=" * 70)

    for label, path in (("RAD", rad_path), ("BET", bet_path), ("ACK", ack_path)):
        if not os.path.isfile(path):
            print(f"ERROR: {label} file not found: {path}")
            return 2

    # --- ALL-only short-circuit -----------------------------------------
    if args.all_only:
        return _generate_all(rad_path, bet_path, ack_path, _resolve(args.all_dir))

    rc = _run_case(args, rad_path, bet_path, ack_path)
    if rc != 0:
        return rc

    if args.all:
        rc = _generate_all(rad_path, bet_path, ack_path, _resolve(args.all_dir))

    print("\nDone.")
    return rc


def _run_case(args, rad_path, bet_path, ack_path) -> int:
    """Process the user inventory: spectroscopic SDEF + binned SDEF + graph."""
    input_path = _resolve(args.input)
    bins_path = _resolve(args.bins)
    case_out = _resolve(args.case_out)
    binned_out = _resolve(args.binned_out)
    graph_dir = _resolve(args.graphdir)

    if not os.path.isfile(input_path):
        print(f"ERROR: input file not found: {input_path}")
        return 2

    # --- 1. Inventory ----------------------------------------------------
    inv = read_inventory(input_path)
    print(f"\nInventory: {len(inv.items)} nuclide(s), total activity "
          f"{inv.total_activity:.6E} Bq, beta points = {inv.num_beta_points}")
    if inv.unparseable:
        print("  ! Unparseable names (checked literally): " + ", ".join(inv.unparseable))
    if not inv.items:
        print("ERROR: inventory is empty.")
        return 2

    # --- 2. Databases ----------------------------------------------------
    wanted = inv.keys
    print("\nReading ICRP-107 databases ...")
    rad = db.load_rad(rad_path, wanted)
    bet = db.load_bet(bet_path, wanted)
    ack = db.load_ack(ack_path, wanted)
    print(f"  RAD: {len(rad)} matched | BET: {len(bet)} matched | ACK: {len(ack)} matched")

    # --- 3. Availability check ------------------------------------------
    missing = check_availability(inv, set(rad.keys()))
    if missing:
        print("\n" + "-" * 70)
        print("UNSUPPORTED NUCLIDES (no emission data in ICRP-07.RAD):")
        for it in missing:
            print(f"   {it.raw_name:<12} {it.activity:.6E} Bq")
        print("-" * 70)
        if not args.force:
            print("\nProcessing halted. Remove the nuclides above from the input file,")
            print("or re-run with --force to drop them automatically and continue.")
            return 1
        print("\n--force given: dropping the unsupported nuclides and continuing.")
        drop = {it.key for it in missing}
        inv.items = [it for it in inv.items if it.key not in drop]
        if not inv.items:
            print("ERROR: no supported nuclides remain.")
            return 1

    # --- 4. Build the spectroscopic source term -------------------------
    print("\nBuilding activity-weighted emission source term ...")
    st = sourceterm.build_source_term(inv, rad, bet, ack)
    _report(st)

    title = f"Mixture of {len(inv.items)} nuclides, total activity {inv.total_activity:.6E} Bq"

    # --- 5. Spectroscopic SDEF (always) ---------------------------------
    _ensure_dir(case_out)
    sdef.write_sdef(case_out, st, title=title + " -- SPECTROSCOPIC")
    print(f"\nSpectroscopic SDEF written to: {case_out}")

    # --- 6. Binned SDEF (always) ----------------------------------------
    bins = energybins.read_energy_bins(bins_path)
    st_binned = energybins.binned_source_term(st, bins)
    _ensure_dir(binned_out)
    sdef.write_sdef(binned_out, st_binned, title=title + " -- DEP/DEE BINNED")
    print(f"Binned SDEF written to:        {binned_out}")
    print(f"  DEP lines {len(st.photon.energy):>6} -> {len(st_binned.photon.energy):>5} binned   "
          f"|  DEE lines {len(st.dee.energy):>6} -> {len(st_binned.dee.energy):>5} binned")

    # --- 7. Graph (binned DEP/DEE + spectroscopic rest) -----------------
    if not args.no_graph:
        print("Rendering source-term graph ...")
        gtitle = f"All-particle emission spectra $Q(E)$ -- {len(inv.items)}-nuclide mixture"
        outs = plotting.plot_source_term(st_binned, graph_dir, stem="SDEF_IREST",
                                         title=gtitle, energy_max=args.graph_energy_max)
        for out in outs:
            print(f"Graph written to:              {out}")

    return 0


def _generate_all(rad_path, bet_path, ack_path, all_dir) -> int:
    """On-demand: write a per-nuclide spectroscopic SDEF for the whole database."""
    print("\n" + "=" * 70)
    print(" ALL mode: per-nuclide spectroscopic SDEF for the entire ICRP-107 database")
    print("=" * 70)
    print("Loading full RAD / BET / ACK databases (this takes a moment) ...")
    rad = db.load_rad(rad_path, None)
    bet = db.load_bet(bet_path, None)
    ack = db.load_ack(ack_path, None)
    print(f"  RAD: {len(rad)} | BET: {len(bet)} | ACK: {len(ack)} isotopes")

    os.makedirs(all_dir, exist_ok=True)
    total = len(rad)
    written = 0
    for n, (key, riso) in enumerate(sorted(rad.items()), start=1):
        nuc = parse_nuclide(riso.name)
        item = InventoryItem(raw_name=riso.name, activity=1.0, nuclide=nuc)
        inv = Inventory(num_beta_points=ALL_BETA_POINTS, items=[item])
        st = sourceterm.build_source_term(inv, {key: riso}, bet, ack)
        out_path = os.path.join(all_dir, f"{riso.name}.txt")
        sdef.write_sdef(out_path, st, title=f"GENERAL SDEF FOR {riso.name} (per decay)")
        written += 1
        if n % 100 == 0 or n == total:
            print(f"  ... {n}/{total} nuclides")
    print(f"\nALL mode complete: {written} per-nuclide SDEF files written to {all_dir}")
    return 0


def _ensure_dir(path: str) -> None:
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)


def _report(st: sourceterm.SourceTerm) -> None:
    """Print a compact summary of per-class yields (per average decay)."""
    rows = [
        ("Neutrons (N)", st.neutron.total, len(st.neutron.energy)),
        ("Discrete photons (DEP)", st.photon.total, len(st.photon.energy)),
        ("Discrete electrons (DEE)", st.dee.total, len(st.dee.energy)),
        ("Beta- (BEN)", st.ben.total, len(st.ben.energy)),
        ("Beta+ (BEP)", st.bep.total, len(st.bep.energy)),
        ("Alpha (DEA)", st.alpha.total, len(st.alpha.energy)),
        ("Recoil atoms (DER)", st.recoil.total, len(st.recoil.energy)),
        ("Heavy fission frag (HFF)", st.hff.total, len(st.hff.energy)),
        ("Light fission frag (LFF)", st.lff.total, len(st.lff.energy)),
    ]
    print(f"  {'class':<28}{'yield/decay':>16}{'points':>9}")
    for name, total, npts in rows:
        print(f"  {name:<28}{total:>16.7E}{npts:>9}")
    if st.dominant_alpha is not None:
        print(f"  Dominant alpha emitter -> recoil ZAID {st.recoil.zaid} "
              f"(from {st.dominant_alpha.name})")
    if st.auger_from_ack:
        shown = ", ".join(st.auger_from_ack[:12]) + (" ..." if len(st.auger_from_ack) > 12 else "")
        print(f"  Detailed ACK Auger spectra used for: {shown}")
    if st.no_beta_data:
        shown = ", ".join(st.no_beta_data[:12]) + (" ..." if len(st.no_beta_data) > 12 else "")
        print(f"  ! Beta emitters lacking BET spectra (omitted from BEN/BEP shape): {shown}")


if __name__ == "__main__":
    raise SystemExit(main())
