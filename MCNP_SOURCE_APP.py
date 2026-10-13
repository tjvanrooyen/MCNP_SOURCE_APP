"""MCNP_SOURCE_APP -- user-facing driver.

Reads a radionuclide activity inventory (``Input/Input.inp``) and the ICRP-107
emission databases, builds the collective ionising-radiation emission source
term, verifies it, and writes:

* ``Output_Spectroscopic_Case/SDEF_MCNP6_spectroscopic.out`` -- the full,
  spectroscopically exact SDEF deck (can be >20 000 lines);
* ``Output_Spectroscopic_Case/QA_Report.txt``                -- the quality
  assurance record for that deck (see below);
* ``Output_Binned/SDEF_MCNP6_binned.out``                    -- the same deck but
  with the discrete photon (DEP) and electron (DEE) classes collapsed onto the
  fixed energy partition in ``data_EnergyBins/EnergyBins_PE.dat``;
* ``GraphBank/SDEF_IREST_ALL.{svg,pdf,png}``                 -- an overview graph
  (binned DEP/DEE plus the spectroscopic remaining classes).

This is the modern Python port of the 2009 FORTRAN ``mcnpsource`` utility.

-------------------------------------------------------------------------------
BETA-MINUS / BETA-PLUS SEPARATION  (corrected in v2.0)
-------------------------------------------------------------------------------
ICRP-107 supplies one continuous beta spectrum per nuclide.  For the 22 nuclides
that emit both beta-minus and beta-plus, that tabulation is the SUM of the two
branches.  Versions before 2.0 handed it to both the BEN and the BEP
distribution, so both branches got the combined yield and the positron branch
got the electron spectrum's shape.  For Cu-64 the old deck wrote
BEN = BEP = 5.649E-01 per decay where the truth is BEN = 3.900E-01 and
BEP = 1.741E-01.

Version 2.0 reconstructs each branch separately and reproduces the ICRP-07.RAD
per-branch yield AND mean energy for all 952 beta branches in the database.  See
``msa_modules/betaspectra.py`` for the method and ``msa_modules/qa.py`` for the
verification.  Nuclides in the shipped example inventory that were affected:
Cl-36, Eu-152, Ag-108 and Mn-54.

-------------------------------------------------------------------------------
FISSION-FRAGMENT PARTITION  (corrected in v2.1)
-------------------------------------------------------------------------------
ICRP-07.RAD gives one FF line per spontaneously fissioning nuclide: the number of
fragments per decay, and the mean kinetic energy of ONE fragment.  Binary fission
gives one light and one heavy fragment, so each class takes HALF the yield, at
its own energy.  Versions before 2.1 gave BOTH classes the full yield at the same
energy -- for Cm-250, 2.96 fragments/decay against a physical 1.48, and twice the
fragment kinetic energy.  See ``msa_modules/fission.py``, which also detects the
Cf-252 yield-convention anomaly in ICRP-107.

-------------------------------------------------------------------------------
QUALITY ASSURANCE
-------------------------------------------------------------------------------
Every run verifies its own output against ICRP-107 and against the deck as
written to disk, prints a verdict, and writes ``QA_Report.txt`` beside the
spectroscopic deck.  The process exit status is 0 on PASS or WARN and 3 on FAIL;
a FAIL means the deck must not be used.  Run ``python -m tests.run_tests`` for
the standing regression suite.

-------------------------------------------------------------------------------
ON-DEMAND "ALL" MODE  (per-nuclide spectroscopic SDEF for the WHOLE database)
-------------------------------------------------------------------------------
The legacy ``ALL`` mode writes one SDEF file per radionuclide for every nuclide
in the ICRP-107 RAD/BET/ACK databases (~1250 files) into
``Output_Spectroscopic_ALL/``.  This is expensive and rarely needs repeating, so
it is **off by default**.  Enable it explicitly:

    python MCNP_SOURCE_APP.py --all              # case run + regenerate ALL
    python MCNP_SOURCE_APP.py --all-only         # regenerate ALL only, skip case

The generated ALL files supersede the older ones in
``Output_Spectroscopic_ALL/Archive_2025_values/``: they include the low-energy
photons and electrons, and they carry the corrected BEN/BEP branches and the
corrected light/heavy fission-fragment partition.
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
from msa_modules import energybins, plotting, qa, sdef, sourceterm      # noqa: E402
from msa_modules.elements import parse_nuclide                          # noqa: E402
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
ALL_BETA_POINTS = 750

# Exit status returned when the QA verdict is FAIL.
EXIT_QA_FAIL = 3


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
    corr_path = os.path.join(emissions_dir, "ICRP-07-CORR.BEP")

    print("=" * 70)
    print(" MCNP_SOURCE_APP  --  ICRP-107 -> MCNP6.3 SDEF source-term builder")
    print("=" * 70)

    for label, path in (("RAD", rad_path), ("BET", bet_path), ("ACK", ack_path)):
        if not os.path.isfile(path):
            print(f"ERROR: {label} file not found: {path}")
            return 2
    if not os.path.isfile(corr_path):
        print("\n" + "!" * 70)
        print("WARNING: corrected positron-spectrum file not found:")
        print(f"         {corr_path}")
        print("         Nuclides emitting BOTH beta-minus and beta-plus will fall back")
        print("         to an allowed-shape (Fermi) positron spectrum.  Yields and mean")
        print("         energies stay correct; the spectral shape is approximated.")
        print("!" * 70)
        corr_path = None

    # --- ALL-only short-circuit -----------------------------------------
    if args.all_only:
        return _generate_all(rad_path, bet_path, ack_path, corr_path, _resolve(args.all_dir))

    rc = _run_case(args, rad_path, bet_path, ack_path, corr_path)
    if rc not in (0, EXIT_QA_FAIL):
        return rc

    if args.all:
        rc2 = _generate_all(rad_path, bet_path, ack_path, corr_path, _resolve(args.all_dir))
        rc = rc or rc2

    print("\nDone.")
    return rc


def _run_case(args, rad_path, bet_path, ack_path, corr_path) -> int:
    """Process the user inventory: spectroscopic SDEF + binned SDEF + graph + QA."""
    input_path = _resolve(args.input)
    bins_path = _resolve(args.bins)
    case_out = _resolve(args.case_out)
    binned_out = _resolve(args.binned_out)
    graph_dir = _resolve(args.graphdir)
    qa_out = os.path.join(os.path.dirname(case_out), "QA_Report.txt")

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
    corr = db.load_bet(corr_path, wanted) if corr_path else {}
    print(f"  RAD: {len(rad)} matched | BET: {len(bet)} matched | "
          f"ACK: {len(ack)} matched | CORR.BEP: {len(corr)} matched")

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
    st = sourceterm.build_source_term(inv, rad, bet, ack, corr)
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

    # --- 7. Quality assurance -------------------------------------------
    rep = qa.QAReport()
    qa.check_source_term(st, rep)
    qa.check_partition(bins, rep)
    qa.check_binning(st, st_binned, bins, rep)
    qa.check_deck(case_out, st, rep, label="spectroscopic")
    qa.check_deck(binned_out, st_binned, rep, label="binned")
    qa.write_report(qa_out, rep, st, title=title)
    print(f"QA report written to:          {qa_out}")

    # --- 8. Graph (binned DEP/DEE + spectroscopic rest) -----------------
    if not args.no_graph:
        print("Rendering source-term graph ...")
        gtitle = f"All-particle emission spectra $Q(E)$ -- {len(inv.items)}-nuclide mixture"
        outs = plotting.plot_source_term(st_binned, graph_dir, stem="SDEF_IREST",
                                         title=gtitle, energy_max=args.graph_energy_max)
        for out in outs:
            print(f"Graph written to:              {out}")

    return _announce_qa(rep)


def _announce_qa(rep: qa.QAReport) -> int:
    """Print the QA verdict and translate it into an exit status."""
    counts = rep.counts
    print("\n" + "=" * 70)
    print(f" QA VERDICT: {rep.status}   ({counts[qa.PASS]} passed, "
          f"{counts[qa.WARN]} warnings, {counts[qa.FAIL]} failures)")
    print("=" * 70)
    for c in rep.failures():
        print(str(c))
    for c in rep.warnings():
        print(str(c))
    if rep.status == qa.FAIL:
        print("\n *** THE GENERATED DECK MUST NOT BE USED. See the QA report. ***")
        return EXIT_QA_FAIL
    return 0


def _generate_all(rad_path, bet_path, ack_path, corr_path, all_dir) -> int:
    """On-demand: write a per-nuclide spectroscopic SDEF for the whole database.

    Each nuclide is also QA-checked; any nuclide whose beta reconstruction misses
    the ICRP-07.RAD mean energy is listed at the end and a summary is written to
    ``<all_dir>/QA_Summary_ALL.txt``.
    """
    print("\n" + "=" * 70)
    print(" ALL mode: per-nuclide spectroscopic SDEF for the entire ICRP-107 database")
    print("=" * 70)
    print("Loading full RAD / BET / ACK / CORR.BEP databases (this takes a moment) ...")
    rad = db.load_rad(rad_path, None)
    bet = db.load_bet(bet_path, None)
    ack = db.load_ack(ack_path, None)
    corr = db.load_bet(corr_path, None) if corr_path else {}
    print(f"  RAD: {len(rad)} | BET: {len(bet)} | ACK: {len(ack)} | CORR.BEP: {len(corr)} isotopes")

    os.makedirs(all_dir, exist_ok=True)
    total = len(rad)
    written = 0
    records: list[sourceterm.BetaRecord] = []
    fissions: list[sourceterm.FissionRecord] = []
    offenders: list[str] = []
    for n, (key, riso) in enumerate(sorted(rad.items()), start=1):
        nuc = parse_nuclide(riso.name)
        item = InventoryItem(raw_name=riso.name, activity=1.0, nuclide=nuc)
        inv = Inventory(num_beta_points=ALL_BETA_POINTS, items=[item])
        st = sourceterm.build_source_term(inv, {key: riso}, bet, ack, corr)
        out_path = os.path.join(all_dir, f"{riso.name}.txt")
        sdef.write_sdef(out_path, st, title=f"GENERAL SDEF FOR {riso.name} (per decay)")
        records.extend(st.beta_records)
        fissions.extend(st.fission_records)
        for r in st.beta_records:
            if abs(r.mean_dev) > qa.TOL_MEAN_FAIL:
                offenders.append(f"{r.nuclide} {r.particle} {r.mean_dev:+.2%}")
        written += 1
        if n % 250 == 0 or n == total:
            print(f"  ... {n}/{total} nuclides")

    summary = os.path.join(all_dir, "QA_Summary_ALL.txt")
    _write_all_summary(summary, records, fissions, offenders, written)
    print(f"\nALL mode complete: {written} per-nuclide SDEF files written to {all_dir}")
    print(f"ALL-mode QA summary written to: {summary}")
    if offenders:
        print(f"  ! {len(offenders)} beta branch(es) exceed the {qa.TOL_MEAN_FAIL:.0%} "
              f"mean-energy tolerance -- see the summary.")
        return EXIT_QA_FAIL
    print(f"  All {len(records)} beta branches reproduce the ICRP-07.RAD mean energy "
          f"within {qa.TOL_MEAN_FAIL:.0%}.")
    return 0


def _write_all_summary(path, records, fissions, offenders, written) -> None:
    """Write the database-wide beta and fission census for ALL mode."""
    from collections import Counter
    from datetime import datetime

    census = Counter(r.origin for r in records)
    worst = sorted(records, key=lambda r: -abs(r.mean_dev))[:25]
    L = [
        "=" * 96,
        " MCNP_SOURCE_APP -- ALL-MODE QUALITY ASSURANCE SUMMARY",
        "=" * 96,
        f" Generated : {datetime.now().isoformat(timespec='seconds')}",
        f" Nuclides  : {written} SDEF files",
        f" Beta branches reconstructed: {len(records)}",
        "",
        " Reconstruction method census",
        " " + "-" * 94,
    ]
    for k, v in sorted(census.items()):
        L.append(f"   {k:<12}{v:>6}")
    L += [
        "",
        f" Branches outside the {qa.TOL_MEAN_FAIL:.0%} mean-energy tolerance: {len(offenders)}",
        " " + "-" * 94,
    ]
    L += [f"   {o}" for o in offenders] or ["   (none)"]
    L += [
        "",
        " Largest 25 mean-energy deviations from ICRP-07.RAD",
        " " + "-" * 94,
        f"   {'nuclide':<12}{'branch':<8}{'Y /decay':>14}{'Ebar /MeV':>13}{'<E> dev':>11}  method",
    ]
    for r in worst:
        L.append(f"   {r.nuclide:<12}{r.particle:<8}{r.yield_:>14.6E}"
                 f"{r.mean_energy:>13.6E}{r.mean_dev:>10.3%}  {r.origin}")

    if fissions:
        L += [
            "",
            " Spontaneous-fission fragment partition",
            " " + "-" * 94,
            "   Each of DELFF (29) and DEHFF (28) carries HALF the ICRP-07.RAD FF yield, because",
            "   binary fission gives one light and one heavy fragment.  Energies follow from",
            "   momentum conservation with TKE = 2 x the tabulated mean fragment energy.",
            "",
            f"   {'nuclide':<10}{'fissions':>14}{'nu':>7}{'TKE/MeV':>10}"
            f"{'LFF zaid':>10}{'E_L/MeV':>10}{'HFF zaid':>10}{'E_H/MeV':>10}",
        ]
        for r in sorted(fissions, key=lambda r: -r.split.fissions):
            s = r.split
            L.append(f"   {s.nuclide:<10}{s.fissions:>14.6E}{s.nu:>7.2f}{s.tke:>10.2f}"
                     f"{s.light.zaid:>10d}{s.light.energy:>10.2f}"
                     f"{s.heavy.zaid:>10d}{s.heavy.energy:>10.2f}")
        flagged = [r.split for r in fissions if r.split.note]
        L += ["", f"   Yield-convention corrections applied: {len(flagged)}"]
        for s in flagged:
            L += ["   " + line for line in _wrap_lines(f"{s.nuclide}: {s.note}", 90)]

    L += ["", "=" * 96, ""]
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L))


def _wrap_lines(text: str, width: int) -> list[str]:
    """Wrap ``text`` to ``width``, indenting continuation lines."""
    out, line = [], ""
    for word in text.split():
        if line and len(line) + len(word) + 1 > width:
            out.append(line)
            line = "  " + word
        else:
            line = f"{line} {word}".strip()
    if line:
        out.append(line)
    return out


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
    if st.dominant_fissioner is not None:
        print(f"  Dominant fissioning nuclide -> fragment ZAIDs "
              f"LFF {st.lff.zaid} / HFF {st.hff.zaid} (from {st.dominant_fissioner.name})")
    if st.auger_from_ack:
        shown = ", ".join(st.auger_from_ack[:12]) + (" ..." if len(st.auger_from_ack) > 12 else "")
        print(f"  Detailed ACK Auger spectra used for: {shown}")

    duals = sorted({r.nuclide for r in st.beta_records
                    if any(o.nuclide == r.nuclide and o.particle != r.particle
                           for o in st.beta_records)})
    if duals:
        print(f"  Dual beta emitters (BEN and BEP separated): {', '.join(duals)}")
    if st.no_beta_data:
        shown = ", ".join(st.no_beta_data[:12]) + (" ..." if len(st.no_beta_data) > 12 else "")
        print(f"  ! Beta emitters lacking BET spectra (allowed shape used): {shown}")


if __name__ == "__main__":
    raise SystemExit(main())
