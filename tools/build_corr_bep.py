"""Rebuild ``data_Emissions_RadioNuclides/ICRP-07-CORR.BEP`` from the reference library.

    python tools/build_corr_bep.py [--library DIR] [--out FILE] [--check]

``ICRP-07-CORR.BEP`` holds the separately evaluated beta-plus (positron)
spectrum of every nuclide that emits *both* beta-minus and beta-plus.  ICRP-107
does not supply those separately -- its single ICRP-07.BET tabulation is the sum
of both branches -- so they are taken from the institutional reference emission
library at ``C:\\Py_Data_Shared\\Source_Emissions``, one ``<nuclide>.xlsx`` per
nuclide, sheet ``BEP``.

The generated file is committed to the project so that ordinary runs never touch
the shared library.  This script exists to make that step reproducible and
auditable: run it to refresh the file after a library update, or with
``--check`` to confirm the committed file still matches the library.

Verification
------------
Every extracted spectrum is checked against ICRP-07.RAD, which tabulates the
beta-plus branch yield and mean energy independently of any spectrum:

    integral( pdf dE )    ==  Y_BEP        from the RAD ``B+`` line
    integral( E pdf dE ) / integral( pdf dE )  ==  Ebar_BEP

Both must agree to better than 0.01 %.  The script refuses to write a file that
fails, and records the achieved agreement in the output file's header so the
evidence travels with the data.

Requires ``openpyxl`` (not needed for ordinary MCNP_SOURCE_APP runs).
"""

from __future__ import annotations

import argparse
import datetime
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from msa_modules import betaspectra                                     # noqa: E402
from msa_modules import databases as db                                 # noqa: E402

DEF_LIBRARY = r"C:\Py_Data_Shared\Source_Emissions"
DEF_OUT = os.path.join(ROOT, "data_Emissions_RadioNuclides", "ICRP-07-CORR.BEP")
RAD_PATH = os.path.join(ROOT, "data_Emissions_RadioNuclides", "ICRP-07.RAD")

# The extracted spectrum must reproduce both ICRP-07.RAD numbers this closely.
TOLERANCE = 1.0e-4


def dual_emitters(rad: dict[str, db.RadIsotope]) -> list[tuple[str, str, float, float]]:
    """``(name, key, Y_BEP, Ebar_BEP)`` for nuclides emitting B- and B+ together."""
    out = []
    for key, iso in rad.items():
        y_ben, _ = betaspectra.rad_branch(iso, "B-")
        y_bep, ebar = betaspectra.rad_branch(iso, "B+")
        if y_ben > 0.0 and y_bep > 0.0:
            out.append((iso.name, key, y_bep, ebar))
    return sorted(out)


def read_bep_sheet(library: str, name: str) -> tuple[np.ndarray, np.ndarray]:
    """Read the ``BEP`` sheet of ``<library>/<nuclide>.xlsx`` as (energy, pdf)."""
    import openpyxl

    path = os.path.join(library, name.replace("-", "") + ".xlsx")
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    try:
        rows = [r for r in wb["BEP"].iter_rows(min_row=3, values_only=True)
                if r[0] is not None]
    finally:
        wb.close()
    e = np.array([float(r[0]) for r in rows])
    p = np.array([float(r[1]) for r in rows])
    order = np.argsort(e)
    e, p = e[order], p[order]

    nz = np.nonzero(p)[0]
    if nz.size == 0:
        raise ValueError(f"{path}: BEP sheet is all zero")
    # Keep one zero anchor on each side of the support; drop the rest.
    return e[max(nz[0] - 1, 0):nz[-1] + 2], p[max(nz[0] - 1, 0):nz[-1] + 2]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--library", default=DEF_LIBRARY, help="reference emission library folder")
    ap.add_argument("--out", default=DEF_OUT, help="data file to write")
    ap.add_argument("--check", action="store_true",
                    help="verify the existing file instead of rewriting it")
    args = ap.parse_args(argv)

    rad = db.load_rad(RAD_PATH, None)
    duals = dual_emitters(rad)
    print(f"ICRP-07.RAD: {len(rad)} nuclides, {len(duals)} emit both B- and B+")

    if args.check:
        return _check(args.out, duals)

    if not os.path.isdir(args.library):
        print(f"ERROR: reference library not found: {args.library}")
        return 2

    blocks, failures = [], []
    for name, key, y_bep, ebar in duals:
        try:
            e, p = read_bep_sheet(args.library, name)
        except (FileNotFoundError, ValueError) as exc:
            failures.append(f"{name}: {exc}")
            continue
        integral = float(np.trapezoid(p, e))
        mean = float(np.trapezoid(p * e, e) / integral)
        dy, dm = (integral - y_bep) / y_bep, (mean - ebar) / ebar
        if abs(dy) > TOLERANCE or abs(dm) > TOLERANCE:
            failures.append(f"{name}: yield dev {dy:+.3%}, mean-energy dev {dm:+.3%}")
            continue
        blocks.append((name, e, p, y_bep, ebar, integral, mean, dy, dm))

    if failures:
        print("\nREFUSING TO WRITE -- these nuclides failed verification:")
        for f in failures:
            print(f"   {f}")
        return 1

    _write(args.out, blocks, args.library)
    print(f"\nWrote {args.out}")
    print(f"  {len(blocks)} positron spectra, all reproducing the ICRP-07.RAD "
          f"yield and mean energy to better than {TOLERANCE:.2%}")
    return 0


def _check(path: str, duals) -> int:
    """Verify the committed file against ICRP-07.RAD without touching the library."""
    if not os.path.isfile(path):
        print(f"ERROR: {path} does not exist")
        return 1
    corr = db.load_bet(path, None)
    problems = []
    for name, key, y_bep, ebar in duals:
        iso = corr.get(key)
        if iso is None:
            problems.append(f"{name}: absent from {os.path.basename(path)}")
            continue
        e, p = np.asarray(iso.energy), np.asarray(iso.pdf)
        integral = float(np.trapezoid(p, e))
        mean = float(np.trapezoid(p * e, e) / integral)
        dy, dm = (integral - y_bep) / y_bep, (mean - ebar) / ebar
        flag = "ok " if abs(dy) <= TOLERANCE and abs(dm) <= TOLERANCE else "FAIL"
        if flag == "FAIL":
            problems.append(f"{name}: yield dev {dy:+.3%}, mean-energy dev {dm:+.3%}")
        print(f"  [{flag}] {name:<10} Y {integral:.6E} vs {y_bep:.6E} ({dy:+.2E})   "
              f"Ebar {mean:.6E} vs {ebar:.6E} ({dm:+.2E})")
    extra = sorted(set(corr) - {k for _, k, _, _ in duals})
    if extra:
        problems.append(f"unexpected nuclides in the file: {extra}")
    print()
    if problems:
        for p in problems:
            print(f"  PROBLEM: {p}")
        return 1
    print(f"  All {len(duals)} positron spectra verify against ICRP-07.RAD.")
    return 0


def _write(path: str, blocks, library: str) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"""c ICRP-07-CORR.BEP -- corrected beta-plus (positron) emission spectra
c
c PURPOSE
c   ICRP-107 (2008) supplies ONE continuous beta spectrum per nuclide in
c   ICRP-07.BET.  For a nuclide that emits BOTH beta-minus and beta-plus that
c   single tabulation is the SUM of the two branches.  Using it for the
c   positron branch (as MCNP_SOURCE_APP did before v2.0) assigns the
c   beta-minus shape and the combined yield to BEP, which is wrong on both
c   counts.  This file supplies the separate, physically correct positron
c   spectrum for every such nuclide.
c
c FORMAT   identical to ICRP-07.BET:
c   <nuclide>  <n_points>
c   <energy/MeV>  <pdf / (positrons per decay per MeV)>      x n_points
c   ...
c   END
c
c NORMALISATION
c   The trapezoidal integral of each spectrum equals the ICRP-07.RAD B+ branch
c   yield Y_BEP, and its first moment equals the ICRP-07.RAD B+ mean energy
c   Ebar_BEP.  Both are reproduced to better than 0.01 % -- see the per-nuclide
c   verification table below.  These are two independent constraints taken from
c   ICRP-107 itself, so the shapes are validated against the same authority
c   that supplies the yields.
c
c PROVENANCE
c   Generated {datetime.date.today().isoformat()} by tools/build_corr_bep.py from the institutional
c   reference emission library {library}\\<nuclide>.xlsx, sheet "BEP".
c   Re-verify at any time with:  python tools/build_corr_bep.py --check
c
c VERIFICATION (integral vs Y_BEP, mean energy vs Ebar_BEP, both from ICRP-07.RAD)
c
c   nuclide      Y_BEP        int(pdf)     dev%      Ebar_BEP     <E>          dev%
c   ------------------------------------------------------------------------------
""")
        for name, e, p, y, ebar, integral, mean, dy, dm in blocks:
            fh.write(f"c   {name:<12} {y:.5E}  {integral:.5E}  {100 * dy:+7.4f}   "
                     f"{ebar:.5E}  {mean:.5E}  {100 * dm:+7.4f}\n")
        fh.write("c\n")
        for name, e, p, *_ in blocks:
            fh.write(f"{name:<12}{len(e):>6}\n")
            for ee, pp in zip(e, p):
                fh.write(f"{ee:.7E} {pp:.7E}\n")
        fh.write("END\n")


if __name__ == "__main__":
    raise SystemExit(main())
