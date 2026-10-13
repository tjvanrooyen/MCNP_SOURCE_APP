# MCNP_SOURCE_APP

**Version 2.1**

A modular Python code that converts a radionuclide **activity inventory** $(i,A_{i})$ into **MCNP6 SDEF**
(source definition) ionising-radiation emission source terms.

## What it does

Given an inventory $(i,A_{i})$ of nuclide identities and their activities in Bq, and the ICRP-107 (2008)
emission data, the code builds the collective ionising radiation emission source term $(E_{g}, Y_{g})_{j}$ for each particle class $j$ and writes, in one run:

| Output | File | Notes |
|--------|------|-------|
| **Spectroscopic SDEF** | `Output_Spectroscopic_Case/SDEF_MCNP6_spectroscopic.out` | full precision; can exceed 20 000 lines |
| **QA report**          | `Output_Spectroscopic_Case/QA_Report.txt` | self-verification of the run; read the verdict first |
| **Binned SDEF**        | `Output_Binned/SDEF_MCNP6_binned.out` | DEP & DEE collapsed onto a fixed energy partition |
| **Overview graph**     | `GraphBank/SDEF_IREST_ALL.{svg,pdf,png}` | binned DEP/DEE + spectroscopic rest; PNG at 600 dpi |

On demand it also writes a **per-nuclide** spectroscopic SDEF for every nuclide
in the database (`Output_Spectroscopic_ALL/<nuclide>.txt`, ~1250 files).

Each nuclide contributes **per average decay of the mixture**: every emission
yield is scaled by the activity fraction `A_i / Σ A`. Discrete lines sharing an
energy are summed, so each class's PAR yield equals the sum of its energy
distribution — the deck is internally consistent.

| Dist. | Class                                | Source                       |
|-------|--------------------------------------|------------------------------|
| 21    | Neutrons (N)                         | `ICRP-07.RAD` `N`            |
| 22    | Discrete-energy photons (DEP)        | `ICRP-07.RAD` `X G PG DG`    |
| 23    | Discrete-energy electrons (DEE)      | `ICRP-07.RAD` `IE` + Auger   |
| 24    | Beta electrons, negative (BEN)       | reconstructed — see below     |
| 25    | Beta positrons (BEP)                 | reconstructed — see below     |
| 26    | Discrete-energy alpha particles (DEA)| `ICRP-07.RAD` `A`            |
| 27    | Discrete-energy recoil atoms (DER)   | `ICRP-07.RAD` `AR`           |
| 28    | Heavy fission fragment (DEHFF)       | partitioned — see below      |
| 29    | Light fission fragment (DELFF)       | partitioned — see below      |

## Energy binning (DEP/DEE)

A spectroscopic source term may carry tens of thousands of photon/electron
lines, which Monte-Carlo codes cannot sample equitably. The binned deck collapses
the DEP and DEE classes onto the partition in
`data_EnergyBins/EnergyBins_PE.dat` (1063 bins, 1e-5 .. 15 MeV): each non-empty
bin becomes one line at the **yield-weighted mean energy** carrying the **summed
yield** (empty bins are dropped). Per-class totals are conserved exactly; all
other classes stay spectroscopically exact.

## Beta-minus / beta-plus separation (corrected in v2.0)

ICRP-107 supplies **one** continuous beta spectrum per nuclide. For the **22
nuclides that emit both beta-minus and beta-plus** that tabulation is the **sum**
of the two branches. This is an approximation; it was corrected in the present code version.

`msa_modules/betaspectra.py` now reconstructs the branches separately, taking the
per-branch yield and mean energy from `ICRP-07.RAD` and the positron spectrum
from the project data file `ICRP-07-CORR.BEP`. **All 952 beta branches in
ICRP-107 reproduce their ICRP-tabulated mean energy to within 2 %.**

The same verification caught a second fault: the `ICRP-07.BET` block for
**Cm-250** is the beta spectrum of its *fission products*, not of Cm-250, and
over-stated that nuclide's beta source ~70-fold.

## Fission-fragment partition (corrected in v2.1)

`ICRP-07.RAD` gives one `FF` line per fissioning nuclide: fragments per decay,
and the mean kinetic energy of **one** fragment. Binary fission gives one light
and one heavy fragment, so each class takes **half** the yield, at its own
energy. Versions before 2.1 gave **both** classes the full yield at the same
energy, with fixed stand-in ZAIDs — for Cm-250, 2.96 fragments/decay against a
physical 1.48, and twice the fragment kinetic energy.

`msa_modules/fission.py` splits the fragments using the near-constant actinide
heavy peak (`A ~ 140`), unchanged charge distribution, and momentum conservation
(`E_L/E_H = A_H/A_L`, `TKE = 2 x Ebar`). The computed fragment energies reproduce
measured values to about 2 % — Cf-252 gives 106.1 / 82.0 MeV against a measured
105.7 / 80.2 — and the ZAIDs are now derived from the dominant fissioning nuclide.

It also detects a yield-convention anomaly: **Cf-252**'s `FF` yield counts
fissions, not fragments, unlike every other nuclide. Read as fragments it implies
a neutron multiplicity of 7.52; read as fissions, 3.760 against the published
3.7676.

## Quality assurance

Every run verifies its own output against ICRP-107 and against the deck as
written to disk, prints a verdict, and writes `QA_Report.txt`. Exit status is 3
on `FAIL`, meaning the deck must not be used.

```
python -m tests.run_tests              # 49-test regression suite
python tools/build_corr_bep.py --check # verify the positron data file
```

## Improvements over the 2009 Necsa RRT code

* **Low-energy completeness.** Electrons and photons are retained down to
  **1.0E-5 MeV (10 eV)**, not the old 1 keV floor. (Validated: emulating the old
  cutoff reproduces the legacy sample's PAR yields exactly.)
* **Detailed Auger electrons.** The `ICRP-07.ACK` Auger / Coster-Kronig spectrum
  *replaces* the coarse RAD `AE` lines where available (per-decay sums are
  identical, so adding them would double-count) — essential for Auger dosimetry.
* **Correct beta-minus / beta-plus separation** (see above) — a defect present
  in the legacy FORTRAN code as well as in the earlier Python versions.
* **Self-verification**: every run writes a QA report and fails loudly.
* **Log-sampled beta spectra**, energy binning, three graph formats.
* Hyphen- and case-insensitive nuclide names (`Na-24`, `Na24`, `na24`, `NA-24`).
* Built-in periodic table (no external `Zlist.txt` needed).

## Usage

```
python MCNP_SOURCE_APP.py                 # spectroscopic + binned SDEF + graph
python MCNP_SOURCE_APP.py --all           # ...and (re)generate per-nuclide ALL output
python MCNP_SOURCE_APP.py --all-only      # ONLY (re)generate per-nuclide ALL output
python MCNP_SOURCE_APP.py --graph-energy-max 100   # widen graph to show fission fragments
python MCNP_SOURCE_APP.py --force         # drop unsupported nuclides and continue
```

The per-nuclide **ALL** output is **off by default** (it writes ~1250 files); run
it on demand with `--all` or `--all-only`. If any inventory nuclide is missing
from `ICRP-07.RAD` the run stops and lists it (manual control is intentional).

Default project layout:

* input    `Input/Input.inp`
* data      `data_Emissions_RadioNuclides/ICRP-07.{RAD,BET,ACK}` + `ICRP-07-CORR.BEP`
* bins      `data_EnergyBins/EnergyBins_PE.dat`
* outputs  `Output_Spectroscopic_Case/`, `Output_Binned/`, `Output_Spectroscopic_ALL/`
* graphs   `GraphBank/`

Requires Python 3.12+, `numpy`, and `matplotlib` (graphs degrade gracefully if
matplotlib is unavailable).

## Package layout (`msa_modules/`)

| Module           | Responsibility                                             |
|------------------|------------------------------------------------------------|
| `elements.py`    | Periodic table; nuclide-name parsing & normalisation       |
| `databases.py`   | Parsers for `ICRP-07.RAD` / `.BET` / `.ACK` / `-CORR.BEP`  |
| `inventory.py`   | Reader for `Input.inp`; availability check                 |
| `betaspectra.py` | BEN / BEP branch separation and verification               |
| `fission.py`     | Light/heavy fission-fragment partition                     |
| `sourceterm.py`  | Activity-weighted source-term construction                 |
| `energybins.py`  | Energy partition + DEP/DEE binning                         |
| `sdef.py`        | MCNP6.3 SDEF writer                                         |
| `qa.py`          | Self-verification of the source term and written decks     |
| `plotting.py`    | Source-term overview graph (SVG/PDF/PNG)                    |

## Documentation

| Document | Contents |
|----------|----------|
| `UserGuide/UserGuide.md` | how to run it, read the outputs, and interpret the QA report |
| `UserGuide/Verification_Record_BEN_BEP.md` | the v2.0 beta correction: cause, fix, and evidence |
| `UserGuide/Verification_Record_Fission_Fragments.md` | the v2.1 fission-fragment correction: cause, fix, and evidence |
