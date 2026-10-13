# MCNP_SOURCE_APP — User Guide

**Version 2.1**

MCNP_SOURCE_APP turns a radionuclide activity inventory into an **MCNP6.3 SDEF**
(source definition) ionising-radiation emission source term, using the ICRP-107
(2008) nuclear decay data.

---

## 1. What the code does

You give it an inventory $(i; A_i)$ — a list of radionuclides and their
activities in Bq. It reads the characteristic emissions of each nuclide from the
ICRP-107 databases and builds the collective emission source term
$(E_g; Y_g)_j$ for each particle class $j$, expressed **per average decay of the
mixture**: every emission yield is scaled by that nuclide's activity fraction
$A_i / \sum A$.

Each run writes:

| Output | File |
|---|---|
| Spectroscopic SDEF | `Output_Spectroscopic_Case/SDEF_MCNP6_spectroscopic.out` |
| Quality-assurance report | `Output_Spectroscopic_Case/QA_Report.txt` |
| Binned SDEF | `Output_Binned/SDEF_MCNP6_binned.out` |
| Overview graph | `GraphBank/SDEF_IREST_ALL.{svg,pdf,png}` |

On demand it also writes one SDEF per radionuclide for the whole ICRP-107
database into `Output_Spectroscopic_ALL/` (~1250 files).

---

## 2. Running it

```
python MCNP_SOURCE_APP.py
```

or double-click `Run.cmd`.

| Command | Effect |
|---|---|
| `python MCNP_SOURCE_APP.py` | spectroscopic + binned SDEF, graph, QA report |
| `python MCNP_SOURCE_APP.py --all` | the above, **and** regenerate the per-nuclide `ALL` library |
| `python MCNP_SOURCE_APP.py --all-only` | regenerate the `ALL` library only |
| `python MCNP_SOURCE_APP.py --no-graph` | skip the graph |
| `python MCNP_SOURCE_APP.py --force` | drop unsupported nuclides and continue |
| `python MCNP_SOURCE_APP.py --graph-energy-max 200` | widen the graph energy axis |

The per-nuclide `ALL` mode is **off by default** because it writes ~1250 files.
It rarely needs repeating — do it after a data-file update or a change to the
emission physics.

**Exit status**

| Code | Meaning |
|---|---|
| 0 | success; QA verdict PASS or WARN |
| 1 | inventory contains nuclides absent from ICRP-07.RAD (re-run with `--force` to drop them) |
| 2 | a required input or data file is missing |
| 3 | **QA verdict FAIL — do not use the deck** |

Requires Python 3.10+, `numpy`, and `matplotlib` (the graph is skipped with a
warning if matplotlib is absent; the SDEF decks are unaffected).

---

## 3. The input file — `Input/Input.inp`

```
Num_Beta_Eng_Points  900
@Start
Co-60    8.63200E+10
Fe-55    5.00600E+10
Ni-63    3.08600E+10
...
@End
```

* `Num_Beta_Eng_Points` sets how many energy points represent the BEN and BEP
  continuous spectra. 500–1000 is ample.
* The inventory sits between `@Start` and `@End`.
* Activities are in **Bq**.
* Nuclide names are flexible: `Co60`, `Co-60`, `CO-60` and `co60` are all
  accepted. Metastable states keep their tag: `Ag-110m`, `Nb-93m`.
* Columns may be separated by spaces, tabs, or both.
* A nuclide listed twice has its activities summed.

If any nuclide is missing from `ICRP-07.RAD` the run **stops** and lists it. That
is deliberate: silently dropping a nuclide from a safety calculation is worse
than stopping. Use `--force` only when you have decided the omission is
acceptable.

---

## 4. The output decks

### 4.1 Particle classes

| Dist. | Class | Source |
|---|---|---|
| 21 | Neutrons (N) | `ICRP-07.RAD` `N` |
| 22 | Discrete-energy photons (DEP) | `ICRP-07.RAD` `X G PG DG` |
| 23 | Discrete-energy electrons (DEE) | `ICRP-07.RAD` `IE` + `ICRP-07.ACK` Auger |
| 24 | Beta electrons, negative (BEN) | reconstructed — see §5 |
| 25 | Beta positrons (BEP) | reconstructed — see §5 |
| 26 | Discrete-energy alpha particles (DEA) | `ICRP-07.RAD` `A` |
| 27 | Discrete-energy recoil atoms (DER) | `ICRP-07.RAD` `AR` |
| 28 | Heavy fission fragment (DEHFF) | partitioned — see §5.1 |
| 29 | Light fission fragment (DELFF) | partitioned — see §5.1 |

Annihilation quanta (`AQ`) are deliberately **not** included as photons. A
positron is emitted instead and MCNP produces the annihilation photons itself
during transport. Including both would double-count them.

### 4.2 Spectroscopic vs binned

The **spectroscopic** deck carries every emission line exactly as ICRP-107 gives
it. It is the reference output and is always written. It can exceed 20 000
lines — the shipped example gives 29 244.

That is more discrete lines than a Monte-Carlo code can sample equitably, so the
**binned** deck collapses the DEP and DEE classes onto the fixed energy partition
in `data_EnergyBins/EnergyBins_PE.dat` (1063 bins, 1E-5 to 15 MeV). Each
non-empty bin becomes **one** line placed at the **yield-weighted mean energy** of
the emissions in it, carrying their **summed yield**; empty bins are dropped.

For the shipped example: 4275 photon lines become 592, and 22 604 electron lines
become 690. **Total yield is conserved exactly** and every other class stays
spectroscopically exact.

Use the binned deck for production transport; keep the spectroscopic deck as the
reference against which the binning can be justified.

### 4.3 Low-energy completeness

Electrons and photons are retained down to **1.0E-5 MeV (10 eV)**. The 2009
FORTRAN code discarded everything below 1 keV, which removed the Auger and
Coster-Kronig electrons entirely. Those matter for targeted radionuclide
therapy, where Auger emitters deposit their energy within nanometres of the DNA.

Set a `CUT:E` card in your MCNP model at around 15 eV to match; the deck writes
suggested cut-off cards as comments.

The detailed `ICRP-07.ACK` Auger spectrum **replaces** the coarse `AE` lines of
the RAD file where it exists — the two describe the same electrons, so adding
them would double-count.

---

## 5. Beta-minus and beta-plus (important)

ICRP-107 supplies **one** continuous beta spectrum per nuclide. For the **22
nuclides that emit both beta-minus and beta-plus**, that single tabulation is the
**sum of the two branches**. It cannot be used for either branch as it stands.

Version 2.0 reconstructs the two branches separately, taking the per-branch yield
and mean energy from `ICRP-07.RAD` (which tabulates them correctly) and the
positron spectrum from the project data file `ICRP-07-CORR.BEP`. Every one of the
952 beta branches in ICRP-107 now reproduces its ICRP-tabulated mean energy to
within 2 %.

**Versions before 2.0 gave both branches the combined yield and the combined
spectrum shape.** For Cu-64 the old deck wrote BEN = BEP = 5.649E-01 per decay,
where the truth is BEN = 3.900E-01 and BEP = 1.741E-01. Regenerate any deck
produced by an older version.

The 22 affected nuclides:

> Ag-106, Ag-108, As-74, Br-78, Br-80, Cl-36, Cs-130, Cs-132, Cu-64, Eu-150m,
> Eu-152, Eu-152m, I-126, I-128, In-112, In-114, K-40, Mn-54, Rb-84, Rh-102,
> Sb-122, Tm-168

The QA report names the method used for every branch:

| Method | Meaning |
|---|---|
| `BET` | ICRP-07.BET used directly (nuclide has one beta branch only) |
| `CORR` | evaluated positron spectrum from `ICRP-07-CORR.BEP` |
| `BET-CORR` | combined ICRP-07.BET minus that positron spectrum |
| `FERMI` | allowed-shape spectrum with its end-point fitted to the ICRP mean energy |

`FERMI` appears where the ICRP tabulation cannot supply the branch — six minor
branches, the largest 0.038 per decay. Yield and mean energy are still exact;
only the spectral shape is an approximation. The QA report lists them explicitly.

Full detail and evidence: **`Verification_Record_BEN_BEP.md`** in this folder.

### 5.1 Fission fragments (corrected in v2.1)

`ICRP-07.RAD` gives one `FF` line per spontaneously fissioning nuclide: the
number of fragments per decay, and the mean kinetic energy of **one** fragment.

Binary fission gives exactly one light and one heavy fragment, so each of the two
classes takes **half** the fragment yield, at its own energy. Fragment masses use
the near-constant actinide heavy peak at `A ~ 140`, sharing `A - nu` between the
two; charges follow `Z_f / A_f = Z / A` and sum to the parent `Z`; energies follow
from momentum conservation, `E_L / E_H = A_H / A_L`, with `TKE = 2 x` the
tabulated mean. The computed fragment energies reproduce measured values to about
2 % (Cf-252: 106.1 / 82.0 MeV computed, 105.7 / 80.2 measured).

**Versions before 2.1 gave both classes the full fragment yield, at the same
energy, with fixed stand-in ZAIDs.** A deck therefore emitted twice the physical
number of fragments and twice their kinetic energy. For Cm-250 that is 2.96
fragments per decay against a physical 1.48.

Fragment ZAIDs are derived from the dominant fissioning nuclide in the inventory
(the QA report names it), the same way the recoil-atom ZAID comes from the
dominant alpha emitter. A class with no emission writes `XX`.

One nuclide needs care: **Cf-252**'s `FF` yield in ICRP-107 counts *fissions*,
not fragments, unlike every other nuclide. The code detects this by testing the
prompt neutron multiplicity the yield implies (7.52 — impossible — against 3.760,
matching the published 3.7676) and corrects it. The QA report records the
correction.

Full detail and evidence: **`Verification_Record_Fission_Fragments.md`**.

---

## 6. Quality assurance

Every run verifies its own output and writes `QA_Report.txt` next to the
spectroscopic deck. Read the verdict line first.

| Verdict | Meaning |
|---|---|
| `PASS` | all checks satisfied |
| `WARN` | usable, but something rests on an approximation — the report says what |
| `FAIL` | **the deck must not be used** |

The checks:

| Group | What it verifies |
|---|---|
| `BETA-BRANCH` | every reconstructed branch reproduces the ICRP-07.RAD mean energy |
| `BETA-MIX` | the mixed BEN/BEP tabulations carry the activity-weighted ICRP yield and mean energy |
| `DUAL-EMITTER` | BEN and BEP are distinct — the signature of the historical defect |
| `FISSION` | one light and one heavy fragment per fission, energy and momentum conserved, implied neutron multiplicity physical |
| `SOURCE` | energies positive, yields non-negative, orderings correct |
| `PARTITION` | the energy partition is contiguous, ascending, strictly positive |
| `BINNING` | binning conserves yield and leaves other classes untouched |
| `DECK` | the written file is read back from disk and re-checked against the source term |

The report also carries a per-nuclide table of every beta branch — activity
fraction, yield, mean energy, reconstruction method and mean-energy error — and,
where the inventory fissions, a per-nuclide fission table giving the fission
rate, implied neutron multiplicity, TKE, both fragment ZAIDs and both fragment
energies. A reviewer can see exactly where each number came from.

`ALL` mode writes its own database-wide summary to
`Output_Spectroscopic_ALL/QA_Summary_ALL.txt`.

### 6.1 Regression suite

```
python -m tests.run_tests
```

49 tests, about 3 seconds. They are posed against ICRP-107's own tabulated
numbers rather than against previously generated output, so they stay valid if
the code is rewritten. Run them after any change to the code or to a data file.

### 6.2 Verifying the positron data file

```
python tools/build_corr_bep.py --check
```

Confirms that every spectrum in `ICRP-07-CORR.BEP` still reproduces the
`ICRP-07.RAD` beta-plus yield and mean energy. This needs no access to the
shared reference library.

---

## 7. The graph

`GraphBank/SDEF_IREST_ALL.{svg,pdf,png}` shows every particle class on one
log-log axis. SVG and PDF are pure vector (sharp at any zoom); PNG is 600 dpi.

* Energy axis **1E-5 to 15 MeV**, logarithmic.
* Yield axis floor **1E-12**.
* The DEP and DEE classes are drawn **binned**. The spectroscopic versions carry
  tens of thousands of overlapping lines that blur into an unreadable smear.
* All other classes are drawn spectroscopically exact.

| Class | Colour |
|---|---|
| Neutrons | green |
| DEP photons | blue (lollipops) |
| DEE electrons | violet (lollipops) |
| BEN beta-minus | light royal blue (continuous line) |
| BEP beta-plus | rouge (continuous line) |
| DEA alpha | brown (lollipops) |
| DER recoil atoms | black (lollipops) |
| Fission fragments | royal purple (lollipops) |

Fission fragments are emitted near 90 MeV, far outside the default window. They
are always plotted; widen the axis to see them:

```
python MCNP_SOURCE_APP.py --graph-energy-max 200
```

A caption beneath the legend names any class that carries yield but cannot appear
in the current window, and gives its yield and energy range — so a class with no
visible mark is still reported rather than silently absent. Note that a class may
be invisible because its yield is below the 1E-12 floor rather than because it is
off the energy axis; the caption distinguishes the two cases.

---

## 8. Project layout

```
MCNP_SOURCE_APP.py               driver
Run.cmd                          double-click launcher
msa_modules/                     the code
  elements.py                    periodic table, nuclide-name parsing
  databases.py                   ICRP-07.RAD / .BET / .ACK / -CORR.BEP parsers
  inventory.py                   Input.inp reader, availability check
  betaspectra.py                 BEN/BEP branch separation
  fission.py                     light/heavy fission-fragment partition
  sourceterm.py                  activity-weighted source-term construction
  energybins.py                  energy partition, DEP/DEE binning
  sdef.py                        MCNP6.3 SDEF writer
  qa.py                          self-verification and reporting
  plotting.py                    overview graph
tests/                           regression suite
tools/build_corr_bep.py          rebuild/verify ICRP-07-CORR.BEP
Input/Input.inp                  the activity inventory
data_Emissions_RadioNuclides/    ICRP-107 data files
data_EnergyBins/                 EnergyBins_PE.dat
Output_Spectroscopic_Case/       case deck + QA report
Output_Binned/                   binned case deck
Output_Spectroscopic_ALL/        per-nuclide library (on demand)
GraphBank/                       graphs
UserGuide/                       this guide, verification record
```

---

## 9. Data files

| File | Contents |
|---|---|
| `ICRP-07.RAD` | discrete radiations; **authoritative** for every branch yield and mean energy |
| `ICRP-07.BET` | continuous beta spectra, one per nuclide (see §5) |
| `ICRP-07.ACK` | detailed Auger / Coster-Kronig electron lines, energies in eV |
| `ICRP-07-CORR.BEP` | evaluated positron spectra for the 22 dual emitters — **project data, generated and verifiable via `tools/build_corr_bep.py`** |
| `EnergyBins_PE.dat` | the DEP/DEE energy partition, `Group LEB UEB` in MeV |

Each of these files carries its own header describing its format and, for
`ICRP-07-CORR.BEP`, its provenance and per-nuclide verification table.

---

## 10. Version history

| Version | Correction |
|---|---|
| 2.1 | Fission-fragment partition (§5.1) — fragment count was doubled, both fragments were at the average energy, ZAIDs were fixed stand-ins. Cf-252 yield anomaly detected. |
| 2.0 | Beta-minus / beta-plus separation (§5) — both branches carried the combined yield and the positron branch carried the electron spectrum. Cm-250 beta anomaly detected. |

Both corrections are documented in full, with evidence, in this folder:
`Verification_Record_BEN_BEP.md` and `Verification_Record_Fission_Fragments.md`.
