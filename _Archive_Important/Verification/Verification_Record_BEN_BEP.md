# Verification Record — Correction of the BEN / BEP Beta Branch Defect

**Code:** MCNP_SOURCE_APP
**Version:** 2.0
**Date:** 2026-09-03
**Applies to:** all decks produced by versions before 2.0

---

## 1. Summary

MCNP_SOURCE_APP gave wrong beta-electron (BEN) and beta-positron (BEP) source
terms for every radionuclide that emits **both** beta-minus and beta-plus. There
are **22** such nuclides in ICRP-107. Four of them — **Cl-36, Eu-152, Ag-108 and
Mn-54** — appear in the shipped example inventory, so the defect affected the
current case as well as the per-nuclide library.

Both the yield and the spectrum shape were wrong:

| Quantity | Before v2.0 | Correct (ICRP-07.RAD) |
|---|---|---|
| Cu-64 BEN yield | 5.649E-01 /decay | **3.900E-01** /decay |
| Cu-64 BEP yield | 5.649E-01 /decay | **1.741E-01** /decay |
| Cu-64 BEN mean energy | 0.2174 MeV | **0.1904 MeV** |
| Cu-64 BEP mean energy | 0.2174 MeV | **0.2782 MeV** |

The defect is corrected in version 2.0. This record states the cause, the
correction, and the evidence that the correction is right.

---

## 2. Cause

ICRP-107 (2008) supplies **one** continuous beta spectrum per nuclide, in
`ICRP-07.BET`. For a nuclide with both branches that single tabulation is the
**sum** of the two:

```
∫ (ICRP-07.BET pdf) dE  =  Y_BEN + Y_BEP
```

This was confirmed for all 22 nuclides. For Cu-64 the integral is 5.649E-01,
against `ICRP-07.RAD` branch yields of 3.900E-01 (B−) and 1.741E-01 (B+),
which sum to 5.641E-01 — agreement to the 0.15 % quadrature noise of the
ICRP tabulation.

The code handed that single tabulation to **both** the BEN and the BEP
distribution. Two independent errors followed:

1. **Yield.** Each branch received the *combined* yield. In a per-nuclide deck
   that doubles the total beta emission and makes the split between the two
   branches arbitrary. In a mixture the error depends on the branching ratios
   and on how much of the activity the dual emitters carry — see §5.1.

2. **Spectrum shape.** The positron branch received the shape of the combined —
   in most cases electron-dominated — spectrum. This is physically wrong, not
   merely imprecise: the daughter nucleus **attracts** the emitted electron and
   **repels** the emitted positron, so the two spectra differ systematically,
   with the positron spectrum displaced to higher energy. For Cu-64 the true
   mean energies differ by 46 %.

The relative size of the yield error depends on the branching ratio and is very
large where the two branches are unequal — see §5.

---

## 3. Correction

`msa_modules/betaspectra.py` reconstructs the two branches separately. The
authority for the per-branch yield and mean energy is `ICRP-07.RAD`, which
tabulates them correctly and separately:

```
 5 3.90000E-01 1.90428E-01 B-       <- Y_BEN, Ebar_BEN
 4 1.74083E-01 2.78220E-01 B+       <- Y_BEP, Ebar_BEP
```

Branch shapes are established in this order of preference.

| Method | Applies to | Count |
|---|---|---|
| `BET` | single-branch nuclide: the ICRP-07.BET tabulation *is* that branch | 907 |
| `CORR` | positron branch of a dual emitter: evaluated spectrum from `ICRP-07-CORR.BEP` | 22 |
| `BET-CORR` | electron branch of a dual emitter: combined BET minus the positron spectrum | 17 |
| `FERMI` | fallback: allowed-shape spectrum with the end-point fitted to `Ebar` | 6 |

Every branch is finally rescaled so that its integral equals the `ICRP-07.RAD`
yield exactly.

### 3.1 The positron spectra

`data_Emissions_RadioNuclides/ICRP-07-CORR.BEP` supplies the separately
evaluated positron spectrum for each of the 22 dual emitters. It is generated
from the institutional reference emission library
`C:\Py_Data_Shared\Source_Emissions\<nuclide>.xlsx` (sheet `BEP`) by
`tools/build_corr_bep.py`, and is committed to the project so that ordinary runs
never depend on the shared folder.

### 3.2 Where subtraction is not usable

`BET − CORR` is exact in principle but is a difference of two nearly equal
quantities. Where `Y_BEN << Y_BEP` the residual is smaller than the quadrature
noise of the tabulations and its shape is meaningless. The reconstruction is
therefore *tested* against `Ebar_BEN` and rejected if it fails, falling back to
an allowed-shape spectrum:

```
N(E) dE  ∝  F(Z_d, E) · p · W · (E0 − E)²  dE
```

with the Fermi function `F` for the daughter charge `Z_d` and the end-point `E0`
solved by bisection so the mean energy equals `Ebar` exactly.

Six branches use this fallback. All are minor branches; the largest is 0.038 per
decay:

| Branch | Y /decay | Reason |
|---|---|---|
| Br-78 BEN | 1.0E-04 | Y_BEN / Y_BEP = 1.1E-04, subtraction ill-conditioned |
| Ag-106 BEN | 1.0E-02 | Y_BEN / Y_BEP = 1.7E-02 |
| Cs-130 BEN | 1.6E-02 | Y_BEN / Y_BEP = 3.7E-02 |
| Rb-84 BEN | 3.8E-02 | Y_BEN / Y_BEP = 1.4E-01 |
| Tm-168 BEN | 1.0E-04 | Y_BEN / Y_BEP = 2.6E+00, but both branches minute |
| Cm-250 BEN | 8.0E-02 | ICRP-07.BET block is not this nuclide's spectrum — see §6 |

---

## 4. Evidence

Four independent lines of evidence support the correction.

### 4.1 Agreement with the independent legacy implementation

The archived 2025 output in `Output_Spectroscopic_ALL/Archive_2025_values/` was
produced by the earlier FORTRAN `mcnpsource` code, which read the branch yields
from `ICRP-07.RAD`. The corrected v2.0 yields reproduce those legacy yields
**exactly (to all 8 written digits) for all 22 dual emitters**. The regression
was introduced by the Python port, and v2.0 restores agreement with the
independent implementation.

The legacy code did *not*, however, separate the spectrum *shapes*. Across all
22 dual emitters its BEN and BEP tabulations agree with each other to within
**0.71 %**, while the true branch mean energies differ by up to **1473 %**
(Ag-106). For Ag-108 the legacy means are 0.62217 MeV (BEN) and 0.62234 MeV
(BEP), against true values of 0.6235 and 0.3972; version 2.0 gives 0.62425 and
0.39723. The spectrum defect was therefore present in **both** the legacy
FORTRAN code and the Python port, and is corrected here for the first time.

### 4.2 Reproduction of the ICRP-107 mean energies

`ICRP-07.RAD` tabulates the mean energy of each beta branch independently of any
spectrum. This is a genuine external check, and it is the check that would have
caught the defect: a positron branch carrying the electron spectrum fails it by
tens of percent.

**All 952 beta branches in ICRP-107 reproduce their tabulated `Ebar` to within
2 %**, worst case 1.51 % (Es-254m BEN, a known ICRP-107 tabulation blemish).

Before the correction the dual emitters failed this check badly, because both
branches were given the mean energy of the combined spectrum:

| | worst error before v2.0 | worst error in v2.0 |
|---|---|---|
| BEP branch mean energy | **+612 %** (I-128) | +0.00 % |
| BEN branch mean energy | **+1448 %** (Ag-106) | +0.79 % (In-112) |

### 4.3 Double validation of the positron data

Each spectrum in `ICRP-07-CORR.BEP` satisfies **two** constraints taken from
ICRP-107 itself — its integral equals `Y_BEP` and its first moment equals
`Ebar_BEP` — to better than 0.001 %. Verify at any time:

```
python tools/build_corr_bep.py --check
```

Independently, the positron spectra agree with an allowed-shape Fermi spectrum
for the daughter charge `Z − 1` to within ±5 % point-by-point, the residual being
the expected difference between a screened relativistic Fermi function and the
non-relativistic approximation.

### 4.4 Conservation of the total beta emission

For every dual emitter, `Y_BEN + Y_BEP` from the corrected reconstruction equals
the integral of the original combined `ICRP-07.BET` tabulation to within 0.5 %.
Nothing has been created or lost — the combined emission has been *partitioned*
correctly instead of duplicated.

---

## 5. Magnitude of the error that has been removed

Per-nuclide SDEF, PAR yields on the `SP7` card. "factor" is old ÷ corrected.

| Nuclide | old BEN | new BEN | factor | old BEP | new BEP | factor |
|---|---|---|---|---|---|---|
| Br-78 | 9.253E-01 | 1.000E-04 | 9253 | 9.253E-01 | 9.245E-01 | 1.0 |
| Ag-106 | 5.965E-01 | 1.000E-02 | 59.6 | 5.965E-01 | 5.856E-01 | 1.0 |
| Cs-130 | 4.523E-01 | 1.600E-02 | 28.3 | 4.523E-01 | 4.355E-01 | 1.0 |
| Rb-84 | 3.090E-01 | 3.800E-02 | 8.1 | 3.090E-01 | 2.710E-01 | 1.1 |
| As-74 | 6.311E-01 | 3.400E-01 | 1.9 | 6.311E-01 | 2.906E-01 | 2.2 |
| Rh-102 | 3.743E-01 | 2.200E-01 | 1.7 | 3.743E-01 | 1.539E-01 | 2.4 |
| In-112 | 6.554E-01 | 4.400E-01 | 1.5 | 6.554E-01 | 2.149E-01 | 3.1 |
| Cu-64 | 5.649E-01 | 3.900E-01 | 1.4 | 5.649E-01 | 1.741E-01 | 3.2 |
| Tm-168 | 1.388E-04 | 1.000E-04 | 1.4 | 1.388E-04 | 3.856E-05 | 3.6 |
| Cs-132 | 2.303E-02 | 1.870E-02 | 1.2 | 2.303E-02 | 4.303E-03 | 5.4 |
| Br-80 | 9.410E-01 | 9.170E-01 | 1.0 | 9.410E-01 | 2.209E-02 | 42.6 |
| I-126 | 4.833E-01 | 4.730E-01 | 1.0 | 4.833E-01 | 1.007E-02 | 48.0 |
| Eu-150m | 8.975E-01 | 8.900E-01 | 1.0 | 8.975E-01 | 6.515E-03 | 138 |
| Ag-108 | 9.746E-01 | 9.715E-01 | 1.0 | 9.746E-01 | 2.869E-03 | 340 |
| Mn-54 | 2.905E-06 | 2.900E-06 | 1.0 | 2.905E-06 | 5.600E-09 | 519 |
| K-40 | 8.923E-01 | 8.914E-01 | 1.0 | 8.923E-01 | 9.000E-04 | 991 |
| Eu-152 | 2.794E-01 | 2.790E-01 | 1.0 | 2.794E-01 | 1.381E-04 | 2023 |
| Cl-36 | 9.813E-01 | 9.810E-01 | 1.0 | 9.813E-01 | 1.390E-04 | 7061 |
| Eu-152m | 7.211E-01 | 7.200E-01 | 1.0 | 7.211E-01 | 7.496E-05 | 9619 |
| Sb-122 | 9.770E-01 | 9.759E-01 | 1.0 | 9.770E-01 | 6.258E-05 | 15613 |
| In-114 | 9.970E-01 | 9.950E-01 | 1.0 | 9.970E-01 | 3.253E-05 | 30645 |
| I-128 | 9.322E-01 | 9.310E-01 | 1.0 | 9.322E-01 | 2.484E-05 | 37529 |

The total beta emission `BEN + BEP` was over-stated by a factor of about **2**
for every one of these nuclides.

### 5.1 Effect on the shipped example inventory

98 nuclides, 1.680E+11 Bq, of which Cl-36, Eu-152, Ag-108 and Mn-54 are dual
emitters. Per average decay of the mixture:

| Quantity | Before v2.0 | v2.0 | factor |
|---|---|---|---|
| BEN yield | 7.005E-01 | 6.999E-01 | 1.001 |
| BEP yield | 3.850E-06 | **8.824E-08** | **43.6** |
| BEN mean energy | 7.568E-02 MeV | 7.569E-02 MeV | 1.000 |
| BEP mean energy | 2.809E-01 MeV | **1.104E-01 MeV** | **2.55** |

The positron source was over-stated **43.6-fold in yield and 2.5-fold in mean
energy**, so the positron *energy* release was over-stated by about 110×. BEN is
essentially unchanged, because the mixture's beta-minus emission is dominated by
pure beta-minus nuclides (Co-60, Ni-63) and because the dual emitters present
are themselves beta-minus dominated.

Positron transport carries 1.022 MeV of annihilation photons per positron, so an
over-stated BEP yield propagates into gamma dose. Any dose or shielding result
produced from a pre-2.0 deck containing one of the 22 nuclides should be
reviewed — most urgently one in which a dual emitter dominates the inventory,
where the per-nuclide factors of §5 apply directly rather than being diluted.

---

## 6. Related defect found and corrected: Cm-250

The verification sweep surfaced a second, unrelated fault.

`ICRP-07.BET` for **Cm-250** does not contain the Cm-250 beta spectrum. It
contains the aggregated **fission-product** beta spectrum from Cm-250's 74 %
spontaneous-fission branch: it integrates to **5.59 betas/decay with a mean of
1.51 MeV**, whereas Cm-250's own B− branch is **0.08 betas/decay with a mean of
9.66 keV**.

Versions before 2.0 took that block at face value and over-stated the Cm-250
beta source by a factor of about **70**, at a mean energy 156× too high. The
mean-energy check now rejects it and an allowed-shape spectrum matched to the
`ICRP-07.RAD` branch is used instead.

The same `ICRP-07.BET` pattern occurs for the other spontaneously fissioning
nuclides (Cf-252, Cf-254, Fm-256, Cm-248 and others), but those have no `B−`
line in `ICRP-07.RAD`, so no beta source was ever built for them and no output
was affected.

---

## 7. Standing verification

Every run of MCNP_SOURCE_APP now verifies its own output and writes
`Output_Spectroscopic_Case/QA_Report.txt`. The process exit status is 3 if the
verdict is FAIL. Checks cover the beta reconstruction against ICRP-107, the
energy-partition integrity, yield conservation under binning, and a read-back of
the written deck from disk.

The regression suite pins this correction:

```
python -m tests.run_tests          # 33 tests
```

The tests that would fail if the defect returned:

* `test_dual_emitter_yields_match_icrp_rad`
* `test_dual_emitter_branches_are_never_identical`
* `test_dual_emitter_total_beta_is_not_doubled`
* `test_positron_spectrum_is_not_the_electron_spectrum`
* `test_cu64_reference_values`
* `test_every_branch_reproduces_the_icrp_mean_energy` (all 952 branches)
* `test_cm250_does_not_inherit_the_fission_product_beta_spectrum`

---

## 8. Actions for the user

1. **Regenerate any deck** produced by a version before 2.0. The per-nuclide
   library in `Output_Spectroscopic_ALL/` has already been regenerated
   (1252 files, 2026-09-03).
2. **Review results** from pre-2.0 decks whose inventory contains any of the 22
   nuclides in §5, or Cm-250.
3. Keep `ICRP-07-CORR.BEP` with the code. Without it the dual emitters fall back
   to the allowed-shape approximation — yields and mean energies stay correct,
   spectral shape is approximated — and the run prints a warning.

---

## 9. Related item — corrected separately in v2.1

The verification sweep for this correction also surfaced a fault in the
fission-fragment classes: DEHFF (dist. 28) and DELFF (dist. 29) were each given
the **full** `ICRP-07.RAD` `FF` yield, so a deck emitted two fragments per
fragment listed — for Cm-250, 2.96 fragments/decay against a physical 1.48.

That was outside the scope of the beta correction and is dealt with separately
in **`Verification_Record_Fission_Fragments.md`** (version 2.1), which also
documents a yield-convention anomaly in ICRP-107's Cf-252 entry.
