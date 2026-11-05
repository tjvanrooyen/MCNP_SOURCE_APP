# Verification Record — Correction of the Fission-Fragment Source Term

**Code:** MCNP_SOURCE_APP
**Version:** 2.1
**Date:** 2026-09-03
**Applies to:** all decks produced by versions before 2.1

Companion record: `Verification_Record_BEN_BEP.md` (the v2.0 beta correction).

---

## 1. Summary

MCNP_SOURCE_APP wrote **twice** the physical number of fission fragments, each
at the wrong energy, for every spontaneously fissioning nuclide. There are **35**
such nuclides in ICRP-107.

The heavy-fragment class (DEHFF, distribution 28) and the light-fragment class
(DELFF, distribution 29) were each given the **full** `ICRP-07.RAD` `FF` yield,
both at the **same** energy, with **fixed stand-in ZAIDs** unrelated to the
fissioning nuclide.

| Cm-250, per decay | Before v2.1 | Correct |
|---|---|---|
| DEHFF yield | 1.480 | **0.740** |
| DELFF yield | 1.480 | **0.740** |
| Fragments emitted | 2.960 | **1.480** |
| DEHFF energy | 90.62 MeV | **78.39 MeV** |
| DELFF energy | 90.62 MeV | **102.86 MeV** |
| Fragment kinetic energy | 268.2 MeV | **134.1 MeV** |
| DEHFF ZAID | 55137 (Cs-137) | **54140 (Xe-140)** |
| DELFF ZAID | 41099 (Mo-99) | **42107 (Mo-107)** |

The behaviour was inherited from the 2009 FORTRAN code, which also called
`totalYield(iso,'FF',yield)` for both classes.

---

## 2. What ICRP-07.RAD actually gives

One `FF` line per fissioning nuclide:

```
10 1.48000E+00 9.06223E+01 FF        <- Cm-250
```

Establishing what the two numbers mean is the whole of the correction, and both
readings are demonstrated from the data rather than assumed.

### 2.1 The energy is per fragment

If the tabulated energy is the mean kinetic energy of **one** fragment, then the
total fission kinetic energy is `TKE = 2 x Ebar`. Against published TKE values:

| Nuclide | 2 × Ebar | Published TKE | dev |
|---|---|---|---|
| Cf-252 | 188.1 MeV | 185.9 MeV | +1.2 % |
| Cm-250 | 181.2 MeV | 182.0 MeV | −0.4 % |
| Cf-254 | 187.6 MeV | 186.9 MeV | +0.4 % |
| U-238 | 169.7 MeV | 169.0 MeV | +0.4 % |
| Pu-240 | 176.4 MeV | 177.9 MeV | −0.9 % |
| Cm-244 | 182.7 MeV | 180.9 MeV | +1.0 % |

Six independent nuclides agree within 1.2 %. The energy is per fragment.

### 2.2 The yield counts fragments

If the tabulated yield counts **fragments**, the fission rate is `Y_FF / 2` and
the prompt neutron multiplicity implied by the RAD `N` line is
`nu = Y_n / (Y_FF / 2)`. That gives physical values — 2.13 to 4.01 — for 34 of
the 35 nuclides, matching published multiplicities:

| Nuclide | implied ν | published ν |
|---|---|---|
| Cf-254 | 3.89 | 3.89 |
| Fm-256 | 4.01 | 3.97 |
| Cm-250 | 3.31 | 3.31 |
| Cm-248 | 3.11 | 3.14 |
| Pu-244 | 2.30 | 2.30 |

---

## 3. The correction

Binary fission produces exactly **one light and one heavy fragment**, so each
class takes **half** the fragment yield. `msa_modules/fission.py` performs the
partition.

**Masses.** In actinide spontaneous fission the heavy-fragment peak sits at an
almost constant `A ≈ 140` while the light peak moves with the fissioning
nucleus. The fragments share `A − ν`, with `ν` taken from the RAD neutron line.

**Charges.** The unchanged-charge-distribution rule `Z_f / A_f = Z / A`, forced
to sum to the parent `Z`.

**Energies.** The fragments recoil against each other, so momentum conservation
gives `E_L / E_H = A_H / A_L`, with `E_L + E_H = TKE = 2 × Ebar`:

```
E_L = TKE . A_H / (A_L + A_H)        E_H = TKE . A_L / (A_L + A_H)
```

**ZAIDs.** Derived from the dominant fissioning nuclide in the inventory, the
same way the recoil-atom ZAID is derived from the dominant alpha emitter. A
class absent from the inventory writes `XX` rather than a misleading real ion.

### 3.1 Validation against measured fragment energies

| Nuclide | E_L computed | E_L measured | E_H computed | E_H measured |
|---|---|---|---|---|
| Cf-252 | 106.1 MeV | 105.7 MeV | 82.0 MeV | 80.2 MeV |
| Cm-250 | 102.9 MeV | 103 MeV | 78.4 MeV | 79 MeV |

Agreement is within about 2 %. The fragment ZAIDs are also chemically sensible:
Cf-252 gives Tc-108 and Cs-140, Cm-250 gives Mo-107 and Xe-140 — all
neutron-rich nuclides in the observed fission-product mass peaks.

---

## 4. A second fault found: the Cf-252 yield anomaly

**Cf-252 is the one nuclide whose `FF` yield is not twice its fission branch.**
ICRP-107 gives 3.0920E-02, which is the fission branch itself.

Read as fragments it implies `ν = 7.52` — roughly twice any multiplicity ever
measured, and impossible. Read as fissions it gives **3.760**, against the
published Cf-252 value of **3.7676**, one of the best-determined constants in
nuclear physics because Cf-252 is the standard neutron-source calibration.

`fission.fission_rate` therefore tests the implied multiplicity against the
physical range and falls back to the fissions-per-decay reading when the
fragments reading is unphysical. The correction is recorded in the QA report
under `FISSION` and in the notes.

**The two errors cancelled for Cf-252.** Because its yield was already per
fission, writing it to *both* classes happened to give the correct fragment
count (6.184E-02) and the correct total energy (5.816 MeV). Cf-252 was the only
nuclide whose totals were accidentally right — but both its fragments were still
placed at the average energy of 94.04 MeV instead of at 106.08 and 82.01 MeV.

Without the anomaly check, halving the yield as for every other nuclide would
have *introduced* a factor-2 error into Cf-252 — the most widely used
spontaneous-fission source there is.

---

## 5. Magnitude of the error that has been removed

Per-nuclide SDEF, per decay. "factor" is old ÷ corrected.

| Nuclide | old fragments | new fragments | factor | old FF energy | new FF energy | factor |
|---|---|---|---|---|---|---|
| Cf-254 | 3.988E+00 | 1.994E+00 | 2.00 | 374.1 MeV | 187.0 MeV | 2.00 |
| Fm-256 | 3.676E+00 | 1.838E+00 | 2.00 | 357.6 MeV | 178.8 MeV | 2.00 |
| Cm-250 | 2.960E+00 | 1.480E+00 | 2.00 | 268.2 MeV | 134.1 MeV | 2.00 |
| Cm-248 | 3.356E-01 | 1.678E-01 | 2.00 | 30.49 MeV | 15.25 MeV | 2.00 |
| **Cf-252** | 6.184E-02 | 6.184E-02 | **1.00** | 5.816 MeV | 5.816 MeV | **1.00** |
| Fm-257 | 8.400E-03 | 4.200E-03 | 2.00 | 0.8161 MeV | 0.4081 MeV | 2.00 |
| Pu-244 | 4.840E-03 | 2.420E-03 | 2.00 | 0.4245 MeV | 0.2123 MeV | 2.00 |
| Cf-250 | 3.080E-03 | 1.540E-03 | 2.00 | 0.2904 MeV | 0.1452 MeV | 2.00 |
| Fm-254 | 2.368E-03 | 1.184E-03 | 2.00 | 0.2309 MeV | 0.1155 MeV | 2.00 |
| Es-254m | 1.800E-03 | 9.000E-04 | 2.00 | 0.1722 MeV | 0.0861 MeV | 2.00 |
| Cm-246 | 1.052E-03 | 5.260E-04 | 2.00 | 0.0958 MeV | 0.0479 MeV | 2.00 |
| Es-255 | 1.800E-04 | 9.000E-05 | 2.00 | 0.0172 MeV | 0.0086 MeV | 2.00 |
| Cf-248 | 1.160E-04 | 5.800E-05 | 2.00 | 0.0110 MeV | 0.0055 MeV | 2.00 |
| Fm-252 | 9.200E-05 | 4.600E-05 | 2.00 | 0.0090 MeV | 0.0045 MeV | 2.00 |

The remaining 21 nuclides follow the same factor of 2; their fission-fragment
yields are below 1E-5 per decay.

**Every nuclide except Cf-252 had twice the physical fragment count and twice
the fission-fragment energy deposition.** All 35, Cf-252 included, had both
fragments at the average energy instead of at the light and heavy energies, and
all 35 carried ZAIDs (Cs-137 / Mo-99) unrelated to the fissioning nuclide.

### 5.1 Who is affected

Fission fragments stop within a few micrometres, so they matter for **local
energy deposition** — self-heating of separated actinides, source-capsule
integrity, alpha/heavy-ion damage studies, and any tally scoring heavy-ion
energy deposition in or immediately around the source region. They do not affect
photon or neutron dose at distance, where the neutron and gamma classes dominate.

The exposure is therefore concentrated in calculations involving **Cf-254,
Fm-256, Cm-250 or Cm-248**, whose fission branches are large. For Cm-250 the
fission-fragment energy deposited per decay was over-stated by 134 MeV — more
than the total energy actually released as fragments.

---

## 6. Standing verification

The `FISSION` check group runs on every execution:

| Check | Verifies |
|---|---|
| `LFF yield == fission rate` | one light fragment per fission |
| `HFF yield == fission rate` | one heavy fragment per fission |
| `LFF + HFF == 2 x fission rate` | no fragment double-count |
| `fragment kinetic energy == fissions x TKE` | energy conservation |
| `light fragment faster than heavy` | momentum conservation |
| `implied neutron multiplicity is physical` | identifies yield-convention anomalies |
| `fragment ZAIDs conserve parent Z and A - nu` | the fragments really are that nucleus |

The QA report carries a per-nuclide fission table (fission rate, ν, TKE, both
ZAIDs, both energies) and records any yield-convention correction in its notes.
`ALL` mode writes the same census for the whole database.

Regression tests (in `tests/test_source_term.py`, class `TestFissionFragments`):

* `test_each_class_carries_the_fission_rate_not_the_fragment_count` — all 35
* `test_cm250_fragment_yield_is_not_doubled`
* `test_fragment_kinetic_energy_is_conserved` — all 35
* `test_light_fragment_is_faster_and_lighter` — all 35
* `test_momentum_conservation` — all 35
* `test_zaids_conserve_charge_and_mass` — all 35
* `test_tke_matches_published_values`
* `test_fragment_energies_match_published_values`
* `test_implied_neutron_multiplicity_is_physical_everywhere`
* `test_cf252_yield_anomaly_is_detected_and_corrected`
* `test_only_cf252_needs_the_yield_correction`
* `test_deck_writes_a_zaid_placeholder_when_there_is_no_fission`

---

## 7. Actions for the user

1. **Regenerate any deck** produced by a version before 2.1 whose inventory
   contains a spontaneously fissioning nuclide. The per-nuclide library in
   `Output_Spectroscopic_ALL/` has already been regenerated (1252 files).
2. **Review results** from pre-2.1 decks that scored heavy-ion energy
   deposition, most urgently for Cf-254, Fm-256, Cm-250 and Cm-248.
3. Note that the fragment ZAIDs are now nuclide-derived. A deck whose inventory
   contains more than one fissioning nuclide takes its ZAIDs from the dominant
   one, exactly as the recoil-atom class does; the QA report names it.

---

## 8. Modelling assumptions, stated plainly

These are approximations, and they are deliberate:

* **Two discrete fragment energies, not a distribution.** Real fission gives a
  broad, correlated mass-and-energy distribution. ICRP-107 supplies one mean
  energy, so one light and one heavy line is all the data supports. The totals —
  fragment count and kinetic energy — are exact.
* **Heavy-fragment peak fixed at A ≈ 140.** True for actinide spontaneous
  fission to within a few mass units. The value is a named constant in
  `fission.py`, clamped so the heavy fragment can never come out lighter than
  the light one.
* **Charge follows unchanged charge distribution.** A standard first
  approximation; real fragments deviate by roughly half a charge unit.
* **Prompt neutrons are already accounted for** in the fragment masses via ν,
  and are emitted separately by the neutron class (distribution 21). No
  double-counting: the neutron class comes from the RAD `N` line, which is
  independent of the `FF` line.
