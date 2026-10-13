"""Regression tests for MCNP_SOURCE_APP.

Run from the project root::

    python -m tests.run_tests          # or: python -m unittest discover -s tests

The tests are posed against ICRP-107's own tabulated numbers, not against
previously generated output, so they stay valid if the code is rewritten.  The
group that matters most is :class:`TestBetaBranchSeparation`: it pins the defect
that motivated version 2.0 -- a nuclide emitting both beta-minus and beta-plus
being given the same yield and the same spectrum for both branches.
"""

from __future__ import annotations

import os
import unittest

import numpy as np

from msa_modules import betaspectra, energybins, fission, qa, sdef, sourceterm
from msa_modules import databases as db
from msa_modules.elements import canonical_key, parse_nuclide
from msa_modules.inventory import Inventory, InventoryItem, read_inventory

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EMISSIONS = os.path.join(ROOT, "data_Emissions_RadioNuclides")
RAD_PATH = os.path.join(EMISSIONS, "ICRP-07.RAD")
BET_PATH = os.path.join(EMISSIONS, "ICRP-07.BET")
ACK_PATH = os.path.join(EMISSIONS, "ICRP-07.ACK")
CORR_PATH = os.path.join(EMISSIONS, "ICRP-07-CORR.BEP")
BINS_PATH = os.path.join(ROOT, "data_EnergyBins", "EnergyBins_PE.dat")

# Loaded once: the RAD file alone is 14 MB.
_DB: dict = {}


def setUpModule() -> None:
    _DB["rad"] = db.load_rad(RAD_PATH, None)
    _DB["bet"] = db.load_bet(BET_PATH, None)
    _DB["ack"] = db.load_ack(ACK_PATH, None)
    _DB["corr"] = db.load_bet(CORR_PATH, None)
    _DB["dual"] = sorted(
        k for k, iso in _DB["rad"].items()
        if betaspectra.rad_branch(iso, "B-")[0] > 0.0
        and betaspectra.rad_branch(iso, "B+")[0] > 0.0
    )


def _single(key: str, points: int = 500) -> sourceterm.SourceTerm:
    """Source term for one nuclide at unit activity, i.e. emissions per decay."""
    riso = _DB["rad"][key]
    item = InventoryItem(raw_name=riso.name, activity=1.0, nuclide=parse_nuclide(riso.name))
    inv = Inventory(num_beta_points=points, items=[item])
    return sourceterm.build_source_term(
        inv, {key: riso}, _DB["bet"], _DB["ack"], _DB["corr"]
    )


class TestBetaBranchSeparation(unittest.TestCase):
    """The defect corrected in v2.0: BEN and BEP must be separate branches."""

    def test_database_has_the_expected_dual_emitters(self):
        self.assertEqual(len(_DB["dual"]), 22)
        self.assertIn("cu64", _DB["dual"])

    def test_corrected_positron_file_covers_every_dual_emitter(self):
        """Without this file the dual emitters fall back to an approximation."""
        missing = [k for k in _DB["dual"] if k not in _DB["corr"]]
        self.assertEqual(missing, [], f"ICRP-07-CORR.BEP is missing {missing}")

    def test_dual_emitter_yields_match_icrp_rad(self):
        """Each branch carries its OWN yield -- not the combined BET integral."""
        for key in _DB["dual"]:
            riso = _DB["rad"][key]
            y_ben, _ = betaspectra.rad_branch(riso, "B-")
            y_bep, _ = betaspectra.rad_branch(riso, "B+")
            with self.subTest(nuclide=riso.name):
                st = _single(key)
                self.assertAlmostEqual(st.ben.total / y_ben, 1.0, places=9)
                self.assertAlmostEqual(st.bep.total / y_bep, 1.0, places=9)

    def test_dual_emitter_branches_are_never_identical(self):
        """The exact signature of the historical defect: BEN == BEP."""
        for key in _DB["dual"]:
            riso = _DB["rad"][key]
            with self.subTest(nuclide=riso.name):
                st = _single(key)
                self.assertNotAlmostEqual(
                    st.ben.total, st.bep.total, places=12,
                    msg=f"{riso.name}: BEN and BEP yields are identical",
                )
                self.assertNotAlmostEqual(
                    st.ben.mean, st.bep.mean, places=12,
                    msg=f"{riso.name}: BEN and BEP mean energies are identical",
                )

    def test_dual_emitter_total_beta_is_not_doubled(self):
        """BEN + BEP must equal the combined ICRP-07.BET integral, not twice it."""
        for key in _DB["dual"]:
            riso = _DB["rad"][key]
            biso = _DB["bet"][key]
            combined = float(np.trapezoid(biso.pdf, biso.energy))
            with self.subTest(nuclide=riso.name):
                st = _single(key)
                self.assertAlmostEqual((st.ben.total + st.bep.total) / combined, 1.0,
                                       delta=0.005)

    def test_cu64_reference_values(self):
        """Worked reference case, straight from the ICRP-07.RAD B-/B+ lines."""
        st = _single("cu64")
        self.assertAlmostEqual(st.ben.total, 3.90000e-01, places=8)
        self.assertAlmostEqual(st.bep.total, 1.74083e-01, places=8)
        self.assertAlmostEqual(st.ben.mean, 1.904280e-01, delta=0.02 * 1.904280e-01)
        self.assertAlmostEqual(st.bep.mean, 2.782200e-01, delta=0.02 * 2.782200e-01)
        # The positron spectrum is harder than the electron spectrum.
        self.assertGreater(st.bep.mean, st.ben.mean)

    def test_positron_spectrum_is_not_the_electron_spectrum(self):
        """Shape check: for a dual emitter the two tabulations must differ."""
        st = _single("cu64", points=900)
        ben = np.asarray(st.ben.pdf)
        bep = np.asarray(st.bep.pdf)
        # Compare normalised shapes on their own grids' overlap.
        grid = np.asarray(st.ben.energy)
        bep_i = np.interp(grid, st.bep.energy, bep, left=0.0, right=0.0)
        shape_ben = ben / max(ben.max(), 1e-300)
        shape_bep = bep_i / max(bep_i.max(), 1e-300)
        self.assertGreater(float(np.abs(shape_ben - shape_bep).max()), 0.1)


class TestWholeDatabaseBeta(unittest.TestCase):
    """Sweep every beta branch in ICRP-107."""

    def test_every_branch_reproduces_the_icrp_mean_energy(self):
        offenders = []
        for key, riso in _DB["rad"].items():
            branches = betaspectra.decompose(
                parse_nuclide(riso.name), riso, _DB["bet"].get(key), _DB["corr"].get(key)
            )
            for b in (branches.ben, branches.bep):
                if b is not None and abs(b.mean_dev) > betaspectra.MEAN_ENERGY_TOL:
                    offenders.append(f"{riso.name} {b.particle} {b.mean_dev:+.2%} ({b.origin})")
        self.assertEqual(offenders, [], f"{len(offenders)} branch(es) off tolerance")

    def test_branch_count_and_method_census(self):
        """Guards against a data file silently going missing or changing."""
        census: dict[str, int] = {}
        for key, riso in _DB["rad"].items():
            branches = betaspectra.decompose(
                parse_nuclide(riso.name), riso, _DB["bet"].get(key), _DB["corr"].get(key)
            )
            for b in (branches.ben, branches.bep):
                if b is not None:
                    census[b.origin] = census.get(b.origin, 0) + 1
        self.assertEqual(sum(census.values()), 952)
        self.assertEqual(census.get(betaspectra.ORIGIN_CORR), 22)
        self.assertEqual(census.get(betaspectra.ORIGIN_SUBTRACT), 17)
        self.assertEqual(census.get(betaspectra.ORIGIN_FERMI), 6)

    def test_single_branch_nuclides_populate_only_their_own_branch(self):
        for key, expect_ben in (("cs137", True), ("na22", False), ("y90", True)):
            riso = _DB["rad"][key]
            with self.subTest(nuclide=riso.name):
                st = _single(key)
                if expect_ben:
                    self.assertGreater(st.ben.total, 0.0)
                    self.assertEqual(st.bep.total, 0.0)
                else:
                    self.assertGreater(st.bep.total, 0.0)
                    self.assertEqual(st.ben.total, 0.0)


class TestBetSpectrumAnomalies(unittest.TestCase):
    """ICRP-07.BET blocks that do not describe the nuclide's own beta branch."""

    def test_cm250_does_not_inherit_the_fission_product_beta_spectrum(self):
        """Cm-250's BET block is the beta spectrum of its fission products.

        It integrates to 5.59 betas/decay with a mean of 1.51 MeV, while Cm-250's
        own B- branch is 0.08 betas/decay with a mean of 9.66 keV.  Taking the
        BET block at face value over-states the beta source ~70-fold.
        """
        riso = _DB["rad"]["cm250"]
        y_ben, ebar = betaspectra.rad_branch(riso, "B-")
        st = _single("cm250")
        self.assertAlmostEqual(st.ben.total / y_ben, 1.0, places=9)
        self.assertAlmostEqual(st.ben.mean / ebar, 1.0, delta=betaspectra.MEAN_ENERGY_TOL)
        self.assertLess(st.ben.total, 0.1)


class TestFermiFallback(unittest.TestCase):
    """The allowed-shape spectrum used when a tabulation cannot be trusted."""

    def test_endpoint_solver_reproduces_the_requested_mean(self):
        for z, ebar, positron in ((29, 0.19043, False), (28, 0.27822, True),
                                  (97, 0.00966, False), (20, 2.5, False)):
            with self.subTest(z=z, ebar=ebar, positron=positron):
                e0 = betaspectra.endpoint_for_mean(z, ebar, positron)
                self.assertGreater(e0, ebar)
                grid, pdf = betaspectra.fermi_spectrum(z, 1.0, ebar, positron)
                area = np.trapezoid(pdf, grid)
                mean = np.trapezoid(pdf * grid, grid) / area
                self.assertAlmostEqual(mean / ebar, 1.0, places=4)

    def test_fermi_function_signs(self):
        """Coulomb attraction enhances electrons; repulsion suppresses positrons."""
        e = np.array([0.05, 0.2, 0.5])
        f_minus = betaspectra.fermi_function(30, e, positron=False)
        f_plus = betaspectra.fermi_function(30, e, positron=True)
        self.assertTrue(np.all(f_minus > 1.0))
        self.assertTrue(np.all(f_plus < 1.0))
        self.assertTrue(np.all(f_minus > f_plus))

    def test_spectrum_normalisation_and_support(self):
        grid, pdf = betaspectra.fermi_spectrum(29, 0.25, 0.19043, positron=False)
        self.assertTrue(np.all(pdf >= 0.0))
        self.assertEqual(pdf[-1], 0.0)


class TestMixture(unittest.TestCase):
    """Activity weighting of a multi-nuclide inventory."""

    def _mix(self, keys_and_activities, points=900):
        items = [
            InventoryItem(raw_name=_DB["rad"][k].name, activity=a,
                          nuclide=parse_nuclide(_DB["rad"][k].name))
            for k, a in keys_and_activities
        ]
        inv = Inventory(num_beta_points=points, items=items)
        rad = {k: _DB["rad"][k] for k, _ in keys_and_activities}
        return sourceterm.build_source_term(inv, rad, _DB["bet"], _DB["ack"], _DB["corr"])

    def test_yields_are_activity_weighted(self):
        st = self._mix([("cu64", 3.0), ("cs137", 1.0)])
        y_ben_cu, _ = betaspectra.rad_branch(_DB["rad"]["cu64"], "B-")
        y_bep_cu, _ = betaspectra.rad_branch(_DB["rad"]["cu64"], "B+")
        y_ben_cs, _ = betaspectra.rad_branch(_DB["rad"]["cs137"], "B-")
        self.assertAlmostEqual(st.bep.total, 0.75 * y_bep_cu, places=8)
        self.assertAlmostEqual(st.ben.total, 0.75 * y_ben_cu + 0.25 * y_ben_cs, places=8)

    def test_tabulated_integral_equals_the_declared_yield(self):
        """The deck must be internally consistent: SP7 entry == integral of SI/SP."""
        st = self._mix([("cu64", 1.0), ("na22", 1.0), ("y90", 1.0)])
        for dist in (st.ben, st.bep):
            self.assertAlmostEqual(dist.total / dist.expected_total, 1.0, places=9)

    def test_mixture_mean_energy_tracks_icrp(self):
        st = self._mix([("cu64", 1.0), ("na22", 1.0), ("cl36", 1.0)])
        for dist in (st.ben, st.bep):
            self.assertAlmostEqual(dist.mean / dist.expected_mean, 1.0,
                                   delta=betaspectra.MEAN_ENERGY_TOL)

    def test_ack_replaces_rad_auger_lines_without_double_counting(self):
        """The ACK breakdown and the RAD AE lines describe the same electrons."""
        key = "cd109"
        riso = _DB["rad"][key]
        rad_ae = sum(e.yield_ for e in riso.entries if e.ptype == "AE")
        ack_total = sum(_DB["ack"][key].yield_)
        st = _single(key)
        rad_ie = sum(e.yield_ for e in riso.entries
                     if e.ptype == "IE" and e.energy >= sourceterm.INCLUDE_THRESHOLD_MEV)
        self.assertAlmostEqual(st.dee.total, ack_total + rad_ie, delta=0.02 * st.dee.total)
        self.assertLess(st.dee.total, rad_ae + ack_total + rad_ie)


class TestFissionFragments(unittest.TestCase):
    """Binary fission gives ONE light and ONE heavy fragment, not two of each."""

    # Published total fission kinetic energy and fragment energies, MeV.
    MEASURED_TKE = {"cf252": 185.9, "cm250": 182.0, "cf254": 186.9,
                    "u238": 169.0, "pu240": 177.9, "cm244": 180.9}
    MEASURED_FRAGMENT_E = {"cf252": (105.7, 80.2), "cm250": (103.0, 79.0)}

    def _fissioners(self):
        return [k for k, iso in _DB["rad"].items()
                if fission.fission_yield(iso)[0] > 0.0]

    def test_database_fissioning_nuclide_count(self):
        self.assertEqual(len(self._fissioners()), 35)

    def test_each_class_carries_the_fission_rate_not_the_fragment_count(self):
        """The historical defect: both classes got the full FF yield."""
        for key in self._fissioners():
            riso = _DB["rad"][key]
            with self.subTest(nuclide=riso.name):
                st = _single(key)
                s = fission.split(parse_nuclide(riso.name), riso)
                self.assertAlmostEqual(st.lff.total / s.fissions, 1.0, places=9)
                self.assertAlmostEqual(st.hff.total / s.fissions, 1.0, places=9)
                self.assertAlmostEqual(
                    (st.lff.total + st.hff.total) / (2.0 * s.fissions), 1.0, places=9)

    def test_cm250_fragment_yield_is_not_doubled(self):
        """Reference case: 1.48 fragments/decay, i.e. 0.74 of each kind."""
        st = _single("cm250")
        self.assertAlmostEqual(st.lff.total, 0.74, places=8)
        self.assertAlmostEqual(st.hff.total, 0.74, places=8)
        self.assertAlmostEqual(st.lff.total + st.hff.total, 1.48, places=8)

    def test_fragment_kinetic_energy_is_conserved(self):
        for key in self._fissioners():
            riso = _DB["rad"][key]
            with self.subTest(nuclide=riso.name):
                st = _single(key)
                s = fission.split(parse_nuclide(riso.name), riso)
                actual = (sum(e * y for e, y in zip(st.lff.energy, st.lff.yield_))
                          + sum(e * y for e, y in zip(st.hff.energy, st.hff.yield_)))
                self.assertAlmostEqual(actual / (s.fissions * s.tke), 1.0, places=9)

    def test_light_fragment_is_faster_and_lighter(self):
        for key in self._fissioners():
            riso = _DB["rad"][key]
            s = fission.split(parse_nuclide(riso.name), riso)
            with self.subTest(nuclide=riso.name):
                self.assertLess(s.light.mass, s.heavy.mass)
                self.assertGreater(s.light.energy, s.heavy.energy)

    def test_momentum_conservation(self):
        """E_L / E_H == A_H / A_L for every fissioning nuclide."""
        for key in self._fissioners():
            riso = _DB["rad"][key]
            s = fission.split(parse_nuclide(riso.name), riso)
            with self.subTest(nuclide=riso.name):
                self.assertAlmostEqual(s.light.energy / s.heavy.energy,
                                       s.heavy.mass / s.light.mass, places=9)
                self.assertAlmostEqual(s.light.energy + s.heavy.energy, s.tke, places=9)

    def test_zaids_conserve_charge_and_mass(self):
        for key in self._fissioners():
            riso = _DB["rad"][key]
            nuc = parse_nuclide(riso.name)
            s = fission.split(nuc, riso)
            with self.subTest(nuclide=riso.name):
                self.assertEqual(s.light.charge + s.heavy.charge, nuc.z)
                self.assertAlmostEqual(s.light.mass + s.heavy.mass, nuc.a - s.nu, places=6)
                self.assertLessEqual(
                    abs((s.light.zaid % 1000) + (s.heavy.zaid % 1000) - (nuc.a - s.nu)), 1.0)

    def test_tke_matches_published_values(self):
        """2 x the RAD mean energy is the total fission kinetic energy."""
        for key, tke in self.MEASURED_TKE.items():
            s = fission.split(parse_nuclide(_DB["rad"][key].name), _DB["rad"][key])
            with self.subTest(nuclide=key):
                self.assertAlmostEqual(s.tke / tke, 1.0, delta=0.03)

    def test_fragment_energies_match_published_values(self):
        for key, (e_l, e_h) in self.MEASURED_FRAGMENT_E.items():
            s = fission.split(parse_nuclide(_DB["rad"][key].name), _DB["rad"][key])
            with self.subTest(nuclide=key):
                self.assertAlmostEqual(s.light.energy / e_l, 1.0, delta=0.04)
                self.assertAlmostEqual(s.heavy.energy / e_h, 1.0, delta=0.04)

    def test_implied_neutron_multiplicity_is_physical_everywhere(self):
        lo, hi = fission.NU_PHYSICAL_RANGE
        for key in self._fissioners():
            riso = _DB["rad"][key]
            s = fission.split(parse_nuclide(riso.name), riso)
            with self.subTest(nuclide=riso.name):
                self.assertTrue(lo <= s.nu <= hi, f"{riso.name}: nu = {s.nu}")

    def test_cf252_yield_anomaly_is_detected_and_corrected(self):
        """ICRP-107 gives Cf-252's FF yield as fissions, not fragments.

        Read as fragments it implies nu = 7.52; read as fissions, nu = 3.760
        against the published 3.7676.  Without the correction the Cf-252
        fission-fragment source is half the truth.
        """
        riso = _DB["rad"]["cf252"]
        y_ff, _, y_n = fission.fission_yield(riso)
        s = fission.split(parse_nuclide("Cf-252"), riso)
        self.assertTrue(s.note, "the Cf-252 yield anomaly was not flagged")
        self.assertAlmostEqual(s.fissions, y_ff, places=12)      # not y_ff / 2
        self.assertAlmostEqual(s.nu, 3.7676, delta=0.05)
        st = _single("cf252")
        self.assertAlmostEqual(st.lff.total, 3.0920e-02, places=9)
        self.assertAlmostEqual(st.hff.total, 3.0920e-02, places=9)

    def test_only_cf252_needs_the_yield_correction(self):
        flagged = [_DB["rad"][k].name for k in self._fissioners()
                   if fission.split(parse_nuclide(_DB["rad"][k].name), _DB["rad"][k]).note]
        self.assertEqual(flagged, ["Cf-252"])

    def test_non_fissioning_nuclide_has_empty_fragment_classes(self):
        st = _single("co60")
        self.assertEqual(st.lff.total, 0.0)
        self.assertEqual(st.hff.total, 0.0)
        self.assertIsNone(st.lff.zaid)
        self.assertIsNone(st.hff.zaid)

    def test_deck_writes_a_zaid_placeholder_when_there_is_no_fission(self):
        import tempfile

        st = _single("co60")
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "deck.out")
            sdef.write_sdef(path, st, title="unit test")
            text = open(path, encoding="utf-8").read()
        self.assertNotIn("None", text)
        self.assertIn("XX         0.0000000E+00       28", text)

    def test_deck_carries_the_dominant_fissioner_zaids(self):
        import tempfile

        st = _single("cm250")
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "deck.out")
            sdef.write_sdef(path, st, title="unit test")
            text = open(path, encoding="utf-8").read()
        self.assertIn(f"{st.hff.zaid}         7.4000000E-01       28", text)
        self.assertIn(f"{st.lff.zaid}         7.4000000E-01       29", text)

    def test_qa_passes_for_every_fissioning_nuclide(self):
        for key in self._fissioners():
            with self.subTest(nuclide=_DB["rad"][key].name):
                rep = qa.QAReport()
                qa.check_source_term(_single(key), rep)
                fis = [c for c in rep.checks if c.group == "FISSION"]
                self.assertTrue(fis)
                self.assertFalse([c for c in fis if c.status == qa.FAIL],
                                 "\n".join(str(c) for c in fis if c.status != qa.PASS))


class TestBinning(unittest.TestCase):
    def setUp(self):
        self.bins = energybins.read_energy_bins(BINS_PATH)
        self.st = _single("cd109")
        self.binned = energybins.binned_source_term(self.st, self.bins)

    def test_partition_is_contiguous_and_ascending(self):
        rep = qa.QAReport()
        qa.check_partition(self.bins, rep)
        self.assertEqual(rep.status, qa.PASS)

    def test_binning_conserves_yield(self):
        self.assertAlmostEqual(self.binned.photon.total, self.st.photon.total, places=12)
        self.assertAlmostEqual(self.binned.dee.total, self.st.dee.total, places=12)

    def test_binning_reduces_line_count(self):
        self.assertLess(len(self.binned.dee.energy), len(self.st.dee.energy))
        self.assertLessEqual(len(self.binned.dee.energy), len(self.bins))

    def test_binning_leaves_other_classes_alone(self):
        self.assertEqual(self.binned.ben.pdf, self.st.ben.pdf)
        self.assertEqual(self.binned.alpha.yield_, self.st.alpha.yield_)

    def test_each_binned_line_sits_inside_its_bin(self):
        import bisect
        for e, _ in zip(self.binned.dee.energy, self.binned.dee.yield_):
            idx = bisect.bisect_right(self.bins.leb, e) - 1
            self.assertGreaterEqual(e, self.bins.leb[idx])
            self.assertLessEqual(e, self.bins.ueb[idx])


class TestDeckWriting(unittest.TestCase):
    """Write a deck to disk and verify it by reading the text back."""

    def test_round_trip_is_self_consistent(self):
        import tempfile

        st = _single("cu64", points=900)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "deck.out")
            sdef.write_sdef(path, st, title="unit test")
            rep = qa.QAReport()
            qa.check_deck(path, st, rep, label="test")
            self.assertEqual(rep.status, qa.PASS,
                             "\n".join(str(c) for c in rep.checks if c.status != qa.PASS))

    def test_ben_and_bep_appear_as_separate_distributions(self):
        import tempfile

        st = _single("cu64", points=900)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "deck.out")
            sdef.write_sdef(path, st, title="unit test")
            text = open(path, encoding="utf-8").read()
        self.assertIn("3.9000000E-01       24", text)     # BEN yield on SP7
        self.assertIn("1.7408300E-01       25", text)     # BEP yield on SP7
        self.assertIn("c Beta-Electron, negative (BEN) emissions:", text)
        self.assertIn("c Beta-positron (BEP) emissions:", text)


class TestInventoryParsing(unittest.TestCase):
    def test_nuclide_name_spellings_are_equivalent(self):
        for spelling in ("Co60", "Co-60", "CO-60", "co60", "  Co-60  "):
            self.assertEqual(canonical_key(spelling), "co60")
        for spelling in ("Ag110m", "AG-110M", "ag-110m"):
            self.assertEqual(canonical_key(spelling), "ag110m")

    def test_metastable_states_are_distinct(self):
        self.assertNotEqual(canonical_key("Ag-110"), canonical_key("Ag-110m"))
        self.assertIn(canonical_key("Ag-110m"), _DB["rad"])

    def test_mixed_tab_and_space_delimiters(self):
        import tempfile

        text = "Num_Beta_Eng_Points  300\n@Start\nCo-60\t1.0E+10\nCS137   2.0E+09\n" \
               "  Ag-110m \t 3.0E+08 \n@End\n"
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "Input.inp")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(text)
            inv = read_inventory(path)
        self.assertEqual(inv.num_beta_points, 300)
        self.assertEqual(len(inv.items), 3)
        self.assertEqual(inv.keys, {"co60", "cs137", "ag110m"})
        self.assertAlmostEqual(inv.total_activity, 1.23e10)

    def test_duplicate_nuclides_are_summed(self):
        import tempfile

        text = "@Start\nCo-60 1.0E+10\nco60 1.0E+10\n@End\n"
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "Input.inp")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(text)
            inv = read_inventory(path)
        self.assertEqual(len(inv.items), 1)
        self.assertAlmostEqual(inv.items[0].activity, 2.0e10)


class TestLowEnergyRetention(unittest.TestCase):
    """The 2009 code discarded everything below 1 keV; v2 keeps it to 10 eV."""

    def test_sub_kev_electrons_are_retained(self):
        st = _single("in113m")
        low = [e for e in st.dee.energy if e < 1.0e-3]
        self.assertGreater(len(low), 100)
        self.assertGreaterEqual(min(st.dee.energy), sourceterm.INCLUDE_THRESHOLD_MEV)

    def test_nothing_below_the_threshold_survives(self):
        st = _single("cd109")
        for d in (st.photon, st.dee, st.alpha, st.neutron):
            if d.energy:
                self.assertGreaterEqual(min(d.energy), sourceterm.INCLUDE_THRESHOLD_MEV)


class TestFullCaseQA(unittest.TestCase):
    """The shipped example inventory must pass every QA check."""

    def test_shipped_inventory_passes_qa(self):
        inv = read_inventory(os.path.join(ROOT, "Input", "Input.inp"))
        wanted = inv.keys
        rad = {k: v for k, v in _DB["rad"].items() if k in wanted}
        st = sourceterm.build_source_term(inv, rad, _DB["bet"], _DB["ack"], _DB["corr"])
        bins = energybins.read_energy_bins(BINS_PATH)
        binned = energybins.binned_source_term(st, bins)
        rep = qa.QAReport()
        qa.check_source_term(st, rep)
        qa.check_partition(bins, rep)
        qa.check_binning(st, binned, bins, rep)
        self.assertEqual(rep.status, qa.PASS,
                         "\n".join(str(c) for c in rep.checks if c.status != qa.PASS))

    def test_shipped_inventory_contains_dual_emitters(self):
        """If this ever fails the example case has stopped exercising the fix."""
        inv = read_inventory(os.path.join(ROOT, "Input", "Input.inp"))
        self.assertTrue(inv.keys & set(_DB["dual"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
