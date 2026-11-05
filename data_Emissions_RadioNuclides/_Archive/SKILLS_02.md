SFN = Spontaneous fission neutron

DEP = Discrete energy photons

DEE = Discrete energy electrons

DEA = Discrete energy alpha

DER = Discrete energy recoils from alpha emissions.

DEFFL = Discrete energy emissions, fission fragment, Light

DEFFH = Discrete energy emissions, fission fragment, High 

BEN = Beta-electron, negative

BEP = Beta positron

The ionising radiation emission data in this folder is that of ICRP-107 (2008), containing augmented SFN and PFG emission data.  PFG = Prompt Fission Gamma.

Problem in ICRP-107: When a nuclide emits both beta-minus electrons and beta-plus electrons, ICRP-107 wrongly allocates the BEN energy spectrum to the BEP positron emission spectrum.  This is not accurate.  This must be corrected when computing emission spectra.

Historical error made by CLAUDE code when generating the output files in “Output_Spectroscopic_ALL” = BEP emission yields are wrong when a radionuclide emits both BEP and BEN.  Yield_BEN was wrongly allocated to Yield_BEP in these cases.  This must be corrected.

Historical error when in output files in “Output_Spectroscopic_ALL” = When a radionuclide emits both BEP and BEN, the BEN energy spectrum was simplistically allocated to BEP energy spectra.  This stems from a systematic error in the ICRP-107 database.  This must be corrected.  Output in subfolder “Output_Spectroscopic_ALL” is wrong for all radionuclides that emit both BEP and BEN.
