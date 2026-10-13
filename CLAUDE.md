## The code “MCNP_SOURCE_APP”

Computational radiation protection analysts are typically faced with a structure that contains a known radionuclide activity inventory, specified as an array $(i;A_{i})$ where $i$ is the identity of a radionuclide and $A_i$ is its activity in Bq.  The collection of radionuclides in the $(i;A_{i})$ array will each emit characteristic ionising radiations.  This gives rise to a ionising radiation emission source term (IREST), $(E_g;Y_g)_j$ where $g$ is the energy index and $j \in \{n,\gamma,e^-,\beta^-,\beta^+,\alpha^{2+},REC,FF\}$.  Here REC refers to the recoiling atom that had emitted an alpha particle, while FF refers to the fission fragment.  There are two fission fragment peaks - the FFL (Fission Fragment, Light) and FFH (Fission Fragment, Heavy) peaks.

The purpose of  MCNP_SOURCE_APP is:

* Read the $(i;A_{i})$ array in the file “Input.inp” in the subfolder “Input_Output”.  Check if all requested radionuclides are indeed present in the ICRP-107 database file “ICRP-07.RAD”.  If not, terminate and supply a list of unsupported radionuclides to the user.
* Read $(E_g;Y_g)_j$ where $g$ is the energy index and $j \in \{n,\gamma,e^-,\beta^-,\beta^+,\alpha^{2+},REC,FF\}$ from the ICRP-107 database files “ICRP-07.RAD” and “ICRP-07.BET” and “ICRP-07.ACK”.
* Produce an an SDEF (source definition) construct in MCNP6.3 format.  This file must have the name “SDEF_MCNP6.out” and must be written in the subfolder “Input_Output”.  An example output is “Output_sample.out” in the subfolder “Input_Output”.

FACT: Good quality radiation emission data with correct beta-plus (positron, i.e. BEP) emission spectra and Yield_BEP and Yield_BEN exist in the shared folder “C:\Py_Data_Shared” on this PC.  You may access and use this data.

## Task 1

Write SDEF radiation emission data into the subfolder “Output_Spectroscopic_ALL” for all nuclides in the ICRP-107 RAD/ACK/BET databases.  Do this on-demand, i.e. do not repeat this each time the code executes.  Add a clear instruction in MCNP_SOURCE_APP.py to enable the on-demand generation of `ALL` radionuclide output.

Note that existing output for nuclides emitting both BEN and BEP, are all wrong.  Yield_BEP is unique.  Yield_BEN is unique.  Energy spectra for BEN and BEP are not identical.. Correct please.

## Task 2

There is an output folder “Output_Spectroscopic_Case” for spectroscopically accurate output for the $(i;A_{i})$ array being analysed.  The code as it is now, writes spectroscopically accurate output.  This can amount to over 20,000 lines of output.  Maintain this ability.  Always write this spectrometrically accurate reference output.

## Task 3

There is a subfolder “data_EnergyBins” which contains the file “EnergyBins_PE.dat”.  This is the energy partition for the binned discrete energy photon and discrete energy electron source term.

Problem: It can easily happen that the spectroscopically accurate SDEF contain, say, 3E4 photon emission lines and 6E4 electron emission lines.  This is excessive and Monte Carlo codes may struggle to sample all emission lines equitably.  The need therefor arises to bin emissions into a fixed energy partition, i.e. energy bins.  The prescribed energy partition is supplied.  The code must write the SDEF with binned DEP and DEE output data in “SDEF_MCNP6_binned.out” in the subfolder “Output_Binned”.

## Task  4

Graphs of the ionising radiation emission source term must be placed in the subfolder GraphBank.  Read the new SKILLS.md in folder GraphBank and refine the graph writing module as needed.  The Graphs must be of BINNED photon emission lines and BINNED discrete electron emission lines.  The Graphs must NOT be of the spectroscopically precise emission lines and discrete electron emission lines, because there can be tens of thousands of these lines, and they overlap and blur excessively.  All other emission lines must still be graphed - Alphas, neutrons, Recoils, beta- and beta+.  Also write the code to put the fission fragment emissions on the graph, even though they will not be visible - they are emitted far above the 15 MeV max energy.

