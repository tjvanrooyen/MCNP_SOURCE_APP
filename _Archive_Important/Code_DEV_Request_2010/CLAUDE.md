Hello CLAUDE

## Request to develop “MCNP_SOURCE_APP”

In this folder and subfolders you will find a legacy FORTRAN code developed in 2009, along with its nuclear data.  Its functionality must be ported to a modular collection of Python Scripts.  The user-facing Python script must have the name “MCNP_SOURCE_APP.py”.  The entire collection of Python scripts will be collectively known as “MCNP_SOURCE_APP”.

MCNP_SOURCE_APP is a utility code that reads the ionising radiation emissions by radionuclide from the files “ICRP-07.BET” and “ICRP-07.RAD”.  There is also a low-energy electron emission file “ICRP-07.ACK”.  All three these files contain text data.  These files are located in the subfolder “Emissions_RN”.  A SKILLS.md file in the latter folder explains the content of the three data files.

The 2009 App Development request is in the file “Request for the development of a code MCNP_SOURCE_APP.pdf” in the subfolder “Code_DEV_Request”.

## The function and logic of the code MCNP_SOURCE_APP

Computational radiation protection analysts are typically faced with a structure that contains a known radionuclide activity inventory, specified as an array $(i;A_{i})$ where $i$ is the identity of a radionuclide and $A_i$ is its activity in Bq.  The collection of radionuclides in the $(i;A_{i})$ array will each emit characteristic ionising radiations.  This gives rise to a ionising radiation emission source term (IREST), $(E_g;Y_g)_j$ where $g$ is the energy index and $j \in \{n,\gamma,e^-,\beta^-,\beta^+,\alpha^{2+},REC,FF\}$.  Here REC refers to the recoiling atom that had emitted an alpha particle, while FF refers to the fission fragment.  There are two fission fragment peaks - the FFL (Fission Fragment, Light) and FFH (Fission Fragment, Heavy) peaks.

The purpose of  MCNP_SOURCE_APP is:

* Read the $(i;A_{i})$ array in the file “Input.inp” in the subfolder “Input_Output”.  Check if all requested radionuclides are indeed present in the ICRP-107 database file “ICRP-07.RAD”.  If not, terminate and supply a list of unsupported radionuclides to the user.
* Read $(E_g;Y_g)_j$ where $g$ is the energy index and $j \in \{n,\gamma,e^-,\beta^-,\beta^+,\alpha^{2+},REC,FF\}$ from the ICRP-107 database files “ICRP-07.RAD” and “ICRP-07.BET” and “ICRP-07.ACK”.
* Produce an an SDEF (source definition) construct in MCNP6.3 format.  This file must have the name “SDEF_MCNP6.out” and must be written in the subfolder “Input_Output”.  An example output is “Output_sample.out” in the subfolder “Input_Output”.

Graphs of the ionising radiation emission source term must be placed in the subfolder GraphBank.



