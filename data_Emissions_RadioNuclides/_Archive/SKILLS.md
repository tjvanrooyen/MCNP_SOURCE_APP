In nuclear medicine radiotherapy, specifically Targeted Radionuclide Therapy (TRT), Auger & Coster-Kroning (ACK) electron emissions are important.  These are low-energy electrons with a high LET and good therapeutic potential if they are carried into the cell nucleus and are emitted close to the DNA of targeted diseased cells.

The ICRP-107 (2008) database file “ICRP-07.ACK” contains ACK low-energy electron emissions for specific radionuclides.  Consider the first few lines of this file.  “In-113m” is the name of the radionuclide.  Its listed half-life is 1.6579 hours (m = minutes, h = hours, d = days, y = years).  In-113m has 682 low-energy electron emissions.  The emission at E = 1.99000E+00 eV has an emission yield of Y = 1.20318E-04.  MCNP_SOURCE_APP must convert the energy unit eV to MeV so that transport codes like MCNP6.3 can make use of the radiation emission source term.  MCNP6 wants to read the radiation emission data as (E,Y) data, i.e. the (Y,E) data in ICRP-107 must be converted to (E,Y) data.

`In-113m      1.6579h         682
1.20318E-04 1.99000E+00 N2 N3 O2
1.73503E-04 2.32000E+00 N2 N3 O3
1.09904E-04 3.13000E+00 M4 M5 O2
3.49403E-04 3.46000E+00 M4 M5 O3
4.39225E-01 3.83000E+00 N4 O1 O1
3.19620E-01 8.84000E+00 N5 O1 O2
1.09408E+00 9.17000E+00 N5 O1 O3`

Information in the last 3 columns is not relevant from a radiation transport perspective.  Just ignore the “N2 N3 O2” and similar stuff.

Important:

* In “ICRP-07.BET” the beta electron emission is listed as $(E,Y_{pdf})$ with the unit of the energy $E$ in MeV.  Here Y_pdf refers to the probability density function of $\beta^-$ (BEN = beta negatron) and $\beta^+$ (BEP = beta positron) emissions, listed at selected energy points.

* In “ICRP-07.RAD” the ionising radiation emissions of the radionuclides are listed in the following format:

  `  n Y E (in MeV)  Z`

  where n is a redundant number corresponding to the radiation particle type Z.

  Z can be:

  `X     = X ray` 
  `G      = Gamma ray = discrete energy photon (DEP)`
  `PG     = Prompt Gamma Ray = discrete energy photon (DEP)`
  `DG     = Delayed Gamma Ray = discrete energy photon (DEP)` 
  `AQ     = Annihilation Quantum -- an annihilation photon from antiparticle annihilation. Ignore and generate a positron, which will be transported in MCNP to eventually produce annihilation photons.` 
  `B-     = Negative Beta particle = BEN = continuous energy negative electron with energy spectrum in the *.BET file.`
  `B+     = Positive Beta particle = Beta positron = BEP = continuous energy positron with emission energy spectrum in the *.BET file.`
  `AE     = Auger electron = Discete Energy Electron (DEE)`
  `IE     = Internal conversion electron  = Discete Energy Electron (DEE)`
  `A      = Alpha particle = Discete Energy Alpha-particle (DEA)`
  `AR     = Recoil atom from alpha particle emission = Discete Energy Recoil-particle (DER)`
  `N = Neutron`
  `FF = Fission Fragment`

  | Particle designator number N | Particle type Z                                           |
  | ---------------------------- | --------------------------------------------------------- |
  | 1                            | PHOTON = DEP                                              |
  | 2                            | PHOTON = DEP                                              |
  | 3                            | PHOTON = DEP                                              |
  | 4                            | Beta Positron = BEP                                       |
  | 5                            | Beta Electron = BEN                                       |
  | 6                            | Discrete energy electron = DEE                            |
  | 7                            | Auger electron = DEE                                      |
  | 8                            | Alpha particle = DEA                                      |
  | 9                            | Alpha recoil atom = Discrete Energy Recoil = DER          |
  | 10                           | FF = Fission Fragment = DEFF.  Partition into LFF and HFF |
  | 11                           | Neutron (SFN).                                            |
  
  
  
  ## Low energy ACK electron emission
  
  The Fortran version of MCNP_SOURCE_APP failed to include electron and photon emission below 1 keV = 0.001 MeV into the MCNP6 format SDEF radiation emission source term.  This must be corrected in the new version of the code.  Include all electron and all photon emissions down to 1.0E-5 MeV = 10 eV into the SDEF.  The MCNP6 user will cut electron transport at circa 15 eV with a CUT:e card.
  
  