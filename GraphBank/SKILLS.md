The code suite MCNP_SOURCE_APP must write data visualisation graphs to this folder.  Ionising radiation emissions must be visualised.  Energy axis must stretch from 1E-5 MeV to 15 MeV, log-scale.

The minimum of the vertical axis (yield axis) must be no smaller than 1E-12.

* Colour code for neutron emissions: Green.
* Colour code for DEP photon emission lollipops: Blue
* Colour code for DEE discrete energy electron lollipops: Violet
* Colour code for BEN negative beta electron continuous line plot: Light Royal Blue
* Colour code for BEP beta positron continuous line plot: Rouge.
* Colour code for DEA = alpha particle emission lollipops: Brown
* Colour for atomic Recoils resulting from alpha-particle emission: Black.
* Colour for Fission fragment emissions: royal purple.  These emissions are at such high energies that they will not normally display.  However, give the user the option to easily increase Eng_max of the energy axis so that even these high energy emissions can be visualised on demand.
