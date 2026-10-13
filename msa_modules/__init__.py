"""MCNP_SOURCE_APP -- modular Python port of the 2009 FORTRAN ``mcnpsource`` utility.

This package converts an ICRP-107 based radionuclide activity inventory into an
MCNP6.3 SDEF (source definition) ionising-radiation emission source term.

Modules
-------
elements     Periodic table and radionuclide-name parsing / normalisation.
databases    Parsers for the ICRP-07.RAD / .BET / .ACK / -CORR.BEP data files.
inventory    Reader for the user ``Input.inp`` activity inventory.
betaspectra  Separation of the beta-minus and beta-plus emission spectra.
fission      Light/heavy partition of the spontaneous-fission fragment source.
sourceterm   Construction of the activity-weighted emission source term (IREST).
energybins   Collapse of the DEP/DEE classes onto a fixed energy partition.
sdef         Writer that emits the MCNP6.3 SDEF construct.
qa           Self-verification of the source term and of the written decks.
plotting     Data-visualisation graphs of the source term.
"""

__all__ = [
    "elements",
    "databases",
    "inventory",
    "betaspectra",
    "fission",
    "sourceterm",
    "energybins",
    "sdef",
    "qa",
    "plotting",
]

__version__ = "2.1.0"
