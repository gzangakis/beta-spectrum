# Beta Decay Energy Spectra for Isotopes and Reactors

A package to model the neutrino energy spectrum for isotopes, single isotope cumulative fission, and mixed-fuel nuclear reactors as well provide a starting point to analyze corresponding experimental datasets. 

This project's code takes as input ENDF and ENSDF formatted nuclear data files and outputs the corresponding beta or neutrino energy spectrum for a chosen physical system. 

## Getting Started

All files are files are written in Python and only require input ENDF or ENSDF files found at the National Nuclear Data Center. If just modeling the beta decay of an isotope, only ENSDF formatted beta decay files are necessary. However, if modeling the beta decay of a fissile isotope or neutron induced fission reactor, fission yield datasets are necessary as well as the ENSDF beta decay files. These are found in ENDF files. 

Familiarizing oneself with ENSDF and ENDF formatted files is suggested as they contain the isotopic nuclear structure and decay informaion. This package can be used as a reference to create your own code for modeling other physical processes or decays.

ENDF: https://www.nndc.bnl.gov/endf/

	-retrieval of datasets using the SIGMA tool is suggested

ENSDF: https://www.nndc.bnl.gov/ensdf/

## INCLUDED

BetaDecay.py - Calculates the electron kinetic energy spectrum or electron-antineutrino energy spectrum of an isotope as specfied by the user.

Output_Data - a directory in which you can store the output files. Contains two examples.

ENSDF_beta.zip - a directory of ENSDF files, for your use. An updated list can be found at the NNDC (link above).

	**ENSDF files are named after the daughter nucleus, therefore to compute spectra one should use the filename corresponding to the daughter nucleus. For example, to calculate a beta spectrum for 16O, one should use the ENSDF file 16N.**

Spontaneous_Fission_Beta_Spectrum.py - Using both ENSDF and ENDF files, this program computes a model beta spectrum for an isotope which undergoes spontaneous fission. The ENDF formatted fission yields for 252Cf is included.

16O - an example ENSDF file for input into BetaDecay.py. It contains the beta decay information of 16N.

### Prerequisites

Python 2.x


## Author: Gabe Zangakis
