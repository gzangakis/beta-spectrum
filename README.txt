# Beta Decay Energy Spectrum for Isotopes and Reactors

A research project to model the neutrino energy spectrum for isotopes, single isotope cumulative fission, and mixed-fuel nuclear reactors as well provide a starting point to analyze corresponding experimental datasets. Supervisor: Professor Jim Napolitano, Temple University.

This project takes as input, ENDF and ENSDF formatted nuclear data files and outputs the corresponding beta or neutrino energy spectrum for a chosen physical system. 

Additionally a file, reactor_neutrino_spectrum.py is included which models the energy spectrum of neutrinos as produced by a fission reactor of pre-determined time-averaged isotope fuel fractions as an ab-initio beta decay calculation. 

## Getting Started

All files are files are written in Python and only require input ENDF or ENSDF files found at the National Nuclear Data Center. If just modeling the beta decay of an isotope, only ENSDF files are necessary. However, if modeling the beta decay of a fissile isotope or neutron induced fission reactor, fission yield datasets are necessary as well as the ENSDF beta decay files. These are found in ENDF files. 

ENDF: https://www.nndc.bnl.gov/endf/

	-retrieval of datasets using SIGMA is suggested

ENSDF: https://www.nndc.bnl.gov/ensdf/

## INCLUDED

Single_Isotope_Beta_Spectrum.py - Calculates the electron kinetic energy spectrum of an isotope
Single_Isotope_Neutrino_Spectrum.py - Calculates the neutrino energy spectrum of an isotope

Output_Data - a directory in which you can store the output files. Contains two examples.

ENSDF_beta.zip - a directory of ENSDF files, for your use. ENSDF files are named after the daughter nucleus, therefore in order to compute the beta spectrum of, for example, 16O, one should use the ENSDF file 16N. 



### Prerequisites

Python 2.x


## Author: Gabe Zangakis
