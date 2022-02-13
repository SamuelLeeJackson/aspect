# aspect
This repository contains a software code to estimate the aspect related uncertainty in individual near-Earth asteroid photometric phase curves.

This code is distributed as a companion to the paper "The Effect of Aspect Changes on Near-Earth Asteroid Phase Curves", (link to be added). 

The releases of this code are archived on zenodo.
For the latest release please see:   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6062631.svg)](https://doi.org/10.5281/zenodo.6062631)

# Compilation and running

To compile the code, we provide a Makefile, so for unix systems with make and g++ installed, just clone the repository and run "make".
This code is provided as-is with no testing across different platforms.
Developed and tested on a CentOS 7 system.

To run the code, you must first have an input geometry file.
This file must contain the geometry information of the phase curve for which you wish to estimate the aspect uncertainty.
The file is tab delimited, with the following columns in the following order:

* Heliocentric distance [AU]
* Heliocentric longitude [degrees]
* Heliocentric latitudes [degrees]
* Geocentric distance [AU]
* Geocentric longitude [degrees]
* Geocentric latitude [degrees]
* Photometric uncertainty [mag]

The photometric uncertainty is required for each position on the phase curve so that the HG fitting for each model is constrained in a similar way to the real phase curve.
We provide example input geometry files for testing.
They must follow the naming convention "{asteroid_num}-geom.txt".

To run the code, two arguments are required.
First argument is the asteroid number (for collecting the geometry and asteroid inclination automatically).
The second is the observed light curve amplitude, to ignore models that do not match the minimum elongation as the amplitude implies.
If the light curve amplitude is unknown, provide a value of 0.

The code is run as follows for (159402) 1999 AP10.

    ./aspect 159402 0.4
