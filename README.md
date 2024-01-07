# openEMS_on_sky130

This repository contains some basic tools and testbenches to run facilitate EM simulation of passive devices designed on the open-source Skywater 130nm CMOS process.

It contains a conversion tool to translate a GDSII layout to geometry in matlab/octave syntax needed to run a simulation, a microstrip testbench and a CPW testbench., both using the matlab/octave interface.

## Installation

Just clone this repository to your working area. It is expected that you already have all necessary tools up and running (openEMS, python, octave). Further instructions can be found on the openEMS page (https://docs.openems.de/install/index.html). Information on the sky130 process can be found at https://skywater-pdk.readthedocs.io/en/main/index.html

## Instructions

To run a simulation just launch octave and run the desired script. It is advisable to read the openEMS instructions on how to set up the environment.
To run the GDS translation script just type

>python gds2matlab.py <path_to_gds>

It is very advisable to read the tutorial _Using OpenEMS with IHP SG13G2 v1.pdf_ found at https://github.com/IHP-GmbH/IHP-Open-PDK/tree/dev/ihp-sg13g2/libs.tech/openems/testcase/SG13_Octagon_L2n0/OpenEMS_Matlab. It contains invaluable information on how to interpret the simulation script and on how to set up a good quality mesh to obtain accurate results.

## Acknowledgments

I would like to thank Volker Muelhaus on teaching me how to set up and understand various openEMS commands, as well as for writting the simulation scripts I based mine on. Also, I would like to thank IHP for putting the effort to facilitate the access for open-source EM simulation on integrated devices.
