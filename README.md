single_particle
===============

Single charged particle in the given electro-magnetic field

At the moment only works for single electron in a plane Gaussian wave travelling along the z direction

Requirements:

C-compiler:
per default gcc is set as default - you can change that in src/Makefile by changing CC appropriately

To run the code:

go to build directory
make

go to run directory
./particle.e

This generates binaries *.dat that can be read (i.e. in fread or similar commands)
