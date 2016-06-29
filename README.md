# flame
Flame data reduction pipeline
-----------------------------

by Sirio Belli (MPE) 2016


To install: download the flame/ directory and add it to your IDL path.
The latest versions of the NASA IDL library, Coyote library, and mpfit library are required.

To run:
1) copy flame_driver.pro to your working directory
2) create the science.txt and the darks.txt files, containing the file names of the relevant science and dark frames
3) edit the top part of flame_driver.pro
4) run flame_driver.pro either line-by-line or all at once by typing "idl flame_driver.pro" at the command line
