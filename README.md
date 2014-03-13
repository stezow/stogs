This is the README of the SToGS Project
=====

SToGS stands for
    Simultation Toolkit fOr Gamma-ray Spectroscopy


The goal is to develop a GEANT4 package in order to fully simulate 
experiment in Nuclear Structure.


Organization: 

csrc: directory containing source files (.hh and .cc)


To configure / compile
cmake -DGeant4_DIR=/where/geant/has/been/built ../

Stucture 
Anything starting by SToGS_G4 is pure geant4 code and the .hh and .cc files could be moved to other G4 packages.
Anything starting by SToGS is complex SToGS code. The .hh and .cc could be moved to other G4 packages.
