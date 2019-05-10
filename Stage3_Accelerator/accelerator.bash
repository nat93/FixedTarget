#!/bin/bash

# Accelerator Beam Optics (AlBO)

make clean; make accelerator;

inputRootFileDir="/media/andrii/F492773C92770302/SPS_SimulationData/"
outputRootFile="accelerator_pion_10mrad.root"

./accelerator $inputRootFileDir $outputRootFile 8 0.010

