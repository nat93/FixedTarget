#!/bin/bash

# AcceLerator Beam Optics (ALBO)

make clean; make accelerator;

inputRootFileDir="/media/andrii/F492773C92770302/SPS_SimulationData/"
outputRootFile="accelerator_10mrad_all.root"

./accelerator $inputRootFileDir $outputRootFile 6 0.010

