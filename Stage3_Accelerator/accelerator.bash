#!/bin/bash

# AcceLerator Beam Optics (ALBO)

make clean; make accelerator;

inputRootFileDir="/media/andrii/F492773C92770302/SPS_SimulationData/"
outputRootFile="accelerator_25mrad_all.root"

./accelerator $inputRootFileDir $outputRootFile 5

