#!/bin/bash

make clean; make accelerator;

inputRootFileDir="/media/andrii/F492773C92770302/SPS_SimulationData/"
outputRootFile="accelerator.root"


./accelerator $inputRootFileDir $outputRootFile 2

