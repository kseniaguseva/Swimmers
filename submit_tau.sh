#!/bin/bash

for stype in NN RR RT; do
for tau in `seq 0.1 1. 12.`; do
python ./align.py 0.5 $tau 0.98 $stype
done
done
