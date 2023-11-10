#!/bin/bash

for stype in NN RR RT; do
for tau in `seq 0.1 0.5 20.`; do 
python ./align.py 0.5 $tau 0.98 $stype
done
done
