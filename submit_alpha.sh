#!/bin/bash

for stype in NN RR RT; do
for alpha in `seq 0.1 0.1 1.`; do 
python ./align.py 0.5 1. $alpha $stype
done
done
