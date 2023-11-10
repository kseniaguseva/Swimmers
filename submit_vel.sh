#!/bin/bash

for stype in NN RR RT; do
for us in `seq 0.1 0.2 5.`; do 
python ./align.py $us 1. 0.98 $stype
done
done
