#!/bin/bash

for a in 0.5 0.25 0.98; do
for us in `seq 0.1 0.3 4.`; do
python ./align.py $us 1. $a NN
done
done
