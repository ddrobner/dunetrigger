#!/bin/bash

idx=100
for file in run_t*.fcl
do
	lar -c $file -n 1 -s test_det.root -o scratch/testboth_$idx.root
	idx=$(($idx+50))
done 
