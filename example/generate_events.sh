#!/bin/bash

output_data_dir=$1
n_events=1

#cosmics generation
lar -c gen_events_cosmiconly.fcl -n $n_events -o $output_data_dir/test_cosmicsonly.root

#g4 stage
lar -c standard_g4_protodunehd.fcl -s $output_data_dir/test_cosmicsonly.root -o $output_data_dir/test_g4_cosmicsonly.root

#detsim stage
lar -c standard_detsim_protodunehd.fcl -s $output_data_dir/test_g4_cosmicsonly.root -o $output_data_dir/test_detsim_cosmicsonly.root

