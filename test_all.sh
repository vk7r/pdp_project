#!/bin/bash

for input in new_coo_input_data/*.mtx
do
  for p in 1 2 4 6
  do
    echo "Running $input with $p process(es)..."
    mpirun -np $p ./power_method $input
  done
done
