#!/bin/bash

#input="coo_input_data/matrix_10000x10000_d0.mtx"
#input="coo_input_data/matrix_10000x10000_d1.mtx"

input="coo_input_data/matrix_10000x10000_d10.mtx"
for p in 1 2 3 4 5 6
do
  echo "Running $input with $p process(es)..."
  mpirun -np $p ./power_method $input
done
