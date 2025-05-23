#!/bin/bash

p=$1

if [ -z "$p" ]; then
  echo "Usage: $0 <number_of_processes>"
  exit 1
fi

input="coo_input_data/matrix_10000x10000_d1.mtx"

#input="coo_input_data/input_5x5.mtx"

#input="coo_input_data/matrix_100000x100000_d0.mtx"
#input="coo_input_data/input_6x6.mtx"
#input="coo_input_data/ue_input_6x6.mtx"
#input="coo_input_data/matrix_100x100_d10.mtx"
#input="coo_input_data/matrix_10000x10000_d0.mtx"

echo "________ Parallel run p=$p - $input ________"
echo ""
mpirun -np $p ./power_method $input
