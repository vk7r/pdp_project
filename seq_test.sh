#!/bin/bash

# Run the program sequentially with 1 MPI process

#input="coo_input_data/matrix_100000x100000_d0.mtx"
#input="coo_input_data/input_6x6.mtx"
#input="coo_input_data/matrix_100x100_d10.mtx"
#input="coo_input_data/matrix_10000x10000_d0.mtx"

input="coo_input_data/input_5x5.mtx"

mpirun -np 1 ./power_method $input
