#!/bin/bash

# Run the program sequentially with 1 MPI process

input="coo_input_data/input_1000x1000_d0.0005.mtx"
#input="coo_input_data/input_6x6.mtx"

mpirun -np 1 ./main $input
