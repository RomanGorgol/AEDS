#!/bin/bash

rm solver
rm main

mpicc solver.c -o solver
mpicc main_final.c -o main

srun solver
srun -n 10 main output
diff output output.txt
