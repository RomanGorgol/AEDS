#!/bin/bash

rm solver
rm main

mpicc solver.c -o solver
mpicc main_final2.c -o main

srun --mem=30G solver
srun --mem-per-cpu=10000 -n 2 main output 100 1000

#for pr in 2 4 8 16 32 64
#do
#	echo -ne "$pr\t"
#	srun --mem=30G -n $pr  main output
#done
	#diff output output.txt
