#!/bin/bash
#PBS -N sample
#PBS -k oe
#PBS -m be


module load amber/18

mpirun -np 56 MMPBSA.py.MPI -O -i mmpbsa.in -cp com.top -rp DNA.top -lp his.top -y sample_MMPBSA.nc
