#!/bin/bash
#PBS -N gb_1aoi
#PBS -k oe
#PBS -m be


module load amber/18

mpirun -np 56 MMPBSA.py.MPI -O -i mmpbsa.in -cp com.top -rp DNA.top -lp his.top -y 1aoi_MMPBSA.nc
#submit job
#sbatch --partition=multinode --constraint=x2680 --ntasks=56 --ntasks-per-core=1 --time=72:00:00 --exclusive sub_job.sh
