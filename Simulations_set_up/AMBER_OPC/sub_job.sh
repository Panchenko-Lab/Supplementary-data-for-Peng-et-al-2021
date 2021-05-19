#!/bin/bash
#PBS -N 1aoi1
#PBS -k oe
#PBS -m be




module load amber/18

mpirun  $AMBERHOME/bin/pmemd.MPI -O -i Min.in -o output/Min.out -p 1kx5_amber.prmtop -c 1kx5_amber.inpcrd -r output/Min.ncrst -inf output/Min.mdinfo 
mpirun  $AMBERHOME/bin/pmemd.MPI -O -i Equil_v.in -o output/Equil_v.out -p 1kx5_amber.prmtop -c output/Min.ncrst -r output/Equil_v.ncrst -x output/Equil_v.nc -inf output/Equil_v.mdinfo
mpirun  $AMBERHOME/bin/pmemd.MPI -O -i Equil_pt.in -o output/Equil_pt.out -p 1kx5_amber.prmtop -c output/Equil_v.ncrst -r output/Equil_pt.ncrst -x output/Equil_pt.nc -inf output/Equil_pt.mdinfo
#mpirun  $AMBERHOME/bin/pmemd.MPI -O -i Prod.in -o output/Prod.out -p 1kx5_amber.prmtop -c output/Equil_pt.ncrst -r output/Prod.ncrst -x output/Prod.nc -inf output/Prod.info 

#submit job
#sbatch --partition=multinode --constraint=x2680 --ntasks=280 --ntasks-per-core=1 --time=168:00:00 --exclusive sub_job.sh
