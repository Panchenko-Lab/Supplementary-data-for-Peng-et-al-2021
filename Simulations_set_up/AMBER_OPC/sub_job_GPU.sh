#!/bin/bash
#PBS -N 1aoi1
#PBS -k oe
#PBS -m be

module load amber/18-gpu



#$AMBERHOME/bin/pmemd.cuda -O -i Min.in -o output/Min.out -p 1kx5_amber.prmtop -c 1kx5_amber.inpcrd -r output/Min.ncrst -inf output/Min.mdinfo 
#$AMBERHOME/bin/pmemd.cuda -O -i Equil_v.in -o output/Equil_v.out -p 1kx5_amber.prmtop -c output/Min.ncrst -r output/Equil_v.ncrst -x output/Equil_v.nc -inf output/Equil_v.mdinfo
#$AMBERHOME/bin/pmemd.cuda -O -i Equil_pt.in -o output/Equil_pt.out -p 1kx5_amber.prmtop -c output/Equil_v.ncrst -r output/Equil_pt.ncrst -x output/Equil_pt.nc -inf output/Equil_pt.mdinfo
#mpirun -np 2 $AMBERHOME/bin/pmemd.cuda.MPI -O -i Prod.in -o output/Prod5.out -p 1kx5_amber.prmtop -c output/Prod4.ncrst -r output/Prod5.ncrst -x output/Prod5.nc -inf output/Prod5.info 
$AMBERHOME/bin/pmemd.cuda -O -i Prod.in -o output/Prod53.out -p 1kx5_amber.prmtop -c output/Prod52.ncrst -r output/Prod53.ncrst -x output/Prod53.nc -inf output/Prod53.info

#submit job
#sbatch --partition=gpu --gres=gpu:p100:2  --time=10:00:00  sub_job_GPU.sh
