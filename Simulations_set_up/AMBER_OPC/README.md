# Files for running simulations using AMBER forcefield and OPC water model

## [Nucleosome_model_PDBs](Nucleosome_model_PDBs)

Different nucleosome models used for simulations

## [PTM_parameters](PTM_parameters)

### Forcefield Parameters for post-translational modifications.

These parameters are originally taken from: Khoury, G. A., Thompson, J. P., Smadbeck, J., Kieslich, C. A. & Floudas, C. A. Forcefield_PTM: Ab Initio Charge and AMBER Forcefield Parameters for Frequently Occurring Post-Translational Modifications. J Chem Theory Comput 9, 5653-5674, doi:10.1021/ct400556v (2013).

## gen_nucl.leap

Script to prepare the residue library and force field parameters for use with LEaP progarm.

To execuate the script:

tleap -f gen_nucl

## Min.in, Equil_v.in, Equil_pt.in, Prod.in

Simulation configration files for energy minimizations, equilibration, and production run.

## sub_job.sh, sub_job_GPU.sh

Sample scripts to run the simulations on CPU or GPU.

## Required Progarms
Amber (https://ambermd.org)

Simulations were performed with the version of Amber18.

## MD trajectories

MD Trajectories are archived at: 
https://zenodo.org/record/4780649#.YKu7Ny9h0mI

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4780649.svg)](https://doi.org/10.5281/zenodo.4780649)
