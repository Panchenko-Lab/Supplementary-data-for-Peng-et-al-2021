# parameter files for running simulations using CHARMM forcefield and TIP3P water model.

## [Nucleosome_model_PDBs](Nucleosome_model_PDBs)
Different nucleosome models used for simulations


##  [toppr](toppr), [par](par)
Topology and parameter files used for perpraing MD files.

CHARMM36m FF used for protein and CHARMM36 FF ued for DNA 

## psfgen_nucl.tcl 
script to prepare the MD files for simulations using psfgen in VMD.
To execuate the script:
vmd -e psfgen_nucl.tcl 

## [example_output](example_output)
Some examples of the generated MD files. 

## min_equil.conf, md.conf
Simulation configration files for energy minimizations, equilibration, and production run.

## sub_job.sh
Sample script to run the simulations on CPU.

## Required Programs (used version)
NAMD 2.12
VMD 1.9.3
