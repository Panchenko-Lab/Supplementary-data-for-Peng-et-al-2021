foreach name  {ModelA ModelB ModelC ModelD} {

mol load parm7 ../${name}_rw.prmtop 
mol addfile ../${name}_rw_run1_4000ns.nc step 10 waitfor all
source rename_chainID_${name}.tcl
set sel [atomselect top all frame last]
$sel writepdb ./All_frames_last/${name}_run1.pdb
mol delete all

mol load parm7 ../${name}_rw.prmtop
mol addfile ../${name}_rw_run2_800ns.nc  step 10 waitfor all
source rename_chainID_${name}.tcl
set sel [atomselect top all frame last]
$sel writepdb ./All_frames_last/${name}_run2.pdb
mol delete all

mol load parm7 ../${name}_rw.prmtop
mol addfile ../${name}_rw_run3_800ns.nc  step 10 waitfor all
source rename_chainID_${name}.tcl
set sel [atomselect top all frame last]
$sel writepdb ./All_frames_last/${name}_run3.pdb
mol delete all

mol load parm7 ../${name}_rw.prmtop
mol addfile ../${name}_rw_run4_800ns.nc step 10 waitfor all
source rename_chainID_${name}.tcl
set sel [atomselect top all frame last]
$sel writepdb ./All_frames_last/${name}_run4.pdb
mol delete all

mol load parm7 ../${name}_rw.prmtop
mol addfile ../${name}_rw_run5_800ns.nc step 10 waitfor all
source rename_chainID_${name}.tcl
set sel [atomselect top all frame last]
$sel writepdb ./All_frames_last/${name}_run5.pdb
mol delete all

}
