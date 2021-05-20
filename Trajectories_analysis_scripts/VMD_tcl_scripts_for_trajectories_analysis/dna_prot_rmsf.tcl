foreach name  {1aoi 1eqz ext2 sym} {

mol delete all
mol load parm7 ../${name}_rw.prmtop
mol addfile ../${name}_rw_run1_4000ns.nc first 2500 last 190000 step 10 waitfor all
mol addfile ../${name}_rw_run2_500ns.nc first 2500 last 24000 step 10 waitfor all
mol addfile ../${name}_rw_run3_500ns.nc first 2500 last 24000 step 10 waitfor all
mol addfile ../${name}_rw_run4_500ns.nc first 2500 last 24000 step 10 waitfor all
mol addfile ../${name}_rw_run5_500ns.nc first 2500 last 24000 step 10 waitfor all


source rename_chainID_${name}.tcl


set nframes [expr  [molinfo top get numframes] - 1 ]

#Align the alpha helix core ti the start structre
set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"


set frame0_alpha_core_ca [atomselect top $alpha_core_ca frame 0]
set sel_alpha_core_ca [atomselect top $alpha_core_ca]

set all [atomselect top "all"]

for { set i 1 } { $i<$nframes } { incr i } {

$sel_alpha_core_ca frame $i
$all frame $i
set trans_mat [measure fit $sel_alpha_core_ca $frame0_alpha_core_ca]
$all move $trans_mat
}

set rmsf [measure rmsf $all first 0 last $nframes]
# set rmsf [measure rmsf $all first 9000 last 10000]

$all frame $nframes
# was 10001
$all set beta $rmsf
$all writepdb ./rmsf/${name}_nucl_rmsf.pdb

#let's do some file output
set chains [list CHA CHE CHB CHF CHC CHG CHD CHH CHI CHJ]
set chains_names [list H3_1 H3_2 H4_1 H4_2 H2A_1 H2A_2 H2B_1 H2B_2 DNA_I DNA_J]


set outfile [open ./rmsf/${name}_rmsf_chains.dat w]
puts $outfile "RMSF of CA or P atoms of nucleosome"
puts $outfile "Generated by VMD rmsf"
puts $outfile "Residue ID"
puts $outfile "RMSF, A"


set ts "Resid"
foreach segname  $chains_names {
append ts "\t$segname"
}
puts $outfile $ts

for  { set i -93 } { $i<=150 } { incr i } {
set ds "$i"
foreach seg  $chains {

set sel [atomselect top "segname $seg and name CA P and resid '$i' "]
if {[$sel num] == 0} { set beta 0
} else {
set  beta [$sel get beta]
}

append ds "\t$beta"

}
puts $outfile $ds
}

close $outfile 


set outfile [open ./rmsf/${name}_rmsf_chains_sch.dat w]
puts $outfile "Average RMSF of side chains and bases atoms of nucleosome"
puts $outfile "Generated by VMD rmsf"
puts $outfile "Residue ID"
puts $outfile "RMSF, A"


set ts "Resid"
foreach segname  $chains_names {
append ts "\t$segname"
}
puts $outfile $ts

for  { set i -93 } { $i<=150 } { incr i } {
set ds "$i"
foreach seg  $chains {
# set sel [atomselect top "segname $seg and (not backbone) and resid '$i' and noh"]
set sel [atomselect top "segname $seg and (not backbone) and (not name P O1P O2P O3' O5' C5' C1' C2' C3' C4' O4') and resid '$i' and noh"]
if {[$sel num] == 0} { set mean 0
} else {
set  beta [$sel get beta]
set bl [llength $beta]
set mean 0
for {set j 0} { $j<$bl } { incr j } {
set mean [expr $mean + [lindex $beta $j]]
}
set mean [expr $mean / $bl]
}


append ds "\t$mean"

}
puts $outfile $ds
}

close $outfile 

}
exit