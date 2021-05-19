#we need to make a heat plot of RMSD evolution of certain bases
#
# 
#
#

foreach {name run} {1aoi run1 1eqz run1 ext2 run1 sym run1} {

mol load parm7 ../${name}_rw.prmtop
mol addfile ../${name}_rw_${run}_1000ns.nc first 2500 step 10 waitfor all
source rename_chainID_${name}.tcl

set nframes [expr  [molinfo top get numframes] - 1 ]

set outfile [open ./rmsd/${name}_dna_rmsd_${run}.dat w]	

#align to the first frame
set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"

set frame0_all [atomselect top "all and noh" frame 0]
set sel_all [atomselect top "all and noh"]

set frame0_alpha_core_ca [atomselect top $alpha_core_ca frame 0]
set sel_alpha_core_ca [atomselect top $alpha_core_ca]


set chi [atomselect top "segname CHI"]	 
set chj [atomselect top "segname CHJ"]	 
# rmsd calculation loop


puts $outfile "DNA RMSD profile"
puts $outfile "Calculated using VMD"
puts $outfile "Time, ps"
puts $outfile "RMSD, A"
puts -nonewline $outfile "Time"
for { set r1 -93 } { $r1<=93 } { incr r1 } {
puts -nonewline $outfile "\t$r1"
}
puts -nonewline $outfile "\n"

for { set i 1 } { $i<=$nframes } { incr i } {
set time [expr 0.2 * $i]
puts -nonewline $outfile [format "%.3f" "$time"]
puts  [format "Time %.3f" "$time"]

for { set r1 -93 } { $r1<=93 } { incr r1 } {

$sel_all frame $i
$sel_alpha_core_ca frame $i

#alignemt to alpha_core_ca the first frame
set trans_mat [measure fit $sel_alpha_core_ca $frame0_alpha_core_ca]
$sel_all move $trans_mat

set r2 [expr $r1 * (-1)]
set sel2 [atomselect top "((segname CHI and resid '$r1') or (segname CHJ and resid '$r2')) and noh" frame $i]
set sel0 [atomselect top "((segname CHI and resid '$r1') or (segname CHJ and resid '$r2')) and noh" frame 0]

set rmsd [measure rmsd $sel2 $sel0]


puts -nonewline $outfile "\t$rmsd"

$sel2 delete
$sel0 delete

}
puts -nonewline $outfile "\n"
}

close $outfile 

}

exit



