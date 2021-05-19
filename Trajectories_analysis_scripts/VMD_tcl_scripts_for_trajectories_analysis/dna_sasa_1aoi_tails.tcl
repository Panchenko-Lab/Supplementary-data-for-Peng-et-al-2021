foreach input  {1aoi} {

mol load parm7 ../${input}_rw.prmtop
mol addfile ../${input}_rw_run1_4500ns.nc first 10000 last 225000 step 250 waitfor all
mol addfile ../${input}_rw_run2_800ns.nc first 2500 last 40000 step 250 waitfor all
mol addfile ../${input}_rw_run3_800ns.nc first 2500 last 40000 step 250 waitfor all
mol addfile ../${input}_rw_run4_800ns.nc first 2500 last 40000 step 250 waitfor all
mol addfile ../${input}_rw_run5_800ns.nc first 2500 last 40000 step 250 waitfor all

source rename_chainID_${input}.tcl

set chi [atomselect top "segname CHI"]	 
set chj [atomselect top "segname CHJ"]	 
# rmsd calculation loop
package require hbonds

set nframes [expr  [molinfo top get numframes] - 1 ]


set outfile0 [open ./sasa/rsasa_all_${input}.dat a]
set outfile1 [open ./sasa/dsasa_all_${input}.dat a]	
set outfile2 [open ./sasa/dsasa_H3_${input}.dat a]
set outfile3 [open ./sasa/dsasa_H4_${input}.dat a]
set outfile4 [open ./sasa/dsasa_H2A_${input}.dat a]
set outfile5 [open ./sasa/dsasa_H2B_${input}.dat a]
#set sel1 [atomselect top "segname CHI and resid '$r1'" frame $i]
puts -nonewline $outfile0 "Time"
puts -nonewline $outfile1 "Time"
puts -nonewline $outfile2 "Time"
puts -nonewline $outfile3 "Time"
puts -nonewline $outfile4 "Time"
puts -nonewline $outfile5 "Time"
for { set r1 -93 } { $r1<=93 } { incr r1 } {
puts -nonewline $outfile0 "\t$r1"
puts -nonewline $outfile1 "\t$r1"
puts -nonewline $outfile2 "\t$r1"
puts -nonewline $outfile3 "\t$r1"
puts -nonewline $outfile4 "\t$r1"
puts -nonewline $outfile5 "\t$r1"
}
puts -nonewline $outfile0 "\n"
puts -nonewline $outfile1 "\n"
puts -nonewline $outfile2 "\n"
puts -nonewline $outfile3 "\n"
puts -nonewline $outfile4 "\n"
puts -nonewline $outfile5 "\n"
for { set i 0 } { $i<=$nframes } { incr i } {
set time [expr 1 * $i]
puts -nonewline $outfile0 [format "%.3f" "$time"]
puts -nonewline $outfile1 [format "%.3f" "$time"]
puts -nonewline $outfile2 [format "%.3f" "$time"]
puts -nonewline $outfile3 [format "%.3f" "$time"]
puts -nonewline $outfile4 [format "%.3f" "$time"]
puts -nonewline $outfile5 [format "%.3f" "$time"]
puts  [format "Time %.3f" "$time"]
for { set r1 -93 } { $r1<=93 } { incr r1 } {
#set sel0 [atomselect top "(nucleic or protein) and not ((segname CHA CHE and resid 1 to 38) or (segname CHB  CHF and resid 1 to 25) or (segname CHC CHG and resid 1 to 15 119 to 128) or (segname CHD CHH and resid 1 to 28))" frame $i]
set sel0 [atomselect top "all" frame $i]
set sel1 [atomselect top "(all not ((segname CHA CHE and resid 1 to 36) or (segname CHB  CHF and resid 1 to 20) or (segname CHC CHG and resid 1 to 13 119 to 128) or (segname CHD CHH and resid 1 to 23)))" frame $i]
set sel2 [atomselect top "all not (segname CHA CHE and resid 1 to 36)"  frame $i]
set sel3 [atomselect top "all not (segname CHB  CHF and resid 1 to 20)"  frame $i]
set sel4 [atomselect top "all not (segname CHC CHG and resid 1 to 13 119 to 128)"  frame $i]
set sel5 [atomselect top "all not (segname CHD CHH and resid 1 to 23)"  frame $i]

set r2 [expr $r1 * (-1)]
set sel6 [atomselect top "((segname CHI and resid '$r1') or (segname CHJ and resid '$r2'))" frame $i]
set sasa0 [measure sasa 1.4  $sel0 -restrict $sel6]
set sasa1 [measure sasa 1.4  $sel1 -restrict $sel6]
set sasa2 [measure sasa 1.4  $sel2 -restrict $sel6]
set sasa3 [measure sasa 1.4  $sel3 -restrict $sel6]
set sasa4 [measure sasa 1.4  $sel4 -restrict $sel6]
set sasa5 [measure sasa 1.4  $sel5 -restrict $sel6]

set rsasa [expr 1 - $sasa0/$sasa1]
set dsasa [expr $sasa1 - $sasa0]
set h3 [expr $sasa2 - $sasa0]
set h4 [expr $sasa3 - $sasa0]
set h2a [expr $sasa4 - $sasa0]
set h2b [expr $sasa5 - $sasa0]

puts -nonewline $outfile0 "\t$rsasa"
puts -nonewline $outfile1 "\t$dsasa"
puts -nonewline $outfile2 "\t$h3"
puts -nonewline $outfile3 "\t$h4"
puts -nonewline $outfile4 "\t$h2a"
puts -nonewline $outfile5 "\t$h2b"

$sel0 delete
$sel1 delete
$sel2 delete
$sel3 delete
$sel4 delete
$sel5 delete
$sel6 delete
unset sasa0 
unset sasa1 
unset sasa2 
unset sasa3 
unset sasa4 
unset sasa5 
unset rsasa 
unset dsasa 
unset h3 
unset h4 
unset h2a 
unset h2b 

}
puts -nonewline $outfile0 "\n"
puts -nonewline $outfile1 "\n"
puts -nonewline $outfile2 "\n"
puts -nonewline $outfile3 "\n"
puts -nonewline $outfile4 "\n"
puts -nonewline $outfile5 "\n"

}

close $outfile0
close $outfile1
close $outfile2
close $outfile3
close $outfile4
close $outfile5


#hbonds -sel1 $chi -sel2 $chj -dist 3.0 -ang 30 -type unique -writefile yes -outfile ../analysis_data/dna_hbonds.dat

mol delete all
}

exit
