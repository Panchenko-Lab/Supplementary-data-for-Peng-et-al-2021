set name 1aoi
mol load parm7 ../${name}_rw.prmtop 
mol addfile ../${name}_rw_run1_4500ns.nc first 10000 last 225000 step 10 waitfor all
mol addfile ../${name}_rw_run2_800ns.nc first 2500 last 40000 step 10 waitfor all
mol addfile ../${name}_rw_run3_800ns.nc first 2500 last 40000 step 10 waitfor all
mol addfile ../${name}_rw_run4_800ns.nc first 2500 last 40000 step 10 waitfor all
mol addfile ../${name}_rw_run5_800ns.nc first 2500 last 40000 step 10 waitfor all

source rename_chainID_${name}.tcl

set chi [atomselect top "segname CHI"]	 
set chj [atomselect top "segname CHJ"]	 
# rmsd calculation loop
package require hbonds

set nframes [expr  [molinfo top get numframes] - 1 ]

set outfile0 [open ./DNA_contacts/dna_prot_contacts_${name}_all_LYS.dat a]
set outfile1 [open ./DNA_contacts/dna_prot_contacts_${name}_h3_LYS.dat a]
set outfile2 [open ./DNA_contacts/dna_prot_contacts_${name}_h4_LYS.dat a]
set outfile3 [open ./DNA_contacts/dna_prot_contacts_${name}_h2a_LYS.dat a]
set outfile4 [open ./DNA_contacts/dna_prot_contacts_${name}_h2b_LYS.dat a]
	
#set sel1 [atomselect top "segname CHI and resid '$r1'" frame $i]

puts -nonewline $outfile0 "Time"
puts -nonewline $outfile1 "Time"
puts -nonewline $outfile2 "Time"
puts -nonewline $outfile3 "Time"
puts -nonewline $outfile4 "Time"

for { set r1 -93 } { $r1<=93 } { incr r1 } {
puts -nonewline $outfile0 "\t$r1"
puts -nonewline $outfile1 "\t$r1"
puts -nonewline $outfile2 "\t$r1"
puts -nonewline $outfile3 "\t$r1"
puts -nonewline $outfile4 "\t$r1"
}
puts -nonewline $outfile0 "\n"
puts -nonewline $outfile1 "\n"
puts -nonewline $outfile2 "\n"
puts -nonewline $outfile3 "\n"
puts -nonewline $outfile4 "\n"
for { set i 0 } { $i<=$nframes } { incr i } {
set time [expr 1 * $i]
puts -nonewline $outfile0 [format "%.3f" "$time"]
puts -nonewline $outfile1 [format "%.3f" "$time"]
puts -nonewline $outfile2 [format "%.3f" "$time"]
puts -nonewline $outfile3 [format "%.3f" "$time"]
puts -nonewline $outfile4 [format "%.3f" "$time"]
puts  [format "Time %.3f" "$time"]
unset time


for { set r1 -93 } { $r1<=93 } { incr r1 } {
#set sel1 [atomselect top "protein and noh" frame $i]
set sel0 [atomselect top "((segname CHA CHE and resid 4 9 14 18 23 27 36) or (segname CHB  CHF and resid 5 8 12 16 20) or (segname CHC CHG and resid 5 9 13 119 124 126 128) or (segname CHD CHH and resid 2 8 9 12 13 17 20 21)) and noh"  frame $i]
set sel1 [atomselect top "(segname CHA CHE and resid 4 9 14 18 23 27 36) and noh"  frame $i]
set sel2 [atomselect top "(segname CHB CHF and resid 5 8 12 16 20) and noh"  frame $i]
set sel3 [atomselect top "(segname CHC CHG and resid 5 9 13 119 124 126 128) and noh"  frame $i]
set sel4 [atomselect top "(segname CHD CHH and resid 2 8 9 12 13 17 20 21) and noh"  frame $i]

set r2 [expr $r1 * (-1)]
set sel5 [atomselect top "((segname CHI and resid '$r1') or (segname CHJ and resid '$r2')) and noh" frame $i]
set contacts0 [lindex [measure contacts  4.0 $sel0 $sel5] 1]
set contacts1 [lindex [measure contacts  4.0 $sel1 $sel5] 1]
set contacts2 [lindex [measure contacts  4.0 $sel2 $sel5] 1]
set contacts3 [lindex [measure contacts  4.0 $sel3 $sel5] 1]
set contacts4 [lindex [measure contacts  4.0 $sel4 $sel5] 1]

set nc0 [llength $contacts0]
set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
set nc3 [llength $contacts3]
set nc4 [llength $contacts4]

#puts "$contacts"
#puts [format "Check: %d" "$delta"]
puts -nonewline $outfile0 "\t$nc0"
puts -nonewline $outfile1 "\t$nc1"
puts -nonewline $outfile2 "\t$nc2"
puts -nonewline $outfile3 "\t$nc3"
puts -nonewline $outfile4 "\t$nc4"
$sel0 delete
$sel1 delete
$sel2 delete
$sel3 delete
$sel4 delete
$sel5 delete
unset contacts0 
unset contacts1 
unset contacts2 
unset contacts3 
unset contacts4 
unset nc0 
unset nc1 
unset nc2 
unset nc3
unset nc4
#hbonds -sel1 $sel1 -sel2 $sel2 -dist 3.0 -ang 30 -frames $i:$i -writefile yes -outfile ../analysis_data/dna_hbonds.dat
#-type unique -writefile yes -outfile ../analysis_data/dna_hbonds.dat

}
puts -nonewline $outfile0 "\n"
puts -nonewline $outfile1 "\n"
puts -nonewline $outfile2 "\n"
puts -nonewline $outfile3 "\n"
puts -nonewline $outfile4 "\n"
}

close $outfile0
close $outfile1 
close $outfile2
close $outfile3
close $outfile4

#hbonds -sel1 $chi -sel2 $chj -dist 3.0 -ang 30 -type unique -writefile yes -outfile ../analysis_data/dna_hbonds.dat

exit


