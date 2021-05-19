foreach name  {1aoi 1eqz ext2 sym} {
mol load parm7 ../${name}_rw.prmtop 
mol addfile ../${name}_rw_run1_4000ns.nc first 2500 last 190000 step 10 waitfor all
mol addfile ../${name}_rw_run2_500ns.nc first 2500 last 24000 step 10 waitfor all
mol addfile ../${name}_rw_run3_500ns.nc first 2500 last 24000 step 10 waitfor all
mol addfile ../${name}_rw_run4_500ns.nc first 2500 last 24000 step 10 waitfor all
mol addfile ../${name}_rw_run5_500ns.nc first 2500 last 24000 step 10 waitfor all

source rename_chainID_${name}.tcl
package require hbonds

set nframes [expr  [molinfo top get numframes] - 1 ]

set outfile1 [open ./histone_contacts/tail_contacts_${name}_h3_A.dat a]
set outfile2 [open ./histone_contacts/tail_contacts_${name}_h3_E.dat a]
set outfile3 [open ./histone_contacts/tail_contacts_${name}_h4_B.dat a]
set outfile4 [open ./histone_contacts/tail_contacts_${name}_h4_F.dat a]
set outfile5 [open ./histone_contacts/tail_contacts_${name}_h2a_C.dat a]
set outfile6 [open ./histone_contacts/tail_contacts_${name}_h2a_G.dat a]
set outfile7 [open ./histone_contacts/tail_contacts_${name}_h2b_D.dat a]
set outfile8 [open ./histone_contacts/tail_contacts_${name}_h2b_H.dat a]
	

puts -nonewline $outfile1 "Time"
puts -nonewline $outfile2 "Time"
puts -nonewline $outfile3 "Time"
puts -nonewline $outfile4 "Time"
puts -nonewline $outfile5 "Time"
puts -nonewline $outfile6 "Time"
puts -nonewline $outfile7 "Time"
puts -nonewline $outfile8 "Time"

#H3
for { set r1 1 } { $r1<=135 } { incr r1 } {
puts -nonewline $outfile1 "\t$r1"
puts -nonewline $outfile2 "\t$r1"
}

#H4
for { set r1 1 } { $r1<=102 } { incr r1 } {
puts -nonewline $outfile3 "\t$r1"
puts -nonewline $outfile4 "\t$r1"
}
#H2A
for { set r1 1 } { $r1<=128 } { incr r1 } {
puts -nonewline $outfile5 "\t$r1"
puts -nonewline $outfile6 "\t$r1"
}
#H2B
for { set r1 1 } { $r1<=122 } { incr r1 } {
puts -nonewline $outfile7 "\t$r1"
puts -nonewline $outfile8 "\t$r1"
}

puts -nonewline $outfile1 "\n"
puts -nonewline $outfile2 "\n"
puts -nonewline $outfile3 "\n"
puts -nonewline $outfile4 "\n"
puts -nonewline $outfile5 "\n"
puts -nonewline $outfile6 "\n"
puts -nonewline $outfile7 "\n"
puts -nonewline $outfile8 "\n"

for { set i 0 } { $i<=$nframes } { incr i } {
set time [expr 1 * $i]
puts -nonewline $outfile1 [format "%.3f" "$time"]
puts -nonewline $outfile2 [format "%.3f" "$time"]
puts -nonewline $outfile3 [format "%.3f" "$time"]
puts -nonewline $outfile4 [format "%.3f" "$time"]
puts -nonewline $outfile5 [format "%.3f" "$time"]
puts -nonewline $outfile6 [format "%.3f" "$time"]
puts -nonewline $outfile7 [format "%.3f" "$time"]
puts -nonewline $outfile8 [format "%.3f" "$time"]
puts  [format "Time %.3f" "$time"]
#contacts for H3
set sel0 [atomselect top "((segname CHA CHE and resid 1 to 36) or (segname CHB  CHF and resid 1 to 15) or (segname CHC CHG and resid 1 to 11 119 to 128) or (segname CHD CHH and resid 1 to 23)) and noh"  frame $i]
for { set r1 1 } { $r1<=135 } { incr r1 } {
set sel1 [atomselect top "(segname CHA and resid '$r1') and noh" frame $i]
set sel2 [atomselect top "(segname CHE and resid '$r1') and noh" frame $i]
set contacts1 [lindex [measure contacts  4.0 $sel0 $sel1] 1]
set contacts2 [lindex [measure contacts  4.0 $sel0 $sel2] 1]
set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
puts -nonewline $outfile1 "\t$nc1"
puts -nonewline $outfile2 "\t$nc2"
$sel1 delete
$sel2 delete
unset contacts1
unset contacts2
unset nc1
unset nc2
}
#contacts for H4
for { set r1 1 } { $r1<=102 } { incr r1 } {
set sel1 [atomselect top "(segname CHB and resid '$r1') and noh" frame $i]
set sel2 [atomselect top "(segname CHF and resid '$r1') and noh" frame $i]
set contacts1 [lindex [measure contacts  4.0 $sel0 $sel1] 1]
set contacts2 [lindex [measure contacts  4.0 $sel0 $sel2] 1]
set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
puts -nonewline $outfile3 "\t$nc1"
puts -nonewline $outfile4 "\t$nc2"
$sel1 delete
$sel2 delete
unset contacts1
unset contacts2
unset nc1
unset nc2
}
#contacts for H2A
for { set r1 1 } { $r1<=128 } { incr r1 } {
set sel1 [atomselect top "(segname CHC and resid '$r1') and noh" frame $i]
set sel2 [atomselect top "(segname CHG and resid '$r1') and noh" frame $i]
set contacts1 [lindex [measure contacts  4.0 $sel0 $sel1] 1]
set contacts2 [lindex [measure contacts  4.0 $sel0 $sel2] 1]
set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
puts -nonewline $outfile5 "\t$nc1"
puts -nonewline $outfile6 "\t$nc2"
$sel1 delete
$sel2 delete
unset contacts1
unset contacts2
unset nc1
unset nc2
}

#contacts for H2B
for { set r1 1 } { $r1<=122 } { incr r1 } {
set sel1 [atomselect top "(segname CHD and resid '$r1') and noh" frame $i]
set sel2 [atomselect top "(segname CHH and resid '$r1') and noh" frame $i]
set contacts1 [lindex [measure contacts  4.0 $sel0 $sel1] 1]
set contacts2 [lindex [measure contacts  4.0 $sel0 $sel2] 1]
set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
puts -nonewline $outfile7 "\t$nc1"
puts -nonewline $outfile8 "\t$nc2"
$sel1 delete
$sel2 delete
unset contacts1
unset contacts2
unset nc1 
unset nc2
}

puts -nonewline $outfile1 "\n"
puts -nonewline $outfile2 "\n"
puts -nonewline $outfile3 "\n"
puts -nonewline $outfile4 "\n"
puts -nonewline $outfile5 "\n"
puts -nonewline $outfile6 "\n"
puts -nonewline $outfile7 "\n"
puts -nonewline $outfile8 "\n"
}

close $outfile1 
close $outfile2
close $outfile3
close $outfile4
close $outfile5
close $outfile6
close $outfile7
close $outfile8

mol delete all
}
exit


