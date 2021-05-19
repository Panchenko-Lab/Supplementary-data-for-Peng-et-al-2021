foreach name  {1aoi} {

mol delete all
mol load parm7 ../${name}_rw.prmtop
mol addfile ../${name}_rw_run1_4500ns.nc first 10000 last 225000 step 10 waitfor all
mol addfile ../${name}_rw_run2_800ns.nc first 2500 last 40000 step 10 waitfor all
mol addfile ../${name}_rw_run3_800ns.nc first 2500 last 40000 step 10 waitfor all
mol addfile ../${name}_rw_run4_800ns.nc first 2500 last 40000 step 10 waitfor all
mol addfile ../${name}_rw_run5_800ns.nc first 2500 last 40000 step 10 waitfor all


source rename_chainID_${name}.tcl
package require hbonds

set nframes [expr  [molinfo top get numframes] - 1 ]

set outfile11 [open ./tail_contacts/tail_contacts_all_${name}_h3_A.dat w]
set outfile12 [open ./tail_contacts/tail_contacts_all_${name}_h3_E.dat w]
set outfile13 [open ./tail_contacts/tail_contacts_backbone_${name}_h3_A.dat w]
set outfile14 [open ./tail_contacts/tail_contacts_backbone_${name}_h3_E.dat w]
set outfile15 [open ./tail_contacts/tail_contacts_base_${name}_h3_A.dat w]
set outfile16 [open ./tail_contacts/tail_contacts_base_${name}_h3_E.dat w]

set outfile21 [open ./tail_contacts/tail_contacts_all_${name}_h4_B.dat w]
set outfile22 [open ./tail_contacts/tail_contacts_all_${name}_h4_F.dat w]
set outfile23 [open ./tail_contacts/tail_contacts_backbone_${name}_h4_B.dat w]
set outfile24 [open ./tail_contacts/tail_contacts_backbone_${name}_h4_F.dat w]
set outfile25 [open ./tail_contacts/tail_contacts_base_${name}_h4_B.dat w]
set outfile26 [open ./tail_contacts/tail_contacts_base_${name}_h4_F.dat w]

set outfile31 [open ./tail_contacts/tail_contacts_all_${name}_h2a_C.dat w]
set outfile32 [open ./tail_contacts/tail_contacts_all_${name}_h2a_G.dat w]
set outfile33 [open ./tail_contacts/tail_contacts_backbone_${name}_h2a_C.dat w]
set outfile34 [open ./tail_contacts/tail_contacts_backbone_${name}_h2a_G.dat w]
set outfile35 [open ./tail_contacts/tail_contacts_base_${name}_h2a_C.dat w]
set outfile36 [open ./tail_contacts/tail_contacts_base_${name}_h2a_G.dat w]

set outfile41 [open ./tail_contacts/tail_contacts_all_${name}_h2b_D.dat w]
set outfile42 [open ./tail_contacts/tail_contacts_all_${name}_h2b_H.dat w]
set outfile43 [open ./tail_contacts/tail_contacts_backbone_${name}_h2b_D.dat w]
set outfile44 [open ./tail_contacts/tail_contacts_backbone_${name}_h2b_H.dat w]
set outfile45 [open ./tail_contacts/tail_contacts_base_${name}_h2b_D.dat w]
set outfile46 [open ./tail_contacts/tail_contacts_base_${name}_h2b_H.dat w]




puts -nonewline $outfile11 "Time"
puts -nonewline $outfile12 "Time"
puts -nonewline $outfile13 "Time"
puts -nonewline $outfile14 "Time"
puts -nonewline $outfile15 "Time"
puts -nonewline $outfile16 "Time"

puts -nonewline $outfile21 "Time"
puts -nonewline $outfile22 "Time"
puts -nonewline $outfile23 "Time"
puts -nonewline $outfile24 "Time"
puts -nonewline $outfile25 "Time"
puts -nonewline $outfile26 "Time"

puts -nonewline $outfile31 "Time"
puts -nonewline $outfile32 "Time"
puts -nonewline $outfile33 "Time"
puts -nonewline $outfile34 "Time"
puts -nonewline $outfile35 "Time"
puts -nonewline $outfile36 "Time"

puts -nonewline $outfile41 "Time"
puts -nonewline $outfile42 "Time"
puts -nonewline $outfile43 "Time"
puts -nonewline $outfile44 "Time"
puts -nonewline $outfile45 "Time"
puts -nonewline $outfile46 "Time"


for { set r1 1 } { $r1<=36 } { incr r1 } {
puts -nonewline $outfile11 "\t$r1"
puts -nonewline $outfile12 "\t$r1"
puts -nonewline $outfile13 "\t$r1"
puts -nonewline $outfile14 "\t$r1"
puts -nonewline $outfile15 "\t$r1"
puts -nonewline $outfile16 "\t$r1"
}
for { set r1 1 } { $r1<=20 } { incr r1 } {
puts -nonewline $outfile21 "\t$r1"
puts -nonewline $outfile22 "\t$r1"
puts -nonewline $outfile23 "\t$r1"
puts -nonewline $outfile24 "\t$r1"
puts -nonewline $outfile25 "\t$r1"
puts -nonewline $outfile26 "\t$r1"
}
for { set r1 1 } { $r1<=13 } { incr r1 } {
puts -nonewline $outfile31 "\t$r1"
puts -nonewline $outfile32 "\t$r1"
puts -nonewline $outfile33 "\t$r1"
puts -nonewline $outfile34 "\t$r1"
puts -nonewline $outfile35 "\t$r1"
puts -nonewline $outfile36 "\t$r1"
}
for { set r1 119 } { $r1<=128 } { incr r1 } {
puts -nonewline $outfile31 "\t$r1"
puts -nonewline $outfile32 "\t$r1"
puts -nonewline $outfile33 "\t$r1"
puts -nonewline $outfile34 "\t$r1"
puts -nonewline $outfile35 "\t$r1"
puts -nonewline $outfile36 "\t$r1"
}
for { set r1 1 } { $r1<=23 } { incr r1 } {
puts -nonewline $outfile41 "\t$r1"
puts -nonewline $outfile42 "\t$r1"
puts -nonewline $outfile43 "\t$r1"
puts -nonewline $outfile44 "\t$r1"
puts -nonewline $outfile45 "\t$r1"
puts -nonewline $outfile46 "\t$r1"
}


puts -nonewline $outfile11 "\n"
puts -nonewline $outfile12 "\n"
puts -nonewline $outfile13 "\n"
puts -nonewline $outfile14 "\n"
puts -nonewline $outfile15 "\n"
puts -nonewline $outfile16 "\n"

puts -nonewline $outfile21 "\n"
puts -nonewline $outfile22 "\n"
puts -nonewline $outfile23 "\n"
puts -nonewline $outfile24 "\n"
puts -nonewline $outfile25 "\n"
puts -nonewline $outfile26 "\n"

puts -nonewline $outfile31 "\n"
puts -nonewline $outfile32 "\n"
puts -nonewline $outfile33 "\n"
puts -nonewline $outfile34 "\n"
puts -nonewline $outfile35 "\n"
puts -nonewline $outfile36 "\n"

puts -nonewline $outfile41 "\n"
puts -nonewline $outfile42 "\n"
puts -nonewline $outfile43 "\n"
puts -nonewline $outfile44 "\n"
puts -nonewline $outfile45 "\n"
puts -nonewline $outfile46 "\n"


for { set i 0 } { $i<=$nframes } { incr i } {
set time [expr 0.1 * $i]
puts -nonewline $outfile11 [format "%.3f" "$time"]
puts -nonewline $outfile12 [format "%.3f" "$time"]
puts -nonewline $outfile13 [format "%.3f" "$time"]
puts -nonewline $outfile14 [format "%.3f" "$time"]
puts -nonewline $outfile15 [format "%.3f" "$time"]
puts -nonewline $outfile16 [format "%.3f" "$time"]

puts -nonewline $outfile21 [format "%.3f" "$time"]
puts -nonewline $outfile22 [format "%.3f" "$time"]
puts -nonewline $outfile23 [format "%.3f" "$time"]
puts -nonewline $outfile24 [format "%.3f" "$time"]
puts -nonewline $outfile25 [format "%.3f" "$time"]
puts -nonewline $outfile26 [format "%.3f" "$time"]

puts -nonewline $outfile31 [format "%.3f" "$time"]
puts -nonewline $outfile32 [format "%.3f" "$time"]
puts -nonewline $outfile33 [format "%.3f" "$time"]
puts -nonewline $outfile34 [format "%.3f" "$time"]
puts -nonewline $outfile35 [format "%.3f" "$time"]
puts -nonewline $outfile36 [format "%.3f" "$time"]

puts -nonewline $outfile41 [format "%.3f" "$time"]
puts -nonewline $outfile42 [format "%.3f" "$time"]
puts -nonewline $outfile43 [format "%.3f" "$time"]
puts -nonewline $outfile44 [format "%.3f" "$time"]
puts -nonewline $outfile45 [format "%.3f" "$time"]
puts -nonewline $outfile46 [format "%.3f" "$time"]

puts  [format "Time %.3f" "$time"]
#contacts for H3
for { set r1 1 } { $r1<=36 } { incr r1 } {
set sel1 [atomselect top "nucleic and noh"  frame $i]
set sel2 [atomselect top "nucleic and backbone and noh"  frame $i]
set sel3 [atomselect top "nucleic not backbone and noh"  frame $i]

set sel4 [atomselect top "(segname CHA and resid '$r1') and noh" frame $i]
set sel5 [atomselect top "(segname CHE and resid '$r1') and noh" frame $i]

set contacts1 [lindex [measure contacts  4.0 $sel1 $sel4] 1]
set contacts2 [lindex [measure contacts  4.0 $sel1 $sel5] 1]

set contacts3 [lindex [measure contacts  4.0 $sel2 $sel4] 1]
set contacts4 [lindex [measure contacts  4.0 $sel2 $sel5] 1]

set contacts5 [lindex [measure contacts  4.0 $sel3 $sel4] 1]
set contacts6 [lindex [measure contacts  4.0 $sel3 $sel5] 1]

set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
set nc3 [llength $contacts3]
set nc4 [llength $contacts4]
set nc5 [llength $contacts5]
set nc6 [llength $contacts6]

puts -nonewline $outfile11 "\t$nc1"
puts -nonewline $outfile12 "\t$nc2"
puts -nonewline $outfile13 "\t$nc3"
puts -nonewline $outfile14 "\t$nc4"
puts -nonewline $outfile15 "\t$nc5"
puts -nonewline $outfile16 "\t$nc6"

$sel1 delete
$sel2 delete
$sel3 delete
$sel4 delete
$sel5 delete

unset contacts1 
unset contacts2
unset contacts3
unset contacts4
unset contacts5
unset contacts6 
unset nc1  
unset nc2 
unset nc3
unset nc4
unset nc5
unset nc6
}
#contacts for H4
for { set r1 1 } { $r1<=20 } { incr r1 } {
set sel1 [atomselect top "nucleic and noh"  frame $i]
set sel2 [atomselect top "nucleic and backbone and noh"  frame $i]
set sel3 [atomselect top "nucleic not backbone and noh"  frame $i]

set sel4 [atomselect top "(segname CHB and resid '$r1') and noh" frame $i]
set sel5 [atomselect top "(segname CHF and resid '$r1') and noh" frame $i]

set contacts1 [lindex [measure contacts  4.0 $sel1 $sel4] 1]
set contacts2 [lindex [measure contacts  4.0 $sel1 $sel5] 1]

set contacts3 [lindex [measure contacts  4.0 $sel2 $sel4] 1]
set contacts4 [lindex [measure contacts  4.0 $sel2 $sel5] 1]

set contacts5 [lindex [measure contacts  4.0 $sel3 $sel4] 1]
set contacts6 [lindex [measure contacts  4.0 $sel3 $sel5] 1]

set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
set nc3 [llength $contacts3]
set nc4 [llength $contacts4]
set nc5 [llength $contacts5]
set nc6 [llength $contacts6]

puts -nonewline $outfile21 "\t$nc1"
puts -nonewline $outfile22 "\t$nc2"
puts -nonewline $outfile23 "\t$nc3"
puts -nonewline $outfile24 "\t$nc4"
puts -nonewline $outfile25 "\t$nc5"
puts -nonewline $outfile26 "\t$nc6"

$sel1 delete
$sel2 delete
$sel3 delete
$sel4 delete
$sel5 delete

unset contacts1
unset contacts2
unset contacts3
unset contacts4
unset contacts5
unset contacts6
unset nc1  
unset nc2
unset nc3
unset nc4
unset nc5
unset nc6

}
#contacts for H2A
for { set r1 1 } { $r1<=13 } { incr r1 } {
set sel1 [atomselect top "nucleic and noh"  frame $i]
set sel2 [atomselect top "nucleic and backbone and noh"  frame $i]
set sel3 [atomselect top "nucleic not backbone and noh"  frame $i]

set sel4 [atomselect top "(segname CHC and resid '$r1') and noh" frame $i]
set sel5 [atomselect top "(segname CHG and resid '$r1') and noh" frame $i]

set contacts1 [lindex [measure contacts  4.0 $sel1 $sel4] 1]
set contacts2 [lindex [measure contacts  4.0 $sel1 $sel5] 1]

set contacts3 [lindex [measure contacts  4.0 $sel2 $sel4] 1]
set contacts4 [lindex [measure contacts  4.0 $sel2 $sel5] 1]

set contacts5 [lindex [measure contacts  4.0 $sel3 $sel4] 1]
set contacts6 [lindex [measure contacts  4.0 $sel3 $sel5] 1]

set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
set nc3 [llength $contacts3]
set nc4 [llength $contacts4]
set nc5 [llength $contacts5]
set nc6 [llength $contacts6]

puts -nonewline $outfile31 "\t$nc1"
puts -nonewline $outfile32 "\t$nc2"
puts -nonewline $outfile33 "\t$nc3"
puts -nonewline $outfile34 "\t$nc4"
puts -nonewline $outfile35 "\t$nc5"
puts -nonewline $outfile36 "\t$nc6"

$sel1 delete
$sel2 delete
$sel3 delete
$sel4 delete
$sel5 delete

unset contacts1
unset contacts2
unset contacts3
unset contacts4
unset contacts5
unset contacts6
unset nc1  
unset nc2
unset nc3
unset nc4
unset nc5
unset nc6

}

for { set r1 119 } { $r1<=128 } { incr r1 } {
set sel1 [atomselect top "nucleic and noh"  frame $i]
set sel2 [atomselect top "nucleic and backbone and noh"  frame $i]
set sel3 [atomselect top "nucleic not backbone and noh"  frame $i]

set sel4 [atomselect top "(segname CHC and resid '$r1') and noh" frame $i]
set sel5 [atomselect top "(segname CHG and resid '$r1') and noh" frame $i]

set contacts1 [lindex [measure contacts  4.0 $sel1 $sel4] 1]
set contacts2 [lindex [measure contacts  4.0 $sel1 $sel5] 1]

set contacts3 [lindex [measure contacts  4.0 $sel2 $sel4] 1]
set contacts4 [lindex [measure contacts  4.0 $sel2 $sel5] 1]

set contacts5 [lindex [measure contacts  4.0 $sel3 $sel4] 1]
set contacts6 [lindex [measure contacts  4.0 $sel3 $sel5] 1]

set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
set nc3 [llength $contacts3]
set nc4 [llength $contacts4]
set nc5 [llength $contacts5]
set nc6 [llength $contacts6]

puts -nonewline $outfile31 "\t$nc1"
puts -nonewline $outfile32 "\t$nc2"
puts -nonewline $outfile33 "\t$nc3"
puts -nonewline $outfile34 "\t$nc4"
puts -nonewline $outfile35 "\t$nc5"
puts -nonewline $outfile36 "\t$nc6"

$sel1 delete
$sel2 delete
$sel3 delete
$sel4 delete
$sel5 delete

unset contacts1
unset contacts2
unset contacts3
unset contacts4
unset contacts5
unset contacts6
unset nc1  
unset nc2
unset nc3
unset nc4
unset nc5
unset nc6

}


#contacts for H2B
for { set r1 1 } { $r1<=23 } { incr r1 } {
set sel1 [atomselect top "nucleic and noh"  frame $i]
set sel2 [atomselect top "nucleic and backbone and noh"  frame $i]
set sel3 [atomselect top "nucleic not backbone and noh"  frame $i]

set sel4 [atomselect top "(segname CHD and resid '$r1') and noh" frame $i]
set sel5 [atomselect top "(segname CHH and resid '$r1') and noh" frame $i]

set contacts1 [lindex [measure contacts  4.0 $sel1 $sel4] 1]
set contacts2 [lindex [measure contacts  4.0 $sel1 $sel5] 1]

set contacts3 [lindex [measure contacts  4.0 $sel2 $sel4] 1]
set contacts4 [lindex [measure contacts  4.0 $sel2 $sel5] 1]

set contacts5 [lindex [measure contacts  4.0 $sel3 $sel4] 1]
set contacts6 [lindex [measure contacts  4.0 $sel3 $sel5] 1]

set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
set nc3 [llength $contacts3]
set nc4 [llength $contacts4]
set nc5 [llength $contacts5]
set nc6 [llength $contacts6]

puts -nonewline $outfile41 "\t$nc1"
puts -nonewline $outfile42 "\t$nc2"
puts -nonewline $outfile43 "\t$nc3"
puts -nonewline $outfile44 "\t$nc4"
puts -nonewline $outfile45 "\t$nc5"
puts -nonewline $outfile46 "\t$nc6"

$sel1 delete
$sel2 delete
$sel3 delete
$sel4 delete
$sel5 delete

unset contacts1
unset contacts2
unset contacts3
unset contacts4
unset contacts5
unset contacts6
unset nc1  
unset nc2
unset nc3
unset nc4
unset nc5
unset nc6

}

puts -nonewline $outfile11 "\n"
puts -nonewline $outfile12 "\n"
puts -nonewline $outfile13 "\n"
puts -nonewline $outfile14 "\n"
puts -nonewline $outfile15 "\n"
puts -nonewline $outfile16 "\n"

puts -nonewline $outfile21 "\n"
puts -nonewline $outfile22 "\n"
puts -nonewline $outfile23 "\n"
puts -nonewline $outfile24 "\n"
puts -nonewline $outfile25 "\n"
puts -nonewline $outfile26 "\n"

puts -nonewline $outfile31 "\n"
puts -nonewline $outfile32 "\n"
puts -nonewline $outfile33 "\n"
puts -nonewline $outfile34 "\n"
puts -nonewline $outfile35 "\n"
puts -nonewline $outfile36 "\n"

puts -nonewline $outfile41 "\n"
puts -nonewline $outfile42 "\n"
puts -nonewline $outfile43 "\n"
puts -nonewline $outfile44 "\n"
puts -nonewline $outfile45 "\n"
puts -nonewline $outfile46 "\n"

}

close $outfile11 
close $outfile12
close $outfile13
close $outfile14
close $outfile15
close $outfile16

close $outfile21
close $outfile22
close $outfile23
close $outfile24
close $outfile25
close $outfile26

close $outfile31
close $outfile32
close $outfile33
close $outfile34
close $outfile35
close $outfile36

close $outfile41
close $outfile42
close $outfile43
close $outfile44
close $outfile45
close $outfile46

mol delete all
}
exit

