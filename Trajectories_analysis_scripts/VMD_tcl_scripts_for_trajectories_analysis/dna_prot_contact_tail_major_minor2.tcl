foreach name  {1eqz} {

mol delete all
mol load parm7 ../${name}_rw.prmtop
mol addfile ../${name}_rw_run1_2000ns.nc first 2500 step 20 waitfor all
mol addfile ../${name}_rw_run2_200ns.nc first 2500 step 20 waitfor all
mol addfile ../${name}_rw_run3_200ns.nc first 2500 step 20 waitfor all
mol addfile ../${name}_rw_run4_200ns.nc first 2500 step 20 waitfor all
mol addfile ../${name}_rw_run5_200ns.nc first 2500 step 20 waitfor all

source rename_chainID_${name}.tcl
package require hbonds

set nframes [expr  [molinfo top get numframes] - 1 ]

set outfile11 [open ./tail_contacts/tail_contacts_minor_base_${name}_h3_A.dat w]
set outfile12 [open ./tail_contacts/tail_contacts_minor_base_${name}_h3_E.dat w]
set outfile13 [open ./tail_contacts/tail_contacts_major_base_${name}_h3_A.dat w]
set outfile14 [open ./tail_contacts/tail_contacts_major_base_${name}_h3_E.dat w]

set outfile21 [open ./tail_contacts/tail_contacts_minor_base_${name}_h4_B.dat w]
set outfile22 [open ./tail_contacts/tail_contacts_minor_base_${name}_h4_F.dat w]
set outfile23 [open ./tail_contacts/tail_contacts_major_base_${name}_h4_B.dat w]
set outfile24 [open ./tail_contacts/tail_contacts_major_base_${name}_h4_F.dat w]

set outfile31 [open ./tail_contacts/tail_contacts_minor_base_${name}_h2a_C.dat w]
set outfile32 [open ./tail_contacts/tail_contacts_minor_base_${name}_h2a_G.dat w]
set outfile33 [open ./tail_contacts/tail_contacts_major_base_${name}_h2a_C.dat w]
set outfile34 [open ./tail_contacts/tail_contacts_major_base_${name}_h2a_G.dat w]

set outfile41 [open ./tail_contacts/tail_contacts_minor_base_${name}_h2b_D.dat w]
set outfile42 [open ./tail_contacts/tail_contacts_minor_base_${name}_h2b_H.dat w]
set outfile43 [open ./tail_contacts/tail_contacts_major_base_${name}_h2b_D.dat w]
set outfile44 [open ./tail_contacts/tail_contacts_major_base_${name}_h2b_H.dat w]




puts -nonewline $outfile11 "Time"
puts -nonewline $outfile12 "Time"
puts -nonewline $outfile13 "Time"
puts -nonewline $outfile14 "Time"

puts -nonewline $outfile21 "Time"
puts -nonewline $outfile22 "Time"
puts -nonewline $outfile23 "Time"
puts -nonewline $outfile24 "Time"

puts -nonewline $outfile31 "Time"
puts -nonewline $outfile32 "Time"
puts -nonewline $outfile33 "Time"
puts -nonewline $outfile34 "Time"

puts -nonewline $outfile41 "Time"
puts -nonewline $outfile42 "Time"
puts -nonewline $outfile43 "Time"
puts -nonewline $outfile44 "Time"


for { set r1 1 } { $r1<=36 } { incr r1 } {
puts -nonewline $outfile11 "\t$r1"
puts -nonewline $outfile12 "\t$r1"
puts -nonewline $outfile13 "\t$r1"
puts -nonewline $outfile14 "\t$r1"
}
for { set r1 1 } { $r1<=20 } { incr r1 } {
puts -nonewline $outfile21 "\t$r1"
puts -nonewline $outfile22 "\t$r1"
puts -nonewline $outfile23 "\t$r1"
puts -nonewline $outfile24 "\t$r1"
}
for { set r1 1 } { $r1<=13 } { incr r1 } {
puts -nonewline $outfile31 "\t$r1"
puts -nonewline $outfile32 "\t$r1"
puts -nonewline $outfile33 "\t$r1"
puts -nonewline $outfile34 "\t$r1"
}
for { set r1 119 } { $r1<=128 } { incr r1 } {
puts -nonewline $outfile31 "\t$r1"
puts -nonewline $outfile32 "\t$r1"
puts -nonewline $outfile33 "\t$r1"
puts -nonewline $outfile34 "\t$r1"
}
for { set r1 1 } { $r1<=23 } { incr r1 } {
puts -nonewline $outfile41 "\t$r1"
puts -nonewline $outfile42 "\t$r1"
puts -nonewline $outfile43 "\t$r1"
puts -nonewline $outfile44 "\t$r1"
}


puts -nonewline $outfile11 "\n"
puts -nonewline $outfile12 "\n"
puts -nonewline $outfile13 "\n"
puts -nonewline $outfile14 "\n"

puts -nonewline $outfile21 "\n"
puts -nonewline $outfile22 "\n"
puts -nonewline $outfile23 "\n"
puts -nonewline $outfile24 "\n"

puts -nonewline $outfile31 "\n"
puts -nonewline $outfile32 "\n"
puts -nonewline $outfile33 "\n"
puts -nonewline $outfile34 "\n"

puts -nonewline $outfile41 "\n"
puts -nonewline $outfile42 "\n"
puts -nonewline $outfile43 "\n"
puts -nonewline $outfile44 "\n"


for { set i 0 } { $i<=$nframes } { incr i } {
set time [expr 1 * $i]
puts -nonewline $outfile11 [format "%.3f" "$time"]
puts -nonewline $outfile12 [format "%.3f" "$time"]
puts -nonewline $outfile13 [format "%.3f" "$time"]
puts -nonewline $outfile14 [format "%.3f" "$time"]

puts -nonewline $outfile21 [format "%.3f" "$time"]
puts -nonewline $outfile22 [format "%.3f" "$time"]
puts -nonewline $outfile23 [format "%.3f" "$time"]
puts -nonewline $outfile24 [format "%.3f" "$time"]

puts -nonewline $outfile31 [format "%.3f" "$time"]
puts -nonewline $outfile32 [format "%.3f" "$time"]
puts -nonewline $outfile33 [format "%.3f" "$time"]
puts -nonewline $outfile34 [format "%.3f" "$time"]

puts -nonewline $outfile41 [format "%.3f" "$time"]
puts -nonewline $outfile42 [format "%.3f" "$time"]
puts -nonewline $outfile43 [format "%.3f" "$time"]
puts -nonewline $outfile44 [format "%.3f" "$time"]

puts  [format "Time %.3f" "$time"]

set DNA_minor "((segname CHI and (resid '-90' to '-85' '-80' to '-75' '-70' to '-65' '-60' to '-55' '-50' to '-45' '-40' to '-35' '-30' to '-25'  '-20' to '-15' '-10' to '-5' 0 to 5 10 to 15 20 to 25 30 to 35 40 to 45 50 to 55 60 to 65 70 to 75 80 to 85)) or (segname CHJ and (resid '-85' to '-80' '-75' to '-70' '-65' to '-60' '-55' to '-50' '-45' to '-40' '-35' to '-30' '-25' to '-20'  '-15' to '-10' '-5' to 0 5 to 10 15 to 20 25 to 30 35 to 40 45 to 50 55 to 60 65 to 70 75 to 80 85 to 90)))"

set DNA_major "((segname CHJ and (resid '-90' to '-85' '-80' to '-75' '-70' to '-65' '-60' to '-55' '-50' to '-45' '-40' to '-35' '-30' to '-25'  '-20' to '-15' '-10' to '-5' 0 to 5 10 to 15 20 to 25 30 to 35 40 to 45 50 to 55 60 to 65 70 to 75 80 to 85)) or (segname CHI and (resid '-85' to '-80' '-75' to '-70' '-65' to '-60' '-55' to '-50' '-45' to '-40' '-35' to '-30' '-25' to '-20'  '-15' to '-10' '-5' to 0 5 to 10 15 to 20 25 to 30 35 to 40 45 to 50 55 to 60 65 to 70 75 to 80 85 to 90)))" 
#contacts for H3
for { set r1 1 } { $r1<=36 } { incr r1 } {
set sel1 [atomselect top "nucleic and noh"  frame $i]
set sel2 [atomselect top "$DNA_minor and (not backbone) and noh"  frame $i]
set sel3 [atomselect top "$DNA_major and (not backbone) and noh"  frame $i]

set sel4 [atomselect top "(segname CHA and resid '$r1') and noh" frame $i]
set sel5 [atomselect top "(segname CHE and resid '$r1') and noh" frame $i]


set contacts1 [lindex [measure contacts  4.0 $sel2 $sel4] 1]
set contacts2 [lindex [measure contacts  4.0 $sel2 $sel5] 1]

set contacts3 [lindex [measure contacts  4.0 $sel3 $sel4] 1]
set contacts4 [lindex [measure contacts  4.0 $sel3 $sel5] 1]

set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
set nc3 [llength $contacts3]
set nc4 [llength $contacts4]

puts -nonewline $outfile11 "\t$nc1"
puts -nonewline $outfile12 "\t$nc2"
puts -nonewline $outfile13 "\t$nc3"
puts -nonewline $outfile14 "\t$nc4"

$sel1 delete
$sel2 delete
$sel3 delete
$sel4 delete
$sel5 delete

unset contacts1 
unset contacts2
unset contacts3
unset contacts4
unset nc1  
unset nc2 
unset nc3
unset nc4
}
#contacts for H4
for { set r1 1 } { $r1<=20 } { incr r1 } {
set sel1 [atomselect top "nucleic and noh"  frame $i]
set sel2 [atomselect top "$DNA_minor and (not backbone) and noh"  frame $i]
set sel3 [atomselect top "$DNA_major and (not backbone) and noh"  frame $i]

set sel4 [atomselect top "(segname CHB and resid '$r1') and noh" frame $i]
set sel5 [atomselect top "(segname CHF and resid '$r1') and noh" frame $i]


set contacts1 [lindex [measure contacts  4.0 $sel2 $sel4] 1]
set contacts2 [lindex [measure contacts  4.0 $sel2 $sel5] 1]

set contacts3 [lindex [measure contacts  4.0 $sel3 $sel4] 1]
set contacts4 [lindex [measure contacts  4.0 $sel3 $sel5] 1]

set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
set nc3 [llength $contacts3]
set nc4 [llength $contacts4]

puts -nonewline $outfile21 "\t$nc1"
puts -nonewline $outfile22 "\t$nc2"
puts -nonewline $outfile23 "\t$nc3"
puts -nonewline $outfile24 "\t$nc4"

$sel1 delete
$sel2 delete
$sel3 delete
$sel4 delete
$sel5 delete

unset contacts1
unset contacts2
unset contacts3
unset contacts4
unset nc1  
unset nc2
unset nc3
unset nc4

}
#contacts for H2A
for { set r1 1 } { $r1<=13 } { incr r1 } {
set sel1 [atomselect top "nucleic and noh"  frame $i]
set sel2 [atomselect top "$DNA_minor and (not backbone) and noh"  frame $i]
set sel3 [atomselect top "$DNA_major and (not backbone) and noh"  frame $i]

set sel4 [atomselect top "(segname CHC and resid '$r1') and noh" frame $i]
set sel5 [atomselect top "(segname CHG and resid '$r1') and noh" frame $i]


set contacts1 [lindex [measure contacts  4.0 $sel2 $sel4] 1]
set contacts2 [lindex [measure contacts  4.0 $sel2 $sel5] 1]

set contacts3 [lindex [measure contacts  4.0 $sel3 $sel4] 1]
set contacts4 [lindex [measure contacts  4.0 $sel3 $sel5] 1]

set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
set nc3 [llength $contacts3]
set nc4 [llength $contacts4]

puts -nonewline $outfile31 "\t$nc1"
puts -nonewline $outfile32 "\t$nc2"
puts -nonewline $outfile33 "\t$nc3"
puts -nonewline $outfile34 "\t$nc4"

$sel1 delete
$sel2 delete
$sel3 delete
$sel4 delete
$sel5 delete

unset contacts1
unset contacts2
unset contacts3
unset contacts4
unset nc1  
unset nc2
unset nc3
unset nc4

}

for { set r1 119 } { $r1<=128 } { incr r1 } {
set sel1 [atomselect top "nucleic and noh"  frame $i]
set sel2 [atomselect top "$DNA_minor and (not backbone) and noh"  frame $i]
set sel3 [atomselect top "$DNA_major and (not backbone) and noh"  frame $i]

set sel4 [atomselect top "(segname CHC and resid '$r1') and noh" frame $i]
set sel5 [atomselect top "(segname CHG and resid '$r1') and noh" frame $i]


set contacts1 [lindex [measure contacts  4.0 $sel2 $sel4] 1]
set contacts2 [lindex [measure contacts  4.0 $sel2 $sel5] 1]

set contacts3 [lindex [measure contacts  4.0 $sel3 $sel4] 1]
set contacts4 [lindex [measure contacts  4.0 $sel3 $sel5] 1]

set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
set nc3 [llength $contacts3]
set nc4 [llength $contacts4]

puts -nonewline $outfile31 "\t$nc1"
puts -nonewline $outfile32 "\t$nc2"
puts -nonewline $outfile33 "\t$nc3"
puts -nonewline $outfile34 "\t$nc4"

$sel1 delete
$sel2 delete
$sel3 delete
$sel4 delete
$sel5 delete

unset contacts1
unset contacts2
unset contacts3
unset contacts4
unset nc1  
unset nc2
unset nc3
unset nc4

}


#contacts for H2B
for { set r1 1 } { $r1<=23 } { incr r1 } {
set sel1 [atomselect top "nucleic and noh"  frame $i]
set sel2 [atomselect top "$DNA_minor and (not backbone) and noh"  frame $i]
set sel3 [atomselect top "$DNA_major and (not backbone) and noh"  frame $i]

set sel4 [atomselect top "(segname CHD and resid '$r1') and noh" frame $i]
set sel5 [atomselect top "(segname CHH and resid '$r1') and noh" frame $i]


set contacts1 [lindex [measure contacts  4.0 $sel2 $sel4] 1]
set contacts2 [lindex [measure contacts  4.0 $sel2 $sel5] 1]

set contacts3 [lindex [measure contacts  4.0 $sel3 $sel4] 1]
set contacts4 [lindex [measure contacts  4.0 $sel3 $sel5] 1]

set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
set nc3 [llength $contacts3]
set nc4 [llength $contacts4]

puts -nonewline $outfile41 "\t$nc1"
puts -nonewline $outfile42 "\t$nc2"
puts -nonewline $outfile43 "\t$nc3"
puts -nonewline $outfile44 "\t$nc4"

$sel1 delete
$sel2 delete
$sel3 delete
$sel4 delete
$sel5 delete

unset contacts1
unset contacts2
unset contacts3
unset contacts4
unset nc1  
unset nc2
unset nc3
unset nc4

}

puts -nonewline $outfile11 "\n"
puts -nonewline $outfile12 "\n"
puts -nonewline $outfile13 "\n"
puts -nonewline $outfile14 "\n"

puts -nonewline $outfile21 "\n"
puts -nonewline $outfile22 "\n"
puts -nonewline $outfile23 "\n"
puts -nonewline $outfile24 "\n"

puts -nonewline $outfile31 "\n"
puts -nonewline $outfile32 "\n"
puts -nonewline $outfile33 "\n"
puts -nonewline $outfile34 "\n"

puts -nonewline $outfile41 "\n"
puts -nonewline $outfile42 "\n"
puts -nonewline $outfile43 "\n"
puts -nonewline $outfile44 "\n"

}

close $outfile11 
close $outfile12
close $outfile13
close $outfile14

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

close $outfile41
close $outfile42
close $outfile43
close $outfile44

mol delete all
}
exit

