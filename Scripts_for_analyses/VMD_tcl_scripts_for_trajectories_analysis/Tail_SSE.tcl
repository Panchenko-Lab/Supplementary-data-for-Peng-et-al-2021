foreach name  {ModelA ModelB ModelC ModelD} {
mol delete all
mol load parm7 ../${name}_rw.prmtop 
mol addfile ../${name}_rw_run1_4500ns.nc first 10000 last 225000 step 10 waitfor all
mol addfile ../${name}_rw_run2_800ns.nc first 2500 last 40000 step 10 waitfor all
mol addfile ../${name}_rw_run3_800ns.nc first 2500 last 40000 step 10 waitfor all
mol addfile ../${name}_rw_run4_800ns.nc first 2500 last 40000 step 10 waitfor all
mol addfile ../${name}_rw_run5_800ns.nc first 2500 last 40000 step 10 waitfor all

source rename_chainID_${name}.tcl
set nframes [expr  [molinfo top get numframes] - 1 ]

set outfile1 [open tail_SSE_${name}_h3_A.dat w]
set outfile2 [open tail_SSE_${name}_h3_E.dat w]
set outfile3 [open tail_SSE_${name}_h4_B.dat w]
set outfile4 [open tail_SSE_${name}_h4_F.dat w]
set outfile5 [open tail_SSE_${name}_h2a_C.dat w]
set outfile6 [open tail_SSE_${name}_h2a_G.dat w]
set outfile7 [open tail_SSE_${name}_h2b_D.dat w]
set outfile8 [open tail_SSE_${name}_h2b_H.dat w]
	

puts -nonewline $outfile1 "Time"
puts -nonewline $outfile2 "Time"
puts -nonewline $outfile3 "Time"
puts -nonewline $outfile4 "Time"
puts -nonewline $outfile5 "Time"
puts -nonewline $outfile6 "Time"
puts -nonewline $outfile7 "Time"
puts -nonewline $outfile8 "Time"


for { set r1 1 } { $r1<=36 } { incr r1 } {
puts -nonewline $outfile1 "\t$r1"
puts -nonewline $outfile2 "\t$r1"
}
for { set r1 1 } { $r1<=20 } { incr r1 } {
puts -nonewline $outfile3 "\t$r1"
puts -nonewline $outfile4 "\t$r1"
}
for { set r1 1 } { $r1<=13 } { incr r1 } {
puts -nonewline $outfile5 "\t$r1"
puts -nonewline $outfile6 "\t$r1"
}
for { set r1 119 } { $r1<=128 } { incr r1 } {
puts -nonewline $outfile5 "\t$r1"
puts -nonewline $outfile6 "\t$r1"
}
for { set r1 1 } { $r1<=23 } { incr r1 } {
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
animate goto $i 
display update ui 
mol ssrecalc top 

set time [expr 0.1 * $i]
puts -nonewline $outfile1 [format "%.3f" "$time"]
puts -nonewline $outfile2 [format "%.3f" "$time"]
puts -nonewline $outfile3 [format "%.3f" "$time"]
puts -nonewline $outfile4 [format "%.3f" "$time"]
puts -nonewline $outfile5 [format "%.3f" "$time"]
puts -nonewline $outfile6 [format "%.3f" "$time"]
puts -nonewline $outfile7 [format "%.3f" "$time"]
puts -nonewline $outfile8 [format "%.3f" "$time"]
puts  [format "Time %.3f" "$time"]
unset time
#sse for H3
for { set r1 1 } { $r1<=36 } { incr r1 } {
set sel1 [atomselect top "(segname CHA and resid '$r1') and name CA"]
set sel2 [atomselect top "(segname CHE and resid '$r1') and name CA"]
$sel1 frame $i 
$sel1 update 
$sel2 frame $i  
$sel2 update  
set sse1 [$sel1 get structure]
set sse2 [$sel2 get structure]
puts -nonewline $outfile1 "\t$sse1"
puts -nonewline $outfile2 "\t$sse2"
$sel1 delete
$sel2 delete
unset sse1 
unset sse2 
}
#sse for H4
for { set r1 1 } { $r1<=20 } { incr r1 } {
set sel1 [atomselect top "(segname CHB and resid '$r1') and name CA" ]
set sel2 [atomselect top "(segname CHF and resid '$r1') and name CA" ]
$sel1 frame $i  
$sel1 update  
$sel2 frame $i
$sel2 update
set sse1 [$sel1 get structure]
set sse2 [$sel2 get structure]
puts -nonewline $outfile3 "\t$sse1"
puts -nonewline $outfile4 "\t$sse2"
$sel1 delete
$sel2 delete
unset sse1 
unset sse2 
}
#sse for H2A
for { set r1 1 } { $r1<=13 } { incr r1 } {
set sel1 [atomselect top "(segname CHC and resid '$r1') and name CA" ]
set sel2 [atomselect top "(segname CHG and resid '$r1') and name CA" ]
$sel1 frame $i  
$sel1 update  
$sel2 frame $i
$sel2 update
set sse1 [$sel1 get structure]
set sse2 [$sel2 get structure]
puts -nonewline $outfile5 "\t$sse1"
puts -nonewline $outfile6 "\t$sse2"
$sel1 delete
$sel2 delete
unset sse1 
unset sse2 
}

for { set r1 119 } { $r1<=128 } { incr r1 } {
set sel1 [atomselect top "(segname CHC and resid '$r1') and name CA" ]
set sel2 [atomselect top "(segname CHG and resid '$r1') and name CA" ]
$sel1 frame $i  
$sel1 update  
$sel2 frame $i
$sel2 update
set sse1 [$sel1 get structure]
set sse2 [$sel2 get structure]
puts -nonewline $outfile5 "\t$sse1"
puts -nonewline $outfile6 "\t$sse2"
$sel1 delete
$sel2 delete
unset sse1 
unset sse2 
}

#contacts for H2B
for { set r1 1 } { $r1<=23 } { incr r1 } {
set sel1 [atomselect top "(segname CHD and resid '$r1') and name CA"]
set sel2 [atomselect top "(segname CHH and resid '$r1') and name CA"]
$sel1 frame $i  
$sel1 update  
$sel2 frame $i
$sel2 update
set sse1 [$sel1 get structure]
set sse2 [$sel2 get structure]
puts -nonewline $outfile7 "\t$sse1"
puts -nonewline $outfile8 "\t$sse2"
$sel1 delete
$sel2 delete
unset sse1 
unset sse2 
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


