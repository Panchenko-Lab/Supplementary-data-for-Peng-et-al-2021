proc renumber { sel start } {
	if { [$sel num] == 0 } {
		puts "Error in renumber: empty selection!"
		return
	}
	set oresid [ $sel get resid ]
	set delta [ expr $start - [ lindex $oresid 0] ]
	set nresid { }
	foreach r $oresid {
		lappend nresid [ expr $r + $delta ]
	}
	$sel set resid $nresid
}


set sel [atomselect top "resid 1 to 187 "]
$sel set chain I
$sel set segname CHI
renumber $sel -93
$sel delete

set sel [atomselect top "resid 188 to 374 "]
$sel set chain J
$sel set segname CHJ
renumber $sel -93
$sel delete

set sel [atomselect top "resid 375 to 509 "]
$sel set chain A
$sel set segname CHA
renumber $sel 1
$sel delete

set sel [atomselect top "resid 510 to 611 "]
$sel set chain B
$sel set segname CHB
renumber $sel 1
$sel delete


set sel [atomselect top "resid 612 to 739 "]
$sel set chain C
$sel set segname CHC
renumber $sel 1
$sel delete

set sel [atomselect top "resid 740 to 841 "]
$sel set chain F
$sel set segname CHF
renumber $sel 1
$sel delete

set sel [atomselect top "resid 842 to 969 "]
$sel set chain G
$sel set segname CHG
renumber $sel 1
$sel delete


set sel [atomselect top "resid 970 to 1091 "]
$sel set chain H
$sel set segname CHH
renumber $sel 1
$sel delete


set sel [atomselect top "resid 1092 to 1226 "]
$sel set chain E
$sel set segname CHE
renumber $sel 1
$sel delete


set sel [atomselect top "resid 1227 to 1348 "]
$sel set chain D
$sel set segname CHD
renumber $sel 1
$sel delete

