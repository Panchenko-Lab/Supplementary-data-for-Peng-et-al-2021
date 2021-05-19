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

set sel [atomselect top "resid 375 to 502 "]
$sel set chain G
$sel set segname CHG
renumber $sel 1
$sel delete

set sel [atomselect top "resid 503 to 604 "]
$sel set chain B
$sel set segname CHB
renumber $sel 1
$sel delete


set sel [atomselect top "resid 605 to 739 "]
$sel set chain E
$sel set segname CHE
renumber $sel 1
$sel delete

set sel [atomselect top "resid 740 to 861 "]
$sel set chain H
$sel set segname CHH
renumber $sel 1
$sel delete

set sel [atomselect top "resid 862 to 983 "]
$sel set chain D
$sel set segname CHD
renumber $sel 1
$sel delete


set sel [atomselect top "resid 984 to 1111 "]
$sel set chain C
$sel set segname CHC
renumber $sel 1
$sel delete


set sel [atomselect top "resid 1112 to 1246 "]
$sel set chain A
$sel set segname CHA
renumber $sel 1
$sel delete


set sel [atomselect top "resid 1247 to 1348 "]
$sel set chain F
$sel set segname CHF
renumber $sel 1
$sel delete

