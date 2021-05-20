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


set sel [atomselect top "resid 1 to 128 "]
$sel set chain G
$sel set segname CHG
renumber $sel 1
$sel delete

set sel [atomselect top "resid 129 to 256 "]
$sel set chain C
$sel set segname CHC
renumber $sel 1
$sel delete

set sel [atomselect top "resid 257 to 391 "]
$sel set chain A
$sel set segname CHA
renumber $sel 1
$sel delete

set sel [atomselect top "resid 392 to 513 "]
$sel set chain D
$sel set segname CHD
renumber $sel 1
$sel delete


set sel [atomselect top "resid 514 to 615 "]
$sel set chain B
$sel set segname CHB
renumber $sel 1
$sel delete

set sel [atomselect top "resid 616 to 737 "]
$sel set chain H
$sel set segname CHH
renumber $sel 1
$sel delete

set sel [atomselect top "resid 738 to 872 "]
$sel set chain E
$sel set segname CHE
renumber $sel 1
$sel delete


set sel [atomselect top "resid 873 to 974 "]
$sel set chain F
$sel set segname CHF
renumber $sel 1
$sel delete


set sel [atomselect top "resid 975 to 1161 "]
$sel set chain I
$sel set segname CHI
renumber $sel -93
$sel delete


set sel [atomselect top "resid 1162 to 1348 "]
$sel set chain J
$sel set segname CHJ
renumber $sel -93
$sel delete

