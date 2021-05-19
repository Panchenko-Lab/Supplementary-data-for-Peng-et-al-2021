##Measure the sasa and contacts for the nusleosomal DNA in current existing PDB structures
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

set PDB [lrange [ls ./PDB_files/structres_contacts_with_DNA/new/] 1 end]

#set outfile1 [open DNA_sasa.dat w]
set outfile2 [open BP_DNA_contact.dat a]

#puts -nonewline $outfile1 "PDB_ID"
puts -nonewline $outfile2 "PDB_ID"

for { set r1 -93 } { $r1<=93 } { incr r1 } {
#puts -nonewline $outfile1 "\t$r1"
puts -nonewline $outfile2 "\t$r1"
}
#puts -nonewline $outfile1 "\n"
puts -nonewline $outfile2 "\n"

foreach input $PDB {
 
mol load pdb ./PDB_files/structres_contacts_with_DNA/new/$input

echo "$input"
#puts -nonewline $outfile1 "$input"
puts -nonewline $outfile2 "$input"

## renumer the residue number in DNA; some structure start from 1 onstead of -73
#set sel1 [atomselect top "chain I and nucleic"]
#set sel2 [atomselect top "chain J and nucleic"]
#set sel3 [atomselect top "(chain I and nucleic) and name P"]
#set len [llength [$sel3 get resid]]
#set start [lindex [$sel3 get resid] 0]
#set end [lindex [$sel3 get resid] end]

#if { ($start >= 1) && ($len <=156) && ($len >=130) } {
#set num_start2 [expr ($len -1)/2]
#set num_start1 [expr $num_start2 * (-1)]
#renumber $sel1 $num_start1
#renumber $sel2 $num_start2
#}

for { set r1 -93 } { $r1<=93 } { incr r1 } {
set r2 [expr $r1 * (-1)]
set sel1 [atomselect top "((chain I and resid '$r1') or (chain J and resid '$r2')) and noh and nucleic"]
set sel2 [atomselect top "protein and (not (chain A B C D E F G H )) and noh"]
set sel3 [atomselect top "(protein or (chain I J and nucleic)) and noh "]
set sel4 [atomselect top "(chain A B C D E F G H I J) and noh and (protein or nucleic)"]

## measure  sasa
#set sasa1 [measure sasa 1.4  $sel3 -restrict $sel1]
#set sasa2 [measure sasa 1.4  $sel4 -restrict $sel1]

#set dsasa [expr $sasa2 - $sasa1]

#puts -nonewline $outfile1 "\t$dsasa"

## measure contacts 
set contacts [lindex [measure contacts  4.0 $sel2 $sel1] 1]
set nc [llength $contacts]
puts -nonewline $outfile2 "\t$nc"


}
#puts -nonewline $outfile1 "\n"
puts -nonewline $outfile2 "\n"

$sel1 delete
$sel2 delete
$sel3 delete
$sel4 delete
unset contacts
unset nc
mol delete all
}
#close $outfile1
close $outfile2
#exit
