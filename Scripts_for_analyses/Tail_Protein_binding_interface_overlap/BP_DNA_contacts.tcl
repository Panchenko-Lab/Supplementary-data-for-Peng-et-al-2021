##Measure nucleosome binding protein-DNA interactions in PDB structures structures
set PDB [lrange [ls ./PDB_files/] 1 end]
set outfile [open BP_DNA_contact.dat a]

puts -nonewline $outfile2 "PDB_ID"

for { set r1 -93 } { $r1<=93 } { incr r1 } {
puts -nonewline $outfile "\t$r1"
}
puts -nonewline $outfile "\n"

foreach $input_PDB $PDB {
mol load pdb ./PDB_files/$input_PDB
puts -nonewline $outfile "$input_PDB"

for { set r1 -93 } { $r1<=93 } { incr r1 } {
set r2 [expr $r1 * (-1)]
set sel1 [atomselect top "((chain I and resid '$r1') or (chain J and resid '$r2')) and noh and nucleic"]
set sel2 [atomselect top "protein and (not (chain A B C D E F G H )) and noh"]

## measure contacts 
set contacts [lindex [measure contacts  4.0 $sel2 $sel1] 1]
set nc [llength $contacts]
puts -nonewline $outfile2 "\t$nc"
}
puts -nonewline $outfile "\n"

$sel1 delete
$sel2 delete
unset contacts
unset nc
mol delete all
}
close $outfile
