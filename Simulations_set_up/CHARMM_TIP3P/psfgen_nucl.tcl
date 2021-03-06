# This file is a VMD/psfgen script to generate topology
# for nucleosome in CHARMM 36 ff
# 
# Will generate a solvated nucleosome with 0.15M ions,
# distance to box 20 A
#

mol load pdb Model_A.pdb

set protchains "A B C D E F G H"
set nuclchains "I J"

#######################


#set water [atomselect top water]
#important
#$water set resid [$water get residue]  
#$water writepdb t_1kx5_water.pdb

set protein [atomselect top protein]
set chains [lsort -unique [$protein get pfrag]]
#sorting is not right!!!
foreach chain $protchains {
  set sel [atomselect top "chain $chain and protein"]
  $sel writepdb t_1kx5_p_${chain}.pdb
}

foreach chain $nuclchains {
  set sel [atomselect top "chain $chain and nucleic"]
  $sel writepdb t_1kx5_n_${chain}.pdb
}

package require psfgen

# the topology files are from MacKerlls package
# These are however in different format and with different patches for DNA

# protein topology is read ok by psfgen
topology top_all36_prot.rtf

# new nucleic acid topology: charmm has problems with autogenerate!!!
# namely it does not understand the PATCH thing, we removed it, hope this is right,
# we well then regenerate everything manually
topology top_all36_na.rtf

# the original ff includes a file for waters in a streamline format
# we had to manually decompose it to par_water_ions_36.prm and top_water_ions_36.rtf
# also deleted the H1-H2 bond
# see http://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l/2666.html
# http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-html/node22.html see end of page
# since it contains NBFIXes for sodium carboxilate interactions, we need to read it last
# in NAMD!!!
topology top_water_ions_36.rtf

pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD  

pdbalias residue DA ADE
pdbalias residue DT THY
pdbalias residue DC CYT
pdbalias residue DG GUA

pdbalias atom DA OP1 O1P
pdbalias atom DA OP2 O2P
pdbalias atom DT OP1 O1P
pdbalias atom DT OP2 O2P
pdbalias atom DC OP1 O1P
pdbalias atom DC OP2 O2P
pdbalias atom DG OP1 O1P
pdbalias atom DG OP2 O2P

pdbalias atom DT C7 C5M

# consult this page also http://www.ks.uiuc.edu/~villa/dna/

# now go with protein segments

 segment chA {
   pdb t_1kx5_p_A.pdb
   # first ACE
  # last CT3
 }
 coordpdb t_1kx5_p_A.pdb chA

segment chB {
  pdb t_1kx5_p_B.pdb
  # first ACE
  # last CT3
 }
 coordpdb t_1kx5_p_B.pdb chB

 segment chC {
  pdb t_1kx5_p_C.pdb
  # first ACE
  # last CT3
 }
 coordpdb t_1kx5_p_C.pdb chC

 segment chD {
  pdb t_1kx5_p_D.pdb
  # first ACE
  # last CT3
 }
 coordpdb t_1kx5_p_D.pdb chD

 segment chE {
  pdb t_1kx5_p_E.pdb
  # first ACE
  # last CT3
 }
 coordpdb t_1kx5_p_E.pdb chE

 segment chF {
  pdb t_1kx5_p_F.pdb
  # first ACE
  # last CT3
 }
 coordpdb t_1kx5_p_F.pdb chF

 segment chG {
  pdb t_1kx5_p_G.pdb
  # first ACE
  # last CT3
 }
 coordpdb t_1kx5_p_G.pdb chG

 segment chH {
  pdb t_1kx5_p_H.pdb
  # first ACE
  # last CT3
 }
 coordpdb t_1kx5_p_H.pdb chH


#now let's do the hard part with nucleic
# see http://www.ks.uiuc.edu/~villa/dna/
# and an example namd biotechnology nanopore tutorial. see make-psf.tcl

# this we will need to loop through resid


 segment chI {
  pdb t_1kx5_n_I.pdb
  first 5TER
  last 3TER
 }
patch DEO5 chI:-93
for {set i -92} {$i <= 93} {incr i 1}  {
patch DEOX chI:$i
}
 coordpdb t_1kx5_n_I.pdb chI

 segment chJ {
  pdb t_1kx5_n_J.pdb
  first 5TER
  last 3TER
 }
patch DEO5 chJ:-93
for {set i -92} {$i <= 93} {incr i 1}  {
patch DEOX chJ:$i
}
 coordpdb t_1kx5_n_J.pdb chJ

#this is important
regenerate angles dihedrals

# # (7) Build water segment
# pdbalias residue HOH TIP3   ; # formerly "alias residue ..."
# segment SOLV {
#  auto none
#  pdb t_1kx5_water.pdb
# }

# # (8) Read water coordinaes from PDB file
# pdbalias atom HOH O OH2     ; # formerly "alias atom ..."
# coordpdb t_1kx5_water.pdb SOLV

# # (9) Guess missing coordinates
 guesscoord

# # (10) Write structure and coordinate files
writepsf t_1kx5_namd.psf
writepdb t_1kx5_namd.pdb


package require solvate	 
#solvate t_1kx5_namd.psf t_1kx5_namd.pdb -minmax {{-125 -125 -125} {125 125 125}}  -o t_1kx5_namd_wb 

solvate t_1kx5_namd.psf t_1kx5_namd.pdb -t 20  -o t_1kx5_namd_wb 
## measure center of the water box
set file [open "centerpro.dat" w]
set sell [atomselect top all]
puts $file [measure center $sell]
close $file

resetpsf

package require autoionize

autoionize -psf t_1kx5_namd_wb.psf -pdb t_1kx5_namd_wb.pdb -sc 0.150 -o 1kx5_ready

mol load pdb 1kx5_ready.pdb
set everyone [atomselect top all]	 
set box [measure minmax $everyone]
echo "DX:"
expr [lindex [lindex $box 1] 0] - [lindex [lindex $box 0] 0]
echo "DY:"
expr [lindex [lindex $box 1] 1] - [lindex [lindex $box 0] 1]
echo "DZ:"
expr [lindex [lindex $box 1] 2] - [lindex [lindex $box 0] 2]
echo "Center:"
measure center $everyone 

#let's mark fixed atoms with 1 beta factor
set all [atomselect top all]
$all set beta 0

set fixed [atomselect top "segname CHA CHB CHC CHD CHE CHF CHF CHH CHI CHJ"]
$fixed set beta 1
$all writepdb 1kx5_ready.pdb

#Now generate file with constraint constants
$all set beta 0
set cons [atomselect top "name CA P"]
$cons set beta 90
$all writepdb 1kx5_cons.pdb

exit

