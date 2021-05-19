AMBER must be installed and the environment variable $AMBERHOME must be set
to the directory AMBER is installed locally. Please see the AMBER manual
regarding how to do this. The zip file ffptm.zip contains all the
files related to the forcefield (ffptm.in and all .frcmod files as described
below).

Before running any simulations, a user must prepare their protein structure
for simulation in AMBER. This is typically done in tleap. The following is the
contents of the file tleapmin.in which contains the syntax needed to
import the parameters to perform calculations using FFPTM's parametrization of
protonated Phosphoserine (SEP) on a protein structure containing one or
more SEP residues. /DirectoryContainingFFPTM/ refers to the path to the
directory containing FFPTM (where ffptm.zip is unzipped). In this
particular example it contains ffptm.in and SEP.frcmod.


File: tleapmin.in
_______________________________________________________
source leaprc.ff03.r1 # imports ff03 parameters by Duan et al.

# loads FFPTM
loadAmberPrep /DirectoryContainingFFPTM/ffptm.in
# loads any residue-specific parameters needed for SEP
frcmod = loadAmberParams /DirectoryContainingFFPTM/SEP.frcmod

protein = loadpdb test.pdb
check protein

saveamberparm protein test.prmtop test.inpcrd
quit

This file produces test.prmtop and test.inpcrd which are needed in sander or
pmemd to perform AMBER simulations.

To call tleapmin.in you must use the command:

COMMAND: Serial call to tleap
_______________________________________________________
$AMBERHOME/bin/tleap -f tleapmin.in


Note, if your structure has multiple post-translational modifications, you may
add more lines to tleapmin.in as per the following command example.

COMMAND: Example of how to load in multiple frcmod files to tleapmin.in
_______________________________________________________
$AMBERHOME/bin/tleap -f tleapmin.in
frcmoda = loadAmberParams /DirectoryContainingFFPTM/SEP.frcmod
frcmodb = loadAmberParams /DirectoryContainingFFPTM/M3L.frcmod


Using the newly created prmtop and inpcrd files generated using
tleapmin.in, you can now perform a local minimization of your protein
containing one or more SEP residues using the following file min1.in.

File: heat1.in
_______________________________________________________
Stage 1 heating of protein 0 to 50K
 &cntrl
  imin=0, irest=0, ntx=1,
  nstlim=100000, dt=0.0005,
  ntc=2, ntf=2,
  ntt=3, gamma_ln = 5,
  tempi=0.0, temp0=50.0,
  ntpr=50, ntwx=500,
  ntb=0, igb=5,ig=-1,gbsa=1,saltcon=0.1,
  cut=16.,rgbmax=16.
 /


The protein can be subsequently heated by adapting the file above and creating
the temperature intervals 50-100, 100-150, 150-200, 200-250, 250-300. Then,
the protein can be submitted for production using the input file below.

File: prod.in
_______________________________________________________
Stage 2 equilibration 0-5.0ns
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=2500000, dt=0.002,
  ntc=2, ntf=2,
  ntt=3, gamma_ln = 5,
  tempi=300.0, temp0=300.0,
  ntpr=500, ntwx=500,
  ntb=0, igb=5,ig=-1,gbsa=1,saltcon=0.1,
  cut=16.,rgbmax=16.
 /


The contents of ffptm.in are the following consistent with AMBER conventions

Dimethylarginine (Anti-symmetric) (name of the residue)
dimethylarg.res  (Name of the output file generated in from Antechamber)
DA2   INT  0     (Unique 3-letter code naming the amino acid)
CORRECT     OMIT DU   BEG
  0.0000
(1) # of atom in tree
(2) Unique atom name for atom I
(3) Symbol for atom I which defines its forcefield atom type
and is used in the module PARM for assigning the forcefield parameters
(4) The topological type (tree symbol) for atom I
(5) The atom number to which atom I is connected
(6) The atom number to which atom I makes an angle along with (5)
(7) The atom number which atom I makes a dihedral along with (5) and (6)
(8) The bond length between atoms I and (5)
(9)The bond angle between the bond angle between (6),(5) and (1)
(10) The dihedral angle between (7), (6), (5), and (1)
(11) The partial atomic charge on atom I generated from the RESP-fitting
procedure.

 (1)  (2)  (3)   (4)  (5) (6) (7)     (8)       (9)       (10)      (11)
   1  DUMM  DU    M    0  -1  -2     0.000      .0        .0      .00000
   2  DUMM  DU    M    1   0  -1     1.449      .0        .0      .00000
   3  DUMM  DU    M    2   1   0     1.522   111.1        .0      .00000
   4  N     N     M    3   2   1     1.540   111.208   180.000 -0.376741
   5  H     H     E    4   3   2     0.994   107.922   -18.167  0.270465
   6  CA    CT    M    4   3   2     1.448    10.429   150.011 -0.119858
   7  CB    CT    3    6   4   3     1.539   110.455     4.603 -0.129377
   8  CG    CT    3    7   6   4     1.540   113.981  -168.737 -0.016346
   9  CD    CT    3    8   7   6     1.534   114.694   -91.332 -0.092919
  10  NE    N2    B    9   8   7     1.465   112.221   171.176 -0.288262

