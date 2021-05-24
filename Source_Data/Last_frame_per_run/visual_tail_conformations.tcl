#load PDB 
mol load pdb ext_sym.pdb
foreach name {1aoi 1eqz ext2 sym} {
mol load  pdb ./${name}_run1.pdb
mol load  pdb ./${name}_run2.pdb
mol load  pdb ./${name}_run3.pdb
mol load  pdb ./${name}_run4.pdb
mol load  pdb ./${name}_run5.pdb
}

mol load pdb sym_gromacs_run1.pdb
mol load pdb sym_gromacs_run2.pdb


display rendermode GLSL
for {set i 0} {$i <= 22} {incr i 1} {
mol delrep 0 $i
} 

package require colorscalebar
axes location off
source input_param.tcl
source add_text_layer.tcl
color Display Background white

#  H3 tail representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color ColorID 0 
mol selection {chain A E and resid 1 to 36}
mol material AOShiny
for {set i 1} {$i <= 22} {incr i 1} {
mol addrep $i
}
mol selupdate 0 top 0
mol colupdate 0 top 0
mol scaleminmax top  0 $minsc $maxsc

# H4 tail representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color ColorID 7
mol selection {chain B F and resid 1 to 20}
mol material AOShiny
for {set i 1} {$i <= 22} {incr i 1} {
mol addrep $i
}
mol selupdate 1 top 0
mol colupdate 1 top 0
mol scaleminmax top  1 $minsc $maxsc

#  H2a tail representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color ColorID 4
mol selection {chain C G and resid 1 to 13 119 to 128}
mol material AOShiny
for {set i 1} {$i <= 22} {incr i 1} {
mol addrep $i
}
mol selupdate 2 top 0
mol colupdate 2 top 0
mol scaleminmax top  2 $minsc $maxsc

# H2b tail representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color ColorID 1
mol selection {chain D H and resid 1 to 23}
mol material AOShiny
for {set i 1} {$i <= 22} {incr i 1} {
mol addrep $i
}
mol selupdate 3 top 0
mol colupdate 3 top 0
mol scaleminmax top  3 $minsc $maxsc


# Dna base  representation
mol representation NewCartoon 1.090000 10.000000 2.070000 1
mol representation Surf
mol color ColorID 6
mol selection {(nucleic not backbone) and noh }
mol material AOEdgy
for {set i 0} {$i <= 0} {incr i 1} {
mol addrep $i
} 
mol selupdate 4 top 0
mol colupdate 4 top 0
mol scaleminmax top  4 $minsc $maxsc


# Dna backbone representation in + strand
mol representation NewCartoon 1.090000 10.000000 2.070000 1
mol representation Surf
mol color ColorID 10
mol selection {chain I and (backbone and nucleic)} 
##name P O1P O2P O3' O5' C5'
mol material AOEdgy
for {set i 0} {$i <= 0} {incr i 1} {
mol addrep $i
}
mol selupdate 5 top 0
mol colupdate 5 top 0
mol scaleminmax top  5 $minsc $maxsc

# Dna backbone representation in - strand
mol representation NewCartoon 1.090000 10.000000 2.070000 1
mol representation Surf
mol color ColorID 13
mol selection {chain J and (backbone and nucleic)}
mol material AOEdgy
for {set i 0} {$i <= 0} {incr i 1} {
mol addrep $i
}
mol selupdate 6 top 0
mol colupdate 6 top 0
mol scaleminmax top 6 $minsc $maxsc


color chain A blue
color chain B green
color chain C yellow
color chain D red
color chain E blue
color chain F green
color chain G yellow
color chain H red
color Highlight Nonback 6
color Highlight Nucback 2
color Display Background white
