#!/bin/bash
while read name1 name2; 
do
{
mv MMPBSA_ddG_${name1}_run1.csv MMPBSA_ddG_${name2}_run1.csv
mv MMPBSA_ddG_${name1}_run2.csv MMPBSA_ddG_${name2}_run2.csv
mv MMPBSA_ddG_${name1}_run3.csv MMPBSA_ddG_${name2}_run3.csv
mv MMPBSA_ddG_${name1}_run4.csv MMPBSA_ddG_${name2}_run4.csv
mv MMPBSA_ddG_${name1}_run5.csv MMPBSA_ddG_${name2}_run5.csv

} 
done < ./name_match.txt
