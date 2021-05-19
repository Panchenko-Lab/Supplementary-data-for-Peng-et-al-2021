#!/bin/bash
while read name1 name2; 
do
{
mv dsasa_H2A_${name1}.dat dsasa_H2A_${name2}.dat
mv dsasa_H2B_${name1}.dat dsasa_H2B_${name2}.dat
mv dsasa_H3_${name1}.dat dsasa_H3_${name2}.dat
mv dsasa_H4_${name1}.dat dsasa_H4_${name2}.dat
mv dsasa_all_${name1}.dat dsasa_all_${name2}.dat
mv rsasa_all_${name1}.dat rsasa_all_${name2}.dat
} 
done < ../name_match.txt
