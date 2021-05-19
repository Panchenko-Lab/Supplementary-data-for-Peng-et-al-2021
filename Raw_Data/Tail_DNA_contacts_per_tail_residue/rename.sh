#!/bin/bash
while read name1 name2; 
do
{
mv tail_contacts_all_${name1}_h2a_C.dat tail_DNA_mean_contacts_${name2}_h2a_C.dat
mv tail_contacts_all_${name1}_h2a_G.dat tail_DNA_mean_contacts_${name2}_h2a_G.dat

mv tail_contacts_all_${name1}_h2b_D.dat tail_DNA_mean_contacts_${name2}_h2b_D.dat
mv tail_contacts_all_${name1}_h2b_H.dat tail_DNA_mean_contacts_${name2}_h2b_H.dat

mv tail_contacts_all_${name1}_h3_A.dat tail_DNA_mean_contacts_${name2}_h3_A.dat
mv tail_contacts_all_${name1}_h3_E.dat tail_DNA_mean_contacts_${name2}_h3_E.dat

mv tail_contacts_all_${name1}_h4_B.dat tail_DNA_mean_contacts_${name2}_h4_B.dat
mv tail_contacts_all_${name1}_h4_F.dat tail_DNA_mean_contacts_${name2}_h4_F.dat
} 
done < ../name_match.txt
