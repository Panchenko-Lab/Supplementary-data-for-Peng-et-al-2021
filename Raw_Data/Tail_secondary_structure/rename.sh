#!/bin/bash
while read name1 name2; 
do
{
mv tail_SSE_${name1}_h2a_C.dat SSE_${name2}_h2a_C.dat
mv tail_SSE_${name1}_h2a_G.dat SSE_${name2}_h2a_G.dat

mv tail_SSE_${name1}_h2b_D.dat SSE_${name2}_h2b_D.dat
mv tail_SSE_${name1}_h2b_H.dat SSE_${name2}_h2b_H.dat

mv tail_SSE_${name1}_h3_A.dat SSE_${name2}_h3_A.dat
mv tail_SSE_${name1}_h3_E.dat SSE_${name2}_h3_E.dat

mv tail_SSE_${name1}_h4_B.dat SSE_${name2}_h4_B.dat
mv tail_SSE_${name1}_h4_F.dat SSE_${name2}_h4_F.dat
} 
done < ../name_match.txt
