#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 12:59:31 2019

@author: pengy10
"""

import pandas as pd
import csv
import collections
import operator

# Protein-DNA contacts from PDB structures 
BP_contacts = pd.read_csv("BP_DNA_contact.dat",sep='\t' ) 
# Tail-DNA contacts from simulations
tail_contacts = pd.read_csv("Tail_DNA_mean_contacts_per_basepair_2_fold_combine.csv",
                       sep=',', header =None,skiprows = 1)
output = dict()
overlap_ratio = dict()
for i in range(0,len(BP_contacts)): # loop over each structure
    count_all_interface = 0 
    count_overlapped_interface = 0
    
    for j in range(1,188): ## loop over each basepairs
        if(BP_contacts.iloc[i,j] !=0 ): ## DNA basepair interact with binding partners in structures
            count_all_interface = count_all_interface +1 ## find the entire protein binding interface on DNA 
            ##  DNA basepair also have mean contacts number with tails > 5 in simulations
            if tail_contacts.iloc[abs(j-94),0] >=5: 
                count_overlapped_interface = count_overlapped_interface +1 ## find the overlapped portion with tails 
    overlap_ratio[BP_contacts.iloc[i,0]] = count_overlapped_interface/count_all_interface ## calculate the overlap ratio
    

sorted_overlap_ratio = sorted(overlap_ratio.items(), key=operator.itemgetter(1), reverse=True)
sorted_overlap_ratio = collections.OrderedDict(sorted_overlap_ratio)

with open('BP_tails_interface_overlap_ratio.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in sorted_overlap_ratio.items():
       writer.writerow([key[0:4], value])
       
