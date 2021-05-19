#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 12:59:31 2019

@author: pengy10
"""

import pandas as pd
import numpy as np
import csv
import collections
import operator
data = pd.read_csv("BP_DNA_contact.dat",sep='\t' )
contacts = pd.read_csv("/Users/yunhuipeng/Desktop/Interactome_paper/tail_paper/figure2/DNA_contacts_4500ns/contact_mean_half.csv",
                       sep=',', header =None,skiprows = 1)
output = dict()
overlap_ratio = dict()
for i in range(0,len(data)):
    output[data.iloc[i,0]] = list()
    count_all = 0
    count_overlap = 0
    for j in range(1,188):
        if(data.iloc[i,j] !=0 ):
            binding = str(j-94) + "_" + str(data.iloc[i,j])
            output[data.iloc[i,0]].append(binding)
            count_all = count_all +1 
            if contacts.iloc[abs(j-94),0] >=5:
                count_overlap = count_overlap +1
    overlap_ratio[data.iloc[i,0]] = count_overlap/count_all
    
f = open('binding_hotspot.csv','w')          
for key in output.keys():
    resid = ','.join(output[key])
    line = str(key) + ":," + resid + "\n"
    f.write(line)
f.close()

sorted_overlap_ratio = sorted(overlap_ratio.items(), key=operator.itemgetter(1), reverse=True)
sorted_overlap_ratio = collections.OrderedDict(sorted_overlap_ratio)

with open('overlap_ratio_5.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in sorted_overlap_ratio.items():
       writer.writerow([key[0:4], value])
       
