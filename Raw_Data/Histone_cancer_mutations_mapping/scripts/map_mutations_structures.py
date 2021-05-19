#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 22:35:51 2020

@author: yunhuipeng
"""
import os
import numpy as np
import pandas as pd 
import requests
import subprocess
import re
import csv

his_mutations = pd.read_csv("histone_mutations_mapped.txt",  sep = ",", header=None)

binding_sites = pd.read_csv("combined_all.txt",  sep = "\t", header=None)

uniprot_structures = []
for lines in binding_sites.iloc[:,0]:
    uniprot_structures.append(lines)
uniprot_structures = list(set(uniprot_structures))

with open("mutation_all.txt", "w") as fh:
    for i in range(0,len(his_mutations)):
        if (his_mutations.iloc[i,1] in uniprot_structures):
            fh.write(his_mutations.iloc[i,1]+"\t" + his_mutations.iloc[i,2] +"\t" +\
                      his_mutations.iloc[i,3]+ "\t" + "mutation" +"\n" )
    
