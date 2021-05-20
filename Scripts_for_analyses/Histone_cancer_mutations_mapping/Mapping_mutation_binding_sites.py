# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
import numpy as np
import pandas as pd 
import requests
import re
import csv

## Retrive histone sequence by uniprotID and perform MSA using clustal-omega
def uniprotID_to_fa_msa(uniprotID_file, fa_file):
    uniprot_ID = pd.read_table(uniprotID_file,header = None)
    with open(fa_file, "w") as fh:
        for i in range(0, len(uniprot_ID)):    
            response = requests.post("https://www.uniprot.org/uniprot/" + uniprot_ID.iloc[i][0] + ".fasta")
            fh.write(response.text) 
    os.system("./clustal-omega-1.2.3-macosx -i " + fa_file + " -o " + fa_file[0:-3] + "_msa.fa --auto -v --force")


def convert_list_to_string(input_list,sting_seq): 
      
    # Converting integer list to string list 
    s = [str(i) for i in input_list] 
      
    # Join list items using join() 
    res = str(sting_seq.join(s)) 
    return(res)


### mapping binding site in the msa and count binding partners per sites
def map_binding_sites_msa(binding_sites_file,msa_file, output_map,output_map_consensus):
#    global binding_sites_align
    binding_sites = pd.read_table(binding_sites_file) 
    binding_sites_align = pd.DataFrame()
    f = open(msa_file, 'r')
    msa = list(set(f.read().split(">")[1:])) ## msa from cluster2 have duplcates ## need check
    f.close
    ## load MSA file
    for lines in msa:
        uniprotID_msa = lines.split("\n")[0].split("|")[1].lower()
        seq = convert_list_to_string(lines.split("\n")[1:],"")
        gap_pos = []
        ## find the gaps in the alignment
        for m in re.finditer(r"-", seq):
            gap_pos.append(int(m.start()))
        ## convert the binding sites number accoring to gaps
        for i in range(0,len(binding_sites)):
            uniprot_ID_sites = binding_sites.iloc[i,0]
            binding_pos = int(binding_sites.iloc[i,2])
            bp_uniprot = binding_sites.iloc[i,4]
            if uniprot_ID_sites == uniprotID_msa:
                for pos in gap_pos:
                    if pos <= binding_pos-1:
                        binding_pos = binding_pos +1
                if binding_sites.iloc[i,3] == "PDB" :
                    binding_sites_align = binding_sites_align.append({"uniprot_ID": uniprotID_msa, "pos": binding_pos, "bp": bp_uniprot}, ignore_index=True)
                    
    with open(output_map, 'w') as myfile:
        uniprot_list = np.unique(binding_sites_align.iloc[:,2])
        for uniprotID in uniprot_list:
            site = []
            binding_sites_map = [0]*len(seq)
            binding_sites_map.insert(0,uniprotID)
            #count binding partner number per sites
            for i in range(0,len(binding_sites_align)):
                if (binding_sites_align.iloc[i,2] == uniprotID):
                   site.append(binding_sites_align.iloc[i,1])
            site_counts_pos = np.unique(site,return_counts=True)[0]
            site_counts_num = np.unique(site,return_counts=True)[1]
            for pos, num  in zip(site_counts_pos, site_counts_num):
                binding_sites_map[int(pos)] = int(num)
            wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
            wr.writerow(binding_sites_map)     
            
    with open(output_map_consensus, 'w') as myfile:
        binding_sites_unique = binding_sites_align.iloc[:,0:2].drop_duplicates()
        for pos in range(1,len(seq)+1):
            count = 0
            for j in range(0,len(binding_sites_unique)):
                if (pos == binding_sites_unique.iloc[j,1]):
                    count +=1
            myfile.write(str(pos) + "," + str(count) + '\n' )    
            
def map_mutation_sites_msa(binding_sites_file,msa_file, output_map,output_map_consensus):
#    global binding_sites_align
    binding_sites = pd.read_table(binding_sites_file) 
    binding_sites_align = pd.DataFrame()
    f = open(msa_file, 'r')
    msa = list(set(f.read().split(">")[1:])) ## msa from cluster2 have duplcates ## need check
    f.close
    ## load MSA file
    for lines in msa:
        uniprotID_msa = lines.split("\n")[0].split("|")[1].lower()
        seq = convert_list_to_string(lines.split("\n")[1:],"")
        gap_pos = []
        ## find the gaps in the alignment
        for m in re.finditer(r"-", seq):
            gap_pos.append(int(m.start()))
        ## convert the binding sites number accoring to gaps
        for i in range(0,len(binding_sites)):
            uniprot_ID_sites = binding_sites.iloc[i,0]
            binding_pos = int(binding_sites.iloc[i,2])
            if uniprot_ID_sites == uniprotID_msa:
                for pos in gap_pos:
                    if pos <= binding_pos-1:
                        binding_pos = binding_pos +1
                if binding_sites.iloc[i,3] == "mutation": ## mutation number should add 1 to match uniprot number 
                    binding_sites_align = binding_sites_align.append({"uniprot_ID": uniprotID_msa, "pos": binding_pos+1}, ignore_index=True)
   
    with open(output_map, 'w') as myfile:
        uniprot_list = np.unique(binding_sites_align.iloc[:,1])
        for uniprotID in uniprot_list:
            site = []
            binding_sites_map = [0]*len(seq)
            binding_sites_map.insert(0,uniprotID)
            #count mutation numbers per sites
            for i in range(0,len(binding_sites_align)):
                if (binding_sites_align.iloc[i,1] == uniprotID):
                   site.append(binding_sites_align.iloc[i,0])
            site_counts_pos = np.unique(site,return_counts=True)[0]
            site_counts_num = np.unique(site,return_counts=True)[1]
            for pos, num  in zip(site_counts_pos, site_counts_num):
                binding_sites_map[int(pos)] = int(num)
            wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
            wr.writerow(binding_sites_map)   
            
    with open(output_map_consensus, 'w') as myfile:
        binding_sites_unique = binding_sites_align.iloc[:,0:2]
        for pos in range(1,len(seq)+1):
            count = 0
            for j in range(0,len(binding_sites_unique)):
                if (pos == binding_sites_unique.iloc[j,0]):
                    count +=1
            myfile.write(str(pos) + "," + str(count) + '\n' )    
            
  
## Perform mutiple sequence alignment for each histone family
uniprotID_to_fa_msa("H2A_uniprot_IDs.txt", "H2A.fa")
uniprotID_to_fa_msa("H2B_uniprot_IDs.txt", "H2B.fa")
uniprotID_to_fa_msa("H3_uniprot_IDs.txt", "H3.fa")
uniprotID_to_fa_msa("H4_uniprot_IDs.txt", "H4.fa")

## mapping protein binding sites onto aligned histone sequences
map_binding_sites_msa("Binding_sites_H2A.txt","H2A_msa.fa", 'H2A_mapping_binding_sites.txt',"H2A_binding_sites_msa.txt")   
map_binding_sites_msa("Binding_sites_H2B.txt","H2B_msa.fa", 'H2B_mapping_binding_sites.txt',"H2B_binding_sites_msa.txt")   
map_binding_sites_msa("Binding_sites_H3.txt","H3_msa.fa", 'H3_mapping_binding_sites.txt',"H3_binding_sites_msa.txt")   
map_binding_sites_msa("Binding_sites_H4.txt","H4_msa.fa", 'H4_mapping_binding_sites.txt',"H4_binding_sites_msa.txt")   

           
## mapping pcancer mutation sites onto aligned histone sequences
map_mutation_sites_msa("Mutation_sites_H2A.txt","H2A_msa.fa", 'H2A_mapping_mutation_sites.txt',"H2A_mutation_sites_msa.txt")   
map_mutation_sites_msa("Mutation_sites_H2B.txt","H2B_msa.fa", 'H2B_mapping_mutation_sites.txt',"H2B_mutation_sites_msa.txt")   
map_mutation_sites_msa("Mutation_sites_H3.txt","H3_msa.fa", 'H3_mapping_mutation_sites.txt',"H3_mutation_sites_msa.txt")   
map_mutation_sites_msa("Mutation_sites_H4.txt","H4_msa.fa", 'H4_mapping_mutation_sites.txt',"H4_mutation_sites_msa.txt")            
