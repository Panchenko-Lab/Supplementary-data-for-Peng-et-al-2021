#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 14:48:45 2019

@author: pengy10
"""

import pandas as pd 

## Different sets of simulation runs
Model_list = ["ModelA", "ModelB", "ModelC", "ModelD", "ModelD_gromacs_run1", "ModelD_gromacs_run2"]
for file_name in Model_list:
    
## unbound state is defined if the percentage of tail residues maintaining contacts with the DNA molecule is no more than the “cut-off”
    for cut_off_ratio in [0,0.1]:
        
        tail_residence_time = { i : [] for i in ["H2A_N", "H2A_C", "H2B", "H3", "H4"] }
        
        
        ## read the measurements of tail-DNA contacts from simulations  
        H2A_contacts1 = pd.read_csv("tail_DNA_mean_contacts_"+ file_name+ "_h2a_C.dat",  sep = "\t", header=0)
        H2A_contacts2 = pd.read_csv("tail_DNA_mean_contacts_"+ file_name+ "_h2a_G.dat",  sep = "\t", header=0)
        H2B_contacts1 = pd.read_csv("tail_DNA_mean_contacts_"+ file_name+ "_h2b_D.dat",  sep = "\t", header=0)
        H2B_contacts2 = pd.read_csv("tail_DNA_mean_contacts_"+ file_name+ "_h2b_H.dat",  sep = "\t", header=0)
        
        H3_contacts1 = pd.read_csv("tail_DNA_mean_contacts_"+ file_name+ "_h3_A.dat",  sep = "\t", header=0)
        H3_contacts2 = pd.read_csv("tail_DNA_mean_contacts_"+ file_name+ "_h3_E.dat",  sep = "\t", header=0)
        H4_contacts1 = pd.read_csv("tail_DNA_mean_contacts_"+ file_name+ "_h4_B.dat",  sep = "\t", header=0)
        H4_contacts2 = pd.read_csv("tail_DNA_mean_contacts_"+ file_name+ "_h4_F.dat",  sep = "\t", header=0)
        

## We only used the measurements from the long runs in calculations
        H2A_contacts = pd.concat([H2A_contacts1.iloc[0:19000,], H2A_contacts2.iloc[0:19000,]]).reset_index().iloc[:,1:]
        H2B_contacts = pd.concat([H2B_contacts1.iloc[0:19000,], H2B_contacts2.iloc[0:19000,]]).reset_index().iloc[:,1:]
        H3_contacts = pd.concat([H3_contacts1.iloc[0:19000,], H3_contacts2.iloc[0:19000,]]).reset_index().iloc[:,1:]
        H4_contacts = pd.concat([H4_contacts1.iloc[0:19000,], H4_contacts2.iloc[0:19000,]]).reset_index().iloc[:,1:]
            
        
        ## timestep for frame intervals
        ts = 0.2  #ns
        
        
        ## H2A_N tails 
        count_bound_time_N = 0
        
        ## H2A_C tails 
        count_bound_time_C = 0
        
        for i in range(0,len(H2A_contacts)-2):   
            contact0_N = filter(lambda x: x >= 1, H2A_contacts.iloc[i, 1:9])
            contact1_N = filter(lambda x: x >= 1, H2A_contacts.iloc[i+1, 1:9])
            
            ratio0_N = len(list(contact0_N))/8 ## percentage of bound tail residues of frame i
            ratio1_N  = len(list(contact1_N))/8 ## percentage of bound tail residues of frame i+i
            

            ## we combined two copeis of tail conformations for each runs run thus they are not continuous in time.
            
            if ratio0_N > cut_off_ratio: ## tail stay in bound state in frame i
                count_bound_time_N += 1                
                if i ==int(len(H2A_contacts))/2 or i ==len(H2A_contacts)-3: ## reach the max number of frames per run
                    tail_residence_time["H2A_N"].append(count_bound_time_N * ts)
                    count_bound_time_N = 0
                    continue
                elif ratio1_N <= cut_off_ratio: ## tail unbind from DNA in frame i+1
                    tail_residence_time["H2A_N"].append(count_bound_time_N * ts)
        
            elif ratio0_N <= cut_off_ratio: ## tail stay in unbound state in frame i
                count_bound_time_N = 0 ## reset residence time

            
        for i in range(0,len(H2A_contacts)-2):  
            contact0_C = filter(lambda x: x >= 1, H2A_contacts.iloc[i, 15:25])
            contact1_C = filter(lambda x: x >= 1, H2A_contacts.iloc[i+1, 15:25])
            
            ratio0_C = len(list(contact0_C))/9
            ratio1_C  = len(list(contact1_C))/9
            
            if ratio0_C > cut_off_ratio: ## tail stay in bound state in frame i
                count_bound_time_C += 1
                if i ==int(len(H2A_contacts))/2 or i ==len(H2A_contacts)-3: ## reach the max number of frames per run
                    tail_residence_time["H2A_C"].append(count_bound_time_C * ts)
                    count_bound_time_C = 0
                    continue
                elif ratio1_C <= cut_off_ratio: ## tail unbind from DNA in frame i+1
                    tail_residence_time["H2A_C"].append(count_bound_time_C * ts)
            elif ratio0_C <= cut_off_ratio:
                count_bound_time_C = 0 ## reset residence time

        
        ## H2B Tails
        count_bound_time = 0
        
        for i in range(0,len(H2A_contacts)-2):    
            contact0 = filter(lambda x: x >= 1, H2B_contacts.iloc[i, 1:19])
            contact1 = filter(lambda x: x >= 1, H2B_contacts.iloc[i+1, 1:19])
            
            ratio0 = len(list(contact0))/18
            ratio1  = len(list(contact1))/18
            
            
            if ratio0 > cut_off_ratio:
                count_bound_time += 1
                if i ==int(len(H2B_contacts))/2 or i ==len(H2B_contacts)-3:
                    tail_residence_time["H2B"].append(count_bound_time * ts)
                    count_bound_time = 0
                    continue
                elif ratio1 <= cut_off_ratio:
                    tail_residence_time["H2B"].append(count_bound_time * ts)
            elif ratio0 <= cut_off_ratio:
                count_bound_time = 0

            
        ## H3 tails
        count_bound_time = 0
        
        for i in range(0,len(H2A_contacts)-2):    
            contact0 = filter(lambda x: x >= 1, H3_contacts.iloc[i, 1:34])
            contact1 = filter(lambda x: x >= 1, H3_contacts.iloc[i+1, 1:34])
            
            ratio0 = len(list(contact0))/33
            ratio1  = len(list(contact1))/33
            
            if ratio0 > cut_off_ratio:
                count_bound_time += 1
                if i ==int(len(H3_contacts))/2 or i ==len(H3_contacts)-3:
                    tail_residence_time["H3"].append(count_bound_time * ts)
                    count_bound_time = 0
                    continue
                elif ratio1 <= cut_off_ratio:
                    tail_residence_time["H3"].append(count_bound_time * ts)
            elif ratio0 <= cut_off_ratio:
                count_bound_time = 0

            
        ## H4
        count_bound_time = 0
        
        for i in range(0,len(H2A_contacts)-2):    
            contact0 = filter(lambda x: x >= 1, H4_contacts.iloc[i, 1:16])
            contact1 = filter(lambda x: x >= 1, H4_contacts.iloc[i+1, 1:16])
            
            ratio0 = len(list(contact0))/15
            ratio1  = len(list(contact1))/15
        
            if ratio0 > cut_off_ratio:
                count_bound_time += 1
                if i ==int(len(H4_contacts))/2 or i ==len(H4_contacts)-3:
                    tail_residence_time["H4"].append(count_bound_time * ts)
                    count_bound_time = 0
                    continue
                elif ratio1 <= cut_off_ratio:
                    tail_residence_time["H4"].append(count_bound_time * ts)
            elif ratio0 <= cut_off_ratio:
                count_bound_time = 0

            
        
        
        df_residence_time = pd.DataFrame({"residence_time": tail_residence_time["H2A_N"],
                                         "tail_type": ["H2A_N"]*len(tail_residence_time["H2A_N"])})
        df_residence_time = df_residence_time.append(pd.DataFrame({"residence_time": tail_residence_time["H2A_C"],
                                         "tail_type": ["H2A_C"]*len(tail_residence_time["H2A_C"])}))
        df_residence_time = df_residence_time.append(pd.DataFrame({"residence_time": tail_residence_time["H2B"],
                                         "tail_type": ["H2B"]*len(tail_residence_time["H2B"])}))
        df_residence_time = df_residence_time.append(pd.DataFrame({"residence_time": tail_residence_time["H3"],
                                         "tail_type": ["H3"]*len(tail_residence_time["H3"])}))
        df_residence_time = df_residence_time.append(pd.DataFrame({"residence_time": tail_residence_time["H4"],
                                         "tail_type": ["H4"]*len(tail_residence_time["H4"])}))
        
        
        
        ## filter the reisdence time < 10 ns, Since full histone tails undergo very rapid fluctuations before retaining stable binding with DNA during the simulations
        cut_off_time = 10 
        df_residence_time_filter = df_residence_time.loc[df_residence_time['residence_time'] >= cut_off_time]
        
        
        with open(file_name+ "_full_tail_residence_time"+ str(cut_off_ratio) + ".csv", "w") as fwh:
            for i in range(0,len(df_residence_time_filter)):
               fwh.write(str(df_residence_time_filter.reset_index().iloc[i,2]) + "," \
                      + str(df_residence_time_filter.reset_index().iloc[i,1]) + "\n")
            
        
                
        
