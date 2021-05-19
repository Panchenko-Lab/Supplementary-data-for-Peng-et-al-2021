#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 14:48:45 2019

@author: pengy10
"""

import pandas as pd 
import numpy as np 
import re
import csv
import statistics 
import seaborn as sns
import matplotlib.pyplot as plt

#["1aoi","1eqz","ext2","sym", "shuxiang_ext_sym_amber_run1", "shuxiang_ext_sym_amber_run2"  ]:
for file_name in ["sym" ]:
    for cut_off in [0,0.1,0.2]:
        
        H2A_contacts1 = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+ file_name+ "_h2a_C.dat",  sep = "\t", header=0)
        H2A_contacts2 = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+ file_name+ "_h2a_G.dat",  sep = "\t", header=0)
        H2B_contacts1 = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+ file_name+ "_h2b_D.dat",  sep = "\t", header=0)
        H2B_contacts2 = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+ file_name+ "_h2b_H.dat",  sep = "\t", header=0)
        
        H3_contacts1 = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+ file_name+ "_h3_A.dat",  sep = "\t", header=0)
        H3_contacts2 = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+ file_name+ "_h3_E.dat",  sep = "\t", header=0)
        H4_contacts1 = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+ file_name+ "_h4_B.dat",  sep = "\t", header=0)
        H4_contacts2 = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+ file_name+ "_h4_F.dat",  sep = "\t", header=0)
        
        
        for H2A_contacts,  H2B_contacts,  H3_contacts, H4_contacts, chains in zip([H2A_contacts1.iloc[0:19000,],H2A_contacts2.iloc[0:19000,]],
                                                                                  [H2B_contacts1.iloc[0:19000,],H2B_contacts2.iloc[0:19000,]],
                                                                          [H3_contacts1.iloc[0:19000,],H3_contacts2.iloc[0:19000,]],
                                                                          [H4_contacts1.iloc[0:19000,],H4_contacts2.iloc[0:19000,]],[1,2]):

            
            tail_residence = { i : [] for i in ["H2A_N", "H2A_C", "H2B", "H3", "H4"] }
            tail_rebinding_count = { i : [] for i in ["H2A_N", "H2A_C", "H2B", "H3", "H4"]}
            unbound_frames_count = { i : [] for i in ["H2A_N", "H2A_C", "H2B", "H3", "H4"]}
            bound_frames_count = { i : [] for i in ["H2A_N", "H2A_C", "H2B", "H3", "H4"]}

        
            ## number of minimum contacts to define binding
            cut_off_contacts = 1
            ## number of minimum ratio to define full tail binding
            cut_off_ratio = cut_off
            ## time step for frames
            ts=0.2 #ns
        
            
            ## H2A 
            count_binding_N = 0
            count_rebinding_N = 0
            count_binding_C = 0
            count_rebinding_C = 0
            bound_frame_N = 0
            unbound_frame_N =0
            bound_frame_C = 0
            unbound_frame_C =0 
            
            for i in range(0,len(H2A_contacts)-2):    
                contact0_N = filter(lambda x: x >= cut_off_contacts, H2A_contacts.iloc[i, 1:9])
                contact1_N = filter(lambda x: x >= cut_off_contacts, H2A_contacts.iloc[i+1, 1:9])
                
                ratio0_N = len(list(contact0_N))/8
                ratio1_N  = len(list(contact1_N))/8
                
    
                ## we combined trajectories from five indenpent run thus they are not continuous in time.
                
                if ratio0_N > cut_off_ratio:
                    bound_frame_N +=1
                    count_binding_N += 1
                    if i ==len(H2A_contacts)-3:
                        tail_residence["H2A_N"].append(count_binding_N * ts)
                        count_binding_N = 0
                        continue
                    elif ratio1_N <= cut_off_ratio:
                        tail_residence["H2A_N"].append(count_binding_N * ts)
            
                elif ratio0_N <= cut_off_ratio:
                    count_binding_N = 0
                    unbound_frame_N +=1
                    if i ==len(H2A_contacts)-3:
                        continue
                    elif ratio1_N > cut_off_ratio:
                        count_rebinding_N += 1
            tail_rebinding_count["H2A_N"] = count_rebinding_N
            bound_frames_count["H2A_N"] = bound_frame_N
            unbound_frames_count["H2A_N"] = unbound_frame_N
            for i in range(0,len(H2A_contacts)-2):  
                contact0_C = filter(lambda x: x >= cut_off_contacts, H2A_contacts.iloc[i, 15:25])
                contact1_C = filter(lambda x: x >= cut_off_contacts, H2A_contacts.iloc[i+1, 15:25])
                
                ratio0_C = len(list(contact0_C))/9
                ratio1_C  = len(list(contact1_C))/9
                
                if ratio0_C > cut_off_ratio:
                    bound_frame_C +=1
                    count_binding_C += 1
                    if  i ==len(H2A_contacts)-3:
                        tail_residence["H2A_C"].append(count_binding_C * ts)
                        count_binding_C = 0
                        continue
                    elif ratio1_C <= cut_off_ratio:
                        tail_residence["H2A_C"].append(count_binding_C * ts)
                elif ratio0_C <= cut_off_ratio:
                    unbound_frame_C +=1
                    count_binding_C = 0
                    if i ==len(H2A_contacts)-3:
                        continue
                    elif ratio1_C > cut_off_ratio:
                        count_rebinding_C += 1
            tail_rebinding_count["H2A_C"] = count_rebinding_C
            bound_frames_count["H2A_C"] = bound_frame_C
            unbound_frames_count["H2A_C"] = unbound_frame_C
            
            ## H2B
            count_binding = 0
            count_rebinding = 0
            bound_frame =0
            unbound_frame=0
            
            for i in range(0,len(H2A_contacts)-2):    
                contact0 = filter(lambda x: x >= cut_off_contacts, H2B_contacts.iloc[i, 1:19])
                contact1 = filter(lambda x: x >= cut_off_contacts, H2B_contacts.iloc[i+1, 1:19])
                
                ratio0 = len(list(contact0))/18
                ratio1  = len(list(contact1))/18
                
                
                if ratio0 > cut_off_ratio:
                    count_binding += 1
                    bound_frame +=1
                    if  i ==len(H2A_contacts)-3:
                        tail_residence["H2B"].append(count_binding * ts)
                        count_binding = 0
                        continue
                    elif ratio1 <= cut_off_ratio:
                        tail_residence["H2B"].append(count_binding * ts)
                elif ratio0 <= cut_off_ratio:
                    count_binding = 0
                    unbound_frame +=1
                    if i ==len(H2A_contacts)-3:
                        continue
                    elif ratio1 > cut_off_ratio:
                        count_rebinding += 1
            tail_rebinding_count["H2B"] = count_rebinding
            bound_frames_count["H2B"] = bound_frame
            unbound_frames_count["H2B"] = unbound_frame
                
            ## H3
            count_binding = 0
            count_rebinding = 0
            bound_frame =0
            unbound_frame=0
            
            for i in range(0,len(H2A_contacts)-2):    
                contact0 = filter(lambda x: x >= cut_off_contacts, H3_contacts.iloc[i, 1:34])
                contact1 = filter(lambda x: x >= cut_off_contacts, H3_contacts.iloc[i+1, 1:34])
                
                ratio0 = len(list(contact0))/33
                ratio1  = len(list(contact1))/33
                
                if ratio0 > cut_off_ratio:
                    count_binding += 1
                    bound_frame +=1
                    if  i ==len(H2A_contacts)-3:
                        tail_residence["H3"].append(count_binding * ts)
                        count_binding = 0
                        continue
                    elif ratio1 <= cut_off_ratio:
                        tail_residence["H3"].append(count_binding * ts)
                elif ratio0 <= cut_off_ratio:
                    count_binding = 0
                    unbound_frame +=1
                    if i ==len(H2A_contacts)-3:
                        continue
                    elif ratio1 > cut_off_ratio:
                        count_rebinding += 1
            tail_rebinding_count["H3"] = count_rebinding
            bound_frames_count["H3"] = bound_frame
            unbound_frames_count["H3"] = unbound_frame
                
            ## H4
            count_binding = 0
            count_rebinding = 0
            bound_frame =0
            unbound_frame =0
            
            for i in range(0,len(H2A_contacts)-2):    
                contact0 = filter(lambda x: x >= cut_off_contacts, H4_contacts.iloc[i, 1:16])
                contact1 = filter(lambda x: x >= cut_off_contacts, H4_contacts.iloc[i+1, 1:16])
                
                ratio0 = len(list(contact0))/15
                ratio1  = len(list(contact1))/15
            
                if ratio0 > cut_off_ratio:
                    count_binding += 1
                    bound_frame +=1
                    if i ==len(H2A_contacts)-3:
                        tail_residence["H4"].append(count_binding * ts)
                        count_binding = 0
                        continue
                    elif ratio1 <= cut_off_ratio:
                        tail_residence["H4"].append(count_binding * ts)
                elif ratio0 <= cut_off_ratio:
                    count_binding = 0
                    unbound_frame +=1
                    if i ==len(H2A_contacts)-3:
                        continue
                    elif ratio1 > cut_off_ratio:
                        count_rebinding += 1
            tail_rebinding_count["H4"] = count_rebinding
            bound_frames_count["H4"] = bound_frame
            unbound_frames_count["H4"] = unbound_frame
                
            
            
            df_residence_time = pd.DataFrame({"residence_time": tail_residence["H2A_N"],
                                             "tail_type": ["H2A_N"]*len(tail_residence["H2A_N"])})
            df_residence_time = df_residence_time.append(pd.DataFrame({"residence_time": tail_residence["H2A_C"],
                                             "tail_type": ["H2A_C"]*len(tail_residence["H2A_C"])}))
            df_residence_time = df_residence_time.append(pd.DataFrame({"residence_time": tail_residence["H2B"],
                                             "tail_type": ["H2B"]*len(tail_residence["H2B"])}))
            df_residence_time = df_residence_time.append(pd.DataFrame({"residence_time": tail_residence["H3"],
                                             "tail_type": ["H3"]*len(tail_residence["H3"])}))
            df_residence_time = df_residence_time.append(pd.DataFrame({"residence_time": tail_residence["H4"],
                                             "tail_type": ["H4"]*len(tail_residence["H4"])}))
            
            
            
            ## filter the reisdence time < 10 ns, maybe due to flucuations
            cut_off_time = 10 
            df_residence_time_filter = df_residence_time.loc[df_residence_time['residence_time'] >= cut_off_time]
            count_H2A_N = len(df_residence_time_filter.loc[df_residence_time_filter['tail_type'] == "H2A_N"])
            count_H2A_C = len(df_residence_time_filter.loc[df_residence_time_filter['tail_type'] == "H2A_C"])
            count_H2B = len(df_residence_time_filter.loc[df_residence_time_filter['tail_type'] == "H2B"])
            count_H3 = len(df_residence_time_filter.loc[df_residence_time_filter['tail_type'] == "H3"])
            count_H4 = len(df_residence_time_filter.loc[df_residence_time_filter['tail_type'] == "H4"])
            unbinding_count_filter = {'H2A_N' :count_H2A_N,'H2A_C':count_H2A_C,'H2B': count_H2B, 'H3':count_H3,'H4': count_H4}
            
            
                
            with open("./Kd/"+ file_name+ "_Kd_"+ str(cut_off) + "_" + str(chains) + ".csv", "w") as fwh:
                 for keys in tail_rebinding_count:
                     fwh.write(keys +"," + str(unbound_frames_count[keys]) + "," + str(bound_frames_count[keys]) + "," +\
                               str(tail_rebinding_count[keys]) + "," + str(unbinding_count_filter[keys]) + "\n")
            
