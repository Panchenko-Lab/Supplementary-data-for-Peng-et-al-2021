#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 14:48:45 2019

@author: pengy10
"""

import pandas as pd 
import os 
import re
import csv
import statistics 
import numpy as np

H2A_residence_mean = { i : [] for i in [*range(1,14),*range(119,129)]}
H2B_residence_mean = { i : [] for i in range(1,24) }
H3_residence_mean = { i : [] for i in range(1,37) }
H4_residence_mean = { i : [] for i in range(1,21) }    

H2A_residence_stdev_all_run = { i : 0 for i in [*range(1,14),*range(119,129)]}
H2B_residence_stdev_all_run = { i : 0 for i in range(1,24) }
H3_residence_stdev_all_run = { i : 0 for i in range(1,37) }
H4_residence_stdev_all_run = { i : 0 for i in range(1,21) } 

H2A_residence_mean_all_run = { i : 0 for i in [*range(1,14),*range(119,129)]}
H2B_residence_mean_all_run = { i : 0 for i in range(1,24) }
H3_residence_mean_all_run = { i : 0 for i in range(1,37) }
H4_residence_mean_all_run = { i : 0 for i in range(1,21) }


## Amber run
for (name,  m) in  zip(["1aoi", "1eqz", "ext2", "sym"],[21500,24000,21500,19000]):
    
#    print(name,m)
        
        H2A_contacts_all_run_C = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h2a_C.dat",  sep = "\t", header=0)
        H2A_contacts_all_run_G = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h2a_G.dat",  sep = "\t", header=0)
        H2B_contacts_all_run_D = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h2b_D.dat",  sep = "\t", header=0)
        H2B_contacts_all_run_H = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h2b_H.dat",  sep = "\t", header=0)
        H3_contacts_all_run_A = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h3_A.dat",  sep = "\t", header=0)
        H3_contacts_all_run_E = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h3_E.dat",  sep = "\t", header=0)
        H4_contacts_all_run_B = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h4_B.dat",  sep = "\t", header=0)
        H4_contacts_all_run_F = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h4_F.dat",  sep = "\t", header=0)
        
        H2A_contacts_run1 = pd.concat([H2A_contacts_all_run_C.iloc[0:m,],H2A_contacts_all_run_G.iloc[0:m,]])
        H2A_contacts_run2 = pd.concat([H2A_contacts_all_run_C.iloc[m:m+3750,],H2A_contacts_all_run_G.iloc[m:m+3750,]])
        H2A_contacts_run3 = pd.concat([H2A_contacts_all_run_C.iloc[m+3750:m+3750*2,],H2A_contacts_all_run_G.iloc[m+3750:m+3750*2,]])
        H2A_contacts_run4 = pd.concat([H2A_contacts_all_run_C.iloc[m+3750*2:m+3750*3,],H2A_contacts_all_run_G.iloc[m+3750*2:m+3750*3,]])
        H2A_contacts_run5 = pd.concat([H2A_contacts_all_run_C.iloc[m+3750*3:m+3750*4,],H2A_contacts_all_run_G.iloc[m+3750*3:m+3750*4,]])
        
        H2B_contacts_run1 = pd.concat([H2B_contacts_all_run_D.iloc[0:m,],H2B_contacts_all_run_H.iloc[0:m,]])
        H2B_contacts_run2 = pd.concat([H2B_contacts_all_run_D.iloc[m:m+3750,],H2B_contacts_all_run_H.iloc[m:m+3750,]])
        H2B_contacts_run3 = pd.concat([H2B_contacts_all_run_D.iloc[m+3750:m+3750*2,],H2B_contacts_all_run_H.iloc[m+3750:m+3750*2,]])
        H2B_contacts_run4 = pd.concat([H2B_contacts_all_run_D.iloc[m+3750*2:m+3750*3,],H2B_contacts_all_run_H.iloc[m+3750*2:m+3750*3,]])
        H2B_contacts_run5 = pd.concat([H2B_contacts_all_run_D.iloc[m+3750*3:m+3750*4,],H2B_contacts_all_run_H.iloc[m+3750*3:m+3750*4,]])

        H3_contacts_run1 = pd.concat([H3_contacts_all_run_A.iloc[0:m,],H3_contacts_all_run_E.iloc[0:m,]])
        H3_contacts_run2 = pd.concat([H3_contacts_all_run_A.iloc[m:m+3750,],H3_contacts_all_run_E.iloc[m:m+3750,]])
        H3_contacts_run3 = pd.concat([H3_contacts_all_run_A.iloc[m+3750:m+3750*2,],H3_contacts_all_run_E.iloc[m+3750:m+3750*2,]])
        H3_contacts_run4 = pd.concat([H3_contacts_all_run_A.iloc[m+3750*2:m+3750*3,],H3_contacts_all_run_E.iloc[m+3750*2:m+3750*3,]])
        H3_contacts_run5 = pd.concat([H3_contacts_all_run_A.iloc[m+3750*3:m+3750*4,],H3_contacts_all_run_E.iloc[m+3750*3:m+3750*4,]])

        H4_contacts_run1 = pd.concat([H4_contacts_all_run_B.iloc[0:m,],H4_contacts_all_run_F.iloc[0:m,]])
        H4_contacts_run2 = pd.concat([H4_contacts_all_run_B.iloc[m:m+3750,],H4_contacts_all_run_F.iloc[m:m+3750,]])
        H4_contacts_run3 = pd.concat([H4_contacts_all_run_B.iloc[m+3750:m+3750*2,],H4_contacts_all_run_F.iloc[m+3750:m+3750*2,]])
        H4_contacts_run4 = pd.concat([H4_contacts_all_run_B.iloc[m+3750*2:m+3750*3,],H4_contacts_all_run_F.iloc[m+3750*2:m+3750*3,]])
        H4_contacts_run5 = pd.concat([H4_contacts_all_run_B.iloc[m+3750*3:m+3750*4,],H4_contacts_all_run_F.iloc[m+3750*3:m+3750*4,]])
       

        for (H2A_contacts,  H2B_contacts, H3_contacts, H4_contacts) in  zip(
                [H2A_contacts_run1,H2A_contacts_run2,H2A_contacts_run3,H2A_contacts_run4,H2A_contacts_run5],
                [H2B_contacts_run1,H2B_contacts_run2,H2B_contacts_run3,H2B_contacts_run4,H2B_contacts_run5],
                [H3_contacts_run1,H3_contacts_run2,H3_contacts_run3,H3_contacts_run4,H3_contacts_run5],
                [H4_contacts_run1,H4_contacts_run2,H4_contacts_run3,H4_contacts_run4,H4_contacts_run5]):
                        
            H2A_residence = { i : [] for i in [*range(1,14),*range(119,129)]}
            H2B_residence = { i : [] for i in range(1,24) }
            H3_residence = { i : [] for i in range(1,37) }
            H4_residence = { i : [] for i in range(1,21) }
            
            H2A_rebinding_count = { i : [] for i in [*range(1,14),*range(119,129)]}
            H2B_rebinding_count = { i : [] for i in range(1,24) }
            H3_rebinding_count = { i : [] for i in range(1,37) }
            H4_rebinding_count = { i : [] for i in range(1,21) }
        
            ## number of minimum contacts to define binding
            cut_off = 1
            ## time step for frames
            ts=0.2 #ns
            
            for key in H2A_residence:
                count_binding = 0
                count_rebinding = 0
            #    print(key)
                for i in range(0,len(H2A_contacts[str(key)])):                    
                ## we combined trajectories from five indenpent run thus they are not continuous in time.
                    if i == len(H2A_contacts[str(key)])-1 or i == int(len(H2A_contacts[str(key)])/2) :
                        H2A_residence[key].append(count_binding*ts)
                        count_binding = 0
            ## tail continuously bound with DNA in frame i
                    if H2A_contacts[str(key)].iloc[i] >= cut_off and i <= len(H2A_contacts[str(key)])-2:
                        count_binding += 1 
                        #tail unbinding in next frame i+1
                        if H2A_contacts[str(key)].iloc[i + 1] < cut_off:
                            #save the binding time
                            H2A_residence[key].append(count_binding*ts)
                    # tail unbind in frame i
                    elif H2A_contacts[str(key)].iloc[i] < cut_off and i <= len(H2A_contacts[str(key)])-2:
                        # restart counting binding time
                        count_binding = 0
                        ## Tail rebound in next frame i+1
                        if H2A_contacts[str(key)].iloc[i+1] >= cut_off:
                            count_rebinding += 1
                H2A_rebinding_count[key] = count_rebinding
            
                            
            for key in H2B_residence:
                count_binding = 0
                count_rebinding = 0
            #    print(key)
                for i in range(0,len(H2B_contacts[str(key)])):
                ## we combined trajectories from five indenpent run thus they are not continuous in time.
                    if i == len(H2B_contacts[str(key)])-1 or i == int(len(H2B_contacts[str(key)])/2) :
                        H2B_residence[key].append(count_binding*ts)
                        count_binding = 0
            ## tail continuously bound with DNA in frame i
                    if H2B_contacts[str(key)].iloc[i] >= cut_off and i <= len(H2B_contacts[str(key)])-2:
                        count_binding += 1 
                        #tail unbinding in next frame i+1
                        if H2B_contacts[str(key)].iloc[i + 1] < cut_off:
                            #save the binding time
                            H2B_residence[key].append(count_binding*ts)
                    # tail unbind in frame i
                    elif H2B_contacts[str(key)].iloc[i] < cut_off and i <= len(H2B_contacts[str(key)])-2:
                        # restart counting binding time
                        count_binding = 0
                        ## Tail rebound in next frame i+1
                        if H2B_contacts[str(key)].iloc[i+1] >= cut_off:
                            count_rebinding += 1
                H2B_rebinding_count[key] = count_rebinding
            
            
            for key in H3_residence:
                count_binding = 0
                count_rebinding = 0
            #    print(key)
                for i in range(0,len(H3_contacts[str(key)])):
                ## we combined trajectories from five indenpent run thus they are not continuous in time.
                    if i == len(H3_contacts[str(key)])-1 or i == int(len(H3_contacts[str(key)])/2) :
                        H3_residence[key].append(count_binding*ts)
                        count_binding = 0
            ## tail continuously bound with DNA in frame i
                    if H3_contacts[str(key)].iloc[i] >= cut_off and i <= len(H3_contacts[str(key)])-2:
                        count_binding += 1 
                        #tail unbinding in next frame i+1
                        if H3_contacts[str(key)].iloc[i + 1] < cut_off:
                            #save the binding time
                            H3_residence[key].append(count_binding*ts)
                    # tail unbind in frame i
                    elif H3_contacts[str(key)].iloc[i] < cut_off and i <= len(H3_contacts[str(key)])-2:
                        # restart counting binding time
                        count_binding = 0
                        ## Tail rebound in next frame i+1
                        if H3_contacts[str(key)].iloc[i+1] >= cut_off:
                            count_rebinding += 1
                H3_rebinding_count[key] = count_rebinding
            
            for key in H4_residence:
                count_binding = 0
                count_rebinding = 0
            #    print(key)
                for i in range(0,len(H4_contacts[str(key)])):
                ## we combined trajectories from five indenpent run thus they are not continuous in time.
                    if i == len(H4_contacts[str(key)])-1 or i == int(len(H4_contacts[str(key)])/2) :
                        H4_residence[key].append(count_binding*ts)
                        count_binding = 0
            ## tail continuously bound with DNA in frame i
                    if H4_contacts[str(key)].iloc[i] >= cut_off and i <= len(H4_contacts[str(key)])-2:
                        count_binding += 1 
                        #tail unbinding in next frame i+1
                        if H4_contacts[str(key)].iloc[i + 1] < cut_off:
                            #save the binding time
                            H4_residence[key].append(count_binding*ts)
                    # tail unbind in frame i
                    elif H4_contacts[str(key)].iloc[i] < cut_off and i <= len(H4_contacts[str(key)])-2:
                        # restart counting binding time
                        count_binding = 0
                        ## Tail rebound in next frame i+1
                        if H4_contacts[str(key)].iloc[i+1] >= cut_off:
                            count_rebinding += 1
                H4_rebinding_count[key] = count_rebinding
                
            
            
            # find the mean residence time for each run.
            cut_off = 740
            for key in H2A_residence_mean:
                if len(sorted(x for x in H2A_residence[key] if x <= cut_off)) >=1:
                   
                    H2A_residence_mean[key].append(statistics.mean(sorted(x for x in H2A_residence[key] if x <= cut_off)))
                
            for key in H2B_residence_mean:
                if len(sorted(x for x in H2B_residence[key] if x <= cut_off)) >=1:
                    H2B_residence_mean[key].append(statistics.mean(sorted(x for x in H2B_residence[key] if x <= cut_off)))
            
            for key in H3_residence_mean:
                if len(sorted(x for x in H3_residence[key] if x <= cut_off)) >=1:
                    H3_residence_mean[key].append(statistics.mean(sorted(x for x in H3_residence[key] if x <= cut_off)))
                
            for key in H4_residence_mean:
                if len(sorted(x for x in H4_residence[key] if x <= cut_off)) >=1:
                    H4_residence_mean[key].append(statistics.mean(sorted(x for x in H4_residence[key] if x <= cut_off)))
            print(H2A_residence_mean)
            print(H2B_residence_mean)
            print(H3_residence_mean)
            print(H4_residence_mean)

#%%
## Gromacs  run
for name in ["shuxiang_ext_sym_amber_run1","shuxiang_ext_sym_amber_run2"]:
        
        H2A_contacts_all_run_C = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h2a_C.dat",  sep = "\t", header=0)
        H2A_contacts_all_run_G = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h2a_G.dat",  sep = "\t", header=0)
        H2B_contacts_all_run_D = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h2b_D.dat",  sep = "\t", header=0)
        H2B_contacts_all_run_H = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h2b_H.dat",  sep = "\t", header=0)
        H3_contacts_all_run_A = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h3_A.dat",  sep = "\t", header=0)
        H3_contacts_all_run_E = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h3_E.dat",  sep = "\t", header=0)
        H4_contacts_all_run_B = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h4_B.dat",  sep = "\t", header=0)
        H4_contacts_all_run_F = pd.read_csv("./tail_DNA_contacts_4500ns/tail_contacts_all_"+name+"_h4_F.dat",  sep = "\t", header=0)
        
        H2A_contacts_run1 = pd.concat([H2A_contacts_all_run_C.iloc[0:24000,],H2A_contacts_all_run_G.iloc[0:24000,]])
        H2B_contacts_run1 = pd.concat([H2B_contacts_all_run_D.iloc[0:24000,],H2B_contacts_all_run_H.iloc[0:24000,]])
        H3_contacts_run1 = pd.concat([H3_contacts_all_run_A.iloc[0:24000,],H3_contacts_all_run_E.iloc[0:24000,]])
        H4_contacts_run1 = pd.concat([H4_contacts_all_run_B.iloc[0:24000,],H4_contacts_all_run_F.iloc[0:24000,]])


        for (H2A_contacts,  H2B_contacts, H3_contacts, H4_contacts) in  zip(
                [H2A_contacts_run1],
                [H2B_contacts_run1],
                [H3_contacts_run1],
                [H4_contacts_run1]):
                        
            H2A_residence = { i : [] for i in [*range(1,14),*range(119,129)]}
            H2B_residence = { i : [] for i in range(1,24) }
            H3_residence = { i : [] for i in range(1,37) }
            H4_residence = { i : [] for i in range(1,21) }
            
            H2A_rebinding_count = { i : [] for i in [*range(1,14),*range(119,129)]}
            H2B_rebinding_count = { i : [] for i in range(1,24) }
            H3_rebinding_count = { i : [] for i in range(1,37) }
            H4_rebinding_count = { i : [] for i in range(1,21) }
        
            ## number of minimum contacts to define binding
            cut_off = 1
            ## time step for frames
            ts=0.2 #ns
            
            for key in H2A_residence:
                count_binding = 0
                count_rebinding = 0
            #    print(key)
                for i in range(0,len(H2A_contacts[str(key)])):                    
                ## we combined trajectories from five indenpent run thus they are not continuous in time.
                    if i == len(H2A_contacts[str(key)])-1 or i == int(len(H2A_contacts[str(key)])/2) :
                        H2A_residence[key].append(count_binding*ts)
                        count_binding = 0
            ## tail continuously bound with DNA in frame i
                    if H2A_contacts[str(key)].iloc[i] >= cut_off and i <= len(H2A_contacts[str(key)])-2:
                        count_binding += 1 
                        #tail unbinding in next frame i+1
                        if H2A_contacts[str(key)].iloc[i + 1] < cut_off:
                            #save the binding time
                            H2A_residence[key].append(count_binding*ts)
                    # tail unbind in frame i
                    elif H2A_contacts[str(key)].iloc[i] < cut_off and i <= len(H2A_contacts[str(key)])-2:
                        # restart counting binding time
                        count_binding = 0
                        ## Tail rebound in next frame i+1
                        if H2A_contacts[str(key)].iloc[i+1] >= cut_off:
                            count_rebinding += 1
                H2A_rebinding_count[key] = count_rebinding
            
                            
            for key in H2B_residence:
                count_binding = 0
                count_rebinding = 0
            #    print(key)
                for i in range(0,len(H2B_contacts[str(key)])):
                ## we combined trajectories from five indenpent run thus they are not continuous in time.
                    if i == len(H2B_contacts[str(key)])-1 or i == int(len(H2B_contacts[str(key)])/2) :
                        H2B_residence[key].append(count_binding*ts)
                        count_binding = 0
            ## tail continuously bound with DNA in frame i
                    if H2B_contacts[str(key)].iloc[i] >= cut_off and i <= len(H2B_contacts[str(key)])-2:
                        count_binding += 1 
                        #tail unbinding in next frame i+1
                        if H2B_contacts[str(key)].iloc[i + 1] < cut_off:
                            #save the binding time
                            H2B_residence[key].append(count_binding*ts)
                    # tail unbind in frame i
                    elif H2B_contacts[str(key)].iloc[i] < cut_off and i <= len(H2B_contacts[str(key)])-2:
                        # restart counting binding time
                        count_binding = 0
                        ## Tail rebound in next frame i+1
                        if H2B_contacts[str(key)].iloc[i+1] >= cut_off:
                            count_rebinding += 1
                H2B_rebinding_count[key] = count_rebinding
            
            
            for key in H3_residence:
                count_binding = 0
                count_rebinding = 0
            #    print(key)
                for i in range(0,len(H3_contacts[str(key)])):
                ## we combined trajectories from five indenpent run thus they are not continuous in time.
                    if i == len(H3_contacts[str(key)])-1 or i == int(len(H3_contacts[str(key)])/2) :
                        H3_residence[key].append(count_binding*ts)
                        count_binding = 0
            ## tail continuously bound with DNA in frame i
                    if H3_contacts[str(key)].iloc[i] >= cut_off and i <= len(H3_contacts[str(key)])-2:
                        count_binding += 1 
                        #tail unbinding in next frame i+1
                        if H3_contacts[str(key)].iloc[i + 1] < cut_off:
                            #save the binding time
                            H3_residence[key].append(count_binding*ts)
                    # tail unbind in frame i
                    elif H3_contacts[str(key)].iloc[i] < cut_off and i <= len(H3_contacts[str(key)])-2:
                        # restart counting binding time
                        count_binding = 0
                        ## Tail rebound in next frame i+1
                        if H3_contacts[str(key)].iloc[i+1] >= cut_off:
                            count_rebinding += 1
                H3_rebinding_count[key] = count_rebinding
            
            for key in H4_residence:
                count_binding = 0
                count_rebinding = 0
            #    print(key)
                for i in range(0,len(H4_contacts[str(key)])):
                ## we combined trajectories from five indenpent run thus they are not continuous in time.
                    if i == len(H4_contacts[str(key)])-1 or i == int(len(H4_contacts[str(key)])/2) :
                        H4_residence[key].append(count_binding*ts)
                        count_binding = 0
            ## tail continuously bound with DNA in frame i
                    if H4_contacts[str(key)].iloc[i] >= cut_off and i <= len(H4_contacts[str(key)])-2:
                        count_binding += 1 
                        #tail unbinding in next frame i+1
                        if H4_contacts[str(key)].iloc[i + 1] < cut_off:
                            #save the binding time
                            H4_residence[key].append(count_binding*ts)
                    # tail unbind in frame i
                    elif H4_contacts[str(key)].iloc[i] < cut_off and i <= len(H4_contacts[str(key)])-2:
                        # restart counting binding time
                        count_binding = 0
                        ## Tail rebound in next frame i+1
                        if H4_contacts[str(key)].iloc[i+1] >= cut_off:
                            count_rebinding += 1
                H4_rebinding_count[key] = count_rebinding
                
            
            
            # find the mean residence time for each run.
            cut_off = 740
            for key in H2A_residence_mean:
                if len(sorted(x for x in H2A_residence[key] if x <= cut_off)) >=1:
                   
                    H2A_residence_mean[key].append(statistics.mean(sorted(x for x in H2A_residence[key] if x <= cut_off)))
                
            for key in H2B_residence_mean:
                if len(sorted(x for x in H2B_residence[key] if x <= cut_off)) >=1:
                    H2B_residence_mean[key].append(statistics.mean(sorted(x for x in H2B_residence[key] if x <= cut_off)))
            
            for key in H3_residence_mean:
                if len(sorted(x for x in H3_residence[key] if x <= cut_off)) >=1:
                    H3_residence_mean[key].append(statistics.mean(sorted(x for x in H3_residence[key] if x <= cut_off)))
                
            for key in H4_residence_mean:
                if len(sorted(x for x in H4_residence[key] if x <= cut_off)) >=1:
                    H4_residence_mean[key].append(statistics.mean(sorted(x for x in H4_residence[key] if x <= cut_off)))
            print(H2A_residence_mean)
            print(H2B_residence_mean)
            print(H3_residence_mean)
            print(H4_residence_mean)



#%%                

for key in H2A_residence_stdev_all_run:
    H2A_residence_stdev_all_run[key] = statistics.stdev(H2A_residence_mean[key])/np.sqrt(len(H2A_residence_mean[key]))
    
for key in H2B_residence_stdev_all_run:
    H2B_residence_stdev_all_run[key] = statistics.stdev(H2B_residence_mean[key])/np.sqrt(len(H2B_residence_mean[key]))

for key in H3_residence_stdev_all_run:
    H3_residence_stdev_all_run[key] = statistics.stdev(H3_residence_mean[key])/np.sqrt(len(H3_residence_mean[key]))
    
for key in H4_residence_stdev_all_run:
    H4_residence_stdev_all_run[key] = statistics.stdev(H4_residence_mean[key])/np.sqrt(len(H4_residence_mean[key]))
    

for key in H2A_residence_mean_all_run:
    H2A_residence_mean_all_run[key] = statistics.mean(H2A_residence_mean[key])
    
for key in H2B_residence_mean_all_run:
    H2B_residence_mean_all_run[key] = statistics.mean(H2B_residence_mean[key])

for key in H3_residence_mean_all_run:
    H3_residence_mean_all_run[key] = statistics.mean(H3_residence_mean[key])
    
for key in H4_residence_mean_all_run:
    H4_residence_mean_all_run[key] = statistics.mean(H4_residence_mean[key])



w = csv.writer(open("./residue_residence_time/h2a_residence_mean_1_all_run_cut_off_4500ns.csv", "w"))
for mean, stdev in zip(H2A_residence_mean_all_run.items(),H2A_residence_stdev_all_run.items()):
    w.writerow([mean[0], mean[1], stdev[1]])
    
w = csv.writer(open("./residue_residence_time/h2b_residence_mean_1_all_run_cut_off_4500ns.csv", "w"))
for mean, stdev in zip(H2B_residence_mean_all_run.items(),H2B_residence_stdev_all_run.items()):
    w.writerow([mean[0], mean[1], stdev[1]])
    
w = csv.writer(open("./residue_residence_time/h3_residence_mean_1_all_run_cut_off_4500ns.csv", "w"))
for mean, stdev in zip(H3_residence_mean_all_run.items(),H3_residence_stdev_all_run.items()):
    w.writerow([mean[0], mean[1], stdev[1]])
    
w = csv.writer(open("./residue_residence_time/h4_residence_mean_1_all_run_cut_off_4500ns.csv", "w"))
for mean, stdev in zip(H4_residence_mean_all_run.items(),H4_residence_stdev_all_run.items()):
    w.writerow([mean[0], mean[1], stdev[1]])
    
            
            
            
                
                     
            
