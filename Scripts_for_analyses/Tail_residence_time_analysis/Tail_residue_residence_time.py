#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 14:48:45 2019

@author: pengy10
"""

import pandas as pd 
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
for (name,  m) in  zip(["ModelA", "ModelB", "ModelC", "ModelD"],[21500,24000,21500,19000]):
            
        H2A_contacts_all_run_C = pd.read_csv("tail_DNA_mean_contacts_"+name+"_h2a_C.dat",  sep = "\t", header=0)
        H2A_contacts_all_run_G = pd.read_csv("tail_DNA_mean_contacts_"+name+"_h2a_G.dat",  sep = "\t", header=0)
        H2B_contacts_all_run_D = pd.read_csv("tail_DNA_mean_contacts_"+name+"_h2b_D.dat",  sep = "\t", header=0)
        H2B_contacts_all_run_H = pd.read_csv("tail_DNA_mean_contacts_"+name+"_h2b_H.dat",  sep = "\t", header=0)
        H3_contacts_all_run_A = pd.read_csv("tail_DNA_mean_contacts_"+name+"_h3_A.dat",  sep = "\t", header=0)
        H3_contacts_all_run_E = pd.read_csv("tail_DNA_mean_contacts_"+name+"_h3_E.dat",  sep = "\t", header=0)
        H4_contacts_all_run_B = pd.read_csv("tail_DNA_mean_contacts_"+name+"_h4_B.dat",  sep = "\t", header=0)
        H4_contacts_all_run_F = pd.read_csv("tail_DNA_mean_contacts_"+name+"_h4_F.dat",  sep = "\t", header=0)
        
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
                        
            H2A_residence_time = { i : [] for i in [*range(1,14),*range(119,129)]}
            H2B_residence_time = { i : [] for i in range(1,24) }
            H3_residence_time = { i : [] for i in range(1,37) }
            H4_residence_time = { i : [] for i in range(1,21) }

        
            ## number of minimum contacts for tail residue bound state with DNA
            cut_off = 1
            ## timestep for frame intervals
            ts=0.2  #ns
            
            ## H2A tails
            for key in H2A_residence_time:
                count_bound_time = 0
                
                for i in range(0,len(H2A_contacts[str(key)])):                    
                ## we combined trajectories from five indenpent run thus they are not continuous in time.
                    if i == len(H2A_contacts[str(key)])-1 or i == int(len(H2A_contacts[str(key)])/2) :
                        H2A_residence_time[key].append(count_bound_time*ts)
                        count_bound_time = 0
            ## tail continuously bound with DNA in frame i
                    if H2A_contacts[str(key)].iloc[i] >= cut_off and i <= len(H2A_contacts[str(key)])-2:
                        count_bound_time += 1 
                        #tail unbinding in next frame i+1
                        if H2A_contacts[str(key)].iloc[i + 1] < cut_off:
                            #save the residence time
                            H2A_residence_time[key].append(count_bound_time*ts)
                    # tail unbind in frame i
                    elif H2A_contacts[str(key)].iloc[i] < cut_off and i <= len(H2A_contacts[str(key)])-2:
                        # restart counting binding time
                        count_bound_time = 0
            
            ## H2B tails                
            for key in H2B_residence_time:
                count_bound_time = 0
                
                for i in range(0,len(H2B_contacts[str(key)])):
                ## we combined trajectories from five indenpent run thus they are not continuous in time.
                    if i == len(H2B_contacts[str(key)])-1 or i == int(len(H2B_contacts[str(key)])/2) :
                        H2B_residence_time[key].append(count_bound_time*ts)
                        count_bound_time = 0
            ## tail continuously bound with DNA in frame i
                    if H2B_contacts[str(key)].iloc[i] >= cut_off and i <= len(H2B_contacts[str(key)])-2:
                        count_bound_time += 1 
                        #tail unbinding in next frame i+1
                        if H2B_contacts[str(key)].iloc[i + 1] < cut_off:
                            #save the residence time
                            H2B_residence_time[key].append(count_bound_time*ts)
                    # tail unbind in frame i
                    elif H2B_contacts[str(key)].iloc[i] < cut_off and i <= len(H2B_contacts[str(key)])-2:
                        # restart counting binding time
                        count_bound_time = 0
            
            ## H3 tails
            for key in H3_residence_time:
                count_bound_time = 0

                for i in range(0,len(H3_contacts[str(key)])):
                ## we combined trajectories from five indenpent run thus they are not continuous in time.
                    if i == len(H3_contacts[str(key)])-1 or i == int(len(H3_contacts[str(key)])/2) :
                        H3_residence_time[key].append(count_bound_time*ts)
                        count_bound_time = 0
            ## tail continuously bound with DNA in frame i
                    if H3_contacts[str(key)].iloc[i] >= cut_off and i <= len(H3_contacts[str(key)])-2:
                        count_bound_time += 1 
                        #tail unbinding in next frame i+1
                        if H3_contacts[str(key)].iloc[i + 1] < cut_off:
                            #save the residence time
                            H3_residence_time[key].append(count_bound_time*ts)
                    # tail unbind in frame i
                    elif H3_contacts[str(key)].iloc[i] < cut_off and i <= len(H3_contacts[str(key)])-2:
                        # restart counting binding time
                        count_bound_time = 0

            ## H4 tails
            for key in H4_residence_time:
                count_bound_time = 0

                for i in range(0,len(H4_contacts[str(key)])):
                ## we combined trajectories from five indenpent run thus they are not continuous in time.
                    if i == len(H4_contacts[str(key)])-1 or i == int(len(H4_contacts[str(key)])/2) :
                        H4_residence_time[key].append(count_bound_time*ts)
                        count_bound_time = 0
            ## tail continuously bound with DNA in frame i
                    if H4_contacts[str(key)].iloc[i] >= cut_off and i <= len(H4_contacts[str(key)])-2:
                        count_bound_time += 1 
                        #tail unbinding in next frame i+1
                        if H4_contacts[str(key)].iloc[i + 1] < cut_off:
                            #save the residence time
                            H4_residence_time[key].append(count_bound_time*ts)
                    # tail unbind in frame i
                    elif H4_contacts[str(key)].iloc[i] < cut_off and i <= len(H4_contacts[str(key)])-2:
                        # restart counting binding time
                        count_bound_time = 0

                            
            
            # calculate the mean residence time for each run.
            cut_off = 740 ##  we ignore the cases that tail residue stays in bound state for entire simulation time
            
            for key in H2A_residence_mean:
                if len(sorted(x for x in H2A_residence_time[key] if x <= cut_off)) >=1:
                   
                    H2A_residence_mean[key].append(statistics.mean(sorted(x for x in H2A_residence_time[key] if x <= cut_off)))
                
            for key in H2B_residence_mean:
                if len(sorted(x for x in H2B_residence_time[key] if x <= cut_off)) >=1:
                    H2B_residence_mean[key].append(statistics.mean(sorted(x for x in H2B_residence_time[key] if x <= cut_off)))
            
            for key in H3_residence_mean:
                if len(sorted(x for x in H3_residence_time[key] if x <= cut_off)) >=1:
                    H3_residence_mean[key].append(statistics.mean(sorted(x for x in H3_residence_time[key] if x <= cut_off)))
                
            for key in H4_residence_mean:
                if len(sorted(x for x in H4_residence_time[key] if x <= cut_off)) >=1:
                    H4_residence_mean[key].append(statistics.mean(sorted(x for x in H4_residence_time[key] if x <= cut_off)))




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



w = csv.writer(open("h2a_tail_residue_mean_residence_time.csv", "w"))
for mean, stdev in zip(H2A_residence_mean_all_run.items(),H2A_residence_stdev_all_run.items()):
    w.writerow([mean[0], mean[1], stdev[1]])
    
w = csv.writer(open("h2b_tail_residue_mean_residence_time.csv", "w"))
for mean, stdev in zip(H2B_residence_mean_all_run.items(),H2B_residence_stdev_all_run.items()):
    w.writerow([mean[0], mean[1], stdev[1]])
    
w = csv.writer(open("h3_tail_residue_mean_residence_time.csv", "w"))
for mean, stdev in zip(H3_residence_mean_all_run.items(),H3_residence_stdev_all_run.items()):
    w.writerow([mean[0], mean[1], stdev[1]])
    
w = csv.writer(open("h4_tail_residue_mean_residence_time.csv", "w"))
for mean, stdev in zip(H4_residence_mean_all_run.items(),H4_residence_stdev_all_run.items()):
    w.writerow([mean[0], mean[1], stdev[1]])
    
            
            
            
                
                     
            
