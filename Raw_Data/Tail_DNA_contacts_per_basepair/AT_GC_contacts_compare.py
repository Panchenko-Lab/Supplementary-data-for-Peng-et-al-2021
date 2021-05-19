#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 12:30:42 2021

@author: pengy10
"""
import pandas as pd 
import numpy as np 
import re
import csv
import statistics 
import seaborn as sns
import matplotlib.pyplot as plt

df_AT_GC_conntacts = pd.read_csv("/Users/pengy10/Desktop/DNA_Methylation/MD/MD_analysis_all/amber/DNA_contacts_4500ns/contact_mean_LYS.csv", \
                                 sep = "\t", header=0)
df_AT_GC_conntacts.rename({'x': 'mean_contacts'}, axis=1, inplace=True)

df_AT_GC_conntacts ["residue"] = list(range(-93,94))

df_AT_GC_conntacts ["AT"] = 0

AT_list = [-84, -76, -73, -58,-57,-56, -50, -46, -42, -32, -28, -26, -24, -23, -7, -5, \
           -1, 1, 7, 11, 14, 20, 22, 24, 27, 28, 32, 37, 40, 45, 48, 56, 64, 69, 73, 75, 76, 77]

for i in range(0,187):
    if(df_AT_GC_conntacts["residue"][i] in AT_list ):
        df_AT_GC_conntacts ["AT"][i] = 1
        

bplot2=sns.boxplot(y='mean_contacts', x='AT', 
                 data=df_AT_GC_conntacts, 
                 width=0.4, boxprops=dict(alpha=0.9),
                 palette=["red","blue"])
# add stripplot to boxplot with Seaborn
bplot2=sns.stripplot(y='mean_contacts', x='AT', 
                   data=df_AT_GC_conntacts, 
                   jitter=True, 
                    marker='o', 
                   alpha=0.8,
                   color='black')
bplot2.set_xlabel("DNA nucleotide content",fontsize=13, weight = "bold")
bplot2.set_ylabel("Mean contacts number" ,fontsize=13, weight = "bold")
bplot2.set_xticklabels(["GC", "AT"], weight = "bold")
bplot2.set_yticklabels(np.asarray(bplot2.get_yticks(), dtype=int), weight = "bold")
bplot2.get_figure().savefig('AT_GC_compare_LYS.png',  dpi = 300 )
plt.clf()

# get ANOVA table as R like output
import statsmodels.api as sm
from statsmodels.formula.api import ols
# reshape the d dataframe suitable for statsmodels package 
#    d_melt = pd.melt(d.reset_index(), id_vars=['index'], value_vars=['A', 'B', 'C', 'D'])
# replace column names
#    d_melt.columns = ['index', 'treatments', 'value']
# Ordinary Least Squares (OLS) model
model = ols('mean_contacts ~ C(AT)', data=df_AT_GC_conntacts).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print(anova_table)

# load packages
from statsmodels.stats.multicomp import pairwise_tukeyhsd
# perform multiple pairwise comparison (Tukey HSD)
m_comp = pairwise_tukeyhsd(endog=df_AT_GC_conntacts['mean_contacts'], groups=df_AT_GC_conntacts['AT'], alpha=0.05)
print(m_comp)
    