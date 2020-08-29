#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 13:48:24 2020

@author: shanedenecke
"""

###Import ideas
import os
import pandas as pd


###Change directory
os.chdir('/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline/')


filtered=pd.read_csv('./Final_outputs/combined_files/Full_counts_long.tsv',sep='\t')
full_taxid=pd.read_csv('/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline/GENERAL_REFERENCE/CAFE/tree_table_full.tsv',sep='\t')



tree_table_filt=full_taxid[full_taxid.Species_name.isin(filtered.Species_name.values)]


for i in set(tree_table_filt.Tree):
    sub=tree_table_filt[tree_table_filt.Tree==i]
    taxids=[str(x)+'_0' for x in sub.taxid_code]
    
    with open('./CAFE/taxid_lists/'+i+'_taxid_codes.txt','w') as f:
        for t in taxids:
            f.write(t+'\n')
    
    