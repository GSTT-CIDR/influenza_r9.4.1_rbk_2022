#!/usr/bin/env python
# coding: utf-8

#IAV WGS epidemiology analysis from consensus sequences

#Input = list of run names to be used in /mnt/flu/irma_output (note takes most recent analysis)

#1) Parse a GISAID download of vaccine strains into usable format
#2) Get list of samples with complete coverage of all segments
#3) mafft align each segment
#4) concatenate sequences for WGS
#5) run snp-dists for WGS
#6) run iqtree for WGS

import pandas as pd
import glob
import os
from datetime import date
import pyfastx

runs_to_use = ['230130FLUTW','230114FLUTW3','230114FLUTW2','230114FLUTW1','221228FLUTW','221222FLUTW40','230103FLUTW','230110FLUTW1','230110FLUTW2','230110FLUTW3','VACCINES']

#Parse background data; format is a GISAID download of sequences with all 8 segments in (cds only), 
#Download stored in /mnt/flu/irma_output/GISAID_BACKGROUND/0/ 
segment_dict_caps = {
    1: 'PB2',
    2: 'PB1',
    3: 'PA',
    4: 'HA',
    5: 'NP',
    6: 'NA' ,
    7: 'MP',
    8: 'NS' }

segment_dict = {
    1: 'pb2',
    2: 'pb1',
    3: 'pa',
    4: 'ha',
    5: 'np',
    6: 'na' ,
    7: 'ma',
    8: 'ns' }

for segment in range(1,9):
    H3_list = []
    H1_list = []
    fasta_list_h1 = []
    fasta_list_h3 = []
    
    #parse background data set
    #this contains a Nothern Hemisphere vaccine strains (modified to contain only coding regions) 
    list_gisaid_files = ['A_Victoria_2570_2019','A_Wisconsin_588_2019','A_Darwin_6_2021','A_Darwin_9_2021']
    
    for file in list_gisaid_files:
        for name, seq, comment in pyfastx.Fastx (f"/mnt/flu/irma_output/GISAID_BACKGROUND/0/{file}.fasta"):

            #make list of sequences of each H_type
            #make files with segement sequence in same format as for local data
            #rename sequences in output fasta file

            if (name.split('|')[2][-4:-2]) == 'H1':
                H1_list.append(name.split('|')[0])
                if (name.split('|')[5]) == segment_dict_caps.get(segment):
                    raw = f">{(name.split('|')[0]) + f'_{segment}'} {comment}\n{seq}\n"
                    fasta_list_h1.append(raw)

            if (name.split('|')[2][-4:-2]) == 'H3':
                H3_list.append(name.split('|')[0])
                if (name.split('|')[5]) == segment_dict_caps.get(segment):
                    raw = f">{(name.split('|')[0])+ f'_{segment}'} {comment}\n{seq}\n"
                    fasta_list_h3.append(raw)

    with open(f"/mnt/flu/irma_output/VACCINES/0/H1_all_consensus_{segment}.fasta", "w") as out:
        for fa in fasta_list_h1:
            out.write(fa)
    
    with open(f"/mnt/flu/irma_output/VACCINES/0/H3_all_consensus_{segment}.fasta", "w") as out:
        for fa in fasta_list_h3:
            out.write(fa)

H1_df = pd.DataFrame (pd.unique(H1_list), columns = ['sample'])
H3_df = pd.DataFrame (pd.unique(H3_list), columns = ['sample'])

H1_df.to_csv(f"/mnt/flu/irma_output/VACCINES/0/H1_wgs_samples.csv")
H3_df.to_csv(f"/mnt/flu/irma_output/VACCINES/0/H3_wgs_samples.csv")

#get list of all samples with > 90% coverage of all 8 segments
#H1
path_to_wgs_samples_H1 = []
for run_name in runs_to_use:
    #pick final analysis for each run
    path_to_wgs_samples_H1.append(glob.glob(f"/mnt/flu/irma_output/{run_name}/*/H1_wgs_samples.csv")[-1])

wgs_all_H1 = pd.DataFrame (columns = ['sample'])

for file in path_to_wgs_samples_H1:
    wgs_from_run_H1 = pd.read_csv(file)
    wgs_all_H1 = pd.concat([wgs_all_H1, wgs_from_run_H1], axis=0, ignore_index=True)

#H3
path_to_wgs_samples_H3 = []
for run_name in runs_to_use:
    #pick final analysis for each run
    path_to_wgs_samples_H3.append(glob.glob(f"/mnt/flu/irma_output/{run_name}/*/H3_wgs_samples.csv")[-1])

wgs_all_H3 = pd.DataFrame (columns = ['sample'])

for file in path_to_wgs_samples_H3:
    wgs_from_run_H3 = pd.read_csv(file)
    wgs_all_H3 = pd.concat([wgs_all_H3, wgs_from_run_H3], axis=0, ignore_index=True)

#align all segments individually with mafft
#align H1 samples
for segment in range(1,9):
    path_to_segment_consensus = []
    for run_name in runs_to_use:
        #pick final analysis for each run
        if len(glob.glob(f"/mnt/flu/irma_output/{run_name}/*/H1_all_consensus_{segment}.fasta"))!= 0:
            path_to_segment_consensus.append(glob.glob(f"/mnt/flu/irma_output/{run_name}/*/H1_all_consensus_{segment}.fasta")[-1])
    
    #remove existing files
    if os.path.exists(f"/mnt/flu/gstt_epi/H1_combined_{segment}.fasta"):
        os.remove(f"/mnt/flu/gstt_epi/H1_combined_{segment}.fasta")
    
    for sequence in path_to_segment_consensus:
        #add consensus sequences for segment to a multifasta
        os.system(f"cat {sequence} >> /mnt/flu/gstt_epi/H1_combined_{segment}.fasta")
        
    #run mafft
    os.system(f"mafft-linsi /mnt/flu/gstt_epi/H1_combined_{segment}.fasta > /mnt/flu/gstt_epi/H1_combined_{segment}.aligned.fasta" )
    
#align H3 samples
for segment in range(1,9):
    path_to_segment_consensus = []
    for run_name in runs_to_use:
        #pick final analysis for each run
        if len(glob.glob(f"/mnt/flu/irma_output/{run_name}/*/H3_all_consensus_{segment}.fasta"))!= 0:
            path_to_segment_consensus.append(glob.glob(f"/mnt/flu/irma_output/{run_name}/*/H3_all_consensus_{segment}.fasta")[-1])
    
    #remove existing files
    if os.path.exists(f"/mnt/flu/gstt_epi/H3_combined_{segment}.fasta"):
        os.remove(f"/mnt/flu/gstt_epi/H3_combined_{segment}.fasta")
    
    for sequence in path_to_segment_consensus:
        #add consensus sequences for segment to a multifasta
        os.system(f"cat {sequence} >> /mnt/flu/gstt_epi/H3_combined_{segment}.fasta")
                
    #run mafft
    os.system(f"mafft-linsi /mnt/flu/gstt_epi/H3_combined_{segment}.fasta > /mnt/flu/gstt_epi/H3_combined_{segment}.aligned.fasta" )

#add mask for N2 segment
#for each sequence, take first 130 bases in alignment, then mask position 131 - 137 with N, then add remaining sequence from base 138 onwards
masked_H3_NA = []
for name, seq, comment in pyfastx.Fastx (f"/mnt/flu/gstt_epi/H3_combined_6.aligned.fasta"):
    masked_sequence = seq[:130] + 'NNNNNNN' + seq[137:]
    raw = f">{name}\n{masked_sequence}\n"
    masked_H3_NA.append(raw)

#overwrite fasta
with open(f"/mnt/flu/gstt_epi/H3_combined_6.aligned.fasta", "w") as out:
    for fa in masked_H3_NA:
        out.write(fa)

#h1
#combine aligned segments for WGS samples
fasta_list = []
for sample in pd.unique(wgs_all_H1['sample']):
    sample_sequence_list = []
    for segment in range(1,9):
        for name, seq, comment in pyfastx.Fastx (f"/mnt/flu/gstt_epi/H1_combined_{segment}.aligned.fasta"):
            if name == str(sample) + '_' + str(segment):
                sample_sequence_list.append(seq)
    #combine segments into a single string             
    final_seq = "".join(sample_sequence_list)
    #write fasta entry to list
    raw = f">{sample}\n{final_seq}\n"
    fasta_list.append(raw)
        
#write fasta
with open(f"/mnt/flu/gstt_epi/H1_combined_wg.aligned.fasta", "w") as out:
    for fa in fasta_list:
        out.write(fa)

#h3
fasta_list = []
for sample in pd.unique(wgs_all_H3['sample']):
    sample_sequence_list = []
    for segment in range(1,9):
        for name, seq, comment in pyfastx.Fastx (f"/mnt/flu/gstt_epi/H3_combined_{segment}.aligned.fasta"):
            if name == str(sample) + '_' + str(segment):
                sample_sequence_list.append(seq)
    #combine segments into a single string             
    final_seq = "".join(sample_sequence_list)
    #write fasta entry to list
    raw = f">{sample}\n{final_seq}\n"
    fasta_list.append(raw)

#write fasta
with open(f"/mnt/flu/gstt_epi/H3_combined_wg.aligned.fasta", "w") as out:
    for fa in fasta_list:
        out.write(fa)

#snp-dists on H1 whole genomes
os.system('/home/themisch/snp-dists/snp-dists -c /mnt/flu/gstt_epi/H1_combined_wg.aligned.fasta > /mnt/flu/gstt_epi/H1_snp_dists_wg.csv')

#snp-dists on H3 whole genomes
os.system('/home/themisch/snp-dists/snp-dists -c /mnt/flu/gstt_epi/H3_combined_wg.aligned.fasta > /mnt/flu/gstt_epi/H3_snp_dists_wg.csv')

#iqtree on H1 whole genomes
os.system(f'iqtree -s /mnt/flu/gstt_epi/H1_combined_wg.aligned.fasta -m TEST -bb 1000')

#iqtree on H3 whole genomes
os.system(f'iqtree -s /mnt/flu/gstt_epi/H3_combined_wg.aligned.fasta -m TEST -bb 1000')
