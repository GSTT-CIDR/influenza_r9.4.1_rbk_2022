#!/usr/bin/env python
# coding: utf-8

from datetime import datetime, timedelta
import pandas as pd
import os
import glob
import pyfastx
import pysam
import argparse

# Create the parser
folder_parser = argparse.ArgumentParser()

# Add the arguments
folder_parser.add_argument( 'run_name', metavar='run_name', type=str, help='Input run name to process')

# Execute the parse_args() method
args = folder_parser.parse_args()

#set run name for the run
run_name = args.run_name

#obtain file path for all fastq files
path_fastq_all = glob.glob(f"/mnt/flu/data/{run_name}/*/*/fastq_pass/", recursive = False)
fastq = glob.glob(f"/mnt/flu/data/{run_name}/*/*/fastq_pass/", recursive = False)[0]

#combine all fastq files into the same directory (this is to avoid issues if experiment stopped and then restarted)
if len(path_fastq_all) > 1:
    for path in path_fastq_all[1:]:
        os.system(f' rsync -rv --ignore-existing {path}* {fastq}')

#generate analysis time; currently gives number of days after sequencing started that analysis performed
analysis_time = int(datetime.today().strftime('%Y%m%d')[2:]) - int(run_name[:6])
        
#import input sheet; run_name.csv
#ANON_ID, SAMPLE_ID, barcode, sample_date, RUN
input_sheet = pd.read_csv(str(f"/mnt/flu/runsheets/{run_name}.csv"))

#convert 3 to 03 etc

for line in input_sheet.index:
    if input_sheet.barcode[line] <10:
        input_sheet.barcode[line] = str("0"+str(input_sheet.barcode[line]))

#create run_id and analysis_id from
input_sheet['RUN_ID'] = (input_sheet.RUN.astype(str)) + (input_sheet.barcode.astype(str))
input_sheet['ANALYSIS_ID'] = input_sheet.RUN_ID + "_" + str(analysis_time)

#count total reads for each barcode
#count number of reads mapping to lambda for each barcode as measure of barcode crosstalk

#make a list of used barcodes in the run and dictionaries to link them with analysis_id
sample_dict = dict(zip(input_sheet.barcode, input_sheet.ANALYSIS_ID))
lambda_list = ["NC_001416.1"]

#build dataframe to insert control processing data
reads_table = pd.DataFrame(columns=list(input_sheet.ANALYSIS_ID), index = lambda_list)
reads_table.iloc[:] = 0

#import references to map reads to
alignment_reference = ('/mnt/flu/control_processing/lamda.fa')

#align and count reads to lambda for each sample
#create dictionaries to store read count for each barcode
total_reads_dict = dict(zip(input_sheet.barcode, [0]*len(sample_dict.keys())))

#create folder to store output if not already made by a previous analysis
if os.path.exists(f"/mnt/flu/control_processing/{run_name}/") == 0:
    os.system(f"mkdir /mnt/flu/control_processing/{run_name}/")

if os.path.exists(f"/mnt/flu/control_processing/{run_name}/{analysis_time}/") == 0:
    os.system(f"mkdir /mnt/flu/control_processing/{run_name}/{analysis_time}/")

#iterate through each barcode in run
for barcode in sample_dict.keys():
    
    #make a new dictionary for each sample
    lamda_dict = dict(zip(lambda_list,[0]*len(lambda_list)))
    #provide path to fastq files for that barcode
    path = fastq + "barcode" + str(barcode) + "/"
        
    #align with minimap2 within singularity
    os.system((f" cat {path}*.fastq.gz | singularity run --bind /mnt/flu/control_processing/:/mnt/flu/control_processing/ /home/tgsw/repositories/singularity/staphb_minimap2.sif  minimap2 -ax map-ont {alignment_reference} - > /mnt/flu/control_processing/{run_name}/{analysis_time}/{barcode}_aligned.sam"))
    
    #sort and index with pysam (samtools)
    pysam.sort('-o', f"/mnt/flu/control_processing/{run_name}/{analysis_time}/{barcode}_aligned_sorted.bam", f"/mnt/flu/control_processing/{run_name}/{analysis_time}/{barcode}_aligned.sam")
    pysam.index(f"/mnt/flu/control_processing/{run_name}/{analysis_time}/{barcode}_aligned_sorted.bam")
    
    #use samtools to count number of mapped reads
    samfile = pysam.AlignmentFile(f"/mnt/flu/control_processing/{run_name}/{analysis_time}/{barcode}_aligned_sorted.bam", "rb")
    
    #count number of reads mapped to lambda
    for read in samfile:
        total_reads_dict[barcode] += 1    
        if read.reference_name in lamda_dict.keys():
            lamda_dict[read.reference_name] += 1         
    #set lamda read counts for that barcode
    reads_table[sample_dict.get(barcode)] = lamda_dict.values()

#set unmapped reads and total reads for all barcodes
for barcode in sample_dict.keys():
    reads_table.loc['Total reads',sample_dict.get(barcode)] = total_reads_dict.get(barcode)

#calculate proportion of lambda reads
reads_table.loc['Lamda reads as proportion of total reads'] = reads_table.apply(lambda x: x['NC_001416.1']/x['Total reads']*100 if x['Total reads'] !=0 else 0)

#output reads_table
reads_table.to_csv(f"/mnt/flu/control_processing/{run_name}/{analysis_time}/control_processing.csv")

#create sample sheet in csv format
sample_sheet = input_sheet[["ANALYSIS_ID", "barcode"]]
sample_sheet = sample_sheet.rename(columns={"ANALYSIS_ID":"sample_id"})
sample_sheet.barcode = sample_sheet.barcode.apply(lambda x: f"barcode{x}")

#write sample sheet
sample_sheet.to_csv(f"/mnt/flu/runsheets/{run_name}_irma.csv")

#dictionary to convert from barcode to sample id
analysis_id_to_sample_id_dict = dict(zip((input_sheet.barcode.astype(str)), input_sheet.SAMPLE_ID))

#set variables for irma command
outdir_higher_folder = str (f"/mnt/flu/irma_output/{run_name}/")
outdir = str (f"/mnt/flu/irma_output/{run_name}/{analysis_time}/")
if os.path.exists(outdir_higher_folder)==False:
    os.mkdir(outdir_higher_folder)
if os.path.exists(outdir)==False:
    os.mkdir(outdir)

#main consensus pipeline
#iterates through each barcode in the run
#filters by read length, trims adaptor/barcode, removes human reads, runs IRMA, then calculates coverage and depth for each segment
segment_dict = {
    1: 'PB2',
    2: 'PB1',
    3: 'PA',
    4: 'HA',
    5: 'NP',
    6: 'NA' ,
    7: 'MP',
    8: 'NS' }

segment_len_dict = {
    1: 2280,
    2: 2274,
    3: 2151,
    4: 1701,
    5: 1497,
    6: 1410 ,
    7: 982,
    8: 838 }

def count_n(text, seq):
    return text.count(seq)

#create coverage df
df_coverage = pd.DataFrame(columns=["barcode","segment_name","coverage","H_type"])

#run pipeline for each barcode
for barcode in sample_sheet.barcode:

    #combine raw fastq into combined input file
    os.system(f"cat {fastq}{barcode}/FA*fastq.gz > {fastq}{barcode}/{barcode}.fastq.gz" )
    os.system(f"gunzip -f {fastq}{barcode}/{barcode}.fastq.gz") 
    
    #trim RBK adapter + barcode (len 120) with nanofilt, remove short reads
    os.system(f"cat {fastq}{barcode}/{barcode}.fastq | singularity run --bind /mnt/flu:/mnt/flu /mnt/flu/nanofilt NanoFilt -l 500 --headcrop 120 > {fastq}{barcode}/{barcode}.trim.fastq")
    
    #remove long headers
    fastq_list = []
    for name,seq,qual,comment in pyfastx.Fastx(f"{fastq}{barcode}/{barcode}.trim.fastq"):
        raw = f"@{name}\n{seq}\n+\n{qual}\n"
        fastq_list.append(raw)
            
    with open(f"{fastq}{barcode}/{barcode}.trim.fastq", "w") as out:
        for fq in fastq_list:
            out.write(fq)
    
    #remove human reads
    #first map to human genome
    os.system(f"singularity run --bind /mnt/flu:/mnt/flu /home/tgsw/repositories/singularity/staphb_minimap2.sif minimap2 -ax map-ont --secondary=no /mnt/flu/hg38_new.mmi {fastq}{barcode}/{barcode}.trim.fastq > {fastq}{barcode}/{barcode}.human.trim.sam")
    
    #then remove reads which map to human genome
    os.system(f"singularity run --bind /mnt/flu:/mnt/flu /mnt/flu/samtools.sif samtools fastq -f4 {fastq}{barcode}/{barcode}.human.trim.sam > {fastq}{barcode}/{barcode}.human_removed.trim.fastq")
    
    #run IRMA on each barcode
    #using local config file (FLU-minion-final)
    os.system(f"/home/tgsw/repositories/irma/flu-amd/IRMA FLU-minion-final {fastq}{barcode}/{barcode}.human_removed.trim.fastq {outdir}{barcode}")
    
    #get coverage for each segment
    for segment in range(1,9):
        segment_name = segment_dict.get(segment)
        
        #get number of Ns, skipping segments with no consensus files
        if len(glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{barcode}/amended_consensus/{barcode}_{segment}*fa"))==0:
            continue
        consensus_file = glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{barcode}/amended_consensus/{barcode}_{segment}*fa")[0]
        
        for name, seq, comment in pyfastx.Fastx (consensus_file):
            n_count = (count_n(seq, 'N'))
        
        #skip segments with no bam files
        if len(glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{barcode}/*_{segment_name}*bam"))==0:
            continue
        
        bam_file = glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{barcode}/*_{segment_name}*bam")[0]

        flu_type = bam_file.split('/')[-1][0]
        
        #get depth at each position
        pysam.sort("-o", f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{barcode}/{flu_type}_{segment_name}.sorted.bam", bam_file)
        read_count = pysam.view("-c", "-F", "260", f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{barcode}/{flu_type}_{segment_name}.sorted.bam")
        depth_df = pd.DataFrame([x.split('\t') for x in pysam.depth(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{barcode}/{flu_type}_{segment_name}.sorted.bam").split('\n')])

        #remove last line which is empty
        depth_df = depth_df[~depth_df[2].isnull()]
        #convert to integer
        depth_df[2] =  depth_df[2].astype(int)
        #count depth at least 20
        cov_20 = (len(depth_df[depth_df[2] > 19]))
        length = len(depth_df)
        coverage = cov_20/length
        
        #mean coverage
        mean_cov = depth_df[2].sum()/length
        
        #add column with sample ID to make downstream analysis easier
        SAMPLE_ID = str(analysis_id_to_sample_id_dict.get(str(barcode[7:9])))
        
        #add H type to coverage df, which is important for later analysis
        #add unknown if no HA for that sample
        if len(glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{barcode}/A_HA_*bam"))==0:
            if len(glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{barcode}/B_HA_*bam"))==1:
                H_type = 'B'
            else:
                H_type = 'Unknown'
            df_row = pd.DataFrame([{"barcode":barcode,
                        "SAMPLE_ID":SAMPLE_ID,
                        "segment_name": segment_name,
                        "coverage":round(coverage*100,1),
                        "mean_coverage" : mean_cov,
                        "number_reads" : read_count[:-1],
                        "length":length,
                        "n_count":n_count,
                        "H_type":H_type,
                        "flu_type":flu_type}])
            df_coverage = pd.concat([df_coverage, df_row], axis=0, ignore_index=True)
            continue
        
        #get H type is there is HA for that sample
        H_type = glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{barcode}/A_HA_*bam")[0]
        H_type = H_type[-6:-4]
        
        df_row = pd.DataFrame([{"barcode":barcode,
                                "SAMPLE_ID":SAMPLE_ID,
                                "segment_name": segment_name,
                                "coverage":round(coverage*100,1),
                                "mean_coverage" : mean_cov,
                                "number_reads" : read_count[:-1],
                                "length":length,
                                "n_count":n_count,
                                "H_type":H_type,
                                "flu_type":flu_type}])
        #add to df
        df_coverage = pd.concat([df_coverage, df_row], axis=0, ignore_index=True)

#write out coverage df
df_coverage.to_csv(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/segment_coverage.csv")

#get all consensus segments above 90% for each sample in the run and place in a file for further analysis
#rename each segment to SAMPLE_ID_SEGMENT
analysis_id_to_sample_id_dict = dict(zip((input_sheet.barcode.astype(str)), input_sheet.SAMPLE_ID))

#clear any previous files
#seperate H1, H3, B files
for segment in range(1,9):
    if os.path.exists(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H1_all_consensus_{segment}.fasta"):
        os.remove(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H1_all_consensus_{segment}.fasta")
    
    if os.path.exists(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H3_all_consensus_{segment}.fasta"):
        os.remove(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H3_all_consensus_{segment}.fasta")
    
    if os.path.exists(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/B_all_consensus_{segment}.fasta"):
        os.remove(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/B_all_consensus_{segment}.fasta")
    
    path_to_segment_consensus = []
    
    #select only segments with > 90% coverage and < 10% ambiguous bases to use for downstream analysis
    for barcode in sample_sheet.barcode:
    #skip samples with no bam for the segment
        if len(df_coverage.coverage[(df_coverage.barcode == barcode) & (df_coverage.segment_name == segment_dict.get(segment))])==0:
            continue
        
        #add path to segment if coverage > 90% at depth x20
        if df_coverage.coverage[(df_coverage.barcode == barcode) & (df_coverage.segment_name == segment_dict.get(segment))].values[0] > 90:
            #and no more than 10% length is Ns
            if df_coverage.n_count[(df_coverage.barcode == barcode) & (df_coverage.segment_name == segment_dict.get(segment))].values[0] < segment_len_dict.get(segment)/10:
                path_to_segment_consensus.append(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{barcode}/amended_consensus/{barcode}_{segment}.fa")
                df_coverage.loc[(df_coverage.barcode == barcode) & (df_coverage.segment_name == segment_dict.get(segment)),'Include']=True

    #combine segments with 90% coverage from the run
    #place into seperate multifasta for different H1/H3/B as these will be aligned seperately
    for sequence in path_to_segment_consensus:
        barcode = sequence.split('/')[6][:9]
        H_type = df_coverage.H_type[(df_coverage.barcode==barcode)].values[0]
        #add consensus sequences for segment to a multifasta
        os.system(f"cat {sequence} >> /mnt/flu/irma_output/{run_name}/{analysis_time}/{H_type}_all_consensus_{segment}.fasta")
    
    #rename segments from barcode to SAMPLE_ID_SEGMENT
    if len(glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H1_all_consensus_{segment}.fasta"))!=0:
        fasta_list = []
        for name, seq, comment in pyfastx.Fastx (f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H1_all_consensus_{segment}.fasta"):
            sample_name = str(analysis_id_to_sample_id_dict.get(str(name[7:9]))) + f"_{segment}" + str(name[11:])
            raw = f">{sample_name} \n{seq}\n"
            fasta_list.append(raw)

        with open(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H1_all_consensus_{segment}.fasta", "w") as out:
            for fa in fasta_list:
                out.write(fa)
    
    if len(glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H3_all_consensus_{segment}.fasta"))!=0:
        fasta_list = []
        for name, seq, comment in pyfastx.Fastx (f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H3_all_consensus_{segment}.fasta"):
            sample_name = str(analysis_id_to_sample_id_dict.get((name[7:9]))) + f"_{segment}" + str(name[11:])
            raw = f">{sample_name} \n{seq}\n"
            fasta_list.append(raw)

        with open(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H3_all_consensus_{segment}.fasta", "w") as out:
            for fa in fasta_list:
                out.write(fa)
    
    if len(glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/B_all_consensus_{segment}.fasta"))!=0:
        fasta_list = []
        for name, seq, comment in pyfastx.Fastx (f"/mnt/flu/irma_output/{run_name}/{analysis_time}/B_all_consensus_{segment}.fasta"):
            sample_name = str(analysis_id_to_sample_id_dict.get((name[7:9]))) + f"_{segment}" + str(name[11:])
            raw = f">{sample_name} \n{seq}\n"
            fasta_list.append(raw)

        with open(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/B_all_consensus_{segment}.fasta", "w") as out:
            for fa in fasta_list:
                out.write(fa)

#Make lists of WGS and majority genome samples, and export to df; split by H1 and H3
df_coverage_90_H1 = df_coverage[(df_coverage.Include == True )&(df_coverage.H_type == 'H1')]
df_coverage_90_H3 = df_coverage[(df_coverage.Include == True)&(df_coverage.H_type == 'H3')]
df_coverage_90_B = df_coverage[(df_coverage.Include == True)&(df_coverage.H_type == 'B')]

df_coverage_90_MG_H1 = df_coverage_90_H1[df_coverage_90_H1.segment_name.isin(['HA','NP','NA' ,'MP','NS','PA','PB2'])]
df_coverage_90_MG_H3 = df_coverage_90_H3[df_coverage_90_H3.segment_name.isin(['HA','NP','NA' ,'MP','NS','PA','PB2'])]
df_coverage_90_MG_B = df_coverage_90_B[df_coverage_90_B.segment_name.isin(['HA','NP','NA' ,'MP','NS','PA','PB2'])]

wgs_samples_H1 = []
wgs_samples_H3 = []
wgs_samples_B = []
MG_samples_H1 = []
MG_samples_H3 = []
MG_samples_B = []

#h1
for barcode in sample_sheet.barcode:
    if len(df_coverage_90_H1[df_coverage_90_H1.barcode==barcode])==8:
        wgs_samples_H1.append(str(barcode[-2:]))

for barcode in sample_sheet.barcode:
    if len(df_coverage_90_MG_H1[df_coverage_90_MG_H1.barcode==barcode])==7:
        MG_samples_H1.append(str(barcode[-2:]))
#h3
for barcode in sample_sheet.barcode:
    if len(df_coverage_90_H3[df_coverage_90_H3.barcode==barcode])==8:
        wgs_samples_H3.append(str(barcode[-2:]))

for barcode in sample_sheet.barcode:
    if len(df_coverage_90_MG_H3[df_coverage_90_MG_H3.barcode==barcode])==7:
        MG_samples_H3.append(str(barcode[-2:]))
        
#b
for barcode in sample_sheet.barcode:
    if len(df_coverage_90_B[df_coverage_90_B.barcode==barcode])==8:
        wgs_samples_B.append(str(barcode[-2:]))

for barcode in sample_sheet.barcode:
    if len(df_coverage_90_MG_B[df_coverage_90_MG_B.barcode==barcode])==7:
        MG_samples_B.append(str(barcode[-2:]))

#convert to sample ID
wgs_samples_H1_sample = []
for item in wgs_samples_H1:
    wgs_samples_H1_sample.append(analysis_id_to_sample_id_dict.get(item))

wgs_samples_H3_sample = []
for item in wgs_samples_H3:
    wgs_samples_H3_sample.append(analysis_id_to_sample_id_dict.get(item))

wgs_samples_B_sample = []
for item in wgs_samples_B:
    wgs_samples_B_sample.append(analysis_id_to_sample_id_dict.get(item))

MG_samples_H1_sample = []
for item in MG_samples_H1:
    MG_samples_H1_sample.append(analysis_id_to_sample_id_dict.get(item))
                                 
MG_samples_H3_sample = []
for item in MG_samples_H3:
    MG_samples_H3_sample.append(analysis_id_to_sample_id_dict.get(item))

MG_samples_B_sample = []
for item in MG_samples_B:
    MG_samples_B_sample.append(analysis_id_to_sample_id_dict.get(item))
    
#output csv file containing sample names
wgs_samples_H1_df = pd.DataFrame (wgs_samples_H1_sample, columns = ['sample'])
MG_samples_H1_df = pd.DataFrame (MG_samples_H1_sample, columns = ['sample'])

MG_samples_H1_df.to_csv(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H1_MG_samples.csv")

wgs_samples_H3_df = pd.DataFrame (wgs_samples_H3_sample, columns = ['sample'])
MG_samples_H3_df = pd.DataFrame (MG_samples_H3_sample, columns = ['sample'])

MG_samples_H3_df.to_csv(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H3_MG_samples.csv")

wgs_samples_B_df = pd.DataFrame (wgs_samples_B_sample, columns = ['sample'])
MG_samples_B_df = pd.DataFrame (MG_samples_B_sample, columns = ['sample'])

MG_samples_B_df.to_csv(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/B_MG_samples.csv")

#align N2 against a reference set of sequences and mask homoploymer pre-resistance calling to prevent frameshift
#run mafft
os.system(f"cat /mnt/flu/irma_output/{run_name}/{analysis_time}/H3_all_consensus_6.fasta /mnt/flu/irma_output/GISAID_BACKGROUND/h3n2/na/sequences.fasta | mafft-linsi - > /mnt/flu/irma_output/{run_name}/{analysis_time}/H3_all_consensus_6.aligned.fasta" )

#NA
#for each sequence, take first 130 bases in alignment, then mask position 131 - 137 with N, then add remaining sequence from base 138
masked_H3_NA = []
for name, seq, comment in pyfastx.Fastx (f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H3_all_consensus_6.aligned.fasta"):
    masked_sequence = seq[:130] + 'NNNNNNN' + seq[137:]
    raw = f">{name}\n{masked_sequence}\n"
    masked_H3_NA.append(raw)

#overwrite fasta
with open(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H3_all_consensus_6.aligned.fasta", "w") as out:
    for fa in masked_H3_NA:
        out.write(fa)

#align NA and PA segments with nextalign to try and predict drug resistance
#done seperately for H1 and H3

segment_dict = {
    1: 'pb2',
    2: 'pb1',
    3: 'pa',
    4: 'ha',
    5: 'np',
    6: 'na' ,
    7: 'ma',
    8: 'ns' }

ha_na_type_dict = {
    'H1':'h1n1pdm',
    'H3':'h3n2',
    'B':'vic'}

#nextalign create protein sequences for N2
for flu_type in ['H3']:
    segment_name = segment_dict.get(6)
    ha_na_type = ha_na_type_dict.get(flu_type)
    
    if len(glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{flu_type}_all_consensus_6.aligned.fasta"))!=0:
    
        os.system(f"singularity run --bind /mnt/flu:/mnt/flu /mnt/flu/nextalign.sif nextalign run --input-ref=/mnt/flu/nextclade-master/data/flu/{ha_na_type}/{segment_name}/reference.fasta --genemap=/mnt/flu/nextclade-master/data/flu/{ha_na_type}/{segment_name}/genemap.gff --output-all=/mnt/flu/irma_output/{run_name}/{analysis_time}/{flu_type}_6  /mnt/flu/irma_output/{run_name}/{analysis_time}/{flu_type}_all_consensus_6.aligned.fasta")

#nextalign create protein sequences for N1 and IBV NA
for flu_type in ['H1','B']:
    segment_name = segment_dict.get(6)
    ha_na_type = ha_na_type_dict.get(flu_type)
    
    if len(glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{flu_type}_all_consensus_6.fasta"))!=0:

        os.system(f"singularity run --bind /mnt/flu:/mnt/flu /mnt/flu/nextalign.sif nextalign run --input-ref=/mnt/flu/nextclade-master/data/flu/{ha_na_type}/{segment_name}/reference.fasta --genemap=/mnt/flu/nextclade-master/data/flu/{ha_na_type}/{segment_name}/genemap.gff --output-all=/mnt/flu/irma_output/{run_name}/{analysis_time}/{flu_type}_6 /mnt/flu/irma_output/{run_name}/{analysis_time}/{flu_type}_all_consensus_6.fasta")

#nextalign create protein sequences for PA
for flu_type in ['H1','H3','B']:
    segment_name = segment_dict.get(3)
    ha_na_type = ha_na_type_dict.get(flu_type)

    if len(glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{flu_type}_all_consensus_3.fasta"))!=0:

        os.system(f"singularity run --bind /mnt/flu:/mnt/flu /mnt/flu/nextalign.sif nextalign run --input-ref=/mnt/flu/nextclade-master/data/flu/{ha_na_type}/{segment_name}/reference.fasta --genemap=/mnt/flu/nextclade-master/data/flu/{ha_na_type}/{segment_name}/genemap.gff --output-all=/mnt/flu/irma_output/{run_name}/{analysis_time}/{flu_type}_3 /mnt/flu/irma_output/{run_name}/{analysis_time}/{flu_type}_all_consensus_3.fasta")

#nextclade for HA gene for clade
for flu_type in ['H1','H3','B']:
    segment=4
    segment_name = segment_dict.get(4)
    ha_na_type = ha_na_type_dict.get(flu_type)
        
    if len(glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{flu_type}_all_consensus_4.fasta"))!=0:

        os.system(f"singularity run --bind /mnt/flu:/mnt/flu /mnt/flu/nextclade.sif nextclade run --input-dataset=/mnt/flu/nextclade-data/{flu_type} --output-all=/mnt/flu/irma_output/{run_name}/{analysis_time}/{flu_type}_{segment} /mnt/flu/irma_output/{run_name}/{analysis_time}/{flu_type}_all_consensus_{segment}.fasta")

#import resistance df
res_df = pd.read_excel('/mnt/flu/nextclade-master/data/flu/IAV_RES_DICT_LOCAL.xlsx',sheet_name='SNP',keep_default_na=False)
res_df['Position'] = res_df.apply(lambda x: int(x[0][1:-1]), axis=1)
res_df['Alternative'] = res_df.apply(lambda x: (x[0][1:]), axis=1)
res_df['Reference'] = res_df.apply(lambda x: (x[0][:-1]), axis=1)

reference_sample_df = pd.DataFrame(columns=['name','flu_type','segment_name','Amino Acid'])
resistance_sample_df = pd.DataFrame(columns=['name','flu_type','segment_name','Amino Acid'])
non_reference_non_resistance_df = pd.DataFrame(columns=['name','flu_type','segment_name','Amino Acid'])

#predict resistance (at consensus level)
for flu_type in ['H1','H3','B']:
    for segment in [3,6]:
        
        if len(glob.glob(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{flu_type}_all_consensus_{segment}.fasta"))!=0:
            
            segment_name = (segment_dict.get(segment)).upper()
            temp_res_df = res_df[(res_df['Type']==flu_type)&(res_df['Gene']==segment_name)]
            
            #for each translated protein sequence, look at each position, and print output to a dataframe of resistant and 
            for name, seq, comment in pyfastx.Fastx (f"/mnt/flu/irma_output/{run_name}/{analysis_time}/{flu_type}_{segment}/nextalign_gene_{segment_name}.translation.fasta"):
                for position in pd.unique(temp_res_df['Position']):
                    if (str(position)+seq[position-1]) in temp_res_df.Alternative.values:
                        res_row = pd.DataFrame([{"name":name,
                                                 "flu_type":flu_type,
                                                 "segment_name":segment_name,
                                                 "Amino Acid":(str(position)+seq[position-1])}])
                        resistance_sample_df = pd.concat([resistance_sample_df, res_row], ignore_index=True)

if len(wgs_samples_H3_df)!=0:

    clade_df = pd.read_csv(f'/mnt/flu/irma_output/{run_name}/{analysis_time}/H3_4/nextclade.csv',sep=';')

    #metadata file
    wgs_samples_H3_df = wgs_samples_H3_df.set_index('sample')
    wgs_samples_H3_df['PA Resistance']='No'
    wgs_samples_H3_df['NA Resistance']='No'

    for line in wgs_samples_H3_df.index:
        wgs_samples_H3_df.loc[line,'ANON_ID'] = line.split('.')[0]
        wgs_samples_H3_df.loc[line,'SAMPLE_ID'] = line
        wgs_samples_H3_df.loc[line,'Location'] =  str(input_sheet.Location[input_sheet.SAMPLE_ID==line].values[0])
        wgs_samples_H3_df.loc[line,'barcode'] = str(input_sheet.barcode[input_sheet.SAMPLE_ID==line].values[0])
        wgs_samples_H3_df.loc[line,'sample_date'] = input_sheet.sample_date[input_sheet.SAMPLE_ID==line].values[0]
        wgs_samples_H3_df.loc[line,'RUN'] = input_sheet.RUN[input_sheet.SAMPLE_ID==line].values[0]
        wgs_samples_H3_df.loc[line,'RUN_ID'] = input_sheet.RUN_ID[input_sheet.SAMPLE_ID==line].values[0]
        wgs_samples_H3_df.loc[line,'ANALYSIS_ID'] = input_sheet.ANALYSIS_ID[input_sheet.SAMPLE_ID==line].values[0]
        wgs_samples_H3_df.loc[line,'Clade'] = clade_df.clade[clade_df.seqName==(line+'_4')].values[0]
        wgs_samples_H3_df.loc[line,'Coverage'] = (13133 - df_coverage.n_count[df_coverage.barcode == 'barcode' + str(input_sheet.barcode[input_sheet.SAMPLE_ID==line].values[0])].sum())/13133*100
        if (line + '_3') in list(resistance_sample_df.name[resistance_sample_df.segment_name == 'PA']):
            wgs_samples_H3_df.loc[line,'PA Resistance'] = resistance_sample_df['Amino Acid'][(resistance_sample_df.segment_name == 'NA')&(resistance_sample_df.name == (line + '_3'))].values[0]
        if (line + '_6') in list(resistance_sample_df.name[resistance_sample_df.segment_name == 'NA']):
            wgs_samples_H3_df.loc[line,'NA Resistance'] = resistance_sample_df['Amino Acid'][(resistance_sample_df.segment_name == 'NA')&(resistance_sample_df.name == (line + '_6'))].values[0]

    wgs_samples_H3_df = wgs_samples_H3_df.round(2)

if len(wgs_samples_H1_df)!=0:

    clade_df = pd.read_csv(f'/mnt/flu/irma_output/{run_name}/{analysis_time}/H1_4/nextclade.csv',sep=';')

    #metadata file
    wgs_samples_H1_df = wgs_samples_H1_df.set_index('sample')
    wgs_samples_H1_df['PA Resistance']='No'
    wgs_samples_H1_df['NA Resistance']='No'

    for line in wgs_samples_H1_df.index:
        wgs_samples_H1_df.loc[line,'ANON_ID'] = line.split('.')[0]
        wgs_samples_H1_df.loc[line,'SAMPLE_ID'] = line
        wgs_samples_H1_df.loc[line,'Location'] =  str(input_sheet.Location[input_sheet.SAMPLE_ID==line].values[0])
        wgs_samples_H1_df.loc[line,'barcode'] = str(input_sheet.barcode[input_sheet.SAMPLE_ID==line].values[0])
        wgs_samples_H1_df.loc[line,'sample_date'] = input_sheet.sample_date[input_sheet.SAMPLE_ID==line].values[0]
        wgs_samples_H1_df.loc[line,'RUN'] = input_sheet.RUN[input_sheet.SAMPLE_ID==line].values[0]
        wgs_samples_H1_df.loc[line,'RUN_ID'] = input_sheet.RUN_ID[input_sheet.SAMPLE_ID==line].values[0]
        wgs_samples_H1_df.loc[line,'ANALYSIS_ID'] = input_sheet.ANALYSIS_ID[input_sheet.SAMPLE_ID==line].values[0]
        wgs_samples_H1_df.loc[line,'Clade'] = clade_df.clade[clade_df.seqName==(line+'_4')].values[0]
        wgs_samples_H1_df.loc[line,'Coverage'] = (13134 - df_coverage.n_count[df_coverage.barcode == 'barcode' + str(input_sheet.barcode[input_sheet.SAMPLE_ID==line].values[0])].sum())/13134*100
        if (line + '_3') in list(resistance_sample_df.name[resistance_sample_df.segment_name == 'PA']):
            wgs_samples_H1_df.loc[line,'PA Resistance'] = resistance_sample_df['Amino Acid'][(resistance_sample_df.segment_name == 'NA')&(resistance_sample_df.name == (line + '_3'))].values[0]
        if (line + '_6') in list(resistance_sample_df.name[resistance_sample_df.segment_name == 'NA']):
            wgs_samples_H1_df.loc[line,'NA Resistance'] = resistance_sample_df['Amino Acid'][(resistance_sample_df.segment_name == 'NA')&(resistance_sample_df.name == (line + '_6'))].values[0]

    wgs_samples_H1_df = wgs_samples_H1_df.round(2)

if len(wgs_samples_B_df)!=0:

    clade_df = pd.read_csv(f'/mnt/flu/irma_output/{run_name}/{analysis_time}/B_4/nextclade.csv',sep=';')

    #metadata file
    wgs_samples_B_df = wgs_samples_B_df.set_index('sample')
    wgs_samples_B_df['PA Resistance']='No'
    wgs_samples_B_df['NA Resistance']='No'

    for line in wgs_samples_B_df.index:
        wgs_samples_B_df.loc[line,'ANON_ID'] = line.split('.')[0]
        wgs_samples_B_df.loc[line,'SAMPLE_ID'] = line
        wgs_samples_B_df.loc[line,'Location'] =  str(input_sheet.Location[input_sheet.SAMPLE_ID==line].values[0])
        wgs_samples_B_df.loc[line,'barcode'] = str(input_sheet.barcode[input_sheet.SAMPLE_ID==line].values[0])
        wgs_samples_B_df.loc[line,'sample_date'] = input_sheet.sample_date[input_sheet.SAMPLE_ID==line].values[0]
        wgs_samples_B_df.loc[line,'RUN'] = input_sheet.RUN[input_sheet.SAMPLE_ID==line].values[0]
        wgs_samples_B_df.loc[line,'RUN_ID'] = input_sheet.RUN_ID[input_sheet.SAMPLE_ID==line].values[0]
        wgs_samples_B_df.loc[line,'ANALYSIS_ID'] = input_sheet.ANALYSIS_ID[input_sheet.SAMPLE_ID==line].values[0]
        wgs_samples_B_df.loc[line,'Clade'] = clade_df.clade[clade_df.seqName==(line+'_4')].values[0]
        wgs_samples_B_df.loc[line,'Coverage'] = (13696 - df_coverage.n_count[df_coverage.barcode == 'barcode' + str(input_sheet.barcode[input_sheet.SAMPLE_ID==line].values[0])].sum())/13696*100
        if (line + '_3') in list(resistance_sample_df.name[resistance_sample_df.segment_name == 'PA']):
            wgs_samples_B_df.loc[line,'PA Resistance'] = resistance_sample_df['Amino Acid'][(resistance_sample_df.segment_name == 'NA')&(resistance_sample_df.name == (line + '_3'))].values[0]
        if (line + '_6') in list(resistance_sample_df.name[resistance_sample_df.segment_name == 'NA']):
            wgs_samples_B_df.loc[line,'NA Resistance'] = resistance_sample_df['Amino Acid'][(resistance_sample_df.segment_name == 'NA')&(resistance_sample_df.name == (line + '_6'))].values[0]

    wgs_samples_B_df = wgs_samples_B_df.round(2)

wgs_samples_H1_df.to_csv(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H1_wgs_samples.csv")
wgs_samples_H3_df.to_csv(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/H3_wgs_samples.csv")
wgs_samples_B_df.to_csv(f"/mnt/flu/irma_output/{run_name}/{analysis_time}/B_wgs_samples.csv")
