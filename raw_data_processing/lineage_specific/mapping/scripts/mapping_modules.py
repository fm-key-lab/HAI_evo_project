## This file contains all necessary python functions for running the mapping snakemake
## Changes to any of these functions should only be done on the dev branch and thoroughly tested prior to merging

## dependencies: genomestats
from Bio import SeqIO
import glob
import numpy as np

## dependencies: pileup_to_diversity_matrix,vcf_to_quals (get_input_output_names,round_half_up,pileup_to_div_matrix_snakemake)
import os
import gzip
from scipy.stats import ttest_ind, fisher_exact,ttest_1samp
from math import log10, floor

## dependencies: read_move_link
from pathlib import Path
import sys
import subprocess

## functions
def genomestats(REFGENOMEFOLDER):
    # parse ref genome to extract relevant stats
    # accepts genome.fasta or genome.fasta.gz (gzip) in refgenomefolder
    fasta_file = glob.glob(REFGENOMEFOLDER + '/genome.fasta')
    if len(fasta_file) != 1:
        fasta_file_gz = glob.glob(REFGENOMEFOLDER + '/genome.fasta.gz')
        if len(fasta_file_gz) != 1:
            raise ValueError('Either no genome.fasta(.gz) or more than 1 genome.fasta(.gz) file found in ' + REFGENOMEFOLDER)
        else: # genome.fasta.gz
            refgenome = SeqIO.parse(gzip.open(fasta_file_gz[0], "rt"),'fasta')
    else: # genome.fasta
        refgenome = SeqIO.parse(fasta_file[0],'fasta')
    Genomelength = 0
    ChrStarts = []
    ScafNames = []
    for record in refgenome:
        ChrStarts.append(Genomelength) # chr1 starts at 0 in analysis.m
        Genomelength = Genomelength + len(record)
        ScafNames.append(record.id)
    # close file
    #refgenome.close() # biopy update SeqIO has no close attribute anymore.
    # turn to np.arrys!
    ChrStarts = np.asarray(ChrStarts,dtype=int)
    Genomelength = np.asarray(Genomelength,dtype=int)
    ScafNames = np.asarray(ScafNames,dtype=object)
    return [ChrStarts,Genomelength,ScafNames]

def get_input_output_names(sample_path_to_pileup,sample_path_to_diversity,sample_path_to_coverage):
    fname_in=os.getcwd() +  '/' + sample_path_to_pileup  # Aro_Change: new; print this
    # Path to diversity.mat (output)
    fname_out=os.getcwd() +  '/' + sample_path_to_diversity  # Aro_Change: new; print this
    # Path to coverage.mat (secondary output)
    fname_out_cov=os.getcwd() +  '/' + sample_path_to_coverage  # Aro_Change: new; print this
    #fname_in='strain.pileup'; # Aro_Change: old
    #fname_out='diversity.mat'; # Aro_Change: old
    return (fname_in,fname_out,fname_out_cov)

def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    return floor(n*multiplier + 0.5) / multiplier    

def parse_multi_genome_smpls(SAMPLE_ls,REF_Genome_ls):
    ## Expand lists if multiple genomes are used within a sample
    REF_Genome_ext_ls = []
    SAMPLE_ext_ls = []
    for sampleID, refgenomes in zip(SAMPLE_ls, REF_Genome_ls):
        for refgenome in refgenomes.split(" "):
            REF_Genome_ext_ls.append(refgenome)
            SAMPLE_ext_ls.append(sampleID)
    return [REF_Genome_ext_ls, SAMPLE_ext_ls]

def get_non_outgroup_bams_for_freebayes(SAMPLE_ls, REF_Genome_ls, outgroup_ls):
    ## note: multiple ref genomes + varied ingroup/outgroup identity might fuck this up if a sample is ingrp for one ref, outgrp for the other
    non_outgroup_sample_ls = {}
    for sampleID,refgenomes,outgroup_bool in zip(SAMPLE_ls,REF_Genome_ls,outgroup_ls):
        for refgenome in refgenomes.split(" "):
            if int(outgroup_bool) == 0:
                if refgenome not in non_outgroup_sample_ls:
                    non_outgroup_sample_ls[refgenome] = [sampleID]
                else: 
                    non_outgroup_sample_ls[refgenome].append(sampleID)
    return non_outgroup_sample_ls

def pileup_to_div_matrix_snakemake(sample_path_to_pileup,sample_path_to_diversity,sample_path_to_coverage,ref_genome_directory, min_reads_on_strand=20): ## default value in tami's code: 
    # Declare constants
    Phred_offset=33; #mpileup is told that the reads are in fastq format and corrects that so its always at Phred+33 when the mpileup comes out
    nts=[ord(nt) for nt in 'ATCGatcg']
    indels=[ord(i) for i in '+-']
    numfields=40; #CHANGED 2019.01
    numtests=50; #how many tests to do in parallel ##NOTE: CURRENTLY NOT IMPLEMENTED
    indelregion=3
    ## get input and output filenames to parse
    (fname_in,fname_out,fname_out_cov)=get_input_output_names(sample_path_to_pileup,sample_path_to_diversity,sample_path_to_coverage)
    # Get reference genome information
    #load alignment_info % Aro_Change: old
    #refgenome=ae.Genome; % Aro_Change: old    
    print(f'\n Reference genome: {ref_genome_directory} \n')
    [ChrStarts,GenomeLength,ScafNames] = genomestats(ref_genome_directory) # Aro_Change: new
    print(f'\n {GenomeLength} \n') # Aro_Change: print this
    #[ChrStarts,GenomeLength,ChromosomeIndicator,ScafNames]=genomestats(['/scratch/mit_lieberman/Reference_Genomes/'
    #refgenome]); %RUN_ON_CLUSTER % Aro_Change: old
    # %initialize
    ## maybe grab Genome length from give genomestats object?
    data=np.zeros((numfields,GenomeLength),dtype=int)
    line_dict={}
    with open(fname_in,'rt') as fid:
        index=0
        for line in fid:
            temp=np.zeros((numfields,1),dtype=int)
            if type(line) == str:
                parsed=line.strip().split('\t') ## split into list and get rid of any EOL 
                chrom=parsed[0]
                if len(ChrStarts) == 1:
                    position=int(parsed[1])-1 ## python is 0 based, not 1 based as in matlab
                else:
                    if chrom not in ScafNames:
                        raise Exception("scaffold name in pileup not found in reference")
                    position=int(ChrStarts[np.where(ScafNames == chrom)] + int(parsed[1])) - 1 ## first part finds start position of this chrom and adds loc of SNP on chrom // MF: changed from .index to np.where and changed int(line[1]) to int(parsed[1])
                if ord(parsed[2]) in nts: ## check first to make sure ref is ATCG or atcg, not N, Y, (other ambiguous calls)
                    ref=nts.index(ord(parsed[2])) ## get numerical value of reference base
                    if ref > 3:
                        ref = ref-4
                else: 
                    ref = None ##may want to be more permissive for parsing these, currently all base calls will not be counted for this position
                calls_ascii_values=np.array([ord(i) for i in parsed[4]])
                bq_ascii_values=np.array([ord(i) for i in parsed[5]]) #base quality, BAQ corrected
                mq_ascii_values=np.array([ord(i) for i in parsed[6]]) #mapping quality
                ## distance from tail, needs loop due to weird encoding of certain chars
                if parsed[7] == "*":
                    p_7=[0]
                else: 
                    p_7=parsed[7].split(',')
                td_ascii_values=np.array([int(i) for i in p_7]) #distance from tail 

                #starts of reads
                start_of_read_vector=np.where(calls_ascii_values==94)[0] ## 94 == '^' == start of read segment ## needs to find all positions
                for i in start_of_read_vector: #(do this to the ^ char and the next)
                    calls_ascii_values[i]=-1 ## convert start of read to -1
                    calls_ascii_values[i+1]=-1 #remove mapping character (mapping qual= ASCII in index after ^ - 33),
                    #absolutely required because the next chracter could be $
                
                #read ends ##TODO: CHECK THIS WITH A LINE WITH END READS
                end_of_read_vector=np.where(calls_ascii_values==36)[0] ## 36 == '$' == end of read segment 
                temp[37]=len(start_of_read_vector) + len(end_of_read_vector)
                calls_ascii_values[end_of_read_vector]=-1 ## convert end of reads to -1

                ## indels and calls from reads supporting indels 
                ## TODO: check indel vector where these appear
                indel_vector=np.concatenate([np.where(calls_ascii_values==43)[0], np.where(calls_ascii_values==45)[0]]) ## all positions that called an indel
                for i in indel_vector:
                    # find size 
                    forward_looking_index_for_indel_size = 1
                    while calls_ascii_values[i + forward_looking_index_for_indel_size] >=48 and calls_ascii_values[i + forward_looking_index_for_indel_size] <58:
                        forward_looking_index_for_indel_size+=1
                    if forward_looking_index_for_indel_size > 1:
                        indel_ascii=list(calls_ascii_values[i+1:i+forward_looking_index_for_indel_size])
                        indel_size=int(''.join(map(chr,indel_ascii)))
                        indeld=forward_looking_index_for_indel_size-1
                    else:
                        indel_size=calls_ascii_values[i+1]-48
                        indeld=1
                    # record that indel was found in nearby region
                    #store directly into data, as it affects lines earlier and later
                    if calls_ascii_values[i]==45: ## deletion
                        if (position - indelregion + 1) > 0 and (position+indel_size+indelregion) <= GenomeLength:
                            data[38:40,position-indelregion+1:position+indel_size+indelregion + 1 ]=data[38:40,position-indelregion+1:position+indel_size+indelregion+1]+1
                        elif (position - indelregion+1) > 0:
                            data[38:40,position-indelregion+1:]=data[38:40,position-indelregion+1:]+1
                        else:
                            data[38:40,0:position+indelregion+indel_size+1]=data[38:40,0:position+indel_size+indelregion+1]+1
                    else: ## insertion 
                        if (position-indelregion+1) > 0 and (position+indelregion) <= GenomeLength:
                            data[38,position-indelregion+1:position+indelregion+1]=data[38,position-indelregion+1:position+indelregion+1]+1
                        elif (position-indelregion+1)>0:
                            data[38,position-indelregion+1:]=data[38,position-indelregion+1:]+1
                        else:
                            data[38,0:position+indel_size+1]=data[38,0:position+indel_size+1]+1
                    #remove these calls from counting
                    calls_ascii_values[i:(i+indeld+indel_size+1)]=-1 # dont remove base that preceeds and indel
                #   ignore '*', as these deletion markers won't be counted towards score
                #   %qualities, etc and the presence of a deletion is recorded by the
                #   upstream '-' indicator 
                ## replace reference with actual calls
                if ref != None:
                    calls_ascii_values[np.where(calls_ascii_values == 46)[0]] = nts[ref] ## '.' forward, store in nts index values ATCGatcg==0:7
                    calls_ascii_values[np.where(calls_ascii_values == 44)[0]] = nts[ref+4] ## ',' reverse, store in nts index values ATCGatcg==0:7         
                ## indx reads for finding scores ## TODO might be bugged
                simple_calls=calls_ascii_values[np.where(calls_ascii_values>0)]## indices of all non-filtered calls (see above)
                for i in range(len(nts)):
                    current_nt=nts[i]
                    current_nt_indices=np.where(simple_calls==current_nt)[0]
                    size_for_this_nt=len(current_nt_indices) ##TODO CHECK IF RIGHT WITH OTHER FILES
                    if size_for_this_nt != 0:
                        temp[i]=size_for_this_nt
                        temp[i+8]=round_half_up(sum(bq_ascii_values[(current_nt_indices)])/size_for_this_nt) - Phred_offset ## average phred score for calls of given base
                        temp[i+16]=round_half_up(sum(mq_ascii_values[(current_nt_indices)])/size_for_this_nt) - 33 ## average mapping qual for this given base
                        temp[i+24]=round_half_up(sum(td_ascii_values[(current_nt_indices)])/size_for_this_nt) ## average tail distance
                ## The following section is critical for metagenomic samples 
                # find major and nextmajor allele
                allele_index_summed_calls=[(x,temp[x]+temp[x+4]) for x in range(4)] # add T and t calls, C and c calls, etc. keep index of major allele as first index in tuple
                allele_index_summed_calls.sort(key=lambda x: x[1][0]) ## sort by number of calls for a given base
                n1 = allele_index_summed_calls[-1][0] # major allele (most calls) (returns actual index of base in nts)
                n2 = allele_index_summed_calls[-2][0] # next most abundant allele (returns actual index of base in nts)
                
                x=np.logical_or(simple_calls==nts[n1], simple_calls==nts[n1+4]) # boolean array of reads with major alelle
                y=np.logical_or(simple_calls==nts[n2], simple_calls==nts[n2+4]) # boolean array of reads with minor alelle
                
                if sum(temp[0:4]) > min_reads_on_strand and sum(temp[4:8])> min_reads_on_strand and sum(y)*200>sum(x): # only calcualte p values of there are greater than 20 reads on each strand and MAF < .995
                    bttests_x=bq_ascii_values[x]-Phred_offset
                    bttests_y=bq_ascii_values[y]-Phred_offset # gets base qualities for major/minor alleles
                    mttests_x=mq_ascii_values[x]-Phred_offset
                    mttests_y=mq_ascii_values[y]-Phred_offset # gets mapping qualities for major/minor alleles
                    fttests_x=td_ascii_values[(np.where(simple_calls==nts[n1])[0])]
                    fttests_y=td_ascii_values[(np.where(simple_calls==nts[n2])[0])] # gets tail distances from fwd reads for major/minor alleles
                    rttests_x=td_ascii_values[(np.where(simple_calls==nts[n1+4])[0])]
                    rttests_y=td_ascii_values[(np.where(simple_calls==nts[n2+4])[0])] # gets tail distances from rev reads for major/minor alleles
                    
                    ## NOTE: matlab ttest will output result if len(y) == 1, treating it as the hypothesized mean of the normal distribution compared to x 
                    bp=ttest_ind(bttests_x,bttests_y)[1]
                    mp=ttest_ind(mttests_x,mttests_y)[1]
                    fp=ttest_ind(fttests_x,fttests_y)[1]
                    rp=ttest_ind(rttests_x,rttests_y)[1]
                    # report pval of above t-tests as -log10
                        ## coverting values of 0 and nan for output as -1:
                    if bp == 0 or np.isnan(bp):
                        temp[33]=-1
                    else: temp[33]=round_half_up(-log10(bp)) 

                    if mp == 0 or np.isnan(mp):
                        temp[34]=-1
                    else: temp[34]=round_half_up(-log10(mp))

                    if fp == 0 or np.isnan(fp):
                        temp[35]=-1
                    else: temp[35]=round_half_up(-log10(fp))
                    
                    if rp == 0 or np.isnan(rp):
                        temp[36]=-1
                    else: temp[36]=round_half_up(-log10(rp))
                    
                    ## currently not broken but testing against matlab script fails since matlab script *is* broken.
                    p=fisher_exact(np.array([[temp[n1][0], temp[n2][0]],[temp[n1+4][0],temp[n2+4][0]]]), alternative='two-sided')[1] # fisher exact test for strand bias (contingency table = # of reads that are fwd/rev and support major/minor allele)
                    if p == 0:
                        temp[32]=-1
                    else:
                        temp[32]=round_half_up(-log10(p))
                ## store the data
                data[0:38,position]=np.ravel(temp[0:38]) 
                ## NOTE: column 38 should possibly be overwritten here, if following matlab script
                ## however, this seems weird and like a bug in their code.
                ## if we want to use indels, we should figure this out thoroughly
            index+=1
    coverage=sum(data[0:8,:])
    np.savez_compressed(fname_out, data = data.astype(int))
    np.savez_compressed(fname_out_cov, coverage = coverage.astype(int))

def vcf_to_quals_snakemake(sample_path_to_vcf, sample_path_to_quals, REF_GENOME_DIRECTORY):

    ## Input file
    fname_in_gz = [os.getcwd() + '/' + sample_path_to_vcf] ## store zipped file names

    [ChrStarts,GenomeLength,ScafNames] = genomestats(REF_GENOME_DIRECTORY);

    print(f'\n {GenomeLength} \n') ## print to stdout or file??

    quals = np.zeros(GenomeLength, dtype=int)

    i = 1
    position = 0

    for file in fname_in_gz:
        with gzip.open(file, 'rt') as fid:

            for line in fid:
                if line[0] != "#":

                    ## print evry 50000 lines a dot (check for running status)
                    if not i%50000:
                        print('.', end = '')

                    ## parsing line
                    l = line.split('\t')

                    ## get single column values
                    chrom = np.where(ScafNames == l[0])
                    position_on_chr = int(l[1])
                    position = ChrStarts[chrom] + position_on_chr

                    ## store nucleotides
                    alt = l[4]
                    ref = l[3]

                    # only store quality for simple calls (not indel, not ambigious)
                    if (alt != '') & (',' not in alt) & (len(alt) == len(ref)) & (len(ref) == 1):
                        ## find and parse quality score in INFO column
                        info_col = l[7].split(';')

                        ## find FQ entry in info column and store value after '=' sign
                        entrywithFQ = [item for item in info_col if item.startswith('FQ')]

                        ## check if FQ was found
                        if not entrywithFQ:
                            print('No FQ entry was found in' + file + ' in row ' + str(i) + '!')
                        else:
                            ## Extract FQ value from string
                            entrywithFQ_str = entrywithFQ[0]
                            fq = float(entrywithFQ_str[entrywithFQ_str.index('=')+1:])

                            ## if already a position with a stronger FQ here, don't include this. (More negative is stronger)
                            if fq < quals[position-1]: ##Note python is 0 based!
                                quals[position-1] = int(round_half_up(fq))

                    ## Needed for progess bar on stdout
                    i += 1

    ## save file
    outpath = os.getcwd() + '/' + sample_path_to_quals
    np.savez_compressed(outpath, quals=quals)

    ## print to stdout
    print(f'\n Saved: {outpath} \n')

def read_samplesCSV(spls):
    # reads in samples.csv file, format: Batch, Sample,Alignments,ProviderName,Patient
    hdr_check = ['Batch', 'Sample', 'Alignments', 'ProviderName', 'Patient', 'Outgroup']
    switch = "on"
    file = open(spls, 'r')
    list_path = []
    list_splID = []
    list_providerNames = []
    list_refG = []
    list_patient = []
    list_outgroup = [] ## 2022.05.05 MF: added outgroup parsing to function
    for line in file:
        line = line.strip('\n').split(',')
        # Test Header. Note: Even when header wrong code continues (w/ warning), but first line not read.
        if switch == "on":
            if (line == hdr_check):
                print("Passed CSV header check")
            else:
                Warning("CSV did NOT pass header check! Code continues, but first line ignored")
            switch = "off"
            continue
        # build lists
        list_path.append(line[0])
        list_splID.append(line[1])
        list_refG.append(line[2])
        list_providerNames.append(line[3])
        list_patient.append(line[4])
        list_outgroup.append(line[5])
    return [list_path,list_splID,list_refG,list_providerNames,list_patient,list_outgroup] # 2019.02.08, Arolyn: removed set(list_patient) provides only unique subject IDs

def split_samplesCSV(PATH_ls,SAMPLE_ls,REF_Genome_ls,PROVIDER_ls,PATIENTID_ls_all,OUTGROUP_ls):
    # Added by Arolyn, 2019.02.11
    # takes info extracted form samples.csv; saves each line of samples.csv as sample_info.csv in data/{sampleID}
    for i, sample in enumerate(SAMPLE_ls):
        # get info for this sample
        sample_info_csv_text = PATH_ls[i] + ',' + SAMPLE_ls[i] + ',' + REF_Genome_ls[i] + ',' + PROVIDER_ls[i] + ',' + PATIENTID_ls_all[i] + ',' + OUTGROUP_ls[i]
        #print( sample )
        #print( sample_info_csv_text )
        # make data directory for this sample if it doesn't already exist
        if not(os.path.isdir('data/' + sample)):
            os.makedirs('data/' + sample, exist_ok=True)
        # check to see if this mini csv with sample info already exists
        if os.path.isfile('data/' + sample + '/sample_info.csv'):
            # if so, read file
            old_file = open('data/' + sample + '/sample_info.csv','r')
            old_info = old_file.readline()
            old_file.close()
            # check to see if the existing file is consistent with samples.csv
            if not(old_info == sample_info_csv_text):
                # if not, remove the old file and save sample info in a new file
                #print('Information file must be updated.')  
                os.remove('data/' + sample + '/sample_info.csv')
                f = open('data/' + sample + '/sample_info.csv','w')
                f.write(sample_info_csv_text) 
                f.close()
            #else:
            #print('Information file already updated.')              
        else: # if mini csv with sample info does not already exist
            # save sample info in mini csv
            #print('Information file must be created.')  
            f = open('data/' + sample + '/sample_info.csv','w')
            f.write(sample_info_csv_text) 
            f.close()

def findfastqfile(dr,ID,filename):
    fwd=[]
    rev=[]
    potentialhits_forward=glob.glob(dr + '/' + filename +'/*1.fastq.gz')
    potentialhits_reverse=glob.glob(dr + '/' + filename +'/*2.fastq.gz')
    if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
        fwd=potentialhits_forward[0]
        rev=potentialhits_reverse[0]
    elif len(potentialhits_forward)==0 or len(potentialhits_reverse)==0: ## need or statement to further screen path if just one file was found i.e. for *R1_001.fastq.gz *R2_001.fastq.gz
        potentialhits_forward=glob.glob(dr + '/' + filename +'/*1_001.fastq.gz')
        potentialhits_reverse=glob.glob(dr + '/' + filename +'/*2_001.fastq.gz')
        if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
            fwd=potentialhits_forward[0]
            rev=potentialhits_reverse[0]
        elif len(potentialhits_forward)==0 or len(potentialhits_reverse)==0:
            potentialhits_forward=glob.glob(dr + '/' + filename +'*1.fastq.gz')
            potentialhits_reverse=glob.glob(dr + '/' + filename +'*2.fastq.gz')
            if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                fwd=potentialhits_forward[0]
                rev=potentialhits_reverse[0]
            elif len(potentialhits_forward)==0 or len(potentialhits_reverse)==0:
                potentialhits_forward=glob.glob(dr + '/' + filename +'*1_001.fastq.gz')
                potentialhits_reverse=glob.glob(dr + '/' + filename +'*2_001.fastq.gz')
                if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                    fwd=potentialhits_forward[0]
                    rev=potentialhits_reverse[0]
                elif len(potentialhits_forward)==0 or len(potentialhits_reverse)==0:
                    potentialhits_forward=glob.glob(dr + '/' + filename +'*1.fastq')
                    potentialhits_reverse=glob.glob(dr + '/' + filename +'*2.fastq')
                    if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                        subprocess.run("gzip " + potentialhits_forward[0], shell=True)
                        subprocess.run("gzip " + potentialhits_reverse[0], shell=True)
                        fwd=potentialhits_forward[0]+'.gz'
                        rev=potentialhits_reverse[0]+'.gz'
                    elif len(potentialhits_forward)==0 or len(potentialhits_reverse)==0:
                        potentialhits_forward=glob.glob(dr + '/' + filename +'*1_001.fastq')
                        potentialhits_reverse=glob.glob(dr + '/' + filename +'*2_001.fastq')
                        if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                            subprocess.run("gzip " + potentialhits_forward[0], shell=True)
                            subprocess.run("gzip " + potentialhits_reverse[0], shell=True)
                            fwd=potentialhits_forward[0]+'.gz'
                            rev=potentialhits_reverse[0]+'.gz'
                    else:
                        foldername=glob.glob(dr + '/' + filename + '*')
                        if foldername and os.path.isdir(foldername[0]):
                            foldername=foldername[0]
                            potentialhits_forward=glob.glob(foldername + '/*' + filename + '*1*.fastq.gz')
                            potentialhits_reverse=glob.glob(foldername + '/*' + filename + '*2*.fastq.gz')
                            if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                                fwd=potentialhits_forward[0]
                                rev=potentialhits_reverse[0]
                            elif len(potentialhits_forward)==0 or len(potentialhits_reverse)==0:
                                print(foldername + '/*' + filename + '*2*.fastq.gz')
                                potentialhits_forward=glob.glob(foldername +  '/*' + filename + '*1*.fastq')
                                potentialhits_reverse=glob.glob(foldername + '/*' + filename + '*2*.fastq')
                                if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                                    subprocess.run("gzip " + potentialhits_forward[0], shell=True)
                                    subprocess.run("gzip " + potentialhits_reverse[0], shell=True)
                                    fwd=potentialhits_forward[0]+'.gz'
                                    rev=potentialhits_reverse[0]+'.gz'
    if not(fwd) or not(rev):
        raise ValueError('Either no file or more than 1 file found in ' + dr + ' for ' + ID)
    ##zip fastq files if they aren't already zipped
    subprocess.run("gzip " + fwd, shell=True)   
    subprocess.run("gzip " + rev, shell=True)   
    return [fwd, rev]

def makelink(path,sample,providername):
    #When sample is run on a single lane
    #Provider name can be either a COMPLETE directory name or a file name in batch(called path in this fx)
    [fwd_file, rev_file]=findfastqfile(path,sample, providername)
    subprocess.run('ln -s -T ' + fwd_file + ' data/' + sample + '/R1.fq.gz', shell=True)    
    subprocess.run('ln -s -T ' + rev_file + ' data/' + sample + '/R2.fq.gz', shell=True)    
        
def cp_append_files(paths,sample,providernames):
    #When sample is run on multiple lanes with same barcode
    fwd_list=''
    rev_list=''
    for path, providername in zip(paths, providernames):
        #Provider name can be either a COMPLETE directory name or a file name in batch(called path in this fx)
        [fwd_file, rev_file]=findfastqfile(path,sample, providername)
        fwd_list=fwd_list+ ' ' +fwd_file
        rev_list=rev_list+ ' ' +rev_file
        print(rev_list)
        print(fwd_list)
    subprocess.run("zcat " + fwd_list + ' | gzip > data/' +  sample + '/R1.fq.gz', shell=True)
    subprocess.run("zcat " + rev_list + ' | gzip > data/' +  sample + '/R2.fq.gz', shell=True)
