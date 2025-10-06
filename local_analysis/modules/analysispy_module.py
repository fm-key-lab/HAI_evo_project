# module accpanying analysis.py
# Module supported and maintaned for analysis of candidate_mutation_table.pickle.gz (created by build_candidate_mutation_table.py; previous matlab-based table not recommended to use due 1-based/0-based conflict!)


## Version History
## 08.04.2022, MF: Added Multiprocesing for plot_coverage_fwd_rev_stacked()
## Initial version from MIT Lieberman Lab version 01.2022


import os
import subprocess
import re
import csv
import glob
import gzip
import collections
import warnings
import time,datetime
import sys

import pickle
import itertools
import numpy as np
import pandas as pd
from collections import OrderedDict, defaultdict

from scipy import stats, sparse
from statsmodels.stats.proportion import proportion_confint
import statsmodels.api as sm
import networkx

from Bio import SeqIO, SeqRecord
from Bio import Phylo
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor #NJ tree
from BCBio import GFF # might be redundant
from functools import partial # for gff_parse() and plot_coverage_fwd_rev_stacked()

from matplotlib import rc
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import concurrent.futures
import upsetplot
import seaborn as sns
import multiprocessing
from multiprocessing import Pool

def dnaparse2lol(file):
    # Each entry will be element in list. Each element in list is a list of length 2, with header [0] and sequence [1]
    # First header line is skipped
    ctr=0
    fa_lol = [] 
    for line in file:
        ctr += 1
        if ctr == 1:
            continue # skip header
        if ctr > 1:
            line = line.strip().split() # read in reference seq
            fa_lol.append( [line[0],line[1]] )
    return fa_lol

def read_samplesCSV(spls):
    # reads in samples.csv file, format: Path,Sample,ReferenceGenome,ProviderName,Subject
    hdr_check = ['Path', 'Sample', 'ReferenceGenome', 'ProviderName', 'Subject']
    switch = "on"
    file = open(spls, 'r')
    list_path = []
    list_splID = []
    list_providerNames = []
    list_refG = []
    list_patient = []
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
    return [list_path,list_splID,list_refG,list_providerNames,list_patient] # set(list_patient) provides only unique subject IDs

def read_candidate_mutation_table_pickle_gzip(file_cmt_pickle_gz, indels = False, isgzip = True):
    # read the candidate_mutation_table.pickle.gz files
    ## p: 0-based index of genomic position (means actual position is p+1)
    ## coverage_stats: each row is sample, col [0-10) == covg bins 1x, 2x... >10x; [10]==median covg; [11]==mean; [12]==stddev 
    ## INDEL ARRAYS:
    # indel_p: (1d, pos) - start position of variant called in joint indel vcf (0-based as genome wide positions)
    # indel_depth: (3d, pos x samples x (cov ref, cov alt (highest GL)) ) -- for multiallelic sites: select alt with highest genotype likelihood
    # indel_support: (3d, pos x samples x [GL ref, GL alt]) -- for multiallelic sites: select alt with highest genotype likelihood
    # indel_size_gp: - (3d, pos x samples x [len(alt)-len(ref)) -- for multiallelic sites: select alt with highest genotype likelihood
    # indel_call: (3d: pos x samples x [ref allele (str), alt allele with highest likelihood (str)) -- for multiallelic sites: select alt with highest genotype likelihood
    if isgzip:
        with gzip.open(file_cmt_pickle_gz, 'rb') as f:
            cmt = pickle.load(f)
            if indels:
                [sampleNames, p, counts, quals, in_outgroup, indel_counter, coverage_stats, indel_p, indel_call, indel_size_gp, indel_depth, indel_support] = [np.array(cit_mat) for cit_mat in cmt]
                return [quals,p,counts,in_outgroup,sampleNames,indel_counter,coverage_stats, ## arrays for SNV calls
                        indel_p,indel_call,indel_size_gp,indel_depth,indel_support] ## arrays for indel calls
            else:
                [sampleNames, p, counts, quals, in_outgroup, indel_counter, coverage_stats] = [np.array(cit_mat) for cit_mat in cmt[:7]]
                return [quals,p,counts,in_outgroup,sampleNames,indel_counter,coverage_stats]
    else:
        with open(file_cmt_pickle_gz, 'rb') as f:
            cmt = pickle.load(f)
            if indels:
                [sampleNames, p, counts, quals, in_outgroup, indel_counter, coverage_stats, indel_p, indel_call, indel_size_gp, indel_depth, indel_support] = [np.array(cit_mat) for cit_mat in cmt]
                return [quals,p,counts,in_outgroup,sampleNames,indel_counter,coverage_stats, ## arrays for SNV calls
                        indel_p,indel_call,indel_size_gp,indel_depth,indel_support] ## arrays for indel calls
            else:
                [sampleNames, p, counts, quals, in_outgroup, indel_counter, coverage_stats] = [np.array(cit_mat) for cit_mat in cmt[:7]]
                return [quals,p,counts,in_outgroup,sampleNames,indel_counter,coverage_stats]

        
def genomestats(REFGENOMEFOLDER):
    # parse ref genome to extract relevant stats
    # accepts genome.fasta or genome.fasta.gz (gzip) in refgenomefolder
    ref_genome_path=str(REFGENOMEFOLDER)+'/genome.fasta'
    fasta_file = glob.glob(ref_genome_path)
    if len(fasta_file) != 1:
        ref_genome_path=str(REFGENOMEFOLDER)+'/genome.fasta.gz'
        fasta_file_gz = glob.glob(ref_genome_path)
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

def p2chrpos(p, ChrStarts):
    '''# return 2col array with chr and pos on chr
    #p...continous, ignores chr
    #pos: like p, 0-based'''
        
    # get chr and pos-on-chr
    chr = np.ones(len(p),dtype=int)
    if len(ChrStarts) > 1:
        for i in ChrStarts[1:]:
            chr = chr + (p >= i) # when (p > i) evaluates 'true' lead to plus 1 in summation. > bcs ChrStarts start with 0...genomestats()
        positions = p - ChrStarts[chr-1] # [chr-1] -1 due to 0based index
        pos = np.column_stack((chr,positions))
    else:
        pos = np.column_stack((chr,p))
    return pos

def extract_outgroup_mutation_positions(REFGENOMEFOLDER,position,CMTpy=True):
    # extracts the ref nucleotide for every position. positions needs to be sorted by chr
    # reads genome.fasta or genome.fasta.gz    
    fasta_file = glob.glob(REFGENOMEFOLDER + '/genome.fasta')
    if len(fasta_file) != 1:
        fasta_file_gz = glob.glob(REFGENOMEFOLDER + '/genome.fasta.gz')
        if len(fasta_file_gz) != 1:
            raise ValueError('Either no genome.fasta(.gz) or more than 1 genome.fasta(.gz) file found in ' + REFGENOMEFOLDER)
        else:
            refgenome=fasta_file_gz[0]
    else:
        refgenome=fasta_file[0]
    refnt = np.zeros(position.shape[0],dtype=object)
    pos_counter=0
    for index,record in enumerate(SeqIO.parse(refgenome, format='fasta')):
        poschr = position[ position[:,0]==index+1 , 1]
        for sglpos in poschr:
            refnt[pos_counter]=record.seq[int(sglpos)]
            pos_counter += 1
    return refnt
    
def extract_seq_from_fasta(refgenomefolder,chrom,start,end,strand):
    # returns seq object between start and end
    # chrom is 0 based integer! start end 1-based (I think!)
    # NOTE: coordinates relative on chromosome, if more than 1 present in fasta
    # reads genome.fasta or genome.fasta.gz    
    fasta_file = glob.glob(refgenomefolder + '/genome.fasta')
    if len(fasta_file) != 1:
        fasta_file_gz = glob.glob(refgenomefolder + '/genome.fasta.gz')
        if len(fasta_file_gz) != 1:
            raise ValueError('Either no genome.fasta(.gz) or more than 1 genome.fasta(.gz) file found in ' + REFGENOMEFOLDER)
        else: # genome.fasta.gz
            refgenome = SeqIO.parse(gzip.open(fasta_file_gz[0], "rt"),'fasta')
    else: # genome.fasta
        refgenome = SeqIO.parse(fasta_file[0],'fasta')
    for c,record in enumerate(refgenome):
        if c == chrom:
            if strand == 1:
                return record.seq[int(start)-1:int(end)]
            elif strand == -1:
                return record.seq[int(start)-1:int(end)].reverse_complement()
            else:
                print('Unresolved strandedness.')


def save_cutoff_vals(timestamp, dict_of_filter_dicts_and_lists):
    with open(timestamp + '_cutoff_values.csv', 'w') as csv_file:  
        writer = csv.writer(csv_file)
        for value_name, list_dict_val in dict_of_filter_dicts_and_lists.items():
            if isinstance(list_dict_val, dict):
                for key, value in list_dict_val.items():
                    writer.writerow([value_name, key, value])
            elif isinstance(list_dict_val, list):
                for value in list_dict_val:
                    writer.writerow([value_name, '', value])
            else: 
                writer.writerow([value_name, '', list_dict_val])

def convert_str_to_value_type(value):
    """
    Converts strings to their value type (int, float or string)
    """
    if value != value:
        return np.nan
    try:
        return int(value)
    except (ValueError, TypeError):
        try:
            return float(value)
        except (ValueError, TypeError):
            return str(value)
        
def read_latest_cutoff_vals(working_directory):
    ## find lates cutoff val file saved and read in values

    ## identify saved latest file by timestamp in name
    latest_cutoff_file = glob.glob(working_directory + '/*/*_cutoff_values.csv')
    latest_cutoff_file = sorted(latest_cutoff_file)[-1]
    print(latest_cutoff_file)
    with open(latest_cutoff_file, 'r') as fo: 
        ## ensure that you read NEW values from file and not leaving old settings 
        saved_filters = {}
        ## loop over lines of file and extract 
        for line in fo:
            line = line.strip().split(',')
            if line[1] != '':
                if not line[0] in saved_filters.keys():
                    saved_filters[line[0]] = {}
                saved_filters[line[0]][line[1]] = convert_str_to_value_type(line[2])
            else:
                if not line[0] in saved_filters.keys(): ## single entry 
                    saved_filters[line[0]] = convert_str_to_value_type(line[2])
                elif isinstance(saved_filters[line[0]], list):
                    saved_filters[line[0]].append(convert_str_to_value_type(line[2]))
                else:
                    saved_filters[line[0]] = [saved_filters[line[0]], convert_str_to_value_type(line[2])] ## a second entry with the same var name identified --> list --> convert single string to list
                    saved_filters[line[0]].append(convert_str_to_value_type(line[2]))
    print(f'Cutoff filter variables read: {", ".join([filter for filter in saved_filters.keys()])}')
    return saved_filters.values() ## just return the values (--> filter_dicts, filter_lists, filter_values)


def div_major_allele_freq(cnts):
    # define matrices with major allele freq; major Nucl. (maNT); minor Nucl.; minor AF
    
    c=cnts[:,0:4,:]+cnts[:,4:8,:]; # flatten frw and rev ATCG counts    

    sorted_arr = np.sort(c,axis=1) # return sorted_arr matrix; axis 0 in 3d array == sort by col
    sortedpositions = np.argsort(c,axis=1) # return matrix indices of sort;axis 0 in 3d array == sort by col
    
    maxcount = sorted_arr[:,3:4:,:] # get allele counts for major allele (4th row); weird "3:4:" indexing required to maintain 3d structure
    minorcount = sorted_arr[:,2:3:,:] # get allele counts for first minor allele (3rd row); tri/quadro-allelic ignored; weird "2:3:" indexing required to maintain 3d structure and get 3rd row
    
    with np.errstate(divide='ignore', invalid='ignore'): # suppress warning for division by 0
        maf = maxcount / sorted_arr.sum(axis=1,keepdims=True)
        minorAF = minorcount / sorted_arr.sum(axis=1,keepdims=True)

    maf = np.squeeze(maf,axis=1) # turn 2d; axis=1 to keep 2d structure when only one position!
    maf[ np.isnan(maf) ] = 0 # set to 0 to indicate no data
    minorAF = np.squeeze(minorAF,axis=1) 
    minorAF[ np.isnan(minorAF) ] = 0 # set to 0 to indicate no data/no minor AF
    
    majorNT = np.squeeze(sortedpositions[:,3:4:,:],axis=1) # index position in sortedpositions represents allele position ATCG; axis=1 to keep 2d structure when only one position!
    minorNT = np.squeeze(sortedpositions[:,2:3:,:],axis=1) # -"-

    # Note: If counts for all bases are zero, then sort won't change the order
    # (since there is nothing to sort), thus majorNT.minorNT put to 4 (ie NA) using maf (REMEMBER: minorAF==0 is a value!)

    majorNT[maf==0] = 4
    minorNT[maf==0] = 4    
    
    return [maf.transpose(), majorNT.transpose(), minorNT.transpose(), minorAF.transpose()]


def major_allele(arr):
    ''' returns 1-dimensional array of size arr.shape[0] with the major allele index [0:3] or NA [4] (if ambigous)'''
    # NOTE: if major NT ambigous (multiple alleles with same number of occurence I report allele with lower numeric value. could be improved as it might lead to negligible (imo) bias). Also loop could prob be removed.
    # input 2d arr (2nd dim can be 1!) with nucleotide indices (0:4)
    nonNA_out = (arr != 4)
    out_idx = []
    for i,row in enumerate(arr):
        if np.any(nonNA_out[i,:]): # any majorNT found
            row = row[nonNA_out[i,:]]
            row_ct = np.unique(row,return_counts=True)
            idx = np.where(row_ct[1] == np.max(row_ct[1]) )[0]
            out_idx.append(row_ct[0][idx][0]) # if multiple outgroup alleles same count I take the first. Could be changed to NA or add more outgroup samples for refined inference.
        else: # no majorNT
            out_idx.append(4)
    out_idx = np.array(out_idx)
    print(np.sum(out_idx == 4),'/', out_idx.size ,'elements of p have no major allele (ie. 4)!'  )
    return out_idx


def ana_mutation_quality(calls,quals,ncpu=1):
    # This functions aims at providing the FQ value for every SNP position 
    # Across all pairwise different allele calls, it reports the best FQ value among the minimum FQ values per pair
    # NOTE: This function requires some more efficient coding!
    [Nmuts, NStrain] = calls.shape ;
    
    # generate template index array to sort out strains gave rise to reported FQ values
    s_template=np.zeros( (NStrain,NStrain) ,dtype=object)
    for i in range(s_template.shape[0]):
        for j in range(s_template.shape[1]):
            s_template[i,j] = str(i)+"_"+str(j)

    if ncpu == 1:
        MutQual = np.zeros((Nmuts,1)) ; 
        MutQualIsolates = np.zeros((Nmuts,2)); 
        for k in range(Nmuts):
            if k%2000 == 0: ## print pattern to stdout to show where in evaluation process runs 
                if k%10000 == 0:
                    print('|', end = '')
                else:
                    print('.', end = '')
            if len(np.unique(np.append(calls[k,:], 4))) <= 2: # if there is only one type of non-N (4) call, skip this location
                MutQual[k] = np.nan ;
                MutQualIsolates[k,:] = 0; 
            else:
                # c = calls[k,:] ; c1 = np.tile(c,(c.shape[0],1)); c2 = c1.transpose() # extract all alleles for pos k and build 2d matrix and a transposed version to make pairwise comparison
                # q = quals[k,:] ; q1 = np.tile(q,(q.shape[0],1)); q2 = q1.transpose() # -"-
                # g = np.all((c1 != c2 , c1 != 4 , c2 != 4) ,axis=0 )  # no data ==4; boolean matrix identifying find pairs of samples where calls disagree (and are not N) at this position
                c = calls[k,:] ; # extract all alleles for pos k 
                q = quals[k,:] ; # extract all quals for pos k 
                c1, c2 = np.meshgrid(c, c, indexing='xy')  # Create pairwise combinations of c
                q1, q2 = np.meshgrid(q, q, indexing='xy')  # Create pairwise combinations of q
                g = (c1 != c2) & (c1 != 4) & (c2 != 4) # no data ==4; boolean matrix identifying find pairs of samples where calls disagree (and are not N) at this position
                # Compute the maximum of the minimum values satisfying the condition
                #positive_pos = find(g); # numpy has no find; only numpy where, which does not flatten 2d array that way
                # get MutQual + logical index for where this occurred
                low_qual_disagree = np.minimum(q1[g],q2[g]) ## gives lower qual for each disagreeing pair of calls
                MutQual[k] = np.max(low_qual_disagree) # np.max(np.minimum(q1[g],q2[g])) find the best of disagreeing pairs; NOTE: np.max > max value in array; np.maximum max element when comparing two arryas
                MutQualIndex = np.argmax(low_qual_disagree) # return index of first encountered maximum!
                # get strain ID of reorted pair (sample number)
                s = s_template
                strainPairIdx = s[g][MutQualIndex]
                MutQualIsolates[k,:] = [strainPairIdx.split("_")[0], strainPairIdx.split("_")[1]]
    else:
        # Set up multiprocessing
        with Pool(processes=ncpu) as pool:
            results = pool.starmap(process_mutqual_row, [(k, calls[k,:], quals[k,:], s_template) for k in range(Nmuts)])
            # results = pool.map(process_mutqual_row, range(Nmuts))
        ## sort results by k
        results_sorted = sorted(results, key=lambda res: res[0])
        MutQual = np.array([res[1] for res in results_sorted]).reshape(-1, 1)
        MutQualIsolates = np.array([res[2] for res in results_sorted]).reshape(-1, 2)
    print('')
    return [MutQual,MutQualIsolates]

def process_mutqual_row(k, c, q, s_template):
    if k%2000 == 0: ## print pattern to stdout to show where in evaluation process runs 
        if k%10000 == 0:
            print('|', end = '')
        else:
            print('.', end = '')
    if len(np.unique(np.append(c, 4))) <= 2: # if there is only one type of non-N (4) call, skip this location
        mutqual_k = np.nan ;
        mutqualisolate_k = [0, 0]; 
    else:
        # c = calls[k,:] ; c1 = np.tile(c,(c.shape[0],1)); c2 = c1.transpose() # extract all alleles for pos k and build 2d matrix and a transposed version to make pairwise comparison
        # q = quals[k,:] ; q1 = np.tile(q,(q.shape[0],1)); q2 = q1.transpose() # -"-
        # g = np.all((c1 != c2 , c1 != 4 , c2 != 4) ,axis=0 )  # no data ==4; boolean matrix identifying find pairs of samples where calls disagree (and are not N) at this position
        c1, c2 = np.meshgrid(c, c, indexing='xy')  # Create pairwise combinations of c
        q1, q2 = np.meshgrid(q, q, indexing='xy')  # Create pairwise combinations of q
        g = (c1 != c2) & (c1 != 4) & (c2 != 4) # no data ==4; boolean matrix identifying find pairs of samples where calls disagree (and are not N) at this position
        # Compute the maximum of the minimum values satisfying the condition
        #positive_pos = find(g); # numpy has no find; only numpy where, which does not flatten 2d array that way
        # get MutQual + logical index for where this occurred
        low_qual_disagree = np.minimum(q1[g],q2[g]) ## gives lower qual for each disagreeing pair of calls
        mutqual_k = np.max(low_qual_disagree) # np.max(np.minimum(q1[g],q2[g])) find the best of disagreeing pairs; NOTE: np.max > max value in array; np.maximum max element when comparing two arryas
        MutQualIndex = np.argmax(low_qual_disagree) # return index of first encountered maximum!
        # get strain ID of reorted pair (sample number)
        s = s_template
        strainPairIdx = s[g][MutQualIndex]
        mutqualisolate_k = [int(strainPairIdx.split("_")[0]), int(strainPairIdx.split("_")[1])]
    return (k, mutqual_k, mutqualisolate_k)  # Include the index `k` in the result for sorting

"""
def ana_mutation_quality(Calls,Quals):
    # This functions aims at providing the FQ value for every SNP position 
    # Across all pairwise different allele calls, it reports the best FQ value among the minimum FQ values per pair
    # NOTE: This function requires some more efficient coding!
    [Nmuts, NStrain] = Calls.shape ;
    MutQual = np.zeros((Nmuts,1)) ; 
    MutQualIsolates = np.zeros((Nmuts,2)); 
    
    # generate template index array to sort out strains gave rise to reported FQ values
    s_template=np.zeros( (len(Calls[0,:]),len(Calls[0,:])) ,dtype=object)
    for i in range(s_template.shape[0]):
        for j in range(s_template.shape[1]):
            s_template[i,j] = str(i)+"_"+str(j)

    for k in range(Nmuts):
        if len(np.unique(np.append(Calls[k,:], 4))) <= 2: # if there is only one type of non-N (4) call, skip this location
            MutQual[k] = np.nan ;
            MutQualIsolates[k,:] = 0; 
        else:
            c = Calls[k,:] ; c1 = np.tile(c,(c.shape[0],1)); c2 = c1.transpose() # extract all alleles for pos k and build 2d matrix and a transposed version to make pairwise comparison
            q = Quals[k,:] ; q1 = np.tile(q,(q.shape[0],1)); q2 = q1.transpose() # -"-
            g = np.all((c1 != c2 , c1 != 4 , c2 != 4) ,axis=0 )  # no data ==4; boolean matrix identifying find pairs of samples where calls disagree (and are not N) at this position
            #positive_pos = find(g); # numpy has no find; only numpy where, which does not flatten 2d array that way
            # get MutQual + logical index for where this occurred
            MutQual[k] = np.max(np.minimum(q1[g],q2[g])) # np.max(np.minimum(q1[g],q2[g])) gives lower qual for each disagreeing pair of calls, we then find the best of these; NOTE: np.max > max value in array; np.maximum max element when comparing two arryas
            MutQualIndex = np.argmax(np.minimum(q1[g],q2[g])) # return index of first encountered maximum!
            # get strain ID of reorted pair (sample number)
            s = s_template
            strainPairIdx = s[g][MutQualIndex]
            MutQualIsolates[k,:] = [strainPairIdx.split("_")[0], strainPairIdx.split("_")[1]]
            
    return [MutQual,MutQualIsolates]
"""

def filter_bed_cov_hist(bed_path,p,scafNames,chrStarts,sampleNames,coverage,cutoff,two_tailed=False,upper=True):
    """
    NOTE: For ancient DNA quality control!!!
    outputs a boolean matrix of index of [p]x[samplenames] to mask from basecall, due to being in the top coverage percentile when using bedtools coverage histogram (or bottom if upper=False, or if two_tailed=True)"""
    # get paths
    bed_histogram_files = glob.glob(bed_path)
    # index accounting
    p_chr_indices=[0]+[np.max(np.where(p < x + 1))+1 for x in chrStarts[1:]]
    # generate output matrix
    to_be_masked_array_covg_percentile=np.full(coverage.shape,False)
    for bed_index,bed_hist in enumerate(bed_histogram_files):
        this_sample_index=-1
        for index,samplename in enumerate(sampleNames):
            if samplename in bed_hist:
                this_sample_index=index
        if this_sample_index==-1:
            print(f'Warning: Bedfile {bed_hist} does not have a corresponding samplename in sampleNames input array, skipping.')
        else:
            this_sample_name_cutoffs=cutoff_bed_covg_histograms(bed_hist,cutoff,two_tailed,upper)
            ## now, mask positions above the percentile, by individual chrom cutoffs
            for index,scaf in enumerate(scafNames):
                this_sample_name_cutoffs_this_scaf=this_sample_name_cutoffs[scaf][1]
                if index<len(scafNames)-1:
                    start=p_chr_indices[index]
                    end=p_chr_indices[index+1]
                    p_to_include_this_chrom=np.array(range(start,end)) ## true/false 
                else:
                    p_to_include_this_chrom=np.array(range(start,len(p)))
                to_mask=p_to_include_this_chrom[np.isin(p_to_include_this_chrom,np.where(coverage[:,this_sample_index] > this_sample_name_cutoffs_this_scaf)[0])]
                to_be_masked_array_covg_percentile[to_mask,this_sample_index]=True
    return to_be_masked_array_covg_percentile

def filter_bed_0_cov_regions(bed_path,p,scafNames,chrStarts,sampleNames,cutoff):
    """outputs a boolean matrix of index of [p]x[samplenames] to mask from basecall, due to being in regions of 0 coverage (eg are coverage islands)"""
    # get paths
    bed_zero_covg_files=glob.glob(bed_path)
    # generate output matrix
    to_be_masked_array_0_covg_regions=np.full((len(p),len(sampleNames)),False)
    for bed_index,bed_zero in enumerate(bed_zero_covg_files):
        this_sample_index=-1
        for index,samplename in enumerate(sampleNames):
            if samplename in bed_zero:
                this_sample_index=index
        if this_sample_index==-1:
            print(f'Warning: Bedfile {bed_zero} does not have a corresponding samplename in sampleNames input array, skipping.')
        else:
            to_mask=parse_bed_zero_covg_regions(bed_zero,p,scafNames,chrStarts,cutoff)
            ## now, mask positions above the percentile 
            to_be_masked_array_0_covg_regions[to_mask,this_sample_index]=True
    return to_be_masked_array_0_covg_regions


def cutoff_bed_covg_histograms(path_to_bed_covg_hist,cutoff,two_tailed=True,upper=True):
    """generates dictionary of scafnames --> cutoff values for coverage, based on bedtools covg histogram"""
    # get cutoffs for coverage
    if cutoff >= 0.5: cutoff=1-cutoff
    with open(path_to_bed_covg_hist) as file:
        lines = csv.reader(file, delimiter="\t")
        output_hist={}
        for line in lines:
            if line[0] not in output_hist:
                output_hist[line[0]] = {line[1]:line[4]}
            else: output_hist[line[0]][line[1]] = line[4]
    if two_tailed:
        cutoffs=(cutoff,1-cutoff)
    else: 
        if upper: cutoffs=(0,1-cutoff)
        else: cutoffs=(cutoff,1)
    cutoffs_scafs={scaf:[0,0] for scaf in output_hist.keys()}
    for scaf in output_hist:
        cumulative_perctile=0
        for cov,value in output_hist[scaf].items():
            cov,value=float(cov),float(value)
            cumulative_perctile+=value
            if cumulative_perctile < cutoffs[0]:
                cutoffs_scafs[scaf][0]=cov
            if cumulative_perctile > cutoffs[1]:
                cutoffs_scafs[scaf][1]=cov
                break
    return cutoffs_scafs

def parse_bed_zero_covg_regions(path_to_bed_zero_covg_covg,p,scafNames,chrStarts,cutoff):
    """generates array of positions that should be masked due to falling within coverage islands. Reach to get to next coverage island defined elsewhere (snakemake, other script) but cutoff defines % of uncovered region which must be uncovered to mask any covered positions within
    
    eg: covg = =======
        noncoverd = _
        genome: ====___________________=_____________
        
        remove the islated = covered position as likely contamination
    """
    with open(path_to_bed_zero_covg_covg) as file:
        lines = csv.reader(file, delimiter="\t")
        output_set=set()
        for line in lines:
            if float(line[6]) > cutoff:
                to_add_to_range=chrStarts[np.where(scafNames==line[0])[0]][0]
                output_set.update([x for x in range(int(line[1])+to_add_to_range,int(line[2])+1+to_add_to_range)])
    should_be_masked=[]
    for index,pos in enumerate(p):
        if pos in output_set:
            should_be_masked.append(index)
    should_be_masked=np.array(should_be_masked)
    return should_be_masked

def combinatorial_filtering_to_goodpos(hasmutation_unfiltered,list_of_filters,list_of_filternames,no_filter_value=0):
    """Generates two lists, data and combination_of_filternames, which store the number of sites that get removed by applying a given combination of filters to candidate mutation table. Generates input for upsetplot.from_memberships(filtercomb_out,data=data_out)"""
    combination_of_filters=[]
    data=[no_filter_value]
    combination_of_filternames=[[]]
    for i in range(1,len(list_of_filters)+1):
        combination_of_filters+=list(itertools.combinations(list(range(len(list_of_filters))),i))
    for combination in combination_of_filters:
        filters_to_apply=[hasmutation_unfiltered]+[list_of_filters[x] for x in combination]
        combination_of_filternames.append([list_of_filternames[x] for x in combination])
        data.append(len(np.setxor1d(np.where(np.sum(hasmutation_unfiltered,axis=1)-np.sum(bool_numpy_combination(filters_to_apply),axis=1)==0)[0],np.where(np.sum(hasmutation_unfiltered,axis=1)==0)[0])))
    return data,combination_of_filternames

def bool_numpy_combination(list_of_bool_arrays):
    """Helper function for combinatorial_filtering_to_goodpos, creates boolean array to sum over to find goodpos (positions with at least 1 believed mutation), after applying a set of filter results (list_of_bool_arrays) """
    bool_array_out=list_of_bool_arrays[0]
    filters=np.full(list_of_bool_arrays[0].shape,False)
    for i in range(1,len(list_of_bool_arrays)):
        filters=(filters | list_of_bool_arrays[i])
    bool_array_out=(bool_array_out & filters)
    return bool_array_out

def combinatorial_filtering_to_goodpos_within_samples(hasnomutation_unfiltered, dict_of_filters, ingroup_bool):
    # Calculate counts for individual filters and combinations of filters
    # Initialize a dictionary to hold the count of each combination
    intersection_counts = {}
    mask_names = list(dict_of_filters.keys())
    for i in range(1, len(dict_of_filters) + 1):
        for comb in itertools.combinations(mask_names, i):
            # Start with a mask of all True values (all rows included)
            combined_mask = hasnomutation_unfiltered[:, ingroup_bool]
            for mask_name in comb:
                combined_mask |= dict_of_filters[mask_name][:, ingroup_bool]
            intersection_counts[comb] = np.sum((np.sum(combined_mask, axis = 1)==np.shape(hasnomutation_unfiltered[:, ingroup_bool])[1]))
    return intersection_counts

def combinatorial_filtering_to_goodpos_across_sites(dict_of_filters, num_mut):
    # Calculate counts for individual filters and combinations of filters
    # Initialize a dictionary to hold the count of each combination
    intersection_counts = {}
    mask_names = list(dict_of_filters.keys())
    for i in range(1, len(dict_of_filters) + 1):
        for comb in itertools.combinations(mask_names, i):
            # Start with a mask of all True values (all rows included)
            combined_mask = np.zeros(num_mut, dtype = bool)
            maskname_tuple = ()
            for mask_name in comb:
                combined_mask |= dict_of_filters[mask_name]
                if isinstance(mask_name, str):
                    maskname_tuple += (mask_name,)
                else:
                    maskname_tuple += tuple([submask for submask in mask_name])
            intersection_counts[maskname_tuple] = np.sum(combined_mask)
    return intersection_counts

def upsetplot_filtering(intersection_counts_dict, timestamp, subject, refgenome, title_suffix=''):
    with warnings.catch_warnings():
        # Cathc warning:  FutureWarning: Downcasting object dtype arrays on .fillna, .ffill, .bfill is deprecated and will change in a future version. Call result.infer_objects(copy=False) instead. To opt-in to the future behavior, set `pd.set_option('future.no_silent_downcasting', True)` df.fillna(False, inplace=True)
        warnings.filterwarnings("ignore", category=FutureWarning, message=".*fillna.*downcasting.*")
        upset_data = upsetplot.from_memberships(intersection_counts_dict.keys(), data = list(intersection_counts_dict.values()))
    # Create UpSet plot
    with warnings.catch_warnings():
        # Cathc warning:  FutureWarning: A value is trying to be set on a copy of a DataFrame or Series through chained assignment using an inplace method. The behavior will change in pandas 3.0. This inplace method will never work because the intermediate object on which we are setting values always behaves as a copy.
        warnings.filterwarnings("ignore", category=FutureWarning, message=".*A value is trying to be set on a copy of a DataFrame or Series through chained assignment using.*")
        fig = plt.figure(figsize=(10, 4))
        upset = upsetplot.UpSet(upset_data, show_counts=True)
        upsetplt = upset.plot(fig = fig)
        upsetplt["intersections"].set_ylabel("SNVs removed")
        upsetplt["totals"].remove()
        plt.ylim(ymin = 0)
    plt.savefig(f'pdf/{timestamp}_{subject}_{refgenome}_removed_sites_by_filter{title_suffix}.pdf')

def parse_gff(REFGENOMEFOLDER,ScafNames,ortholog_info_series=pd.Series(dtype = 'int64'),forceReDo=False):
    # parse gff file (tested with version3) with genome annotation (gff or gff.gz)
    # requires one gff file in REFGENOMEFOLDER, which is detected automatically
    # provide ScafNames to maintain same order of contigs/chr/scaffolds in output list as in previously generated ref-genome-based variables, eg. contig_positions,ChrStarts, GenomeLength
    # some columns potentially contain multiple entries. Those are separated by ";"
    # no data is always reported as "."
    # https://biopython.org/wiki/GFF_Parsing 
    # NOTE: annotation is read from gff file!
    # only execute if dataframes not yet made
    short_aa_sequence_counter = 0
    if os.path.isfile(REFGENOMEFOLDER+"/annotation_genes.pandas.py.pk1") and (forceReDo == False):
        afile = open(REFGENOMEFOLDER+"/annotation_genes.pandas.py.pk1", 'rb')
        list_of_dataframes = pickle.load(afile)
        afile.close()
        return list_of_dataframes
    else: # search for gff or gff.gz
        examiner = GFF.GFFExaminer()        
        gff_file = glob.glob(REFGENOMEFOLDER + '/*.gff*') # search for gff or gff.gz or gff3
        if len(gff_file) != 1:
            raise ValueError('Either no gff(.gz) file or more than 1 *gff(.gz) file found in ' + REFGENOMEFOLDER)
        if gff_file[0][-2:] == 'gz':
            encoding = 'gzip'
        else:
            encoding = 'unzip'
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open # define opening procedure gzip.open or open
        with _open(gff_file[0]) as gff_handle:
            possible_limits = examiner.available_limits(gff_handle)
        # start processing
        chromosomes = possible_limits["gff_id"].keys() # list of tuples containing contig ID; do not loop over those, as I only loop over contigs observed (ScafNames)...might need a fix if !=
        tagnumber_counter = 0 # unique, continous numerical id for all features across all chromosomes
        list_of_dataframes = [] # output: each element contains a dtaframe w/ all annotation info for chr. ordered as in ScafNames.
        # define gff_type to extracted. Use all but gene, gap and region. gene annotated in NCBI but lacks most info that extra CDS has. Region just describes chromosome extensions, gaps describes gaps in scaffolds annotated by bakta
        annotation_types_to_exclude = ['gene', 'region', 'gap']
        limits = dict(gff_type = [i[0] for i in possible_limits['gff_type'].keys() if not i[0] in annotation_types_to_exclude] )
        # loop over chr with established order
        for chrom in ScafNames:
            limits["gff_id"] = [chrom]
            # limits = dict(gff_id=[chrom])
            with _open(gff_file[0]) as gff_handle:
                for rec in GFF.parse(gff_handle, limit_info=limits):
                    # for loop over every chr but only defined [limits] has data. Memory efficient!
                    if rec.id == chrom:
                        # if chr has any feature build list of dicts and append to list_of_dataframes, else append empty dataframe
                        if len(rec.features) > 0:
                            # test if seq object part of gff (prokka-based yes, but NCBI-based no >> then load ref genome.fasta)
                            # if len(rec.seq) == rec.seq.count('?'):
                            for seq_record in SeqIO.parse(REFGENOMEFOLDER+"/genome.fasta", "fasta"):
                                if seq_record.id == rec.id:
                                    rec.seq = seq_record.seq
                            if len(rec.seq) == rec.seq.count('?'): # test if succesful
                                print('Warning: No reference genome found that matches chromosome:' + rec.id)
                            lod_genes = [] # list-of-dictionary. Easy to turn to pd.dataframe
                            for gene_feature in rec.features:
                                gene_dict = {}
                                tagnumber_counter += 1
                                gene_dict['type'] = gene_feature.type
                                gene_dict['locustag'] = gene_feature.id
                                # add ortholog info if locustag (eg. repeat region has none)
                                if gene_feature.id != "" and gene_feature.type == 'CDS' and not ortholog_info_series.empty:
                                    try: 
                                        gene_dict['orthologtag'] = ortholog_info_series[ortholog_info_series.str.findall(gene_feature.id).str.len() == 1].index[0]
                                    except:
                                        short_aa_sequence_counter +=1
                                        print('Gene ID:')
                                        print(gene_feature.id)
                                        print('Length of sequence:')
                                        print(len(rec.seq[gene_feature.location.start:gene_feature.location.end].translate(table = "Bacterial")))
                                        gene_dict['orthologtag'] = np.nan ## NOTE: Bakta annotates some sequences < 10 aa which cd-hit excludes by default!
                                #print(rec.id+"   "+gene_dict['locustag'])
                                if 'gene' in gene_feature.qualifiers.keys():
                                    gene_dict['gene'] = ";".join(gene_feature.qualifiers['gene'])
                                else:
                                    gene_dict['gene'] = "." # add "." instead of []
                                gene_dict['type'] = gene_feature.type
                                gene_dict['locustag'] = gene_feature.id
                                if gene_dict['type'] == "CDS" or gene_dict['type'] == "gene":
                                    gene_dict['tagnumber'] = tagnumber_counter
                                else:
                                    gene_dict['tagnumber'] = 0
                                if 'product' in gene_feature.qualifiers.keys():
                                    gene_dict['product'] = ";".join(gene_feature.qualifiers['product']) 
                                else:
                                    gene_dict['product'] = "."
                                if 'protein_id' in gene_feature.qualifiers.keys():
                                    gene_dict['protein_id'] = gene_feature.qualifiers['protein_id']
                                else:
                                    gene_dict['protein_id'] = "."
                                if "Dbxref" in gene_feature.qualifiers.keys():
                                    gene_dict['db_xref'] = ";".join(gene_feature.qualifiers['Dbxref'])
                                else:
                                    gene_dict['db_xref'] = "."
                                if "note" in gene_feature.qualifiers.keys():
                                    gene_dict['note'] = ";".join(gene_feature.qualifiers['note'])
                                elif "Note" in gene_feature.qualifiers.keys():
                                    gene_dict['note'] = ";".join(gene_feature.qualifiers['Note'])
                                else:
                                    gene_dict['note'] = "."
                                if 'phase' in gene_feature.qualifiers.keys():
                                    gene_dict['codon_start'] = int(gene_feature.qualifiers['phase'][0])
                                else:
                                    gene_dict['codon_start'] = "."
                                gene_dict['indices'] = [int(gene_feature.location.start),int(gene_feature.location.end)]
                                gene_dict['loc1'] = int(gene_feature.location.start) # automatically 0-based
                                gene_dict['loc2'] = int(gene_feature.location.end)
                                gene_dict['location'] = gene_feature.location
                                gene_dict['strand'] = gene_feature.location.strand
                                dna_seq = rec.seq[gene_feature.location.start:gene_feature.location.end]
                                if gene_dict['strand'] == 1:
                                    gene_dict['sequence'] = dna_seq
                                elif gene_dict['strand'] == -1:
                                    gene_dict['sequence'] = dna_seq.reverse_complement()
                                else:
                                    gene_dict['sequence'] = dna_seq # eg. repeat region
                                # translation, add info where codon starts if info was available. Usually it starts at 0
                                if isinstance( gene_dict['codon_start'] , int):
                                    sequence2translate = gene_dict['sequence'][gene_dict['codon_start']:]
                                    gene_dict['translation'] = sequence2translate.translate(table="Bacterial") # bacterial genetic code GTG is a valid start codon, and while it does normally encode Valine, if used as a start codon it should be translated as methionine. http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:translation
                                elif gene_dict['type'] == "CDS":
                                    sequence2translate = gene_dict['sequence']
                                    gene_dict['translation'] = sequence2translate.translate(table="Bacterial")
                                else:
                                    gene_dict['translation'] = "." # all non-CDS (RNA's or repeat regions) not translated (as those are sometimes also off-frame)
                                lod_genes.append(gene_dict)
                            # after all features parsed > add as dataframe to output list
                            # ATTENTION: SORT pandas DF (annotation not necessarily sorted)
                            df_sort = pd.DataFrame(lod_genes)
                            df_sort = df_sort.sort_values(by=['loc1'])
                            list_of_dataframes.append(df_sort)
                            
                        else:
                            list_of_dataframes.append( pd.DataFrame([{}]) )
        print(f'Bakta annotates some sequences < 10 aa which cd-hit excludes by default! Identified {short_aa_sequence_counter} such instances in gff file.')
        afile = open(REFGENOMEFOLDER+"/annotation_genes.pandas.py.pk1", 'wb')
        pickle.dump(list_of_dataframes, afile)
        afile.close()
        return list_of_dataframes

def idx2nts(calls,missingdata="?"):
    # translate index array to array containing nucleotides
    nucl = np.array(['A','T','C','G',missingdata],dtype=object) # add 5th element --> no data! == index 4
    palette = [0,1,2,3,4] # values present in index-array
    index = np.digitize(calls.ravel(), palette, right=True)
    return nucl[index].reshape(calls.shape)

def nts2idx(nts_array,missingdata="?"):
    # translate nucl array to array containing nucleotide indices
    nucl = np.array(['A','T','C','G',missingdata],dtype=object) # add 5th element --> no data! == index 4
    indices = [0,1,2,3,4] # values present in index-array
    for i in indices:
        nts_array[ nts_array == nucl[i] ] = i
    return nts_array.astype(int)


def genomic_position_all(anno_genes_ls,genomelength,chrstarts):
    # get initial placeholder for all genes (differentiate by chr and genic/nongenic)  
    # UPDATE: gene_num > 0 means CDS. changed 'tagnumber' assignment in parse_gff()
    # 2024/02 MF: generate dict to assign for each pos all possible genetic regions, if genetic regions overlap
    # 2024/02 MF: select assignment to longest gene, if regulatory region is longer than CDS, select CDS!

    locus_tag_nums = np.ones(genomelength,dtype=float)*0.5 # genome-wide unique CDS tag ('tagnumber'), and intergenic 0.5
    cds_indices = np.ones(genomelength,dtype=float)*0.5 # per chromosome unique; intergenic #chr.5-1; intragenic #gene
    pos_to_genenum_dict = {pos: [] for pos in range(genomelength)}
    
    for i,this_chr_df in enumerate(anno_genes_ls):
        if this_chr_df.empty: #Skip contigs w/o CDS annotation.
            continue
        
        gene_num = this_chr_df['tagnumber'].values # tdl called cds_num but actually it was not always cds. e.g tRNA; .values  turn pd.DF to 1D numpy array
        genestarts = chrstarts[i] + this_chr_df['loc1'].values ; # .values  turn pd.DF to 1D numpy array
        geneends = chrstarts[i] + this_chr_df['loc2'].values ; # .values  turn pd.DF to 1D numpy array
        genetype = this_chr_df['type'].values

        for j in range(len(genestarts)-1): # loop over all but last element
            for genomepos in range(genestarts[j], geneends[j]):
                if len(pos_to_genenum_dict[genomepos]) > 0:
                    if len(pos_to_genenum_dict[genomepos][0]) == 1:
                        pos_to_genenum_dict[genomepos] = pos_to_genenum_dict[genomepos][1:] ## remove intergenic if there is one stored
                pos_to_genenum_dict[genomepos].append([gene_num[j], j+1, abs((genestarts[j]) - geneends[j]), genetype[j]]) # set gene number, chromosome and gene length

            for genomepos in range(geneends[j],(genestarts[j+1])):
                if pos_to_genenum_dict[genomepos] == []: ## check if intergenic is still possible 
                    pos_to_genenum_dict[genomepos].append([j+1+0.5])# intergenic are #chr+0.5; starts with 0.5 prior chr1 which is set when array created
                    
        # mark positions for last gene of chr
        for genomepos in range(genestarts[-1], geneends[-1]):
            if genomepos > genomelength:
                print('Annotation exceeding the genome length boundaries!') ## NOTE: Can happen e.g.: in genome ASM1413175v1
                #break
            pos_to_genenum_dict[genomepos].append([gene_num[-1], len(gene_num), abs((genestarts[-1]) - geneends[-1]), genetype[-1]]) # set gene number, chromosome and gene length
        # mark trailing cds_intergenic from last gene till end of chromosome (if last chr >> till end of genome)
        if ((i+1) < len(chrstarts)):
            for genomepos in range(geneends[-1],chrstarts[i+1]):
                if pos_to_genenum_dict[genomepos] == []: ## check if intergenic is still possible 
                    pos_to_genenum_dict[genomepos].append([len(gene_num)+0.5])
        else:
            for genomepos in range(geneends[-1],genomelength):
                if pos_to_genenum_dict[genomepos] == []: ## check if intergenic is still possible 
                    pos_to_genenum_dict[genomepos].append([len(gene_num)+0.5])

    ## select associations based on longest gene 
    for key, entry in pos_to_genenum_dict.items():
        if len(entry) == 1: 
            if (len(entry[0]) == 1) & (abs(cds_indices[key-1]-entry[0][0]) < 1): ## check if entry is length of 1 --> intergenic, and the difference genenumber between this and the pos before is less then 1 (intergenics between genes have gene number!)
                cds_indices[key] = entry[0][0]
            elif (len(entry[0]) == 1): ## check if intergenic but gene number is further away than that from previos position --> there is a gene within another gene --> update intergenic position to correct association!
                if cds_indices[key-1] == int(cds_indices[key-1]):
                    cds_indices[key] = cds_indices[key-1] + 0.5
                else:
                    cds_indices[key] = cds_indices[key-1]
            else:
                [locus_tag_nums[key], cds_indices[key]] = entry[0][:2]
        elif len(entry) > 1:
            where_cds = np.where(np.asanyarray(entry) == 'CDS')[0] ## get if there are coding sequences --> favour them! 
            if (len(where_cds) > 0) & (len(where_cds) < len(entry)): ## select longest coding sequence at position (if not all are CDS but at least one --> else directly subset!)
                subset_of_entries = [entry[idx] for idx in where_cds]
                max_entry = max(subset_of_entries, key=lambda x: x[2])    
            else: ## select longest genetic region 
                max_entry = max(entry, key=lambda x: x[2])
            [locus_tag_nums[key], cds_indices[key]] = max_entry[:2]
    return [locus_tag_nums,cds_indices]

def annotate_mutations(annotation_genes , p_gp , refnti_gp , ancnti_gp , calls_gp , counts_gp , hasmutation_gp , mutQual_gp, promotersize, ref_genome_folder, goodpos2use=[]):
    ''' produce dataframe with annotation info for each goodposition used for tree
    all positions are 1-based! (pos, nt_pos, loc1/2, distance1/2)
    # function combines TDL: annotate_mutations_gb.m and append_annotations.m 
    CHANGELOG:
    17.10.2022 IL: Added parsing of ancestral AA state by getting annotation info from ancnti_gp
    02/2024 MF: changed AA gt return to AA ref, anc and alt and evaluation of nonsyn mutations'''
    # p_gp: genomic position of goodpos
    # p_gp=p[goodpos2use]
    # refnti_gp=refnti_m[goodpos2use,:]
    # calls_gp=calls[goodpos2use,:]
    # counts_gp=counts[:,:,goodpos2use]
    # hasmutation_gp=hasmutation[goodpos2use,:]
    # mutQual_gp = mutQual[goodpos2use,].flatten() 
    # ref_genome_folder=REFGENOMEFOLDER
    if len(goodpos2use) > 0:
        p_gp = p_gp[goodpos2use]
        refnti_gp = refnti_gp[goodpos2use,:]
        ancnti_gp = ancnti_gp[goodpos2use,:]
        calls_gp = calls_gp[goodpos2use,:]
        counts_gp = counts_gp[:,:,goodpos2use]
        hasmutation_gp = hasmutation_gp[goodpos2use,:]
        mutQual_gp = mutQual_gp[goodpos2use,].flatten()

    [maf_gp, maNT_gp, minorNT_gp, minorAF_gp] = div_major_allele_freq(counts_gp)
    nts = ['A','T','C','G'] #'atcg';
    rc = ['T','A','G','C']  # 'tagc';
    [chrstarts, genomelength,scafnames]= genomestats(ref_genome_folder);
    [locus_tag,chr_tag] = genomic_position_all(annotation_genes, genomelength, chrstarts); #loc_tag numbers all genes across genome (intergenic:0.5); chr_tag: numbers genes per chr (intergenic: #chr+0.5)
    mut_genenum = chr_tag[p_gp]
    loc_tag_p_gp = locus_tag[p_gp]
    pos_chr = p2chrpos(p_gp,chrstarts) # get two-col table chr,pos   
    lod_mutAnno = []
    for i,pos in enumerate(p_gp):

        mut_annotations = {}
        mut_annotations['gene_num'] = mut_genenum[i]
        mut_annotations['gene_num_global'] = loc_tag_p_gp[i]
        mut_annotations['chr'] = pos_chr[i][0]
        mut_annotations['pos'] = pos_chr[i][1] + 1 # turn 0-based pos into 1-based pos.
        mut_annotations['p'] = p_gp[i]
        mut_annotations['quals'] = mutQual_gp[i]
        mut_annotations['nts'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_ref'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_anc'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_alt'] = "." # default. overwritten below if NT defined.
        p_chr = pos_chr[i][0]
        if mut_genenum[i] == int(mut_genenum[i]): # intragenic
            p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])-1] # both -1 necessary bcs list of df and df 0-based
            mut_annotations['product'] = p_anno.loc['product']
            mut_annotations['gene'] = p_anno.loc['gene']
            mut_annotations['protein_id'] = p_anno.loc['protein_id']
            mut_annotations['strand'] = p_anno.loc['strand']
            mut_annotations['loc1'] = p_anno.loc['loc1'] + 1 # 1-based first position
            mut_annotations['loc2'] = p_anno.loc['loc2'] # last position of gene (inclusive)
            mut_annotations['sequence'] = p_anno.loc['sequence']
            mut_annotations['protein_id'] = p_anno.loc['protein_id']
            mut_annotations['note'] = p_anno.loc['note']
            mut_annotations['locustag'] = p_anno.loc['locustag']
            if 'orthologtag' in p_anno:
                mut_annotations['orthologtag'] = p_anno.loc['orthologtag']
            if p_anno.loc['strand'] == 1: # get position within gene, consider strandedness
                mut_annotations['nt_pos'] = mut_annotations['pos'] - mut_annotations['loc1']+1 # pos/loc1 1-based. nt_pos 1-based. +1 to get 1-based nt_pos in gene (checked)!
            elif p_anno.loc['strand'] == -1:
                mut_annotations['nt_pos'] = p_anno.loc['loc2'] - mut_annotations['pos'] +1 # pos/loc2 1-based. +1 to get 1-based nt_pos in gene (checked)!
            else:
                mut_annotations['nt_pos'] = "." # can happen. eg. Crispr
            
            if p_anno.loc['type'] == 'CDS': ## check if CDS --> otherwise it is ncRNA or regulatory region (i.e. rRNA, tRNA, tmRNA, ncRNA, CRISPR, oriT, oriC, regulatory region,...)
                mut_annotations['translation'] = p_anno.loc['translation']
                if mut_annotations['nt_pos'] == ".": # Observed with special 'type's like crispr annotated. rare! leads to no AA inference.
                    aan = 9999999
                else:
                    aan = int( (mut_annotations['nt_pos']-1 )/3 ) + 1 # get #codon that harbours mutation. 1-based.
                    ncn = mut_annotations['nt_pos'] - ((aan-1)*3) # get #nucl within triplett. 1-based
                    mut_annotations['aa_pos'] = aan
                codons_ls = []
                aa_ls = [] #potential AA changes at given codon
                if len(mut_annotations['sequence']) >= (aan*3) and mut_annotations['translation'] != ".":
                    codon = mut_annotations['sequence'][aan*3-2-1:aan*3] # -1 bcs seq-object 0-based; but positional argument aan not 
                    codon = [n for n in codon  ] # turn seq.object to list, seq object does not allow reassignment
                    for n in range(4): # test all four nucleotide options for resulting AA change
                        if p_anno.loc['strand'] == 1:
                            codon[ncn-1] = nts[n]
                        else:
                            codon[ncn-1] = rc[n]
                        codonSeqO = Seq( "".join(codon)) 
                        codons_ls.append(codonSeqO)
                        aa_ls.append(codonSeqO.translate())
                mut_annotations['codons'] = codons_ls
                mut_annotations['AA'] = [aa[0] for aa in aa_ls]
                # append_annotations intragenic
                if len(mut_annotations['AA']) < 4:
                    mut_annotations['type'] = 'U'
                mut_annotations['AA_gt_ref'] = '.' # Fill in annotations with whether NS or Syn mutation
                mut_annotations['AA_gt_anc'] = '' 
                mut_annotations['nt_alt'] = ''
                mut_annotations['AA_gt_alt'] = ''
                if np.unique(refnti_gp[i,:])[0] != 4: # only if ref defined; [0] ok...see next line comment
                    mut_annotations['nt_ref'] = nts[np.unique(refnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                    mut_annotations['nts'] = mut_annotations['nt_ref']
                if np.unique(ancnti_gp[i,:])[0] != 4:
                    mut_annotations['nt_anc'] = nts[np.unique(ancnti_gp[i,:])[0]] # ancnti_gp should by definition never have > 1 unique element per row (one outgroup only!!! >> [0] apropriate
                if len(mut_annotations['AA']) == 4:
                    refAA = mut_annotations['AA'][ np.unique(refnti_gp[i,:])[0] ]
                    mut_annotations['AA_gt_ref'] = refAA # the AA list order corresponds to NTS list!
                # extract derived genotype(s) and according AA across all samples
                # call Nonsynonymous if >1 amino acid type exists across all called samples (eg ref AA and alt )
                for j,_ in enumerate(calls_gp[i,:]):
                    #fixed mutation
                    # i is index of current pos, j is index of samples 
                    if hasmutation_gp[i,j] == True :
                        if calls_gp[i,j] < 4:
                            # add each NT called, later only keep uniq, but kept order
                            mut_annotations['nts'] = mut_annotations['nts'] + nts[calls_gp[i,j]]
                            mut_annotations['nt_alt'] = mut_annotations['nt_alt'] + nts[calls_gp[i,j]]
                            if len(mut_annotations['AA']) == 4:
                                # mut_annotations['AA_gt'] = mut_annotations['AA_gt'] + mut_annotations['AA'][ maNT_gp[i,j] ]
                                mut_annotations['AA_gt_alt'] = mut_annotations['AA_gt_alt'] + mut_annotations['AA'][ maNT_gp[i,j] ]
                                if np.unique(ancnti_gp[i,:])[0] != 4:
                                    mut_annotations['AA_gt_anc'] = mut_annotations['AA_gt_anc'] + mut_annotations['AA'][ np.unique(ancnti_gp[i,:])[0] ]
                # remove duplicates
                mut_annotations['AA_gt_anc'] = ''.join(OrderedDict.fromkeys( mut_annotations['AA_gt_anc'] ).keys()) # get only unique AA and keep order
                mut_annotations['nts'] = ''.join(OrderedDict.fromkeys( mut_annotations['nts'] ).keys()) # get only unique Nuc and keep order
                unique_nts_alt = ''.join(OrderedDict.fromkeys( mut_annotations['nt_alt'] ).keys()) # get only unique Nuc and keep order; keep in separate variable and save in dict just after subsetting for the codons used for AA_gt_alt!
                mut_annotations['AA_gt_ref'] = ''.join(OrderedDict.fromkeys( mut_annotations['AA_gt_ref'] ).keys()) # get only unique AA and keep order
                mut_annotations['AA_gt_alt'] = ''.join([mut_annotations['AA_gt_alt'][mut_annotations['nt_alt'].index(entry)] for entry in unique_nts_alt]) # get only AA encoded from different codons!
                mut_annotations['nt_alt'] = unique_nts_alt
                # record if nonsynonymous/synonymous mutation
                mut_annotations['type'] = ''.join(['S' if aa_alt == refAA else 'N' for aa_alt in mut_annotations['AA_gt_alt']])
                # Record all observed mutations across all isolates; E.g. K134Y, W47A, etc.
                if len(mut_annotations['AA_gt_alt'])>0:
                    mut_annotations['muts'] = [refAA+str(mut_annotations['aa_pos'])+derAA for derAA in mut_annotations['AA_gt_alt'] ]
                    if mut_annotations['muts'] == []:
                        mut_annotations['muts'] = '.'
                else:
                    mut_annotations['muts'] = "."
                mut_annotations['NonSyn_rel_ref'] = [(altAA != mut_annotations['AA_gt_ref']) for altAA in mut_annotations['AA_gt_alt']]
                if len(mut_annotations['AA_gt_anc'])>1:
                    ancestralAA = mut_annotations['AA_gt_anc'][0]
                    derivedAAs = mut_annotations['AA_gt_anc'][1:]
                    mut_annotations['muts_anc'] = [ancestralAA+str(mut_annotations['aa_pos'])+derivedAA for derivedAA in derivedAAs]
                else:
                    mut_annotations['muts_anc'] = "."
                mut_annotations['NonSyn_rel_anc'] = [(altAA != mut_annotations['AA_gt_anc']) for altAA in mut_annotations['AA_gt_alt']]
                if any(mut_annotations['NonSyn_rel_anc']) or any(mut_annotations['NonSyn_rel_ref']): ## need to change to something different --> NonSyn_rel_anc and NonSyn_rel_ref do not need to be same length!
                    mut_annotations['NonSyn'] = True
                else:
                    mut_annotations['NonSyn'] = False
            else: # it is ncRNA or regulatory region (i.e. rRNA, tRNA, tmRNA, ncRNA, CRISPR, oriT, oriC, regulatory region,...)
                mut_annotations['type'] = 'R'
        else: #intergenic
            if int(mut_genenum[i])>0: # get info for gene prior SNP (if any)
                p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])-1] # both -1 necessary bcs list of df and df 0-based
                mut_annotations['gene1'] = p_anno.loc['gene']
                mut_annotations['locustag1'] = p_anno.loc['locustag']
                if 'orthologtag' in p_anno:
                    mut_annotations['orthologtag1'] = p_anno.loc['orthologtag']
                mut_annotations['product1'] = p_anno.loc['product']
                mut_annotations['distance1'] = mut_annotations['pos'] - p_anno.loc['loc2']
                if p_anno.loc['strand'] == -1:
                    mut_annotations['distance1'] = mut_annotations['distance1'] * -1

            if int(mut_genenum[i]+0.5) <= annotation_genes[p_chr-1].shape[0] and annotation_genes[p_chr-1].shape[1] !=0 and (mut_annotations['pos'] < annotation_genes[p_chr-1]['loc2'].max()): # get info gene after SNP (if any); second conditional to evade empty chr
                idx_next_gene_dwnstream = 0
                while annotation_genes[p_chr-1].loc[int(mut_genenum[i])+idx_next_gene_dwnstream, 'loc1'] < mut_annotations['pos']:
                    idx_next_gene_dwnstream += 1
                p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])+idx_next_gene_dwnstream] # -1 necessary bcs list of df 0-based; gene_id 0-based by we want following

                mut_annotations['gene2'] = p_anno.loc['gene']
                mut_annotations['locustag2'] = p_anno.loc['locustag']
                if 'orthologtag' in p_anno:
                    mut_annotations['orthologtag2'] = p_anno.loc['orthologtag']
                mut_annotations['product2'] = p_anno.loc['product']
                mut_annotations['distance2'] = p_anno.loc['loc1'] - mut_annotations['pos'] +1 # +1 to get correct bp distance
                if p_anno.loc['strand'] == 1:
                    mut_annotations['distance2'] = mut_annotations['distance2'] * -1
            # append_annotations intragenic
            #print(mut_annotations['distance1'])
            if ( 'distance1' in mut_annotations and mut_annotations['distance1'] > (-1*promotersize) and mut_annotations['distance1'] < 0) or ( 'distance2' in mut_annotations and mut_annotations['distance2'] > (-1*promotersize) and mut_annotations['distance2'] < 0):
                mut_annotations['type'] = 'P'
            else:
                mut_annotations['type'] = 'I'
            # define nts (repeat of intragenic code)
            if np.unique(refnti_gp[i,:])[0] != 4: # only if ref defined; [0] ok...see next line comment
                mut_annotations['nt_ref'] = nts[np.unique(refnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                mut_annotations['nts'] = mut_annotations['nt_ref']
                if np.unique(ancnti_gp[i,:])[0] != 4:
                    mut_annotations['nt_anc'] = nts[np.unique(ancnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
            # extract derived genotype(s) across all samples
            for j,callidx in enumerate(calls_gp[i,:]):
                #print(str(j)+" "+str(callidx))
                #fixed mutation
                if hasmutation_gp[i,j] == True:
                    if calls_gp[i,j] < 4:
                        mut_annotations['nts'] = mut_annotations['nts'] + nts[calls_gp[i,j]]
                        mut_annotations['nt_alt'] = nts[calls_gp[i,j]]
            # remove duplicates
            mut_annotations['nts'] = ''.join(OrderedDict.fromkeys( mut_annotations['nts'] ).keys()) # get only unique Nuc and keep order
        lod_mutAnno.append(mut_annotations)
    dataframe_mut = pd.DataFrame(lod_mutAnno)
    return dataframe_mut
"""
def genomic_position_all_v2(anno_genes_ls,genomelength,chrstarts, subselect_to_cds_if_present=False):
    # get initial placeholder for all genes (differentiate by chr and genic/nongenic)  
    # UPDATE: gene_num > 0 means CDS. changed 'tagnumber' assignment in parse_gff()
    # 2024/02 MF: generate dict to assign for each pos all possible genetic regions, if genetic regions overlap
    # 2024/02 MF: select assignment by returning multiple CDS (to select the genetic region based on mutation type order: N > S > R > P > I > U)

    locus_tag_nums = np.ones(genomelength,dtype=float)*0.5 # genome-wide unique CDS tag ('tagnumber'), and intergenic 0.5
    cds_indices = np.ones(genomelength,dtype=float)*0.5 # per chromosome unique; intergenic #chr.5-1; intragenic #gene
    pos_to_genenum_dict = {pos: [] for pos in range(genomelength)} ## store per position and overlapping genetic region: (1) locustag_num: genome-wide unique CDS tag ('tagnumber'), and intergenic 0.5; (2) cds_indices: per chromosome unique , intergenic #chr 0.5-1, intragenic #gene; (3) genetype: genetype from annotation
    
    for i,this_chr_df in enumerate(anno_genes_ls):
        if this_chr_df.empty: #Skip contigs w/o CDS annotation.
            continue
        
        gene_num = this_chr_df['tagnumber'].values # tdl called cds_num but actually it was not always cds. e.g tRNA; .values  turn pd.DF to 1D numpy array
        genestarts = chrstarts[i] + this_chr_df['loc1'].values ; # .values  turn pd.DF to 1D numpy array
        geneends = chrstarts[i] + this_chr_df['loc2'].values ; # .values  turn pd.DF to 1D numpy array
        genetype = this_chr_df['type'].values

        for j in range(len(genestarts)-1): # loop over all but last element
            for genomepos in range(genestarts[j], geneends[j]):
                pos_to_genenum_dict[genomepos].append([gene_num[j], j+1, genetype[j]]) # set gene number, chromosome and gene length
            for genomepos in range(geneends[j],(genestarts[j+1])):
                if pos_to_genenum_dict[genomepos] == []: ## check if intergenic is still possible 
                    pos_to_genenum_dict[genomepos].append([j+1+0.5])# intergenic are #chr+0.5; starts with 0.5 prior chr1 which is set when array created
                    
        # mark positions for last gene of chr
        for genomepos in range(genestarts[-1], geneends[-1]):
            pos_to_genenum_dict[genomepos].append([gene_num[-1], len(gene_num), genetype[-1]]) # set gene number, chromosome and gene length
        # mark trailing cds_intergenic from last gene till end of chromosome (if last chr >> till end of genome)
        if ((i+1) < len(chrstarts)):
            for genomepos in range(geneends[-1],chrstarts[i+1]):
                if pos_to_genenum_dict[genomepos] == []: ## check if intergenic is still possible 
                    pos_to_genenum_dict[genomepos].append([len(gene_num)+0.5])
        else:
            for genomepos in range(geneends[-1],genomelength):
                if pos_to_genenum_dict[genomepos] == []: ## check if intergenic is still possible 
                    pos_to_genenum_dict[genomepos].append([len(gene_num)+0.5])
    if subselect_to_cds_if_present:
        ## if one want to prefer CDS (if some are annotated) for N/S mutations rather then annotate everything
        for pos, annotation in pos_to_genenum_dict.items():
            if len(annotation) > 1:
                where_cds = np.where(np.asanyarray(annotation) == 'CDS')[0] ## get if there are coding sequences --> favour them! 
                if (len(where_cds) > 0) & (len(where_cds) < len(annotation)): ## select longest coding sequence at position (if not all are CDS but at least one --> else directly subset!)
                    pos_to_genenum_dict[pos] = [annotation[idx] for idx in where_cds]
    return pos_to_genenum_dict
"""

def genomic_position_all_v3(anno_genes_ls,genomelength,chrstarts, p_to_subset_to=[]):
    # get initial placeholder for all genes (differentiate by chr and genic/nongenic)  
    # UPDATE: gene_num > 0 means CDS. changed 'tagnumber' assignment in parse_gff()
    # 2024/02 MF: generate dict to assign for each pos all possible genetic regions, if genetic regions overlap
    # 2024/02 MF: select assignment by returning multiple CDS (to select the genetic region based on mutation type order: premature stop > N > S > R > P > I > U)

    pos_to_genenum_l = [[] for pos in range(genomelength)] ## store per position and overlapping genetic region: (1) locustag_num: genome-wide unique CDS tag ('tagnumber'), and intergenic 0.5; (2) cds_indices: per chromosome unique , intergenic #chr 0.5-1, intragenic #gene; (3) genetype: genetype from annotation
    stop_loop = False
    for i,this_chr_df in enumerate(anno_genes_ls):
        if this_chr_df.empty: #Skip contigs w/o CDS annotation.
            continue
        gene_num = this_chr_df['tagnumber'].values # tdl called cds_num but actually it was not always cds. e.g tRNA; .values  turn pd.DF to 1D numpy array
        genestarts = chrstarts[i] + this_chr_df['loc1'].values ; # .values  turn pd.DF to 1D numpy array
        geneends = chrstarts[i] + this_chr_df['loc2'].values ; # .values  turn pd.DF to 1D numpy array

        for j in range(len(genestarts)-1): # loop over all but last element
            for genomepos in range(genestarts[j], geneends[j]):
                pos_to_genenum_l[genomepos].append([gene_num[j], j+1]) # set gene number, chromosome and gene length
            for genomepos in range(geneends[j],(genestarts[j+1])):
                if pos_to_genenum_l[genomepos] == []: ## check if intergenic is still possible 
                    pos_to_genenum_l[genomepos].append([0.5, j+1+0.5])# intergenic are #chr+0.5; starts with 0.5 prior chr1 which is set when array created        
        # mark positions for last gene of chr
        print_exceeding_annotation = True
        for genomepos in range(genestarts[-1], geneends[-1]):
            if ((i+1) < len(chrstarts)):
                pos_in_chr_range = (genomepos < chrstarts[i+1])  ## pos on any non-last chr and not exceeding chr boundaries
            else:
                pos_in_chr_range = (genomepos < genomelength) ## pos on last chr and not exceeding genome boundaries
            if pos_in_chr_range: ## note some refseq annotations can exceed the fasta sequence! --> we map to the fasta --> we just need the region until the sequence ends and cannot do much with the info after
                pos_to_genenum_l[genomepos].append([gene_num[-1], len(gene_num)]) # set gene number, chromosome and gene length
            elif print_exceeding_annotation:
                print(f'Regions of annotation exceed contig {i+1}')
                print_exceeding_annotation = False
        # mark trailing cds_intergenic from last gene till end of chromosome (if last chr >> till end of genome)
        if ((i+1) < len(chrstarts)):
            for genomepos in range(geneends[-1],chrstarts[i+1]):
                if pos_to_genenum_l[genomepos] == []: ## check if intergenic is still possible 
                    pos_to_genenum_l[genomepos].append([0.5, len(gene_num)+0.5])
        else:
            for genomepos in range(geneends[-1],genomelength):
                if pos_to_genenum_l[genomepos] == []: ## check if intergenic is still possible 
                    pos_to_genenum_l[genomepos].append([0.5, len(gene_num)+0.5])
    if len(p_to_subset_to) > 0: 
        ## subset to positions of interest (e.g. goodpos)
        pos_to_genenum_l = np.array(pos_to_genenum_l, dtype=object)[p_to_subset_to]
    ## replace all empty lists with 0.5 --> intergenic 
    pos_to_genenum_l = [entries if entries != [] else [[0.5, 0.5]] for entries in pos_to_genenum_l]
    return pos_to_genenum_l


def extract_nth_instance(row, n, columns):
    """
    Function to return the nth letter of a string from multiple columns of a pandas df 
    """
    pdseries_dict = {}
    for col in columns:
        if row[col] != row[col]: ## echek for nans
            pdseries_dict[col] = np.nan
        elif n < len(row[col]):
            if isinstance(row[col], list):
                pdseries_dict[col] = [row[col][n]]
            else:
                pdseries_dict[col] = str(row[col][n])
    return pd.Series(pdseries_dict) # pd.Series({col: str(row[col][n-1]) if () else '' for col in columns})

def deduplicate_annotation_mutations(annotation_mutations, timestamp = '', association_order = ['N', 'S', 'R', 'P', 'I', 'U'], columns_to_reduce = ['nt_alt', 'AA_gt_alt', 'type', 'muts', 'NonSyn_rel_ref', 'NonSyn_rel_anc']):
    """
    Deduplicate annotation mutations if on a given position mutliple genes overlap to annotation_mutations df with one genetic region per position and allele
    """
    association_order_dict = {mut_type: idx for idx, mut_type in enumerate(reversed(association_order))}
    unique_pos_anno_mut = annotation_mutations[~annotation_mutations[['chr', 'pos']].duplicated(keep = False)].copy()
    dudup_p_anno_mut = annotation_mutations[annotation_mutations[['chr', 'pos']].duplicated(keep = False)].copy()
    
    removed_duplicates = pd.DataFrame()
    if dudup_p_anno_mut.empty:
        print('All positions are unique. No overlapping genetic regions contain any SNV')
        return unique_pos_anno_mut
    else:
        ## check if all columns are present
        if not all([col in dudup_p_anno_mut.columns for col in columns_to_reduce]):
            columns_to_reduce = [col for col in columns_to_reduce if col in dudup_p_anno_mut.columns]
            print('Not all columns which should be reduced based on deduplication are present')
            print(f'Proceeding with {"; ".join([col for col in columns_to_reduce])}')
        ## translate mutation types 
        dudup_p_anno_mut['type_prio'] = dudup_p_anno_mut['type'].apply(lambda type_col: [association_order_dict[mut_type] for mut_type in type_col])
        ## loop over each position 
        for chr, pos in dudup_p_anno_mut[['chr', 'pos']].drop_duplicates().values:
            dedup_p_anno_mut_sub = dudup_p_anno_mut[(dudup_p_anno_mut['chr'] == chr) & (dudup_p_anno_mut['pos'] == pos)]
            num_mutation_types = dedup_p_anno_mut_sub['type_prio'].str.len()
            max_num_mutation_types_idx = np.where(num_mutation_types == max(num_mutation_types))[0]
            non_max_num_mutation_types_idx = np.where(num_mutation_types != max(num_mutation_types))[0]
            ## subset to entries with the maximal number of mutation types --> just higher prioritized mutation types (N, S) should have multiple mutation types 
            dedup_p_anno_mut_sub_prio = dedup_p_anno_mut_sub.iloc[max_num_mutation_types_idx].reset_index(drop = True)
            removed_duplicates = pd.concat([removed_duplicates, dedup_p_anno_mut_sub.iloc[non_max_num_mutation_types_idx]]) ## save removed row
            if sum(num_mutation_types > 1) == 1: ## only one entry has more than 1 mutational type --> select
                unique_pos_anno_mut = pd.concat([unique_pos_anno_mut, dedup_p_anno_mut_sub_prio]) ## save mutation with higher prio to anno mutations
            else: ## multiple entries have same length of mutational type --> select higher prios
                dedup_p_anno_mut_sub_prio_clean = dedup_p_anno_mut_sub_prio.copy()
                dedup_p_anno_mut_sub_prio_toremove_clean = dedup_p_anno_mut_sub_prio.copy()
                for col in columns_to_reduce:
                    dedup_p_anno_mut_sub_prio_clean[col] = dedup_p_anno_mut_sub_prio_clean[col].apply(lambda x: [] if isinstance(x, list) else '')
                    dedup_p_anno_mut_sub_prio_toremove_clean[col] = dedup_p_anno_mut_sub_prio_toremove_clean[col].apply(lambda x: [] if isinstance(x, list) else '')
                for i in range(max(num_mutation_types)): ## loop over every mutation type and select the one with the higher association
                    highest_prio = max(dedup_p_anno_mut_sub_prio['type_prio'].str[i])
                    highest_prio_idx = np.where(dedup_p_anno_mut_sub_prio['type_prio'].str[i] == highest_prio)[0]
                    low_prio_idx = np.where(dedup_p_anno_mut_sub_prio['type_prio'].str[i] != highest_prio)[0]
                    if len(highest_prio_idx) > 1: ## multiple genetic regions have same mutational type priority 
                        print(f'Genes involved: {"; ".join(dedup_p_anno_mut_sub_prio.loc[highest_prio_idx, "gene"].values)}')
                        if all(dedup_p_anno_mut_sub_prio.loc[highest_prio_idx, 'type'].str[i] == 'N'): ## we have just nonsyn mutations --> try to select premature stops over other nonsyn mutations
                            highest_prio_premature_stop = np.where(dedup_p_anno_mut_sub_prio.loc[highest_prio_idx, 'muts'].str[i].str[-1] == '*')[0]
                            non_premature_stop = np.where(dedup_p_anno_mut_sub_prio.loc[highest_prio_idx, 'muts'].str[i].str[-1] != '*')[0]
                            non_premature_stop = np.append(non_premature_stop, values = low_prio_idx) ## append non-premature stop idx to low_prios 
                            if len(highest_prio_premature_stop) == 0: ## no premature stop --> select first instance of highest priority
                                low_prio_idx = np.append(low_prio_idx, values = highest_prio_idx[1:]) ## append non first highest prio to low prio idx 
                                highest_prio_idx = highest_prio_idx[0]
                            elif len(highest_prio_premature_stop) > 0: ## select gene which might have a premature stop at first (if multiple, select the first of multiple!) 
                                low_prio_idx = np.append(non_premature_stop, values = highest_prio_premature_stop[1:]) ## append non first highest prio to low prio idx 
                                highest_prio_idx = highest_prio_premature_stop[0]
                        else: ## multiple genes have same mutation type and are no nonsyn mutations --> select the first one 
                            low_prio_idx = np.append(low_prio_idx, values = highest_prio_idx[1:]) ## append non first highest prio to low prio idx 
                            highest_prio_idx = highest_prio_idx[0]
                    else: ## just one genetic region has a higher priority
                        highest_prio_idx = highest_prio_idx[0]
                    ## store based on associations 
                    if not isinstance(highest_prio_idx, int):
                        SystemError(f'Multiple highest priority genes would have been selected for the deduplicated annotation mutation df.\nPlease check code for chromosome {chr}, position {pos}, allele id {str(i+1)}')
                    dedup_p_anno_mut_sub_prio_clean.loc[highest_prio_idx, columns_to_reduce] = dedup_p_anno_mut_sub_prio_clean.loc[highest_prio_idx, columns_to_reduce] + extract_nth_instance(dedup_p_anno_mut_sub_prio.iloc[highest_prio_idx], i, columns_to_reduce)
                    for idx in low_prio_idx: ## loop over each 
                        dedup_p_anno_mut_sub_prio_toremove_clean.loc[idx, columns_to_reduce] = dedup_p_anno_mut_sub_prio_toremove_clean.loc[idx, columns_to_reduce] + extract_nth_instance(dedup_p_anno_mut_sub_prio.iloc[idx], i, columns_to_reduce)
                dedup_p_anno_mut_sub_prio_clean = dedup_p_anno_mut_sub_prio_clean[(dedup_p_anno_mut_sub_prio_clean[columns_to_reduce] != '').all(axis = 1)] ## remove empty lines 
                dedup_p_anno_mut_sub_prio_toremove_clean = dedup_p_anno_mut_sub_prio_toremove_clean[(dedup_p_anno_mut_sub_prio_toremove_clean[columns_to_reduce] != '').all(axis = 1)] ## remove empty lines 
                ## modify the nts column to keep propriate amount of nts per gene identified
                dedup_p_anno_mut_sub_prio_clean['nts'] = dedup_p_anno_mut_sub_prio_clean['nt_ref'] + dedup_p_anno_mut_sub_prio_clean['nt_alt']
                dedup_p_anno_mut_sub_prio_toremove_clean['nts'] = dedup_p_anno_mut_sub_prio_toremove_clean['nt_ref'] + dedup_p_anno_mut_sub_prio_toremove_clean['nt_alt']
                unique_pos_anno_mut = pd.concat([unique_pos_anno_mut, dedup_p_anno_mut_sub_prio_clean]) ## save mutation with higher prio to anno mutations
                removed_duplicates = pd.concat([removed_duplicates, dedup_p_anno_mut_sub_prio_toremove_clean]) ## save removed row
        timestamp = timestamp + '_'
        removed_duplicates.to_csv(f'{timestamp}annotation_mutations_deduplicated_unselected_regions_pos_alleles.csv')
    unique_pos_anno_mut = unique_pos_anno_mut.sort_values(['chr', 'pos']).reset_index(drop = True)
    return unique_pos_anno_mut


def annotate_mutations_v2(annotation_genes , p_gp , refnti_gp , ancnti_gp , calls_gp , counts_gp , hasmutation_gp , mutQual_gp, promotersize, ref_genome_folder, goodpos2use=[]):
    ''' produce dataframe with annotation info for each goodposition used for tree
    all positions are 1-based! (pos, nt_pos, loc1/2, distance1/2)
    # function combines TDL: annotate_mutations_gb.m and append_annotations.m 
    CHANGELOG:
    17.10.2022 IL: Added parsing of ancestral AA state by getting annotation info from ancnti_gp
    02/2024 MF: changed AA gt return to AA ref, anc and alt and evaluation of nonsyn mutations
    03/2024 MF: return annotation mutations as a dataframe with all possible associations (==> Some positions are duplicated!); selection to unique positions and alleles can be done via deduplicate_annotation_mutations()
    '''
    # ref_genome_folder = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/metadata/reference_genomes/P07_Ehormaechei-c1_231009/')
    # p_gp: genomic position of goodpos
    # p_gp=p[goodpos2use]
    # refnti_gp=refnti_m[goodpos2use,:]
    # calls_gp=calls[goodpos2use,:]
    # counts_gp=counts[:,:,goodpos2use]
    # hasmutation_gp=hasmutation[goodpos2use,:]
    # mutQual_gp = mutQual[goodpos2use,].flatten() 
    # ref_genome_folder=REFGENOMEFOLDER

    if len(goodpos2use) > 0:
        p_gp = p_gp[goodpos2use]
        refnti_gp = refnti_gp[goodpos2use,:]
        ancnti_gp = ancnti_gp[goodpos2use,:]
        calls_gp = calls_gp[goodpos2use,:]
        counts_gp = counts_gp[:,:,goodpos2use]
        hasmutation_gp = hasmutation_gp[goodpos2use,:]
        mutQual_gp = mutQual_gp[goodpos2use,].flatten()

    [maf_gp, maNT_gp, minorNT_gp, minorAF_gp] = div_major_allele_freq(counts_gp)
    nts = ['A','T','C','G'] #'atcg';
    rc = ['T','A','G','C']  # 'tagc';
    [chrstarts, genomelength,scafnames]= genomestats(ref_genome_folder);
    gp_pos_to_genenum_l = genomic_position_all_v3(annotation_genes,genomelength,chrstarts, p_to_subset_to=p_gp); #loc_tag numbers all genes across genome (intergenic:0.5); chr_tag: numbers genes per chr (intergenic: #chr+0.5)
    # mut_genenum = chr_tag[p_gp]
    # loc_tag_p_gp = locus_tag[p_gp]
    pos_chr = p2chrpos(p_gp,chrstarts) # get two-col table chr,pos   
    lod_mutAnno = []
    for i,annotations in enumerate(gp_pos_to_genenum_l): ## iterate over goodpos positions which are associated to none, one or multiple genetic regions
        for p_gene_assoc in annotations: ## iterate over all genetic regions associated with current position 
            mut_annotations = {}
            mut_annotations['gene_num'] = p_gene_assoc[1]
            mut_annotations['gene_num_global'] = p_gene_assoc[0]
            mut_annotations['chr'] = pos_chr[i][0]
            mut_annotations['pos'] = pos_chr[i][1] + 1 # turn 0-based pos into 1-based pos.
            mut_annotations['p'] = p_gp[i]
            mut_annotations['quals'] = mutQual_gp[i]
            mut_annotations['nts'] = "." # default. overwritten below if NT defined.
            mut_annotations['nt_ref'] = "." # default. overwritten below if NT defined.
            mut_annotations['nt_anc'] = "." # default. overwritten below if NT defined.
            mut_annotations['nt_alt'] = "." # default. overwritten below if NT defined.
            p_chr = pos_chr[i][0]
            if mut_annotations['gene_num'] == int(mut_annotations['gene_num']): # intragenic
                p_anno = annotation_genes[p_chr-1].iloc[int(mut_annotations['gene_num'])-1] # both -1 necessary bcs list of df and df 0-based
                mut_annotations['product'] = p_anno.loc['product']
                mut_annotations['gene'] = p_anno.loc['gene']
                mut_annotations['protein_id'] = p_anno.loc['protein_id']
                mut_annotations['strand'] = p_anno.loc['strand']
                mut_annotations['loc1'] = p_anno.loc['loc1'] + 1 # 1-based first position
                mut_annotations['loc2'] = p_anno.loc['loc2'] # last position of gene (inclusive)
                mut_annotations['sequence'] = p_anno.loc['sequence']
                mut_annotations['protein_id'] = p_anno.loc['protein_id']
                mut_annotations['note'] = p_anno.loc['note']
                mut_annotations['locustag'] = p_anno.loc['locustag']
                if 'orthologtag' in p_anno:
                    mut_annotations['orthologtag'] = p_anno.loc['orthologtag']
                if p_anno.loc['strand'] == 1: # get position within gene, consider strandedness
                    mut_annotations['nt_pos'] = mut_annotations['pos'] - mut_annotations['loc1']+1 # pos/loc1 1-based. nt_pos 1-based. +1 to get 1-based nt_pos in gene (checked)!
                elif p_anno.loc['strand'] == -1:
                    mut_annotations['nt_pos'] = p_anno.loc['loc2'] - mut_annotations['pos'] +1 # pos/loc2 1-based. +1 to get 1-based nt_pos in gene (checked)!
                else:
                    mut_annotations['nt_pos'] = "." # can happen. eg. Crispr
                
                if p_anno.loc['type'] == 'CDS': ## check if CDS --> otherwise it is ncRNA or regulatory region (i.e. rRNA, tRNA, tmRNA, ncRNA, CRISPR, oriT, oriC, regulatory region,...)
                    mut_annotations['translation'] = p_anno.loc['translation']
                    if mut_annotations['nt_pos'] == ".": 
                        aan = 9999999
                    else:
                        aan = int( (mut_annotations['nt_pos']-1 )/3 ) + 1 # get #codon that harbours mutation. 1-based.
                        ncn = mut_annotations['nt_pos'] - ((aan-1)*3) # get #nucl within triplett. 1-based
                        mut_annotations['aa_pos'] = aan
                    codons_ls = []
                    aa_ls = [] #potential AA changes at given codon
                    if len(mut_annotations['sequence']) >= (aan*3) and mut_annotations['translation'] != ".":
                        codon = mut_annotations['sequence'][aan*3-2-1:aan*3] # -1 bcs seq-object 0-based; but positional argument aan not 
                        codon = [n for n in codon  ] # turn seq.object to list, seq object does not allow reassignment
                        for n in range(4): # test all four nucleotide options for resulting AA change
                            if p_anno.loc['strand'] == 1:
                                codon[ncn-1] = nts[n]
                            else:
                                codon[ncn-1] = rc[n]
                            codonSeqO = Seq( "".join(codon)) 
                            codons_ls.append(codonSeqO)
                            aa_ls.append(codonSeqO.translate())
                    mut_annotations['codons'] = codons_ls
                    mut_annotations['AA'] = [aa[0] for aa in aa_ls]
                    # append_annotations intragenic
                    if len(mut_annotations['AA']) < 4:
                        mut_annotations['type'] = 'U'
                    mut_annotations['AA_gt_ref'] = '.' # Fill in annotations with whether NS or Syn mutation
                    mut_annotations['AA_gt_anc'] = '' 
                    mut_annotations['nt_alt'] = ''
                    mut_annotations['AA_gt_alt'] = ''
                    if np.unique(refnti_gp[i,:])[0] != 4: # only if ref defined; [0] ok...see next line comment
                        mut_annotations['nt_ref'] = nts[np.unique(refnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                        mut_annotations['nts'] = mut_annotations['nt_ref']
                    if np.unique(ancnti_gp[i,:])[0] != 4:
                        mut_annotations['nt_anc'] = nts[np.unique(ancnti_gp[i,:])[0]] # ancnti_gp should by definition never have > 1 unique element per row (one outgroup only!!! >> [0] apropriate
                    if len(mut_annotations['AA']) == 4:
                        refAA = mut_annotations['AA'][ np.unique(refnti_gp[i,:])[0] ]
                        mut_annotations['AA_gt_ref'] = refAA # the AA list order corresponds to NTS list!
                    # extract derived genotype(s) and according AA across all samples
                    # call Nonsynonymous if >1 amino acid type exists across all called samples (eg ref AA and alt )
                    for j,_ in enumerate(calls_gp[i,:]):
                        #fixed mutation
                        # i is index of current pos, j is index of samples 
                        if hasmutation_gp[i,j] == True :
                            if calls_gp[i,j] < 4:
                                # add each NT called, later only keep uniq, but kept order
                                mut_annotations['nts'] = mut_annotations['nts'] + nts[calls_gp[i,j]]
                                mut_annotations['nt_alt'] = mut_annotations['nt_alt'] + nts[calls_gp[i,j]]
                                if len(mut_annotations['AA']) == 4:
                                    # mut_annotations['AA_gt'] = mut_annotations['AA_gt'] + mut_annotations['AA'][ maNT_gp[i,j] ]
                                    mut_annotations['AA_gt_alt'] = mut_annotations['AA_gt_alt'] + mut_annotations['AA'][ maNT_gp[i,j] ]
                                    if np.unique(ancnti_gp[i,:])[0] != 4:
                                        mut_annotations['AA_gt_anc'] = mut_annotations['AA_gt_anc'] + mut_annotations['AA'][ np.unique(ancnti_gp[i,:])[0] ]
                    # remove duplicates
                    mut_annotations['AA_gt_anc'] = ''.join(OrderedDict.fromkeys( mut_annotations['AA_gt_anc'] ).keys()) # get only unique AA and keep order
                    mut_annotations['nts'] = ''.join(OrderedDict.fromkeys( mut_annotations['nts'] ).keys()) # get only unique Nuc and keep order
                    unique_nts_alt = ''.join(OrderedDict.fromkeys( mut_annotations['nt_alt'] ).keys()) # get only unique Nuc and keep order; keep in separate variable and save in dict just after subsetting for the codons used for AA_gt_alt!
                    mut_annotations['AA_gt_ref'] = ''.join(OrderedDict.fromkeys( mut_annotations['AA_gt_ref'] ).keys()) # get only unique AA and keep order
                    mut_annotations['AA_gt_alt'] = ''.join([mut_annotations['AA_gt_alt'][mut_annotations['nt_alt'].index(entry)] for entry in unique_nts_alt]) # get only AA encoded from different codons!
                    mut_annotations['nt_alt'] = unique_nts_alt
                    # record if nonsynonymous/synonymous mutation
                    mut_annotations['type'] = ''.join(['S' if aa_alt == refAA else 'N' for aa_alt in mut_annotations['AA_gt_alt']])
                    # Record all observed mutations across all isolates; E.g. K134Y, W47A, etc.
                    if len(mut_annotations['AA_gt_alt'])>0:
                        mut_annotations['muts'] = [refAA+str(mut_annotations['aa_pos'])+derAA if (derAA != refAA) else '.' for derAA in mut_annotations['AA_gt_alt']]
                        if mut_annotations['muts'] == []:
                            mut_annotations['muts'] = ['.']
                    else:
                        mut_annotations['muts'] = ["."]
                    mut_annotations['NonSyn_rel_ref'] = [(altAA != mut_annotations['AA_gt_ref']) for altAA in mut_annotations['AA_gt_alt']]
                    if len(mut_annotations['AA_gt_anc'])>1:
                        ancestralAA = mut_annotations['AA_gt_anc'][0]
                        derivedAAs = mut_annotations['AA_gt_anc'][1:]
                        mut_annotations['muts_anc'] = [ancestralAA+str(mut_annotations['aa_pos'])+derivedAA if (derivedAA != ancestralAA) else '.' for derivedAA in derivedAAs]
                    else:
                        mut_annotations['muts_anc'] = ["."]
                    mut_annotations['NonSyn_rel_anc'] = [(altAA != mut_annotations['AA_gt_anc']) for altAA in mut_annotations['AA_gt_alt']]
                    if any(mut_annotations['NonSyn_rel_anc']) or any(mut_annotations['NonSyn_rel_ref']): ## need to change to something different --> NonSyn_rel_anc and NonSyn_rel_ref do not need to be same length!
                        mut_annotations['NonSyn'] = True
                    else:
                        mut_annotations['NonSyn'] = False
                else: # it is ncRNA or regulatory region (i.e. rRNA, tRNA, tmRNA, ncRNA, CRISPR, oriT, oriC, regulatory region,...)
                    mut_annotations['nt_alt'] = ''
                    if np.unique(refnti_gp[i,:])[0] != 4: # only if ref defined; [0] ok...see next line comment
                        mut_annotations['nt_ref'] = nts[np.unique(refnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                        mut_annotations['nts'] = mut_annotations['nt_ref']
                    if np.unique(ancnti_gp[i,:])[0] != 4:
                        mut_annotations['nt_anc'] = nts[np.unique(ancnti_gp[i,:])[0]] # ancnti_gp should by definition never have > 1 unique element per row (one outgroup only!!! >> [0] apropriate
                    # extract derived genotype(s)
                    for j,_ in enumerate(calls_gp[i,:]):
                        #fixed mutation
                        # i is index of current pos, j is index of samples 
                        if hasmutation_gp[i,j] == True :
                            if calls_gp[i,j] < 4:
                                # add each NT called, later only keep uniq, but kept order
                                mut_annotations['nts'] = mut_annotations['nts'] + nts[calls_gp[i,j]]
                                mut_annotations['nt_alt'] = mut_annotations['nt_alt'] + nts[calls_gp[i,j]]
                    # remove duplicates
                    mut_annotations['nts'] = ''.join(OrderedDict.fromkeys( mut_annotations['nts'] ).keys()) # get only unique Nuc and keep order
                    mut_annotations['nt_alt'] = ''.join(OrderedDict.fromkeys( mut_annotations['nt_alt'] ).keys()) # get only unique Nuc and keep order
                    mut_annotations['type'] = 'R'
            else: #intergenic
                if int(mut_annotations['gene_num'])>0: # get info for gene prior SNP (if any)
                    p_anno = annotation_genes[p_chr-1].iloc[int(mut_annotations['gene_num'])-1] # both -1 necessary bcs list of df and df 0-based
                    mut_annotations['gene1'] = p_anno.loc['gene']
                    mut_annotations['locustag1'] = p_anno.loc['locustag']
                    if 'orthologtag' in p_anno:
                        mut_annotations['orthologtag1'] = p_anno.loc['orthologtag']
                    mut_annotations['product1'] = p_anno.loc['product']
                    mut_annotations['distance1'] = mut_annotations['pos'] - p_anno.loc['loc2']
                    if p_anno.loc['strand'] == -1:
                        mut_annotations['distance1'] = mut_annotations['distance1'] * -1

                if (int(mut_annotations['gene_num']+0.5) <= annotation_genes[p_chr-1].shape[0]) and (annotation_genes[p_chr-1].shape[1] !=0) and (mut_annotations['pos'] < annotation_genes[p_chr-1]['loc2'].max()): # get info gene after SNP (if any); second conditional to evade empty chr
                    idx_next_gene_dwnstream = 0
                    while annotation_genes[p_chr-1].loc[int(mut_annotations['gene_num'])+idx_next_gene_dwnstream, 'loc1'] < mut_annotations['pos']:
                        idx_next_gene_dwnstream += 1
                    p_anno = annotation_genes[p_chr-1].iloc[int(mut_annotations['gene_num'])+idx_next_gene_dwnstream] # -1 necessary bcs list of df 0-based; gene_id 0-based by we want following

                    mut_annotations['gene2'] = p_anno.loc['gene']
                    mut_annotations['locustag2'] = p_anno.loc['locustag']
                    if 'orthologtag' in p_anno:
                        mut_annotations['orthologtag2'] = p_anno.loc['orthologtag']
                    mut_annotations['product2'] = p_anno.loc['product']
                    mut_annotations['distance2'] = p_anno.loc['loc1'] - mut_annotations['pos'] +1 # +1 to get correct bp distance
                    if p_anno.loc['strand'] == 1:
                        mut_annotations['distance2'] = mut_annotations['distance2'] * -1
                # append_annotations intragenic
                if ( 'distance1' in mut_annotations and mut_annotations['distance1'] > (-1*promotersize) and mut_annotations['distance1'] < 0) or ( 'distance2' in mut_annotations and mut_annotations['distance2'] > (-1*promotersize) and mut_annotations['distance2'] < 0):
                    mut_annotations['type'] = 'P'
                else:
                    mut_annotations['type'] = 'I'
                # define nts (repeat of intragenic code)
                if np.unique(refnti_gp[i,:])[0] != 4: # only if ref defined; [0] ok...see next line comment
                    mut_annotations['nt_ref'] = nts[np.unique(refnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                    mut_annotations['nts'] = mut_annotations['nt_ref']
                if np.unique(ancnti_gp[i,:])[0] != 4:
                    mut_annotations['nt_anc'] = nts[np.unique(ancnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                # extract derived genotype(s) across all samples
                mut_annotations['nt_alt'] = ''
                for j,callidx in enumerate(calls_gp[i,:]):
                    #fixed mutation
                    if hasmutation_gp[i,j] == True:
                        if calls_gp[i,j] < 4:
                            mut_annotations['nts'] = mut_annotations['nts'] + nts[calls_gp[i,j]]
                            mut_annotations['nt_alt'] = mut_annotations['nt_alt'] + nts[calls_gp[i,j]]
                # remove duplicates
                mut_annotations['nts'] = ''.join(OrderedDict.fromkeys( mut_annotations['nts'] ).keys()) # get only unique Nuc and keep order
                mut_annotations['nt_alt'] = ''.join(OrderedDict.fromkeys( mut_annotations['nt_alt'] ).keys()) # get only unique Nuc and keep order
            lod_mutAnno.append(mut_annotations)
    dataframe_mut = pd.DataFrame(lod_mutAnno)
    return dataframe_mut

def repolarize_snvs(anno_table):
    ## function to repolarize samples based on ancestral sequence 
    ## NOTE: This is a helper function until annotation mutations is cleaned as the anc annotation is not working as expected!
    ## NOTE: To polarize based on anc, one should repolarize hasmutations as well!
    if any(anno_table['nt_anc'] == '.'):
        print('repolarization impossible due to incompletely resolved ancestral sequence. Please provide full ancestral sequence!')
        return 
    wrong_pol = anno_table['nt_ref'] != anno_table['nt_anc']
    for rowidx, row in anno_table[wrong_pol].iterrows():
        anno_table.loc[rowidx, 'nt_alt'] = anno_table.loc[rowidx, 'nt_alt'].replace(row['nt_anc'], row['nt_ref'])
        if (('N' in row['type']) or ('S' in row['type'])) and (row['nt_anc'] in 'ATCG'):
            codons = np.array(row['codons'])
            aminoacids = np.array(row['AA'])
            pos_changing = [i for i in range(3) if any([s[i] != codons[0][i] for s in codons])][0]
            if row['strand'] == -1:
                anc_codon = aminoacids[[str(Seq(row['nt_anc']).reverse_complement()) == codon[pos_changing] for codon in codons]][0]
                alt_codons = [aminoacids[[str(Seq(nt_alt).reverse_complement()) == codon[pos_changing] for codon in codons]][0] for nt_alt in row['nt_alt']]
            else:
                anc_codon = aminoacids[[row['nt_anc'] == codon[pos_changing] for codon in codons]][0]
                alt_codons = [aminoacids[[nt_alt == codon[pos_changing] for codon in codons]][0] for nt_alt in row['nt_alt']]
            polarized_muts = [f'{anc_codon}{int(row["aa_pos"])}{alt_codon}' if (alt_codon != anc_codon) else '.' for alt_codon in alt_codons]
            anno_table.loc[rowidx, 'muts'] = ', '.join(polarized_muts)
            anno_table.loc[rowidx, 'type'] = ''.join(['N' if mut != '.' else 'S' for mut in polarized_muts])
    return anno_table

def annotate_mutations_single_sample_vs_custom_nt_call(annotation_genes , p_goodpos , custom_nt_calls_goodpos , calls_single_sample_goodpos , hasmutation_goodpos , promotersize, ref_genome_folder, goodpos=[]):
    ''' 
    produces annotation_mutations-like dataframe with annotation info for each goodposition used for tree for a given sample
    Comparison and type called based on "calls" from user provided goodpos
    all positions are 1-based! (pos, nt_pos, loc1/2, distance1/2)
    CHANGELOG:
    17.10.2022 IL: Added parsing of ancestral AA state by getting annotation info from ancnti_gp
    23.03.2023 IL: Removed reference nt submission for clarity, all reliance on provided custom_nt_call call'''
    # p_gp: genomic position of goodpos
    # p_gp=p[goodpos2use]
    # refnti_gp=refnti_m[goodpos2use,:]
    # calls_gp=calls[goodpos2use,:]
    # counts_gp=counts[:,:,goodpos2use]
    # hasmutation_gp=hasmutation[goodpos2use,:]
    # mutQual_gp = mutQual[goodpos2use,].flatten() 
    # ref_genome_folder=REFGENOMEFOLDER
    if len(goodpos) > 0:
        p_goodpos = p_goodpos[goodpos]
        custom_nt_calls_goodpos = custom_nt_calls_goodpos[goodpos,:]
        calls_single_sample_goodpos = calls_single_sample_goodpos[goodpos,:]
        hasmutation_goodpos = hasmutation_goodpos[goodpos,:]
    nts = ['A','T','C','G'] #'atcg';
    rc = ['T','A','G','C']  # 'tagc';
    [chrstarts, genomelength,scafnames]= genomestats(ref_genome_folder);
    [locus_tag,chr_tag] = genomic_position_all(annotation_genes, genomelength, chrstarts); #loc_tag numbers all genes across genome (intergenic:0.5); chr_tag: numbers genes per chr (intergenic: #chr+0.5)
    mut_genenum = chr_tag[p_goodpos]
    loc_tag_p_gp = locus_tag[p_goodpos]
    pos_chr = p2chrpos(p_goodpos,chrstarts) # get two-col table chr,pos   
    lod_mutAnno = []
    for i,pos in enumerate(p_goodpos):
        #print(i,pos)
        mut_annotations = {}
        mut_annotations['gene_num'] = mut_genenum[i]
        mut_annotations['gene_num_global'] = loc_tag_p_gp[i]
        mut_annotations['chr'] = pos_chr[i][0]
        mut_annotations['pos'] = pos_chr[i][1] + 1 # turn 0-based pos into 1-based pos.
        mut_annotations['p'] = p_goodpos[i]
        mut_annotations['nts'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_ref'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_anc'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_alt'] = "." # default. overwritten below if NT defined.
        p_chr = pos_chr[i][0]
        if mut_genenum[i] == int(mut_genenum[i]): # intragenic
            p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])-1] # both -1 necessary bcs list of df and df 0-based
            mut_annotations['product'] = p_anno.loc['product']
            mut_annotations['gene'] = p_anno.loc['gene']
            mut_annotations['protein_id'] = p_anno.loc['protein_id']
            mut_annotations['strand'] = p_anno.loc['strand']
            mut_annotations['loc1'] = p_anno.loc['loc1'] + 1 # 1-based first position
            mut_annotations['loc2'] = p_anno.loc['loc2'] # last position of gene (inclusive)
            mut_annotations['sequence'] = p_anno.loc['sequence']
            mut_annotations['protein_id'] = p_anno.loc['protein_id']
            mut_annotations['note'] = p_anno.loc['note']
            mut_annotations['locustag'] = p_anno.loc['locustag']
            mut_annotations['chr_locustag'] = str(pos_chr[i][0]) + '_' + p_anno.loc['locustag']
            if 'orthologtag' in p_anno:
                mut_annotations['orthologtag'] = p_anno.loc['orthologtag']
            mut_annotations['translation'] = p_anno.loc['translation']
            if p_anno.loc['strand'] == 1: # get position within gene, consider strandedness
                mut_annotations['nt_pos'] = mut_annotations['pos'] - mut_annotations['loc1']+1 # pos/loc1 1-based. nt_pos 1-based. +1 to get 1-based nt_pos in gene (checked)!
            elif p_anno.loc['strand'] == -1:
                mut_annotations['nt_pos'] = p_anno.loc['loc2'] - mut_annotations['pos'] +1 # pos/loc2 1-based. +1 to get 1-based nt_pos in gene (checked)!
            else:
                mut_annotations['nt_pos'] = "." # can happen. eg. Crispr
            if mut_annotations['nt_pos'] == ".": # Observed with special 'type's like crispr annotated. rare! leads to no AA inference.
                aan = 9999999
            else:
                aan = int( (mut_annotations['nt_pos']-1 )/3 ) + 1 # get #codon that harbours mutation. 1-based.
                ncn = mut_annotations['nt_pos'] - ((aan-1)*3) # get #nucl within triplett. 1-based
                mut_annotations['aa_pos'] = aan
            codons_ls = []
            aa_ls = []
            if len(mut_annotations['sequence']) >= (aan*3) and mut_annotations['translation'] != ".":
                codon = mut_annotations['sequence'][aan*3-2-1:aan*3] # -1 bcs seq-object 0-based; but positional argument aan not 
                codon = [n for n in codon  ] # turn seq.object to list, seq object does not allow reassignment
                for n in range(4): # test all four nucleotide options for resulting AA change
                    if p_anno.loc['strand'] == 1:
                        codon[ncn-1] = nts[n]
                    else:
                        codon[ncn-1] = rc[n]
                    codonSeqO = Seq( "".join(codon)) 
                    codons_ls.append(codonSeqO)
                    aa_ls.append(codonSeqO.translate())
            mut_annotations['codons'] = codons_ls
            mut_annotations['AA'] = [aa[0] for aa in aa_ls]
            # append_annotations intragenic
            if len(mut_annotations['AA']) < 4:
                mut_annotations['type'] = 'U'
            mut_annotations['AA_gt_custom'] = '' 
            if custom_nt_calls_goodpos[i] != 4: # only if ref defined; [0] ok...see next line comment
                mut_annotations['nt_custom'] = nts[custom_nt_calls_goodpos[i]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                mut_annotations['nts'] = mut_annotations['nt_custom']
                if custom_nt_calls_goodpos[i] != 4:
                    mut_annotations['nt_custom'] = nts[custom_nt_calls_goodpos[i]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
            # extract derived genotype(s) and according AA across all samples
            # call Nonsynonymous if >1 amino acid type exists across all called samples (eg ref AA and alt )
            #print(str(j)+" "+str(callidx))
            #fixed mutation
            # i is index of current pos, j is index of samples 
            if hasmutation_goodpos[i] == True :
                if calls_single_sample_goodpos[i] < 4:
                    # add each NT called, later only keep uniq, but kept order
                    mut_annotations['nts'] = mut_annotations['nts'] + nts[calls_single_sample_goodpos[i]]
                    mut_annotations['nt_alt'] = nts[calls_single_sample_goodpos[i]]
                    if len(mut_annotations['AA']) == 4:
                        mut_annotations['AA_gt'] = mut_annotations['AA_gt'] + mut_annotations['AA'][ calls_single_sample_goodpos[i] ]
                        if np.unique(custom_nt_calls_goodpos[i])[0] != 4:
                            mut_annotations['AA_gt_custom'] = mut_annotations['AA_gt_custom'] + mut_annotations['AA'][ custom_nt_calls_goodpos[i] ]
            if len(mut_annotations['AA']) == 4:
                mut_annotations['type'] = 'S' # eventually overwritten below if N
            # remove duplicates
            mut_annotations['AA_gt_custom'] = ''.join(OrderedDict.fromkeys( mut_annotations['AA_gt_custom'] ).keys()) # get only unique AA and keep order
            mut_annotations['nts'] = ''.join(OrderedDict.fromkeys( mut_annotations['nts'] ).keys()) # get only unique Nuc and keep order
            # record if nonsynonymous mutation
            if len(mut_annotations['AA_gt_custom'])>1:
                mut_annotations['type'] = 'N'
            # Record all observed mutations across all isolates; E.g. K134Y, W47A, etc.
            if len(mut_annotations['AA_gt_custom'])>1:
                mut_annotations['muts'] = []
                ancAA = mut_annotations['AA_gt_custom'][0]
                derAAs = mut_annotations['AA_gt_custom'][1:]
                for j,derAA in enumerate(derAAs):
                    mut_annotations['muts'].append( ancAA+str(mut_annotations['aa_pos'])+derAA )
            else:
                mut_annotations['muts'] = "."
        else: #intergenic
            if int(mut_genenum[i])>0: # get info for gene prior SNP (if any)
                p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])-1] # both -1 necessary bcs list of df and df 0-based
                mut_annotations['gene1'] = p_anno.loc['gene']
                mut_annotations['locustag1'] = p_anno.loc['locustag']
                if 'orthologtag' in p_anno:
                    mut_annotations['orthologtag1'] = p_anno.loc['orthologtag']
                mut_annotations['product1'] = p_anno.loc['product']
                mut_annotations['distance1'] = mut_annotations['pos'] - p_anno.loc['loc2']
                if p_anno.loc['strand'] == -1:
                    mut_annotations['distance1'] = mut_annotations['distance1'] * -1

            if int(mut_genenum[i]+0.5) <= annotation_genes[p_chr-1].shape[0] and annotation_genes[p_chr-1].shape[1] !=0: # get info gene after SNP (if any); second conditional to evade empty chr
                p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])] # -1 necessary bcs list of df 0-based; gene_id 0-based by we want following
                mut_annotations['gene2'] = p_anno.loc['gene']
                mut_annotations['locustag2'] = p_anno.loc['locustag']
                if 'orthologtag' in p_anno:
                    mut_annotations['orthologtag2'] = p_anno.loc['orthologtag']
                mut_annotations['product2'] = p_anno.loc['product']
                mut_annotations['distance2'] = p_anno.loc['loc1'] - mut_annotations['pos'] +1 # +1 to get correct bp distance
                if p_anno.loc['strand'] == 1:
                    mut_annotations['distance2'] = mut_annotations['distance2'] * -1
            # append_annotations intragenic
            #print(mut_annotations['distance1'])
            if ( 'distance1' in mut_annotations and mut_annotations['distance1'] > (-1*promotersize) and mut_annotations['distance1'] < 0) or ( 'distance2' in mut_annotations and mut_annotations['distance2'] > (-1*promotersize) and mut_annotations['distance2'] < 0):
                mut_annotations['type'] = 'P'
            else:
                mut_annotations['type'] = 'I'
            # define nts (repeat of intragenic code)
            if custom_nt_calls_goodpos[i] != 4: # only if ref defined; [0] ok...see next line comment
                mut_annotations['nt_custom'] = nts[custom_nt_calls_goodpos[i]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                mut_annotations['nts'] = mut_annotations['nt_custom']
            # extract derived genotype(s) across all samples
            if hasmutation_goodpos[i] == True:
                if calls_single_sample_goodpos[i] < 4:
                    mut_annotations['nts'] = mut_annotations['nts'] + nts[calls_single_sample_goodpos[i]]
                    mut_annotations['nt_alt'] = nts[calls_single_sample_goodpos[i]]
            # remove duplicates
            mut_annotations['nts'] = ''.join(OrderedDict.fromkeys( mut_annotations['nts'] ).keys()) # get only unique Nuc and keep order
        lod_mutAnno.append(mut_annotations)
    dataframe_mut = pd.DataFrame(lod_mutAnno)
    return dataframe_mut

def write_snp_table(calls_gp_ingroup, sampleNames_ingroup, annotation_mutations_dedup):
        
        ## get NT,sampleName df
        # build sampleNames w/ metainfo
        # snp_table.csv contains 'nt_anc' > nt of explicit single outgroup isolate. '.' == NA
        # sampleNamesLong = apy.annotate_sampleNames(sampleNames,locations_long_names,patients,timepoints,locations) # all sample names translated
        ## apy.annotate_sampleNames commented and replaced by command below 
        
        # translate index to nucleotide
        calls_for_tree = idx2nts(calls_gp_ingroup) # ATCGN translation
        snp_data = pd.DataFrame(calls_for_tree,columns=sampleNames_ingroup)

        ## get snp metadata
        # Contig	Location	Cluster ID	NCTC9343 homolog	Assembly locus	Gene	Protein annotation (Prokka)	Type*	Ancestor allele	Amino acid changes	Nucleotide position in the gene**
        snp_metadata = annotation_mutations_dedup.copy()
        # for I/P mutations: turn all to I (P somewhat random); add up/downstream info to 0tag and gene description to product
        for i,row in snp_metadata.iterrows(): # I/P SNV
            if row.isnull()['locustag']:
                if row.isnull()['locustag1'] and not row.isnull()['locustag2']: # no preceding gene. start of contig
                    snp_metadata.at[i,'locustag'] = "NA;"+str(int(row['distance2']))+":"+row['locustag2']
                    snp_metadata.at[i,'gene'] = "NA;"+row['gene2']
                    snp_metadata.at[i,'product'] = "NA;"+row['product2']
                    snp_metadata.at[i,'type'] = 'I'
                elif row.isnull()['locustag2'] and not row.isnull()['locustag1']: # no subsequent gene. end of contig
                    snp_metadata.at[i,'locustag'] = str(row['distance1'])+":"+str(row['locustag1'])+";NA"
                    snp_metadata.at[i,'gene'] = row['gene1']+";NA"
                    snp_metadata.at[i,'product'] = row['product1']+";NA"
                    snp_metadata.at[i,'type'] = 'I'
                elif row.isnull()['locustag1'] and row.isnull()['locustag2']: # no annotated gene on contig
                    snp_metadata.at[i,'locustag'] = "NA"
                    snp_metadata.at[i,'gene'] = "NA"
                    snp_metadata.at[i,'product'] = "NA"
                    snp_metadata.at[i,'type'] = 'I'
                else: # intergenic SNV with preceding/subsequent gene                
                    snp_metadata.at[i,'locustag'] = str(row['distance1'])+":"+str(row['locustag1'])+";"+str(int(row['distance2']))+":"+row['locustag2']
                    snp_metadata.at[i,'gene'] = row['gene1']+";"+row['gene2']
                    snp_metadata.at[i,'product'] = row['product1']+";"+row['product2']
                    snp_metadata.at[i,'type'] = 'I'
        
        #snp_metadata = snp_metadata[['chr','pos','type','muts','locustag','gene','loc1','loc2','strand','product','nt_pos','aa_pos','nt_ref','nt_alt','nt_anc']]
        #snp_metadata[['nt_anc']] = snp_metadata[['nt_anc']].replace('.','?') # turn NA to '?', similar to NT data 
        snp_metadata = snp_metadata[['chr','pos','type','muts','locustag','gene','loc1','loc2','strand','product','nt_pos','aa_pos','nt_ref', 'nt_anc', 'nt_alt']]
        snp_metadata[['nt_anc']] = snp_metadata[['nt_anc']].replace('.','?') # turn NA to '?', similar to NT data 
        snp_table = pd.concat([snp_metadata,snp_data.reset_index(drop=True)], axis=1)
        
        ## store
        with open('snp_table.csv', 'w') as file:
            snp_table.to_csv(file, header=True, index=False)

def annotate_indels(annotation_genes , indel_p_gp , indel_call_gp, indel_size_gp, hasindel_gp , indel_GL_diff_to_ref , promotersize, ref_genome_folder):
    ''' produce dataframe with annotation info for each indel
    all positions are 1-based! (pos, nt_pos, loc1/2, distance1/2)
    populates with a naive guess at impact of indel, no parsing of other mutation information on the gene (which could impact inference)
    eg 3bp in/del == inframe
    eg 2bp in/del == out of frame, checks is premature stop or no

    NOTE: to be extended using indel impact prediction based on other mutations in this sample (eg if Framshift BUT also compensatory mutation/indel, then OK)
    '''
    # indel_p_gp: genomic position of candidate indels
    # annotation_genes=annotation_genes
    # indel_p_gp=indel_p[cand_indel]
    # indel_call_gp=indel_call[cand_indel, :, :]
    # indel_size_gp=indel_size[cand_indel, :]
    # hasindel_gp=has_indel[cand_indel, :]
    # indel_GL_diff_to_ref=indel_support[cand_indel, :, :]
    # promotersize=promotersize
    # ref_genome_folder=ref_genome_folder

    [chrstarts, genomelength,_]= genomestats(ref_genome_folder);
    gp_pos_to_genenum_l = genomic_position_all_v3(annotation_genes,genomelength,chrstarts); #loc_tag numbers all genes across genome (intergenic:0.5); chr_tag: numbers genes per chr (intergenic: #chr+0.5) // NOTE: cannot subset to indel_p as we need all to identify overlaps of indels vs genes!
    # [locus_tag,chr_tag] = genomic_position_all_v(annotation_genes, genomelength, chrstarts); #loc_tag numbers all genes across genome (intergenic:0.5); chr_tag: numbers genes per chr (intergenic: #chr+0.5)
    # mut_genenum = chr_tag[indel_p_gp]
    # loc_tag_p_gp = locus_tag[indel_p_gp]
    pos_chr = p2chrpos(indel_p_gp,chrstarts) # get two-col table chr,pos   
    lod_indelAnno = []
    for i,indel_pos in enumerate(indel_p_gp):
        indel_annotations = {}
        indel_ref_call = indel_call_gp[i, :, 0]
        indel_ref_call = indel_ref_call[indel_ref_call==indel_ref_call][0]## ref is same across all samples! --> take first (non nan) instance
        indel_ref_end = indel_pos+len(indel_ref_call) ## end of indel on reference 
        alt_calls = indel_call_gp[i, :, 1]
        unique_alt_calls = np.unique(alt_calls[(alt_calls == alt_calls) & (alt_calls != None)]) ## remove nans and get unique list
        gene_num_overlapping = gp_pos_to_genenum_l[indel_pos:indel_ref_end]
        indel_annotations['gene_num'] = np.unique([geneid[1] for pos in gene_num_overlapping for geneid in pos]) ## list per position with list of per gene overlapping with position of [chr_specific_geneid, globals_geneid]
        indel_annotations['gene_num_global'] = np.unique([geneid[0] for pos in gene_num_overlapping for geneid in pos])
        if len(indel_annotations['gene_num']) > 1: ## multiple genetic annotations overlap with site
            indel_annotations['gene_num'] = [geneid for geneid in indel_annotations['gene_num'] if geneid == int(geneid)] ## remove intragenic sites
            indel_annotations['gene_num_global'] = [geneid for geneid in indel_annotations['gene_num_global'] if geneid == int(geneid)] ## remove intragenic sites
        if len(indel_annotations['gene_num']) > 1: ## multiple genes overlap with site
            print('WARNING MUTLIPLE GENES ARE OVERLAPPING WITH INDEL. HANDLING THIS IS NOT YET IMPLEMENTED!')
        indel_annotations['chr'] = pos_chr[i][0]
        indel_annotations['pos'] = pos_chr[i][1] + 1 # turn 0-based pos into 1-based pos.
        indel_annotations['chr_boundary'] = (indel_pos in chrstarts) | (indel_ref_end in chrstarts) ## check if indel overlaps with chr boundaries --> poor alignments
        indel_annotations['indel_ref'] = indel_ref_call
        indel_annotations['indel_anc'] = "." # default. should be called during CMT generation!
        indel_annotations['indel_alt'] = unique_alt_calls
        # min value in indel_support which is in has indel
        indel_annotations['max_indel_GL_diff_to_ref'] = min(indel_GL_diff_to_ref[i,:])
        p_chr = pos_chr[i][0]
        indel_annotations['indel_size_gp']=np.unique(indel_size_gp[i,hasindel_gp[i,:]])
        # indel_size_gp_gt=[]
        # for j,indel_callidx in enumerate(indel_size_gp[i,:]):
        #     # i is index of current pos, j is index of samples 
        #     if hasindel_gp[i,j] == True :
        #         indel_size_gp_gt.append(indel_size_gp[i,j])
        # indel_annotations['indel_size_gp'] = list(dict.fromkeys( indel_size_gp_gt ))
        for indel_genenum in indel_annotations['gene_num']:
            if indel_genenum == int(indel_genenum): # intragenic
                p_anno = annotation_genes[p_chr-1].iloc[int(indel_genenum)-1] # both -1 necessary bcs list of df and df 0-based
                indel_annotations['type'] = 'G'
                indel_annotations['product'] = p_anno.loc['product']
                indel_annotations['gene'] = p_anno.loc['gene']
                indel_annotations['protein_id'] = p_anno.loc['protein_id']
                indel_annotations['strand'] = p_anno.loc['strand']
                indel_annotations['loc1'] = p_anno.loc['loc1'] + 1 # 1-based first position
                indel_annotations['loc2'] = p_anno.loc['loc2'] # last position of gene (inclusive)
                indel_annotations['sequence'] = p_anno.loc['sequence']
                indel_annotations['protein_id'] = p_anno.loc['protein_id']
                indel_annotations['note'] = p_anno.loc['note']
                indel_annotations['locustag'] = p_anno.loc['locustag']
                if 'orthologtag' in p_anno:
                    indel_annotations['orthologtag'] = p_anno.loc['orthologtag']
                indel_annotations['translation'] = p_anno.loc['translation']
                gene_length = (indel_annotations['loc2']-indel_annotations['loc1'])
                for key in ['indel_pos_start', 'indel_pos_end', 'alt_start', 'alt_stop', 'alt_stop_pos', 'frameshift', 'mut_translation', 'mut', 'num_isolates']:
                    indel_annotations[key] = []
                # impact of indel
                for unique_call_idx, unique_call in enumerate(unique_alt_calls):
                    indel_annotations['num_isolates'].append(np.sum([call == unique_call for call in indel_call_gp[i, hasindel_gp[i, :], 1]]))
                    if p_anno.loc['strand'] == 1: # get position within gene, consider strandedness
                        indel_annotations['indel_pos_start'].append(indel_annotations['pos'] - indel_annotations['loc1']+1) # pos/loc1 1-based. nt_pos 1-based. +1 to get 1-based nt_pos in gene (checked)!
                        indel_annotations['indel_pos_end'].append(indel_annotations['indel_pos_start'][unique_call_idx]+len(indel_ref_call))
                        unique_call = Seq(unique_call)
                    elif p_anno.loc['strand'] == -1:
                        indel_annotations['indel_pos_end'].append(indel_annotations['loc2'] - indel_annotations['pos'] +1) # pos/loc2 1-based. +1 to get 1-based nt_pos in gene (checked)!
                        indel_annotations['indel_pos_start'].append(indel_annotations['indel_pos_end'][unique_call_idx]-len(indel_ref_call)) ## subtraction as we are on negative strand!
                        unique_call = Seq(unique_call).reverse_complement() ## note genes are already rev complement for antisense strands!
                    else:
                        indel_annotations['indel_pos_end'].append('.')
                        indel_annotations['indel_pos_start'].append('.')
                        unique_call = Seq(unique_call)
                    start_of_indel = indel_annotations['indel_pos_start'][unique_call_idx]
                    end_of_indel = indel_annotations['indel_pos_end'][unique_call_idx]
                    ## check for gene boundaries
                    if start_of_indel < 0: ## indel overlaps with gene start
                        ## identify if the start codon is still present or a start codon is present at any other site
                        alt_seq = unique_call+indel_annotations['sequence'][end_of_indel:]
                        start_codons = [alt_seq.find(string) for string in ('ATG', 'GTG', 'TTG')]
                        if start_codons != []:
                            alt_seq = alt_seq[start_codons[0]:]
                            indel_annotations['alt_start'].append(start_codons)
                        else:
                            indel_annotations['alt_start'].append('-')
                    elif (start_of_indel > 0) & (end_of_indel < gene_length): ## indel is in gene
                        alt_seq = indel_annotations['sequence'][:start_of_indel] + unique_call + indel_annotations['sequence'][end_of_indel:]
                        indel_annotations['alt_start'].append('-')
                    elif end_of_indel > gene_length: ## indel overlaps with gene end
                        alt_seq = indel_annotations['sequence'][:start_of_indel] + unique_call
                    else: 
                        print(f'WARNING: Something went wrong in the annotation of the indels for indel {indel_annotations}')
                        
                    frameshift = len(alt_seq) % 3
                    indel_annotations['frameshift'].append(frameshift)
                    if frameshift > 0:
                        alt_seq_pad = alt_seq[:-frameshift] ## remove nt which do not generate entire codon due to frameshift
                        indel_annotations['mut_translation'].append(alt_seq_pad.translate())
                    else:
                        indel_annotations['mut_translation'].append(alt_seq.translate())
                    indel_annotations['mut'].append(indel_annotations['translation'] != indel_annotations['mut_translation'][unique_call_idx]) ## check if indel caused actual alteration of amino acid code
                    alt_stop_codon_num = indel_annotations['mut_translation'][unique_call_idx].find("*")+1 ## check if ALT creates an alternative stop
                    indel_annotations['alt_stop'].append(str(np.where((alt_stop_codon_num == -1), 'extended+', ## protein is extended beyond previous stop
                                                                       np.where(alt_stop_codon_num > len(indel_annotations['translation']), 'extended',  ## insertion extended protein but stop is still in sequence of previous gene
                                                                                np.where(alt_stop_codon_num == len(indel_annotations['translation']), 'same', ## length of gene stayed the same
                                                                                        'truncated'))))
                                                        )
                    indel_annotations[f'alt_stop_pos'].append( len(indel_annotations['mut_translation'][unique_call_idx].split("*", 1)[0]) +1)
            else: #intergenic
                indel_annotations['num_isolates'] = []
                for unique_call in unique_alt_calls:
                    indel_annotations['num_isolates'].append([np.sum([call == unique_call for call in indel_call_gp[i, hasindel_gp[i, :], 1]])])
                if int(indel_genenum)>0: # get info for gene prior SNP (if any)
                    p_anno = annotation_genes[p_chr-1].iloc[int(indel_genenum)-1] # both -1 necessary bcs list of df and df 0-based
                    indel_annotations['gene1'] = p_anno.loc['gene']
                    indel_annotations['locustag1'] = p_anno.loc['locustag']
                    if 'orthologtag' in p_anno:
                        indel_annotations['orthologtag1'] = p_anno.loc['orthologtag']
                    indel_annotations['product1'] = p_anno.loc['product']
                    indel_annotations['distance1'] = indel_annotations['pos'] - p_anno.loc['loc2']
                    if p_anno.loc['strand'] == -1:
                        indel_annotations['distance1'] = indel_annotations['distance1'] * -1
                if int(indel_genenum+0.5) <= annotation_genes[p_chr-1].shape[0] and annotation_genes[p_chr-1].shape[1] !=0: # get info gene after Indel (if any); second conditional to evade empty chr
                    p_anno = annotation_genes[p_chr-1].iloc[int(indel_genenum)] # -1 necessary bcs list of df 0-based; gene_id 0-based by we want following
                    indel_annotations['gene2'] = p_anno.loc['gene']
                    indel_annotations['locustag2'] = p_anno.loc['locustag']
                    if 'orthologtag' in p_anno:
                        indel_annotations['orthologtag2'] = p_anno.loc['orthologtag']
                    indel_annotations['product2'] = p_anno.loc['product']
                    indel_annotations['distance2'] = p_anno.loc['loc1'] - indel_annotations['pos'] +1 # +1 to get correct bp distance
                    if p_anno.loc['strand'] == 1:
                        indel_annotations['distance2'] = indel_annotations['distance2'] * -1
                if ( 'distance1' in indel_annotations and indel_annotations['distance1'] > (-1*promotersize) and indel_annotations['distance1'] < 0) or ( 'distance2' in indel_annotations and indel_annotations['distance2'] > (-1*promotersize) and indel_annotations['distance2'] < 0):
                    indel_annotations['type'] = 'P'
                else:
                    indel_annotations['type'] = 'I'
        lod_indelAnno.append(indel_annotations)
    dataframe_indel = pd.DataFrame(lod_indelAnno)
    return dataframe_indel

def write_indel_table(annotation_indels, hasindel_gp, indel_call_gp, sampleNames):
    call_mat = np.where(hasindel_gp, indel_call_gp[:, :, 1], indel_call_gp[:, :, 0])
    call_df = pd.DataFrame(data = call_mat, columns = sampleNames)
    indel_table = pd.merge(annotation_indels, call_df, how = 'outer', left_index=True, right_index=True)
    indel_table.to_csv('indel_table.csv', header=True, index=False)
    

def annotate_sampleNames(samplenames,locations_long_names,patients_sbj,visits_sbj,locations_sbj):
    # extend sample name with patient/visit/location identifier. all in same order!
    extendend_sampleNames = np.copy(samplenames)
    for i,name in enumerate(extendend_sampleNames):
        extendend_sampleNames[i] = "S"+patients_sbj[i]+"_V"+str(visits_sbj[i])+"_"+locations_long_names[locations_sbj[i]]+"_"+name
    return extendend_sampleNames
       
def write_calls_sampleName_to_fasta(calls_for_tree,treeSampleNames,timestamp):
    fa_file = open(timestamp+".fa", "w")
    for i,name in enumerate(treeSampleNames):
        nucl_string = "".join(list(calls_for_tree[:,i]))
        fa_file.write(">" + name + "\n" + nucl_string + "\n")
    fa_file.close()

def parse_tree(tree_path):
    """helper function to parse tree file, with error catching. Returns first tree in a file, if multiple"""
    try:
        parsed_tree=Phylo.read(tree_path,'newick')
    except ValueError:
        parsed_tree=None
        print("Multiple trees in file, trying to parse first (filled) tree in file.")
        trees= Phylo.parse(tree_path,'newick')
        for tree in trees:
            if tree.count_terminals() > 1:
                parsed_tree=tree
                print("Found filled tree, with length", tree.count_terminals())
                break
    return parsed_tree

def find_clade_terminals(clade_name, tree):
    """Finds a sample name in a tree and returns all sample names in a list at or below this branch in the clade"""
    parsed_tree = parse_tree(tree)
    clades_found=[]
    for clade in parsed_tree.find_elements(clade_name):
        clades_found.append(clade)
    if len(clades_found) == 0:
        print("Error, no terminal branch found with input clade_name:",clade_name)
    else:
        parent=parsed_tree.get_path(clades_found[0])[-2] ## find parent node
        ancestor_clades=parent.get_terminals()
        output=[]
        for ancestor in ancestor_clades:
            output.append(ancestor.name)
        return output


def count_number_mutational_events(ancestral_reconstruction_labelled_tree,ancestral_reconstruction_fasta, skip_preterminals=True, track_nodes=False):
    """
    Counts the number of mutational events at all positions.
    Uses a treetime output ancestral reconstruction tree (with labelled internal nodes) and calls for all the goodpos that went into building the tree and ancestral reconstruction. 
    These outputs are from treetime ancestral, see function create_ancestral_reconstruction.
    
    Homoplasic locations are all cells with value >1

    Options to count or skip terminal node transitions, (eg only count internal nodes).

    Default output names from treetime:
    annotated_tree.nexus == ancestral_reconstruction_labelled_tree
    ancestral_sequences.fasta == ancestral_reconstruction_fasta

    Example call:
    count_number_mutational_events(f'{analysis_params_output_name}_ancestral_reconstruction/annotated_tree.nexus',f'{analysis_params_output_name}_ancestral_reconstruction/ancestral_sequences.fasta'"""
    # parse fasta
    calls_pos={}
    for record in SeqIO.parse(ancestral_reconstruction_fasta, 'fasta'):
        calls_pos[record.id] = str(record.seq)
    # parse tree structure into adjacency graph to traverse
    parsed_tree = parse_tree(ancestral_reconstruction_labelled_tree) 
    root=parsed_tree.common_ancestor(parsed_tree.get_terminals())
    net = Phylo.to_networkx(parsed_tree)
    tree_as_dict_of_lists=networkx.to_dict_of_lists(net)
    # start tree traversal, checking if each internal node has same call as parent, if not iterate the val for that base +1
    transitions_goodpos=np.zeros(len(str(record.seq)))
    to_visit=[(x, root) for x in tree_as_dict_of_lists[root]]
    visited=set()
    visited.add(root)
    nodes_output={i:[] for i in range(len(transitions_goodpos))}
    while len(to_visit)>0:
        currently_processing=to_visit[0]
        parent=currently_processing[1]
        current_node=currently_processing[0]
        ## check if parent is really parent of node --> otherwise at root it will be 
        if not parsed_tree.find_any(parent.name).is_parent_of(current_node.name):
            to_visit=to_visit[1:]
            continue
        visited.add(current_node)
        is_preterminal=len(tree_as_dict_of_lists[current_node])==1
        if skip_preterminals and is_preterminal:
            pass
        else:
            for index in range(len(transitions_goodpos)):
                if calls_pos[parent.name][index] != calls_pos[current_node.name][index] and calls_pos[current_node.name][index]!= 'N' and calls_pos[parent.name][index]!= 'N':
                    transitions_goodpos[index]+=1
                    if track_nodes:
                        nodes_output[index].append((current_node.name,calls_pos[current_node.name][index]))
            for children in tree_as_dict_of_lists[current_node]:
                if children not in visited:
                    to_visit.append((children,current_node))
        to_visit=to_visit[1:]
    if track_nodes:
        return transitions_goodpos,nodes_output
    return transitions_goodpos

def check_multiply_mutated_pos_bootraps(bootstrap_tree, ancestral_reconstruction_labelled_tree ,multiply_mutated_pos_nodes_to_basecalls,minimum_bootstrap_value_at_MRCA,return_nodes_to_compare=False,return_all_bs_values=False):
    """Confirms that homoplasies have MRCA node with bootstrap above given threshold"""
    # for each pos which is multiply mutated,
    #   find all clades which got mutated (number of mutational events)
    #   find all pairwise MRCAs, for each clade, use MRCA comparison which is smallest amount of descendants (eg compare closeby preferrentially ) 
    #   for each pariwise smallest MRCA, check MRCA-node bootstrap, if >thresh then OK, we believe the independent mutation event, if not then it is not independent mutational event, subtract 1 from the counts of multiply mutated pos
    bs_tree=parse_tree(bootstrap_tree)
    named_internal_node_tree=parse_tree(ancestral_reconstruction_labelled_tree)
    max=named_internal_node_tree.count_terminals()+len(named_internal_node_tree.get_nonterminals())
    collector_of_node_to_compare=[]
    for pos_nodes_calls in multiply_mutated_pos_nodes_to_basecalls.values():
        node_to_compare={}
        if len(pos_nodes_calls) > 1:
            node_combos_to_check={}
            for index_in_node_call_list,node_call in enumerate(pos_nodes_calls):
                shortest=max
                node_combos_to_check[node_call] = [pos_nodes_calls[x] for x in range(len(pos_nodes_calls)) if x != index_in_node_call_list]
                if len(node_combos_to_check[node_call])>0:
                    for node,value in node_combos_to_check[node_call]:
                        if len(named_internal_node_tree.trace(node_call[0],node)) < shortest:
                            node_to_compare[node_call[0]]=node
                            shortest=len(named_internal_node_tree.trace(node_call[0],node))
        collector_of_node_to_compare.append(node_to_compare)
    if return_nodes_to_compare: 
        return collector_of_node_to_compare
    elif return_all_bs_values:
        return return_all_bs_values_homoplasy(bs_tree,named_internal_node_tree,collector_of_node_to_compare)
    else: 
        return helper_for_homoplasy_boostrap_check(bs_tree,named_internal_node_tree,minimum_bootstrap_value_at_MRCA,collector_of_node_to_compare)

def helper_for_homoplasy_boostrap_check(parsed_bs_tree,parsed_anc_internal_node_tree,minimum_bootstrap_value_at_MRCA,nodes_to_check):
    correction_for_nodes_to_check=[]
    for node_to_check_dict in nodes_to_check:
        correction_for_nodes_to_check.append(0)
        if len(node_to_check_dict)>1:
            for key,value in node_to_check_dict.items():
                key_descendant = list(parsed_anc_internal_node_tree.find_clades(key))[0].get_terminals()[0].name
                value_descendant = list(parsed_anc_internal_node_tree.find_clades(value))[0].get_terminals()[0].name
                value_of_MRCA = parsed_bs_tree.common_ancestor(list(parsed_bs_tree.find_clades(value_descendant))+list(parsed_bs_tree.find_clades(key_descendant))).confidence
                descendants_of_parent=[list(parsed_bs_tree.find_clades(x.name))[0] for x in parsed_anc_internal_node_tree.get_path(list(parsed_anc_internal_node_tree.find_clades(key))[0])[-2].get_terminals() if len(list(parsed_bs_tree.find_clades(x.name)))>0]
                value_of_parent_node = parsed_bs_tree.common_ancestor(descendants_of_parent).confidence
                if value_of_parent_node < minimum_bootstrap_value_at_MRCA:
                    correction_for_nodes_to_check[-1]-=1
    return correction_for_nodes_to_check

def return_all_bs_values_homoplasy(parsed_bs_tree,parsed_anc_internal_node_tree,nodes_to_check):
    bs_values_closest_intermediates=[]
    print(len(nodes_to_check))
    for node_to_check_dict in nodes_to_check:
        if len(node_to_check_dict)>1:
            bs_collector_for_this_pos=[]
            for key,value in node_to_check_dict.items():
                key_descendants = [x.name for x in list(parsed_anc_internal_node_tree.find_clades(key))[0].get_terminals() if x.name != 'Sref']
                value_descendants = [x.name for x in list(parsed_anc_internal_node_tree.find_clades(value))[0].get_terminals() if x.name != 'Sref']
                key_node_in_bs_tree = parsed_bs_tree.common_ancestor(key_descendants)
                value_node_in_bs_tree = parsed_bs_tree.common_ancestor(value_descendants)
                nodes_between_appearance_of_homoplasy = parsed_bs_tree.trace(key_node_in_bs_tree,value_node_in_bs_tree)
                bs_values_between_appearance_of_homoplasy = [x.confidence for x in nodes_between_appearance_of_homoplasy]
                bs_collector_for_this_pos.append(bs_values_between_appearance_of_homoplasy)
            bs_values_closest_intermediates.append(bs_collector_for_this_pos)
        else:
            bs_values_closest_intermediates.append([])
    return bs_values_closest_intermediates

def count_homoplasy_snp_differences(ancestral_reconstruction_labelled_tree,ancestral_reconstruction_fasta,multiply_mutated_pos_nodes_to_basecalls, output_all_pairwise=True, output_one_each=False):
    """
        Function outputs number of SNP differences between nodes where homoplasic mutation reappears
        This functions hould plot a histogram of branch lengths separating homoplasmic internal nodes
        Hopefully this helps ID a reasonable cutoff for removing homoplasic loci (if biomodal)
    """
    named_internal_node_tree=parse_tree(ancestral_reconstruction_labelled_tree)
    calls_pos={}
    for record in SeqIO.parse(ancestral_reconstruction_fasta, 'fasta'):
        calls_pos[record.id] = str(record.seq)
    collector_of_node_to_compare=[]
    max=len(calls_pos)
    if not output_one_each:
        for idx,pos_nodes_calls in multiply_mutated_pos_nodes_to_basecalls.items():
            node_to_compare={}
            if len(pos_nodes_calls) > 1:
                node_combos_to_check={}
                for index_in_node_call_list,node_call in enumerate(pos_nodes_calls):
                    shortest=max
                    node_combos_to_check[node_call] = [pos_nodes_calls[x] for x in range(len(pos_nodes_calls)) if x != index_in_node_call_list and pos_nodes_calls[x][1] == node_call[1]]
                    if len(node_combos_to_check[node_call])>0:
                        for node,value in node_combos_to_check[node_call]:
                            if len(named_internal_node_tree.trace(node_call[0],node)) < shortest:
                                node_to_compare[idx]=[node_call[0],node]
                                if not output_all_pairwise:
                                    shortest=len(named_internal_node_tree.trace(node_call[0],node))
                                else:
                                    collector_of_node_to_compare.append(np.sum(np.array([x for x in calls_pos[node_call[0]]]) != np.array([x for x in calls_pos[node]])))
            if not output_all_pairwise:
                for node_list in node_to_compare.values():
                    collector_of_node_to_compare.append(np.sum(np.array([x for x in calls_pos[node_list[0]]]) != np.array([x for x in calls_pos[node_list[1]]])))
        return collector_of_node_to_compare
    else:
        node_to_compare={}
        for idx,pos_nodes_calls in multiply_mutated_pos_nodes_to_basecalls.items():
            if len(pos_nodes_calls) > 1:
                node_combos_to_check={}
                for index_in_node_call_list,node_call in enumerate(pos_nodes_calls):
                    shortest=max
                    node_combos_to_check[node_call] = [pos_nodes_calls[x] for x in range(len(pos_nodes_calls)) if x != index_in_node_call_list]
                    if len(node_combos_to_check[node_call])>0:
                        for node,value in node_combos_to_check[node_call]:
                            if len(named_internal_node_tree.trace(node_call[0],node)) < shortest:
                                node_to_compare[idx]=[node_call[0],node]
                                shortest=len(named_internal_node_tree.trace(node_call[0],node))
        for node_list in node_to_compare.values():
            collector_of_node_to_compare.append(np.sum(np.array([x for x in calls_pos[node_list[0]]]) != np.array([x for x in calls_pos[node_list[1]]])))
        return collector_of_node_to_compare

def count_all_pairwise_snp_differences(ancestral_reconstruction_fasta,include_terminal=True):
    """
        Function outputs number of SNP differences between all pairwise node combinations
        This functions hould plot a histogram of branch lengths separating homoplasmic internal nodes
        Hopefully this helps ID a reasonable cutoff for removing homoplasic loci (if biomodal)
    """
    calls_pos={}
    for record in SeqIO.parse(ancestral_reconstruction_fasta, 'fasta'):
        if not include_terminal:
            if 'NODE' not in record.id:
                calls_pos[record.id] = str(record.seq)
        else:
            calls_pos[record.id] = str(record.seq)
    collector_of_distances=[]
    counted_nodes=set()
    for self_node in calls_pos:
        for comparison_node in calls_pos:
            if self_node != comparison_node and comparison_node not in counted_nodes:
                collector_of_distances.append(np.sum(np.array([x for x in calls_pos[self_node]]) != np.array([x for x in calls_pos[comparison_node]])))
        counted_nodes.add(self_node)
    return collector_of_distances


def count_number_mutational_events_basal(ancestral_reconstruction_labelled_tree,calls_nodes_to_goodpos_fasta, skip_preterminals=True, node_to_start_basal_to=''):
    """Counts the number of mutational events using an ancestral reconstruction tree (with labelled internal nodes) and calls for all the goodpos that went into building the tree and ancestral reconstruction. 
    These outputs are from treetime ancestral, see function create_ancestral_reconstruction.
    
    Default output names from treetime:
    annotated_tree.nexus == ancestral_reconstruction_labelled_tree
    ancestral_sequences.fasta == calls_nodes_to_goodpos_fasta"""
    ## NOTE: not sure if it 100% works right now.
    # parse fasta
    calls_pos={}
    for record in SeqIO.parse(calls_nodes_to_goodpos_fasta, 'fasta'):
        calls_pos[record.id] = str(record.seq)
    # parse tree structure into adjacency graph to traverse
    parsed_tree = parse_tree(ancestral_reconstruction_labelled_tree) 
    if node_to_start_basal_to != '':
        node_to_start_basal_to_output=[]
        for clade in parsed_tree.find_elements(node_to_start_basal_to):
            node_to_start_basal_to_output.append(clade)
        if len(node_to_start_basal_to_output) == 0:
            print("Error, no terminal branch found with input clade_name:",node_to_start_basal_to)
        else:
            root=parsed_tree.get_path(node_to_start_basal_to_output[0])[-2] ## find parent node
    else:
        root=parsed_tree.common_ancestor(parsed_tree.get_terminals())
    downstream_from_root=root.get_terminals()+root.get_nonterminals()
    net = Phylo.to_networkx(parsed_tree)
    tree_as_dict_of_lists=networkx.to_dict_of_lists(net)
    # start tree traversal, checking if each internal node has same call as parent, if not iterate the val for that base +1
    transitions_goodpos=np.zeros(len(str(record.seq)))
    to_visit=[(x, root) for x in tree_as_dict_of_lists[root] if x in downstream_from_root and x not in [root]]
    visited=set()
    visited.add(root)
    while len(to_visit)>0:
        currently_processing=to_visit[0]
        parent=currently_processing[1]
        current_node=currently_processing[0]
        visited.add(current_node)
        is_terminal=len(tree_as_dict_of_lists[current_node])==1
        if is_terminal and not skip_preterminals:
            for index in range(len(transitions_goodpos)):
                if calls_pos[parent.name][index] != calls_pos[current_node.name][index] and calls_pos[current_node.name][index]!= 'N' and calls_pos[parent.name][index]!= 'N':
                    transitions_goodpos[index]+=1
        elif not is_terminal:
            for index in range(len(transitions_goodpos)):
                if calls_pos[parent.name][index] != calls_pos[current_node.name][index] and calls_pos[current_node.name][index]!= 'N' and calls_pos[parent.name][index]!= 'N':
                    transitions_goodpos[index]+=1
        for children in tree_as_dict_of_lists[current_node]:
            if children not in visited:
                to_visit.append((children,current_node))
        to_visit=to_visit[1:]
    return transitions_goodpos

def depth_first_search_order(tree,node_to_start_basal_to):
    """Function returns breadth first search order of all samples downstream of provided node_to_star_basal_to"""
    parsed_tree = parse_tree(tree) 
    node_to_start_basal_to_output=[]
    for clade in parsed_tree.find_elements(node_to_start_basal_to):
        node_to_start_basal_to_output.append(clade)
    if len(node_to_start_basal_to_output) == 0:
        print("Error, no terminal branch found with input clade_name:",node_to_start_basal_to)
    else:
        root=parsed_tree.get_path(node_to_start_basal_to_output[0])[-2] ## find parent node
    downstream_from_root=root.get_terminals()+root.get_nonterminals()
    net = Phylo.to_networkx(parsed_tree)
    tree_as_dict_of_lists=networkx.to_dict_of_lists(net)
    # start tree traversal, checking if each internal node has same call as parent, if not iterate the val for that base +1
    to_visit=[x for x in tree_as_dict_of_lists[root] if x in downstream_from_root and x not in [root]]
    order=[]
    visited=set()
    while len(to_visit)>0:
        visiting=to_visit[0]
        visited.add(visiting)
        if visiting.is_terminal():
            order.append(visiting.name)
        for children in tree_as_dict_of_lists[visiting]:
            if children not in visited and children in downstream_from_root:
                to_visit.append(children)
        to_visit=to_visit[1:]
    return order

def get_internal_node_nt_calls(ancestral_reconstruction_tree,ancestral_reconstruction_fasta,treminals_to_find_MRCA_of=[]):
    """
    Returns the (goodpos) sequence of the internal nodes which is an MRCA of the input list of sample names
    Use case: generating analysis sets for dn/ds, muts, etc. or finding sequence of MRCA for comparissons.
    """
    # parse fasta
    calls_pos={}
    for record in SeqIO.parse(ancestral_reconstruction_fasta, 'fasta'):
        calls_pos[record.id] = str(record.seq)
    # parse tree structure into adjacency graph to traverse
    parsed_tree = parse_tree(ancestral_reconstruction_tree)
    found_nodes=[x for x in parsed_tree.get_terminals() if x.name in treminals_to_find_MRCA_of]    
    internal_node=parsed_tree.common_ancestor(found_nodes)
    return calls_pos[internal_node.name]


def generate_tree(calls_for_tree,treeSampleNamesLong,sampleNamesDnapars,filetag='',buildTree=False,writeDnaparsAlignment=False,root_on_smpl='',coloring_pattern='', raxml_model = 'GTR+G', additional_raxml_parameters='--threads 2 --seed 123', bootstrapping = False):
    # Write alignment file (as fasta)
    # calc NJ or Parsimonous tree or None
    # writeDnaparsAlignment==True for writing dnapars input for usage on cluster
    # can use --ancestral with additional raxml params to reconstruct internal nodes
    ts = time.time()
    timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
        
    # write alignment fasta,read alignment
    write_calls_sampleName_to_fasta(calls_for_tree,treeSampleNamesLong,timestamp) #timestamp.fa
    if writeDnaparsAlignment:
        # change tip labels and write phylip
        write_calls_sampleName_to_fasta(calls_for_tree,sampleNamesDnapars,timestamp+"_"+filetag+"_dnapars") #timestamp_dnapars.fa > for dnapars...deleted later
        # turn fa to phylip and delete fasta with short tip labels    
        aln = AlignIO.read(timestamp+"_"+filetag+"_dnapars.fa", 'fasta')
        AlignIO.write(aln, timestamp+"_"+filetag+".phylip", "phylip")
        subprocess.run(["rm -f " + timestamp+"_"+filetag+"_dnapars.fa"],shell=True)
        print("Written file: " + timestamp+"_"+filetag+".phylip")
        # write parameter file
        with open(timestamp+"_"+filetag+"_options.txt",'w') as file:
            file.write(timestamp+"_"+filetag+".phylip"+"\n")
            file.write("f"+"\n")
            file.write(timestamp+"_"+filetag+"_out.dnapars"+"\n")
            file.write("5"+"\n")
            file.write("V"+"\n")
            file.write("1"+"\n")
            file.write("y"+"\n")
            file.write("f"+"\n")
            file.write(timestamp+"_"+filetag+".tree"+"\n"+"\n")



    if buildTree=='PS':      
        # write phylip file with dnaparse compatible 10c samplenames
        write_calls_sampleName_to_fasta(calls_for_tree,sampleNamesDnapars,timestamp+"_dnapars") #timestamp_dnapars.fa > for dnapars...deleted later
        # turn fa to phylip and delete fasta with short tip labels    
        aln = AlignIO.read(timestamp+"_dnapars.fa", 'fasta')
        AlignIO.write(aln, timestamp+".phylip", "phylip")
        subprocess.run(["rm -f " + timestamp+"_dnapars.fa"],shell=True)
            
        # find dnapars executable
        dnapars_path = glob.glob('dnapars')
        path_extension = "../"
        backstop = 0
        while len(dnapars_path) == 0 and backstop <= 5:
            dnapars_path = glob.glob(path_extension+'dnapars')
            path_extension = path_extension + "../"
            backstop = backstop + 1
        if len(dnapars_path) == 0:
            raise ValueError('dnapars executable could not be located.')
        elif dnapars_path[0]=='dnapars':
            dnapars_path[0] = './dnapars'
        # write parameter file
        with open(timestamp+"_options.txt",'w') as file:
            file.write(timestamp+".phylip"+"\n")
            file.write("f"+"\n")
            file.write(timestamp+"_out.dnapars"+"\n")
            file.write("5"+"\n")
            file.write("V"+"\n")
            file.write("1"+"\n")
            file.write("y"+"\n")
            file.write("f"+"\n")
            file.write(timestamp+".tree"+"\n"+"\n")

        # run dnapars
        print("Build parsimony tree...")
        #print( dnapars_path[0] + " < " + timestamp+"_options.txt > " + timestamp+"_dnapars.log")
        subprocess.run([ "touch outtree"  ],shell=True)
        subprocess.run([ dnapars_path[0] + " < " + timestamp+"_options.txt > " + timestamp+"_dnapars.log"  ],shell=True)
        # print('done')
        # re-write tree with new long tip labels        
        tree = Phylo.read(timestamp+".tree", "newick")
        for leaf in tree.get_terminals():
            # print(leaf.name)
            idx = np.where(sampleNamesDnapars==leaf.name) 
            if len(idx[0]) > 1:
                warnings.warn("Warning: dnapars 10c limit leads to ambigous re-naming for "+leaf.name)
                idx = idx[0][0] #np.where returns: tuple with array with index
            else:
                idx = idx[0][0] #np.where returns: tuple with array with index
            leaf.name = treeSampleNamesLong[idx]
        
        Phylo.write(tree, timestamp+".tree", 'nexus')

        if root_on_smpl != '': 
            tree_path = reroot_tree(timestamp+".tree", root_on_smpl, tree_type = 'nexus')
            print('done')
        else:
            tree_path = timestamp+".tree"

        ## section to color treebranches in nexus file
        if coloring_pattern:
            ## get unique patterns
            patterns = []
            for i, sampleName in enumerate(treeSampleNamesLong):
                if i != 0: ## skip first name entry which is reference
                    sampleName_substr = sampleName.split('_')
                    pattern = [pattern for pattern in sampleName_substr if coloring_pattern in pattern]
                    if not pattern or len(pattern) > 1: ## if pattern was not found or multiple patterns have been found
                        print('Multiple entries of coloring pattern were found! Please specify something specific for coloring or do not include coloring in function!')
                        ## stop function to do not mess things up
                        return 
                    else:
                        patterns.append(pattern[0]) ## index is important as a list was returned via list comprehension!
                    
            # patterns = [int(sampleName.split('_')[1][1:]) for i, sampleName in enumerate(treeSampleNamesLong) if i != 0] ## patterns are in second entry, and remove string (first character); mind, that reference genome has no pattern entry! --> Skip
            unique_patterns = np.sort(np.unique(patterns))

            minval = 0.25; maxval = 0.9; num_colors = len(unique_patterns) ## set borders in which colors should be taken from color map 
            cmap = plt.get_cmap('Greens', num_colors) 
            cmap_subset = mcolors.LinearSegmentedColormap.from_list('Custom cmap', cmap(np.linspace(minval, maxval, num_colors)), cmap.N)  ## subset cmap to remove unreadable bright colors 
            color = {}
            for color_entry, pattern in zip(range(cmap_subset.N), unique_patterns):
                rgba = cmap(color_entry)
                color[pattern] = (int(rgba[0]*255) << 16) + (int(rgba[1]*255) << 8) + int(rgba[2]*255) ## set key back to pattern value (python is 0-based!) including leading character; get decimal color code by bit shifting 
            ## read tree again 
            f = open(tree_path, 'r').readlines()
            ## open outfile
            fo = open(tree_path.replace('_latest', '_latest_colored'), 'w')
            ## get beginning of sample names 
            start_manipulating = False
            for line in f: 
                if start_manipulating:
                    if ';' not in line:
                        print('Warning! Taxlabels are on separate lines and coloring of samples will not be successful!')
                    taxlabels = line.strip().split(' ')
                    new_taxlabels_all = ''
                    for taxlabel in taxlabels: 
                        if ('taxlabels' in taxlabel.lower()) or (taxlabel == ''): 
                            new_taxlabels_all = '\t' + taxlabel + '\n' ## get taxlabel string (first entry of line)
                        else:
                            try: ## try splitting
                                curr_pattern = [pattern for pattern in taxlabel.split('_') if coloring_pattern in pattern][0] ## no checks for multiple pattern necessary as it was checked above, extract first (and only) one
                                new_taxlabel = '\'' + taxlabel.strip(';') + '\'[&!color=#-' + str(color[curr_pattern]) + ']'
                            except: ## no pattern available --> black labelling
                                new_taxlabel = '\'' + taxlabel.strip(';') + '\'[&!color=#-16777216]'
                            ## append all to one line again
                            new_taxlabels_all = new_taxlabels_all + '\t' + new_taxlabel  + '\n'
                            if taxlabel[-1] == ';': ## last string in taxlabels
                                start_manipulating = False
                                new_taxlabels_all = new_taxlabels_all + ';' ## add ending, too
                    fo.write(new_taxlabels_all + '\n') ## write all taxlabels to outfile 
                else: ## write all other lines to outfile 
                    fo.write(line)
                if 'dimensions ntax=' in line.lower():
                    start_manipulating = True
            fo.close()
        Phylo.write(tree, filetag+"_latest.nwk.tree", 'newick')
    
        if root_on_smpl != '': 
            reroot_tree(filetag+"_latest.nwk.tree", root_on_smpl)

        # clean up
        subprocess.run(["rm -f " + timestamp+".phylip " + timestamp+"_options.txt " + timestamp+"_dnapars.log"],shell=True)

    elif buildTree == 'NJ':
        ## biopython tree build
        print("Build NJ tree...")
        # build starting tree (NJ)
        aln = AlignIO.read(timestamp+'.fa', 'fasta')
        calculator = DistanceCalculator('identity')
        constructor = DistanceTreeConstructor(calculator, 'nj')
        treeNJ = constructor.build_tree(aln)
        # Phylo.draw(treeNJ)
        Phylo.write(treeNJ,timestamp+"_NJ.tree","nexus")
        # build parsimonous tree
        #scorer = ParsimonyScorer()
        #searcher = NNITreeSearcher(scorer)
        #constructor = ParsimonyTreeConstructor(searcher, treeNJ)
        #treePS = constructor.build_tree(aln)
        #Phylo.write(treePS,timestamp+"_PS.tree","nexus")
    elif buildTree=='ML':
        # raxml-ng run. # set up: conda install -c bioconda raxml-ng 
        # raxml_parameters='': Pass additional parameters to the raxml command.
        if not root_on_smpl == '': 
            additional_raxml_parameters = additional_raxml_parameters + ' --outgroup ' + root_on_smpl
        ml_prefix = f"{timestamp}_ML_{raxml_model.replace('+', '').replace('{', '').replace('}', '')}"
        subprocess.run([f"raxml-ng --msa {timestamp}.fa --model {raxml_model} --prefix {ml_prefix} {additional_raxml_parameters}"],shell=True)
        if bootstrapping:
            subprocess.run([f"raxml-ng --support --tree {timestamp}_ML_{raxml_model.replace('+', '')}.raxml.bestTreeCollapsed --bs-trees {timestamp}_ML_{raxml_model.replace('+', '')}.raxml.bootstraps --prefix {timestamp}_ML_{raxml_model.replace('+', '')}_with_bootstraps"], shell = True)
    elif buildTree=='ML-basic':
        # raxml-ng run. # set up: conda install -c bioconda raxml 
        # raxml_parameters='': Pass additional parameters to the raxml command.
        if not root_on_smpl == '': 
            additional_raxml_parameters = additional_raxml_parameters.replace('--seed', '-p').replace('--threads', '-T') + ' -o ' + root_on_smpl ## should change to a parameter dict
        subprocess.run([f"raxmlHPC -s {timestamp}.fa -m {raxml_model} -n {timestamp}_ML_{raxml_model.replace('+', '')}.tree {additional_raxml_parameters}"],shell=True)
    return timestamp

def convert_indel_for_trees(size_for_indels,missingdata="?"):
    # translate index array to array containing nucleotides
    all = size_for_indels.ravel().astype(int).astype(str)
    all[np.where(all=='-9223372036854775808')]=missingdata
    return all.reshape(size_for_indels.shape)


def build_table_for_tree_labeling(p_chr_table,treeSampleNamesLong,calls_for_tree,patient="",indels=False):
    """Purges any previous per-SNP trees and recreates folder structure to fill. Creates for_tree_labelling.csv in appropriate locaiton (which is just a SNP table) to color trees in later functions."""
    # make new folder ('tree_counting'), wipe all previous content, add table
    if calls_for_tree.shape[1] != treeSampleNamesLong.shape[0]:
        print("WARNING: inputs for treeSampleNamesLong and calls_for_tree inputs have different dimensions (differing number of samples), downstream errors will occur with indexing")
        print("Aborting function execution.")
        return
    if patient != "":
        subprocess.run(["rm -fr tree_counting/subject"+patient+" ; mkdir tree_counting/subject"+patient ],shell=True)
        with open("tree_counting/subject"+patient+"/for_tree_labeling.csv",'w') as csv_file:
            csv_file.write(",".join(np.append(np.array(['chr','pos']),treeSampleNamesLong))+"\n") # write header
            for i in range(p_chr_table.shape[0]):
                csv_file.write(",".join(np.append( np.array([str(p_chr_table[i,0]),str(p_chr_table[i,1])]) ,calls_for_tree[i,]))+"\n")
    else:
        if indels:
            subprocess.run(["rm -fr tree_counting/indels ; mkdir -p tree_counting/indels" ],shell=True)
        else: 
            subprocess.run(["rm -fr tree_counting/snps ; mkdir -p tree_counting/snps" ],shell=True)
    # build table    
    if indels:
        with open("tree_counting/indels/for_tree_labeling.csv",'w') as csv_file:
            csv_file.write(",".join(np.append(np.array(['chr','pos']),treeSampleNamesLong))+"\n") # write header
            for i in range(p_chr_table.shape[0]):
                csv_file.write(",".join(np.append( np.array([str(p_chr_table[i,0]),str(p_chr_table[i,1])]) ,calls_for_tree[i,]))+"\n")
    else:
        with open("tree_counting/snps/for_tree_labeling.csv",'w') as csv_file:
            csv_file.write(",".join(np.append(np.array(['chr','pos']),treeSampleNamesLong))+"\n") # write header
            for i in range(p_chr_table.shape[0]):
                csv_file.write(",".join(np.append( np.array([str(p_chr_table[i,0]),str(p_chr_table[i,1])]) ,calls_for_tree[i,]))+"\n")

###########################################################################
##################### ancestral reconstruction module #####################
###########################################################################
def create_ancestral_reconstruction(treename,treesampleNamesLong_ref_outgroup,calls_for_tree_ref_outgroup, outdir="ancestral_reconstruction"):
    """Creates internal node fasta SNP alignments from input tree, using treetime"""
    subprocess.run([f"rm -fr {outdir} ; mkdir -p {outdir}" ],shell=True)
    with open(f'{outdir}/snp_table_treetime.fasta','w') as file:
        for index,sample in enumerate(treesampleNamesLong_ref_outgroup):
            str_alignment=''.join(calls_for_tree_ref_outgroup[:,index])
            for_output=str_alignment.replace('?','N')
            file.write(f'>{sample}')
            file.write('\n')
            file.write(f'{for_output}')
            file.write('\n')
    # run treetime
    subprocess.run([f"treetime ancestral --aln {outdir}/snp_table_treetime.fasta --tree {treename} --outdir {outdir}"],shell=True)

def parse_ancestral_fasta_treetime(fasta, node_label_to_use=None, convert_to_ntsidx = True):
    """
    Loads in ancestral sequence fasta output from treetime. Note this output is only the indices that went into the calls for tree (usually treetime)
    """
    if (type(fasta) == str) and node_label_to_use:
        ref=SeqIO.parse(fasta,'fasta')
        for record in ref:
            if (record.id == node_label_to_use) and convert_to_ntsidx:
                records = nts2idx(np.array([x for x in record.seq]))
            elif (record.id == node_label_to_use):
                records = [nt for nt in str(record.seq)]
    return records

def create_full_genome_ancestral_reconstruction_fasta(ancestral_nts_all_p,p,ChrStarts,REFGENOMEFOLDER,outdir="ancestral_reconstruction"):
    """
    Creates new (ref genome sized) fasta with full sequence record of reference genome, with ancestral base calls from above treetime inference.
    
    example running:

    apyt.create_ancestral_reconstruction('Yersinia_pestis_CO92_latest_rerooted_YAC.nwk.tree',treesampleNamesLong_ref_outgroup,calls_for_tree_ref_outgroup)
    ancestral_nti_goodpos=apyt.parse_ancestra_fasta_treetime(f'ancestral_reconstruction/ancestral_sequences.fasta','NODE_0000002')

    ancestral_nti=refnti.copy()
    ancestral_nti[goodpos]=ancestral_nti_goodpos

    apyt.create_ancestral_reconstruction_fasta(ancestral_nti,p,chrStarts,ref_genome_folder)

    """
    # parse ref fasta
    ref_genome_path=str(REFGENOMEFOLDER)+'/genome.fasta'
    fasta_file = glob.glob(ref_genome_path)
    if len(fasta_file) != 1:
        ref_genome_path=str(REFGENOMEFOLDER)+'/genome.fasta.gz'
        fasta_file_gz = glob.glob(ref_genome_path)
        if len(fasta_file_gz) != 1:
            raise ValueError('Either no genome.fasta(.gz) or more than 1 genome.fasta(.gz) file found in ' + REFGENOMEFOLDER)
        else: # genome.fasta.gz
            refgenome = SeqIO.parse(gzip.open(fasta_file_gz[0], "rt"),'fasta')
    else: # genome.fasta
        refgenome = SeqIO.parse(fasta_file[0],'fasta')
    ref_genome_seqs={}
    ref_genome_chr_names={}
    for index_in_record,record in enumerate(refgenome):
        ref_genome_seqs[index_in_record] = str(record.seq)
        ref_genome_chr_names[index_in_record] = [record.id,str(record.description)+' ancestral reconstruction']
    chrpos=p2chrpos(p, ChrStarts) ## 1-indexed chr, 0-indexed pos
    for index_in_p,variant_to_reannotate in enumerate(idx2nts(ancestral_nts_all_p)):
        chr=chrpos[index_in_p][0]-1
        pos=chrpos[index_in_p][1]
        ref_genome_seqs[chr]=ref_genome_seqs[chr][:pos] + variant_to_reannotate + ref_genome_seqs[chr][pos+1:]
    seqs_to_output=[]
    for index_in_dicts in range(len(ref_genome_seqs)):
        seq=ref_genome_seqs[index_in_dicts]
        name=ref_genome_chr_names[index_in_dicts][0]
        description=ref_genome_chr_names[index_in_dicts][1]
        seqs_to_output.append(SeqRecord.SeqRecord(Seq(seq), id=name, description=description))
    SeqIO.write(seqs_to_output, f'{outdir}/ancestral_genome.fasta', "fasta")

###################################################################################################
## Generate trees colored by basecall for each SNP functions ## START ##
###################################################################################################

def generate_mutation_colored_tree_helper_generate_temp_tree(tree, samples_basecalls_csv, outgroup_name, color_lca, label_lca, indels, mutational_events_for_name=[]):
    """
    Generates new trees for each mutation (position) with coloring based on base call at the position
    iterates over all the mutations in samples_basecalls_csv (csv generated in build_table_for_tree_labling)
    samples_basecalls_csv is csv of chrom,pos,samples.... on first line, following lines are chrom,pos, SNP call for each sample


    Changelog:
    06.2022 IL: added docstring, parameter options for: 1) designating outgroup, 2) counting mutation events for naming, 3) coloring reference based on basecall, 4) labelling reference name with basecall
    09.2022 IL: added parsing option for indels
    03.2023 IL: added mutational_events_for_name param which is num_mutational_events :)
    """
    if len(mutational_events_for_name)>0:
        print("Note: Output name of tree will include number of mutational events for this mutation (position), based on current tree structure")
    f=open(samples_basecalls_csv,'r').readlines()
    lca_index_found=False ## flag to find outgroup name only once
    for i, line in enumerate(f[1:]):
        ## start iterating over each mutation (position), skip header
        l=line.strip().split(',')
        if len(l) < 5:
            continue
        chromosome=l[0]
        pos=str(int(l[1])) # turn 0-based pos into 1-based pos.
        if not lca_index_found:
            if len(outgroup_name) > 1:
                lca_index=f[0].strip().split(',').index(outgroup_name)
            else:
                #use first strain as lca
                lca_index=2
            lca_index_found=True
            print("lca index =",lca_index, "corresponding to sample name", f[0].strip().split(',')[lca_index])
        lca=l[lca_index]
        #count mutations
        newtree_nexus_snp_colored = generate_mutation_colored_tree_helper_add_snp_call(tree, samples_basecalls_csv, chromosome, pos,outgroup_name, color_lca, label_lca,indels)
        pos_for_output=str(int(pos)+1)
        if len(mutational_events_for_name)>0:
            num_mutational_events=mutational_events_for_name[i]
        #if a==0:
            #print('NO MUTS:')
            #print(chromosome, pos)
        #save trees
            f1=open(f'{chromosome}_{pos_for_output}_{num_mutational_events}.tree','w')
        else:
            f1=open(chromosome+'_'+pos_for_output+'.tree','w')
        t=open('tempnexus.txt','r').readlines()
        for line in t:
            f1.write(line)

def generate_mutation_colored_tree_helper_add_snp_call(tree, bascalls_per_sample_csv, chromosome, pos, reference_name_in_tree, color_lca, label_lca,indels):
    """
    creates temptree.txt tsv file w/ basecall for each sample, then calls generate_mutation_colored_tree_helper_add_colors_to_basecall_tree to generate new tree file with basecall as color labels

    changelog:
    6.2022 IL: added docstring, removed limit to sample numbers and naming convention required to add sample names to temptree.txt, renamed dict --> bascalls_per_sample_csv for clarity, added parsing for generate_mutation_colored_tree_helper_add_colors_to_basecall_tree
    """
    f=open(bascalls_per_sample_csv).readlines()
    fo=open('temptree.txt','w')
    header=f[0].strip().split(',')
    locations=[]
    for n,i in enumerate(header):
    #        if (i.startswith('S') or i.startswith('D') or i.startswith('R') or i.startswith('P')) and len(i)<100:
        locations.append(n)
    for line in f:
        l=line.strip('\n').split(',')
        if l[0]==chromosome and l[1]==pos:
            for i in locations:
                if len(l) > i and len(l[i])>0:
                    fo.write(header[i]+'\t'+l[i]+'\n')
                else:
                    fo.write(header[i]+'\t?\n')
                #if i > len(l):
                    #print(line, l, i)
            break
    fo.close()
    newtree_nexus_snp_colored = generate_mutation_colored_tree_helper_add_colors_to_basecall_tree(tree,'temptree.txt',reference_name_in_tree, color_lca, label_lca,indels,chromosome,pos)
    return newtree_nexus_snp_colored

def generate_mutation_colored_tree_helper_add_colors_to_basecall_tree(tree, dictionary, reference_name_in_tree, color_lca, label_lca, indels,chromosome,pos):
    """
    generates tempnexus.txt with color labels and tiplabels based on basecall
    
    Changelog:
    09.2022 IL: Added parsing of indels and coloring based on classes, coloring classes will repeat if >6 indel classes exist (unlikely)
    """
    f=open(dictionary).readlines()
    numStrains=len(f)
    annotation={}
    intree_labels_to_outtree_labels_with_snp={}
    #print header for tempnexus
    fo=open('tempnexus.txt','w')
    fo.write('#NEXUS\nbegin taxa;\n\tdimensions ntax='+str(numStrains-2)+';\n\ttaxlabels\n')
    if not indels:
        colors={'A':'[&!color=#-16776961]', 'C':'[&!color=#-16725916]','G':'[&!color=#-3670016]','T':'[&!color=#-3618816]','?':'[&!color=#-16777216]',
                'a':'[&!color=#-16776961]', 'c':'[&!color=#-16725916]','g':'[&!color=#-3670016]','t':'[&!color=#-3618816]'} ## used for pseudogenized_SNVs
    else:
        colors={'0':'[&!color=#-10223416]','?':'[&!color=#-16777216]'}
        ## not all colors will be used, repeats to ensure it has enough anyways
        colors_multiple_indels=['[&!color=#-3670016]','[&!color=#-32768]','[&!color=#-3618816]','[&!color=#-65408]','[&!color=#-16725916]','[&!color=#-16776961]']
    ambiguous=['R','Y','M','K','S','W']
    for each in ambiguous:
        colors[each]='[&!color=#-16777216]'
    indels_colored={}
    #get annotations
    f=open(dictionary, 'r').readlines()
    for line in f:
        if not line.startswith('#'):
            l=line.split()
            annotation[l[0]]=l[1]
    # combine names and annotations, dict sample --> sample--*A/T/C/G
    for i in annotation.keys():
        if not label_lca and i == reference_name_in_tree:
            intree_labels_to_outtree_labels_with_snp[i] = i
        else: 
            intree_labels_to_outtree_labels_with_snp[i] = i + '--*' + annotation[i]
    #make new tree
    f=open(tree,'r').readlines()
    for line in f:
        ## first need order of samples as they appear in tree, to put color labels in correct order
        sample_order_in_tree=[]
        for sample in line.strip().replace('(','').replace(')','').split(','): ## each sample will have own string, with other info from previous tree
            sample_name=sample.split(':')[0] ## separate sample_name from info from tree (eg relatedness, not important for us right now)
            sample_order_in_tree.append(intree_labels_to_outtree_labels_with_snp[sample_name]) ## add to list of sample order in tree
            samplename_with_basecall=sample_order_in_tree[-1] ## samplename just added
            if not indels:
                if not color_lca and sample_name == reference_name_in_tree: 
                    fo.write('\t\''+samplename_with_basecall+'\'[&!color=#-16777216]\n')
                else:
                    fo.write('\t\''+samplename_with_basecall+'\''+colors[samplename_with_basecall[-1]]+'\n') #write down new colors to nexus file
            else:
                if not color_lca and sample_name == reference_name_in_tree: 
                    fo.write('\t\''+samplename_with_basecall+'\'[&!color=#-16777216]\n')
                else:
                    if samplename_with_basecall[-1] not in colors:
                        if samplename_with_basecall[-1] in indels_colored:
                            color = indels_colored[samplename_with_basecall[-1]]
                        else:
                            index_to_add=len(indels_colored) % len(colors_multiple_indels)
                            if len(indels_colored) > len(colors_multiple_indels): print(f"Coloring classes of indels will repeat for tree {chromosome}_{pos}")
                            indels_colored[samplename_with_basecall[-1]]=colors_multiple_indels[index_to_add]
                            color = indels_colored[samplename_with_basecall[-1]]
                        fo.write('\t\''+samplename_with_basecall+'\''+color+'\n')
                    else:
                        fo.write('\t\''+samplename_with_basecall+'\''+colors[samplename_with_basecall[-1]]+'\n')
        fo.write(';\nend\n\nbegin trees;\n\ttree tree_1=[&R] ')
        ## next, need to replace all samplenames in tree file with newly updated labels with colors, as above
        for intree_labels,outtree_labels in intree_labels_to_outtree_labels_with_snp.items():
            line=line.replace(intree_labels,outtree_labels)
        newtree_nexus_snp_colored=line
    for line in newtree_nexus_snp_colored:
        fo.write(line)
    fo.write('end;\n\nbegin figtree;\n')
    fo.write('\tset tipLabels.fontSize=10;\n')
    fo.write('end;')
    fo.close()
    return newtree_nexus_snp_colored

def reroot_tree(tree_path, outgroup_name_for_lca = "", tree_type = 'newick', verbose = False):
    if verbose: print('Note: Uses first strain as LCA, unless outgroup_name_for_lca provided as a string.')
    if outgroup_name_for_lca != "": print("Outgroup name provided, will attempt to use sample", outgroup_name_for_lca, "as outgroup.")
    rc('font', **{'sans-serif':'Arial'})
    if verbose: print("Attempting to reroot tree, based on provided outgroup_name_for_lca",outgroup_name_for_lca)
    input_tree=Phylo.read(tree_path, tree_type)
    if outgroup_name_for_lca:
        input_tree.root_with_outgroup(outgroup_name_for_lca)
        tree_path=tree_path.replace("_latest","_latest_rerooted_"+outgroup_name_for_lca) 
    else:
        input_tree.root_at_midpoint()
        tree_path=tree_path.replace("_latest","_latest_rerooted_midpoint")        
    Phylo.write(input_tree, tree_path, tree_type)
    if verbose: 
        print(f"\nSuccessfully rerooted tree to filename",tree_path)
    return tree_path

def generate_mutation_colored_tree(tree,samples_basecalls_csv, outgroup_name_for_lca="", count_mutations_for_name=[], color_lca=False, label_lca=True, indels=False):
    """
    Main function for generating trees for every SNP
    Calls generate_mutation_colored_tree_helper_generate_temp_tree and sets font. 
    Set values for downstream helper function calls in this function. 
    
    Changelog:
    06.2022 IL: added additional print statements for flagging information, expected output behavior, added additional flags for modifying output
    06.2022 IL: added optional rerooting on all SNP trees to user input LCA, using biopython's Phylo, modified parameter names to lca instead of reference, outputting parameters used for run
    09.2022 IL: added parsing option for indel calls.
    03.2023 IL: externalized rerooting of tree to separate function, to be called prior to coloring tree.
    """
    print("Note: Coloring last common ancestor (lca) tip by basecall option is currently",color_lca, "... set color_lca to True/False to modify")
    print("Note: Labelling last common ancestor with basecall is currently", label_lca, "... set label_lca to True/False to modify")
    # save parameters into output directory
    paramters_df=pd.DataFrame({"parameter_name":list(locals().keys()),"values_input":list(locals().values())})
    paramters_df.to_csv("parameters.tsv",sep="\t",index=False)
    ## print parameters and log for start, end at runtime
    if indels:
        print("Generating per-indel trees, uncalled black, indel as red, no indel call as green.")
    else:
        print("Generating per-snp trees, using basecalls as color.")
    generate_mutation_colored_tree_helper_generate_temp_tree(tree,samples_basecalls_csv,outgroup_name_for_lca, color_lca, label_lca, indels, count_mutations_for_name)
    if indels:
        print("Successfully generated per-indel trees, uncalled black, Reference (no indel) as purple, other indel classes as other colors.")
    else:
        print('Successfully generated per-snp trees, colored by mutation basecall')

###################################################################################################
## Generate trees colored by basecall for each SNP functions ## END ##
###################################################################################################

def generate_colored_tree(input_tree,samples_to_color, default_color='[&!color=#-3670016]',reroot=''):
    """will recolor tree on samples_to_color (either list/nparray of samplenames or dict mapping samplenames --> colors"""
    #NOTE: not fully implemented / bugging out atm, TODO for hackathon
    colors_dict={}
    if not isinstance(samples_to_color,dict):
        colors_dict={sample: default_color for sample in samples_to_color}
    if reroot!='':
        phylo_tree=Phylo.read(input_tree, "newick")
        phylo_tree.root_with_outgroup(reroot)
        input_tree=input_tree.replace('_latest','_latest_rerooted')
        Phylo.write(phylo_tree, input_tree, "newick")
        
    f=open(input_tree,'r').readlines()
    numStrains=f[0].count(',')+1
    #print header for colored_tree
    output_tree=input_tree.replace('_latest','_colored')
    fo=open(output_tree,'w')
    fo.write('#NEXUS\nbegin taxa;\n\tdimensions ntax='+str(numStrains)+';\n\ttaxlabels\n')

    for line in f:
        ## first need order of samples as they appear in tree, to put color labels in correct order
        sample_order_in_tree=[]
        for sample in line.strip().replace('(','').replace(')','').split(','): ## each sample will have own string, with other info from previous tree
            sample_name=sample.split(':')[0] ## separate sample_name from info from tree (eg relatedness, not important for us right now)
            sample_order_in_tree.append(sample_name) ## add to list of sample order in tree
            samplename_with_basecall=sample_order_in_tree[-1] ## samplename just added
            if sample_name in colors_dict: 
                fo.write('\t\''+samplename_with_basecall+'\''+colors_dict[sample_name]+'\n')
            else:  
                fo.write('\t\''+samplename_with_basecall+'\'[&!color=#-16777216]\n')
        fo.write(';\nend\n\nbegin trees;\n\ttree tree_1=[&R] ')
        ## next, need to replace all samplenames in tree file with newly updated labels with colors, as above
        newtree_nexus_snp_colored=line
    for line in newtree_nexus_snp_colored:
        fo.write(line)
    fo.write('end;\n\nbegin figtree;\n')
    fo.write('\tset tipLabels.fontSize=10;\n')
    fo.write('end;')
    fo.close()


def build_dataframe_coverage_info(goodpos2useTree,NTs,SampleNames,maNT,minorNT,coverage_forward_strand,coverage_reverse_strand,maf,minorAF):
    ## Very slow, numpy integration would be great!

    listOfDF_fwd = []
    listOfDF_rev = []
    NTs_wN = np.append(NTs,'N') # there are N in maf/ninorNT
    for i in goodpos2useTree: # array contains index of p
        lod_fwd = []
        lod_rev = []
        # extract for each pos-sample pair the read_count for ACTG (we record maf & minor >> thus max 2 nucl per pos)
        for j,sample in enumerate(SampleNames):
            spl_fwd = {}
            spl_rev = {}
            spl_fwd['sample'] = sample
            spl_rev['sample'] = sample
            for k,nuc in enumerate(NTs_wN):
                # dict had to be defined within if else, otherwise get overwritten
                if k == maNT[i,j]:
                    # print('maNT',k,nuc)
                    spl_fwd[nuc] = ( coverage_forward_strand[i,j] * maf[i,j] ) # recalc #Nucleotide calls based on coverage and major allele freq
                    spl_rev[nuc] = ( coverage_reverse_strand[i,j] * maf[i,j] ) # recalc #Nucleotide calls based on coverage and major allele freq
                elif k == minorNT[i,j]:
                    # only report minorAF when majorNT not set 4 during filtering!
                    if maNT[i,j] != 4:
                        spl_fwd[nuc] = ( coverage_forward_strand[i,j] * minorAF[i,j] ) # recalc #Nucleotide calls based on coverage and major allele freq
                        spl_rev[nuc] = ( coverage_reverse_strand[i,j] * minorAF[i,j] ) # recalc #Nucleotide calls based on coverage and major allele freq
                else:
                    spl_fwd[nuc] = 0
                    spl_rev[nuc] = 0
            lod_fwd.append(spl_fwd)
            lod_rev.append(spl_rev)
            df_fwd = pd.DataFrame(lod_fwd)
            df_fwd = df_fwd.set_index('sample')
            df_rev = pd.DataFrame(lod_rev)
            df_rev = df_rev.set_index('sample')
        listOfDF_fwd.append(df_fwd)
        listOfDF_rev.append(df_rev)
    return [listOfDF_fwd,listOfDF_rev]

def plot_coverage_fwd_rev_stacked_red(chr_pos_gp,lod_fwd_cov,lod_rev_cov,timestamp,subject=""):
    ## BUG FIX REQUEST: Use counts for barplot!!! More informative about nucleotide info at site!!!
    # reconstruct TDL coverage plot for every SNP across all samples
    # pdfs: pdf/coverage_snp_fwd_rev/timestamp
    subprocess.run(["mkdir -p pdf/coverage_snp_fwd_rev/" + timestamp + "_"+subject+ " ; rm -f pdf/coverage_snp_fwd_rev/" + timestamp + "_"+subject + "/* "],shell=True)
    os.chdir('pdf/coverage_snp_fwd_rev/' + timestamp+ "_"+subject)

    ## 2022.04.08: old version running as one process
    for i in range(chr_pos_gp.shape[0]):
        # get pos,chr,locID,N/S/I/P, FQ for each positon and plot
        chr = chr_pos_gp[i,0]
        pos = chr_pos_gp[i,1] + 1 # turn 1-based bcs annotation_mutation is 1-based
        identifier_string = str(chr) + "_" + str(pos)
        print("Plot coverage fwd/rev: "+identifier_string )
        plot_stacked_paired( [lod_fwd_cov[i],lod_rev_cov[i]] , labels=["fwd", "rev"] , title=identifier_string )
    os.chdir('../../../')


def plot_coverage_fwd_rev_stacked(chr_pos_gp,annotation_mutations,lod_fwd_cov,lod_rev_cov,timestamp,subject=""):
    ## BUG FIX REQUEST: Use counts for barplot!!! More informative about nucleotide info at site!!!
    # reconstruct TDL coverage plot for every SNP across all samples
    # pdfs: pdf/coverage_snp_fwd_rev/timestamp
    subprocess.run(["mkdir -p pdf/coverage_snp_fwd_rev/" + timestamp + "_"+subject+ " ; rm -f pdf/coverage_snp_fwd_rev/" + timestamp + "_"+subject + "/* "],shell=True)
    os.chdir('pdf/coverage_snp_fwd_rev/' + timestamp+ "_"+subject)

    ## 2022.04.08: old version running as one process
    for i in range(chr_pos_gp.shape[0]):
        # get pos,chr,locID,N/S/I/P, FQ for each positon and plot
        chr = chr_pos_gp[i,0]
        pos = chr_pos_gp[i,1] + 1 # turn 1-based bcs annotation_mutation is 1-based
        bool_pos_anno_df = (annotation_mutations['chr'] == chr) &  (annotation_mutations['pos'] == pos)
        if annotation_mutations.loc[ bool_pos_anno_df , 'locustag' ].isnull().any(axis=0):
            locID = "nan"
        else:
            locID = annotation_mutations.loc[ bool_pos_anno_df , 'locustag' ].values[0] # locustag
        anno = annotation_mutations.loc[ bool_pos_anno_df , 'type' ].values[0] # N/S/P/I
        qual = str(int(annotation_mutations.loc[ bool_pos_anno_df , 'quals' ].values[0])) # single qual value (based on FQ: best of lowest scoring pair)
        identifier_string = str(chr) + "_" + str(pos) + "_" + locID + "_" + anno + "_" + qual        
        print("Plot coverage fwd/rev: "+identifier_string )
        plot_stacked_paired( [lod_fwd_cov[i],lod_rev_cov[i]] , labels=["fwd", "rev"] , title=identifier_string )
    os.chdir('../../../')

def plot_coverage_fwd_rev_stacked_multi(chr_pos_gp,annotation_mutations,lod_fwd_cov,lod_rev_cov,timestamp,subject="",no_processes=multiprocessing.cpu_count()-1):
    ## BUG FIX REQUEST: Use counts for barplot!!! More informative about nucleotide info at site!!!
    # reconstruct TDL coverage plot for every SNP across all samples
    # pdfs: pdf/coverage_snp_fwd_rev/timestamp
    subprocess.run(["mkdir -p pdf/coverage_snp_fwd_rev/" + timestamp + "_"+subject+ " ; rm -f pdf/coverage_snp_fwd_rev/" + timestamp + "_"+subject + "/* "],shell=True)
    os.chdir('pdf/coverage_snp_fwd_rev/' + timestamp+ "_"+subject)

    tasks = []
    for i in range(chr_pos_gp.shape[0]):
        # get pos,chr,locID,N/S/I/P, FQ for each positon and plot
        chr = chr_pos_gp[i,0]
        pos = chr_pos_gp[i,1] + 1 # turn 1-based bcs annotation_mutation is 1-based
        bool_pos_anno_df = (annotation_mutations['chr'] == chr) &  (annotation_mutations['pos'] == pos)
        if annotation_mutations.loc[ bool_pos_anno_df , 'locustag' ].isnull().any(axis=0):
            locID = "nan"
        else:
            locID = annotation_mutations.loc[ bool_pos_anno_df , 'locustag' ].values[0] # locustag
        anno = annotation_mutations.loc[ bool_pos_anno_df , 'type' ].values[0] # N/S/P/I
        qual = str(int(annotation_mutations.loc[ bool_pos_anno_df , 'quals' ].values[0])) # single qual value (based on FQ: best of lowest scoring pair)
        identifier_string = str(chr) + "_" + str(pos) + "_" + locID + "_" + anno + "_" + qual        
        tasks.append([ [lod_fwd_cov[i],lod_rev_cov[i]],["fwd", "rev"],identifier_string ])
    # Create a pool of workers
    with concurrent.futures.ThreadPoolExecutor(max_workers=no_processes) as executor:
        # Submit tasks to the executor
        futures = [executor.submit(plot_stacked_paired, *task) for task in tasks]
        
        # Optionally, wait for all futures to complete
        concurrent.futures.wait(futures)
    os.chdir('../../../')

def plot_coverage_indel_ref_alt_stacked(chr_pos_gp,annotation_mutations,lod_fwd_cov,lod_rev_cov,timestamp,subject=""):
    ## BUG FIX REQUEST: Use counts for barplot!!! More informative about nucleotide info at site!!!
    # reconstruct TDL coverage plot for every SNP across all samples
    # pdfs: pdf/coverage_snp_fwd_rev/timestamp
    subprocess.run(["mkdir -p pdf/coverage_indel_ref_alt/" + timestamp + "_"+subject+ " ; rm -f pdf/coverage_indel_ref_alt/" + timestamp + "_"+subject + "/* "],shell=True)
    os.chdir('pdf/coverage_indel_ref_alt/' + timestamp+ "_"+subject)

    ## 2022.04.08: old version running as one process
    for i in range(chr_pos_gp.shape[0]):
        # get pos,chr,locID,N/S/I/P, FQ for each positon and plot
        chr = chr_pos_gp[i,0]
        pos = chr_pos_gp[i,1] + 1 # turn 1-based bcs annotation_mutation is 1-based
        bool_pos_anno_df = (annotation_mutations['chr'] == chr) &  (annotation_mutations['pos'] == pos)
        if annotation_mutations.loc[ bool_pos_anno_df , 'locustag' ].isnull().any(axis=0):
            locID = "nan"
        else:
            locID = annotation_mutations.loc[ bool_pos_anno_df , 'locustag' ].values[0] # locustag
        anno = annotation_mutations.loc[ bool_pos_anno_df , 'type' ].values[0] # N/S/P/I
        qual = str(int(annotation_mutations.loc[ bool_pos_anno_df , 'quals' ].values[0])) # single qual value (based on FQ: best of lowest scoring pair)
        identifier_string = str(chr) + "_" + str(pos) + "_" + locID + "_" + anno + "_" + qual        
        print("Plot coverage fwd/rev: "+identifier_string )
        plot_stacked_paired( [lod_fwd_cov[i],lod_rev_cov[i]] , labels=["fwd", "rev"] , title=identifier_string )
    os.chdir('../../../')


def plot_input(chr_pos_gp, annotation_mutations, lod_fwd_cov, lod_rev_cov, i):
    # get pos,chr,locID,N/S/I/P, FQ for each positon and plot
    chr = chr_pos_gp[i,0]
    pos = chr_pos_gp[i,1] + 1 # turn 1-based bcs annotation_mutation is 1-based
    bool_pos_anno_df = (annotation_mutations['chr'] == chr) &  (annotation_mutations['pos'] == pos)
    if annotation_mutations.loc[ bool_pos_anno_df , 'locustag' ].isnull().any(axis=0):
        locID = "nan"
    else:
        locID = annotation_mutations.loc[ bool_pos_anno_df , 'locustag' ].values[0] # locustag
    anno = annotation_mutations.loc[ bool_pos_anno_df , 'type' ].values[0] # N/S/P/I
    qual = str(int(annotation_mutations.loc[ bool_pos_anno_df , 'quals' ].values[0])) # single qual value (based on FQ: best of lowest scoring pair)
    identifier_string = str(chr) + "_" + str(pos) + "_" + locID + "_" + anno + "_" + qual        
    print("Plot coverage fwd/rev: "+identifier_string )
    plot_stacked_paired( [lod_fwd_cov[i],lod_rev_cov[i]] , labels=["fwd", "rev"] , title=identifier_string )


def plot_stacked_paired(dfall, labels=None, title="notitle",  H="/", **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
        labels is a list of the names of the dataframe, used for the legend
        title is a string for the title of the plot
        H is the hatch used for identification of the different dataframe"""

    """Please check if columns are always the same! Make this as a condition maybe!
    Make also multithreading!"""

    #NOTE: when major NT == 'N', minor allele frequency not shown. Column only represents major AF in this specific case. Otherwise column presents majorAF plus minorAF
    n_df = len(dfall)
    n_col = len(dfall[0].columns) 
    n_ind = len(dfall[0].index)
    
    # define plot width, min = 7
    plot_width = int(n_ind/5)
    if plot_width < 7:
        plot_width = 7
    elif plot_width > 200:
        plot_width = 200
    fig = plt.figure(figsize=(plot_width, 8))
    axe = plt.subplot(111)

    for df in dfall : # for each data frame

        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ylim=(0,100),
                      ax=axe,
                      legend=False,
                      grid=False,
                      **kwargs)  # make bar plots
        
    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                rect.set_hatch(H * int(i / n_col)) #edited part     
                rect.set_width(1 / float(n_df + 1))

    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    axe.set_xticklabels(df.index, rotation = 90)
    axe.set_title(title)

    # Add invisible data to add another legend
    n=[]        
    for i in range(n_df):
        n.append(axe.bar(0, 0, color="gray", hatch=H * i))

    l1 = axe.legend(h[:n_col], l[:n_col], loc=[1.01, 0.5])
    if labels is not None:
        l2 = plt.legend(n, labels, loc=[1.01, 0.1]) 
    axe.add_artist(l1)
    plt.tight_layout()

    fig.savefig(title + ".pdf")
    plt.close()


def plot_stacked_paired_multi(dfall, labels=None, title="notitle",  H="/", **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
        labels is a list of the names of the dataframe, used for the legend
        title is a string for the title of the plot
        H is the hatch used for identification of the different dataframe"""

    """Please check if columns are always the same! Make this as a condition maybe!
    Make also multithreading!"""

    #NOTE: when major NT == 'N', minor allele frequency not shown. Column only represents major AF in this specific case. Otherwise column presents majorAF plus minorAF
    n_df = len(dfall)
    n_col = len(dfall[0].columns) 
    n_ind = len(dfall[0].index)
    
    # define plot width, min = 7
    plot_width = int(n_ind/5)
    if plot_width < 7:
        plot_width = 7
    elif plot_width > 200:
        plot_width = 200
    fig = plt.figure(figsize=(plot_width, 8))
    axe = plt.subplot(111)

    for df in dfall : # for each data frame

        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ylim=(0,100),
                      ax=axe,
                      legend=False,
                      grid=False,
                      **kwargs)  # make bar plots
        
    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                rect.set_hatch(H * int(i / n_col)) #edited part     
                rect.set_width(1 / float(n_df + 1))

    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    axe.set_xticklabels(df.index, rotation = 90)
    axe.set_title(title)

    # Add invisible data to add another legend
    n=[]        
    for i in range(n_df):
        n.append(axe.bar(0, 0, color="gray", hatch=H * i))

    l1 = axe.legend(h[:n_col], l[:n_col], loc=[1.01, 0.5])
    if labels is not None:
        l2 = plt.legend(n, labels, loc=[1.01, 0.1]) 
    axe.add_artist(l1)

    fig.savefig(title + "_multi.pdf")
    plt.close()
    pause(0.1)

def plot_indels_cov_iaf_ref_alt(sampleNames, cand_indel, indel_p, indel_depth, indel_call, annotation_indels, chrStarts, path_ext):
        subprocess.run(["mkdir -p pdf/cov_iaf_ref_alt_indel/" + path_ext + " ; rm -f pdf/cov_iaf_ref_alt_indel/" + path_ext + "/* "],shell=True)
        os.chdir('pdf/cov_iaf_ref_alt_indel/' + path_ext)

        plot_width = int(len(sampleNames)/5)
        if plot_width < 7:
            plot_width = 7
            
        for indel_idx in cand_indel:
            ## get unique alternative calls
            alt_indel_sample = indel_call[indel_idx, :, 1]
            unique_alt_indel = np.unique([alt for alt in alt_indel_sample if (alt == alt) and (alt != None)])
            unique_alt_indel_dict = {alt: idx+1 for idx, alt in enumerate(unique_alt_indel)}
            ## generate unique colormap for number of alternative alleles
            alt_hsv_color = mcolors.rgb_to_hsv(sns.color_palette()[1])
            saturation_interspace = 1/(len(unique_alt_indel) - 0.5)
            unique_ab_class_color = {}
            unique_ab_class_color['ref'] = np.array(sns.color_palette()[0])
            for alt, alt_id in unique_alt_indel_dict.items():
                unique_ab_class_color[f'alt_{alt_id}'] = mcolors.hsv_to_rgb((alt_hsv_color[0], alt_hsv_color[1] - alt_hsv_color[1]*saturation_interspace*(alt_id-1), alt_hsv_color[2]))
            unique_ab_class_color['undefined'] = np.array([0.8, 0.8, 0.8]) ## lightgrey
            ## prepare df for plotting
            cov_df = pd.DataFrame(indel_depth[indel_idx, :, :], columns = ['ref', 'alt', 'total'])
            cov_df['undefined'] = cov_df['total'] - (cov_df['ref']+cov_df['alt'])
            cov_df['sample'] = sampleNames
            ## inset alternatve alleles and spread df
            cov_df['alt_allele'] = alt_indel_sample
            cov_df['alt_allele'] = cov_df['alt_allele'].map(unique_alt_indel_dict)
            cov_df = cov_df.pivot(index=['ref', 'total', 'undefined', 'sample'], columns='alt_allele', values='alt')
            alt_cols = [f'alt_{int(subset)}' if (subset == subset) else 'alt' for subset in cov_df.columns ]
            cov_df.columns = alt_cols
            cov_df = cov_df.reset_index()
            ## sort df
            cov_df = cov_df[['ref'] + alt_cols + ['undefined', 'total', 'sample']]
            cov_df['sample'] = pd.Categorical(cov_df['sample'], categories=sampleNames, ordered=True)
            cov_df = cov_df.sort_values('sample').reset_index(drop = True)

            ## calculate frequencies
            cov_df_freq = cov_df.copy()
            for col in ['ref', 'undefined'] + alt_cols:
                cov_df_freq[col] = cov_df_freq[col] / cov_df_freq['total']
            cov_df = cov_df.drop(['total', 'alt'], axis = 1)
            cov_df_freq = cov_df_freq.drop(['total', 'alt'], axis = 1)
            
            cov_df_long = cov_df.melt('sample')
            
            indel_chr_pos = p2chrpos(np.array([indel_p[indel_idx]]), chrStarts)
            for col in ['gene', 'gene1', 'gene2']:
                if not col in annotation_indels.columns:
                    annotation_indels[col] = ''
            anno = annotation_indels[(annotation_indels['chr'] == indel_chr_pos[0][0]) & (annotation_indels['pos'] == indel_chr_pos[0][1]+1)]
            anno.loc[anno['gene'].isna(), 'gene'] = anno.loc[anno['gene'].isna(), 'gene1'] + ';' + anno.loc[anno['gene'].isna(), 'gene2']
            titlestring = f'{indel_chr_pos[0][0]}_{indel_chr_pos[0][1]+1}_{anno["gene"].values[0]}_{anno["max_indel_GL_diff_to_ref"].values[0]}'
            
            ## plot absolute values (coverage)
            plot_cov_per_sample_bar(cov_df_long, plot_width, unique_ab_class_color, titlestring)

            ## plot frequencies (allele frequencies)
            plot_iaf_per_sample_bar(cov_df_freq, sampleNames, plot_width, unique_ab_class_color, titlestring)

            ## plot alt allele frequency against alt coverage
            plot_iaf_alt_cov(indel_depth, indel_idx, titlestring)

        os.chdir('../../../')

def plot_cov_per_sample_bar(cov_df_long, plot_width, unique_ab_class_color, titlestring):
    fig, ax = plt.subplots(figsize = (plot_width, 8))
    sns.barplot(data = cov_df_long, x = 'sample', y = 'value', hue = 'variable', palette = unique_ab_class_color, linewidth = 0, dodge=True, ax = ax)
    plt.xticks(rotation = 90)
    ax.set_ylim(0, 100)
    ax.set_ylabel('Allele Coverage')
    ax.legend(title = 'Allele', loc='center left', bbox_to_anchor=(1, 0.5), fancybox=False, shadow=False)
    ax.set_title(titlestring)
    plt.tight_layout()
    fig.savefig(titlestring + "_cov.pdf")
    plt.close('all')

def plot_iaf_per_sample_bar(cov_df_freq, sampleNames, plot_width, unique_ab_class_color, titlestring):
    fig, ax = plt.subplots(figsize = (plot_width, 8))
    ax = cov_df_freq.plot(kind="bar",
                            stacked=True,
                            width = 0.8,
                            edgecolor = 'k', 
                            color = [unique_ab_class_color[col] for col in cov_df_freq.columns if col in unique_ab_class_color.keys()], ## important that df and dict has same order --> loop to extract proper colors!
                            linewidth = 0.2,
                            ylim=(0,1),
                            ax=ax,
                            grid=False)  # make stacked bar plots of allele frequencies
    ax.set_xticklabels(sampleNames, rotation = 90)
    ax.set_ylim(0, 1)
    ax.set_ylabel('Allele frequency')
    ax.legend(title = 'Allele', loc='center left', bbox_to_anchor=(1, 0.5), fancybox=False, shadow=False)
    ax.set_title(titlestring)
    plt.tight_layout()
    fig.savefig(titlestring + "_iaf.pdf")
    plt.close('all')

def plot_iaf_alt_cov(indel_depth, indel_idx, titlestring):
    fig, ax = plt.subplots(figsize = (4, 4))
    sns.scatterplot(x = indel_depth[indel_idx, :, 1], y = indel_depth[indel_idx, :, 1]/indel_depth[indel_idx, :, 2], alpha = 0.25, edgecolor = 'k', linewidth = 0.2, ax = ax)
    ax.set_ylabel('Alternative allele frequency')
    ax.set_xlabel('Alternative allele coverage')
    ax.set_title(titlestring)
    plt.tight_layout()
    fig.savefig(titlestring + "_cov_vs_iaf.pdf")
    plt.close('all')
        

def homoplasy_simulation(annotation_mutations,numtrials,ref_genome_folder,plot=True,plot_title=None,output_name=None,num_mutations=None):
    """Calculates the expected number of homoplasmic positions and probability of observing the actual number of homoplasmic positions, based on simulation of all pos equally likely to be mutated. Plotting optional."""
    [chrstarts, genomelength,scafnames]= genomestats(ref_genome_folder);
    expectedHomoplasmicPos = np.zeros(numtrials)
    num_homoplasies=0
    if not num_mutations:
        if 'num_mutational_events' in annotation_mutations.columns:
            num_mutations=np.sum(annotation_mutations['num_mutational_events'])
            num_homoplasies=len(np.where(annotation_mutations['num_mutational_events']>1)[0])
        else:
            num_mutations=len(annotation_mutations)
    for i in range(numtrials):
        # Pick random positions on the genome to mutate
        randpos = np.random.randint(1, int(genomelength), int(num_mutations) , dtype=int)
        expectedHomoplasmicPos[i]=np.sum(np.unique(randpos,return_counts=True)[1]>1)
    simProbForObsCand = 1-sum(expectedHomoplasmicPos < num_homoplasies)/len(expectedHomoplasmicPos)
    if plot:
        plt.rcParams.update({'font.size': 14}) # all label size
        f = figure()
        ax=f.add_subplot(111)
        # histogram depicting simulated distribution of expected candidate genes, filter criteria for candidates, P for observation cand count
        mybin = np.arange(max(np.unique(expectedHomoplasmicPos))+2)-0.5 # +2 in order to now loose last bin due to arange; -0.5 needed for bars align at center with x axis ticks
        plt.hist(expectedHomoplasmicPos,bins=mybin,rwidth=0.8,color='#607c8e')
        max_on_x_axis=max(max(expectedHomoplasmicPos),num_homoplasies)
        if max_on_x_axis>100:
            step=20
        elif max_on_x_axis>50:
            step=10
        elif max_on_x_axis>20:
            step=5
        else:
            step=1
        plt.xlim(0,max(max(expectedHomoplasmicPos),num_homoplasies+1))
        xticks(np.arange(0,max(max(expectedHomoplasmicPos),num_homoplasies)+step,step))
        plt.ylabel('Simulated counts, N='+str(numtrials))
        plt.xlabel('Number of homoplasies')
        if plot_title:
            plt.title(plot_title)
        plt.axvline(x=num_homoplasies,color='violet') 
        text(0.98, 0.88, "P("+ str(num_homoplasies) + ") = " + str(np.around(simProbForObsCand,3)), fontsize=12,horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
        #save pdf
        if output_name:
            subprocess.run(["mkdir -p pdf/adaptive_evo/ "],shell=True)
            f.savefig('pdf/adaptive_evo/' + output_name + ".pdf",bbox_inches="tight")
            plt.close()
            print('Plotted: pdf/adaptive_evo/' + output_name + ".pdf")
        else:
            subprocess.run(["mkdir -p pdf/adaptive_evo/ "],shell=True)
            f.savefig('pdf/adaptive_evo/' + 'homoplasies' + ".pdf",bbox_inches="tight")
            plt.close()
            print('Plotted: pdf/adaptive_evo/' + 'homoplasies' + ".pdf")
    return [expectedHomoplasmicPos,simProbForObsCand]


def mutation_distribution_simulation(num_mutations_genome,numtrials, ref_genome_folder, annotation_genes, prob_N, promotersize=500):
    """
    gets simulated distribution of expected numbers of genic, nongenic, promoter (defined by promotersize bases upstream), and nonsyn mutations (defined by genic * prob_N) 
    """
    [chrstarts, genomelength,scafnames]= genomestats(ref_genome_folder);
    [locus_tag,chr_tag] = genomic_position_all(annotation_genes, genomelength, chrstarts) #loc_tag numbers all genes(CDS) across genome (intergenic/rna:0.5); chr_tag: numbers genes per chr (intergenic: #chr+0.5)
    #    locus_tag >> cds_number_by_position_tags
    #    chr_tag >> cds_number_by_position_indices
    # extract vector genelengths that follows order
    
    annotation_genes_sglDF = pd.concat(annotation_genes,sort=False)
    start_pos_gene = annotation_genes_sglDF['loc1'].values
    start_pos_gene = start_pos_gene[~np.isnan(start_pos_gene)] # remove nan, caused by empty chr
    end_pos_gene = annotation_genes_sglDF['loc2'].values
    end_pos_gene = end_pos_gene[~np.isnan(end_pos_gene)] # remove nan, caused by empty chr
    # gene lengths. incl rna genes (however when used those indexes not present in locus_tag (0.5))
    genelengths = (end_pos_gene - start_pos_gene)
    
    expectedNumberPromoterMutations = np.zeros(numtrials)
    expectedNumberNonsynMutations = np.zeros(numtrials)
    expectedNumberIntergenicMutations_nonpromoter = np.zeros(numtrials)
    expectedNumberIntragenicMutations = np.zeros(numtrials)
    ## get p location of promoters
    promoter_indices=set()
    for chr_num in range(len(chrstarts)):
        p_anno = annotation_genes[chr_num-1]
        neg_strand_locs_starts=p_anno[p_anno['strand']==-1]['loc2'] + 1 ## since 1 range starts at start index, promoter range must be pos+1 .. pos+1+promotersize
        neg_strand_locs_end=neg_strand_locs_starts+promotersize
        ## create dict ranges for all promoter lengths
        for s,e in zip(neg_strand_locs_starts,neg_strand_locs_end):
            for pos_to_add in range(s,e):
                promoter_indices.add(pos_to_add)
        pos_strand_locs_ends=p_anno[p_anno['strand']==1]['loc1'] ## since range ends end index, promoter range must be pos+1 .. pos+1+promotersize
        pos_strand_locs_starts=pos_strand_locs_ends-promotersize
        for s,e in zip(pos_strand_locs_starts,pos_strand_locs_ends):
            for pos_to_add in range(s,e):
                promoter_indices.add(pos_to_add)
    # start sim
    for i in range(numtrials):
        # Pick random positions on the genome to mutate
        randpos = np.random.randint(1, genomelength, num_mutations_genome , dtype=int)
        # Does NOT assume any mutational spectrum!!!
            
        # Find out in which genes, promoter, intergenic regions these mutations occurred
        genenums = locus_tag[randpos];
        genenums_intragenic = genenums[genenums>0.5]; #remove intergenic/rna-gene mutations
        expectedNumberIntragenicMutations[i]=len(genenums_intragenic)

        expectedNumberNonsynMutations[i]=stats.binom.rvs(len(genenums_intragenic),prob_N)
        # Find out promoter and other intergenic
        total_num_promoter_muts=0
        intergenic_pos=randpos[np.where(genenums==0.5)[0]]
        for pos in intergenic_pos:
            if pos in promoter_indices:
                total_num_promoter_muts+=1
        expectedNumberPromoterMutations[i]=total_num_promoter_muts
        expectedNumberIntergenicMutations_nonpromoter[i]=len(intergenic_pos)-total_num_promoter_muts
    return [expectedNumberPromoterMutations,expectedNumberNonsynMutations,expectedNumberIntergenicMutations_nonpromoter,expectedNumberIntragenicMutations]

    
def parallel_evolution_counting_and_simulation(num_mutations_genome, num_mutations_genic , mutation_number_threshold , mutation_density_threshold , numtrials, max_muts_per_gene_to_track, chr_pos_gp, ref_genome_folder, annotation_genes):
    ## Parallel evolution counting and simulation
    # Inputs: genome, number of expected mutations in the genome
    # Output: number of genes with a mutation density (mutations per length of
    # gene) over some threshold

    # num_mutations_genome=number_of_mutations_on_genome
    # num_mutations_genic=number_of_mutations_genic
    # mutation_number_threshold=params_dict['Min_num_mutations_cand']
    # mutation_density_threshold=params_dict['Min_mutation_density_cand']
    # numtrials=params_dict['NumTrialsSim']
    # max_muts_per_gene_to_track = params_dict['max_muts_per_gene_to_track']
    # chr_pos_gp=chr_pos_gp
    # ref_genome_folder=params_dict['ref_genome_folder']
    


    # Get data structures for analysis
    [chrstarts, genomelength,scafnames]= genomestats(ref_genome_folder);
    [locus_tag,chr_tag] = genomic_position_all(annotation_genes, genomelength, chrstarts) #loc_tag numbers all genes(CDS) across genome (intergenic/rna:0.5); chr_tag: numbers genes per chr (intergenic: #chr+0.5)
    #    locus_tag >> cds_number_by_position_tags
    #    chr_tag >> cds_number_by_position_indices
    # extract vector genelengths that follows order
    
    annotation_genes_sglDF = pd.concat(annotation_genes,sort=False)
    start_pos_gene = annotation_genes_sglDF['loc1'].values
    start_pos_gene = start_pos_gene[~np.isnan(start_pos_gene)] # remove nan, caused by empty chr
    end_pos_gene = annotation_genes_sglDF['loc2'].values
    end_pos_gene = end_pos_gene[~np.isnan(end_pos_gene)] # remove nan, caused by empty chr
    # gene lengths. incl rna genes (however when used those indexes not present in locus_tag (0.5))
    genelengths = (end_pos_gene - start_pos_gene)
    
    expectedNumberGenesMultipleMutations = np.zeros(numtrials) #initialize vector to store simulation results
    expectedNumberOfGenesWithNmutations = np.zeros((numtrials,max_muts_per_gene_to_track));
    # start sim
    for i in range(numtrials):
        # Pick random positions on the genome to mutate
        # Initially get 10x the number of positions you need --> later filter only for num_mutations_genic mutations
        randpos = np.random.randint(1, genomelength, 10*num_mutations_genome , dtype=int)
        # Does NOT assume any mutational spectrum!!!
        
        #TDL has scripts for doing the same thing for operon and pathway
        #levels, but the inputs they are generated on may not be relevant for
        #your genome
    
        # Find out in which genes these mutations occurred
        genenums = locus_tag[randpos];
        genenums_intragenic = genenums[genenums>0.5]; #remove intergenic/rna-gene mutations
        genenums_intragenic = genenums_intragenic[0:num_mutations_genic]; # only take as many as you need (generated extra because some removed due to not being on a gene)  
        # NOTE: may be overestimating the number of gene mutations in genic genes vs total num mutations.....
        expectedNumberIntragenicMutations=len(genenums_intragenic)
        # Get a list of the genes mutated, along with the number of times each
        [sim_unique_genes , sim_mutspergene] = np.unique(genenums_intragenic,return_counts=True)
        
        # This calculates the number of mutations per gene length for each gene
        # on which there were simulated mutations
        sim_mutspergenelength = sim_mutspergene/genelengths[sim_unique_genes.astype(int)-1]; #-1 bcs sim_unique_genes 1-based tagnumber > 1st gene is position 0 in genelengths
        sim_mutspergene_above_threshold=sim_mutspergene[(sim_mutspergenelength > mutation_density_threshold)]

        # The final step finds the number of genes that were mutated multiple 
        # times and above the threshold mutation density (per unit length). 
        # Result saved, indexed by trial number.
        expectedNumberGenesMultipleMutations[i] = sum( (sim_mutspergenelength >= mutation_density_threshold) & (sim_mutspergene >= mutation_number_threshold)); # Arolyn, 2018.12.14 > to >= to match input from function
    
        # The second piece of information this script returns is the number of
        # genes with >m mutations
        for j in range(max_muts_per_gene_to_track):
            expectedNumberOfGenesWithNmutations[i,j] = sum( sim_mutspergene_above_threshold >= j+1 ) # +1 due to 0-based
    return [expectedNumberGenesMultipleMutations, expectedNumberOfGenesWithNmutations ]

def codon_composition_table( allbases, allcodons ):
    ''' Build table containing base counts (col) for each codon (row) '''
    codoncompositiontable = np.zeros((len(allcodons),len(allbases)),dtype=int)
    for i,codon in enumerate(allcodons):
        for j,base in enumerate(allbases):
            codoncompositiontable[i,j] = codon.count(base)
    return codoncompositiontable

def codon_mutation_table(allmuts , allcodons , codon_all_dict):
    ''' Generates table of probabilities that a given mutation on a codon is nonsynonymous '''    
    table = np.zeros( (len(allcodons), len(allmuts) ) ); # Initialize table
    for i,codon in enumerate(allcodons):
        for j,mut in enumerate(allmuts):
            # Calculates the probability that a given mutation is nonsynonymous on a given codon and then updates the table
            table[i,j] = prob_nonsyn_codon_mutation( codon, mut , codon_all_dict);
    return table

def prob_nonsyn_codon_mutation(codon,mut,codon_all_dict):
    ''' Calculates the probability that a given mutation leads to a nonsynonymous change across a given codon '''
    aa0 = codon_all_dict[codon] # AA of codon
            
    # Find the positions on the codon at which mutation could occur
    possiblemuts=[i for (i, base) in enumerate(codon) if base == mut[0]]
    ctr = 0
    if len(possiblemuts) == 0: # if the mutation cannot occur on this codon
        probability = float('nan');
    else: # mut can occur at least once
        for pos in possiblemuts:
            newcodon = list(codon)
            newcodon[pos] = mut[1] # mutate codon position that carries mut[0]
            aa1 = codon_all_dict[ "".join(newcodon) ]
            if aa0 != aa1:
                ctr += 1
        probability = ctr/len(possiblemuts) # fraction of mutations that were nonsynonymous
    return probability
     
def codons_in_genome(annotation_genes,allcodons,genewise=False):
    ''' Get probability for each codon across all CDS annotated in the reference genome (annotation_genes), optional output as a per-gene dictionary of codon proportions (NOTE: corrected later if investigating subset of genes)'''
    # possibility to add a flag in order to restrict analysis to genomic region
    annotation_genes_sglDF = pd.concat(annotation_genes,sort=False)
    annotation_genes_sglDF_CDS = annotation_genes_sglDF.loc[ (annotation_genes_sglDF['type'] ==  'CDS') | (annotation_genes_sglDF['type'] ==  'gene') ]
    # Tally codon occurrences over all proteins in genome
    codonCounts = np.zeros( len(allcodons) ,dtype=int ); # for storing tally of codons
    if not genewise:
        for i, row in annotation_genes_sglDF_CDS.iterrows():
            seq = str(row['sequence'])
            startCodon = row['codon_start']
            startPos = startCodon * 3
            codons_gene = [seq[i:i+3] for i in range(startPos, len(seq), 3)]
            codons_gene = collections.Counter(codons_gene) # builds dictionary with all codons (key) in list and counts (value)
            for i,codon in enumerate(allcodons):
                if codon in codons_gene.keys():
                    codonCounts[i] += codons_gene[codon]
        return codonCounts/sum(codonCounts) # probabilities of codons sorted by allcodons order
    else:
        codon_counts_genewise={}
        for i, row in annotation_genes_sglDF_CDS.iterrows():
            codonCounts = np.zeros( len(allcodons) ,dtype=int )
            locustag = str(row['locustag'])
            seq = str(row['sequence'])
            startCodon = row['codon_start']
            startPos = startCodon * 3
            codons_gene = [seq[i:i+3] for i in range(startPos, len(seq), 3)]
            codons_gene = collections.Counter(codons_gene) # builds dictionary with all codons (key) in list and counts (value)
            for i,codon in enumerate(allcodons):
                if codon in codons_gene.keys():
                    codonCounts[i] += codons_gene[codon]
            codon_counts_genewise[locustag]=codonCounts/sum(codonCounts)
        return codon_counts_genewise # probabilities of codons per gene

def cacluate_expected_dn_ds_per_gene(params_dict,annotation_genes):
    ''' This script examines the probability of a nonsynonymous mutation on a 
    reference genome given some mutational spectrum. '''
    # HACKATHON
    # based on arolyn's m script: @Feb 2019
    ## Define DNA bases, possible mutations, and codons
    # All possible mutations to be considered:
    allbases = np.array(['A','T','G','C']) #fourDNAbases
    allmuts = np.array(['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG']); # 'AT' denotes A to T
    
    # All standard codons
    standard_codon_table = CodonTable.unambiguous_dna_by_id[1]
    allcodons = np.array([c for c in standard_codon_table.forward_table.keys()],dtype=object) # standard codon table, but miss stop codosn
    allcodons = np.append(allcodons , np.array([c for c in standard_codon_table.stop_codons],dtype=object) ) # 64 codons
    # build a combined dictionary (codon>>AA) containing regular and stop codon AA (stop:*)
    codon_all_dict = {}
    for c in allcodons:
        if c in standard_codon_table.forward_table.keys():
            codon_all_dict[c] = standard_codon_table.forward_table[c]
        else:
            codon_all_dict[c] = "*"

    # Generate table of codon composition by base
    codoncompositiontable = codon_composition_table( allbases, allcodons );
    # Rows (64) = all possible codons
    # Columns (4) = number of A's/T's/G's/C's in each codon
    
    ## Generate table of probabilities of nonsynonymous mutations
    # Rows (64) = all possible codons
    # Columns (12) = all possible mutations
    # Entries = probability of nonsynonymous mutation (given a mutation and a codon)
    # Note: If the given mutation cannot be applied to a given codon, the entry is nan. 

    codonnonsyntable = codon_mutation_table( allmuts, allcodons , codon_all_dict );

    ## Calculate mutational spectrum
    # Assuming all types of mutations are equally likely to occur:
    # Get mutational spectrum from experiments, but ignore if fwd or rev strand
    # allmuts = ['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG']
    #AT, TA  0
    #AC, TG  1
    #AG, TC  2
    #GC, CG  3
    #GT, CA  4
    #GA, CT  5
    
    if isinstance(params_dict['substitution_spectrum'], str):
        afile = open(params_dict['substitution_spectrum'], 'rb')
        mutationalspectrum = pickle.load(afile)
        afile.close()
    elif isinstance(params_dict['substitution_spectrum'], dict):
        mutationalspectrum = [params_dict['substitution_spectrum'][mut] for mut in allmuts] ## extract from dictionary in the order of interest the probabilities
    else:
        mutationalspectrum = [1/12] * 12 # uniform distribution. replace with time stratified observation based on all samples

    ## Calculate codon distribution in reference genome ordered as in allcodons
    codondistribution = codons_in_genome( annotation_genes , allcodons, True );

    # Calculate probability of nonsynonymous mutations per gene and output as dict for gene locustag --> prob    
    locustag_based_nonsyn_prob={}
    for genelocustag in codondistribution:
        probnonsyn = 0; # probability of nonsynonymous mutations over all possible mutations

        for i,mut in enumerate(allmuts): # loop through all possible mutations
            # Probability that this mutation occurs
            prob_base_base = mutationalspectrum[i]
            # Find the codons that can undergo this mutation:
            base = mut[0] # base that gets mutated; ex. A
            baseindex = np.where(allbases==base) # ex. A is indexed in position 1
            basecodonoccurrences = codoncompositiontable[:,baseindex].flatten()  # ex. how many A's in each codon
            # Indices of codons that have the relevant initial base:
            basecodonoccurrences_bool = basecodonoccurrences > 0; # ex. AAT has A's but GGC does not
            # Probability that this mutation occurs on a given relevant codon
            # Take into account base composition of codons
            basecountincodon = basecodonoccurrences[ basecodonoccurrences_bool ];
            # Take into account codon abundance on reference genome
            probcodonongenome = codondistribution[genelocustag][ basecodonoccurrences_bool ]
            # Combine these two probabilities
            probmutoncodon = basecountincodon*probcodonongenome
            # Renormalize (sum = 1 over all relevant codons)
            probmutoncodon = probmutoncodon/sum(probmutoncodon);

            # Probability that this mutation is nonsynonymous at each relevant codon
            thismutnonsynoncodon = codonnonsyntable[:,i];
            probmutnonsynoncodon = thismutnonsynoncodon[basecodonoccurrences_bool];
            # Overall probability that this mutation is nonsynonymous over all possible codons
            probmutnonsyn=prob_base_base*sum(probmutoncodon*probmutnonsynoncodon);

            # Add contribution of this mutation to the total probability of a nonsynonymous mutation
            probnonsyn += probmutnonsyn  
        locustag_based_nonsyn_prob[genelocustag]=probnonsyn
    return locustag_based_nonsyn_prob # Prob for N occuring 

def calculate_dn_ds_per_gene(params_dict,annotation_genes,annotation_mutation,calc_proportion_confint=False):
    """returns probability of every gene excess nonsyn and synonymous based on mutations provided in annotation mutation, optional output of 95 CI using wilson"""
    # HACKATHON
    dn_prob_per_gene=cacluate_expected_dn_ds_per_gene(params_dict,annotation_genes)
    output_probability_enrichment_nonsyn_or_syn={}
    for locustag in np.unique(annotation_mutation['locustag'][annotation_mutation['type'].isin(['N','S'])]):
        this_locustag_annotation_mutations=annotation_mutation[annotation_mutation['locustag'].isin([locustag])]
        cand_muts_N=0
        cand_muts_S=0
        for i,row in this_locustag_annotation_mutations.iterrows():
            ## ensures only count if in genic regions
            cand_muts_N += sum([muttype=='N' for muttype in row['type']])
            cand_muts_S += sum([muttype=='S' for muttype in row['type']])

        prob_obsN_excess = stats.binomtest(cand_muts_N, n=(cand_muts_N+cand_muts_S), p=dn_prob_per_gene[locustag], alternative='greater')
        prob_obsS_excess = stats.binomtest(cand_muts_S, n=(cand_muts_N+cand_muts_S), p=(1-dn_prob_per_gene[locustag]), alternative='greater')
        if calc_proportion_confint:
            N_prop_conf=proportion_confint(cand_muts_N,cand_muts_N+cand_muts_S,0.05,method='wilson')
            S_prop_conf=proportion_confint(cand_muts_S,cand_muts_N+cand_muts_S,0.05,method='wilson')
            output_probability_enrichment_nonsyn_or_syn[locustag]=(prob_obsN_excess,prob_obsS_excess,N_prop_conf,S_prop_conf)
        else: output_probability_enrichment_nonsyn_or_syn[locustag]=(prob_obsN_excess,prob_obsS_excess)
    return output_probability_enrichment_nonsyn_or_syn

def calculate_dn_ds_across_genes(params_dict,annotation_genes,annotation_mutation,calc_proportion_confint=False,precalculated_gene_dict=None):
    if precalculated_gene_dict:
        dn_prob_per_gene=precalculated_gene_dict
    else: dn_prob_per_gene=cacluate_expected_dn_ds_per_gene(params_dict,annotation_genes)
    locustags_investigated=[]
    annotation_genes_sglDF = pd.concat(annotation_genes,sort=False)
    cand_muts_N=0
    cand_muts_S=0
    for locustag in np.unique(annotation_mutation['locustag'][annotation_mutation['type'].isin(['N','S'])]):
        locustags_investigated.append(locustag)
        this_locustag_annotation_mutations=annotation_mutation[annotation_mutation['locustag'].isin([locustag])]
        for i,row in this_locustag_annotation_mutations.iterrows():
            ## ensures only count if in genic regions
            ## ensures only count if in genic regions
            cand_muts_N += sum([muttype=='N' for muttype in row['type']])
            cand_muts_S += sum([muttype=='S' for muttype in row['type']])

    # get proportion weighting each gene should get for contribution to expected nonsyn muts (weight by proportion of total gene lengths)
    seq_lengths=[]
    for locustag in locustags_investigated:
        seq_lengths.append(len(annotation_genes_sglDF[annotation_genes_sglDF['locustag'].isin([locustag])]['sequence']))
    weighting=np.array(seq_lengths)/np.sum(np.array(seq_lengths))
    overall_prob=0
    for index,locustag in enumerate(locustags_investigated):
        overall_prob+=dn_prob_per_gene[locustag]*weighting[index]
    sum_all_genic_mutations=(cand_muts_N+cand_muts_S)
    if sum_all_genic_mutations >0:
        prob_obsN_excess = stats.binomtest(cand_muts_N, n=sum_all_genic_mutations, p=overall_prob, alternative='greater')
        prob_obsS_excess = stats.binomtest(cand_muts_S, n=sum_all_genic_mutations, p=(1-overall_prob), alternative='greater')
        # calculate binom
        if calc_proportion_confint:
            N_prop_conf=proportion_confint(cand_muts_N,sum_all_genic_mutations,0.05,method='wilson')
            S_prop_conf=proportion_confint(cand_muts_S,sum_all_genic_mutations,0.05,method='wilson')
            return (cand_muts_N/sum_all_genic_mutations,cand_muts_S/sum_all_genic_mutations,N_prop_conf,S_prop_conf,overall_prob,sum_all_genic_mutations)
        else: return (prob_obsN_excess,prob_obsS_excess,overall_prob)
    else: 
        if calc_proportion_confint: return (np.NaN,np.NaN,np.NaN,np.NaN,overall_prob)
        else:
            return (np.NaN,np.NaN,overall_prob)


def mutation_probability(params_dict,annotation_genes,print_statements=True):
    ''' This script examines the probability of a nonsynonymous mutation on a 
        reference genome given some mutational spectrum. '''
    # based on arolyn's m script: @Feb 2019
    ## Define DNA bases, possible mutations, and codons
    # All possible mutations to be considered:
    allbases = np.array(['A','T','G','C']) #fourDNAbases
    allmuts = np.array(['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG']); # 'AT' denotes A to T
    
    # All standard codons
    standard_codon_table = CodonTable.unambiguous_dna_by_id[1]
    allcodons = np.array([c for c in standard_codon_table.forward_table.keys()],dtype=object) # standard codon table, but miss stop codosn
    allcodons = np.append(allcodons , np.array([c for c in standard_codon_table.stop_codons],dtype=object) ) # 64 codons
    # build a combined dictionary (codon>>AA) containing regular and stop codon AA (stop:*)
    codon_all_dict = {}
    for c in allcodons:
        if c in standard_codon_table.forward_table.keys():
            codon_all_dict[c] = standard_codon_table.forward_table[c]
        else:
            codon_all_dict[c] = "*"

    # Generate table of codon composition by base
    codoncompositiontable = codon_composition_table( allbases, allcodons );
    # Rows (64) = all possible codons
    # Columns (4) = number of A's/T's/G's/C's in each codon
    
    ## Generate table of probabilities of nonsynonymous mutations
    # Rows (64) = all possible codons
    # Columns (12) = all possible mutations
    # Entries = probability of nonsynonymous mutation (given a mutation and a codon)
    # Note: If the given mutation cannot be applied to a given codon, the entry is nan. 

    codonnonsyntable = codon_mutation_table( allmuts, allcodons , codon_all_dict );

    ## Calculate mutational spectrum
    # Assuming all types of mutations are equally likely to occur:
    # Get mutational spectrum from experiments, but ignore if fwd or rev strand
    # allmuts = ['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG']
    #AT, TA  0
    #AC, TG  1
    #AG, TC  2
    #GC, CG  3
    #GT, CA  4
    #GA, CT  5
    if isinstance(params_dict['substitution_spectrum'], str):
        afile = open(params_dict['substitution_spectrum'], 'rb')
        mutationalspectrum = pickle.load(afile)
        afile.close()
    elif isinstance(params_dict['substitution_spectrum'], dict):
        mutationalspectrum = [params_dict['substitution_spectrum'][mut] for mut in allmuts] ## extract from dictionary in the order of interest the probabilities
    else:
        mutationalspectrum = [1/12] * 12 # uniform distribution. replace with time stratified observation based on all samples

    ## Calculate codon distribution in reference genome ordered as in allcodons
    codondistribution = codons_in_genome( annotation_genes , allcodons, False );
    
    ## Calculate probability of nonsynonymous mutation
    # Takes into account: mutation spectrum, abundance of codons on genome,
    # abundance of bases in codons, probability of a given mutation being
    # nonsynonymous on a givne codon...
    probnonsyn = 0; # probability of nonsynonymous mutations over all possible mutations

    for i,mut in enumerate(allmuts): # loop through all possible mutations
        # Probability that this mutation occurs
        prob_base_base = mutationalspectrum[i]
        # Find the codons that can undergo this mutation:
        base = mut[0] # base that gets mutated; ex. A
        baseindex = np.where(allbases==base) # ex. A is indexed in position 1
        basecodonoccurrences = codoncompositiontable[:,baseindex].flatten()  # ex. how many A's in each codon
        # Indices of codons that have the relevant initial base:
        basecodonoccurrences_bool = basecodonoccurrences > 0; # ex. AAT has A's but GGC does not
        # Probability that this mutation occurs on a given relevant codon
        # Take into account base composition of codons
        basecountincodon = basecodonoccurrences[ basecodonoccurrences_bool ];
        # Take into account codon abundance on reference genome
        probcodonongenome = codondistribution[ basecodonoccurrences_bool ]
        # Combine these two probabilities
        probmutoncodon = basecountincodon*probcodonongenome
        # Renormalize (sum = 1 over all relevant codons)
        probmutoncodon = probmutoncodon/sum(probmutoncodon);

        # Probability that this mutation is nonsynonymous at each relevant codon
        thismutnonsynoncodon = codonnonsyntable[:,i];
        probmutnonsynoncodon = thismutnonsynoncodon[basecodonoccurrences_bool];
        # Overall probability that this mutation is nonsynonymous over all possible codons
        probmutnonsyn=prob_base_base*sum(probmutoncodon*probmutnonsynoncodon);

        # Add contribution of this mutation to the total probability of a nonsynonymous mutation
        probnonsyn += probmutnonsyn  
        
    if print_statements: print('Probability of nonsynonymous mutation across genome: ' + str(probnonsyn) )
    return probnonsyn # Prob for N occuring 

def parallel_evo_module(goodpos2use,contig_positions,annotation_mutations, annotation_genes, params_dict ,plot=True,ortholog=True,plot_title=None, count_multiply_mut_pos=True, dn_ds_only_cand_genes=True,return_prob_sim_cand_genes=False, testing=False, quiet=False, exclude_adj_muts = False):
    # goodpos2use = goodpos
    # contig_position = contig_positions
    # annotation_mutations = annotation_mutations_dedup
    # annotation_genes = annotation_genes
    # params_dict = parameters
    # plot=True
    # ortholog=True
    # plot_title=None
    # count_multiply_mut_pos=True
    # dn_ds_only_cand_genes=True
    # return_prob_sim_cand_genes=False
    # testing=False
    # quiet=False
    # exclude_adj_muts = False
    
    '''
    Changelog:
    06/2022 IL: added option to not parse orthologs and still output probabilities and plots
    06/2022 IL: added optional plotting output
    09/2022 IL: added optional parsing of 'output_name' to params dictionary to name plot output
    09/2022 IL: added optional plot_title parameter for having plot title
    09/2022 IL: added dynamic xtick plotting
    10/2022 IL: added tight layout for plot export
    11/2022 IL: added parsing of multiple mutation entry in annotation_mutations input (num_mutational_events)
    11/2022 IL: moved dN/dS calculation to helper function, added option to calc dN/dS on all genes or just cand genes (default)
    
    ## Module to calculate parallel evolution 
    ### Test excess of genes with multiple mutations (thresholds for candidates defined in parameters dict ), recommendation is to use 2 to find all parallely mutated genes above threshold
    ### Test excess of NonSyn (dNdS) in candidate genes compared to expectation given reference genome 
    ### NOTE: dNdS uses simple substitution model (equal P per mut) unless given in params_dict, which should be updated based on observations '''
    #HACKATHON -- decide what to do
    if not quiet: print('Parallel evolution inference.')
    # Find mutations that are adjacent to each other.
    # True for any mutation that is followed by an adjacent mutation, False if not
    # keep trailing bases of adjacent sets (2+)
    if len(contig_positions) == len(goodpos2use):
        chr_pos_gp = contig_positions
    else:
        chr_pos_gp = contig_positions[goodpos2use,]
    bool_adjacent_mut = np.full( chr_pos_gp.shape[0] , False, dtype=bool) #mutated_genes_pos = []; # keep track of the positions of each event in the table
    if exclude_adj_muts:
        for i in range(chr_pos_gp.shape[0]):
            chr = chr_pos_gp[i,0]
            pos = chr_pos_gp[i,1]
            if (i+1) <= (chr_pos_gp.shape[0]-1): # (chr_pos_gp.shape[0]-1) bcs shape[0] not 0-based
                if chr == chr_pos_gp[i+1,0] and (pos+1 == chr_pos_gp[i+1,1]):
                    bool_adjacent_mut[i] = True
        
    annotation_mutations['chr_locustag'] = annotation_mutations['chr'].apply(str)+'_'+annotation_mutations['locustag']
    
    # get info for candidate genes
    mutated_genes = np.zeros(0,dtype=int); # keep track of gene number
    mutated_genes_tally = np.zeros(0,dtype=int); # keep track of how many SNPs there are on this gene
    mutated_genes_lengths = np.zeros(0,dtype=int); # keep track of the length of this gene
    locustags_all = np.zeros(0,dtype=object) # record the locustag
    orthologtags_all = np.zeros(0,dtype=object) # record the orthologtags
    gene_tally_dict = {} # key=locustag; value=list(first entry tally, second is len, third is nonsyn, 4th is syn)
    for i, row in annotation_mutations.iterrows():
        if bool_adjacent_mut[i] == False: # ignore leading adjacent mutations
            gene_num_global = row['gene_num_global']
            chr_locustag = row['chr_locustag']
            if gene_num_global != 0.5 and gene_num_global != 0: # non-genic 0.5
                if 'num_mutational_events' in annotation_mutations and count_multiply_mut_pos: ## mutational events counted
                    if chr_locustag in gene_tally_dict:
                        gene_tally_dict[chr_locustag][0] += row['num_mutational_events']
                        gene_tally_dict[chr_locustag][2] += sum([nonsyn for nonsyn in row['NonSyn_rel_anc']])
                        gene_tally_dict[chr_locustag][3] += sum( np.invert([nonsyn for nonsyn in row['NonSyn_rel_anc']]) )
                    else:
                        mutated_genes_lengths =  (row['loc2']-row['loc1']+1) # +1 bcs loc1/2 1-based, thus underestimate L by 1bp
                        num_Nonsyns = sum([nonsyn for nonsyn in row['NonSyn_rel_anc']])
                        num_Syns = sum([nonsyn for nonsyn in row['NonSyn_rel_anc']])
                        gene_tally_dict[chr_locustag] = [row['num_mutational_events'],mutated_genes_lengths, num_Nonsyns,num_Syns]
                        locustags_all = np.append(locustags_all , chr_locustag)
                        if ortholog:
                            orthologtags_all = np.append(orthologtags_all , row['orthologtag'])
                        if testing:
                            print(mutated_genes_tally,mutated_genes)
                else:
                    if chr_locustag in gene_tally_dict:
                        gene_tally_dict[chr_locustag][0] += len(row['NonSyn_rel_anc'])
                        gene_tally_dict[chr_locustag][2] += sum([nonsyn for nonsyn in row['NonSyn_rel_anc']])
                        gene_tally_dict[chr_locustag][3] += sum( np.invert([nonsyn for nonsyn in row['NonSyn_rel_anc']]) )
                    else:
                        mutated_genes_lengths =  (row['loc2']-row['loc1']+1) # +1 bcs loc1/2 1-based, thus underestimate L by 1bp
                        num_Nonsyns = sum([nonsyn for nonsyn in row['NonSyn_rel_anc']])
                        num_Syns = sum([nonsyn for nonsyn in row['NonSyn_rel_anc']])
                        gene_tally_dict[chr_locustag] = [row['num_mutational_events'],mutated_genes_lengths, num_Nonsyns,num_Syns]
                        locustags_all = np.append(locustags_all , chr_locustag)
                        if ortholog:
                            orthologtags_all = np.append(orthologtags_all , row['orthologtag'])
                        if testing:
                            print(mutated_genes_tally,mutated_genes)
    mutated_genes_tally=np.array([tally_len[0] for tally_len in gene_tally_dict.values()])
    mutated_genes_lengths=np.array([tally_len[1] for tally_len in gene_tally_dict.values()])
    mutated_genes_nonsyn=np.array([tally_len[2] for tally_len in gene_tally_dict.values()])
    mutated_genes_syn=np.array([tally_len[3] for tally_len in gene_tally_dict.values()])
    mutated_genes=np.array([gene_locustag for gene_locustag in gene_tally_dict.keys()])
    mutated_genes_tally_perGeneLen = mutated_genes_tally/mutated_genes_lengths;
    ## get all gene lengths
    total_CDS_lengths=0
    for genomic_element in annotation_genes:
        for i,gene in genomic_element.iterrows():
            if not gene.empty: ## check if contig withou annotated gene
                if gene['type'] == 'CDS': ## get only cds as we count mutations only in CDS!
                    total_CDS_lengths+=gene['loc2']-gene['loc1']

    #% Define candidates for selection
    mutation_number_threshold = params_dict['Min_num_mutations_cand']; # minimum number of mutations per gene
    mutation_density_threshold = params_dict['Min_mutation_density_cand']; # minimum number of mutations per 1000 bp
    
    mutated_genes_of_interest = ( mutated_genes_tally >= mutation_number_threshold) & (mutated_genes_tally_perGeneLen >= mutation_density_threshold );
    num_mutated_genes_of_interest = sum(mutated_genes_of_interest);
    if not quiet:
        print('Number of genes with multiple mutations: ' + str(num_mutated_genes_of_interest))
    with open("parallel_evo_module_results" + params_dict['subjectID'] + ".txt", "w") as myfile:
        myfile.write('Number of genes with multiple mutations: ' + str(num_mutated_genes_of_interest) + "\n")
        myfile.write('Minimum number mutations required: ' + str(mutation_number_threshold) + "\n")
        myfile.write('Minimum SNP density for candidate genes: ' + str(mutation_density_threshold) + "\n")        
    
    # break if no candidate genes present
    if num_mutated_genes_of_interest == 0:
        if not dn_ds_only_cand_genes:
            if not quiet:
                print('Calculating dN/dS on all genes, mutations')
                calculate_genome_wise_dn_ds(params_dict,annotation_genes,annotation_mutations)
            else:
                calculate_genome_wise_dn_ds(params_dict,annotation_genes,annotation_mutations, print_output=False)
        if not quiet:
            print('No genes with multiple mutation found! >> skip adaptive evolution analysis')
        return [np.array([]),pd.DataFrame(),np.NaN] # return empty array and dataframe
    # get annotation_mutation for SNPs with signature of parallel evolution
    annotation_mutation_paraSignal = annotation_mutations.loc[ annotation_mutations['chr_locustag'].isin( mutated_genes[mutated_genes_of_interest] ) ]

    # =============================================================================
    # Test candidate genes for parallel mutation enrichment
    # =============================================================================
    # NOTE: Simulation based on random genic SNPs in genome and does not use a specific substitution model
    number_of_mutations_on_genome = len(goodpos2use)
    if 'num_mutational_events' in annotation_mutations.columns:
        number_of_mutations_on_genome = sum(annotation_mutations['num_mutational_events'])
    else:
        number_of_mutations_on_genome = sum(annotation_mutations['nts'].str.len()-1)
    number_of_mutations_genic = int(sum(mutated_genes_tally)); # only genic/CDS!

    [expectedNumberGenesMultipleMutations, expectedNumberOfGenesWithNmutations] = parallel_evolution_counting_and_simulation(number_of_mutations_on_genome,number_of_mutations_genic,params_dict['Min_num_mutations_cand'],params_dict['Min_mutation_density_cand'],params_dict['NumTrialsSim'],params_dict['max_muts_per_gene_to_track'],chr_pos_gp , params_dict['ref_genome_folder'],annotation_genes)
    # expectedNumberGenesMultipleMutations: num candidate genes per sim for parallel evolution based on params_dict['Min_num_mutations_cand'],params_dict['Min_mutation_density_cand']
    # expectedNumberOfGenesWithNmutations: row==sim; col==genes with mutations. col[0] == 1 mutation!, col[1]==2mut...
    
    simProbForObsCand = 1-sum(expectedNumberGenesMultipleMutations < num_mutated_genes_of_interest)/len(expectedNumberGenesMultipleMutations)
    if not quiet:
        print('Simulation-based probability to observe ' + str(num_mutated_genes_of_interest) + ' canddiate genes for parallel evolution is: '+str(simProbForObsCand))
        print('Simulated a mean number of candidate genes is ' + str(np.mean(expectedNumberGenesMultipleMutations)))
    with open("parallel_evo_module_results" + params_dict['subjectID'] + ".txt", "a") as myfile:
        myfile.write('Simulation-based probability to observe ' + str(num_mutated_genes_of_interest) + ' canddiate genes for parallel evolution is: '+str(simProbForObsCand) + "\n")
        myfile.write('Simulated a mean number of candidate genes is ' + str(np.mean(expectedNumberGenesMultipleMutations)) + "\n")

    # calc prob to observe candidate mut counts
    locustags_cand = locustags_all[mutated_genes_of_interest]
    gene_lengths_cand = mutated_genes_lengths[mutated_genes_of_interest].astype(int)
    if ortholog: orthologtags_cand = orthologtags_all[mutated_genes_of_interest]
    mut_cand_tally = mutated_genes_tally[mutated_genes_of_interest]
    prob_cand_nummut = np.ones(mut_cand_tally.shape)
    prob_cand_nummut_poisson = np.ones(mut_cand_tally.shape)
    prob_cand_nummut_binom = np.ones(mut_cand_tally.shape)
    syn_cand=mutated_genes_syn[mutated_genes_of_interest].astype(int)
    nonsyn_cand=mutated_genes_nonsyn[mutated_genes_of_interest].astype(int)
    for i,nummut in enumerate(mut_cand_tally):
        # original implemntation
        prob_cand_nummut[i] = 1- (np.sum(expectedNumberOfGenesWithNmutations[:,0:int(nummut-1)])/np.sum(expectedNumberOfGenesWithNmutations))
        #poisson
        p=number_of_mutations_genic/total_CDS_lengths
        prob_cand_nummut_poisson[i] = stats.poisson.sf(nummut,mutated_genes_lengths[mutated_genes_of_interest][i]*p)
        # binomial
        prob_cand_nummut_binom[i] = stats.binom.sf(nummut,mutated_genes_lengths[mutated_genes_of_interest][i],p)
    prob_cand_nummut_poisson=prob_cand_nummut_poisson.astype('float64')
    prob_cand_nummut_binom=prob_cand_nummut_binom.astype('float64')
    if ortholog:
        res_cand_nummut = np.vstack((orthologtags_cand,locustags_cand,mut_cand_tally,nonsyn_cand,syn_cand,gene_lengths_cand,prob_cand_nummut, prob_cand_nummut_poisson, prob_cand_nummut_binom)).T
        #res_cand_nummut = np.vstack((orthologtags_cand,locustags_cand,mut_cand_tally,prob_cand_nummut)).T
    
        if not quiet: 
            print('Orthologtag Locustag NumMuts NumMutsNonsyn NumMutsSyn GeneLen PercentileSim ProbPoisson ProbBinom')
            print(res_cand_nummut)    
        with open("parallel_evo_module_results" + params_dict['subjectID'] + ".txt", "a") as myfile:
            myfile.write('Orthologtag Locustag NumMuts NumMutsNonsyn NumMutsSyn GeneLen PercentileSim ProbPoisson ProbBinom'  + "\n")
            np.savetxt(myfile,  res_cand_nummut,fmt="%s")
    else:
        res_cand_nummut = np.vstack((locustags_cand,mut_cand_tally,nonsyn_cand,syn_cand,gene_lengths_cand,prob_cand_nummut,prob_cand_nummut_poisson, prob_cand_nummut_binom)).T
        if not quiet:
            print('Locustag NumMuts NumMutsNonsyn NumMutsSyn GeneLen PercentileSim ProbPoisson ProbBinom')
            print(res_cand_nummut)
        with open("parallel_evo_module_results" + params_dict['subjectID'] + ".txt", "a") as myfile:
            myfile.write('Locustag NumMuts NumMutsNonsyn NumMutsSyn GeneLen PercentileSim ProbPoisson ProbBinom'  + "\n")
            np.savetxt(myfile,  res_cand_nummut,fmt="%s")
    # =============================================================================
    # dN/dS calculation
    # =============================================================================
    # calc if we observe more NS than expected (simulated),
    
    if dn_ds_only_cand_genes:
        if not quiet:
            print('Calculating dN/dS on only candidate genes')
            calculate_genome_wise_dn_ds(params_dict,annotation_genes,annotation_mutations, mutated_genes, mutated_genes_of_interest)
        else:
            calculate_genome_wise_dn_ds(params_dict,annotation_genes,annotation_mutations, mutated_genes, mutated_genes_of_interest,print_output=False)
    else:
        if not quiet:
            print('Calculating dN/dS on all genes, mutations')
            calculate_genome_wise_dn_ds(params_dict,annotation_genes,annotation_mutations)
        else:
            calculate_genome_wise_dn_ds(params_dict,annotation_genes,annotation_mutations,print_output=False)
    # =============================================================================
    #     Plot
    # =============================================================================
    if plot:
        plt.rcParams.update({'font.size': 14}) # all label size
        f = plt.figure()
        ax=f.add_subplot(111)
    
        # histogram depicting simulated distribution of expected candidate genes, filter criteria for candidates, P for observation cand count
        mybin = np.arange(max(np.unique(expectedNumberGenesMultipleMutations))+2)-0.5 # +2 in order to now loose last bin due to arange; -0.5 needed for bars align at center with x axis ticks
        plt.hist(expectedNumberGenesMultipleMutations,bins=mybin,rwidth=0.8,color='#607c8e', edgecolor = 'k', lw = 0.4)
        max_on_x_axis=max(len(mut_cand_tally),max(np.unique(expectedNumberGenesMultipleMutations)))
        if max_on_x_axis>800:
            step=100
        elif max_on_x_axis>500:
            step=50
        elif max_on_x_axis>100:
            step=20
        elif max_on_x_axis>50:
            step=10
        elif max_on_x_axis>20:
            step=5
        else:
            step=1
        plt.xlim(0,max(len(mut_cand_tally),max(np.unique(expectedNumberGenesMultipleMutations)))+1)
        plt.xticks(np.arange(0,max(len(mut_cand_tally),max(np.unique(expectedNumberGenesMultipleMutations)))+step,step))
        plt.ylabel('Simulated counts, N='+str(params_dict['NumTrialsSim']))
        plt.xlabel('Number of genes with ' + str(params_dict['Min_num_mutations_cand']) + ' or more mutations')
        if plot_title:
            plt.title(plot_title)
        plt.axvline(x=len(mut_cand_tally),color='violet') 
        plt.text(0.98, 0.95, "min #mut:"+str(params_dict['Min_num_mutations_cand'])+"; min density:"+str(params_dict['Min_mutation_density_cand']), fontsize=12,horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
        plt.text(0.98, 0.88, "P("+ str(len(mut_cand_tally)) + ") = " + str(np.around(simProbForObsCand,3)), fontsize=12,horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
        #save pdf
        if 'analsysis_params_output_name_folder' in params_dict:
            output_folder=params_dict['analsysis_params_output_name_folder']+'/'
        else: output_folder=''
        if 'output_name' in params_dict:
            subprocess.run([f"mkdir -p pdf/adaptive_evo/{output_folder} "],shell=True)
            f.savefig(f'pdf/adaptive_evo/{output_folder}' + params_dict['output_name'] + "_" + params_dict['subjectID'] + ".pdf",bbox_inches="tight")
            plt.close()
            if not quiet:
                print(f'Plotted: pdf/adaptive_evo/{output_folder}' + params_dict['output_name'] + "_" + params_dict['subjectID'] + ".pdf")
        else:
            subprocess.run([f"mkdir -p pdf/adaptive_evo/{output_folder} "],shell=True)
            f.savefig(f'pdf/adaptive_evo/{output_folder}' + params_dict['timestamp'] + "_" + params_dict['subjectID'] + ".pdf",bbox_inches="tight")
            plt.close()
            if not quiet: print(f'Plotted: pdf/adaptive_evo/{output_folder}' + params_dict['timestamp'] + "_" + params_dict['subjectID'] + ".pdf")
    if return_prob_sim_cand_genes:
        return [res_cand_nummut,annotation_mutation_paraSignal,simProbForObsCand]
    return [res_cand_nummut,annotation_mutation_paraSignal]

def calculate_genome_wise_dn_ds(params_dict,annotation_genes,annotation_mutations, mutated_genes=[], genes_of_interest=[], save_results=True, print_output=True,return_values=False, calc_proportion_confint=False, quiet = False):
    """
    Function to calculate dN/dS ratio and compare to theoretical expectations at a genome wide level
    Set of total mutations can be subset to specific genes providing a set of chr_locustags for annotation_mutations
    
    was previously implemented in parallel evo module, which now calls this as a helper function
    IL 14.02.2023 added optional outputting of probs and CI, rather than just printing/saving results
    """
    if print_output: print('Investigate dN/dS')

    ## Gene numbers of genes of interest, if provided use subset, if not use all mutations (genic filtered later during counting)
    if len(genes_of_interest) > 0 and len(mutated_genes) > 0:
        cand_genes_chr_locustags = mutated_genes[genes_of_interest]; # these_gene_nums
    else:
        cand_genes_chr_locustags = annotation_mutations['chr_locustag'].unique()

    # Observed dNdS
    cand_muts_N = 0;
    cand_muts_S = 0;
    for chr_locustag in cand_genes_chr_locustags:
        cand_mut_anno = annotation_mutations.loc[ annotation_mutations['chr_locustag'] == chr_locustag ]
        print(cand_mut_anno[['NonSyn', 'type']])
        for i,row in cand_mut_anno.iterrows():
            ## ensures only count if in genic regions
            cand_muts_N += sum([muttype=='N' for muttype in row['type']])
            cand_muts_S += sum([muttype=='S' for muttype in row['type']])

    if print_output:
        if (cand_muts_N+cand_muts_S) == 0: # necessary in case no S/N mutations but all P/U
            print("No dN/dS calculated >> no S/N mutations.")
        else:
            probN_obs = cand_muts_N/(cand_muts_N+cand_muts_S)        
            # Simulate genome specific expected probability for NonSyn to occur
            # simulate only if value not yet calculated. > stored in regenome/probNsim.pk     
            # probN_sim = 0.6829; 
            probN_sim = mutation_probability(params_dict,annotation_genes,print_output)
            # NOTE: mutation_probability uses a pre-calculted substitution spectrum (mutationalspectrum) when specified in params or a uniform mutation spectrum.
            
            # Binomial test
            prob_obsN_excess = stats.binomtest(cand_muts_N, n=(cand_muts_N+cand_muts_S), p=probN_sim, alternative='greater')
            prob_obsS_excess = stats.binomtest(cand_muts_S, n=(cand_muts_N+cand_muts_S), p=(1-probN_sim), alternative='greater')
            print('Observed P(N): ' + str(probN_obs) + '; N: ' + str(cand_muts_N) + ', S: ' + str(cand_muts_S))
            print('Expected P(N): ' + str(probN_sim) )
            print('Probability of enrichment in observed count of N: ' + str(prob_obsN_excess) )
            print('Probability of enrichment in observed count of S: ' + str(prob_obsS_excess) )
            if save_results:
                with open("dn_ds_results_" + params_dict['subjectID'] + ".txt", "a") as myfile:
                    myfile.write('Observed P(N): ' + str(probN_obs) + '; N: ' + str(cand_muts_N) + ', S: ' + str(cand_muts_S) + "\n")
                    myfile.write('Expected P(N): ' + str(probN_sim) + "\n")
                    myfile.write('Probability of enrichment in observed count of N: ' + str(prob_obsN_excess) + "\n")
                    myfile.write('Probability of enrichment in observed count of S: ' + str(prob_obsS_excess) + "\n")
                print('Saved results in file ' + "dn_ds_results_" + params_dict['subjectID'] + ".txt")
    else:
        if (cand_muts_N+cand_muts_S) == 0:
            return (np.NaN,np.NaN)
        probN_obs = cand_muts_N/(cand_muts_N+cand_muts_S)  
        probN_sim = mutation_probability(params_dict,annotation_genes,print_output)
        # NOTE: mutation_probability uses a pre-calculted substitution spectrum (mutationalspectrum) when specified in params or a uniform mutation spectrum.
        # Binomial test
        prob_obsN_excess = stats.binomtest(cand_muts_N, n=(cand_muts_N+cand_muts_S), p=probN_sim, alternative='greater')
        prob_obsS_excess = stats.binomtest(cand_muts_S, n=(cand_muts_N+cand_muts_S), p=(1-probN_sim), alternative='greater')
        if save_results:
            subprocess.run([f"mkdir -p dn_ds_results"], shell=True)
            with open("dn_ds_results/dn_ds_results" + params_dict['subjectID'] + ".txt", "a") as myfile:
                myfile.write('Observed P(N): ' + str(probN_obs) + '; N: ' + str(cand_muts_N) + ', S: ' + str(cand_muts_S) + "\n")
                myfile.write('Expected P(N): ' + str(probN_sim) + "\n")
                myfile.write('Probability of enrichment in observed count of N: ' + str(prob_obsN_excess) + "\n")
                myfile.write('Probability of enrichment in observed count of S: ' + str(prob_obsS_excess) + "\n")
        if return_values==False:
            if calc_proportion_confint:
                N_prop_conf=proportion_confint(cand_muts_N,cand_muts_N+cand_muts_S,0.05,method='wilson')
                S_prop_conf=proportion_confint(cand_muts_S,cand_muts_N+cand_muts_S,0.05,method='wilson')
                return (prob_obsN_excess,prob_obsS_excess,N_prop_conf,S_prop_conf)
            return (prob_obsN_excess,prob_obsS_excess)
        else:
            return (prob_obsN_excess,prob_obsS_excess,probN_obs)



def get_index_visit_sampleNames(spl_names_long,treesampleNamesLong=True):
    # get idx V1,2,3,4,5 based on treesampleNamesLong
    # later extended funtion to also use visits variable instead of treesampleNamesLong
    list_visit_idx = []
    if treesampleNamesLong:
        for v in ['TP' + str(i) for i in range(1, 36)]:
            match_visit = np.array([i-1 for i,x in enumerate(spl_names_long) if re.search('_'+v+'_',x)])
            list_visit_idx.append(match_visit)
    else:
        for v in [1,2,3,4,5]:
            match_visit = np.array([i for i,x in enumerate(spl_names_long) if x==v])
            list_visit_idx.append(match_visit)        
    return list_visit_idx

def calc_regression(mean_count_per_visit,subject_fld_label,num_samples,num_visits):
    # calculate regression 
    # return dictionary with results and some additional metainfo
    slope, intercept, r_value, p_value, std_err = stats.linregress(mean_count_per_visit[:,0] , mean_count_per_visit[:,1] )
    regress_dict = {}
    regress_dict['slope'] = slope
    regress_dict['intercept'] = intercept
    regress_dict['r_sq'] = r_value**2
    regress_dict['p_val'] = p_value
    regress_dict['std_err'] = std_err
    regress_dict['subject'] = subject_fld_label
    regress_dict['num_samples'] = num_samples
    regress_dict['num_visits'] = num_visits
    return regress_dict

def plot_molecular_clock(input_data,hdr,mean_count_per_visit,regress_dict,basescale,subject_identifier, genome_size = 1e6, jitter_x = 0, jitter_y = 0, title_suffix = '', plot_confidence_interval=False, plot_predictor_interval=False):
    ''' plot molecular clock for genome-size corrected inference. rescale to genome-wide clock values using eg. 'basescale'=100000 '''
    subprocess.run(["mkdir -p pdf/molclock/ "],shell=True) # build output folder
    # turn to df for seaborn
    # basescale 
    data_df = pd.DataFrame(data=input_data,    # values
                 index=np.arange(input_data.shape[0]),    # 1st column as index
                 columns=hdr)  # 1st row as the column names        
    data_df['is_even'] = (data_df['colonization_days'] % 2) == 0 # add T/F for coloring light/darkgrey
    data_df_scaled = data_df.copy()
    data_df_scaled['rate'] = data_df_scaled['rate']*basescale
    xpos = np.unique(data_df['colonization_days']) # time in months for adding line
    slope = regress_dict['slope']
    std_err = regress_dict['std_err']
    r_value_sq = regress_dict['r_sq']
    p_value = regress_dict['p_val']
    subject_id = regress_dict['subject']
    intercept = regress_dict['intercept']
    mut_rate_per_day = slope*basescale
    mean_count_per_visit_d = mean_count_per_visit.keys()
    mean_count_per_visit_vals = [cnt * basescale for cnt in mean_count_per_visit.values()]
    basescale_mb = basescale/1e6
    if int(basescale_mb) == basescale_mb:
        basescale_mb = int(basescale_mb)
    dof = len(data_df)-2 ## degrees of freedome for linear regression = len(df)-2
    alpha = 0.05
    ## calculate 95% CI values for mutation rate
    ts = abs(stats.t.ppf( 1-(alpha / 2), dof)) ## (confidence interval and the degrees of freedom minus num of parameters to estimate)
    margin_of_error = std_err*basescale*365 * ts ## moe = std_err * t-score
    ci_lower = mut_rate_per_day*365 - margin_of_error
    ci_upper = mut_rate_per_day*365 + margin_of_error
    
    annotation_text = f'Mutation rate (mut/y/{basescale_mb}Mb) = {np.round(365*mut_rate_per_day ,3)} [{np.round(ci_lower, 3)}, {np.round(ci_upper, 3)}]\nr^2 = {np.round(r_value_sq ,3)}\np-value = {p_value:.2e}'
    useylim = max(data_df['rate']*basescale)*1.2
    
    ## to include confidence intervals
    if plot_confidence_interval:
        x_intercept = -intercept / slope
        time_since_xintercept = xpos-x_intercept
        
        ## get variance to calculate shape and scale (https://library.virginia.edu/data/articles/getting-started-with-gamma-regression)
        y_pred = intercept*basescale + mut_rate_per_day*data_df['colonization_days']
        residuals = data_df_scaled['rate'] - y_pred
        variance = np.sum(residuals**2) / dof

        ## estimate scale and shape ot gamma distribution
        scale = variance / (time_since_xintercept*mut_rate_per_day)
        scale = np.clip(scale, a_min = 1e-32, a_max = None) ## clip to avoid negative values
        shape = (time_since_xintercept*mut_rate_per_day)**2 / variance
        shape = np.clip(shape, a_min = 1e-32, a_max = None) ## clip to avoid negative values
        
        ## get 95% intervals
        lower_gam_ci = stats.gamma.ppf(alpha/2, a=shape, scale=scale)
        upper_gam_ci = stats.gamma.ppf(1-(alpha/2), a=shape, scale=scale)
        useylim = max(useylim, max(upper_gam_ci))

    ## include predictor intervals given the assumption of a strict clock model (fixing scale to 1 (->variance == mean) in respect to its genome size (to account for basescale!))
    if plot_predictor_interval:
        genome_size_mb = genome_size/1e6
        x_intercept = -intercept / slope
        time_since_xintercept = xpos-x_intercept

        ## estimate scale and shape of gamma distribution for a strict clock model
        scale = 1/genome_size_mb
        scale = np.clip(scale, a_min = 1e-32, a_max = None) ## clip to avoid negative values
        shape = time_since_xintercept*mut_rate_per_day*genome_size_mb
        shape = np.clip(shape, a_min = 1e-32, a_max = None) ## clip to avoid negative values

        ## get 95% intervals
        lower_gam_pi = stats.gamma.ppf(alpha/2, a=shape, scale=scale) 
        upper_gam_pi = stats.gamma.ppf(1-(alpha/2), a=shape, scale=scale) 
        useylim = max(useylim, max(upper_gam_pi))

    ## plot
    if useylim != useylim:
        useylim = basescale
    jitter_x_array = np.random.uniform(low=0-(jitter_x/2), high=0+(jitter_x/2), size=(len(data_df),)) ## generate random array for jitter on x axis
    jitter_y_array = np.random.uniform(low=0-(jitter_y/2), high=0+(jitter_y/2), size=(len(data_df),)) ## generate random array for jitter on x axis

    fig, ax = plt.subplots(figsize=(4.125,4))
    ax.scatter(data_df['colonization_days'] - jitter_x_array,
                (data_df['rate']*basescale) - jitter_y_array,
                facecolors='none', edgecolors='darkgrey',alpha=1, s = 40, lw = 0.5)
                #facecolors='none', c='darkgrey',alpha=0.8, s = 50, marker = 'x')
    # sns.regplot(data=data_df_scaled, x='colonization_days', y='rate', ci=95, scatter=False, color='k', ax=ax)
    ax.plot(xpos,(intercept*basescale) + slope*xpos*basescale,'k-', lw = 2)
    if plot_confidence_interval:
        ## plots 
        plt.fill_between(xpos, lower_gam_ci, upper_gam_ci,
                         color='darkgrey', alpha=0.2, label='95% CI')
    if plot_predictor_interval:
        ## plots expected stochastic variability in mutation accumulation over time
        plt.fill_between(xpos, lower_gam_pi, upper_gam_pi,
                         color='darkgrey', alpha=0.2, label='95% PI')
    ax.scatter(mean_count_per_visit_d, mean_count_per_visit_vals, marker="x", s=150, c='k', lw = 2.5)
    ax.text(0.04, 0.96, annotation_text, transform=ax.transAxes, 
            ha='left', va = 'top', fontsize=10)
    ax.set_ylim((0, max(ax.get_ylim()[-1], useylim)))
    ax.set_xlabel("Time (days)",fontsize=14)
    ax.set_title(subject_identifier, fontsize=12)
    ax.set_ylabel("Observed mutations per "+str(basescale_mb)+'Mb',fontsize=14)    
    plt.tight_layout()
    fig.savefig(f'pdf/molclock/{subject_id}_mol_clock_assembleBased{title_suffix}.pdf')
    plt.close()
    print('Molecular clock analysis: slope: ' + str(np.round(365*mut_rate_per_day ,3)) + ", r^2: " + str( np.round(r_value_sq ,3) ) )
    print("Done. pdf/molclock/" + subject_id + "_mol_clock_assembleBased{title_suffix}.pdf")


def infer_non_snp_events(mutantAF,p,distance_for_nonsnp=500,fraction_covariance_nonsnp=0.98, print_outputs=True):
    # NOTE: resolve which function to keep --> hackathon, see findrecominantSNPs
    ''' Mutations that  covary in close physical proximity (unexpected by chance) are likely non-independent mutations due to HGT event '''
    # return index of p with evidence for recombination (nonsnp) , boolean of length p with nonsnp == True
    p_nonsnps_events = np.array([])
    cov_m = np.cov(mutantAF)
    dist_p = p[1:]-p[:-1] # list of bp differences between elements of p (len(p)-1)
    proxy_pi = dist_p < distance_for_nonsnp
    proxy_pi = np.argwhere(proxy_pi) # index p and other p in close proximity
    for i,p_i in enumerate(proxy_pi):
        roi_min , roi_max = p[p_i] - distance_for_nonsnp , p[p_i] + distance_for_nonsnp # define region of interest with mutations, where we know at least 1 other mutation is in close proximity
        p_bool_roi = np.all( (p>roi_min , p<roi_max ),axis=0) # bool of p in proximity
        cov_p_i = cov_m[p_i,:][0,:] # cov of focal SNP p_i with all other p
        cov_p_bool = cov_p_i > (np.max(cov_p_i) * fraction_covariance_nonsnp) # bool of all SNPs that have x*max(cov) with SNP. Note: The focal SNP has ~highest cov with self
        if np.sum(cov_p_bool[p_bool_roi]) > 1: # if high cov of focal SNP with at least one more SNP in close proximity
            p_nonsnps_events = np.append( p_nonsnps_events, p[p_bool_roi][cov_p_bool[p_bool_roi]])
    p_nonsnps_events = np.unique(p_nonsnps_events)
    p_idx_nonsnp = np.array([np.where(p==i)[0][0] for i in p_nonsnps_events]) # get index of for nonsnps
    p_bool_nonsnp = np.zeros(p.shape) # build bool array of length p with nonsnps == True
    if p_idx_nonsnp.size > 0:
        p_bool_nonsnp[p_idx_nonsnp] = 1
    p_bool_nonsnp = p_bool_nonsnp.astype(bool)
    if print_outputs:
        print(p_idx_nonsnp.size,'/',p.size,'elements of candpos evidence for recombination event.')
    return [p_idx_nonsnp,p_bool_nonsnp]


"""
def findrecombinantSNPs(p, mutantAF, distance_for_nonsnp, corr_threshold_recombination, ncpu = 1):
    ''' Mutations that  covary in close physical proximity (unexpected by chance) are likely non-independent mutations due to HGT event '''
    print(f'{len(p)} sites will be evaluated for indications of recombination events.')

    #look for recombination regions
    nonsnp_idx = np.zeros(0,dtype='int')

    for i in range(len(p)):
        if i%2000 == 0: ## print pattern to stdout to show where in evaluation process runs 
            if i%10000 == 0:
                print('|', end = '')
            else:
                print('.', end = '')
        gp = p[i]
        #find nearby snps
        if gp > distance_for_nonsnp:
            region = np.array(np.where((p > gp - distance_for_nonsnp) & (p < gp + distance_for_nonsnp)) ).flatten()
            if len(region)>1: 
                r = mutantAF[region,:]
                corrmatrix = np.corrcoef(r)
                [a,b]=np.where(corrmatrix > corr_threshold_recombination)
                nonsnp_idx=np.concatenate((nonsnp_idx,region[a[np.where(a!=b)]]))

    nonsnp_idx=np.unique(nonsnp_idx)
    ## print number of sites which have indication for recombination
    print('\n' + str(nonsnp_idx.shape[0]) + ' of a total ' + str(p.shape[0]) + ' ('  + str(nonsnp_idx.shape[0]/p.shape[0]*100) + '%) positions in goodpos were found to be recombinant.')
    nonsnp_bool = np.zeros(p.shape) # build bool array of length p with nonsnps == True
    if nonsnp_idx.size > 0:
        nonsnp_bool[nonsnp_idx] = 1
    nonsnp_bool = nonsnp_bool.astype(bool)

    return [nonsnp_idx, nonsnp_bool]
"""
def expand_arr_by_nt(snv_p, mutantAF, mant, nts_idx):
    """Expand mutantAF and snv_p in a way that for every unique call at the given site, a separate row in mutantAF is given"""
    expanded_mutantAF = []
    expanded_snv_p = []
    p_mapper = []
    pos_idx = []

    ## duplicate rows for every nt idx
    for single_nt in nts_idx:
        # mask for positions matching current single_nt and keep only those
        mask = (mant == single_nt)  
        new_row = np.where(mask, mutantAF, 0) ## matrix either mutantAF or filled bt 0 for unmasked (different) nts
        expanded_mutantAF.append(new_row)
        expanded_snv_p.append(snv_p.copy()) # duplicate snv_p
        p_mapper.append([single_nt]*np.shape(new_row)[0]) ## generate a mapper to map positions back in the end
        pos_idx.append(np.arange(len(snv_p))) ## generate a mapper for every position idx to map back in the end

    # Convert to numpy arrays
    expanded_mutantAF = np.vstack(expanded_mutantAF)
    expanded_snv_p = np.concatenate(expanded_snv_p)
    p_mapper = np.concatenate(p_mapper)
    pos_idx = np.concatenate(pos_idx)

    ## remove uninofrmative rows --> all have no mutantAF as this allele does not appear at given position
    mask_uninfo_sites = np.any(expanded_mutantAF != 0, axis = 1)
    expanded_mutantAF = expanded_mutantAF[mask_uninfo_sites, :]
    expanded_snv_p = expanded_snv_p[mask_uninfo_sites]
    p_mapper = p_mapper[mask_uninfo_sites]
    pos_idx = pos_idx[mask_uninfo_sites]
    return [expanded_mutantAF, expanded_snv_p, p_mapper, pos_idx]

def process_recombinant_per_pos(i, snv_p, mutantAF, distance_for_nonsnp, corr_threshold_recombination):
    if i%2000 == 0: ## print pattern to stdout to show where in evaluation process runs 
        if i%10000 == 0:
            print('|', end = '')
        else:
            print('.', end = '')
    gp = snv_p[i]
    #find nearby snps
    if gp > distance_for_nonsnp:
        region = np.array(np.where((snv_p > gp - distance_for_nonsnp) & (snv_p < gp + distance_for_nonsnp)) ).flatten()
        if len(region)>1: 
            r = mutantAF[region,:]
            corrmatrix = np.corrcoef(r)
            [a,b]=np.where(corrmatrix > corr_threshold_recombination)
            nonsnp_idx=region[a[np.where(a!=b)]]
            return nonsnp_idx
        else:
            return np.array([]) ## return empty array for concatenation
    else:
        return np.array([]) ## return empty array for concatenation

def findrecombinantSNPs_new(snv_p, mutantAF, distance_for_nonsnp, corr_threshold_recombination, eval_recom_p_allele = False, nts_idx = np.arange(0, 5), mant = np.array([]), ncpu = 1):
    ''' Mutations that covary in close physical proximity (unexpected by chance) are likely non-independent mutations due to HGT event '''
    # up to 4 cores are recommended, but above that the I/O gets the restricting factor and can make the improvement worse!
    # NOTE: if using eval_recom_p_allele one should implement unfiltered (!) maNT instead of calls, otherwise it can happen 
    # that some recombinant sites are masked and separated into two groups (initial, unfiltered allele & filtered, unsupported (=N) allele)
    num_p = len(snv_p)
    print(f'{num_p} sites will be evaluated for indications of recombination events.')

    #look for recombination regions
    if eval_recom_p_allele:
        [expanded_mutantAF, expanded_snv_p, p_mapper, pos_idx] = expand_arr_by_nt(snv_p, mutantAF, mant, nts_idx)
        ## identify recombination sites which correlate with other nts
        if ncpu == 1:
            nonsnp_idx = []
            for i in range(len(expanded_snv_p)):
                nonsnp_idx.append(process_recombinant_per_pos(i, expanded_snv_p, expanded_mutantAF, distance_for_nonsnp, corr_threshold_recombination))
        else:
            # Set up multiprocessing
            with Pool(processes=ncpu) as pool:
                nonsnp_idx = pool.starmap(process_recombinant_per_pos, [(i, expanded_snv_p, expanded_mutantAF, distance_for_nonsnp, corr_threshold_recombination) for i in range(len(expanded_snv_p))])
        if nonsnp_idx != []:
            nonsnp_idx = np.concatenate(nonsnp_idx).astype(int)
        ## mask all mutantAF with evidence for recombination site as sites with no support
        expanded_mutantAF[nonsnp_idx, :] = 0 ## note every nucleotide per position has its own row --> just set the nt's AF to 0 where there is evidence for a recombination event
        nonrecomb_mutantAF = np.zeros_like(mutantAF)
        collapsed_mant = mant.copy()

        for single_nt in nts_idx:
            is_maNT_in_sample = (mant == single_nt)
            is_nt_in_mutantAF = (p_mapper == single_nt)
            pos_has_allele = pos_idx[is_nt_in_mutantAF]
            nonrecomb_mutantAF[pos_has_allele, :] += np.where(is_maNT_in_sample[pos_has_allele, :], expanded_mutantAF[is_nt_in_mutantAF, :], 0) ## expanded_mutantAF = for every nt a separate row per position exists --> mask to just the current nucleotide and set all AF to 0 where 
            is_recombinant_site = ((is_maNT_in_sample[pos_has_allele, :]) & (nonrecomb_mutantAF[pos_has_allele, :] == 0))
            collapsed_mant[pos_has_allele, :] = np.where(~is_recombinant_site, collapsed_mant[pos_has_allele, :], 4)
            
        ## identify differents stats
        all_nonsnp_idx = np.where(np.all(nonrecomb_mutantAF == 0, axis = 1))[0]
        evidence_of_recomb = (mant != collapsed_mant)
        mixed_snp_nonsnp_idx = np.where(np.any(nonrecomb_mutantAF > 0, axis = 1) & (np.any(evidence_of_recomb, axis = 1)))[0]
        print(f'\nIdentified {len(all_nonsnp_idx)}/{str(num_p)} ({str(len(all_nonsnp_idx)/num_p*100)}%) sites with exclusive evidence for recombination and {len(mixed_snp_nonsnp_idx)}/{str(num_p)} ({str(len(mixed_snp_nonsnp_idx)/num_p*100)}%) sites with evidence for recombination within a subset of samples')
        all_nonsnp_bool = np.zeros(num_p) # build bool array of length p with nonsnps == True
        if all_nonsnp_bool.size > 0:
            all_nonsnp_bool[all_nonsnp_idx] = 1
        all_nonsnp_bool = all_nonsnp_bool.astype(bool)
        
        return [all_nonsnp_idx, all_nonsnp_bool, evidence_of_recomb, mixed_snp_nonsnp_idx, nonrecomb_mutantAF] ## note: all_nonsnp_idx and bool just reports sites where all sites have evidence for recombination, sites with mixed evidence remain but are flagged in the 3rd matrix of bools!
    else:
        if ncpu == 1:
            nonsnp_idx = []
            for i in range(num_p):
                nonsnp_idx.append(process_recombinant_per_pos(i, snv_p, mutantAF, distance_for_nonsnp, corr_threshold_recombination))
        else:
            # Set up multiprocessing
            with Pool(processes=ncpu) as pool:
                nonsnp_idx = pool.starmap(process_recombinant_per_pos, [(i, snv_p, mutantAF, distance_for_nonsnp, corr_threshold_recombination) for i in range(num_p)])
        if nonsnp_idx != []:
            nonsnp_idx = np.concatenate(nonsnp_idx).astype(int)
            
        nonsnp_idx=np.unique(nonsnp_idx)
        ## print number of sites which have indication for recombination
        print(f'\n{str(nonsnp_idx.shape[0])} of a total {str(num_p)} ({str(len(nonsnp_idx)/num_p*100)}%) positions in goodpos were found to be recombinant.')
        nonsnp_bool = np.zeros(num_p) # build bool array of length p with nonsnps == True
        if nonsnp_idx.size > 0:
            nonsnp_bool[nonsnp_idx] = 1
        nonsnp_bool = nonsnp_bool.astype(bool)
        return [nonsnp_idx, nonsnp_bool]
    
def findrecombinantSNPs(p, mutantAF, distance_for_nonsnp, corr_threshold_recombination, ncpu = 1):
    ''' Mutations that  covary in close physical proximity (unexpected by chance) are likely non-independent mutations due to HGT event '''
    # up to 4 cores are recommended, but above that the I/O gets the restricting factor and can make the improvement worse!
    print(f'{len(p)} sites will be evaluated for indications of recombination events.')

    #look for recombination regions
    if ncpu == 1:
        #nonsnp_idx = np.zeros(0,dtype='int')
        nonsnp_idx = []
        for i in range(len(p)):
            if i%2000 == 0: ## print pattern to stdout to show where in evaluation process runs 
                if i%10000 == 0:
                    print('|', end = '')
                else:
                    print('.', end = '')
            gp = p[i]
            #find nearby snps
            if gp > distance_for_nonsnp:
                region = np.array(np.where((p > gp - distance_for_nonsnp) & (p < gp + distance_for_nonsnp)) ).flatten()
                if len(region)>1: 
                    r = mutantAF[region,:]
                    corrmatrix = np.corrcoef(r)
                    [a,b]=np.where(corrmatrix > corr_threshold_recombination)
                    #[a,b]=np.where(np.abs(corrmatrix)  > 0.5) #> corr_threshold_recombination)
                    #nonsnp_idx=np.concatenate((nonsnp_idx,region[a[np.where(a!=b)]]))
                    nonsnp_idx.append(region[a[np.where(a!=b)]])
    else:
        # Set up multiprocessing
        with Pool(processes=ncpu) as pool:
            nonsnp_idx = pool.starmap(process_recombinant_per_pos, [(i, p, mutantAF, distance_for_nonsnp, corr_threshold_recombination) for i in range(len(p))])
    if nonsnp_idx != []:
        nonsnp_idx = np.concatenate(nonsnp_idx).astype(int)
        
    nonsnp_idx=np.unique(nonsnp_idx)
    ## print number of sites which have indication for recombination
    print('\n' + str(nonsnp_idx.shape[0]) + ' of a total ' + str(p.shape[0]) + ' ('  + str(nonsnp_idx.shape[0]/p.shape[0]*100) + '%) positions in goodpos were found to be recombinant.')
    nonsnp_bool = np.zeros(p.shape) # build bool array of length p with nonsnps == True
    if nonsnp_idx.size > 0:
        nonsnp_bool[nonsnp_idx] = 1
    nonsnp_bool = nonsnp_bool.astype(bool)

    return [nonsnp_idx, nonsnp_bool]


def read_genome(ref_genome_folder):
        ## read reference in
        with open(f'{ref_genome_folder}/genome.fasta', 'r') as fasta:
            genome_dict = {record.name: list(str(record.seq)) for record in SeqIO.parse(fasta, "fasta")}
        return genome_dict

def generate_pseudogenome(genome_dict, calls_ingroup, sampleNames_ingroup, hasmutation_ingroup, chrpos_all, scafNames, lineage_dict, contig_padding = 0):
    ## generate pseudogenomes using SNVs only and write an alignment file used for Gubbins as input to identify recombination sites
    alignment_records = []
    
    for lineage_name, lineage in lineage_dict.items():
        for sample_idx, sample in enumerate(sampleNames_ingroup):
            if not sample in lineage:
                continue
            genome_dict_sample = genome_dict.copy()
            smpl_snv_l = np.where(hasmutation_ingroup[:, sample_idx])[0]
            for snv_idx in smpl_snv_l:
                chr_idx = chrpos_all[snv_idx][0] -1
                pos_idx = chrpos_all[snv_idx][1]
                genome_dict_sample[scafNames[chr_idx]][pos_idx] = str(idx2nts(calls_ingroup[snv_idx, sample_idx])) # ATCGN translation

            ## align sequences --> Note: all sequences have the same lengths --> generate dict of alignments of pseudogenomes with contigs concatenated and padded by the max recombination size
            padding_seq = 'N'*contig_padding
            concat_genome = Seq(''.join(["".join(seq) + padding_seq for seq in genome_dict_sample.values()]))
            alignment_records.append(SeqRecord.SeqRecord(concat_genome, id=sample, description=""))

        # Write to output
        #SeqIO.write(alignment_records, f"{lineage_name}_aligned_pseudogenomes.fasta", "fasta")
        with open(f"{lineage_name}_aligned_pseudogenomes.fasta", "w") as out_handle:
            writer = SeqIO.FastaIO.FastaWriter(out_handle, wrap=None)  # disable wrapping to have one sequence per line (important for clonalframe ml!)
            writer.write_file(alignment_records)

def run_gubbins(lineage_dict, gubbins_params):
    ## get the conda path were gubbins is installed
    ## gubbins installation: `conda create -n "gubbins" -c bioconda -c conda-forge -c r gubbins`

    conda_path = os.environ.get('CONDA_EXE').replace('bin/conda', 'bin/activate') ## get conda path
    activate_conda = f" source {conda_path} {gubbins_params['conda_name']} ; "
    for lineage_name, _ in lineage_dict.items():
        tree_builder_cmds = f"--first-tree-builder {gubbins_params['first-tree-builder']} --first-model {gubbins_params['first-model']} --tree-builder {gubbins_params['tree-builder']} --model {gubbins_params['model']}"
        recomb_detection_options = f"--min-snps {gubbins_params['min-snps']} --min-window-size {gubbins_params['min-window-size']} --max-window-size {gubbins_params['max-window-size']}"
        additional_params = f"--threads {gubbins_params['threads']} --iterations {gubbins_params['iterations']} --seed {gubbins_params['seed']} {gubbins_params['kwargs']}"
        gubbins_call = f" run_gubbins.py --prefix {lineage_name} {tree_builder_cmds} {recomb_detection_options} {additional_params} {lineage_name}_aligned_pseudogenomes.fasta ; "
    
        print(f'Running gubbins with command:\n{activate_conda} {gubbins_call}\n\n\n')
        subprocess.run(f'{activate_conda} {gubbins_call}', shell = True, executable='/bin/bash')


def count_sample_per_pos(ranges):
    ## function uses ranges from recombination site detection and generates a list of ranges with the number of samples masked as recombinant thereafter
    ## The function takes care of overlapping regions even those which are overlapping from the same sample (e.g sample_a has region 8-12 and 10-15)
    
    # build lookup dict for star/stop events
    events = defaultdict(list)
    for start, stop, ids in ranges:
        events[int(start)].append(("add", ids.split(',')))
        events[int(stop)].append(("remove", ids.split(',')))
    
    active = set()
    last_pos = None
    segments = []
    # loop over every postion wgere range begin/end
    for pos in sorted(events.keys()):
        ## if segment overlaps with current one save previous segment until the overlap 
        if (last_pos is not None) and (last_pos < pos - 1) and active:
            segments.append((last_pos, pos - 1, sorted(active)))
        # update samples which overlap with the current position and thereafter
        for action, ids in events[pos]:
            if action == "add":
                active.update(ids)
            elif action == "remove":
                active.difference_update(ids)
        # if ids are stored, save the current position to check if overlap exists
        if active:
            last_pos = pos
        else:
            last_pos = None
    return segments


def plot_genomic_site_freq_masked_by_gubbins(recomb_features, hasmutation_ingroup, calls_ingroup, sampleNames_ingroup, genomeLength, ambiguity_filter = ''):

    segments = count_sample_per_pos(recomb_features)
    ## calculate the frequency per genomic position
    frac_recombinant_pos = np.zeros(genomeLength, dtype = float)
    for seg in segments:
        start, stop, ids = seg
        freq = len(ids) / len(sampleNames_ingroup)
        frac_recombinant_pos[start:stop+1] = freq

    ## plot the frequency of recmobinant sites identified acorss the entire genome per position
    fig, ax = plt.subplots(figsize = (6, 4))
    ax.hist(frac_recombinant_pos, bins = np.arange(0, 1.01, 0.01), cumulative = True, density = True, color = 'darkgrey', edgecolor = 'k', lw = 0.2, zorder = 10)
    if ambiguity_filter != '':
        ax.axvline(ambiguity_filter, lw = 1, ls = '--', c = 'k', zorder = 15)
        fraction_masked = f'Ambiguity ≥{ambiguity_filter}\nGenomic sites masked: {np.round(np.mean(frac_recombinant_pos >= ambiguity_filter)*100, 2)}%'
        ax.text(x = 0.98, y = 0.02, s = fraction_masked, va = 'bottom', ha = 'right', zorder = 20)
    ax.set_xlim((0, 1))
    ax.set_ylim((0, 1))
    ax.set_xticks(np.arange(0, 1.1, 0.1))
    ax.grid(which = 'major', axis = 'both', zorder = 0)
    ax.set_xlabel('Fraction of samples with recombinant calls')
    ax.set_ylabel('Fraction of genomic sites')
    ax.set_title('Cumultative fraction of genomic sites identified with\nrecombinant regions across all ingroup samples')
    plt.tight_layout()
    fig.savefig('pdf/gubbins_cumultative_fraction_samples_maskedrecombinant_per_genomic_location.pdf')

    fig, ax = plt.subplots(figsize = (6, 4))
    ax.hist(frac_recombinant_pos, bins = np.arange(0, 1.01, 0.01), density = True, color = 'darkgrey', edgecolor = 'k', lw = 0.2, zorder = 10)
    if ambiguity_filter != '':
        ax.axvline(ambiguity_filter, lw = 1, ls = '--', c = 'k', zorder = 15)
        fraction_masked = f'Ambiguity ≥{ambiguity_filter}\nGenomic sites masked: {np.round(np.mean(frac_recombinant_pos >= ambiguity_filter)*100, 2)}%'
        ax.text(x = 0.98, y = 0.98, s = fraction_masked, va = 'top', ha = 'right', zorder = 20, transform=ax.transAxes)
    ax.set_xlim((0, 1))
    ax.set_xticks(np.arange(0, 1.1, 0.1))
    ax.grid(which = 'major', axis = 'both', zorder = 0)
    ax.set_xlabel('Fraction of samples with recombinant calls')
    ax.set_ylabel('Fraction of genomic sites [%]')
    ax.set_title('Fraction of genomic sites identified with\nrecombinant regions across all ingroup samples')
    plt.tight_layout()
    fig.savefig('pdf/gubbins_fraction_samples_maskedrecombinant_per_genomic_location.pdf')

    ## plot the frequency of ambigous sites
    candpos = np.where( np.sum(hasmutation_ingroup, axis=1) > 0 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!
    fig, ax = plt.subplots(figsize = (6, 4))
    ax.hist(np.sum(calls_ingroup[candpos, :] == 4, axis = 1)/len(sampleNames_ingroup), bins = np.arange(0, 1.01, 0.01), color = 'darkgrey', edgecolor = 'k', lw = 0.2, zorder = 10)
    if ambiguity_filter != '':
        ax.axvline(ambiguity_filter, lw = 1, ls = '--', c = 'k', zorder = 15)
        num_sites_ambig = np.sum(calls_ingroup[candpos, :] == 4, axis = 1)  >= (len(sampleNames_ingroup)*ambiguity_filter)
        fraction_masked = f'Ambiguity ≥{ambiguity_filter}\nCandidate SNVs masked: {np.sum(num_sites_ambig)}'
        ax.text(x = 0.98, y = 0.98, s = fraction_masked, va = 'top', ha = 'right', zorder = 20, transform=ax.transAxes)
    ax.set_xlim((0, 1))
    ax.set_xticks(np.arange(0, 1.1, 0.1))
    ax.grid(which = 'major', axis = 'both', zorder = 0)
    ax.set_xlabel('Fraction of ambiguous candpos sites after recombination site masking\n(without following ambiguity filter)')
    ax.set_ylabel('Number of SNVs')
    plt.tight_layout()
    fig.savefig('pdf/candpos_ambiguity_frequency_post_gubbins_masking.pdf')


def mask_recombinant_sites_gubbins(calls, p, sampleNames, chr_starts, padding_size, lineage_dict, mask_recomb_per_sample = True):
    ## mask recombination sites in calls and set them to 4 (= N)
    ## if mask_per_sample = True: mask only the sample which has the recombination site identified
    ## if mask_per_sample = False: mask the identified positions within identified recombination locus across all sites
    padding_length = (p2chrpos(p, chrStarts)[:, 0]-1)*padding_size ## get the amount of nucleotides included to pad the contigs
    p_adj = p + padding_length ## adjust p for padding done for keeping contigs separated
    chr_starts_padded = chr_starts + ( padding_size*np.arange(len(chr_starts)) ) 

    calls_masked_recomb = calls.copy()
    recomb_features = []
    ## read in gubbin's flagged sites 
    if mask_recomb_per_sample:
        masked_num_sites = {lineage_name: {sid: 0 for sid in lineage_samples} for lineage_name, lineage_samples in lineage_dict.items()}
    else:
        masked_num_sites = {lineage_name: 0 for lineage_name, _ in lineage_dict.items()}
    for lineage_name, lineage_samples in lineage_dict.items():
        file_has_content = False
        try:
            with open(f'{lineage_name}.recombination_predictions.gff') as gff_file:
                for record in GFF.parse(gff_file):
                    if not file_has_content:
                        file_has_content = True
                    for feature in record.features:
                        calls_to_flag = np.where((p_adj >= int(feature.location.start)) & (p_adj <= int(feature.location.end)-1))[0]
                        if mask_recomb_per_sample:
                            sampleid_with_recombination = [sid for sid in feature.qualifiers['taxa'][0].split(' ') if sid in lineage_samples]
                            sampleid_with_recombination_idx = np.array([idx for idx, sid in enumerate(sampleNames) if sid in sampleid_with_recombination])
                            calls_masked_recomb[ np.ix_(calls_to_flag, sampleid_with_recombination_idx) ] = 4 ## flag all recombinant sites
                        else:
                            calls_masked_recomb[ calls_to_flag, : ] = 4 ## flag all recombinant sites
                        ## store all information to identiy frequency of recombinant sites
                        recomb_chr = np.where(chr_starts_padded > int(feature.location.start))[0][0] ## chromosome id is 1-based, here we get the chr start of the next in 0-based --> chr-id carrying the recomb site reported as 1-based
                        recomb_pos_start = int(feature.location.start) - (padding_size * (recomb_chr-1)) ## -1 to apply 0-based
                        recomb_pos_end = int(feature.location.end) - (padding_size * (recomb_chr-1)) ## -1 to apply 0-based
                        taxa = ','.join([t for t in feature.qualifiers['taxa'][0].strip().split(' ') if t != ''])
                        recomb_features.append([int(recomb_pos_start), int(recomb_pos_end), taxa])
                ## get number of positions which are different from before
                if not file_has_content:
                    print(f'No recombinant sites identified for lineage cluster {lineage_name}')
        except:
            print(f'No recombinant sites identified for lineage cluster {lineage_name}')
        if mask_recomb_per_sample:
            masked_num_sites[lineage_name] = {smpl_id: sum(calls_masked_recomb[:, sidx] != calls[:, sidx]) for sidx, smpl_id in enumerate(sampleNames) if smpl_id in lineage_samples}
        else:
            first_sample_of_lineage = min([sidx for sidx, smpl_id in enumerate(sampleNames) if smpl_id in lineage_samples])
            masked_num_sites[lineage_name] = sum((calls_masked_recomb != calls)[first_sample_of_lineage])
        ## sort the recombinant features by start and stop 
        recomb_features = np.array(recomb_features)
        if len(recomb_features) > 0:
            recomb_features = recomb_features[np.lexsort((recomb_features[:,1].astype(int), recomb_features[:,0].astype(int)))]
    return calls_masked_recomb, masked_num_sites, recomb_features

def identify_recombination_sites_via_gubbins(calls_ingroup, p, sampleNames_ingroup, hasmutation_ingroup, ref_genome_folder, chrStarts, scafNames, genomeLength, gubbins_params, lineage_dict, max_fraction_ambigious_samples, mask_recomb_per_sample = True):
    ## NOTE: Please incoporare to subset to lineage prior to call of function to simplify!!!
    ## read reference genome as dict with sequence per chromososme
    genome_dict = read_genome(ref_genome_folder)
    chrpos_all = p2chrpos(p, chrStarts)
    padding_size = gubbins_params['max-window-size']

    ## insert all identified and filtered SNVs to a pseudogenome using SNVs only
    generate_pseudogenome(genome_dict, calls_ingroup, sampleNames_ingroup, hasmutation_ingroup, chrpos_all, scafNames, lineage_dict, padding_size)
    print('Whole genome alignment of pseudogenomes written')

    ## run gubbins per lineage within lineage_dict
    run_gubbins(lineage_dict, gubbins_params)
    
    ## get the identified sites and mask them in the calls matrix as ambigous (=4)
    calls_ingroup, masked_num_sites, recomb_features = mask_recombinant_sites_gubbins(calls_ingroup, p, sampleNames_ingroup, chrStarts, padding_size, lineage_dict, mask_recomb_per_sample)

    ## plot the frequencies of sites identified to be recombinant by gubbins
    if len(recomb_features) > 0:
        plot_genomic_site_freq_masked_by_gubbins(recomb_features, hasmutation_ingroup, calls_ingroup, sampleNames_ingroup, genomeLength, ambiguity_filter = max_fraction_ambigious_samples)
    
    ## flag sites which are now below the ambiguity filter
    siteFilt = (calls_ingroup>3).sum(axis=1) >= (len(sampleNames_ingroup) * max_fraction_ambigious_samples)
    calls_ingroup[ siteFilt ,:] = 4 # sites that fail qc -> 4, for all samples incl. outgroup   

    hasmutation_ingroup = hasmutation_ingroup & (calls_ingroup != 4) ## update hasmutation to excldue all recombination sites identified

    for lineage_name in lineage_dict.keys():
        if mask_recomb_per_sample:
            print(f'Identified between {min(list(masked_num_sites[lineage_name].values()))} and {max(list(masked_num_sites[lineage_name].values()))} sites (median: {np.median(list(masked_num_sites[lineage_name].values()))}) per sample')
        else:
            print(f'Identified {masked_num_sites[lineage_name]} within samples for lineage cluster {lineage_name}')
    print(f'After masking, {sum(siteFilt)} are now below the ambiguity filter and have been removed')
    print('Masked identified recombination sites')
    
    return calls_ingroup, hasmutation_ingroup

def read_trans_trav_rate_from_tree(tree):
    ## get kappa (transistions / transversion rate) from bestmodel 
    with open(tree, 'r') as fid:
        line = fid.read().strip()
        gtr_vals = np.array(line.split('+')[0].split('{')[1].split('}')[0].split('/')).astype(float)
        ## GTR frequency order: A-C, A-G, A-T, C-G, C-T, G-T; see https://github.com/amkozlov/raxml-ng/wiki/Input-data#state-encoding--order
        transitions_bool = np.array([False, True, False, False, True, False])
        kappa = np.mean(gtr_vals[transitions_bool]) / np.mean(gtr_vals[~transitions_bool])
    return kappa

def run_clonalframeml(alignment, kappa, clonalframeml_params, ml_prefix = ''):
    
    conda_path = os.environ.get('CONDA_EXE').replace('bin/conda', 'bin/activate') ## get conda path
    activate_conda = f" source {conda_path} {clonalframeml_params['conda_name']} ; "
    clonalframeml_cmd = f" ClonalFrameML {ml_prefix}.raxml.bestTree {alignment} {ml_prefix}_clonalframe.output -kappa {kappa} -emsim 100 > {ml_prefix}_clonalframe.log ; " ## Note: do not feed the collapsed tree as clonalframe expects a fully bifurcating tree!
    subprocess.run([f"{activate_conda} {clonalframeml_cmd}"], shell=True, executable='/bin/bash')
    return file_prefix

def get_tips(clade):
    if clade.is_terminal():
        return [clade.name]
    tips = []
    for subclade in clade.clades:
        tips.extend(get_tips(subclade))
    return tips


def mask_recombinant_sites_clonalframeml(clonalframe_ml_prefix, sampleNames_ingroup, p_candpos, calls_candpos_ingroup, hasmutation_candpos_ingroup, max_fraction_ambigious_samples):
    tree = Phylo.read(f'{os.getcwd()}/{clonalframe_ml_prefix}.labelled_tree.newick', 'newick')
    node_to_tips = {}

    # get for each node the tips to flag sequences in all samples individually
    for clade in tree.find_clades():
        tips = get_tips(clade)
        node_to_tips[clade.name] = tips

    masked_sites_p_sample = np.zeros((len(sampleNames_ingroup)))
    ## flag for the respective samples the calls
    with open(f'{clonalframe_ml_prefix}.importation_status.txt', 'r') as fid:
        for lid, line in enumerate(fid):
            if lid == 0:
                continue
            recomb_smpl, recomb_start, recomb_end = line.strip().split()
            smpl_have_recomb_bool = np.isin(sampleNames_ingroup, node_to_tips[recomb_smpl])
            calls_to_flag = np.where((p_candpos >= int(recomb_start)) & (p_candpos <= int(recomb_end)))[0]
            calls_candpos_ingroup[np.ix_(calls_to_flag, smpl_have_recomb_bool)] = 4

    masked_sites_p_sample = np.sum(hasmutation_candpos_ingroup & (calls_candpos_ingroup == 4), axis = 0)

    siteFilt = (calls_candpos_ingroup>3).sum(axis=1) >= (len(sampleNames_ingroup) * max_fraction_ambigious_samples)
    calls_candpos_ingroup[ siteFilt ,:] = 4 # sites that fail qc -> 4, for all samples incl. outgroup   

    hasmutation_candpos_ingroup = hasmutation_candpos_ingroup & (calls_candpos_ingroup != 4) ## update hasmutation to excldue all recombination sites identified
    
    print(f'Identified between {min(masked_sites_p_sample)} and {max(masked_sites_p_sample)} sites (median: {np.median(masked_sites_p_sample)}) per sample')
    print(f'After masking, {sum(siteFilt)} are now below the ambiguity filter and have been removed')
    print('Masked identified recombination sites')
    return calls_candpos_ingroup, hasmutation_candpos_ingroup

def identify_recombination_sites_via_clonalframeml(calls_candpos_ingroup, p_candpos, sampleNames_ingroup, hasmutation_candpos_ingroup, ref_genome_folder, chrStarts, scafNames, clonalframeml_params, lineage_dict, max_fraction_ambigious_samples):
    ## NOTE: Please incoporare to subset to lineage prior to call of function to simplify!!!
    ## read reference genome as dict with sequence per chromososme
    genome_dict = read_genome(ref_genome_folder)
    chrpos_all = p2chrpos(p, chrStarts)
    
    ## insert all identified and filtered SNVs to a pseudogenome using SNVs only
    generate_pseudogenome(genome_dict, calls_candpos_ingroup, sampleNames_ingroup, hasmutation_candpos_ingroup, chrpos_all, scafNames, lineage_dict, contig_padding = 0)
    pseudogenome_filename = f"{list(lineage_dict.keys())[0]}_aligned_pseudogenomes.fasta"
    print('Whole genome alignment of pseudogenomes written')
    
    ## generate maximum likelihood tree from pseudogenomes 
    ml_prefix = f"{timestamp}_ML_{clonalframeml_params['model'].replace('+', '')}_candpos"
    raxml_ng_call = f"raxml-ng --msa {pseudogenome_filename} --model {clonalframeml_params['model']} --prefix {ml_prefix} --threads {clonalframeml_params['threads']} --seed {clonalframeml_params['seed']}"
    subprocess.run([raxml_ng_call],shell=True)
    print('Generated tree')
    
    ## get kappa from tree
    kappa = read_trans_trav_rate_from_tree(f'{ml_prefix}.raxml.bestModel')
    print('Got kappa')

    ## run gubbins per lineage within lineage_dict
    results_prefix = run_clonalframeml(pseudogenome_filename, kappa, clonalframeml_params, ml_prefix)
    
    ## get the identified sites and mask them in the calls matrix as ambigous (=4)
    mask_recombinant_sites_clonalframeml(results_prefix, sampleNames_ingroup, p_candpos, calls_candpos_ingroup, hasmutation_candpos_ingroup, max_fraction_ambigious_samples)
    
    return calls_candpos_ingroup, hasmutation_candpos_ingroup

def identify_overlapping_recomb_regions_across_two_lineages(masked_regions_lineage1, masked_regions_lineage2):
    region_lineage1_idx, region_lineage2_idx = 0, 0
    overlapping_recomb_regions_2lineages = []
    while region_lineage1_idx < len(masked_regions_lineage1) and region_lineage2_idx < len(masked_regions_lineage2):
        start = max(masked_regions_lineage1[region_lineage1_idx, 0], masked_regions_lineage2[region_lineage2_idx, 0])
        end   = min(masked_regions_lineage1[region_lineage1_idx, 1], masked_regions_lineage2[region_lineage2_idx, 1])
        if start < end:  # real overlap
            overlapping_recomb_regions_2lineages.append([start, end])
        # move the one that ends earlier
        if masked_regions_lineage1[region_lineage1_idx, 1] < masked_regions_lineage2[region_lineage2_idx, 1]:
            region_lineage1_idx += 1
        else:
            region_lineage2_idx += 1
    return np.array(overlapping_recomb_regions_2lineages)

def identify_overlapping_recomb_regions_across_lineages(merged_recomb_regions_dict):
    lineage_names = list(merged_recomb_regions_dict.keys())
    # set first entry as all overlaps 
    overlapping_recomb_regions = merged_recomb_regions_dict[lineage_names[0]]
    for lineage_name in lineage_names[1:]:
        if np.shape(merged_recomb_regions_dict[lineage_name])[0] == 0:
            overlapping_recomb_regions = np.empty((0, 2), dtype=int)
            break
        overlapping_recomb_regions = identify_overlapping_recomb_regions_across_two_lineages(overlapping_recomb_regions, merged_recomb_regions_dict[lineage_name])
        if len(overlapping_recomb_regions) == 0:  # no overlap left
            break
    return overlapping_recomb_regions

def mask_genes_w_recombinant_sites(annotation_genes, annotation_mutations, lineage_dict, gubbins_params, chrStarts, genomeLength):
    ## mask recombination regions identified by gubbins from annotation mutations and annotation genes
    ## NOTE: This function does only identify regions which are masked from the entire dataset and does 
    ## not work if one has masked recombination sites per sample!
    padding_size = gubbins_params['max-window-size']

    pos_w_mut = (annotation_mutations['p'] + (annotation_mutations['chr']-1)*padding_size).to_numpy()
    
    ## read in gubbin's flagged sites 
    merged_recomb_regions = {lineage_name: [] for lineage_name in lineage_dict.keys()}
    for lineage_name in lineage_dict.keys():
        try:
            recomb_regions = []
            ## read in gff with start and end values
            with open(f'{lineage_name}.recombination_predictions.gff') as gff_file:
                for record in GFF.parse(gff_file):
                    for feature in record.features:
                        if not ((pos_w_mut >= int(feature.location.start)) & (pos_w_mut <= int(feature.location.end)-1)).any(): ## check if the recombintion region has no called snv --> otherwise, region was not removed, but just flagged in a subset of samples, and should not be masked from the remaining genome!!!
                            recomb_regions.append([int(feature.location.start), int(feature.location.end)-1])
            recomb_regions = np.array(recomb_regions)
            ## sort regions
            recomb_regions = recomb_regions[np.argsort(recomb_regions[:, 0])]
            ## merge regions
            current_start, current_end = recomb_regions[0]
            for start, end in recomb_regions[1:]:
                if start <= current_end:  # overlap
                    current_end = max(current_end, end)
                else:
                    merged_recomb_regions[lineage_name].append([current_start, current_end])
                    current_start, current_end = start, end
            merged_recomb_regions[lineage_name].append([current_start, current_end]) ## append last region
            merged_recomb_regions[lineage_name] = np.array(merged_recomb_regions[lineage_name]) ## convert to numpy array
        except:
            merged_recomb_regions[lineage_name] = np.empty((0, 2), dtype=int)
    ## identify only regions which are found across all lineages --> flag only those sites which are masked from entire dataset!!!
    if len(lineage_dict.keys()) > 1:
        overlapping_recomb_regions = identify_overlapping_recomb_regions_across_lineages(merged_recomb_regions)
    else:
        overlapping_recomb_regions = list(merged_recomb_regions.values())[0]
    ## mask regions
    latest_tag_num = 0
    tag_num_dict = {}
    for chr, chr_start in enumerate(chrStarts):
        ## identify the correct region start and end by removing padded regions
        padded_region_chr_start = chr_start + chr*padding_size
        if chr < len(chrStarts)-1:
            padded_region_chr_end = chrStarts[chr+1]-1 + chr*padding_size
        else:
            padded_region_chr_end = genomeLength-1 + chr*padding_size
        ## identify regions on chr
        region_on_chr_idx = np.where((overlapping_recomb_regions[:, 0] >= padded_region_chr_start) & (overlapping_recomb_regions[:, 1] <= padded_region_chr_end))[0]
        masked_region_on_chr = overlapping_recomb_regions[region_on_chr_idx, ]
        masked_region_on_chr_p = masked_region_on_chr - padded_region_chr_start ## convert to positions on chromosome!
        for region in masked_region_on_chr_p:
            ## mask annotation_genes which overlap with recombinant region
            if not annotation_genes[chr].empty:
                gene_start = annotation_genes[chr]['loc1'].values
                gene_end = annotation_genes[chr]['loc2'].values
                gene_overlaps_region_start = np.where((gene_start < region[0]) & (region[0] < gene_end))[0]
                gene_overlaps_region_end = np.where((gene_start < region[1]) & (region[1] < gene_end))[0]
                gene_overlaps_region_full = (gene_start >= region[0]) & (gene_end <= region[1])
                ## truncate genes which are overlapping
                for overlap in gene_overlaps_region_start:
                    annotation_genes[chr].loc[overlap, 'loc2'] = region[0]
                    annotation_genes[chr].loc[overlap, 'indices'][1] = region[0]
                    seq_len = annotation_genes[chr].loc[overlap, 'loc2'] - annotation_genes[chr].loc[overlap, 'loc1']
                    annotation_genes[chr].loc[overlap, 'sequence'] = str(annotation_genes[chr].loc[overlap, 'sequence'][:int(seq_len)]) ## pandas dislikes iterable seq objects --> convert to string and later back to seq objects
                    annotation_genes[chr].loc[overlap, 'translation'] = str(annotation_genes[chr].loc[overlap, 'translation'][:int(seq_len/3)]) ## get only the full codons covered! ; pandas dislikes iterable seq objects --> convert to string and later back to seq objects
                for overlap in gene_overlaps_region_end:
                    annotation_genes[chr].loc[overlap, 'loc1'] = region[1]
                    annotation_genes[chr].loc[overlap, 'indices'][0] = region[1]
                    seq_len = annotation_genes[chr].loc[overlap, 'loc2'] - annotation_genes[chr].loc[overlap, 'loc1']
                    annotation_genes[chr].loc[overlap, 'sequence'] = str(annotation_genes[chr].loc[overlap, 'sequence'][-int(seq_len):]) ## pandas dislikes iterable seq objects --> convert to string and later back to seq objects
                    annotation_genes[chr].loc[overlap, 'translation'] = str(annotation_genes[chr].loc[overlap, 'translation'][-int(seq_len/3):]) ## get only the full codons covered! ; pandas dislikes iterable seq objects --> convert to string and later back to seq objects
                annotation_genes[chr] = annotation_genes[chr][~gene_overlaps_region_full].reset_index(drop = True) ## drop all which overlap with region
                ## convert seq columns back to seq objects 
                annotation_genes[chr]['sequence'] = annotation_genes[chr]['sequence'].apply(Seq)
                annotation_genes[chr]['translation'] = annotation_genes[chr]['translation'].apply(Seq)
                ## translate gene_num_global
            ## mask annotation_mutations which overlap with recombinant region (shorten their sequence and check if masking is valid (no mutation is allowed to be masked!!!))
            if not annotation_mutations.empty:
                mut_on_chr = (annotation_mutations['chr'] == chr+1).values
                if 'loc1' in annotation_mutations.columns:
                    mut_intragenic = ~annotation_mutations['loc1'].isna().values ## just select those which are intragenic!
                    mut_cds = ~annotation_mutations['translation'].isna().values ## get for translation also info which are CDS
                else:
                    mut_intragenic = np.zeros(len(annotation_mutations), dtype = bool) ## all are intergenic
                if (sum(mut_on_chr & mut_intragenic) > 0): ## at least one mutation on current chr and intragenic
                    gene_start = annotation_mutations.loc[mut_on_chr & mut_intragenic]['loc1'].values
                    gene_end = annotation_mutations.loc[mut_on_chr & mut_intragenic]['loc2'].values
                    mut_pos = annotation_mutations.loc[mut_on_chr & mut_intragenic]['p'].values
                    gene_overlaps_region_start = np.where((gene_start < region[0]) & (region[0] < gene_end))[0]
                    gene_overlaps_region_end = np.where((gene_start < region[1]) & (region[1] < gene_end))[0]
                    mut_overlaps_region = (mut_pos >= region[0]) & (mut_pos <= region[1])
                    ## convert idx back to df idx to remove unecessary double subsetting 
                    gene_overlaps_region_start_idx_full_df = np.where(mut_on_chr & mut_intragenic)[0][gene_overlaps_region_start]
                    gene_overlaps_region_end_idx_full_df = np.where(mut_on_chr & mut_intragenic)[0][gene_overlaps_region_end]
                    ## truncate genes which are overlapping
                    for overlap in gene_overlaps_region_start_idx_full_df:
                        annotation_mutations.loc[overlap, 'loc2'] = region[0]
                        seq_len = annotation_mutations.loc[overlap, 'loc2'] - annotation_mutations.loc[overlap, 'loc1']
                        annotation_mutations.loc[overlap, 'sequence'] = str(annotation_mutations.loc[overlap, 'sequence'][:int(seq_len)]) ## pandas dislikes iterable seq objects --> convert to string and later back to seq objects
                        if annotation_mutations.loc[overlap, 'translation'] == annotation_mutations.loc[overlap, 'translation']: ## check if it has seq
                            annotation_mutations.loc[overlap, 'translation'] = str(annotation_mutations.loc[overlap, 'translation'][:int(seq_len/3)]) ## get only the full codons covered! ; pandas dislikes iterable seq objects --> convert to string and later back to seq objects
                    for overlap in gene_overlaps_region_end_idx_full_df:
                        annotation_mutations.loc[overlap, 'loc1'] = region[1]
                        seq_len = annotation_mutations.loc[overlap, 'loc2'] - annotation_mutations.loc[overlap, 'loc1']
                        annotation_mutations.loc[overlap, 'sequence'] = str(annotation_mutations.loc[overlap, 'sequence'][-int(seq_len):]) ## pandas dislikes iterable seq objects --> convert to string and later back to seq objects
                        if annotation_mutations.loc[overlap, 'translation'] == annotation_mutations.loc[overlap, 'translation']: ## check if it has seq
                            annotation_mutations.loc[overlap, 'translation'] = str(annotation_mutations.loc[overlap, 'translation'][-int(seq_len/3):]) ## get only the full codons covered! ; pandas dislikes iterable seq objects --> convert to string and later back to seq objects
                    annotation_mutations.loc[mut_intragenic, 'sequence'] = annotation_mutations.loc[mut_intragenic, 'sequence'].apply(Seq)
                    annotation_mutations.loc[mut_intragenic & mut_cds, 'translation'] = annotation_mutations.loc[mut_intragenic & mut_cds, 'translation'].apply(Seq)
                    if sum(mut_overlaps_region):
                        print('WARNING: A mutation overlaps with recombination regions! SOMETHING WENT WRONG.')
                        print(f'Recombinant region: {region[0]}-{region[1]} on chromosome {chr+1}')
                        print(f'Mutations in region: {annotation_mutations[mut_on_chr & mut_intragenic][mut_overlaps_region].to_string(index = False)}')
                        sys.exit()
        if not annotation_genes[chr].empty:
            gene_has_tagnum = annotation_genes[chr]['tagnumber'] > 0
            tags_old = annotation_genes[chr].loc[gene_has_tagnum, 'tagnumber'].values
            tags_new = annotation_genes[chr].loc[gene_has_tagnum, 'tagnumber'].index+1 + latest_tag_num
            tag_num_dict_tmp = {tag_old: tag_new for tag_old, tag_new in zip(tags_old, tags_new)}
            annotation_genes[chr].loc[gene_has_tagnum, 'tagnumber'] = annotation_genes[chr].loc[gene_has_tagnum, 'tagnumber'].map(tag_num_dict_tmp)
            latest_tag_num = annotation_genes[chr].loc[gene_has_tagnum, 'tagnumber'].max()
            tag_num_dict = {**tag_num_dict, **tag_num_dict_tmp}
    ## translate gene nums in annotation_mutation 
    annotation_mutations['gene_num_global'] = annotation_mutations['gene_num_global'].map(tag_num_dict).fillna(annotation_mutations['gene_num_global'])
    return annotation_genes, annotation_mutations     


def median_above_zero(row):
    positive_values = row[row > 0]  # Filter values > 0
    if positive_values.size == 0:
        return -1  # Return -1 if no values are greater than 0
    else:
        return np.median(positive_values)  # Calculate median for positive values


def save_num_muts_p_sample_num_mutsamples_per_p(hasmutation_bool, sampleNames, p, contig_positions, suffix=''):
    """
    This function safes 2 csv files to save the number of mutated samples per position and number of mutations per sample
    hasmutation_bool : 2d numpy bollean matrix of of samples*positions stating where mutations appear (e.g. (calls[:,ingroup_bool] != refnti_m[:,ingroup_bool]) & (calls[:,ingroup_bool] < 4))
    """
    ## get number of SNPs per sample and number of samples per covered SNP 
    mutsamples_per_pos = np.sum(hasmutation_bool, axis = 1)
    muts_per_sample = np.sum(hasmutation_bool, axis = 0)
    with open(f'num_muts_per_sample{suffix}.csv', 'w') as fid:
        fid.write('sampleid,num_muts\n')
        for sample, num_muts in zip(sampleNames, muts_per_sample):
            fid.write(f'{sample},{num_muts}\n')
    with open(f'num_mutated_smpls_per_pos{suffix}.csv', 'w') as fid:
        fid.write('genome_pos,chr,pos,num_muts\n')
        for pos, chrpos, num_smpls in zip(p, contig_positions, mutsamples_per_pos):
            fid.write(f'{pos},{chrpos[0]},{chrpos[1]},{num_smpls}\n')

def get_num_core_genome_pos(cmtdir, smpl_of_interest_bool, p, goodpos2use, genomeLength, filter_parameter_site_per_sample, filter_parameter_site_across_samples, median_coverage_position_with_cov = False):
    ## get number of positions acorss samples which can be considered to be within the core genome given the SNV thresholds 
    ## The number is required for Felsenstein's ascertainment correction
    ## NOTE: subsetting of the cov_matrix to itself allows good memory usage as this matrix can be huge!
    ########################################
    # LOAD COV MATRIX FOR ML TREES!
    cov_mat_csr = sparse.load_npz(f'{cmtdir}cov_raw_sparsecsr_mat.npz')
    cov_matrix = np.array(cov_mat_csr.todense())
    del cov_mat_csr ## delete variable to save memory
    
    p_masked_idx = np.array([cand_pos for cand_pos_idx, cand_pos in enumerate(p) if cand_pos_idx not in goodpos2use]) ## mask all positions which have been removed by filtering!
    p_mask_bool = np.ones(genomeLength, dtype=bool)
    p_mask_bool[p_masked_idx] = False

    cov_matrix = cov_matrix[smpl_of_interest_bool, :] ## remove all samples which are not of interest
    cov_matrix = cov_matrix[:, p_mask_bool] ## note cov_matrix is not aligned with genome anymore!!!
    min_cov_required = filter_parameter_site_per_sample['min_cov_per_strand_for_call'] * 2 ## acount for fwd and reverse strand!
    max_smpls_ambiguity = len(smpl_of_interest_bool) * filter_parameter_site_across_samples['max_fraction_ambigious_samples']
    non_ambig_pos_bool = np.sum(cov_matrix < min_cov_required, axis = 0) < max_smpls_ambiguity ## flag all positions which would be flagged by ambiguity (0 coverage --> N --> would be masked if SNV identified in at least one sample, but exceeds ambiguity threshold)
    cov_matrix = cov_matrix[:, non_ambig_pos_bool]

    number_pos = np.shape(cov_matrix)[1]
    number_chunks = np.ceil(number_pos / 1e6).astype(int)
    pos_has_suff_cov_across_smpls_l = []
    if median_coverage_position_with_cov:
        for i in range(number_chunks): ## workaround for large datasets which would fail due to memory limits
                start_mask = int(i * 1e6)
                end_mask = int((i+1)*1e6)
                if end_mask > number_pos:
                        end_mask = number_pos
                masked_arr = np.ma.masked_equal(cov_matrix[:, start_mask:end_mask], 0) ## mask 0
                chunk_pos_has_suff_cov_across_smpls = np.ma.median(masked_arr, axis=0).filled(np.nan) >= filter_parameter_site_across_samples['min_median_coverage_position_with_cov']
                pos_has_suff_cov_across_smpls_l.append(chunk_pos_has_suff_cov_across_smpls)
        pos_has_suff_cov_across_smpls = np.concatenate(pos_has_suff_cov_across_smpls_l)
    else:
        pos_has_suff_cov_across_smpls = np.median(cov_matrix, axis = 0) >= filter_parameter_site_across_samples['min_median_coverage_position']

    number_core_pos = np.sum(pos_has_suff_cov_across_smpls)
    number_invariable_pos = number_core_pos - len(goodpos2use)
    number_variable_pos = len(goodpos2use)
    number_unresolved_pos = genomeLength - number_core_pos
    with open('num_pos_coregenome.csv', 'w') as fid:
        fid.write('number_core_pos,number_invariable_pos,number_variable_pos,number_unresolved_pos\n')
        fid.write(f'{number_core_pos},{number_invariable_pos},{number_variable_pos},{number_unresolved_pos}\n')
    return number_core_pos, number_invariable_pos

##############################
## QC plot functions
##############################

def plot_minorAF_rates(minorAF, hasmutation, sampleNames, subject, refgenome, timestamp, titlelabel = ''):
    ## get amount of snp
    # no_snps = np.shape(minorAF)[0]
    no_samples = np.shape(minorAF)[1]
    ## get frequencies of minor AF rates where sample has a true snp
    lower_bounds = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1]
    upper_bounds = [0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 1]
    labels = ['<= ' + str(int(u_bound * 100)) + '%' for u_bound in upper_bounds]

    freq_minorAF = {}
    for l_bound, u_bound in zip(lower_bounds, upper_bounds):
        ## note need to use np.divide to remove via where cases which would divide by 0
        freq_minorAF[u_bound] = np.divide(np.sum(((minorAF > l_bound) & (minorAF <= u_bound)) & hasmutation, axis = 0), 
                                        np.sum(hasmutation, axis = 0), 
                                        where = np.sum(hasmutation, axis = 0) != 0)
    ## plot barplot
    plot_width = no_samples / 5
    if plot_width < 6:
        plot_width = 6
    elif plot_width > 200:
        plot_width = 200
    fig = plt.figure(figsize=(plot_width, 8))
    axe = plt.subplot(111)

    bottom_tracker = np.zeros(np.shape(minorAF)[1])
    for freq_steps in freq_minorAF.values():
        axe.bar(height = freq_steps, bottom = bottom_tracker, x = sampleNames)
        bottom_tracker = bottom_tracker + freq_steps
    
    ## modify plot 
    axe.set_xticks(sampleNames)
    axe.set_xticklabels(sampleNames, rotation=90,fontsize=6)
    axe.legend(title = 'minor allele frequency', labels = labels) #['< 5%', '< 10%', '< 25%', '>= 25%'])
    axe.set_ylabel('rate of observed minor allel frequencies in sample')
    axe.set_xlabel('samples')
    axe.set_title(f'Rates of observed minor allel frequencies per sample in {subject} {refgenome}')
    plt.tight_layout()  
    ## save plot
    fig.savefig(f'pdf/{timestamp}_{subject}_{refgenome}_obs_minorAF_rates{titlelabel}.pdf')
    plt.close('all')


def plot_quals_rates(quals, hasmutation, sampleNames, subject, refgenome, timestamp, titlelabel = ''):
    ## get amount of snp
    no_snps = np.shape(quals)[0]
    no_samples = np.shape(quals)[1]
    ## get frequencies of minor AF rates where sample has a true snp
    upper_bounds = [300, 200, 100, 75, 50, 30]
    lower_bounds = [200, 100, 75, 50, 30, 0]
    labels = ['>= ' + str(l_bound) if l_bound != 0 else '< ' + str(upper_bounds[-1]) for l_bound in lower_bounds]

    freq_quals = {}
    for l_bound, u_bound in zip(lower_bounds, upper_bounds):
        ## note need to use np.divide to remove via where cases which would divide by 0
        freq_quals[u_bound] = np.divide(np.sum(((quals < u_bound) & (quals >= l_bound)) & hasmutation, axis = 0), 
                                        np.sum(hasmutation, axis = 0), 
                                        where = np.sum(hasmutation, axis = 0) != 0)
    ## plot barplot
    plot_width = no_samples / 5
    if plot_width < 6:
        plot_width = 6
    elif plot_width > 200:
        plot_width = 200
    fig = plt.figure(figsize=(plot_width, 8))
    axe = plt.subplot(111)

    bottom_tracker = np.zeros(np.shape(quals)[1])
    for freq_steps in freq_quals.values():
        axe.bar(height = freq_steps, bottom = bottom_tracker, x = sampleNames)
        bottom_tracker = bottom_tracker + freq_steps
    
    ## modify plot 
    axe.set_xticks(sampleNames)
    axe.set_xticklabels(sampleNames, rotation=90,fontsize=6)
    axe.legend(title = 'FQ score', labels = labels) #['>= 200', '>= 100', '>= 30', '< 30'])
    axe.set_ylabel('rate of FQ score in sample')
    axe.set_xlabel('samples')
    axe.set_title(f'Rates of FQ score per sample in {subject} {refgenome}')
    plt.tight_layout()  
    ## save plot
    fig.savefig(f'pdf/{timestamp}_{subject}_{refgenome}_FQ_rates{titlelabel}.pdf')
    plt.close('all')


def plot_count_rates(counts, hasmutation, sampleNames, subject, refgenome, timestamp, titlelabel = ''):
    ## get amount of snp
    no_samples = np.shape(counts)[0]
    ## get frequencies of minor AF rates where sample has a true snp
    upper_bounds = [8, 10, 15, 20, 25, 9999]
    lower_bounds = [0, 8, 10, 15, 20, 25]
    labels = ['>= ' + str(l_bound) if l_bound != lower_bounds[0] else '< ' + str(lower_bounds[0]) for l_bound in lower_bounds]

    ## get coverage on major allel frequency:
    counts_both_strands = counts[:, :4, :] + counts[:, 4:, :] ## calculate coverage per snp on both strands
    counts_maf = np.max(counts_both_strands, axis = 1).T ## extract major allel frequency (note: it is not cleaned for ref nt coverage yet!)
    
    ## subset matrices to only samples which have mutations
    freq_counts = {}
    for l_bound, u_bound in zip(lower_bounds, upper_bounds):
        ## calculate how many snps are between bounds _and_ have a SNP on pos (hasmutation)np.sum(hasmutation, axis = 0) ## note need to use np.divide to remove via where cases which would divide by 0
        freq_counts[l_bound] = np.divide(np.sum(((counts_maf < u_bound) & (counts_maf >= l_bound)) & hasmutation, axis = 0), 
                                        np.sum(hasmutation, axis = 0), 
                                        where = np.sum(hasmutation, axis = 0) != 0)
    
    ## plot barplot
    plot_width = no_samples / 5
    if plot_width < 6:
        plot_width = 6
    elif plot_width > 200:
        plot_width = 200
    fig = plt.figure(figsize=(plot_width, 8))
    axe = plt.subplot(111)

    bottom_tracker = np.zeros(no_samples)
    for freq_steps in freq_counts.values():
        axe.bar(height = freq_steps, bottom = bottom_tracker, x = sampleNames)
        bottom_tracker = bottom_tracker + freq_steps
    
    ## modify plot 
    axe.set_xticks(sampleNames)
    axe.set_xticklabels(sampleNames, rotation=90,fontsize=6)
    axe.legend(title = 'Coverage', labels = labels) #['>= 200', '>= 100', '>= 30', '< 30'])
    axe.set_ylabel('rates of coverage on SNP in sample')
    axe.set_xlabel('samples')
    axe.set_title(f'Rates of coverage on SNP per sample in {subject} {refgenome}')
    plt.tight_layout()  
    ## save plot
    fig.savefig(f'pdf/{timestamp}_{subject}_{refgenome}_cov_rates{titlelabel}.pdf')
    plt.close('all')


def plot_Ns_rate(maNT, sampleNames, pos, genomeLength, num_frames, subject, refgenome, timestamp, titlelabel = ''):
    ## get number of samples
    no_samples = np.shape(maNT)[1]
    no_snps = np.shape(maNT)[0]
    
    ## plot barplot
    plot_width = no_samples / 5
    if plot_width < 6:
        plot_width = 6
    elif plot_width > 200:
        plot_width = 200
    fig = plt.figure(figsize=(plot_width, 8))
    axe = plt.subplot(111)
    frame_bp_Length = genomeLength / num_frames
    ## get number of Ns
    # rate_Ns_sample = np.sum(maNT > 3, axis = 0) / no_snps ## store rate of Ns per sample

    bottom_tracker = np.zeros(no_samples) ## create bottom tracker for bar plot
    
    ## plot 
    for genomic_frame in range(num_frames):
        pos_consider = np.where((pos >= frame_bp_Length * genomic_frame) & (pos < frame_bp_Length * (genomic_frame + 1)))[0] ## get idx of positions to consider 
        rate_snps_genomic_frame = np.sum(maNT[pos_consider, :] > 3, axis = 0) / no_snps ## ratio of mutations in genomic frame
        axe.bar(height = rate_snps_genomic_frame, bottom = bottom_tracker, x = sampleNames) ## plot 
        bottom_tracker = bottom_tracker + rate_snps_genomic_frame ## update bottom for barplot

    ## add description to plot 
    legend_labels = [('< ' + str( np.trunc( ((genomic_frame+1) * frame_bp_Length/1000)*100) /100) ) for genomic_frame in range(num_frames)] ## generate list of genomic frames in kb with two decimals (np.trunc(genomic_loc *100/100))
    axe.legend(title = 'Genomic region [kbp]', labels = legend_labels)
    #axe.bar(height = rate_Ns_sample, x = sampleNames)
    
    ## modify plot 
    axe.set_xticks(sampleNames)
    axe.set_xticklabels(sampleNames, rotation=90,fontsize=6)
    axe.set_ylabel('rate of Ns called on pos with identified SNPs')
    axe.set_xlabel('samples')
    axe.set_title(f'Rates of Ns on positions with an identified SNP {subject} {refgenome}')
    plt.tight_layout()  
    ## save plot
    fig.savefig(f'pdf/{timestamp}_{subject}_{refgenome}_Ns_rate{titlelabel}.pdf')
    plt.close('all')


def number_snps_per_genomic_region(pos, hasmutation, pos_subselect_idx, sampleNames, genomeLength, num_frames, timestamp, subject, refgenome, ratio, titlelabel = ''):
    
    ## plot barplot
    num_samples = len(sampleNames)
    plot_width = num_samples / 5
    if plot_width < 6:
        plot_width = 6
    elif plot_width > 200:
        plot_width = 200
    fig = plt.figure(figsize=(plot_width, 8))
    axe = plt.subplot(111)
    ## divide genome into n frames 
    frame_bp_Length = genomeLength / num_frames
    bottom_tracker = np.zeros(num_samples) ## create bottom tracker for bar plot
    if ratio:
        title_append = 'Ratio'
        no_snps_sample = np.sum(hasmutation[pos_subselect_idx, :], axis = 0) ## get number of all snps in sample 
    else: 
        title_append = 'Number'
    
    ## plot 
    for genomic_frame in range(num_frames):
        pos_consider = np.where((pos >= frame_bp_Length * genomic_frame) & (pos < frame_bp_Length * (genomic_frame + 1)))[0] ## get idx of positions to consider 
        snps_consider = np.intersect1d(pos_subselect_idx, pos_consider) ## intersect goodpos with positions in frame 
        if ratio:
            rate_snps_genomic_frame = np.sum(hasmutation[snps_consider, : ], axis = 0) / no_snps_sample ## ratio of mutations in genomic frame
        else: 
            rate_snps_genomic_frame = np.sum(hasmutation[snps_consider, : ], axis = 0) ## number of mutations in genomic frame
        axe.bar(height = rate_snps_genomic_frame, bottom = bottom_tracker, x = sampleNames) ## plot 
        bottom_tracker = bottom_tracker + rate_snps_genomic_frame ## update bottom for barplot

    ## add description to plot 
    legend_labels = [('< ' + str( np.trunc( ((genomic_frame+1) * frame_bp_Length/1000)*100) /100) ) for genomic_frame in range(num_frames)] ## generate list of genomic frames in kb with two decimals (np.trunc(genomic_loc *100/100))
    axe.legend(title = 'Genomic region [kbp]', labels = legend_labels)
    axe.set_xticks(sampleNames)
    axe.set_xticklabels(sampleNames, rotation=90,fontsize=6)
    axe.set_xlabel('samples')
    axe.set_ylabel(f'{title_append} of SNPs in genomic region')
    axe.set_title(f'{title_append} of SNPs per genomic region per sample in {subject} {refgenome}')
    ## save plot
    plt.tight_layout()
    fig.savefig(f'pdf/{timestamp}_{subject}_{refgenome}_{title_append}_SNPs_per_genomic_region{titlelabel}.pdf')
    plt.close('all')


def onpick(event, countdata, pos, smplnames, nsample, saveplots = False, timestamp='', subject='', refgenome='', cutoff_value='', title_append=''):
    i=event.ind[0]
    plot_width = nsample / 10
    
    drilldownfig = plt.figure(i, clear = True)
        
    formated_count_data=np.stack((countdata[0,:4,i], countdata[0,4:,i]))
    for s in range(1,nsample):
        new=np.stack((np.zeros((4)),countdata[s,:4,i], countdata[s,4:,i])) 
        formated_count_data= np.concatenate((formated_count_data,new))
            
    
    axsub = drilldownfig.subplots()
    stackedbar = pd.DataFrame(formated_count_data )
    stackedbar.plot(kind='bar',stacked=True, width=1, ax=axsub)
    legend(['A','T','C','G'])
    axsub.set_ylabel('counts')
    axsub.set_xticks(ticks=range(0,nsample*3,3), labels=smplnames)
    axsub.set_xticklabels(smplnames, rotation=90, fontsize=6) 
    axsub.set_title('SNP at pos ' + str(pos[i]))
    drilldownfig.tight_layout()
    drilldownfig.show()

    if saveplots:
        plt.gcf().set_size_inches(plot_width, 8)
        plt.tight_layout()
        drilldownfig.savefig(f'pdf/coverage_snp_fwd_rev/{timestamp}_{subject}_{refgenome}_snp{str(pos[i])}_quals_per_sample_{cutoff_value}{title_append}.pdf')


def plot_interactive_scatter_barplots(xcoor,ycoor,xlabel,ylabel,samplenames,countdata,timestamp,cutoff_value,subject,refgenome, min_qual=0, max_qual=None, min_pos = 0, max_pos = None,saveplots = False, title_append=''):
    ''' Interactive scatterplot for troubleshooting classify step 
    Args:
        Ctscov (list): counts_coverage object output from counts_CSS.
        expected_freq (df): sxc df giving true clade frequencies.
        measured_freq (df): sxc df giving frequencies output from counts_CSS.
        strain_freq (TYPE): sxstrains df giving true strain frequencies .

    Returns:
        fig (TYPE): DESCRIPTION.

    Example call:
    apy.plot_interactive_scatter_barplots(p[goodpos],mutQual[goodpos],'pos','qual', sampleNames,counts[:,:,goodpos],timestamp, filter_parameter_site_across_samples['corr_threshold_recombination'], subject_fld_label, refgenome, saveplots = True)    
    '''
    if saveplots:
        if not os.path.exists(os.getcwd() + "/pdf/coverage_snp_fwd_rev"):
        ## generate dir 
            os.mkdir(os.getcwd() + "/pdf/coverage_snp_fwd_rev")

    ## set maxquals and maxpos to default (max values) for plotting, if not set otherwise
    if max_qual is None:
        max_qual == max(ycoor)
    if max_pos is None:
        max_pos == max(xcoor)

    nsample=np.size(countdata,axis=0)

    fig, axs = plt.subplots()
    axs.scatter(xcoor,  ycoor, alpha=0.6, s=12, picker=True)
    axs.set_xlim(xmin=min_pos, xmax=max_pos)
    axs.set_ylim(ymin=min_qual, ymax=max_qual)

    plt.xlabel(xlabel,fontsize=16)
    plt.ylabel(ylabel,fontsize=16)
    plt.title(str(len(xcoor)) + ' goodpos sites defined with an ' + str(cutoff_value) + ' cutoff value')
    plt.tight_layout()
    if saveplots:
        fig.canvas.mpl_connect('pick_event', lambda event: onpick(event, countdata = countdata, pos = xcoor, smplnames = samplenames, nsample = nsample, saveplots = True, timestamp=timestamp, subject=subject, refgenome=refgenome, cutoff_value=cutoff_value, title_append=title_append))
        fig.savefig(f'pdf/coverage_snp_fwd_rev/{timestamp}_{subject}_{refgenome}_quals_vs_pos_{cutoff_value}{title_append}.pdf')
    else:
        fig.canvas.mpl_connect('pick_event', lambda event: onpick(event, countdata = countdata, pos = xcoor, smplnames = samplenames, nsample = nsample))
    return fig


## plot SNP distribution across genome
def plot_snp_distribution(positions, refgenome, chrStarts):
    """
    Analysis and QC Plotting function, plots p matrix across genome, with lines separating chroms/plasmids
    """
    fig = plt.figure(figsize=(12, 5))
    axe = plt.subplot(111)
    ## mark begin of new chromosomes by vertical line
    plt.vlines(chrStarts, ymin = 0, ymax = 2, color = 'black', alpha = 0.5)
    ## plot positions scattered around 1 
    plt.scatter(x = positions, y = np.random.random(len(positions))+0.5, s = 10, alpha = 0.4)
    ## add labels
    plt.xlabel('Genomic position [Mb]')
    axe.axes.yaxis.set_visible(False)
    axe.set_title('SNP distribution across ' + refgenome + ' genome')
    fig.savefig('pdf/snp_distribution_' + refgenome + ".pdf")
    plt.close('all')


## plot which outgroup sample is polarizing which SNP and which SNPs are not polarized at all 
def plot_outgroup_snv_polarization(p_gp, outgroup_sampleNames, outsplsnti, outchimerasplancnti, goodpos, chrStarts, genomeLength, refgenome, plot_fileprefix = ''):      
        ## plot which outgroup sample is polarizing which SNP and which SNPs are not polarized at all 
        outsplsnti_chim = np.concatenate((outsplsnti, np.reshape(outchimerasplancnti, (-1, 1))), axis = 1) ## include the major (chimeric) calls
        outgrpspl_idx_sort = np.argsort(np.sum(outsplsnti_chim[goodpos]==4, axis = 0)) ## idx to sort samples based on the number of SNPs they polarize 
        outgrp_sampleNames = np.append(outgroup_sampleNames, 'Chimeric outgroup') ## extract sample names and add chimeric sample name
        outgrp_sampleNames = outgrp_sampleNames + '\n(unpol. SNPs: ' + (np.sum(outsplsnti_chim[goodpos]==4, axis = 0)).astype(str) + ')' ## add information how many snps not resolved by outgroup
        outgrp_sampleNames_s = outgrp_sampleNames[outgrpspl_idx_sort] ## sort sample names 
        outgrp_sampleNames_sext = np.tile(outgrp_sampleNames_s, len(goodpos)) ## repeat samplenames n times 

        gnome_pos_snps_ext = np.repeat(p_gp, len(outgrpspl_idx_sort)) ## repeat each pos * outgroup sample number 
        
        outgrpspl_pol_pos = outsplsnti_chim[np.ix_(goodpos, outgrpspl_idx_sort)] ## get calls of outgroups and sort by idx 
        #outgrpspl_pol_pos = (outgrpspl_pol_pos != 4)
        outgrpspl_pol_pos_1d = outgrpspl_pol_pos.ravel() ## convert to 1d array for plotting
        outgrpspl_pol_pos_1d[outgrpspl_pol_pos_1d != 4] = True
        outgrpspl_pol_pos_1d[outgrpspl_pol_pos_1d == 4] = False

        ## get cmap   
        cmap = mcolors.LinearSegmentedColormap.from_list("", ['black', 'darkorange'])
        
        ## plot
        plot_width = 8 + len(outgrp_sampleNames)/2
        if plot_width < 6:
            plot_width = 6
        elif plot_width > 200:
            plot_width = 200
        fig, ax = plt.subplots(figsize = (8 + len(outgrp_sampleNames)/2, 10))
        ax.scatter(x = gnome_pos_snps_ext, y = outgrp_sampleNames_sext, c = outgrpspl_pol_pos_1d, alpha = 0.35, s = 60, norm = plt.Normalize(vmin=0, vmax=1), cmap = cmap, zorder = 10)
        ax.vlines(x = chrStarts, ymin=-1, ymax=len(outgrp_sampleNames_s), linestyles='--', linewidth = 0.75, colors='k', alpha = 0.2, zorder = 5)
        legend_info = []
        legend_info += [Line2D([0], [0], markerfacecolor='darkorange', color = 'white', alpha = 0.35, marker = 'o', label='covered', markersize = 10)]
        legend_info += [Line2D([0], [0], markerfacecolor='black', color = 'white', alpha = 0.35, marker = 'o', label='uncovered', markersize = 10)]
        ax.legend(title = 'SNP', handles = legend_info, loc='lower center', bbox_to_anchor=(0.5, -0.35), fancybox=False, shadow=False)
        ax.set_ylim(-0.5, len(outgrp_sampleNames_s)-0.5)
        ax.set_xlim(0, genomeLength)
        ax.yaxis.grid(True, zorder = 1, linewidth = 0.2)
        ax.set_xlabel('Genomic position [Mb]')
        ax.set_ylabel('Outgroup sample name')
        ax.set_title(f'Polarized SNPs per outgroup in {refgenome}')
        plt.tight_layout()
        fig.savefig(f'pdf/{plot_fileprefix}_{refgenome}_outgroup_snp_pol.pdf')
        plt.close('all')


##############################
