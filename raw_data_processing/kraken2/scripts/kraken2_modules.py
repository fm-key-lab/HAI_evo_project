## This file contains all necessary python functions for running the KRAKEN2 snakemake
## Changes to any of these functions should only be done on the dev branch and thoroughly tested prior to merging

## dependencies: read_samples
import os
import subprocess
import csv
import gzip
import glob
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
import urllib.request


################################################
## functions used to identify input/output files for rules 

def get_ST_species():
    ## get always newest state of MLST DB to type all possible species
    mlst_db_url = "https://pubmlst.org/data/dbases.xml"
    response = urllib.request.urlopen(mlst_db_url).read()
    tree = ET.fromstring(response)
    st_species = []

    for species in tree.iter('species'):
        st_species.append(species.text.strip('\n').replace(' ', '_'))
    st_species_dict = dict(zip(st_species, [st_sp.split('#')[0] for st_sp in st_species]))
    return st_species_dict


################################################
## functions used in snakemake rules
def read_samplesCSV(spls):
    # reads in samples.csv file, format: Path,Sample,ReferenceGenome,ProviderName,Subject 
    header_check = ['Path','Sample','ReferenceGenome','ProviderName','Subject']
    switch = True
    file = open(spls, 'r')
    list_path = []
    list_splID = []
    list_providerNames = []
    list_refG = []
    list_subject = []
    for line in file:
        line = line.strip('\n').split(',')
        # Test Header. Note: Even when header wrong code continues (w/ warning), but first line not read.
        if switch:
            if (line == header_check):
                print("Passed CSV header check")
            else:
                Warning("CSV did NOT pass header check! Code continues, but first line ignored")
            switch = False
            continue
        # build lists
        list_path.append(line[0])
        list_splID.append(line[1])
        list_refG.append(line[2])
        list_providerNames.append(line[3])
        list_subject.append(line[4])
    return [list_path,list_splID,list_refG,list_providerNames,list_subject]

def split_samplesCSV(PATH_ls,SAMPLE_ls,REF_Genome_ls,PROVIDER_ls,PATIENTID_ls_all):
    # Added by Arolyn, 2019.02.11
    # takes info extracted form samples.csv; saves each line of samples.csv as sample_info.csv in data/{sampleID}
    for i, sample in enumerate(SAMPLE_ls):
        # get info for this sample
        sample_info_csv_text = PATH_ls[i] + ',' + SAMPLE_ls[i] + ',' + REF_Genome_ls[i] + ',' + PROVIDER_ls[i] + ',' + PATIENTID_ls_all[i]
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


## functions for taxoomic rank parsing 
def map_dirname_to_species(ref_genome_dir, ref_genome):
    ## get current reference genomes 
    ## get just species names from unprocessed ref_genomes 
    ## if file does already exists, no need for extraction of full species name
    if not os.path.isfile(ref_genome_dir + ref_genome + '/taxonomic_rank.txt'):
        if os.path.isfile(ref_genome_dir + ref_genome + '/genome.fasta'):
            with open(ref_genome_dir + ref_genome + '/genome.fasta', 'r') as fid:
                header = fid.readline().split(' ')
                ref_genome_spec = header[1] + ' ' + header[2] ## extract genus (second entry) and species (third entry) and save to ref_genome folder name (dict.key())
                return ref_genome_spec
        else:
            print(f'WARNING! genome.fasta for {ref_genome} is missing which - if that species is used in the analysis will crash the run!')
            print('Else, the pipeline will still run through!')
    

def taxrank_parser(ref_genome_dir, ncbi_tax_rank, ref_genome):

    ## convert set of ref_genome_dir to string (snakemake specific)
    if type(ref_genome_dir) == set:
        ref_genome_dir = list(ref_genome_dir)[0]
    if type(ncbi_tax_rank) == set:
        ncbi_tax_rank = list(ncbi_tax_rank)[0]

    ## get dirnames to species mapping:
    if ref_genome_dir[-1] != '/':
        ref_genome_dir = ref_genome_dir + '/'

    ## check if ref_genome is set --> get first value 
    if isinstance(ref_genome, set):
        ref_genome = str(list(ref_genome)[0])
    elif isinstance(ref_genome, list):
        ref_genome = str(ref_genome[0])

    ref_genome_spec = map_dirname_to_species(ref_genome_dir, ref_genome)
    if ref_genome_spec == None: ## species was not found or nothing to be done!
        return
    ## parse rankedlineage.dmp file, downloaded from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/
    ## and get taxonomic ranks of species of interest
    header = ['ref_genome_dir', 'species','genus','family','order','class','phylum','kingdom','domain']
    with gzip.open(ncbi_tax_rank, 'rt') as f:
        for line in f:
            l = line.strip().replace('\t', '').split('|')
            if (l[2] == '') and (l[3] != '') and (l[1] in ref_genome_spec): ## line is species entry (we do not want to resolve supspecies here!)
                ## search entry in dict to append to
                if ref_genome_spec in l[1]: ## account for different expansions of species which might be seen? If we need to refine this part 
                    line_out = np.delete(l, (0, 2, -1)) ## remove empty entry
                    line_out = np.concatenate(([ref_genome], line_out)) ## add ref_genome folder name
                    print(ref_genome)
                    print(ref_genome_dir + ref_genome + '/taxonomic_rank.txt')
                    if not os.access(ref_genome_dir + ref_genome, os.W_OK):
                        print(f'No write permissions for {ref_genome_dir + ref_genome}. Script will not break but also not write taxonomic_rank.txt file!')
                    else:
                        ## check if file has been written to already
                        if not os.path.isfile(ref_genome_dir + ref_genome + '/taxonomic_rank.txt'):
                            with open(ref_genome_dir + ref_genome + '/taxonomic_rank.txt', 'w') as outfile:
                                ### Add header to file
                                writer = csv.writer(outfile)
                                writer.writerows([header, line_out])
                        else: ## Write to outfile without adding header
                            with open(ref_genome_dir + ref_genome + '/taxonomic_rank.txt', 'a') as outfile:
                                ## Write header at first
                                writer = csv.writer(outfile)
                                writer.writerows([line_out])
