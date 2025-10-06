## This file contains all necessary python functions for running the CASE snakemake
## Changes to any of these functions should only be done on the dev branch and thoroughly tested prior to merging

## dependencies: genomestats
from Bio import SeqIO
import glob
import numpy as np

## dependencies: generate_positions_snakemake
import os

## dependencies: generate_positions_single_sample_snakemake
import gzip


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
    # turn to np.arrys!
    ChrStarts = np.asarray(ChrStarts,dtype=int)
    Genomelength = np.asarray(Genomelength,dtype=int)
    ScafNames = np.asarray(ScafNames,dtype=object)
    return [ChrStarts,Genomelength,ScafNames]

def read_samples_caseCSV(spls):
	# reads in samples_case.csv file, format: Path,Sample,ReferenceGenome,Outgroup
	hdr_check = ['Path','Sample','ReferenceGenome','Outgroup']
	switch = "on"
	file = open(spls, 'r')
	list_path = []
	list_splID = []
	list_refG = []
	list_outgroup = []
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
		list_outgroup.append(line[3])
	return [list_path,list_splID,list_refG,list_outgroup]

def unique_key_generator(key_ls, value_ls):
    ## function to allow processing of multiple genomes
	## note key_ls and value_ls must have same lengths
	## for each unique key, the values are stored in a dict entry
	d = {key_val: [] for key_val in set(key_ls)}
	for i, key in enumerate(key_ls):
		d[key].append(value_ls[i])
	return d


def chrpos2index(chrpos,chr_starts,genomelength=np.nan):
    '''
    Args:
        chrpos (arr): 2D array of chromsome ID (1-based) and position on chromosome (0-based).
        chr_starts (arr): Vector of chromosome starts (begins at 0).
        genomelength (float, optional): Total genome length. If specified, checks if any positions exceed this length. Default is np.nan.

    Returns:
        p (arr): Vector of position indexes. (0-based)

    '''
    if (np.size(chrpos,1) != 2) & (np.size(chrpos,0) == 2):
        chrpos=chrpos.T
        print('Reversed orientation of chrpos')
    elif (np.size(chrpos,1) != 2) & (np.size(chrpos,0) != 2):
        print('Wrong input format given. Provide np.array with shape (p, 2)')
        return 
    
    if len(chr_starts) == 1:
        p=chrpos[:,1]
    else:
        p=chr_starts[chrpos[:,0]-1]+chrpos[:,1]
    if not np.isnan(genomelength) and any(p>=genomelength):
        print('Invalid genome positions given with positions larger than genome length')
        return sys.exit()

    return p

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

def generate_positions_snakemake(paths_to_input_p_files, REF_GENOME_DIRECTORY):
    ## Get positions on reference genome
    [ChrStarts,GenomeLength] = genomestats(REF_GENOME_DIRECTORY)[0:2]
    ## Get variant positions
    timesvariant = np.zeros(GenomeLength, dtype=int) # initialize vector to count occurrances of variants across samples
    ## Load files 
    for path in paths_to_input_p_files:
        pos = np.load(os.getcwd() + '/' + path)
        if np.size(pos['Positions']) >= 2: # HC 9/13/2013; MF 09.02.23: need to be >= to check if at least one entry (which consists of two values)! 
            x = chrpos2index(pos['Positions'],ChrStarts)
            timesvariant[x] = timesvariant[x] + 1   
    # p = np.where((timesvariant > 0) & (timesvariant < len(paths_to_input_p_files)))[0] ## get first entry from tuple as there is just one array   # NOTE:old code, can be removed
    p = np.where(timesvariant > 0)[0] ## get first entry from tuple as there is just one array; MF 10.02.23: excluded filtering for fixed mutations (those could be multiallelic and there is no filtering based on nt information here!)
    print(f'Not considering {sum(timesvariant == len(paths_to_input_p_files))} positions where all samples have a variant compared to the reference...\n')
    return p

def generate_positions_single_sample_snakemake_withoutbool(sample_path_to_variant_vcf, sample_path_to_variant_positions, maxFQ, REF_GENOME_DIRECTORY, outgroup_boolean):
    ## ensure type of variable
    maxFQ = int(maxFQ)
    outgroup_boolean = int(outgroup_boolean)
    ## For debugging
    print(f'\nCurrently examining the following vcf file: {sample_path_to_variant_vcf}\n')
    print(f'\nFQ threshold: {maxFQ}\n')
    ## Get reference genome information
    [ChrStarts,GenomeLength,ScafNames] = genomestats(REF_GENOME_DIRECTORY)
    # Initialize boolean vector for positions to include as candidate SNPs that vary from the reference genome
    include = np.zeros(GenomeLength, dtype=int)
    ## For outgroup samples only:
    if outgroup_boolean:
        # Then save no positions
        Positions = p2chrpos(np.where(include), ChrStarts)  ## Positions=p2chrpos(find(include),ChrStarts); find function does not make sense to me?!?
        np.savez_compressed(os.getcwd() + '/' + sample_path_to_variant_positions, Positions = Positions)
        print("Sample is in outgroup")
        return # end script, no need to look at vcf file
    ## Get candidate SNP positions from variant vcf file
    # Input file
    fname_in_gz = os.getcwd() + '/' + sample_path_to_variant_vcf ## store zipped file name
    with gzip.open(fname_in_gz, 'rt') as fid:
        for line in fid:
            if line[0] != "#":
                l = line.split('\t')
                position_on_chr = int(l[1])
                position = ChrStarts[l[0] == ScafNames] + position_on_chr 
                alt = l[4]
                ref = l[3]
                #only consider for simple calls (not indel, not ambigious)
                if (alt != '') & (',' not in alt) & (len(alt) == len(ref)) & (len(ref) == 1):
                    ## find and parse quality score
                    info_col = l[7].split(';')
                    ## find FQ entry in info column and store value after '=' sign
                    entrywithFQ = [item for item in info_col if item.startswith('FQ')]
                    ## Extract FQ value from string
                    entrywithFQ_str = entrywithFQ[0]
                    fq = float(entrywithFQ_str[entrywithFQ_str.index('=')+1:])
                    if fq < maxFQ: #better than maxFQ??
                        include[position] = 1 
    Positions = p2chrpos(np.where(include)[0], ChrStarts)
    np.savez_compressed(os.getcwd() + '/' + sample_path_to_variant_positions, Positions = Positions) 

def combine_positions_snakemake_out(path_to_list_of_input_p_files, path_to_other_p_file, path_to_output_p_file, path_to_outgroup_boolean_file, REF_GENOME_DIRECTORY, looseFQmax):
    ## in_outgroup: booleans for whether or not each sample is in the outgroup
    print('Processing outgroup booleans...\n')
    # Input is space separated text file??
    fname = os.getcwd() + '/' + path_to_outgroup_boolean_file
    with open(fname, 'r') as fid:
        in_outgroup_string = fid.readline().strip() ## read just one line       
        in_outgroup_cell = in_outgroup_string.split(' ')   
    in_outgroup = np.array(in_outgroup_cell, dtype = int)
    print('Outgroup booleans:') # Print for troubleshooting
    print(in_outgroup) # Print for troubleshooting
    ## Get positions on reference genome
    ## ChrStarts = genomestats(REF_GENOME_DIRECTORY)[0]
    ## 1. Find positions with at least one fixed mutation relative to the reference genome
    print('\n\nFinding positions with at least 1 fixed mutation...\n')
    # Import list of directories for where to find variant positions for each
    # sample; turn this into a cell array
    fname = os.getcwd() + '/' + path_to_list_of_input_p_files
    with open(fname, 'r') as fid:
        fstring = fid.readline().strip() # from echo, just one line
        paths_to_input_p_files = fstring.split(' ') # currently printing out
    print('Paths used in generate positions:') # Print for troubleshooting
    paths_to_input_p_files_ingroup = np.delete(paths_to_input_p_files, np.where(in_outgroup != 0)[0])
    print(paths_to_input_p_files_ingroup) # Print for troubleshooting
    cp = generate_positions_snakemake(paths_to_input_p_files_ingroup, REF_GENOME_DIRECTORY) # Changed 2019.02.21
    #cp = generate_positions_snakemake( paths_to_input_p_files, REF_GENOME_DIRECTORY ); # Changed 2019.02.21
    print(f'Found {len(cp)} positions where provided vcf called a fixed variant in at least one in-group sample with FQ score < {looseFQmax}\n')
    ##TODO
    ## 2. Find positions with within-sample polymorphisms
    #fprintf(1,'\nFinding single nucleotide positions with within-sample polymorphisms...\n');
    # PLACEHOLDER
    print('\nWARNING! Finding within-sample polymorphisms has NOT been implemented in the snakemake pipeline!!!!!\n')
    dp = np.array([], dtype = int)
    #[dp, coveragethresholds] = find_diverse_positions_no_control(loose_parameters, {SampleDirs{~in_outgroup}}, {SampleNames{~in_outgroup}}', RUN_ON_COMPUTING_CLUSTER, jobsubmitoptions_short,TEMPORARYFOLDER,SCRIPTSDIRECTORY);
    #fprintf(1,'Found #i positions with within-sample polymorphism that meets loose parameters in at least 1 in-group sample \n',length(dp)) ;
    ## 3. Add candidate positions manually
    #fprintf(1,'\nAdding other positions to consider manually...\n');
    if os.path.exists(path_to_other_p_file):
        other_p_file = np.load(path_to_other_p_file)
        if 'op' not in other_p_file:
            SystemError('ERROR! File for other positions (' + path_to_other_p_file + ') does not contain variable op!')
        else:
            op = other_p_file['op']
            SystemError('Please check the output of op before going further on! This was not checked and compared to the matlab output as there was no such file up to now!')
    else:
        op = np.array([], dtype=int) # empty vector if no positions specified manually
    # possibly change to if statement with if nargin < N...     
    print(f'\nConsidering {len(op)} positions previously specified \n')
    ## Combine all three types of positions
    print('\nCombining all three types of positions into one list...\n') # Aro_Change: print this
    allp = np.unique(np.concatenate([dp, op, cp], axis = 0)) 
    p = np.sort(allp)
    p = p[p > 0] #remove any 0s
    #positions=p2chrpos(p,ChrStarts); # Aro_Change: not clear where this is used??
    ## Save positions
    print('\nSaving list of all positions...\n') # Aro_Change: print this
    print(f'\nSaving: {os.getcwd()}/{path_to_output_p_file}\n') # Aro_Change: print this
    ## store in matlab shape:
    p = np.atleast_2d(p).T
    np.savez_compressed(os.getcwd() + '/' + path_to_output_p_file, p=p)

