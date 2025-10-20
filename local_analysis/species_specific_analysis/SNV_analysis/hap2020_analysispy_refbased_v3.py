# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Started 2022/05/02

refbased analysis with SM runs from 2024/08

@author: mf
"""

# %%
# =============================================================================
# SNP-based analyses:
# =============================================================================
# - parsimony tree
# - SNP-specific tree coloring
# - barplots fwd/rev for SNPs
# - parallel evolution analysis
# - molecular clock analysis
# - read in unified ortholog table (creation NOT yet part of assembly snakemake)
# - allele frequency change analysis
# - read-in fucntions from apy (analysis.py) module

## To ADD:
# - tree branch coloring according to majority visit

# %% apy module (contains all functions)

import sys,os,re,subprocess
import glob,pickle,time,datetime,warnings
import numpy as np
from scipy import stats
import pandas as pd
import seaborn as sns

from pylab import * #what is that

from matplotlib import rc
# enable plotting in separate windows
%matplotlib qt

## Directory to analysispy_module.py
SCRIPTS_DIRECTORY = os.getcwd() + "/../../modules/"
sys.path.insert(0, SCRIPTS_DIRECTORY)

#import analysispy_module_2022 as apy
import analysispy_module as apy

# %% Filter Parameter -- IMPORTANT TO ADJUST !!!!
# Much of this will need to vary according to your reference genome,
# coverage, and particular samples
# consider to add SD filter in addition to coverage. Isolates of other species that harbor plasmid present in pangenome at very high cov can otherwise meet the coverage threshold even though FP

# Remove samples that are not high quality

filter_parameter_sample_across_sites = {\
                                        'min_median_coverage_to_include_sample': 8, 
                                        'min_basecalls_to_include_sample': 0.1, # remove samples that have too many undefined base (ie. N). added this filter.
                                        }

# Remove sites within samples that are not high quality
filter_parameter_site_per_sample = {\
                                    'min_maf_for_call' : 0.95, #on individual samples, calls
                                    'min_cov_per_strand_for_call' : 2,  # on individual samples, calls
                                    'min_qual_for_call' : 30,  #on individual samples, calls
                                    }

# Remove sites across samples that are not high quality or indications for recombination event
filter_parameter_site_across_samples = {\
                                        'max_fraction_ambigious_samples' : 0.01, #across samples per position
                                        'min_median_coverage_position' : 4, #across samples per position
                                        'distance_for_nonsnp' : 500, #region in bp on either side of goodpos that is considered for recombination
                                        'corr_threshold_recombination' : 0.75 #minimum threshold for correlation to remove site 
                                        }

## how far upstream of the nearest gene to annotate something a promoter
## mutation (not used if no annotation)
promotersize=250;

reuse_latest_cutoff_vals = False
on_the_fly_eval = False ## reduces max_fraction_ambigious_samples for samples with extensive SNPs to lower fraction (> 1000 SNPs)
## if old values shall be taken, all old shall be used and this need to be hard coded to set off!
if reuse_latest_cutoff_vals == True: 
    on_the_fly_eval = False

# %% set samples which do not place well in tree
samplesNames_offSamples = {}

# %%  cd(workingdir)
analysis_run = 'refbased_2024_08'
#analysis_run = 'refbased_2022_04'
#analysis_run = 'denovo_2022_06'

main_workingdir = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/' + analysis_run + '/analysis/')
os.chdir(os.path.expanduser(main_workingdir))

#cand_mut_table_dirs = '/../case/2-candidate_mutation_table_2022_06_11'
cand_mut_table_dirs = '/../case/2-candidate_mutation_table_2024_08_23'

## get species and their reference genomes (stored in subdir name) + 
## remove hidden files from list
refgenomes_ls = [dir for dir in os.listdir(os.getcwd() + cand_mut_table_dirs) if not dir.startswith('.')]


# %% Other variables
# TIMESTAMP
# used in output files
ts = time.time()
timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S') # WARN: fixed timestamp not yet build in

## Nucleotides
NTs = np.array(['A','T','C','G'],dtype=object) # NTs='ATCG'

# %% Load in specimen information
# parse specimen info to csv
## grep newest specimen logfile (done via date in filename in order YYYY_MM_DD_specimenlog.csv)
specimenLog_path = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/Drylab/metadata/*_specimenlog.csv')
specimenLog_file = sorted(glob.glob(specimenLog_path), key = lambda x: x.rsplit('/')[-1])[-1]

SpecimenLog = {}
SpecimenLog["Patient"] = [];SpecimenLog["date"] = [];SpecimenLog["Kit"] = [];SpecimenLog['sampleids'] = [];SpecimenLog["timepoint"] = [];SpecimenLog["plate"] = [];
with open(specimenLog_file,"r") as file_specimen:
    for line in file_specimen:
        if not line.startswith("patient"):
            line = line.strip().split(",")
            SpecimenLog["Patient"].append(line[0])
            SpecimenLog["date"].append(line[1])
            SpecimenLog["Kit"].append(line[2])
            SpecimenLog["sampleids"].append(line[3].replace('-','_'))
            SpecimenLog["timepoint"].append(line[5])
            SpecimenLog["plate"].append(line[6])

# %% make dict of species as keys and patients in which those species have been found
## Read in samples csv fot that
samples_csv = glob.glob(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/' + analysis_run + '/mapping/*samples.csv'))
samples_csv = sorted(samples_csv)[-1] ## return newest samples_case.csv

species_in_patients_dict = {}
with open(samples_csv, 'r') as fo:
    for line in fo: 
        if not 'ReferenceGenome' in line: ## skip first line 
            line = line.strip().split(',')
            ## include up to 3 leading 0 (5 chars for patient ID each)
            if line[4] != 'Key':
                pat_id = line[4][0] + str(line[4][1:]).zfill(4)

                ## samples csv is ordered by Path,Sample,ReferenceGenome,ProviderName,Subject
                for refgenome in line[2].split(' '): ## account for multiple ref genomes stated in samples csv:
                    if not refgenome in species_in_patients_dict.keys(): ## if reference genome not in keys, then generate new key
                        species_in_patients_dict[refgenome] = [pat_id]
                    elif pat_id not in species_in_patients_dict[refgenome]: ## if patient id not in values of key, then store 
                        species_in_patients_dict[refgenome].append(pat_id)
species_in_patients_dict = dict(sorted(species_in_patients_dict.items()))

## list of samples identified to be contaminated 
contaminated_samples = ['P07_3453', 'P07_2730', 'L00071_07_1', 'L00071_02_1', 'P21_0043', 'P21_MIBI03']

# %% Define variables to store subject-specific results 
para_evo_cand = np.array([],dtype=object)
mut_genes_df = {}
snp_freq_shifts_anno_lod = []
regression_res_lod = [] # stores all the results of the regression, and turned to pandas DF and saved
annotation_mutation_allParaSignal = {} # store all annotation_mutation for all detected candidates
mol_clock_data_dc = {} # dc to store for each subject the calculated data for mol clock estimation
allP_gene_paraevo = {} # store all para evo candidates that fullfill bonferroni

# %% Initiate loop over all subjects
# subject_fld_name = ["P0021"]
# refgenome = refgenomes_ls[-1]
for refgenome, subject_fld_name in species_in_patients_dict.items():
    
    for subj_idx,subject_fld_label in enumerate(subject_fld_name):
        print('\n'+'Process ' + subject_fld_label + ' and ' + refgenome)

        dir_timestamp = timestamp.split('_')[0].replace('-', '_')
        ## make per species a separate folder
        if not os.path.exists(main_workingdir + subject_fld_label + '_' + refgenome):
            ## generate dir 
            os.mkdir(main_workingdir + subject_fld_label + '_' + refgenome)
        if not os.path.exists(main_workingdir + subject_fld_label + '_' + refgenome + '/' + dir_timestamp + "_" + refgenome):
            ## generate dir 
            os.mkdir(main_workingdir + subject_fld_label + '_' + refgenome + '/' + dir_timestamp + "_" + refgenome)
            os.mkdir(main_workingdir + subject_fld_label + '_'  + refgenome + '/' + dir_timestamp + "_" + refgenome + '/pdf')
        os.chdir(main_workingdir + subject_fld_label + '_' + refgenome + '/' + dir_timestamp + "_" + refgenome)

        if 'refbased' in analysis_run:
            ref_genome_folder = os.path.expanduser('~/Nextcloud/keylab/reference_genomes/') + refgenome
        else:
            ref_genome_folder = os.path.expanduser(main_workingdir + '../assembly/assembled_genomes/') + refgenome

        # %% 
        # =============================================================================
        # Load data from candidate_mutation_
        
        # =============================================================================
        # load('candidate_mutation_table')
        # 
        # Import candidate_mutation_table
        cmtdir = main_workingdir + cand_mut_table_dirs + '/' + refgenome + '/'
        cmtFile = 'candidate_mutation_table.pickle.gz'
        [quals,p,counts,in_outgroup,sampleNames,indel_counter,coverage_stats] = apy.read_candidate_mutation_table_pickle_gzip(cmtdir + cmtFile)


        # %%
        # =============================================================================

        # Before doing any data parsing, save cutoff values to csv or read in latest cutoff vals if wanted
        if reuse_latest_cutoff_vals:
            try:
                ## load old data
                [filter_parameter_sample_across_sites, filter_parameter_site_per_sample, filter_parameter_site_across_samples, promotersize] = apy.read_latest_cutoff_vals(main_workingdir + subject_fld_label + '_' + refgenome)
                print('Old cut-offs are used again!')
            except:
                ## if no old cut off values are stored, set the new ones
                set_filters = {'filter_parameter_sample_across_sites': filter_parameter_sample_across_sites, 
                              'filter_parameter_site_per_sample': filter_parameter_site_per_sample, 
                              'filter_parameter_site_across_samples': filter_parameter_site_across_samples,
                              'promotersize': promotersize}
                apy.save_cutoff_vals(timestamp, set_filters)
                print('Old cut-offs not found. Setting new ones!')
        else:
            if on_the_fly_eval & (len(p) > 1000):
                filter_parameter_site_across_samples['max_fraction_ambigious_samples'] = 0.01
            elif on_the_fly_eval:
                filter_parameter_site_across_samples['max_fraction_ambigious_samples'] = 0.99
            set_filters = {'filter_parameter_sample_across_sites': filter_parameter_sample_across_sites, 
                              'filter_parameter_site_per_sample': filter_parameter_site_per_sample, 
                              'filter_parameter_site_across_samples': filter_parameter_site_across_samples,
                              'promotersize': promotersize}
            apy.save_cutoff_vals(timestamp, set_filters)
            print('Setting new cut-off values!')
        
        ## extract only samples (= entries) which are of the patient(s) of interest 
        subject_smpls_speclog = [sampleid for subject, sampleid in zip(SpecimenLog['Patient'], SpecimenLog['sampleids']) if (subject == subject_fld_label)] ## extract which sampleids are derived from current patient
              
        ## rename samples in sampleNames if patient sample and it only has a 3 digits ID
        sampleNames = np.array([smplname.split('_')[0] + '_' + smplname.split('_')[1].zfill(4) if ( (smplname[0] == 'P') & (smplname.split('_')[-1].isdigit()) ) else smplname for smplname in sampleNames])

        ## subset matrices of cand_mutation_table:
        subj_specific_samples = [idx for idx, sampleid in enumerate(sampleNames) if (((sampleid in subject_smpls_speclog) & (sampleid not in contaminated_samples)) | ('_GCA_' in sampleid) | ('_GCF_' in sampleid) | ('_ASM_' in sampleid))] ## extract indices which include samples of specific patient
        
        quals = quals[:, subj_specific_samples]
        counts = counts[subj_specific_samples, :, :]
        in_outgroup = in_outgroup[:, subj_specific_samples]
        sampleNames = sampleNames[subj_specific_samples]
        indel_counter = indel_counter[subj_specific_samples, :, :]
        coverage_stats = coverage_stats[subj_specific_samples, :]

        # %% 
        # get median cov per samples based on all positions in counts (aka p)
        median_cov_p_per_sample = coverage_stats[:, 10] ## each row is sample, col [0-10) == covg bins 1x, 2x... >10x; [10]==median covg; [11]==mean; [12]==stddev 

        # % indel counter reduction to indel count
        # indel_counter: The first statistic is the number of reads (at this position and in this sample) that 
        # support an indel. The second statistics is the number of reads (at this position and in this sample) that support a deletion.
        # for now we need only the count for indels overall 
        
        indel_counter = indel_counter[:,0,:].transpose()
        # indel counter >> 50% indel row:1; row2 >> deletion
        
        # %% Save everything using more permanent sample names to avoid overwriting
        
        sampleNames_all=np.asarray(sampleNames,dtype=object)
        quals_all=-quals;
        counts_all=counts;
        coverage_all = counts_all.sum(axis=1).transpose() # axis=1 == rows; transpose needed > rows: pos and col: samples

        # %% 
        # Assign patient & visit number to each sample; Parse specimen name to human readable
        patients_all = np.empty(len(sampleNames_all), dtype=object) #np.zeros(shape=(1, len(sampleNames_all)), dtype=np.int) # row w/ 0's 
        timepoints_all = np.zeros(len(sampleNames_all), dtype=np.int64) # row w/ 0's
        locations_all = np.zeros(len(sampleNames_all), dtype=np.int64) # row w/ 0's
        kit_date_all = np.zeros(len(sampleNames_all), dtype='object') # row w/ 0's
        original_plate_all = np.zeros(len(sampleNames_all), dtype='object') # row w/ 0's

        locations_abbreviations = ['N','O','R','S','B','L','G','U']; 
        smpl_type_abbreviations = ['R', 'R', 'R', 'R', 'D', 'D', 'D', 'D']; ## routine = R, diagnostic = D
        locations_long_names = ['Nasal', 'Oral', 'Rectal', 'Skin', 'Blood', 'Lung', 'Gastric', 'Urine'];

        for i in range(0,len(sampleNames_all)):
            if ('_GCA_' in sampleNames_all[i]) | ('_GCF_' in sampleNames_all[i]) | ('_ASM_' in sampleNames_all[i]):
                continue
            if sampleNames_all[i].endswith('_S2') or sampleNames_all[i].endswith('_S3'):
                kit_idx = SpecimenLog['sampleids'].index( sampleNames_all[i][:-3] )
            else:
                kit_idx = SpecimenLog['sampleids'].index( sampleNames_all[i] )
            patients_all[i] =  SpecimenLog['Patient'][kit_idx]
            timepoints_all[i] =  SpecimenLog['timepoint'][kit_idx]
            kit_date_all[i] =  SpecimenLog['date'][kit_idx].replace('-', '')[2:8] ## concat YYYYMMDD and remove time (HH:MM)
            original_plate_all[i] = SpecimenLog['plate'][kit_idx]
            if SpecimenLog['Kit'][kit_idx][0] == 'K':
                currID = SpecimenLog['Kit'][kit_idx].split('-')[-1] # replace everything but location identifier
            else:
                currID = SpecimenLog['Kit'][kit_idx][0] ## Diagnostic samples
            for l in locations_abbreviations:
                if currID in l:
                    locationmatch = locations_abbreviations.index(l) # get location of sublist that contains currID
            locations_all[i] = locationmatch
        
        patients_sampled = np.unique(patients_all[in_outgroup[0] == '0'])
        
        ## generate dict with days since first isolate and tp mapping 
        kit_date_ingroup_datetime = [datetime.datetime.strptime(date, '%y%m%d').date() for date in kit_date_all[in_outgroup[0] == '0']] 
        first_occurence = min(kit_date_ingroup_datetime)
        timepoint_to_colonization_day_dict = {tp: (sorted(kit_date_ingroup_datetime)[idx] - first_occurence).days for idx, tp in enumerate(sorted(timepoints_all[in_outgroup[0] == '0']))} ## NOTE: extraction is done on sorted lists for TP and dates to have sorted dict output


        # %%
        # =============================================================================
        #     Read in genome information
        # =============================================================================
        [chrStarts, genomeLength, scafNames] = apy.genomestats(ref_genome_folder);

        
        # %% Display samples the NOT fullfill min_median_coverage_to_include_sample: 
        lowcovsamples = sampleNames_all[ median_cov_p_per_sample < filter_parameter_sample_across_sites['min_median_coverage_to_include_sample'] ]
        print(lowcovsamples)
        np.savetxt('sampleNames_lowCovSamples.txt',lowcovsamples,fmt='%s')

        # %% 
        # =============================================================================
        #     Remove undesired samples based on genome-wide coverage
        # =============================================================================

        # %% Define goodsamples and filter data

        goodsamples =  (median_cov_p_per_sample >= filter_parameter_sample_across_sites['min_median_coverage_to_include_sample'])

        sampleNames = sampleNames_all[goodsamples]
        counts = counts_all[goodsamples , : , : ] # keep only level (samples) that fullfil filter!
        quals = quals_all[ : , goodsamples ]
        coverage = coverage_all[ : ,goodsamples]
        patients = patients_all[goodsamples]
        timepoints = timepoints_all[goodsamples]
        kit_date = kit_date_all[goodsamples]
        locations = locations_all[goodsamples]
        original_plates = original_plate_all[goodsamples]
        indels = indel_counter[:,goodsamples]
        
        num_samples = len(sampleNames)
        
        coverage_forward_strand = counts[:,0:4,:].sum(axis=1).transpose()
        coverage_reverse_strand = counts[:,4:8,:].sum(axis=1).transpose()

        # %% Breakpoint: Too few samples passed filter
        if np.sum(goodsamples) < 2:
            print("Too few samples fullfill filter criteria! >> skip: " + refgenome)
            continue

        # %%
        # =============================================================================
        # Extract refnt and define out/in-group bools
        # =============================================================================
        ## Note ancnti/outs_nti defined below after filtered calls has been generated!
        ## get reference allele for all p; NOTE: analysis.m stored ATCG as double-digit-numeric
        # use ref nt for mutation calls. important if multiple outgroups called 
        refnt = apy.extract_outgroup_mutation_positions(ref_genome_folder, apy.p2chrpos(p,chrStarts));
        refnti = apy.nts2idx(refnt)
        refnti_m = np.tile(refnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele

        ## Estimate outgroup (ancestral allele) from ALL samples added as outgroup 
        pattern = re.compile('_GCA_|_GCF_|_ASM_') # string tag in sample Name to identify outgroup in staphAD data
        outgroup_name = np.array(list(filter(pattern.search, list(sampleNames))))
        outgroup_bool = np.isin(sampleNames , outgroup_name)
        
        # ingroup array (bool, idx) used later
        ingroup_bool = np.invert(outgroup_bool)
        ingroup_idx = np.nonzero(ingroup_bool)[0]

        ## bool for all single outgroup sample AND ingroup-samples with dedicated single outgroup sample
        outgroup_currSub = np.array(["GCA", 'GCF', 'ASM'])
        outgroup_spl_idx = [i for i,item in enumerate(sampleNames) if (outgroup_currSub[0] in item) | (outgroup_currSub[1] in item) | (outgroup_currSub[2] in item)]
        ingroup_wOutSpl_bool = np.copy(ingroup_bool)
        ingroup_wOutSpl_bool[outgroup_spl_idx] = True
        ingroup_wOutSpl_onlyIngroup_bool = np.array([False if ('GCA' in s) | ('GCF' in s) | ('ASM' in s) else True for s in sampleNames[ingroup_wOutSpl_bool]])

        # %%
        # =============================================================================
        #     Extract allele and frequencies
        # =============================================================================
        contig_positions = apy.p2chrpos(p,chrStarts) # 1col: chr, 2col: pos on chr; for all p
        [maf, maNT, minorNT, minorAF] = apy.div_major_allele_freq(counts)     #generates data structures of allele frequencies
        minorAF[np.isnan(minorAF)]=0     #set nan values to zero
        # NOTE: function assumes first 8 rows in counts == 4nucl fwd&rev! watch out if extended counts used!
        # NOTE: maf==0 -> no data;minorAF==0 -> no minor allele/ or no major allele; NT number corresponds to index in NTs [ATCG] or if maf==0 > NA == 4  

        # %%
        # =============================================================================
        # Make some basic structures for finding mutations
        # =============================================================================
        #create mutantAF
        mutantAF=np.zeros((maNT.shape))
        mutantAF[np.where(maNT!=refnti_m)] = maf[np.where(maNT!=refnti_m)]

        # %%
        # =============================================================================
        #     Find positions with fixed mutations
        # =============================================================================
        # Define mutations we do not trust in each and across samples.
        # goodpos are indices of p that we trust
        
        ## Filter per mutation
        # added indel filter! (less than 50% of reads had site are allowed to support indel)
        def_base_samples = np.sum(maNT < 4, axis = 0) ## extract how many SNPs per samples are not N 
        ## print how many SNPs at least, at max and at median have been removed by filter --> control to see if some filter are too strict/loose
        ## Note: this must be calculated before maNT is changed by subsetting!
        snps_Ns = np.sum(maNT > 3, axis = 0) / np.sum(maNT <= 4, axis = 0) * 100
        snps_removed_by_quals = np.sum(((quals < filter_parameter_site_per_sample['min_qual_for_call']) & (maNT < 4)), axis = 0) / def_base_samples * 100
        snps_removed_by_maf = np.sum(((maf < filter_parameter_site_per_sample['min_maf_for_call']) & (maNT < 4)), axis = 0) / def_base_samples * 100
        snps_removed_by_fwdstrand = np.sum(((coverage_forward_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call']) & (maNT < 4)), axis = 0) / def_base_samples * 100
        snps_removed_by_revstrand = np.sum(((coverage_reverse_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call']) & (maNT < 4)), axis = 0) / def_base_samples * 100
        snps_removed_by_indel = np.sum(((indels > (0.5*coverage)) & (maNT < 4)), axis = 0) / def_base_samples * 100
        
        calls = np.copy(maNT) # NOTE: reassigning calls to 4 changes also maNT! --> with copy statment not any more!
        calls[ (quals < filter_parameter_site_per_sample['min_qual_for_call']) ] = 4
        calls[ (maf < filter_parameter_site_per_sample['min_maf_for_call']) ] = 4
        calls[ (coverage_forward_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call']) ] = 4
        calls[ (coverage_reverse_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call']) ] = 4
        calls[ (indels > (0.5*coverage) ) ] = 4 

        print(f'Rate of undefined bases (Ns) in samples: \tMinimum {min(snps_Ns):.2f}%, \tmaximum {max(snps_Ns):.2f}% and at \tmedian {median(snps_Ns):.2f}% of SNPs per sample')
        print(f'Filtering by FQ scores: \t\tRemoved minimum {min(snps_removed_by_quals):.2f}%, \tmaximum {max(snps_removed_by_quals):.2f}% and at \tmedian {median(snps_removed_by_quals):.2f}% of SNPs per sample')
        print(f'Filtering by major allele frequency: \tRemoved minimum {min(snps_removed_by_maf):.2f}%, \tmaximum {max(snps_removed_by_maf):.2f}% and at \tmedian {median(snps_removed_by_maf):.2f}% of SNPs per sample')
        print(f'Filtering by forward strand coverage: \tRemoved minimum {min(snps_removed_by_fwdstrand):.2f}%, \tmaximum {max(snps_removed_by_fwdstrand):.2f}% and at \tmedian {median(snps_removed_by_fwdstrand):.2f}% of SNPs per sample')
        print(f'Filtering by reverse strand coverage: \tRemoved minimum {min(snps_removed_by_revstrand):.2f}%, \tmaximum {max(snps_removed_by_revstrand):.2f}% and at \tmedian {median(snps_removed_by_revstrand):.2f}% of SNPs per sample')
        print(f'Filtering by indel support: \t\tRemoved minimum {min(snps_removed_by_indel):.2f}%, \tmaximum {max(snps_removed_by_indel):.2f}% and at \tmedian {median(snps_removed_by_indel):.2f}% of SNPs per sample')

        ## Filter per site across samples
        # Ignore here outgroup samples!
        siteFilt = np.any(( (calls[:,ingroup_bool]>3).sum(axis=1) >= ((num_samples-np.sum(outgroup_bool)) * filter_parameter_site_across_samples['max_fraction_ambigious_samples']) \
                            ,np.median( coverage[:,ingroup_bool], axis=1) < filter_parameter_site_across_samples['min_median_coverage_position'] ),axis=0)
        calls[ siteFilt ,:] = 4 # sites that fail qc -> 4, for all samples incl. outgroup     
            
        # NOTE: func below takes forever with many SNPs...saved below
        [mutQual, mutQualIsolates] = apy.ana_mutation_quality(calls[:,ingroup_bool],quals[:,ingroup_bool]) # get FQ value for SNP across samples. mutQualIsolates contains sample indices for sample pair FQ based on. 
        mutQual = np.nan_to_num(mutQual, nan=-1) # turn mutQual nan's to -1; necessary to avoid later warning
        
        ## translate filtered calls of ingroup into goodpos. mutations we believe. fixedmutation part removed in v6.
        hasmutation = (calls != refnti_m) & (calls < 4) & (np.tile(mutQual,(1,num_samples)) >= 1) # consider only ingroup samples; mutQual >= 1 is very loose. Important filter with low qual data! refnt not ancnt!!!
        hasmutation[:,outgroup_bool] = False # put outgroup samples 4 in order to identify ingroup mutations only
        
        candpos = np.where( np.sum(hasmutation, axis=1) > 0 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!

        print(f'Identified {len(candpos)} candidate positions for SNPs')
        if len(candpos) == 0:
            ## save filtered data to npz file 
            save_candpos_cmt_dict = {'calls': calls[candpos, :], 
                                    'mutantAF': mutantAF[candpos, :], 
                                    'p': p[candpos], 
                                    'goodpos': candpos, 
                                    'refnti': refnti[candpos],
                                    'sampleNames': sampleNames,
                                    'counts': counts[:,:, candpos],
                                    'quals': quals[candpos,:],
                                    'coverage': coverage[candpos,:],
                                    'patients': patients,
                                    'timepoints': timepoints,
                                    'kit_date': kit_date,
                                    'locations': locations,
                                    'original_plates': original_plates,
                                    'indels': indels[candpos,:]}
            np.savez_compressed(f'{timestamp}_goodpos_cmt.npz', **save_candpos_cmt_dict)
            continue

        # %%
        # =============================================================================
        #  Check for recombination in p and remove positions from goodpos
        # =============================================================================
        
        #candpos = candpos[1:10000]
        if len(candpos) > 0:
            [recombpos, recombpos_bool, evidence_of_recomb, mixed_snp_nonsnp_idx, nonrecomb_mutantAF] = apy.findrecombinantSNPs_new(p[candpos], 
                                                                                        mutantAF[ np.ix_(candpos, ingroup_bool) ],
                                                                                        filter_parameter_site_across_samples['distance_for_nonsnp'],
                                                                                        filter_parameter_site_across_samples['corr_threshold_recombination'],
                                                                                        eval_recom_p_allele = True, 
                                                                                        nts_idx = apy.nts2idx(np.array([nt for nt in 'ATCG?'])), 
                                                                                        mant = maNT[ np.ix_(candpos, ingroup_bool) ]) 
            #These are the positions in p that are likely recombinant that we will remove from goodpos
            goodpos = candpos[ np.invert(recombpos_bool) ] # NOTE: candpos/goodpos is INDEX of good positions for p! 
            ## add to norecomb_calls the outgroup samples again
            calls[np.ix_(candpos, ingroup_bool)] = np.where(evidence_of_recomb, 4, calls[np.ix_(candpos, ingroup_bool)])
        else:
            goodpos = candpos
        print('======================================================')
        print(f'{np.shape(goodpos)[0]} good positions have been identified.')
        print('======================================================\n')
        
        ## save filtered data to npz file 
        save_goodpos_cmt_dict = {'calls': calls[goodpos, :], 
                                 'mutantAF': mutantAF[goodpos, :], 
                                 'p': p[goodpos], 
                                 'goodpos': goodpos, 
                                 'refnti': refnti[goodpos],
                                 'sampleNames': sampleNames,
                                 'counts': counts[:,:, goodpos],
                                 'quals': quals[goodpos,:],
                                 'coverage': coverage[goodpos,:],
                                 'patients': patients,
                                 'timepoints': timepoints,
                                 'kit_date': kit_date,
                                 'locations': locations,
                                 'original_plates': original_plates,
                                 'indels': indels[goodpos,:]}
        np.savez_compressed(f'{timestamp}_goodpos_cmt.npz', **save_goodpos_cmt_dict)
        
        # %% 
        # Define ancnti and outsplnti based on filtered calls
        ## ancnti == all outgroup samples regardless of phylogeny
        calls_outgroup = calls[:,outgroup_bool] 
        print(f'Checking outgroup samples (n = {np.shape(calls_outgroup)[1]}) for major allele:')
        ancnti = apy.major_allele(calls_outgroup) # NOTE: the filter criteria (cov,qual etc.) are not applied before major allele call
        ancnti_m = np.tile(ancnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele    

        ## outsplnti == phylogenetically closest high-cov outgroup
        outsplsnti = calls[:,outgroup_spl_idx]
        ## for identification of good outgroup samples:
        polarized_site_file_txt = '' ## save number of polarized sites to store in file 
        for i in range(np.shape(outsplsnti)[1]): ## loop over all outgroup samples
            unpol_snps = np.sum(outsplsnti[:,i][goodpos]==4)
            outgrp_splname = sampleNames[outgroup_spl_idx][i]
            polarized_site_file_txt = '\n'.join([polarized_site_file_txt, f'{outgrp_splname}:\t\t{unpol_snps} unpolarized SNPs'])
        ## select the outgroup which is phylogenetically the closest
        min_outgroupsplnti = np.argmin(np.sum(outsplsnti[goodpos]==4, axis = 0))
        outsplnti = outsplsnti[:,min_outgroupsplnti]
        outsplnti_m = np.tile(outsplnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele    
        polarized_site_file_txt = '\n'.join([polarized_site_file_txt, f'\nUnresolved goodpos in outsplnti outgroup:\t\t{np.sum(outsplnti[goodpos]==4)}'])
        polarized_site_file_txt = '\n'.join([polarized_site_file_txt, f'Closest outgroup sample name:\t\t{sampleNames[outgroup_spl_idx][min_outgroupsplnti]}'])
                
        ## out chimera. based on outsplnti but all NA (==4) are replaced by major allele call in ancnti
        # we use this for molecular clock and derived allele freq (but not tree!)
        outchimerasplancnti = np.array([ancnti[i] if nti==4 else nti for i,nti in enumerate(outsplnti) ])
        polarized_site_file_txt = '\n'.join([polarized_site_file_txt, f'Unresolved goodpos in chimera outgroup:\t\t{np.sum(outchimerasplancnti[goodpos]==4)}'])
        polarized_site_file_txt = '\n'.join([polarized_site_file_txt, f'Resolved goodpos in chimera outgroup:\t\t{len(goodpos) - np.sum(outchimerasplancnti[goodpos]==4)} ({((len(goodpos) - np.sum(outchimerasplancnti[goodpos]==4)) / len(goodpos)) * 100}%)'])
        print(polarized_site_file_txt)

        with open(f'{timestamp}_polarized_sites_by_outgroup.txt', 'w') as fo_pol:
            fo_pol.write(polarized_site_file_txt)

        ## plot which outgroup sample is polarizing which SNP and which SNPs are not polarized at all 
        outsplsnti_chim = np.concatenate((outsplsnti, np.reshape(outchimerasplancnti, (-1, 1))), axis = 1) ## include the major (chimeric) calls
        outgrpspl_idx_sort = np.argsort(np.sum(outsplsnti_chim[goodpos]==4, axis = 0)) ## idx to sort samples based on the number of SNPs they polarize 
        outgrp_sampleNames = np.append(sampleNames[outgroup_spl_idx], 'Chimeric outgroup') ## extract sample names and add chimeric sample name
        outgrp_sampleNames = outgrp_sampleNames + '\n(unpol. SNPs: ' + (np.sum(outsplsnti_chim[goodpos]==4, axis = 0)).astype(str) + ')' ## add information how many snps not resolved by outgroup
        outgrp_sampleNames_s = outgrp_sampleNames[outgrpspl_idx_sort] ## sort sample names 
        outgrp_sampleNames_sext = np.tile(outgrp_sampleNames_s, len(goodpos)) ## repeat samplenames n times 

        gnome_pos_snps = p[goodpos] ## get genomic positions of goodpos
        gnome_pos_snps_ext = np.repeat(gnome_pos_snps, len(outgrpspl_idx_sort)) ## repeat each pos * outgroup sample number 
        
        outgrpspl_pol_pos = outsplsnti_chim[np.ix_(goodpos, outgrpspl_idx_sort)] ## get calls of outgroups and sort by idx 
        #outgrpspl_pol_pos = (outgrpspl_pol_pos != 4)
        outgrpspl_pol_pos_1d = outgrpspl_pol_pos.ravel() ## convert to 1d array for plotting
        outgrpspl_pol_pos_1d[outgrpspl_pol_pos_1d != 4] = True
        outgrpspl_pol_pos_1d[outgrpspl_pol_pos_1d == 4] = False

        ## get contig cutoffs 
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ['black', 'darkorange'])
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
        fig.savefig(f'pdf/{timestamp}_{subject_fld_label}_{refgenome}_outgroup_snp_pol.pdf')

        # %%
        # =============================================================================
        #     Annotate genome
        # =============================================================================
        # Uncomment the following if there is an annotated genome
        # NOTE: annotation is read from *.gff file in ref_folder!
        annotation_genes = apy.parse_gff(ref_genome_folder,scafNames, forceReDo=True) # ref_genome_folder+"/annotation_genes.pandas.py.pk1"

        # %%
        # =============================================================================
        #  Make table with annotations
        # =============================================================================

        order = np.arange(sampleNames.shape[0])
        num_contigs = np.max(contig_positions[:,0]);
        contig_lengths = np.append(chrStarts[1:], genomeLength) - chrStarts
        
        # Uncomment the following if there is an annotated genome
        # NOTE: annotation is read from *.gff file in ref_folder!
        # annotation_genes = apy.parse_gff(ref_genome_folder,scafNames,ortholog_df[subject_fld_label],forceReDo=False) # ref_genome_folder+"/annotation_genes.pandas.py.pk1"
        ## If no ortholog_df is there just leave it out
        annotation_mutations = apy.annotate_mutations_v2(annotation_genes , p[goodpos] , refnti_m[goodpos,:] , ancnti_m[goodpos,:], calls[goodpos,:] , counts[:,:,goodpos] , hasmutation[goodpos,:], mutQual[goodpos,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP
        mutational_priorities = ['N', 'S', 'R', 'P', 'I', 'U']
        annotation_mutations_dedup = apy.deduplicate_annotation_mutations(annotation_mutations, timestamp, mutational_priorities)
        pd.concat(annotation_genes).to_csv(timestamp + '_annotation_genes.csv')
        
        annotation_mutations.to_csv(timestamp + '_annotation_SNVs.csv')

        # %% 
        # =============================================================================
        #     Breakpoint: Too few positions pass filter
        # =============================================================================
        if len(goodpos) < 2:
            print("Too few positions after filter! >> skip: " + refgenome)
            continue

        # %% 
        # =============================================================================
        #     parsimony tree. (incl. samples filter)
        # =============================================================================    
        # define pos to use in tree
        goodpos2useTree = goodpos #[0::6] #(1:1000); %trim for easier tree view; TDL called quality_positions
        
        # get data and filter for goodpos
        calls_for_treei = calls; 
        calls_for_treei = calls_for_treei[ goodpos2useTree, : ]
        calls_for_treei = calls_for_treei[ : , ingroup_wOutSpl_bool ]
        
        seq_types = dict()
        with open(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/kraken2/5-srst2/sumSRST2_results.csv'), 'r') as fid:
            header = fid.readline()
            for line in fid:
                line = line.strip().split(',')
                seq_types[line[1]] = line[2].replace('*', 'x')
        seq_types.pop('L00071_07_01', None) ## Enterobacter sample defined as P aeruginosa by kraken 

        # build sampleNames w/ metainfo
        treesampleNamesLong = np.copy(sampleNames)
        for i,name in enumerate(sampleNames):
            if ('_GCA_' in treesampleNamesLong[i]) | ('_GCF_' in treesampleNamesLong[i]): ## need to check if subsetting in for i, name in enumerate(sampleNames[ingroup_bool]) would be successful as well! 
                treesampleNamesLong[i] = 'OG_GC' + treesampleNamesLong[i].split('_GC')[1] + '_TP0'
            elif ('_ASM_' in treesampleNamesLong[i]): ## need to check if subsetting in for i, name in enumerate(sampleNames[ingroup_bool]) would be successful as well! 
                treesampleNamesLong[i] = 'OG_ASM' + treesampleNamesLong[i].split('_ASM')[1] + '_TP0'
            else:
                treesampleNamesLong[i] = f"{patients[i]}_TP{str(timepoints[i]).zfill(2)}_{kit_date[i]}_{smpl_type_abbreviations[locations[i]]}_{locations_long_names[locations[i]]}"
                if name in seq_types.keys():
                    treesampleNamesLong[i] = f"{treesampleNamesLong[i]}_ST{seq_types[name]}" ## add ST if available
                if original_plates[i] != '':
                    treesampleNamesLong[i] = f"{treesampleNamesLong[i]}_p{original_plates[i]}" ## add original plate of sample if available
                treesampleNamesLong[i] = f"{treesampleNamesLong[i]}__{name}" ## add name

        treesampleNamesLong = treesampleNamesLong[ingroup_wOutSpl_bool] # remove outgroup samples

        # sampleNamesDnapars : max 10c. Use numeric with 10c (works for up to 10^10 samples! )
        sampleNamesDnapars = [ "{:010d}".format(i) for i in range(len(treesampleNamesLong))]

        # translate index to nucleotide
        calls_for_tree = apy.idx2nts(calls_for_treei) # ATCGN translation
        # add reference nucleotide for all positions
        refgenome_nts = apy.extract_outgroup_mutation_positions(ref_genome_folder, apy.p2chrpos(p[goodpos2useTree],chrStarts));
        calls_for_tree = np.concatenate((refgenome_nts[:, None],calls_for_tree),axis=1) # first column now refgenome_nts; refgenome_nts[:, None] to make ndims (2) same for both

        reference_label = refgenome[:5] + '_ref'
        sampleNamesDnapars_wRef = np.append(reference_label,sampleNamesDnapars) # add name for outgroup // if this is changed it will cause issues with the per SNP tree as this has the Reference genome hardcoded!
        treesampleNamesLong_wRef = np.append(reference_label,treesampleNamesLong) # add name for outgroup // if this is changed it will cause issues with the per SNP tree as this has the Reference genome hardcoded!

        # ## build parsimony tree; added flag "buildTree": if false > only build dnaparse/fasta input files but do not calc tree 
        # Uncomment to get PS tree
        if np.shape(calls_for_tree)[0] > 0:
            outgroup_for_root = treesampleNamesLong_wRef[1:][outgroup_spl_idx][min_outgroupsplnti] ## note, reference is excluded to get right outgroup
            ## remove all other outgroups for calculations:
            redundant_outgroups_bool = np.invert(np.append(np.array(False), outgroup_bool)) + (treesampleNamesLong_wRef == outgroup_for_root)
            calls_for_tree_1outgrp = calls_for_tree[:, redundant_outgroups_bool]
            treesampleNamesLong_wRef_1outgrp = treesampleNamesLong_wRef[redundant_outgroups_bool]
            sampleNamesDnapars_wRef_1outgrp = sampleNamesDnapars_wRef[redundant_outgroups_bool]
            apy.generate_tree(calls_for_tree_1outgrp,treesampleNamesLong_wRef_1outgrp,sampleNamesDnapars_wRef_1outgrp,filetag=refgenome,buildTree='PS',root_on_smpl = outgroup_for_root, coloring_pattern = 'TP')
            
            # Uncomment to get ML tree
            ## get number of positions within core genome (NOTE: This takes some time as it handles the coverage matrix!)
            good_smpls_of_subj = np.array(subj_specific_samples)[goodsamples]
            ingrp_1outgrp_smpl = good_smpls_of_subj[redundant_outgroups_bool[1:]] ## remove from bool the reference!
            num_pos_coregenome, num_pos_invar = apy.get_num_core_genome_pos(cmtdir, ingrp_1outgrp_smpl, p, goodpos2useTree, genomeLength, filter_parameter_site_per_sample, filter_parameter_site_across_samples, median_coverage_position_with_cov = False)

            ml_model = 'GTR+G+ASC_FELS{' + str(num_pos_invar) + '}'
            ## NOTE: --force model_lh_impr is required to allow raxml to run on even if the branch length optimization algorithm allows a worse likelihood (required for clonal trees with low diversity to let raxml run on!), However, when using validate the generated tree output!
            apy.generate_tree(calls_for_tree_1outgrp,treesampleNamesLong_wRef_1outgrp,sampleNamesDnapars_wRef,buildTree='ML',root_on_smpl = outgroup_for_root, raxml_model=ml_model, additional_raxml_parameters = '--threads 4 --seed 123 --force model_lh_impr --blmin 1e-10 --precision 16') ## set precision above blmin to collspase tree upon clonal nodes
            
        # %% Write sampleNames of goodSamples w/o outgroup to file
        # Necessary to filter coverage matrix and remove false positives
        np.savetxt('sampleNames_goodsamples_noOut.txt',sampleNames[ingroup_bool],fmt='%s')
        np.save('sampleNames_goodsamples_noOut', sampleNames[ingroup_bool]) # appends .npy
        np.savetxt('sampleNames_goodsamples_treeNamesLongAll.txt',treesampleNamesLong,fmt='%s')

        # %%
        # =============================================================================
        #     Plot barchart with fwd/rev read count for each position for each allele
        # =============================================================================
        # get dataframes that carry fwd/rev coverage for all bases (incl. N)      
        [lod_fwd_cov,lod_rev_cov] = apy.build_dataframe_coverage_info(goodpos2useTree,NTs,sampleNames,maNT,minorNT,coverage_forward_strand,coverage_reverse_strand,maf,minorAF)
        
        # get chr and pos for reelvant SNPs
        chr_pos_gp = contig_positions[goodpos2useTree,]
        
        # loop over SNPs and plot. 
        # All results in: pdf/coverage_snp_fwd_rev/chr_poos_locID_anno_qual.pdf
        matplotlib.use('agg')
        apy.plot_coverage_fwd_rev_stacked(chr_pos_gp,annotation_mutations_dedup,lod_fwd_cov,lod_rev_cov,timestamp, f'{subject_fld_label}_{refgenome}')

        # %%
        # =============================================================================
        #     Make a tree for each SNP location 
        # =============================================================================
        
        ## create empty tree_counting folder and add for_tree_labeling.csv
        if np.shape(calls_for_tree)[0] > 0:
            ## create empty tree_counting folder and add for_tree_labeling.csv
            apy.build_table_for_tree_labeling(apy.p2chrpos(p[goodpos], chrStarts),treesampleNamesLong_wRef_1outgrp,calls_for_tree_1outgrp)

            ## add tree for each mutation dsiplaying mutation (countMutations.py)
            os.chdir('tree_counting/snps')
            tree_path = f'../../{refgenome}_latest_rerooted_{outgroup_for_root}.nwk.tree'
            apy.generate_mutation_colored_tree(tree_path, "for_tree_labeling.csv", outgroup_name_for_lca = outgroup_for_root, color_lca = False)
            os.chdir('../../')
        
        # =============================================================================
        # %% 
        # Store SNP table
        # =============================================================================
        ## SOM SNP Table. Should be a function.

        apy.write_snp_table(calls[ np.ix_(goodpos2useTree, ingroup_bool) ], treesampleNamesLong[ingroup_bool], annotation_mutations_dedup)
        
