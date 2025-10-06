# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Analysis started 2023/04

Public data analysis from Ehormaechei isolates 2024/01
analysis module from 2023/04

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
# - tree branch coloring according to majority timepoint

# %% apy module (contains all functions)

import sys,os,re,subprocess
import glob,pickle,time,datetime,warnings
import numpy as np
from scipy import stats
import pandas as pd
import seaborn as sns
import gzip

from pylab import * #what is that // Maybe change to matplotlib.pyplot?

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
print(f'Starting the script from {os.getcwd()}')

## Directory to analysispy_module.py
SCRIPTS_DIRECTORY = os.getcwd() + '/modules/'
sys.path.insert(0, SCRIPTS_DIRECTORY)

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
                                        'max_fraction_ambigious_samples_pregubbins' : 0.01, #across samples per position, be as restrict as you want to be in the end, the additional post gubbins filter is just there to remove sites which have high ambiguity and interfer with the tree construction!
                                        'max_fraction_ambigious_samples' : 0.1, #across samples per position
                                        'min_median_coverage_position' : 3, #across samples per position
                                        'distance_for_nonsnp' : 500, #region in bp on either side of goodpos that is considered for recombination
                                        'corr_threshold_recombination' : 0.75 #minimum threshold for correlation to remove site 
                                        }
gubbins_params = {\
                  'conda_name': 'gubbins', ## to install gubbins: `conda create -n "gubbins" -c bioconda -c conda-forge -c r gubbins`
                  'min-snps': 3,
                  'min-window-size': 100,
                  'max-window-size': 10000,
                  'iterations': 10,
                  'first-tree-builder': 'rapidnj',
                  'first-model': 'JC',
                  'tree-builder': 'raxmlng',
                  'model': 'GTRGAMMA',
                  'kwargs': '--extensive-search',
                  'threads': 64,
                  'seed': 123
                }
gubbins_cond = 'gubbinsdefault_ambig10'

## how far upstream of the nearest gene to annotate something a promoter
## mutation (not used if no annotation)
promotersize=250;

reuse_latest_cutoff_vals = False
on_the_fly_eval = False ## reduces max_fraction_ambigious_samples for samples with extensive SNPs to lower fraction (> 1000 SNPs)
## if old values shall be taken, all old shall be used and this need to be hard coded to set off!
if reuse_latest_cutoff_vals == True: 
    on_the_fly_eval = False


# %%  cd(workingdir)
analysis_run = 'public_data_analysis'

main_workingdir = os.path.expanduser('~/data/mf_2020_hap/data_analysis/2023/' + analysis_run + '/analysis/')
os.chdir(os.path.expanduser(main_workingdir))

case_run='case'
cand_mut_table_dirs = f'/../{case_run}/2-candidate_mutation_table'

## get species and their reference genomes (stored in subdir name) + 
## remove hidden files from list
refgenomes_ls = [dir for dir in os.listdir(os.getcwd() + cand_mut_table_dirs) if not dir.startswith('.')]
refgenomes_ls = list(sort(refgenomes_ls)) ## sort refgenomes 

# %% Other variables
# TIMESTAMP
# used in output files
ts = time.time()
timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S') # WARN: fixed timestamp not yet build in

## Nucleotides
NTs = np.array(['A','T','C','G'],dtype=object) # NTs='ATCG'


# %% make dict of species as keys and patients in which those species have been found
## Read in samples csv fot that
samples_csv = glob.glob(os.path.expanduser('~/data/mf_2020_hap/data_analysis/2023/' + analysis_run + '/mapping/samples.csv'))
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

# %% Define variables to store subject-specific results 
para_evo_cand = np.array([],dtype=object)
mut_genes_df = {}
snp_freq_shifts_anno_lod = []
regression_res_lod = {} # stores all the results of the regression, and turned to pandas DF and saved
annotation_mutation_allParaSignal = {} # store all annotation_mutation for all detected candidates
mol_clock_data_dc = {} # dc to store for each subject the calculated data for mol clock estimation
allP_gene_paraevo = {} # store all para evo candidates that fullfill bonferroni

# %% Initiate loop over all subjects
# subject_fld_name = ["P0007"]
# refgenome = refgenomes_ls[0]

for refgenome, subject_fld_name in species_in_patients_dict.items():
    
    for subj_idx,subject_fld_label in enumerate(subject_fld_name):
        print('\n'+'Process ' + subject_fld_label + ' and ' + refgenome)

        dir_timestamp = timestamp.split('_')[0].replace('-', '_')
        ## make per species a separate folder
        if not os.path.exists(main_workingdir + case_run + '_' + subject_fld_label + '_' + refgenome):
            ## generate dir 
            os.mkdir(main_workingdir + case_run + '_' + subject_fld_label + '_' + refgenome)
        if not os.path.exists(main_workingdir + case_run + '_' + subject_fld_label + '_' + refgenome + '/' + dir_timestamp + "_" + refgenome + "_" + gubbins_cond):
            ## generate dir 
            os.mkdir(main_workingdir + case_run + '_' + subject_fld_label + '_' + refgenome + '/' + dir_timestamp + "_" + refgenome + "_" + gubbins_cond)
            os.mkdir(main_workingdir + case_run + '_' + subject_fld_label + '_'  + refgenome + '/' + dir_timestamp + "_" + refgenome + "_" + gubbins_cond + '/pdf')
        os.chdir(main_workingdir + case_run + '_' + subject_fld_label + '_' + refgenome + '/' + dir_timestamp + "_" + refgenome + "_" + gubbins_cond)

        ref_genome_folder = '/home/fenk/data/mf_2020_hap/data_analysis/metadata/reference_genomes/' + refgenome

        # %% 
        # =============================================================================
        # Load data from candidate_mutation
        # =============================================================================
        # Import candidate_mutation_table
        cmtdir = main_workingdir + cand_mut_table_dirs + '/' + refgenome + '/'
        cmtFile = 'candidate_mutation_table_refgenomeoutgroups_dedup.pickle.gz'
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
                            'promotersize': promotersize,
                            'gubbins_params': gubbins_params}
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
                            'promotersize': promotersize,
                            'gubbins_params': gubbins_params}
            apy.save_cutoff_vals(timestamp, set_filters)
            print('Setting new cut-off values!')
        
        # %% 
        # get median cov per samples based on all positions in counts (aka p)
        median_cov_p_per_sample = coverage_stats[:, 10] ## each row is sample, col [0-10) == covg bins 1x, 2x... >10x; [10]==median covg; [11]==mean; [12]==stddev 

        # % indel counter reduction to indel count
        # indel_counter: The first statistic is the number of reads (at this position and in this sample) that 
        # support an indel. The second statistics is the number of reads (at this position and in this sample) that support a deletion.
        # for now we need only the count for indels overall 
        
        ## NEED TO BE FIXED AGAIN!
        indel_counter = indel_counter[:,0,:].transpose()
        # indel counter >> 50% indel row:1; row2 >> deletion
        
        # %% Save everything using more permanent sample names to avoid overwriting
        
        sampleNames_all=np.asarray(sampleNames,dtype=object)
        quals_all=-quals;
        counts_all=counts;
        coverage_all = counts_all.sum(axis=1).transpose() # axis=1 == rows; transpose needed > rows: pos and col: samples

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
        goodsamples = (median_cov_p_per_sample >= filter_parameter_sample_across_sites['min_median_coverage_to_include_sample']) | (in_outgroup == '1')[0]
        
        sampleNames = sampleNames_all[goodsamples]
        counts = counts_all[goodsamples , : , : ] # keep only level (samples) that fullfil filter!
        quals = quals_all[ : , goodsamples ]
        coverage = coverage_all[ : ,goodsamples]
        coverage_stats = coverage_stats[goodsamples, :]
        indels = indel_counter[:,goodsamples]
        in_outgroup = in_outgroup[:, goodsamples]
        
        ## open up memory
        del globals()['sampleNames_all']
        del globals()['counts_all']
        del globals()['quals_all']
        del globals()['coverage_all']
        del globals()['indel_counter']

        num_samples = len(sampleNames)
        
        coverage_forward_strand = counts[:,0:4,:].sum(axis=1).transpose()
        coverage_reverse_strand = counts[:,4:8,:].sum(axis=1).transpose()

        # %% Breakpoint: Too few samples passed filter
        if np.sum(goodsamples) < 2:
            print("Too few samples fullfill filter criteria! >> skip: " + refgenome)
            continue

        # %%
        # =============================================================================
        # Read in number of positions covered across genome
        # =============================================================================
        # obtain isolate specific-count of covered based required to correct observed mutation count (ie. for molecular clock inference and paraevo-poisson)
        # NOTE: This data has to be generated by user as it dependes on the assembly build. Here I put a placeholder
        numBaseGenome_covThreshold = coverage_stats[:, 8] ## bases covered > 8 per sample 
        numBaseGenome_min_cov_snv_filter = coverage_stats[:, int(filter_parameter_site_per_sample['min_cov_per_strand_for_call']*2)-1] ## bases covered ≥ min cov per strand
        ## open up memory
        del globals()['coverage_stats']

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
        
        ## Estimate outgroup (ancestral allele) from ALL samples added as outgroup to SM pipeline (ancnti* == major allele)
        pattern = re.compile('_GCA_|_GCF_') # string tag in sample Name to identify outgroup in staphAD data
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
        siteFilt = np.any(( (calls[:,ingroup_bool]>3).sum(axis=1) >= ((num_samples-np.sum(outgroup_bool)) * filter_parameter_site_across_samples['max_fraction_ambigious_samples_pregubbins']) \
                            ,np.median( coverage[:,ingroup_bool], axis=1) < filter_parameter_site_across_samples['min_median_coverage_position'] ),axis=0)
        calls[ siteFilt ,:] = 4 # sites that fail qc -> 4, for all samples incl. outgroup     
        
        # NOTE: func below takes forever with many SNPs...saved below
        [mutQual, mutQualIsolates] = apy.ana_mutation_quality(calls[:,ingroup_bool],quals[:,ingroup_bool]) # get FQ value for SNP across samples. mutQualIsolates contains sample indices for sample pair FQ based on. 
        mutQual = np.nan_to_num(mutQual, nan=-1) # turn mutQual nan's to -1; necessary to avoid later warning
        
        ## translate filtered calls of ingroup into goodpos. mutations we believe. fixedmutation part removed in v6.
        hasmutation = (calls != refnti_m) & (calls < 4) & (np.tile(mutQual,(1,num_samples)) >= 1) # consider only ingroup samples; mutQual >= 1 is very loose. Important filter with low qual data! refnt not ancnt (as we want to find snps _between_ samples not between samples and outgroup!)
        hasmutation[:,outgroup_bool] = False # put outgroup samples 4 in order to identify ingroup mutations only
        
        candpos = np.where( np.sum(hasmutation, axis=1) > 0 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!

        print(f'Identified {len(candpos)} candidate positions for SNPs')
        if len(candpos) == 0:
            with open(f'{timestamp}_no_goodpos_identified.txt', 'w') as fid:
                print('No variable positions identified. End of analysis.')
            continue

        cmtdir = main_workingdir + cand_mut_table_dirs + '/' + refgenome + '/'        
        # cmtFile_sub_gp = 'candidate_mutation_table_processed_gp.gubbins10reps_pad.pickle'
        cmtFile_sub_gp = f'candidate_mutation_table_processed_candpos.acrossallsmpls_{gubbins_cond}.pickle'
        with open(cmtdir + cmtFile_sub_gp, 'wb') as pickle_file:
            pickle.dump({
                'counts': counts,
                'quals': quals,
                'coverage_forward_strand': coverage_forward_strand,
                'coverage_reverse_strand': coverage_reverse_strand,
                'refnti_m': refnti_m,
                'p': p,
                'refgenome': refgenome,  
                'sampleNames': sampleNames,  
                'outgroup_bool': outgroup_bool,  
                'numBaseGenome_covThreshold': numBaseGenome_covThreshold,  
                'contig_positions': contig_positions,
                'mutantAF': mutantAF,
                'maf': maf,
                'maNT': maNT,
                'minorNT': minorNT,
                'minorAF': minorAF,
                'calls': calls,
                'hasmutation': hasmutation,
                'candpos': candpos,
                'ingroup_idx': ingroup_idx, 
                'ref_genome_folder': ref_genome_folder,
                'chrStarts': chrStarts,
                'scafNames': scafNames, 
                }, pickle_file)
        subprocess.run(f'gzip --force {cmtdir}{cmtFile_sub_gp}', shell = True) ## Note: CLI gzip is faster than pythons gzip.open for larger files!
        
        ########################################################
        ########################################################
        ########################################################

        lineage_dict = {refgenome: [sample for sample in sampleNames[ingroup_bool]]}
        if gubbins_params != {}:
            calls[:, ingroup_idx], hasmutation[:, ingroup_idx] = apy.identify_recombination_sites_via_gubbins(calls[:, ingroup_idx], p, sampleNames[ingroup_idx], hasmutation[:, ingroup_idx], ref_genome_folder, chrStarts, scafNames, genomeLength, gubbins_params, lineage_dict, filter_parameter_site_across_samples['max_fraction_ambigious_samples'], mask_recomb_per_sample=True)

        goodpos = np.where( np.sum(hasmutation, axis=1) > 0 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!
        
        print(f'{np.shape(goodpos)[0]} good positions have been identified.\n')
        
        counts = counts[:, :, goodpos]
        quals = quals[goodpos, :]
        coverage_forward_strand = coverage_forward_strand[goodpos, :]
        coverage_reverse_strand = coverage_reverse_strand[goodpos, :]
        refnti_m = refnti_m[goodpos, :]
        p = p[goodpos]
        contig_positions = contig_positions[goodpos, :]
        mutantAF = mutantAF[goodpos, :]
        maf = maf[goodpos, :]
        maNT = maNT[goodpos, :]
        minorNT = minorNT[goodpos, :]
        minorAF = minorAF[goodpos, :]
        calls = calls[goodpos, :]
        hasmutation = hasmutation[goodpos, :]
        mutQual = mutQual[goodpos]

        try:
            ## open up memory 
            del globals()['coverage']
            del globals()['refnt']
            del globals()['refnti']
        except:
            print('in troubleshoot mode')

        ## need to be set to plot in new window!!!
        if np.shape(goodpos)[0] > 0: 
            ##%matplotlib qt
            apy.plot_interactive_scatter_barplots(p,mutQual,'pos','qual', sampleNames,counts,timestamp, filter_parameter_site_across_samples['corr_threshold_recombination'], subject_fld_label, refgenome, saveplots = True)

            ## need to be set to ensure saving of plot in right size without blocking a figure window
            matplotlib.use('agg')
            apy.plot_minorAF_rates(minorAF, hasmutation, sampleNames, subject_fld_label, refgenome, timestamp, '_goodpos')
            apy.plot_quals_rates(quals, hasmutation, sampleNames, subject_fld_label, refgenome, timestamp, '_goodpos')
            apy.plot_count_rates(counts, hasmutation, sampleNames, subject_fld_label, refgenome, timestamp, '_goodpos')
            apy.plot_Ns_rate(maNT, sampleNames, p, genomeLength, 10, subject_fld_label, refgenome, timestamp, '_goodpos')
        
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
            unpol_snps = np.sum(outsplsnti[:,i]==4)
            outgrp_splname = sampleNames[outgroup_spl_idx][i]
            polarized_site_file_txt = '\n'.join([polarized_site_file_txt, f'{outgrp_splname}:\t\t{unpol_snps} unpolarized SNPs'])
        ## select the outgroup which is phylogenetically the closest
        min_outgroupsplnti = np.where(sampleNames[outgroup_spl_idx] == 'Equasihormaechei_GCF_004331385_1')[0][0]
        outsplnti = outsplsnti[:,min_outgroupsplnti]
        outsplnti_m = np.tile(outsplnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele    
        polarized_site_file_txt = '\n'.join([polarized_site_file_txt, f'\nUnresolved goodpos in outsplnti outgroup:\t\t{np.sum(outsplnti==4)}'])
        polarized_site_file_txt = '\n'.join([polarized_site_file_txt, f'Closest outgroup sample name:\t\t{sampleNames[outgroup_spl_idx][min_outgroupsplnti]}'])
                
        ## out chimera. based on outsplnti but all NA (==4) are replaced by major allele call in ancnti
        # we use this for molecular clock and derived allele freq (but not tree!)
        outchimerasplancnti = np.array([ancnti[i] if nti==4 else nti for i,nti in enumerate(outsplnti) ])
        polarized_site_file_txt = '\n'.join([polarized_site_file_txt, f'Unresolved goodpos in chimera outgroup:\t\t{np.sum(outchimerasplancnti==4)}'])
        polarized_site_file_txt = '\n'.join([polarized_site_file_txt, f'Resolved goodpos in chimera outgroup:\t\t{len(goodpos) - np.sum(outchimerasplancnti==4)} ({((len(goodpos) - np.sum(outchimerasplancnti==4)) / len(goodpos)) * 100}%)'])
        print(polarized_site_file_txt)

        with open(f'{timestamp}_polarized_sites_by_outgroup.txt', 'w') as fo_pol:
            fo_pol.write(polarized_site_file_txt)


        ## plot which outgroup sample is polarizing which SNP and which SNPs are not polarized at all 
        outsplsnti_chim = np.concatenate((outsplsnti, np.reshape(outchimerasplancnti, (-1, 1))), axis = 1) ## include the major (chimeric) calls
        outgrpspl_idx_sort = np.argsort(np.sum(outsplsnti_chim==4, axis = 0)) ## idx to sort samples based on the number of SNPs they polarize 
        outgrp_sampleNames = np.append(sampleNames[outgroup_spl_idx], 'Chimeric outgroup') ## extract sample names and add chimeric sample name
        outgrp_sampleNames = outgrp_sampleNames + '\n(unpol. SNPs: ' + (np.sum(outsplsnti_chim==4, axis = 0)).astype(str) + ')' ## add information how many snps not resolved by outgroup
        outgrp_sampleNames_s = outgrp_sampleNames[outgrpspl_idx_sort] ## sort sample names 
        outgrp_sampleNames_sext = np.tile(outgrp_sampleNames_s, len(goodpos)) ## repeat samplenames n times 

        gnome_pos_snps = p ## get genomic positions of goodpos
        gnome_pos_snps_ext = np.repeat(gnome_pos_snps, len(outgrpspl_idx_sort)) ## repeat each pos * outgroup sample number 
        
        outgrpspl_pol_pos = outsplsnti_chim[:, outgrpspl_idx_sort] ## get calls of outgroups and sort by idx 
        outgrpspl_pol_pos_1d = outgrpspl_pol_pos.ravel() ## convert to 1d array for plotting
        outgrpspl_pol_pos_1d[outgrpspl_pol_pos_1d != 4] = True
        outgrpspl_pol_pos_1d[outgrpspl_pol_pos_1d == 4] = False
        
        # %%
        # =============================================================================
        #  Make table with annotations
        # =============================================================================

        order = np.arange(sampleNames.shape[0])
        num_contigs = np.max(contig_positions[:,0]);
        contig_lengths = np.append(chrStarts[1:], genomeLength) - chrStarts
        
        # Uncomment the following if there is an annotated genome
        # NOTE: annotation is read from *.gff file in ref_folder!
        ## If no ortholog_df is there just leave it out
        ortholog_column_tag = 'HAP2020_mf_' + refgenome[:-7]
        annotation_genes = apy.parse_gff(ref_genome_folder,scafNames, forceReDo=False) # ref_genome_folder+"/annotation_genes.pandas.py.pk1"
        annotation_mutations = apy.annotate_mutations_v2(annotation_genes , p , refnti_m , ancnti_m, calls , counts , hasmutation, mutQual.flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP
        mutational_priorities = ['N', 'S', 'R', 'P', 'I', 'U']
        annotation_mutations_dedup = apy.deduplicate_annotation_mutations(annotation_mutations, timestamp, mutational_priorities)
        pd.concat(annotation_genes).to_csv(timestamp + 'annotation_genes.csv')
        
        annotation_mutations.to_csv(timestamp + 'annotation_mutations.csv')

        ## save everything
        cmtdir = main_workingdir + cand_mut_table_dirs + '/' + refgenome + '/'        
        # cmtFile_sub_gp = 'candidate_mutation_table_processed_gp.gubbins10reps_pad.pickle'
        cmtFile_sub_gp = f'candidate_mutation_table_processed_gp.acrossallsmpls_{gubbins_cond}.pickle'
        with open(cmtdir + cmtFile_sub_gp, 'wb') as pickle_file:
            pickle.dump({
                'counts': counts,
                'quals': quals,
                'coverage_forward_strand': coverage_forward_strand,
                'coverage_reverse_strand': coverage_reverse_strand,
                'refnti_m': refnti_m,
                'p': p,
                'refgenome': refgenome,  
                'sampleNames': sampleNames,  
                'outgroup_bool': outgroup_bool,  
                'numBaseGenome_covThreshold': numBaseGenome_covThreshold,  
                'contig_positions': contig_positions,
                'mutantAF': mutantAF,
                'maf': maf,
                'maNT': maNT,
                'minorNT': minorNT,
                'minorAF': minorAF,
                'calls': calls,
                'hasmutation': hasmutation,
                'goodpos': goodpos,
                'annotation_mutations': annotation_mutations,
                'annotation_genes': annotation_genes
                }, pickle_file)
        subprocess.run(f'gzip --force {cmtdir}{cmtFile_sub_gp}', shell = True) ## Note: CLI gzip is faster than pythons gzip.open for larger files!

        # %% 
        # =============================================================================
        #     parsimony tree. (incl. samples filter)
        # =============================================================================    
        # define pos to use in tree
        # goodpos2useTree = goodpos[0::1000] #(1:1000); %trim for easier tree view; TDL called quality_positions
        goodpos2useTree = np.array([i for i in range(len(goodpos))])

        # get data and filter for goodpos
        calls_for_treei = calls; 
        calls_for_treei = calls_for_treei[ goodpos2useTree, : ]
        calls_for_treei = calls_for_treei[ : , ingroup_wOutSpl_bool ]
        
        # build sampleNames w/ metainfo
        treesampleNamesLong = np.copy(sampleNames)

        # sampleNamesDnapars : max 10c. Use numeric with 10c (works for up to 10^10 samples! )
        sampleNamesDnapars = [ "{:010d}".format(i) for i in range(len(treesampleNamesLong))]

        # translate index to nucleotide
        calls_for_tree = apy.idx2nts(calls_for_treei) # ATCGN translation
        # add reference nucleotide for all positions
        refgenome_nts = apy.extract_outgroup_mutation_positions(ref_genome_folder, apy.p2chrpos(p[goodpos2useTree],chrStarts));
        calls_for_tree = np.concatenate((refgenome_nts[:, None],calls_for_tree),axis=1) # first column now refgenome_nts; refgenome_nts[:, None] to make ndims (2) same for both

        reference_label = refgenome.split('_')[1][:5] + '_ref'
        sampleNamesDnapars_wRef = np.append(reference_label,sampleNamesDnapars) # add name for outgroup // if this is changed it will cause issues with the per SNP tree as this has the Reference genome hardcoded!
        treesampleNamesLong_wRef = np.append(reference_label,treesampleNamesLong) # add name for outgroup // if this is changed it will cause issues with the per SNP tree as this has the Reference genome hardcoded!

        # ## build parsimony tree; added flag "buildTree": if false > only build dnaparse/fasta input files but do not calc tree 
        # Uncomment to get PS tree
        ## note: rerooting is just implemented for PS tree yet!
        if np.shape(calls_for_tree)[0] > 0:
            outgroup_for_root = treesampleNamesLong_wRef[1:][outgroup_spl_idx][min_outgroupsplnti] ## note, reference is excluded to get right outgroup
            ## remove all other outgroups for calculations:
            redundant_outgroups_bool = np.invert(np.append(np.array(False), outgroup_bool)) + (treesampleNamesLong_wRef == outgroup_for_root)
            calls_for_tree_1outgrp = calls_for_tree[:, redundant_outgroups_bool]
            treesampleNamesLong_wRef_1outgrp = treesampleNamesLong_wRef[redundant_outgroups_bool]
            sampleNamesDnapars_wRef = sampleNamesDnapars_wRef[redundant_outgroups_bool]
            
            ## write alignment file for ML tree:
            apy.write_calls_sampleName_to_fasta(calls_for_tree_1outgrp,treesampleNamesLong_wRef_1outgrp,'snp_msa')
                
            # Uncomment to get ML tree
            raxml_model = 'GTR+G'
            ml_tree_timestamp = apy.generate_tree(calls_for_tree_1outgrp,treesampleNamesLong_wRef_1outgrp,sampleNamesDnapars_wRef,filetag=refgenome,buildTree='ML',raxml_model=raxml_model,root_on_smpl = outgroup_for_root, additional_raxml_parameters='--threads 64 --seed 123')
            
            ## infer sequence on each node (needed for molecular clock calculation!)
            ml_suffix = raxml_model.replace('+', '').replace('{', '').replace('}', '')
            apy.create_ancestral_reconstruction(f'{ml_tree_timestamp}_ML_{ml_suffix}.raxml.bestTreeCollapsed',treesampleNamesLong_wRef_1outgrp,calls_for_tree_1outgrp, outdir="node_sequence_inference")

        
        # %% Write sampleNames of goodSamples w/o outgroup to file
        # Necessary to filter coverage matrix and remove false positives
        np.savetxt('sampleNames_goodsamples_noOut.txt',sampleNames[ingroup_bool],fmt='%s')
        np.save('sampleNames_goodsamples_noOut', sampleNames[ingroup_bool]) # appends .npy
        np.savetxt('sampleNames_goodsamples_treeNamesLongAll.txt',treesampleNamesLong_wRef,fmt='%s')

        # =============================================================================
        # %% 
        # Repolarize SNPs to the ancestral sequence
        # =============================================================================
        ## get the ancestral sequence
        anc_seq = apy.parse_ancestral_fasta_treetime(f'{os.getcwd()}/node_sequence_inference/ancestral_sequences.fasta', node_label_to_use='NODE_0000000', convert_to_ntsidx = False)
        annotation_mutations_dedup['nt_anc'] = anc_seq
        
        annotation_mutations_dedup = apy.repolarize_snvs(annotation_mutations_dedup)
        annotation_mutations_dedup.to_csv(timestamp + 'annotation_mutations_dedup_repolarized.csv')

        # =============================================================================
        # %% 
        # Store SNP table
        # =============================================================================
        ## SOM SNP Table. Should be a function.
        ## get NT,sampleName df
        # build sampleNames w/ metainfo
        if 'ref' in treesampleNamesLong_wRef[0]:
            sampleNamesLong = treesampleNamesLong_wRef[1:] ## remove reference 
        else:
            sampleNamesLong = treesampleNamesLong_wRef
        sampleNamesLong = sampleNamesLong[ingroup_wOutSpl_bool] # remove outgroup samples
        sampleNamesLong = sampleNamesLong[ingroup_wOutSpl_onlyIngroup_bool] # remove outgroup
        
        # translate index to nucleotide
        calls_for_treei = calls; 
        calls_for_treei = calls_for_treei[ goodpos2useTree, : ]
        calls_for_treei = calls_for_treei[ : , ingroup_wOutSpl_bool ]
        calls_for_treei_ingroup = calls_for_treei[ : , ingroup_wOutSpl_onlyIngroup_bool ]
        calls_for_tree = apy.idx2nts(calls_for_treei_ingroup) # ATCGN translation
        snp_data = pd.DataFrame(calls_for_tree,columns=sampleNamesLong)

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
        
        snp_metadata = snp_metadata[['chr','pos','type','muts','locustag','gene','loc1','loc2','strand','product','nt_pos','aa_pos','nt_ref', 'nt_anc', 'nt_alt']]
        snp_metadata[['nt_anc']] = snp_metadata[['nt_anc']].replace('.','?') # turn NA to '?', similar to NT data 
        snp_table = pd.concat([snp_metadata,snp_data.reset_index(drop=True)], axis=1)
        
        ## store
        with open('snp_table.csv', 'w') as file:
            snp_table.to_csv(file, header=True, index=False)

        # %%
        # =============================================================================
        # Estimate substitution rate
        # =============================================================================
        ## DONE using all subjects mapped to USA300

        allmuts = np.array(['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG'],dtype=object)
        # AT, TA  0
        # AC, TG  1
        # AG, TC  2
        # GC, CG  3
        # GT, CA  4
        # GA, CT  5
        allmuts_types = np.array([0,2,1,0,1,2,5,4,3,4,5,3])
        mutationalspectrum = [1/12] * 12 # uniform distribution. 

        mrca_seq = apy.get_internal_node_nt_calls('node_sequence_inference/annotated_tree.nexus','node_sequence_inference/ancestral_sequences.fasta',treesampleNamesLong[ingroup_bool])
        mrca_seq_ntidx = apy.nts2idx(np.array([x for x in mrca_seq]))

        ## loop over all ref/anc nucleotides and count mutation types. count twice for each direction
        allmuts = np.array(['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG'],dtype=object)
        allmuts_counts = np.zeros( allmuts.shape,dtype=int )
        for i,anc in enumerate(mrca_seq): 
            obs_allele = np.unique(calls_for_tree[i,1:]) # Note: 1st col is reference!
            obs_allele = obs_allele[ ~(obs_allele == '?') & ~(obs_allele == anc) ]
            for j in obs_allele:
                idx = np.where( (anc+j)==allmuts)
                allmuts_counts[idx] += 1
                idx = np.where( (j+anc)==allmuts)
                allmuts_counts[idx] += 1
        mutationalspectrum = allmuts_counts/np.sum(allmuts_counts) # true mutationalspectrum // np.sum(allmuts_counts) == number of all observed and resolved SNVs (no ambiguities!) *2; mutational spectrum = frequency of mutational shift
        # store
        afile = open("mutationalspectrum.py.pk1", 'wb')
        pickle.dump(mutationalspectrum, afile)
        afile.close()
        

        # %%
        # =============================================================================
        # Parallel evolution module
        # =============================================================================
        # define dict with parameters for inference
        
        parameters = {
                'NumTrialsSim':1000,
                'Min_num_mutations_cand':2, # minimum number of mutations per gene candidate
                'Min_mutation_density_cand':1/1000, # minimum number of mutations per 1000 bp per gene candidate
                'ref_genome_folder':ref_genome_folder,
                'subjectID': refgenome, # used for naming pdf in pdf/adaptive_evo
                'substitution_spectrum':'mutationalspectrum.py.pk1', # put None if not calculated. Alternatively path
                'max_muts_per_gene_to_track': 250,
                'timestamp':timestamp
                }
        
        print('Hard fix for num_mutations --> might not always be correct (see homoplasies!)')
        print('need to test "homoplasy_simulation" function!')
        annotation_mutations_dedup['num_mutational_events'] = annotation_mutations_dedup['nts'].str.len()-1

        parameters['analsysis_params_output_name_folder'] = 'nomasked_sites_para_evo'
        [res_cand_nummut,annotation_mutation_paraSignal] = apy.parallel_evo_module( goodpos, contig_positions , annotation_mutations_dedup , annotation_genes , parameters, plot = True, ortholog = False )

        # %%
        # =============================================================================
        # Report para evo results w/ Bonferroni Poisson P val 
        # =============================================================================
        # get poisson p value for observed candidates of para evo 

        ## variables for lambda (expectation) calculation
        tot_num_mutations = len(goodpos)
        median_genome_length = np.median(numBaseGenome_covThreshold)

        ## bonferroni corrected p for subject (based on total number of annotated genes)
        num_annotated_genes = sum([t.shape[0] for t in annotation_genes if not t.empty])
        alpha_sig_level = 0.05 # pval considered significant for multiple testing correction
        p_bonf_corrected = alpha_sig_level/num_annotated_genes

        
        lod_gene_paraevo = []    
        for i in range(res_cand_nummut.shape[0]):
            para_gene_dict = {}
            ### calculate poisson lambda - expectation
            ## get gene of para evo candidate gene
            contig_paraevo_candmut = res_cand_nummut[i,0].split('_', 1)[0]
            locustag_paraevo_candmut = res_cand_nummut[i,0].split('_', 1)[1]
            idxGeneAll = np.where((annotation_mutations_dedup['locustag']==locustag_paraevo_candmut) & 
                                (annotation_mutations_dedup['chr'].astype(str)==contig_paraevo_candmut))[0] 
            idxGene = idxGeneAll[0] # only take idx of first instance of gene
            lengthGene = len(annotation_mutations_dedup['sequence'][idxGene])
            my_lambda = tot_num_mutations * (lengthGene/median_genome_length)
            ### poisson cdf for num-obs_mut or more
            obs_mut = res_cand_nummut[i,1]
            p_poisson_obs = 1-stats.poisson.cdf(obs_mut-1,my_lambda)
            ### record if below bonferroni corrected p for that patient: p/#genes
            # record all candidates
            para_gene_dict['Patient'] = subject_fld_label
            para_gene_dict['Contig'] = annotation_mutations_dedup['chr'][idxGene]
            para_gene_dict['Gene_Start'] = annotation_mutations_dedup['loc1'][idxGene]
            para_gene_dict['Gene_End'] = annotation_mutations_dedup['loc2'][idxGene]
            para_gene_dict['Gene_length'] = annotation_mutations_dedup['loc2'][idxGene]- annotation_mutations_dedup['loc1'][idxGene]
            para_gene_dict['Strand'] = annotation_mutations_dedup['strand'][idxGene]
            para_gene_dict['Gene_id'] = annotation_mutations_dedup['locustag'][idxGene]
            para_gene_dict['Gene_name'] = annotation_mutations_dedup['gene'][idxGene]
            para_gene_dict['Product_description'] = annotation_mutations_dedup['product'][idxGene]
            para_gene_dict['Poisson_p_value'] = p_poisson_obs
            para_gene_dict['Bonferroni_p_value'] = p_bonf_corrected
            para_gene_dict['Number_mutations'] = res_cand_nummut[i,1]
            para_gene_dict['ProbSim'] = res_cand_nummut[i,6]
            para_gene_dict['Mutations'] = [annotation_mutations_dedup['muts'][m][0] for m in idxGeneAll]
            para_gene_dict['Mutation_type'] = [annotation_mutations_dedup['type'][m] for m in idxGeneAll]
            para_gene_dict['Mutation_position'] = [annotation_mutations_dedup['pos'][m] for m in idxGeneAll]
            para_gene_dict['Translation'] = annotation_mutations_dedup['translation'][idxGene]
            para_gene_dict['Sequence'] = annotation_mutations_dedup['sequence'][idxGene]            
            lod_gene_paraevo.append(para_gene_dict)
        df_gene_paraevo = pd.DataFrame(lod_gene_paraevo)  # turn patient results to dataframe and save
        df_gene_paraevo.to_csv('para_evo_snp_table.csv', index = False)
        