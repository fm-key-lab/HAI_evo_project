# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Analysis started 2024/08

denovo analysis with SM runs from 2024/08
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

import sys,os,re
import glob,pickle,time,datetime
import numpy as np
from scipy import stats
import pandas as pd
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", message=".*timestamp seems very low.*")

from pylab import * #what is that // Maybe change to matplotlib.pyplot?

# enable plotting in separate windows
matplotlib.use('Qt5Agg')
plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

## Directory to analysispy_module.py
SCRIPTS_DIRECTORY = os.getcwd() + '/../../modules/'
#SCRIPTS_DIRECTORY = os.getcwd() + "/modules/"
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
                                        'max_fraction_ambigious_samples' : 1, #across samples per position // 1 == filter off (all sites can be (in theory) ambigous)
                                        'min_median_coverage_position' : 0, #across samples per position; NOTE: just set if you want to look at core genome, otherwise, set threshold min_median_coverage_position_with_cov
                                        'min_median_coverage_position_with_cov' : 4, #across samples per position with any coverage (this accounts for deletions between samples!)
                                        'distance_for_nonsnp' : 500, #region in bp on either side of goodpos that is considered for recombination
                                        'corr_threshold_recombination' : 0.75 #minimum threshold for correlation to remove site 
                                        }

# Remove indels within samples that are not high quality
filter_parameter_indel_per_sample = {\
                                    'min_iaf_for_call' : 0.8, # minimum indel allele frequency on individual samples, calls; NOTE: Indels are more variable and tend to have higher variations in indel allele frequencies due to misalignments despite high reliability!
                                    'min_cov_at_pos' : 3,  
                                    'min_indel_length' : 1,  # minimum indel length --> 0 == ref and alt are same length --> SNV!; on individual samples, calls
                                    'max_indel_support' : -10,  # the smaller the value the better the support for the best alternative indel; on individual samples, calls; GL is not PHRED-like (-10*log10(p)) but (log10(p))
                                    'max_fraction_ambigious_samples': 0.2 # set to lower values compared to SNVs to avoid highly variable sites (e.g. transposons)
                                    }

## how far upstream of the nearest gene to annotate something a promoter
## mutation (not used if no annotation)
promotersize=250;

reuse_latest_cutoff_vals = False
on_the_fly_eval = False ## reduces max_fraction_ambigious_samples for samples with extensive SNPs to lower fraction (> 1000 SNPs)
## if old values shall be taken, all old shall be used and this need to be hard coded to set off!
if reuse_latest_cutoff_vals == True: 
    on_the_fly_eval = False

# %%  cd(workingdir)
analysis_run = 'denovo_2024_08'

main_workingdir = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/' + analysis_run + '/analysis/backup/')
os.chdir(os.path.expanduser(main_workingdir))

cand_mut_table_dirs = '/../../case/2-candidate_mutation_table_2024_10_18'

## get species and their reference genomes (stored in subdir name) + 
## remove hidden files from list
refgenomes_ls = [dir for dir in os.listdir(os.getcwd() + cand_mut_table_dirs) if not dir.startswith('.')]
refgenomes_ls = sort(refgenomes_ls).tolist() ## sort refgenomes 

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
specimenLog_path = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/Drylab/metadata/*specimenlog.csv')
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
            SpecimenLog["sampleids"].append(line[3].replace('-', '_'))
            SpecimenLog["timepoint"].append(line[5])
            SpecimenLog["plate"].append(line[6])

# %% make dict of species as keys and patients in which those species have been found
## Read in samples csv fot that
samples_csv = glob.glob(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/' + analysis_run + '/mapping/*samples*.csv'))
samples_csv = sorted(samples_csv)[-1] ## return newest samples_case.csv

species_in_patients_dict = {}
with open(samples_csv, 'r') as fo:
    for line in fo: 
        if not 'ReferenceGenome' in line: ## skip first line 
            line = line.strip().split(',')
            ## include up to 3 leading 0 (5 chars for patient ID each)
            if not line[4] in ['Key', 'Outgroup', 'P17']:
                pat_id = line[4][0] + str(line[4][1:]).zfill(4)
                ## samples csv is ordered by Path,Sample,ReferenceGenome,ProviderName,Subject
                for refgenome in line[2].split(' '): ## account for multiple ref genomes stated in samples csv:
                    if not refgenome in species_in_patients_dict.keys(): ## if reference genome not in keys, then generate new key
                        species_in_patients_dict[refgenome] = [pat_id]
                    elif pat_id not in species_in_patients_dict[refgenome]: ## if patient id not in values of key, then store 
                        species_in_patients_dict[refgenome].append(pat_id)
## order dictionary
species_in_patients_dict = dict(sorted(species_in_patients_dict.items()))


## Load ortholog df 
ortholog_data = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/assembly/7-ortholog_identification/annotation_orthologs.tsv') 
ortholog_df = pd.read_csv(ortholog_data,sep="\t",index_col=0)

## list of samples identified to be contaminated 
contaminated_samples = ['P07_3453', 'P07_2730', 'L00071_07_1', 'L00071_02_1', 'P21_0043']


# %% Define variables to store subject-specific results 
para_evo_cand = np.array([],dtype=object)
snp_freq_shifts_anno_lod = []
regression_res_lod = {} # stores all the results of the regression, and turned to pandas DF and saved
annotation_mutation_allParaSignal = {} # store all annotation_mutation for all detected candidates
mol_clock_data_dc = {} # dc to store for each subject the calculated data for mol clock estimation
allP_gene_paraevo = {} # store all para evo candidates that fullfill bonferroni
mut_sig_species_dict = {} # store all unfiltered mutational signatures 

# %% Initiate loop over all subjects
subject_fld_name = ["P0021"]
refgenome = refgenomes_ls[12]

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
        elif 'denovo' in analysis_run:
            ref_genome_folder = os.path.expanduser(main_workingdir + '../../../../metadata/reference_genomes/') + refgenome

        # %% 
        # =============================================================================
        # Import candidate_mutation_table
        cmtdir = main_workingdir + cand_mut_table_dirs + '/' + refgenome + '/'
        cmtFile = 'candidate_mutation_table.pickle.gz'
        [quals,p,counts,in_outgroup,sampleNames,indel_counter,coverage_stats,indel_p,indel_call,indel_size,indel_depth,indel_support] = apy.read_candidate_mutation_table_pickle_gzip(cmtdir + cmtFile, indels = True)

        # %%
        # =============================================================================
        # Before doing any data parsing, save cutoff values to csv or read in latest cutoff vals if wanted
        if reuse_latest_cutoff_vals:
            try:
                ## load old data
                [filter_parameter_sample_across_sites, filter_parameter_site_per_sample, filter_parameter_site_across_samples, filter_parameter_indel_per_sample, promotersize] = apy.read_latest_cutoff_vals(main_workingdir + subject_fld_label + '_' + refgenome)
                print('Old cut-offs are used again!')
            except:
                ## if no old cut off values are stored, set the new ones
                set_filters = {'filter_parameter_sample_across_sites': filter_parameter_sample_across_sites, 
                              'filter_parameter_site_per_sample': filter_parameter_site_per_sample, 
                              'filter_parameter_site_across_samples': filter_parameter_site_across_samples, 
                              'filter_parameter_indel_per_sample': filter_parameter_indel_per_sample,
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
                              'filter_parameter_indel_per_sample': filter_parameter_indel_per_sample,
                              'promotersize': promotersize}
            apy.save_cutoff_vals(timestamp, set_filters)
            print('Setting new cut-off values!')
      
        ## rename samples in sampleNames if patient sample and it only has a 3 digits ID
        sampleNames = np.array([f"{smplname.split('_')[0]}_{smplname.split('_')[1].zfill(4)}_{'_'.join(smplname.split('_')[2:])}".strip("_") if ( (smplname[0] == 'P') & (smplname.split('_')[1].isdigit()) ) else smplname for smplname in sampleNames])
        ## rename 2 sample names which have lineage appendix 
        if refgenome == 'P21_Ecoli-ST23-c2_240825':
            sampleNames = np.array([smpl.replace('-c2', '') if smpl.endswith('-c2') else smpl for smpl in sampleNames])
        
        ## extract only samples (= entries) which are of the patient(s) of interest 
        subject_smpls_speclog = [sampleid for subject, sampleid in zip(SpecimenLog['Patient'], SpecimenLog['sampleids']) if subject == subject_fld_label] ## extract which sampleids are derived from current patient
        repeated_samples_bool = [(sampleid[-3:-1] == '_S') and sampleid[-1].isdigit() and sampleid.split('_S')[0] not in contaminated_samples for sampleid in sampleNames] ## extracts samples which are resequenced and just keeps those which have been identified to be contaminated previously
        subj_specific_samples = [idx for idx, sampleid in enumerate(sampleNames) if (((sampleid.split('_S')[0] in subject_smpls_speclog) & (sampleid not in contaminated_samples)) | ('_GCA_' in sampleid) | ('_GCF_' in sampleid) | ('_ASM_' in sampleid)) & (~repeated_samples_bool[idx])] ## extract indices which include samples of specific patient    
        
        ## subset matrices of cand_mutation_table:
        quals = quals[:, subj_specific_samples]
        counts = counts[subj_specific_samples, :, :]
        in_outgroup = in_outgroup[:, subj_specific_samples]
        sampleNames = sampleNames[subj_specific_samples]
        indel_counter = indel_counter[subj_specific_samples, :, :]
        coverage_stats = coverage_stats[subj_specific_samples, :]
        indel_call = indel_call[:, subj_specific_samples,:]
        indel_size = indel_size[:, subj_specific_samples]
        indel_depth = indel_depth[:, subj_specific_samples, :]
        indel_support = indel_support[:, subj_specific_samples, :]

        # %% 
        # get median cov per samples based on all positions in counts (aka p)
        # --> INDICES WILL CHANGE for new data!!!
        median_cov_p_per_sample = coverage_stats[:, 11] ## each row is sample, col [0-10) == covg bins 0x, 1x, 2x... >10x; [11]==median covg; [12]==mean; [13]==stddev 

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
        # Assign patient & timepoint number to each sample; Parse specimen name to human readable
        patients_all = np.empty(len(sampleNames_all), dtype=object) #np.zeros(shape=(1, len(sampleNames_all)), dtype=np.int) # row w/ 0's 
        timepoints_all = np.zeros(len(sampleNames_all), dtype=np.int64) # row w/ 0's
        kit_date_all = np.zeros(len(sampleNames_all), dtype='object') # row w/ 0's
        locations_all = np.zeros(len(sampleNames_all), dtype=np.int64) # row w/ 0's
        original_plate_all = np.zeros(len(sampleNames_all), dtype='object') # row w/ 0's

        locations_abbreviations = ['N','O','R','S','B','L','G','U']; 
        smpl_type_abbreviations = ['R', 'R', 'R', 'R', 'D', 'D', 'D', 'D']; ## routine = R, diagnostic = D
        locations_long_names = ['Nasal', 'Oral', 'Rectal', 'Skin', 'Blood', 'Lung', 'Gastric', 'Urine'];

        for i in range(len(sampleNames_all)):
            current_smpl = sampleNames_all[i]
            if current_smpl[-3:-1] == '_S':
                current_smpl = current_smpl[:-3]
            if ('_GCA_' in current_smpl) | ('_GCF_' in current_smpl):
                continue
            kit_idx = SpecimenLog['sampleids'].index( current_smpl )
            patients_all[i] =  SpecimenLog['Patient'][kit_idx]
            timepoints_all[i] =  SpecimenLog['timepoint'][kit_idx]
            kit_date_all[i] =  SpecimenLog['date'][kit_idx].replace('-', '')[0:8] ## concat YYYYMMDD and remove time (HH:MM)
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
        kit_date_ingroup_datetime = [datetime.datetime.strptime(date, '%Y%m%d').date() for date in kit_date_all[in_outgroup[0] == '0']] 
        first_occurence = min(kit_date_ingroup_datetime)
        timepoint_to_colonization_day_dict = {tp: (sorted(kit_date_ingroup_datetime)[idx] - first_occurence).days for idx, tp in enumerate(sorted(timepoints_all[in_outgroup[0] == '0']))} ## NOTE: extraction is done on sorted lists for TP and dates to have sorted dict output

        # %%
        # =============================================================================
        #     Read in genome information
        # =============================================================================
        [chrStarts, genomeLength, scafNames] = apy.genomestats(ref_genome_folder);
        contig_positions = apy.p2chrpos(p,chrStarts) # 1col: chr, 2col: pos on chr; for all p
        
        # %%
        # =============================================================================
        #     Annotate genome
        # =============================================================================
        # NOTE: annotation is read from *.gff file in ref_folder!
        ## If no ortholog_df is there just leave it out
        ortholog_column_tag = 'HAP2020_mf_' + refgenome[:-7] ## [:-7] removes the date of the reference genome
        annotation_genes = apy.parse_gff(ref_genome_folder,scafNames,ortholog_df[ortholog_column_tag], forceReDo=True) # ref_genome_folder+"/annotation_genes.pandas.py.pk1"

        # %% 
        # =============================================================================
        #     Remove undesired samples based on genome-wide coverage
        # =============================================================================
        # %% Define goodsamples and filter data
        lowcovsamples = sampleNames_all[ median_cov_p_per_sample < filter_parameter_sample_across_sites['min_median_coverage_to_include_sample'] ]
        print(lowcovsamples)
        np.savetxt('sampleNames_lowCovSamples.txt',lowcovsamples,fmt='%s')

        goodsamples =  (median_cov_p_per_sample >= filter_parameter_sample_across_sites['min_median_coverage_to_include_sample'])
        
        sampleNames = sampleNames_all[goodsamples]
        counts = counts_all[goodsamples , : , : ] # keep only level (samples) that fullfil filter!
        quals = quals_all[ : , goodsamples ]
        coverage = coverage_all[ : ,goodsamples]
        coverage_stats = coverage_stats[goodsamples, :]
        patients = patients_all[goodsamples]
        timepoints = timepoints_all[goodsamples]
        kit_date = kit_date_all[goodsamples]
        locations = locations_all[goodsamples]
        original_plates = original_plate_all[goodsamples]
        indels = indel_counter[:,goodsamples]
        indel_call = indel_call[:, goodsamples,:]
        indel_size = indel_size[:, goodsamples]
        indel_depth = indel_depth[:, goodsamples, :]
        indel_support = indel_support[:, goodsamples, :]
        
        num_samples = len(sampleNames)
        
        coverage_forward_strand = counts[:,0:4,:].sum(axis=1).transpose()
        coverage_reverse_strand = counts[:,4:8,:].sum(axis=1).transpose()

        # %% Breakpoint: Too few samples passed filter
        if np.sum(goodsamples) < 2:
            print("Too few samples fullfill filter criteria! >> skip: " + refgenome)
            #continue

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
        pattern = re.compile('_GCA_|_GCF_') # string tag in sample Name to identify outgroup in data
        outgroup_name = np.array(list(filter(pattern.search, list(sampleNames))))
        if outgroup_name.size == 0:
            print('No outgroup is amongst good samples. Proceeding without outgroup')
        outgroup_bool = np.isin(sampleNames , outgroup_name)
        
        # ingroup array (bool, idx) used later
        ingroup_bool = np.invert(outgroup_bool)
        ingroup_idx = np.nonzero(ingroup_bool)[0]

        ## bool for all single outgroup sample AND ingroup-samples with dedicated single outgroup sample
        outgroup_currSub = np.array(["GCA", 'GCF', 'ASM'])
        outgroup_spl_idx = [i for i,item in enumerate(sampleNames) if (outgroup_currSub[0] in item) | (outgroup_currSub[1] in item) | (outgroup_currSub[2] in item)]
        ingroup_wOutSpl_bool = np.copy(ingroup_bool)
        ingroup_wOutSpl_bool[outgroup_spl_idx] = True

        # %%
        # =============================================================================
        #     Extract allele and frequencies
        # =============================================================================
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
        
        ## get unfiltered calls per site 
        ## ancnti == all outgroup samples regardless of phylogeny; outsplnti == phylogenetically closest high-cov outgroup
        calls_outgroup = maNT[:,outgroup_bool] 
        outsplsnti = maNT[:,outgroup_spl_idx]
        ## for identification of good outgroup samples:
        polarized_site_file_txt = '' ## save number of polarized sites to store in file 
        for i in range(np.shape(outsplsnti)[1]): ## loop over all outgroup samples
            unpol_snps = np.sum(outsplsnti[:,i]==4)
            outgrp_splname = sampleNames[outgroup_spl_idx][i]
            polarized_site_file_txt = '\n'.join([polarized_site_file_txt, f'{outgrp_splname}:\t\t{unpol_snps} unpolarized SNPs'])
        ## select the outgroup which is phylogenetically the closest
        min_outgroupsplnti = np.argmin(np.sum(outsplsnti==4, axis = 0))
        outsplnti = outsplsnti[:,min_outgroupsplnti]
        outsplnti_m = np.tile(outsplnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele    
        ## get ancestral sequence from outgroups 
        outgroup_maNT = apy.major_allele(calls_outgroup) # NOTE: the filter criteria (cov,qual etc.) are not applied before major allele call
        ancnti = np.array([closest_outgroup_nti if closest_outgroup_nti != 4 else maNT_outgroups_nti for closest_outgroup_nti, maNT_outgroups_nti in zip(outsplnti, outgroup_maNT)])
        
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
        
        apy.save_num_muts_p_sample_num_mutsamples_per_p((calls[:,ingroup_bool] != refnti_m[:,ingroup_bool]) & (calls[:,ingroup_bool] < 4), sampleNames[ingroup_bool], p, contig_positions, suffix='_post_sitefilter')

        ## Filter per site across samples
        # Ignore here outgroup samples!
        calls_with_var_bool = (calls != refnti_m) & (calls < 4)
        siteFilt = np.any( ((calls[:,ingroup_bool]>3).sum(axis=1) >= ((num_samples-np.sum(outgroup_bool)) * filter_parameter_site_across_samples['max_fraction_ambigious_samples']),
                            # np.median( coverage[:,ingroup_bool], axis=1) < filter_parameter_site_across_samples['min_median_coverage_position'],
                            np.apply_along_axis(apy.median_above_zero, 1, coverage[:,ingroup_bool]) < filter_parameter_site_across_samples['min_median_coverage_position_with_cov']),
                            axis=0)
        calls[ siteFilt ,:] = 4 # sites that fail qc -> 4, for all samples incl. outgroup     
            
        # NOTE: func below takes forever with many SNPs...saved below
        [mutQual, mutQualIsolates] = apy.ana_mutation_quality(calls[:,ingroup_bool],quals[:,ingroup_bool]) # get FQ value for SNP across samples. mutQualIsolates contains sample indices for sample pair FQ based on. 
        mutQual = np.nan_to_num(mutQual, nan=-1) # turn mutQual nan's to -1; necessary to avoid later warning
        
        apy.save_num_muts_p_sample_num_mutsamples_per_p((calls[:,ingroup_bool] != refnti_m[:,ingroup_bool]) & (calls[:,ingroup_bool] < 4), sampleNames[ingroup_bool], p, contig_positions, suffix='_post_sitesmplfilter')

        ## NOTE: The removed sites by across samples filter migth overlap between each other/within sample filters --> should add the combinations!
        dict_of_per_smpl_filters = { 'Ns': (maNT == 4),
                                     'minQual': (quals < filter_parameter_site_per_sample['min_qual_for_call']), 
                                     'mAF': (maf < filter_parameter_site_per_sample['min_maf_for_call']), 
                                     'fwd': (coverage_forward_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call']), 
                                     'rev': (coverage_reverse_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call']), 
                                     'indel': (indels > (0.5*coverage) ) }
        hasnomutation_unfiltered = ~(maNT != refnti_m)
        intersection_counts_within = apy.combinatorial_filtering_to_goodpos_within_samples(hasnomutation_unfiltered, dict_of_per_smpl_filters, ingroup_bool)
        dict_of_across_smpl_filters = {
            tuple(dict_of_per_smpl_filters.keys()): np.sum(hasnomutation_unfiltered[:, ingroup_bool] | np.logical_or.reduce(list(dict_of_per_smpl_filters.values()))[:, ingroup_bool], axis=1) == np.sum(ingroup_bool), ## the sum of which sites have been removed by all within sample filters including unvar pos (hasnomut)
            'site_filter': siteFilt,
            'mutQual': (mutQual < 1).flatten()
            }
        intersection_counts_across = apy.combinatorial_filtering_to_goodpos_across_sites(dict_of_across_smpl_filters, np.shape(hasnomutation_unfiltered)[0])
        intersection_counts = intersection_counts_within
        intersection_counts.update(intersection_counts_across)
        matplotlib.use('agg')
        apy.upsetplot_filtering(intersection_counts, timestamp, subject_fld_label, refgenome)
        
        ## translate filtered calls of ingroup into goodpos. mutations we believe. fixedmutation part removed in v6.
        hasmutation = (calls != refnti_m) & (calls < 4) & (np.tile(mutQual,(1,num_samples)) >= 1) # consider only ingroup samples; mutQual >= 1 is very loose. Important filter with low qual data! refnt not ancnt (as we want to find snps _between_ samples not between samples and outgroup!)
        hasmutation[:,outgroup_bool] = False # put outgroup samples 4 in order to identify ingroup mutations only
        
        candpos = np.where( np.sum(hasmutation, axis=1) > 0 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!

        print('======================================================')
        print(f'Identified {len(candpos)} candidate positions for SNPs')
        print('======================================================\n')
        
        # %%
        # =============================================================================
        #  Find positions with candidate indels
        # =============================================================================
        ############
        ## filter indels 
        ############
        ## indel has high enough coverage on position
        pos_has_cov = indel_depth[:, :, 2] >= filter_parameter_indel_per_sample['min_cov_at_pos']
        pos_has_len = indel_size != 0 ## Indel need to generate sequence longer or shorter than before, else it is a multi SNV site and should be called within SNV part!
        ## indel has indel allele frequency above threshold
        pos_is_alt = indel_depth[:, :, 1] > (indel_depth[:, :, 2]*filter_parameter_indel_per_sample['min_iaf_for_call'])
        pos_is_ref = indel_depth[:, :, 0] > indel_depth[:, :, 2]*filter_parameter_indel_per_sample['min_iaf_for_call'] ## sample is ref coverage
        ## genotype likelihood difference of ref-indel is below threshold --> the smaller the higher the support for the indel 
        indel_has_support = (indel_support[:, :, 0] - indel_support[:, :, 1]) <= filter_parameter_indel_per_sample['max_indel_support'] ## the smaller the value the better the support for the best alternative indel
        ref_has_support = (indel_support[:, :, 1] - indel_support[:, :, 0]) <= filter_parameter_indel_per_sample['max_indel_support'] ## the smaller the value the better the support for the best alternative indel

        ## create matrix which supports indels and remove invariant positions (all samples have indel)
        has_indel = (pos_has_cov & pos_has_len & pos_is_alt & indel_has_support)
        has_ref = (pos_has_cov & pos_has_len & pos_is_ref & ref_has_support)
        
        ## remove sites which have support for ref and alt (e.g. transposons, tRNAs,...) --> misalignments 
        num_samples_w_intermediat_af = np.sum(~(pos_is_alt[:, ~outgroup_bool]) & ~(pos_is_ref[:, ~outgroup_bool]) & pos_has_cov[:, ~outgroup_bool], axis = 1) ## positions which have sufficient coverage, but are neither ref nor alt
        high_variablitity_pos = num_samples_w_intermediat_af >= (np.sum(pos_has_cov[:, ~outgroup_bool], axis = 1) * filter_parameter_indel_per_sample['max_fraction_ambigious_samples'])
        has_indel[high_variablitity_pos, :] = False ## set all invar pos to false
        
        ## remove sites which are eiter alt or ref only
        variable_pos = (np.sum(has_indel[:, ~outgroup_bool], axis = 1) > 0) & (np.sum(has_ref[:, ~outgroup_bool], axis = 1) > 0) ## Remove invariant sites of ingroup samples (one sample must have ref support and one alt support!)
        has_indel[~variable_pos, :] = False ## set all invar pos to false

        cand_indel = np.where( np.sum(has_indel, axis=1) > 0 )[0]
        print('======================================================')
        print(f'Identified {len(cand_indel)} candidate positions for Indels')
        print('======================================================\n')

        ## extract major allele
        maIndel = np.where(has_indel[cand_indel, :], indel_call[cand_indel, :, 1], ## set all indels if indel identified
                            np.where(has_ref[cand_indel, :], indel_call[cand_indel, :, 0], '?')) ## set all reference calls if high allele frequency, else set to undefined call
        indel_ref = np.full((len(cand_indel)), '', dtype = 'object')
        for i, indel_idx in enumerate(cand_indel):
            indel_ref[i] = np.unique([ref for ref in indel_call[indel_idx, :, 0] if (ref == ref) and (ref != None)])[0]

        # %%
        # =============================================================================
        #  Check if there are any candidate mutations identified
        # =============================================================================
        if (len(candpos) == 0) & (len(cand_indel) == 0):
            with open(f'{timestamp}_no_candpos_identified.txt', 'w') as fid:
                print('No variable positions identified. End of analysis.')
            #continue
        elif (len(cand_indel) == 0):
            with open(f'{timestamp}_no_candpos_identified.txt', 'w') as fid:
                print('No variable Indel positions identified. End of indel analysis.')
        elif (len(candpos) == 0):
            with open(f'{timestamp}_no_candpos_identified.txt', 'w') as fid:
                print('No variable SNV positions identified. End of SNV analysis.')

        # %%
        # =============================================================================
        #  Check for recombination in p and remove positions from goodpos
        # =============================================================================
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

        # %%
        # =============================================================================
        #  Save all goodpos/good_indels to new cmt 
        # =============================================================================
        ## save filtered data to npz file 
        save_goodpos_cmt_dict = {'calls': calls[goodpos, :], 
                                 'maNT': maNT[goodpos, :], 
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
                                 'indels': indels[goodpos,:],
                                 'indel_p': indel_p[cand_indel],
                                 'indel_call': indel_call[cand_indel, :, :],
                                 'indel_size': indel_size[cand_indel, :],
                                 'indel_depth': indel_depth[cand_indel, :, :],
                                 'indel_support': indel_support[cand_indel, :, :]}
        np.savez_compressed(f'{timestamp}_goodpos_candindel_cmt.npz', **save_goodpos_cmt_dict)

        # %%
        # =============================================================================
        #  Make some QC plots of SNVs
        # =============================================================================
        ## need to be set to plot in new window!!!
        if np.shape(goodpos)[0] > 0: 
            matplotlib.use('Qt5Agg')
            apy.plot_interactive_scatter_barplots(p[goodpos],mutQual[goodpos],'pos','qual', sampleNames,counts[:,:,goodpos],timestamp, filter_parameter_site_across_samples['corr_threshold_recombination'], subject_fld_label, refgenome, saveplots = True)

            ## need to be set to ensure saving of plot in right size without blocking a figure window
            matplotlib.use('agg')
            apy.plot_minorAF_rates(minorAF[goodpos, :], hasmutation[goodpos, :], sampleNames, subject_fld_label, refgenome, timestamp, '_goodpos')
            apy.plot_quals_rates(quals[goodpos, :], hasmutation[goodpos, :], sampleNames, subject_fld_label, refgenome, timestamp, '_goodpos')
            apy.plot_count_rates(counts[:, :, goodpos], hasmutation[goodpos, :], sampleNames, subject_fld_label, refgenome, timestamp, '_goodpos')
            apy.plot_Ns_rate(maNT[goodpos, :], sampleNames, p[goodpos], genomeLength, 10, subject_fld_label, refgenome, timestamp, '_goodpos')
        
        print(goodpos.size,'goodpos found.')


        # %% 
        # =============================================================================
        #  Define ancnti and outsplnti based on filtered calls
        # =============================================================================

        ## ancnti == all outgroup samples regardless of phylogeny
        calls_outgroup = calls[:,outgroup_bool] 
        ## NOTE: pushed down and we select the outgroup closest to the sample which polarizes the most snps, unpolarized snps will be exchanged by major allele of all!
        
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
                
        ## get ancestral sequence from outgroups 
        print(f'Checking outgroup samples (n = {np.shape(calls_outgroup)[1]}) for major allele:')
        outgroup_maNT = apy.major_allele(calls_outgroup) # NOTE: the filter criteria (cov,qual etc.) are not applied before major allele call
        ancnti = np.array([closest_outgroup_nti if closest_outgroup_nti != 4 else maNT_outgroups_nti for closest_outgroup_nti, maNT_outgroups_nti in zip(outsplnti, outgroup_maNT)])
        ancnti_m = np.tile(ancnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele    

        ## out chimera. based on outsplnti but all NA (==4) are replaced by major allele call in ancnti
        # we use this for molecular clock and derived allele freq (but not tree!)
        outchimerasplancnti = np.array([ancnti[i] if nti==4 else nti for i,nti in enumerate(outsplnti) ])
        polarized_site_file_txt = '\n'.join([polarized_site_file_txt, f'Unresolved goodpos in chimera outgroup:\t\t{np.sum(outchimerasplancnti[goodpos]==4)}'])
        polarized_site_file_txt = '\n'.join([polarized_site_file_txt, f'Resolved goodpos in chimera outgroup:\t\t{len(goodpos) - np.sum(outchimerasplancnti[goodpos]==4)} ({((len(goodpos) - np.sum(outchimerasplancnti[goodpos]==4)) / len(goodpos)) * 100}%)'])
        print(polarized_site_file_txt)

        with open(f'{timestamp}_polarized_sites_by_outgroup.txt', 'w') as fo_pol:
            fo_pol.write(polarized_site_file_txt)

        apy.plot_outgroup_snv_polarization(p[goodpos], sampleNames[outgroup_spl_idx], outsplsnti, outchimerasplancnti, goodpos, chrStarts, genomeLength, refgenome, f'{timestamp}_{subject_fld_label}')

        # %%
        # =============================================================================
        #  Make table with annotations
        # =============================================================================
        
        ## SNVs
        if len(goodpos) > 0:
            annotation_mutations = apy.annotate_mutations_v2(annotation_genes , p[goodpos] , refnti_m[goodpos,:] , ancnti_m[goodpos,:], calls[goodpos,:] , counts[:,:,goodpos] , hasmutation[goodpos,:], mutQual[goodpos,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP
            mutational_priorities = ['N', 'S', 'R', 'P', 'I', 'U']
            annotation_mutations_dedup = apy.deduplicate_annotation_mutations(annotation_mutations, timestamp, mutational_priorities)
            pd.concat(annotation_genes).to_csv(timestamp + '_annotation_genes.csv')
            
            annotation_mutations.to_csv(timestamp + '_annotation_SNVs.csv')
        
        ## Indels
        if len(cand_indel) > 0:
            annotation_indels = apy.annotate_indels(annotation_genes , indel_p[cand_indel], indel_call[cand_indel, :, :], indel_size[cand_indel, :], has_indel[cand_indel, :], indel_support[cand_indel, :, 0]-indel_support[cand_indel, :, 1], promotersize, ref_genome_folder)
            annotation_indels.to_csv(timestamp + '_annotation_indels.csv')
            ## write indel table
            apy.write_indel_table(annotation_indels, has_indel[cand_indel, :], indel_call[cand_indel, :, :], sampleNames)

        # %% 
        # =============================================================================
        #     Breakpoint: Too few positions pass filter
        # =============================================================================
        if len(goodpos) < 2:
            print("Too few positions after filter! >> skip: " + refgenome)
            with open(f'{timestamp}_unsufficient_num_goodpos.txt', 'w') as fid:
                print('Less than 2 goodpos identifed. End of analysis.')
            continue

        # %% 
        # =============================================================================
        #     Prepare calls for tree
        # =============================================================================    
        # define pos to use in tree
        goodpos2useTree = goodpos #[0::6] #(1:1000); %trim for easier tree view; TDL called quality_positions
        
        # get data and filter for goodpos
        calls_for_treei = calls[ np.ix_(goodpos2useTree, ingroup_wOutSpl_bool) ]
        calls_for_tree = apy.idx2nts(calls_for_treei) # ATCGN translation

        # add reference nucleotide for all positions
        refgenome_nts = apy.extract_outgroup_mutation_positions(ref_genome_folder, apy.p2chrpos(p[goodpos2useTree],chrStarts));
        calls_for_tree = np.concatenate((refgenome_nts[:, None],calls_for_tree),axis=1) # first column now refgenome_nts; refgenome_nts[:, None] to make ndims (2) same for both


        # %% 
        # =============================================================================
        #     Prepare sample Names for tree
        # =============================================================================    
        # build sampleNames w/ metainfo
        seq_types = dict()
        with open(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/kraken2/5-srst2/sumSRST2_results.csv'), 'r') as fid:
            header = fid.readline()
            for line in fid:
                line = line.strip().split(',')
                if (line[1] != 'L00071_07_01') and (refgenome != 'P07_Ehormaechei-c1_240825'): ## skip Enterobacter sample defined as P aeruginosa by kraken 
                    seq_types[line[0]] = line[1].replace('*', 'x') ## asterisks must be replaced for SNP trees (there asterisks separate label and nt and therefore those would mess up the count mutation function!)
        
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

        if 'refbased' in analysis_run:
            reference_label = refgenome[:5] + '_ref'
        elif 'denovo' in analysis_run:
            reference_label = refgenome.split('_')[1][:5] + '_ref'
        sampleNamesDnapars_wRef = np.append(reference_label,sampleNamesDnapars) # add name for outgroup // if this is changed it will cause issues with the per SNP tree as this has the Reference genome hardcoded!
        treesampleNamesLong_wRef = np.append(reference_label,treesampleNamesLong) # add name for outgroup // if this is changed it will cause issues with the per SNP tree as this has the Reference genome hardcoded!

        # %% 
        # =============================================================================
        #     Infer tree (parsimonious/ML/NJ)
        # =============================================================================    
        # ## build parsimony tree; added flag "buildTree": if false > only build dnaparse/fasta input files but do not calc tree 
        # Uncomment to get PS tree
        ## note: rerooting is just implemented for PS tree yet!
        if np.shape(calls_for_tree)[0] > 0:
            outgroup_for_root = treesampleNamesLong[outgroup_spl_idx][min_outgroupsplnti] ## note, reference is excluded to get right outgroup
            ## remove all other outgroups for calculations:
            redundant_outgroups_bool = np.invert(np.append(np.array(False), outgroup_bool)) + (treesampleNamesLong_wRef == outgroup_for_root)
            calls_for_tree_1outgrp = calls_for_tree[:, redundant_outgroups_bool]
            treesampleNamesLong_wRef_1outgrp = treesampleNamesLong_wRef[redundant_outgroups_bool]
            sampleNamesDnapars_wRef = sampleNamesDnapars_wRef[redundant_outgroups_bool]

            apy.generate_tree(calls_for_tree_1outgrp,treesampleNamesLong_wRef_1outgrp,sampleNamesDnapars_wRef,filetag=refgenome,buildTree='PS',root_on_smpl = outgroup_for_root, coloring_pattern = 'TP')
            
            ## infer sequence on each node (needed for molecular clock calculation!)
            apy.create_ancestral_reconstruction(f'{refgenome}_latest_rerooted_{outgroup_for_root}.nwk.tree',treesampleNamesLong_wRef_1outgrp,calls_for_tree_1outgrp, outdir="node_sequence_inference")

            ## get number of positions within core genome (NOTE: This takes some time as it handles the coverage matrix!)
            good_smpls_of_subj = np.array(subj_specific_samples)[goodsamples]
            ingrp_1outgrp_smpl = good_smpls_of_subj[redundant_outgroups_bool[1:]] ## remove from bool the reference!
            num_pos_coregenome, num_pos_invar = apy.get_num_core_genome_pos(cmtdir, ingrp_1outgrp_smpl, p, goodpos2useTree, genomeLength, filter_parameter_site_per_sample, filter_parameter_site_across_samples, median_coverage_position_with_cov = True)

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
        #     Make a tree for each Mutation location 
        # =============================================================================
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
        ## SOM SNP Table 
        apy.write_snp_table(calls[ np.ix_(goodpos2useTree, ingroup_bool) ], treesampleNamesLong[ingroup_bool], annotation_mutations_dedup)

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
        ## Note: This is like the clickable plot but slower as it runs on own cpu! uncomment if you want to run it!
        apy.plot_coverage_fwd_rev_stacked(chr_pos_gp,annotation_mutations_dedup,lod_fwd_cov,lod_rev_cov,timestamp, f'{subject_fld_label}_{refgenome}')

        # %%
        # =============================================================================
        # Estimate substitution rate
        # =============================================================================
        allmuts = np.array(['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG'],dtype=object)
        allmuts_types = np.array([0,2,1,0,1,2,5,4,3,4,5,3])
        mutationalspectrum = [1/12] * 12 # uniform distribution. 
        
        ## get calls for MRCA 
        ## NOTE: Outgroups not necessarily represent MRCAs!!!
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
        # Read in number of positions covered across genome
        # =============================================================================
        # obtain isolate specific-count of covered based required to correct observed mutation count (ie. for molecular clock inference and paraevo-poisson)
        # NOTE: This data has to be generated by user as it dependes on the assembly build. Here I put a placeholder
        numBaseGenome_covThreshold = np.sum(coverage_stats[:, filter_parameter_site_across_samples['min_median_coverage_position_with_cov']:11], axis = 1) ## bases covered > 4 per sample to be consistent with filters (filtering sites ≥4, therefore normalize on the same!)

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
                'max_muts_per_gene_to_track': 100,
                'timestamp':timestamp
                }
        
        print('Hard fix for num_mutations --> might not always be correct (see homoplasies!)')
        print('need to test "homoplasy_simulation" function!')
        annotation_mutations_dedup['num_mutational_events'] = annotation_mutations_dedup['nts'].str.len()-1

        try:
            [res_cand_nummut,annotation_mutation_paraSignal] = apy.parallel_evo_module( goodpos, contig_positions , annotation_mutations_dedup , annotation_genes , parameters )
        except:
            res_cand_nummut = np.array([])
            annotation_mutation_paraSignal = pd.DataFrame()

        ## report number of mutational events across indel and snv table 
        annotation_indels['num_mutational_events'] = annotation_indels['indel_alt'].str.len()
        mut_table_indel = annotation_indels.explode('gene_num_global')
        gene_mut_dict = {}
        mut_table = pd.concat([annotation_mutations_dedup, mut_table_indel])
        cols_oi = ['gene_num', 'gene_num_global', 'chr', 'pos', 'product', 'gene', 'protein_id', 'strand', 'type', 'loc1', 'loc2', 'sequence', 'note', 'locustag', 'orthologtag', 'translation', 'gene1', 'locustag1', 'orthologtag1', 'product1', 'distance1', 'gene2', 'locustag2', 'orthologtag2', 'product2', 'distance2', 'num_mutational_events']
        mut_table = mut_table[cols_oi].sort_values(['chr', 'pos'])
        
        mut_table_genic = mut_table.groupby(['gene_num_global'], as_index = False).agg({
                                                        'num_mutational_events': 'sum', 
                                                        'gene': 'first', 
                                                        'product': 'first'
                                                    })
        mut_table_genic = mut_table_genic[mut_table_genic['gene_num_global'] > 0.5] ## remove non_genic 
        num_mutations_genic = mut_table_genic['num_mutational_events'].sum()
        num_muts_genome = mut_table['num_mutational_events'].sum()
        [expectedNumberGenesMultipleMutations, expectedNumberOfGenesWithNmutations] = apy.parallel_evolution_counting_and_simulation(num_muts_genome,num_mutations_genic, parameters['Min_num_mutations_cand'],parameters['Min_mutation_density_cand'],parameters['NumTrialsSim'],parameters['max_muts_per_gene_to_track'],chr_pos_gp , parameters['ref_genome_folder'],annotation_genes)
        mut_table_genic['num_simulations'] = parameters['NumTrialsSim']
        for num_mut in mut_table_genic['num_mutational_events'].unique():
            mut_table_genic.loc[mut_table_genic['num_mutational_events'] == num_mut, 'p_val'] = 1-sum(expectedNumberGenesMultipleMutations < num_mut)/len(expectedNumberGenesMultipleMutations)

        mut_table_genic_multiple = mut_table_genic[mut_table_genic['num_mutational_events'] > 1].sort_values('num_mutational_events', ascending = False)
        print(mut_table_genic_multiple.to_string(index = False))
        mut_table_genic_multiple.to_csv(f"parallel_evo_multimut_genes_combined_SNVindel_results_{refgenome}.csv", index = False)

        mybin = np.arange(max(np.unique(expectedNumberGenesMultipleMutations))+2)-0.5 # +2 in order to now loose last bin due to arange; -0.5 needed for bars align at center with x axis ticks
        max_on_x_axis=max(mut_table_genic['num_mutational_events'].max(),max(np.unique(expectedNumberGenesMultipleMutations)))
        multiple_mut_counts = mut_table_genic_multiple.groupby(['num_mutational_events', 'p_val'], as_index = False).size()

        import matplotlib.patches as mpatches
        import subprocess
        fig, ax = plt.subplots(figsize = (7.5, 5))
        ax.hist(expectedNumberGenesMultipleMutations,bins=mybin,rwidth=0.8,color='#607c8e', edgecolor = 'k', lw = 0.4)
        ax.set_xlim(-0.5,max(mut_table_genic['num_mutational_events'].max(),np.max(expectedNumberGenesMultipleMutations))*1.25)
        ax.set_xticks(np.arange(0,ax.get_xlim()[1], 1))
        ax.set_ylabel(f"Simulated counts, N={parameters['NumTrialsSim']}")
        ax.set_xlabel(f"Number of genes with {parameters['Min_num_mutations_cand']} or more mutations")
        plt.title(f"Expectation vs observation of multiple mutated genes\nwithin {refgenome}\nincluding SNVs and indels", fontsize = 14)
        ax.text(0.98, 0.96, f"min #mut: {parameters['Min_num_mutations_cand']}", fontsize=10,horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
        ax.text(0.98, 0.92, f"min density: {parameters['Min_mutation_density_cand']}", fontsize=10,horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
        for idx, (num_mut, p_val, num_occ) in enumerate(multiple_mut_counts.values[::-1]):
            for i in range(int(num_occ)):
                y_pos = ax.get_ylim()[1]*(0.01+(0.025*i))
                ax.scatter(x = int(num_mut), y = y_pos, edgecolor = 'k', lw = 0.4, s = 50, c = 'violet')
            ax.text(0.98, 0.88-(0.04*idx), f"P({int(num_mut)}) = {round(p_val,3)}", fontsize=10,horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
            if idx == 0:
                legend_h = [mpatches.Patch(facecolor='#607c8e', edgecolor='black', lw = 0.4, label='Simulation')]
                legend_h += [plt.Line2D([0], [0], marker='o', color='w', lw = 0.2, markeredgecolor = 'k', markerfacecolor='violet', markersize=8, label='Observation')]
                ax.legend(handles=legend_h, loc = 'center left', bbox_to_anchor = (1, 0.5))
        plt.tight_layout()
        subprocess.run(f"mkdir -p pdf/adaptive_evo/",shell=True)
        fig.savefig(f'pdf/adaptive_evo/{timestamp}_{refgenome}_exp_vs_obs_multi_SNVindel_genes.pdf')

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
            contig_paraevo_candmut = res_cand_nummut[i,1].split('_', 1)[0]
            locustag_paraevo_candmut = res_cand_nummut[i,1].split('_', 1)[1]
            idxGeneAll = np.where((annotation_mutations_dedup['locustag']==locustag_paraevo_candmut) & 
                                (annotation_mutations_dedup['chr'].astype(str)==contig_paraevo_candmut))[0] 
            idxGene = idxGeneAll[0] # only take idx of first instance of gene
            lengthGene = len(annotation_mutations_dedup['sequence'][idxGene])
            my_lambda = tot_num_mutations * (lengthGene/median_genome_length)
            ### poisson cdf for num-obs_mut or more
            obs_mut = res_cand_nummut[i,2]
            p_poisson_obs = 1-stats.poisson.cdf(obs_mut-1,my_lambda)
            ### record if below bonferroni corrected p for that patient: p/#genes
            # record all candidates
            para_gene_dict['Patient'] = subject_fld_label
            para_gene_dict['Contig'] = annotation_mutations_dedup['chr'][idxGene]
            para_gene_dict['Gene_Start'] = annotation_mutations_dedup['loc1'][idxGene]
            para_gene_dict['Gene_End'] = annotation_mutations_dedup['loc2'][idxGene]
            para_gene_dict['Gene_length'] = annotation_mutations_dedup['loc2'][idxGene]- annotation_mutations_dedup['loc1'][idxGene]
            para_gene_dict['Strand'] = annotation_mutations_dedup['strand'][idxGene]
            para_gene_dict['Gene_id_prokka'] = annotation_mutations_dedup['locustag'][idxGene]
            para_gene_dict['Gene_id'] = annotation_mutations_dedup['gene'][idxGene]
            para_gene_dict['Orthologtag'] = annotation_mutations_dedup['orthologtag'][idxGene]
            para_gene_dict['Product_description_prokka'] = annotation_mutations_dedup['product'][idxGene]
            para_gene_dict['Gene_id_orthologtag'] = annotation_mutations_dedup['orthologtag'][idxGene]
            para_gene_dict['Poisson_p_value'] = p_poisson_obs
            para_gene_dict['Bonferroni_p_value'] = p_bonf_corrected
            para_gene_dict['Numer_mutations'] = res_cand_nummut[i,2]
            para_gene_dict['ProbSim'] = res_cand_nummut[i,6]
            para_gene_dict['Mutations'] = [annotation_mutations_dedup['muts'][m][0] for m in idxGeneAll]
            para_gene_dict['Mutation_type'] = [annotation_mutations_dedup['type'][m] for m in idxGeneAll]
            para_gene_dict['Mutation_position'] = [annotation_mutations_dedup['pos'][m] for m in idxGeneAll]
            para_gene_dict['Translation'] = annotation_mutations_dedup['translation'][idxGene]
            para_gene_dict['Sequence'] = annotation_mutations_dedup['sequence'][idxGene]            
            lod_gene_paraevo.append(para_gene_dict)
        df_gene_paraevo = pd.DataFrame(lod_gene_paraevo)  # turn patient results to dataframe and save
        df_gene_paraevo.to_csv('para_evo_snp_table.csv', index = False)
        if not df_gene_paraevo.empty:
            allP_gene_paraevo[subject_fld_label] = df_gene_paraevo # ... and store
        
        if (len(treesampleNamesLong[ingroup_bool]) < 2) or (len(goodpos) == 0):
            continue

        # %%
        # =============================================================================
        #  Molecular Clock Module
        # =============================================================================
        # mutations polarized using outchimerasplanc which is based on phylogenetically close outgroup sample, with possible unknown bases (==4) filled up by major allele across all outgroups
        
        ## read in hopsital admission 
        redcap_data = f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/RedCap_Data/Proc/{subject_fld_label}/{subject_fld_label}_redcap_instance_timepoints.csv'
        if os.path.exists(os.path.expanduser(redcap_data)):
            redcap_tps_df = pd.read_csv(os.path.expanduser(redcap_data))
            try:
                hosp_admission_date = redcap_tps_df.loc[redcap_tps_df['entry_type'] == 'basedata_hospadm', 'date'].values[0]
                hosp_admission_date = hosp_admission_date.split(' ')[0] ## remove hours of day
                hosp_admission_date = datetime.datetime.strptime(hosp_admission_date, '%Y-%m-%d')
            except:
                print('No hospital admission date found!')    
                hosp_admission_date = np.nan
        else:
            print('No hospital admission data could have been read as the provided path does not exist')

        list_timepoint_idx = []
        list_diagnostic_idx = []
        list_routine_idx = []
        for tp_str in ['TP' + str(i).zfill(2) for i in unique(timepoints)]:
            if tp_str != 'TP00': ## skip reference
                match_timepoint = np.array([i for i,x in enumerate(treesampleNamesLong[ingroup_bool]) if re.search('_'+tp_str+'_',x)])
                list_timepoint_idx.append(match_timepoint)

                ## match routine separately 
                match_routine = np.array([i for i,x in enumerate(treesampleNamesLong[ingroup_bool]) if re.search('_'+tp_str+'_',x) and re.search('_R_',x)])
                list_routine_idx.append(match_routine)

        outchimerasplancnti_gp = outchimerasplancnti[goodpos2useTree] # ATCGN translation; outchimerasplancnti is chimera of patient-specific outgroup and goodpos that remained unresolved (ie. 4) are substituted by major allele of all outgroup samples.
        numBaseGenome_covThresholdIngroup = numBaseGenome_covThreshold[ingroup_bool]
        patient_time = [datetime.datetime.strptime(date, '%Y%m%d') for date in kit_date[ingroup_bool]]
        patient_time_since_v1 = np.array([(date - min(patient_time)).days for date in patient_time])
        patient_time_since_admission = np.array([(date - hosp_admission_date).days for date in patient_time])


        # loop over each sample and assess #diffs to ancestral sequence
        exclude_region = False ## set to true if you want to exlde a region
        diff_count = np.zeros(timepoints[ingroup_bool].shape,dtype=float64)
        for sample_idx in range(len(treesampleNamesLong[ingroup_bool])):
            counter = 0
            for i,mrca_nt in enumerate( mrca_seq_ntidx ):
                if (mrca_nt != 4) and (calls_for_treei[:, ingroup_bool][i,sample_idx] != 4) and (mrca_nt != calls_for_treei[:, ingroup_bool][i,sample_idx]): ## sum(hasmutation[goodpos]) == ? or via calls? 
                    counter += 1 # counts nt distance sample and ancnt)
            diff_count[sample_idx] = counter

        # get counts 
        hdr = np.array(['hosp_days', 'colonization_days','rate','numBaseGenome'])
        num_mut_data_noinf = np.zeros( shape=(len(diff_count),4),dtype=float)
        num_mut_data_winf = np.zeros( shape=(len(diff_count),4),dtype=float)
        mean_count_per_timepoint_noinf = {}
        mean_count_per_timepoint_winf = {}
        for tp, r_tp in zip(list_timepoint_idx, list_routine_idx):
            if (r_tp.size != 0): ## only non infectious isolates (gastric samples excluded!)
                num_mut_data_noinf[r_tp,0] = patient_time_since_admission[r_tp]
                num_mut_data_noinf[r_tp,1] = patient_time_since_v1[r_tp]
                num_mut_data_noinf[r_tp,2] = diff_count[r_tp]
                num_mut_data_noinf[r_tp,3] = numBaseGenome_covThresholdIngroup[r_tp]
                mean_count_per_timepoint_noinf[np.unique(patient_time_since_v1[r_tp])[0]] = mean(num_mut_data_noinf[r_tp,2]/num_mut_data_noinf[r_tp,3])
            if (tp.size != 0): ## store including infection tp
                num_mut_data_winf[tp,0] = patient_time_since_admission[tp]
                num_mut_data_winf[tp,1] = patient_time_since_v1[tp]
                num_mut_data_winf[tp,2] = diff_count[tp]
                num_mut_data_winf[tp,3] = numBaseGenome_covThresholdIngroup[tp]
                mean_count_per_timepoint_winf[np.unique(patient_time_since_v1[tp])[0]] = mean(num_mut_data_winf[tp,2]/num_mut_data_winf[tp,3])
                
        ## remove entries from infection tps to save df 
        num_mut_data_noinf_save = num_mut_data_noinf[~np.all(num_mut_data_noinf == 0, axis=1)]
        
        num_samples_noinf = num_mut_data_noinf_save.shape[0]
        num_timepoints_noinf = len(np.unique(num_mut_data_noinf_save[:,1]))
        
        # save num_mut_data_noinf to csv and dict 
        np.savetxt(f'{timestamp}_mol_clock_data_noinf.csv', num_mut_data_noinf_save, header = ','.join(hdr), delimiter=",")
        np.savetxt(f'{timestamp}_mol_clock_data_winf.csv', num_mut_data_winf, header = ','.join(hdr), delimiter=",")
        mol_clock_data_dc[f'{subject_fld_label}_{refgenome}'] = num_mut_data_noinf_save # store data for subsequent regression analysis that includes correction for genome-wide base calls

        # %% PLOT seaborn
        # calc regression
        # Warning: make conditional for data with > 1 timepoint! adjust plotting function, too.
        if num_timepoints_noinf > 1:
            regress_dict = apy.calc_regression(np.column_stack((num_mut_data_noinf_save[:,1],num_mut_data_noinf_save[:,2]/num_mut_data_noinf_save[:,3])),subject_fld_label,num_samples_noinf,num_timepoints_noinf)
            regression_res_lod[f'{subject_fld_label}_{refgenome}'] = regress_dict # store, turn to pd dataframe for all subjects

            basescale = 1e6 # rescale genome-wide molecular rate to human-readable numbers
            apy.plot_molecular_clock(np.column_stack( (num_mut_data_noinf_save[:,1],(num_mut_data_noinf_save[:,2]/num_mut_data_noinf_save[:,3]) )),hdr[1:3],mean_count_per_timepoint_noinf,regress_dict,basescale,f'{subject_fld_label} {refgenome}', genome_size = genomeLength, jitter_x = 0.4, plot_confidence_interval=True, plot_predictor_interval=False) # plot in pdf/molclock/

            ## calculate tMRCA 
            ## select all tps with >= 10 isolates 
            confidence=0.95
            req_isolate_num = 10

            mut_rate_obs = regress_dict['slope'] * 365 * basescale ## mutation rate based on regression over all isolates 
            mut_rate_dict = {}
            species = refgenome.split('-')[0].split('_')[1]
            
            if not any([len(tp_idx) >= req_isolate_num for tp_idx in list_routine_idx]):
                print(f'Did not found any tp with more than {req_isolate_num} isolates')
                req_isolate_num = max([len(tp_idx) for tp_idx in list_routine_idx])
                print(f'Changed to {req_isolate_num} isolates which is the maximum present')
            if req_isolate_num < 2:
                print('with less than 2 isolates the calculation will be skipped')
                continue
            for tp, tp_idx in enumerate(list_routine_idx):
                if len(tp_idx) >= req_isolate_num:
                    scaled_num_muts = num_mut_data_noinf[tp_idx,2]/num_mut_data_noinf[tp_idx,3] ## note reference to old matrix with inf as 0 in matrix to ensure proper indexing!
                    # scaled_num_muts_perday = scaled_num_muts / num_mut_data_noinf[tp_idx,0]
                    mean_mut = np.mean(scaled_num_muts)
                    std_err = stats.sem(scaled_num_muts)
                    moe = std_err * stats.t.ppf((1 + confidence) / 2., len(tp_idx)-1)
                    ci95 = [mean_mut - moe, mean_mut + moe]
                    tmrca_pat_spec_mutrate = mean_mut * 365 * basescale / mut_rate_obs
                    ci95_tmrca_pat_spec_mutrate = [ci_bound * 365 * basescale / mut_rate_obs for ci_bound in ci95]

                    ## get how many days since hosp admission to tp:
                    d_since_hospadmission = np.unique(patient_time_since_admission[tp_idx])[0]
                    tmrca_pre_hospital = d_since_hospadmission - tmrca_pat_spec_mutrate 
                    ci95_tmrca_pre_hospital = [d_since_hospadmission - ci_bound for ci_bound in ci95_tmrca_pat_spec_mutrate]

                    mut_rate_dict[tp + 1] = {'days_in_hosp': np.unique(num_mut_data_noinf[tp_idx,0])[0],
                                            'days_to_root': np.unique(num_mut_data_noinf[tp_idx,1])[0],
                                            'num_isolates': len(tp_idx),
                                            'mut_rate_obs': mut_rate_obs,
                                            'mean_mut': mean_mut, 
                                            'ci95_mean_mut': ci95,
                                            'tmrca_pat_spec': tmrca_pat_spec_mutrate,
                                            'ci95_tmrca_pat_spec': ci95_tmrca_pat_spec_mutrate,
                                            'tmrca_pre_hospadmission': tmrca_pre_hospital,
                                            'ci95_tmrca_pre_hospadmission': ci95_tmrca_pre_hospital}
            mutation_rate_df = pd.DataFrame(mut_rate_dict).T.reset_index()
            mutation_rate_df.rename(columns={'index':'timepoint'})
            mutation_rate_df.to_csv(f'{subject_fld_label}_{refgenome}_mutation_rate.csv')
