## Script to summarise Bracken results 

######################################################
## load modules
######################################################

import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

######################################################
## positional and optional argument parser
######################################################

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''Generating summary of Bracken results''')

parser.add_argument("-b", dest="brackenres", help="Summarized bracken results with species resolution.", default = '4-bracken/sumBrackTest_species.txt')
parser.add_argument("-s", dest="samplecsv", help="Path to sample.csv", default = 'samples.csv')
parser.add_argument("-c", dest="cutoffPercRead", help="Samples need at least cutoffPercRead % reads assigned to species of interest to be included in filtered output.", default = 80, type=int)
parser.add_argument("-g", dest="refgenomedir", help="Folder to reference genomes with trailing '/' which has a subdirectory for each species including a taxonomic_rank.txt file")

args = parser.parse_args()

######################################################
## Sanity check and output dir generation
######################################################

## Sanity check of reference genome directory
refgenomedir = args.refgenomedir
if refgenomedir[-1] != '/':
    refgenomedir = refgenomedir + '/'

## generate outdirs if not existent
outdirs = ['pdf/target_spec_histogram', 'pdf/contaminated_samples', 'pdf/spec_vs_taxnodes', '3-samplesFiltPerSubject']
[os.makedirs(outdir, exist_ok=True) for outdir in outdirs]


######################################################
## Load sample csv
######################################################

## read sample csv which includes subject identifier
samplecsv = pd.read_csv(args.samplecsv)

## List of patterns to search for control wells
ctr_sample_names = ['ctr_well', 'ctr-well', 'ctr well', 'ctr_sample', 'ctr-sample', 'ctr sample', 
                    'ctr_smpl', 'ctr-smpl', 'ctr smpl', 'control_well', 'control-well', 'control well', 
                    'control_sample', 'control-sample', 'control sample', 'control_smpl', 'control-smpl', 
                    'control smpl']
ctr_sample_names_regex = '|'.join(ctr_sample_names)


## separate control data and samplescsv data 
ctr_samples = samplecsv[samplecsv['Sample'].str.contains(ctr_sample_names_regex, regex = True)]
samples = samplecsv[~samplecsv['Sample'].str.contains(ctr_sample_names_regex, regex = True)]

######################################################
## Load taxonomic ranks
######################################################

## get taxonomic_rank files for all species in samplecsv
unique_refgenomes = samples['ReferenceGenome'].unique()
taxonomic_ranks = pd.DataFrame()
for refgenome in unique_refgenomes:
    ## load tax rank file 
    taxfile = f'{refgenomedir}{refgenome}/taxonomic_rank.txt'
    if os.path.exists(taxfile):
        tmp_df = pd.read_csv(taxfile)
        taxonomic_ranks = pd.concat([taxonomic_ranks, tmp_df])
    else:
        SystemError(f'Reference genome directory or taxonomic_rank.txt file is missing! Error occured for {refgenome}.')

## clean taxonomic rank file 
taxonomic_ranks = taxonomic_ranks.drop('kingdom', axis = 1)
tax_levels = taxonomic_ranks.columns[1:]
taxonomic_ranks[tax_levels] = taxonomic_ranks[tax_levels].replace(' ', '_', regex = True) ## remove whitespaces across all taxlevels

## add taxonomic rank information to sample csv
samples_taxinfo = pd.merge(samples, taxonomic_ranks, how = 'left', left_on = 'ReferenceGenome', right_on = 'ref_genome_dir')

######################################################
## Load bracken results
######################################################

bracken_results = pd.read_csv(args.brackenres, sep = ' ', header = None)
bracken_results.columns = ['frag_perc_assigned', 'frag_num_assigned', 'species', 'Sample']

## subset info for control wells, if ctr well df is not empty and select first 3 entries per sample
if not ctr_samples.empty:
    ctr_hits = bracken_results.loc[bracken_results['Sample'].isin(ctr_samples['Sample'])]
    ctr_tophits = ctr_hits.sort_values('frag_perc_assigned', ascending=False).groupby('Sample').head(3)
else:
    ctr_tophits = pd.DataFrame() 


######################################################
## Filter samples for species-of-interest 
######################################################

## Get IDs of samples that fullfill cutoff 
## > output one file per Subject ID that includes seq IDs which fullfill quality criteria (ie: %reads on species >= cutoffPercRead )

## select rows with species which has been expected for a given sample
b_results_expected = pd.merge(bracken_results, samples_taxinfo, how = 'inner', on = ['Sample', 'species']) 
b_results_expected.loc[b_results_expected['Subject'].isnull(), 'Subject'] = 'unknown'

## get ratios and counts of species identified with assignment rates above given threshold 
expected_observed_spec_ratio = b_results_expected.groupby(['Subject', 'species']).agg(ratio=('frag_perc_assigned', lambda x: sum(x >= args.cutoffPercRead) / len(x)),
                                                                                      cnt=('frag_perc_assigned', lambda x: sum(x >= args.cutoffPercRead)))
expected_observed_spec_ratio = expected_observed_spec_ratio.sort_values(['Subject', 'species']).reset_index()
## order by subject and species for plotting 
expected_observed_spec_ratio['subj_spec'] = expected_observed_spec_ratio['Subject'] + '\n' + expected_observed_spec_ratio['species']
expected_observed_spec_ratio['subj_spec'] = pd.Categorical(expected_observed_spec_ratio['subj_spec'], 
                                                            expected_observed_spec_ratio.sort_values('subj_spec')['subj_spec'].unique(), 
                                                            ordered = True)

## filter for samples which are above threshold 
b_results_expected_thr = b_results_expected[b_results_expected['frag_perc_assigned'] >= args.cutoffPercRead]
unique_subjects = b_results_expected_thr['Subject'].unique()

## get all samples which have been failed given the cutoff
smpls_failed = pd.merge(bracken_results, samples_taxinfo[['Sample', 'Subject', 'species']], how = 'left', on = ['Sample'], suffixes = ['', '_exp']) 
low_classification_freq = ( (smpls_failed['frag_perc_assigned'] < args.cutoffPercRead) & (smpls_failed['species'] == smpls_failed['species_exp']) ) 
unexp_classification = ( (smpls_failed['species'] != smpls_failed['species_exp']) ) 
ctrls = ( smpls_failed['Sample'].isin(ctr_samples['Sample']) )  
smpls_failed = smpls_failed.loc[(low_classification_freq | unexp_classification) & ~ctrls]
smpls_failed_tophits = smpls_failed.sort_values('frag_perc_assigned', ascending=False).groupby('Sample').head(3)
smpls_failed_tophits.columns = ['frag_perc_assigned', 'frag_num_assigned', 'obs_species', 'Sample', 'Subject', 'expected_spec']

######################################################
## Save files
######################################################

# write sample IDs that had sufficient %reads on target into per Subject file. NOTE: subjects with 0 successful samples will not be saved!
for subject in unique_subjects:
    tmp_subj_df = b_results_expected_thr[b_results_expected_thr['Subject'] == subject]
    unique_species_in_subj = tmp_subj_df['species'].unique()
    for subj_species in unique_species_in_subj:
        ## Subselect the df for samples per patient and per species which is expected (and true) and save
        smpl_subj_spec_hit = tmp_subj_df.loc[tmp_subj_df['species'] == subj_species, 'Sample']
        ## generate a shorter species name i.e. Escherichia_coli get Ecoli (first letter of genus name is extracted and merged with species name)
        spec_abbrev = subj_species[0] + subj_species.split('_', 2)[1]
        smpl_subj_spec_hit.to_csv(f'3-samplesFiltPerSubject/subject{subject}_{spec_abbrev}_samples.txt', 
                                    sep='\t', index=False, header=False, quoting=0)

        ## print some statements for logging:
        num_input_samples_subj_spec = len(samples_taxinfo[(samples_taxinfo['species'] == subj_species) & (samples_taxinfo['Subject'] == subject)])
        percentage_kept = round((len(smpl_subj_spec_hit) / num_input_samples_subj_spec) * 100, 1)
        print(f'IMPORTANT: {len(smpl_subj_spec_hit)}/{num_input_samples_subj_spec} ({percentage_kept}%) samples are positive for {subj_species} and exceed >= {args.cutoffPercRead}% species-specific reads which will be used for assembly.')


# write file with all patients that fullfill criteria
with open('3-samplesFiltPerSubject/SubjectsPassBracken.txt', 'w') as fo:
    fo.write('\n'.join(unique_subjects))


## write failed samples and their contaminants
smpls_failed.to_csv('3-samplesFiltPerSubject/SamplesFailedBracken.txt', sep='\t', index = False)

## write top head of control wells if df exists 
if not ctr_tophits.empty:
    ctr_tophits.to_csv('3-samplesFiltPerSubject/ControlSamplesBracken.txt', sep='\t', index = False)


######################################################
## Plot summaries to assess quality
######################################################

####
## Barplot %samples that fullfilled %reads on target cutoff
plt_width = len(unique_subjects)/2+2
if plt_width < 8:
    plt_width = 8

fig, ax = plt.subplots(figsize = (plt_width, 8))
sns.barplot(data = expected_observed_spec_ratio, x = 'subj_spec', y = 'ratio', hue = 'species', dodge = False, alpha = 0.9, zorder = 3, ax = ax)
## annotate bars (note, subject and species are ordered )
idx = 0
sorted_patches = sorted(ax.patches, key=lambda p: p.get_x()) ## sort patches
for p in sorted_patches:
    ## skip all pathces which have no height (every color, defined by hue, has an x axis with a nan height value)
    if not np.isnan(p._height):
        cnt = expected_observed_spec_ratio.loc[idx, 'cnt']
        ax.annotate(f'n={str(cnt)}', (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha='center', va='bottom', rotation=45, fontsize=8)
        idx += 1
ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fancybox=False, shadow=False)
# Set plot title and axis labels
ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, rotation_mode = 'anchor', ha="right")
ax.set_xlabel("Subject")
ax.set_ylabel("%Samples usable for Assembly")
ax.set_title("%Samples w/ sufficient %reads-on-target per Subject (N top of bar)")
ax.grid(visible = True, axis = 'y', lw = 0.3, zorder = 0)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.tight_layout()
fig.savefig('pdf/BrackenAssess_percSamplePerSubject.pdf')
plt.close()

####
## Histograms of %reads on target species and homo_sapiens
def plot_hist_target_assignment(bracken_result_df, target_species, threshold):
    b_results_target_species = bracken_result_df[bracken_result_df['species'] == target_species]
    spec_abbr = target_species[0] + target_species.split('_', 2)[1]
    fig, ax = plt.subplots(figsize = (8, 6))
    sns.histplot(data = b_results_target_species, x = 'frag_perc_assigned', binwidth = 1, color = 'grey', ax = ax)
    ax.axvline(x = threshold, ls = '--', color = 'k', lw = 1)
    ax.annotate('classification\ncut-off', xy = (threshold - 3, ax.get_ylim()[1] * 0.825), rotation = 90, size = 8, horizontalalignment = 'center')
    ax.set_xlabel('% reads on tax node')
    ax.set_xlim(-1, 101)
    ax.set_title(f'% reads assigned per sample for species {target_species}')
    ax.grid(visible = True, axis = 'y', lw = 0.3, zorder = 0)
    plt.tight_layout()
    fig.savefig(f'pdf/target_spec_histogram/BrackenAssess_percReadHist_{spec_abbr}.pdf')
    plt.close()

## get all unique expected target species + Homo sapiens 
species_of_interest = list(samples_taxinfo['species'].dropna().unique())
## plot for each target species a own plot
for species in species_of_interest:
    plot_hist_target_assignment(b_results_expected, species, args.cutoffPercRead)
## plot human contamination per sample
plot_hist_target_assignment(bracken_results, 'Homo_sapiens', args.cutoffPercRead)

####
## plot which species is contaminated by which other 
def plot_failed_samples(failed_df, species):
    failed_target_species_df = failed_df[failed_df['expected_spec'] == species].copy()
    failed_target_species_df['cnt'] = failed_target_species_df.groupby('obs_species')['Sample'].transform('count') ## generate label with number of samples in which observed species has been identified
    failed_target_species_df['label'] = failed_target_species_df['obs_species']  + '\nn=' + failed_target_species_df['cnt'].astype(str)
    ## sort values for boxplot by median 
    failed_target_species_df['meds'] = failed_target_species_df.groupby('obs_species')['frag_perc_assigned'].transform('median')
    failed_target_species_df = failed_target_species_df.sort_values('meds', ascending = False).reset_index(drop = True) 

    spec_abbr = species[0] + species.split('_', 2)[1]
    
    ## plot boxplot 
    fig, ax = plt.subplots(figsize = (9, 6))
    sns.boxplot(data = failed_target_species_df, x = 'label', y = 'frag_perc_assigned', color= 'lightgrey', ax = ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, rotation_mode = 'anchor', ha="right")
    ax.set_ylim(-1, 101)
    ax.set_xlabel('Identified species')
    ax.set_ylabel('% reads on tax node')
    ax.set_title(f'Contamination of failed samples for expected {species}')
    ax.grid(visible = True, axis = 'y', lw = 0.3, zorder = 0)
    plt.tight_layout()
    fig.savefig(f'pdf/contaminated_samples/BrackenAssess_ContaminatingSpeciesOfFailed_boxplt_{spec_abbr}.pdf')
    plt.close()
    
    ## plot barplot with top 3 contaminating species sorted by most prominent contamination

    ## sort values for barplot by species which have highest assignment rate 
    failed_target_species_df = failed_target_species_df.sort_values(['frag_perc_assigned', 'cnt'], ascending = False)
    failed_target_species_df['Sample'] = pd.Categorical(failed_target_species_df['Sample'], failed_target_species_df['Sample'].unique(), ordered = True)

    plt_width = len(failed_target_species_df['Sample'].unique()) / 5
    if plt_width < 8:
        plt_width = 8
    fig, ax = plt.subplots(figsize = (plt_width, 6))
    sns.barplot(data = failed_target_species_df, x = 'Sample', y = 'frag_perc_assigned', hue = 'obs_species', dodge = False, ax = ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, rotation_mode = 'anchor', ha="right", size = 8)
    ax.set_ylim(-0, 100)
    ax.set_xlabel('Identified species')
    ax.set_ylabel('% reads on tax node')
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fancybox=False, shadow=False)
    ax.set_title(f'% read of top 3 contaminating species\nof failed samples for expected {species}')
    ax.grid(visible = True, axis = 'y', lw = 0.3, zorder = 0)
    plt.tight_layout()
    fig.savefig(f'pdf/contaminated_samples/BrackenAssess_ContaminatingSpeciesOfFailed_barplt_{spec_abbr}.pdf')
    plt.close()

uniq_expected_species_failed = smpls_failed_tophits['expected_spec'].dropna().unique()
for species in uniq_expected_species_failed:
    plot_failed_samples(smpls_failed_tophits, species)

####
## Scatter plot species vs. higher order tax levels

## loop over all tax nodes, exclude species taxlevel 
for tax_node in tax_levels:
    if tax_node in ['species', 'phylum']: ## we don't want to compare species to species and phylum has old naming scheme compared to taxonomic class file
        continue
    ## get bracken result for given tax node
    bracken_results_taxnode = pd.read_csv(args.brackenres.replace('_species.txt', f'_{tax_node}.txt'), sep = ' ', header = None)
    bracken_results_taxnode.columns = ['frag_perc_assigned', 'frag_num_assigned', 'tax_name', 'Sample']
    species_vs_other_taxnode = pd.merge(b_results_expected[['frag_perc_assigned', 'Sample', 'species', tax_node]], 
                                        bracken_results_taxnode[['frag_perc_assigned', 'Sample', 'tax_name']], 
                                        how = 'inner',
                                        left_on = ['Sample', tax_node],
                                        right_on = ['Sample', 'tax_name'],
                                        suffixes = ['_spec', '_othertaxnode'])

    species_vs_other_taxnode['label'] = species_vs_other_taxnode['species'] + '/\n' + species_vs_other_taxnode['tax_name']
    
    fig, ax = plt.subplots(figsize = (8, 6))
    sns.scatterplot(data = species_vs_other_taxnode, x = 'frag_perc_assigned_spec', y = 'frag_perc_assigned_othertaxnode', hue = 'label', alpha = 0.6, ax = ax)
    ax.set_xlabel('% reads assigned to species')
    ax.set_ylabel(f'% reads assigned to {tax_node}')
    ax.legend(title = f'species/{tax_node}', loc='center left', bbox_to_anchor=(1.02, 0.5), fancybox=False, shadow=False)
    ax.set_title(f'{tax_node} vs. species')
    ax.grid(visible = True, axis = 'both', lw = 0.3, zorder = 0)
    plt.tight_layout()
    fig.savefig(f'pdf/spec_vs_taxnodes/BrackenAssess_percRead_Species_vs_{tax_node}.pdf')
    plt.close()


