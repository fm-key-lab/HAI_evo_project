
import os
import numpy as np
import glob
import itertools
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import ScalarFormatter, FixedLocator
import seaborn as sns
from Bio import Phylo
import networkx as nx
from datetime import datetime

# Get current date and format it
current_date = datetime.now().strftime("%Y_%m_%d")

## 
## Updated to calculate the pairwise snp distance, the dMRCA (to root, not internal node as we want to have a measure how much comb like tree structre we have) 
# and the time difference as well as the time of the first isolate  

plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

## functions
def load_goodpos_cmt_npz_file(filepath):
    ## Load goodpos cmt npz file
    with np.load(os.path.expanduser(filepath), allow_pickle = True) as npfile:
        goodpos_cmt = {np_archive: npfile[np_archive] for np_archive in npfile}
    
    return goodpos_cmt

def get_sublists(df, label_order_in_tree=[], samples_analysed = [], lineage_cutoff=100):
    # filter based on given snp cutoff for lineages
    df_filtered = df[df['snp_diff'] <= lineage_cutoff]
    
    # create a graph and add edges
    G = nx.Graph()
    G.add_edges_from(zip(df_filtered['sample_1'], df_filtered['sample_2']))
    # find connected components (lineages)
    lineage_clusters = [list(component) for component in nx.connected_components(G)]
    samples_in_clusters = [sample for cluster in lineage_clusters for sample in cluster]

    # sort the lineages according to the order in tree
    if len(label_order_in_tree) >= 1:
        for sample in label_order_in_tree:
            if (sample not in samples_in_clusters) and (sample != ''):
                lineage_clusters.append([sample]) ## append missing samples which have a cluster by themselves
        sample_to_index = {sample: i for i, sample in enumerate(label_order_in_tree) if sample != ''}
        lineage_clusters = sorted(lineage_clusters, key=lambda group: sample_to_index[group[0]])
        lineage_clusters = [sorted(group, key=lambda s: sample_to_index[s]) for group in lineage_clusters]

    ## filter clusters to just those with samples (remove ref and outgroup)
    if len(samples_analysed) >=1:
        lineage_clusters = [cluster for cluster in lineage_clusters if any(item in samples_analysed for item in cluster)]
    
    return lineage_clusters

def read_tree(filepath, sampleNames):
    tree_paths = glob.glob(os.path.expanduser(filepath))
    if len(tree_paths) > 0: ## should just be one instance !
        if len(tree_paths) > 1:
            print('Mulitple trees have been identified. Please check!')
        mytree=sorted(tree_paths, key = lambda x: x.rsplit('/')[-1])[-1] ## extract latest tree
        tree = Phylo.read(mytree, "newick")
        tip_labels = np.array([i.name for i in tree.get_terminals()]) # tip labels ordered as tree
        tip_labels = ['_'.join([sl for sl in l.split('_')[5:] if not sl.startswith(('p', 'ST'))]) for l in tip_labels]
    else:
        print('No tree file was found. proceeding with the sampleNames as tiplabels')
        tip_labels = sampleNames.copy()
    return tip_labels

## get specimen log file
specimenLog_path = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/Drylab/metadata/*specimenlog.csv')
specimenLog_file = sorted(glob.glob(specimenLog_path), key = lambda x: x.rsplit('/')[-1])[-1]
specimenLog = pd.read_csv(specimenLog_file)

## get goodpos based CMT 
files = glob.glob(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/refbased_2024_08/analysis/P0*/*/2025-08-12*_goodpos_cmt.npz'))

## sanity check 
subject_species =  sorted([file.split('/')[-3] for file in files])
if len(subject_species) != len(set(subject_species)):
    SystemExit('There are duplicated files! Specify data you want to use for analysis!')

sample_name_dict = {}
snp_diff_df = pd.DataFrame()
snp_diff_pm_df = pd.DataFrame()
cum_sum_smpls = 0
for file in files:
    subject_species = file.split('/')[-3]
    goodpos_cmt_dict = load_goodpos_cmt_npz_file(file)
    calls = goodpos_cmt_dict['calls']
    sampleNames = goodpos_cmt_dict['sampleNames']
    kitdate = np.array([f'20{str(date)[:2]}-{str(date)[-4:-2]}-{str(date)[-2:]}' if date != 0 else '' for date in goodpos_cmt_dict['kit_date']], dtype=np.datetime64) ## convert to time
    ## exclude outgroups and ref if in cmt:
    ingroup_idx = [idx for idx, smpl in enumerate(sampleNames) if not '_GCA_' in smpl and not '_GCF_' in smpl and not '_ASM_' in smpl and not '_ref' in smpl]
    sampleNames = sampleNames[ingroup_idx]
    calls = calls[:, ingroup_idx]
    if 'P0010_Kmichiganensis' in subject_species:
        rm_gastric_idx = [idx for idx, smpl in enumerate(sampleNames) if not 'G00040_' in smpl]
        sampleNames = sampleNames[rm_gastric_idx]
        calls = calls[:, rm_gastric_idx]
    print(subject_species)
    print(f'Observed {len(sampleNames)} samples in data set')
    cum_sum_smpls += len(sampleNames)
    sample_name_dict[subject_species] = sampleNames

    ###############
    ## calculate pairwise SNP diff
    snp_diff_l = []
    smpl_i = []
    smpl_j = []
    t_diff = []

    for i, j in itertools.combinations(range(np.shape(calls)[1]), 2):
        diff_nt = (calls[:, i] != calls[:, j])
        no_undef_nt_i = (calls[:, i] != 4)
        no_undef_nt_j = (calls[:, j] != 4)
        snp_diff_l.append(np.sum(diff_nt & no_undef_nt_i & no_undef_nt_j))
        t_diff.append(abs((kitdate[i] - kitdate[j]) / np.timedelta64(1, 'D')))
        smpl_i.append(sampleNames[i])
        smpl_j.append(sampleNames[j])

    snp_diff_tmp_df = pd.DataFrame({'subject_species': [subject_species]*len(snp_diff_l), 'sample_1': smpl_i, 'sample_2': smpl_j, 'snp_diff': snp_diff_l, 'time_diff': t_diff})
    
    ###############
    ## identify pathogenic and nonpathogenic isolates
     
    ## get pathogenic and non pathogenic isolates
    ## exclude gastric juice samples as they are not pathogenic nor are considered in pilot study to colonize consistently the stomach!
    smpl_pathogenic = []
    smpl_mbiome = []

    for i, smpl in enumerate(sampleNames):
        if not '_GCA_' in smpl and not '_GCF_' in smpl and not '_ASM_' in smpl:
            if specimenLog.loc[specimenLog['ID'] == smpl, 'kit'].values[0][0] in ['B', 'L', 'U']:
                smpl_pathogenic.append(i)
            elif specimenLog.loc[specimenLog['ID'] == smpl, 'kit'].values[0][0] == 'K':
                smpl_mbiome.append(i)
            elif specimenLog.loc[specimenLog['ID'] == smpl, 'kit'].values[0][0] == 'G':
                print('Gastric juice sample. Excluding this sample for SNP differences in pathogen to mbiome pair')
            else:
                print(f'{smpl} not in specimenLog file')
    
    ###############
    ## add metadata about pathogenic/non-clinical isolates to df
    s1_is_pathogen = snp_diff_tmp_df['sample_1'].isin(sampleNames[smpl_pathogenic])
    s2_is_pathogen = snp_diff_tmp_df['sample_2'].isin(sampleNames[smpl_pathogenic])
    snp_diff_tmp_df['s1_type'] = np.where(s1_is_pathogen, 'pathogen', 'mbiome')
    snp_diff_tmp_df['s2_type'] = np.where(s2_is_pathogen, 'pathogen', 'mbiome')
    
    snp_diff_df = pd.concat([snp_diff_df, snp_diff_tmp_df])
    
    ###############
    ## calculate minimal SNP difference for all pathogenic vs all mbiome isolates

    snp_diff_pm_l_min = []
    smpl_p_min = []
    smpl_m_min = []
    t_diff_min = []

    for path_idx in smpl_pathogenic:
        
        snp_diff_pm_l = [] ## append all data, and extract for each pathogen then the minimum distance to any of the mbiome samples 
        smpl_m = [] ## append all data, and extract for each pathogen the sample name with minimum distance at the end
        for mbiome_idx in smpl_mbiome:
            diff_nt = (calls[:, path_idx] != calls[:, mbiome_idx])
            no_undef_nt_i = (calls[:, path_idx] != 4)
            no_undef_nt_j = (calls[:, mbiome_idx] != 4)
            snp_diff_pm_l.append(np.sum(diff_nt & no_undef_nt_i & no_undef_nt_j))
            smpl_m.append(sampleNames[mbiome_idx])
        ## get the id with the minimum distance
        if not subject_species.startswith(('P0007_Paeruginosa', 'P0021_Paeruginosa')):
            min_dist_idx = np.argmin(snp_diff_pm_l)
            snp_diff_pm_l_min.append(snp_diff_pm_l[min_dist_idx])
            t_diff_min.append(abs((kitdate[path_idx] - kitdate[min_dist_idx]) / np.timedelta64(1, 'D')))
            smpl_p_min.append(sampleNames[path_idx])
            smpl_m_min.append(smpl_m[min_dist_idx])
        else:
            snp_diff_pm_l_min.append(np.nan)
            t_diff_min.append(np.nan)
            smpl_p_min.append(np.nan)
            smpl_m_min.append(np.nan)

    tmp_pm_df = pd.DataFrame({'subject_species': [subject_species]*len(snp_diff_pm_l_min), 'sample_pathogen': smpl_p_min, 'sample_mbiome': smpl_m_min, 'snp_diff': snp_diff_pm_l_min, 'time_diff': t_diff_min})
    snp_diff_pm_df = pd.concat([snp_diff_pm_df, tmp_pm_df])
    
print(f'Calculated pairwise SNP difference on {cum_sum_smpls} samples in total.')

snp_diff_df = snp_diff_df.reset_index(drop = True)
snp_diff_pm_df = snp_diff_pm_df.reset_index(drop = True)
snp_diff_df['color'] = np.where(snp_diff_df['snp_diff'] < 100, 'close', 'distant')

## plt SNP differences between all isolates and on top all pathogenic isolates and their closest mbiome sample (across all pathogenic isolates!)
## NOTE: Plotting cannot illustrate well the isolates with a lot of difference!
name_include = ''

small_stepsize = 1
large_stepsize = 50
breakpoint_between_stepsizes = 100
max_of_bins = int((np.floor(max(snp_diff_df['snp_diff']) / large_stepsize)+1) * large_stepsize)
bins = np.concatenate([np.arange(0, breakpoint_between_stepsizes, small_stepsize), 
                       np.arange(breakpoint_between_stepsizes, max_of_bins, large_stepsize)]) ## add different bin sizes for the small (<= 100) and large SNV (> 100) differences

for legend in ['nolegend', 'wlegend']:
    figsize = (4.75, 3)
    fig, ax = plt.subplots(figsize = figsize)
    sns.histplot(data = snp_diff_df, x = 'snp_diff', bins = bins, color = 'lightgrey', alpha = 1,
                edgecolor = 'none', ax = ax)
    sns.histplot(data = snp_diff_pm_df, x = 'snp_diff', bins = bins, color = '#eb4824', alpha = 1, edgecolor = 'none', ax = ax)
    ax.set_yscale('symlog', linthresh=5)
    ax.set_xscale('symlog', linthresh=1)
    ax.axvline(x = 30, lw = 0.5, ls = '--', c = 'black')
    ax.set_xlim((0.75, max_of_bins*1.2))
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_formatter(ScalarFormatter())
        axis.get_major_formatter().set_scientific(False)
        axis.get_major_formatter().set_useOffset(False)
        max_val_axis = int(np.ceil(np.log10(axis.get_view_interval()[1]))) ## the max log base of the axis to set minor ticks in that range
        axis.set_minor_locator(FixedLocator([k * 10**n for n in range(max_val_axis) for k in range(1, 10)])) ## set minor ticks for all integer values (remove minor ticks < 1!)
    ax.set_xlabel('Intrahost pairwise SNV differences\nper species')
    if legend == 'wlegend':
        legend_info = [mpatches.Patch( facecolor='lightgrey', alpha = 0.6, edgecolor = 'none', label='All isolates per HAI')]
        legend_info += [mpatches.Patch( facecolor='#eb4824', edgecolor = 'none', label='Isolates from site of\ninfection to their closest\nmicrobiome isolate')]
        ax.legend(title='Pairwise SNV differences', handles = legend_info, loc='upper right', bbox_to_anchor=(1.02, 1.05), fancybox=False, shadow=False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2025_06_09_pairwise_SNP_diff_path_vs_mbiome_refbased_all_specs_log{name_include}_{legend}.pdf'))
    plt.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2025_06_09_pairwise_SNP_diff_path_vs_mbiome_refbased_all_specs_log{name_include}_{legend}.svg'))
        
## plot the distance within each timepoint vs across timepoints for pathogens
s1_is_pathogen = (snp_diff_df['s1_type'] == 'pathogen')
s2_is_pathogen = (snp_diff_df['s2_type'] == 'pathogen')
snp_diff_pathogen = snp_diff_df[s1_is_pathogen & s2_is_pathogen]
snp_diff_pathogen_collapsed = snp_diff_pathogen[['subject_species', 'snp_diff', 'time_diff']].groupby(['subject_species', 'snp_diff', 'time_diff']).size().reset_index(name="Number of pairwise comparisons")
snp_diff_pathogen_collapsed['Pathogen-Patient pair'] = snp_diff_pathogen_collapsed['subject_species'].str.split('_').str[:2].str.join(' ')
snp_diff_pathogen_collapsed['Pathogen-Patient pair'] = snp_diff_pathogen_collapsed['subject_species'].str.split('_').str[:2]
snp_diff_pathogen_collapsed['Pathogen-Patient pair'] = 'Patient ' + snp_diff_pathogen_collapsed['Pathogen-Patient pair'].str[0].str[1:].astype(int).astype(str) + ' ' + snp_diff_pathogen_collapsed['Pathogen-Patient pair'].str[1].str[0] + '. ' + snp_diff_pathogen_collapsed['Pathogen-Patient pair'].str[1].str[1:]
snp_diff_pathogen_collapsed['Pathogen-Patient pair'] = snp_diff_pathogen_collapsed['Pathogen-Patient pair'].str.replace('hormaecheiYT3', 'hormaechei')

snp_diff_all_collapsed = snp_diff_df[['subject_species', 'snp_diff', 'time_diff']].groupby(['subject_species', 'snp_diff', 'time_diff']).size().reset_index(name="Number of pairwise comparisons")
snp_diff_all_collapsed['Pathogen-Patient pair'] = snp_diff_all_collapsed['subject_species'].str.split('_').str[:2]
snp_diff_all_collapsed['Pathogen-Patient pair'] = 'Patient ' + snp_diff_all_collapsed['Pathogen-Patient pair'].str[0].str[1:].astype(int).astype(str) + ' ' + snp_diff_all_collapsed['Pathogen-Patient pair'].str[1].str[0] + '. ' + snp_diff_all_collapsed['Pathogen-Patient pair'].str[1].str[1:]
snp_diff_all_collapsed['Pathogen-Patient pair'] = snp_diff_all_collapsed['Pathogen-Patient pair'].str.replace('hormaecheiYT3', 'hormaechei')

uniq_pat_patho_pairs = sorted(snp_diff_pathogen_collapsed['Pathogen-Patient pair'].unique(), key = lambda x: int(x.split(' ')[1])) ## sort by the patient ID
num_cols = 3
num_rows = int(np.ceil(len(uniq_pat_patho_pairs)/num_cols))
min_num_pairwise_comp = 1
max_num_pairwise_comp = max(snp_diff_pathogen_collapsed['Number of pairwise comparisons'])
min_scatter_size = 20
max_scatter_size = 200
fig, axs = plt.subplots(figsize = (2.5*num_cols+1.5, 2.5*num_rows), ncols = num_cols, nrows = num_rows, sharex = True, sharey = True)
axs_l = axs.flatten()
for axid, (ax, ppp) in enumerate(zip(axs_l, uniq_pat_patho_pairs)):
    plt_df = snp_diff_pathogen_collapsed[snp_diff_pathogen_collapsed['Pathogen-Patient pair'] == ppp]
    sns.scatterplot(data = plt_df, 
                    x = 'time_diff', 
                    y = 'snp_diff', 
                    size = 'Number of pairwise comparisons', 
                    sizes=(min_scatter_size, max_scatter_size), 
                    size_norm=(min_num_pairwise_comp, max_num_pairwise_comp), 
                    color = 'grey',
                    edgecolor = 'k', 
                    lw = 0.25, 
                    alpha = 0.4, 
                    ax = ax)
    ax.axhline(y = 30, ls = '--', lw = 0.5, c = 'k')
    ax.annotate(text='Lineage\nthreshold', xy = (snp_diff_pathogen_collapsed['time_diff'].max(), 30), va = 'bottom', ha = 'right', fontsize = 9)
    ax.set_yscale('symlog', linthresh=5)
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.get_major_formatter().set_scientific(False)
    ax.yaxis.get_major_formatter().set_useOffset(False)
    max_val_axis = int(np.ceil(np.log10(ax.yaxis.get_view_interval()[1]))) ## the max log base of the axis to set minor ticks in that range
    ax.yaxis.set_minor_locator(FixedLocator([k * 10**n for n in range(max_val_axis) for k in range(1, 10)])) ## set minor ticks for all integer values (remove minor ticks < 1!)
    if axid >= (num_cols*num_rows)-num_cols:
        ax.set_xlabel('Days between isolates')
    if axid % num_cols == 0:
        ax.set_ylabel('Intrahost\npairwise SNV differences')
    if axid == (round(num_rows/2, 0) * num_cols)-1: 
        legend_sizes = [1, 40, 80, 120, 160, 200]
        legend_marker_sizes = np.interp(legend_sizes, (min_num_pairwise_comp, max_num_pairwise_comp), (min_scatter_size, max_scatter_size))
        legend_marker_sizes = np.sqrt(legend_marker_sizes) ## matplotlib uses diameter based scaling, while seaborn uses area based scaling
        legend_hanles = [plt.Line2D([0], [0], 
                                marker='o', 
                                color = 'grey',
                                markeredgecolor = 'k',
                                markeredgewidth = 0.25, 
                                linestyle = '',
                                alpha = 0.4,
                                markersize=s, 
                                label=l) 
                                for l, s in zip(legend_sizes, legend_marker_sizes)]
        ax.legend(title = 'Number of\npairwise comparisons\namongst isolates\nfrom site of infection', handles = legend_hanles, loc='center left', bbox_to_anchor=(1, 0.5), fancybox=False, shadow=False)
    else:
        ax.legend().remove()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title(f"{' '.join(ppp.split(' ')[:2])}\n{' '.join(ppp.split(' ')[2:])}")
plt.tight_layout()
plt.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2025_04_14_pairwise_SNP_diff_pathogens_vs_time_refbased_all_specs_log{name_include}.pdf'))

print('Unique SNP differences below 100 snps all vs all:')
print(snp_diff_df[snp_diff_df['snp_diff'] < 100].sort_values('snp_diff')['snp_diff'].unique())
snp_diff_df.to_csv(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/refbased_2024_08/analysis/2025_06_09_pairwise_SNP_diff.csv'), index = False)

print('Minimal SNP differences of pathogen species to closest microbiome isolate:')
print(snp_diff_pm_df.sort_values('snp_diff').groupby('subject_species')['snp_diff'].min())
print(snp_diff_pm_df.sort_values('snp_diff').groupby('subject_species')['snp_diff'].unique())
snp_diff_pm_df.to_csv(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/refbased_2024_08/analysis/2025_06_09_min_pairwise_SNP_diff_per_pathogenisolate_vs_mbiome.csv'), index = False)

print('Unique SNP differences of pathogens vs. closest microbiome isolate:')
min_pwSNP_diff_per_pathogenic_iso = snp_diff_pm_df.sort_values('snp_diff').groupby(['subject_species', 'snp_diff'])['sample_pathogen'].count().reset_index()
min_pwSNP_diff_per_pathogenic_iso.columns = ['subject_species', 'snp_diff', 'cnt']
print(min_pwSNP_diff_per_pathogenic_iso)
min_pwSNP_diff_per_pathogenic_iso.to_csv(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/refbased_2024_08/analysis/2025_06_09_min_pairwise_SNP_diff_per_pathogen_vs_mbiome.csv'), index = True)


########################################################################
########################################################################
########################################################################

## identify the lineages given a cutoff and write them to a file to use for the assembly
lineage_cutoff = 30
kraken_file_path = '~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/kraken2/3-samplesFiltPerSubject'

for file in files:
    analysis_path = '/'.join(file.split('/')[:-1])
    subject_species = file.split('/')[-3]
    subject_species_short = f'{subject_species[0]}{subject_species[3:5]}_{subject_species.split("_")[1]}'
    if 'P07_Ehorm' in subject_species_short:
        subject_species_short = subject_species_short.replace('YT3', '') ## remove for Ehormaechei the subspecies label

    tree_path = f'{analysis_path}/*_latest_rerooted_OG_*_TP0.nwk.tree'
    tree_labels = read_tree(tree_path, sample_name_dict[subject_species])
    lineage_cluster = get_sublists(snp_diff_df[snp_diff_df['subject_species'] == subject_species], tree_labels, sample_name_dict[subject_species], lineage_cutoff)

    ## load samples which have been identified to be good via bracken (cutoff >= 80%)
    good_samples = []
    with open(os.path.expanduser(f'{kraken_file_path}/subject{subject_species_short}_samples.txt'), 'r') as fid:
        for line in fid:
            good_samples.append(line.strip())
    
    ## save clusters to independent files 
    for lineage_id, lineage in enumerate(lineage_cluster):
        ## all samples of a cluster
        with open(os.path.expanduser(f'{analysis_path}/{current_date}_{subject_species_short}-c{lineage_id+1}.txt'), 'w') as fo:
            for sample in lineage:
                fo.write(f'{sample}\n')
        ## good samples of a cluster
        
        with open(os.path.expanduser(f'{analysis_path}/{current_date}_{subject_species_short}-c{lineage_id+1}_bracken80.txt'), 'w') as fo:
            for sample in lineage:
                if (sample in good_samples):
                    fo.write(f'{sample}\n')

## to inlcude ST in assembly names curate them manually as this is not yet implemented!
