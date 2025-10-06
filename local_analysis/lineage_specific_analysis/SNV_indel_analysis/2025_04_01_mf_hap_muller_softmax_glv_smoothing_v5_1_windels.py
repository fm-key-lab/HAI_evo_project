
############################################################################################################################################
############################################################################################################################################
########## Muller plots - allele freq change > 0.3
############################################################################################################################################
############################################################################################################################################

# =============================================================================
# %% Modules
# =============================================================================

import os
import sys
import glob
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from datetime import datetime
from Bio import Phylo
from scipy.interpolate import PchipInterpolator

# =============================================================================
# %% Functions
# =============================================================================

########################################
## get isolates with mutation of interest

def get_tips_w_mutation(mut_table, mut_table_isolate_cols, muller_tbl):
    ## get the childs of the branch 
    genotype_dict = {}
    for _, muller_row in muller_tbl.iterrows():
        gt = muller_row['Mutation'] # 'genotype-' + str(muller_row['Trajectory'])
        is_chr = (mut_table['chr'] == muller_row['Chromosome'])
        is_pos = (mut_table['pos'] == muller_row['Position'])
        mut_oi = mut_table[is_chr & is_pos]
        if not mut_oi.empty:
            isolates_w_gt = [col for col in mut_table_isolate_cols if all(mut_oi[col] == muller_row['Alt_allele'])]
            genotype_dict[gt] = isolates_w_gt
    return genotype_dict

def collapse_linked_mutations(muller_tbl, gt_to_isolate_dict, muller_table_columns_to_combine = []):
    ## if all isolates within the lineage_dict are the same, the mutations are linked and should be collapsed as we cannot polarize what occured first
    mutational_order = [mut for mut in muller_tbl['Mutation']]
    
    linked_mut = {}
    seen = set()
    for mut, isolates in gt_to_isolate_dict.items():
        sorted_isolates_tup = tuple(sorted(isolates))
        if sorted_isolates_tup not in seen:
            seen.add(sorted_isolates_tup)
            linked_mut[sorted_isolates_tup] = [mut]
        else:
            linked_mut[sorted_isolates_tup].append(mut)
    gt_to_isolate_dict_dedup = {'; '.join(sorted(muts, key=lambda x: mutational_order.index(x))): list(isolates) for isolates, muts in linked_mut.items()}
    
    ## deduplicate the muller table
    muller_tbl[muller_table_columns_to_combine] = muller_tbl[muller_table_columns_to_combine].astype(str)
    for mut_on_branch in gt_to_isolate_dict_dedup.keys():
        muts = mut_on_branch.split('; ')
        if len(muts) > 1:
            for col in muller_table_columns_to_combine:
                new_val = '; '.join([muller_tbl.loc[(muller_tbl['Mutation'] == mut), col].astype(str).values[0] for mut in muts])  ## keep order of annotation!
                muller_tbl.loc[(muller_tbl['Mutation'] == muts[0]), col] = new_val
            muller_tbl = muller_tbl[~muller_tbl['Mutation'].isin(muts[1:])]
    
    return muller_tbl, gt_to_isolate_dict_dedup

########################################
## rename locustags

def rename_locustags(mut_df, locustag_str = ''):
    for row_idx, row in mut_df.iterrows():
        new_gene_tag = ''
        if row['gene'] != row['gene']:
            continue
        for gene_idx, gene in enumerate(row['gene'].split(';')):
            if gene == '.':
                lt = row["locustag"].split(';')[gene_idx].split(':')[-1]
                lt = lt.replace(locustag_str, '').replace('_', '')
                new_gene_tag += f';LT{lt}'
            else:
                new_gene_tag += f';{gene}'
        mut_df.loc[row_idx, 'gene'] = new_gene_tag.strip('; ')
    return mut_df

########################################
## Tree functions

def load_tree(tree_oi_file):
    # Load the tree from a Newick file
    tree = Phylo.read(tree_oi_file, "newick")
    ## set unique internal node names
    for i, node in enumerate(tree.get_nonterminals()):
        if not node.name:  # Only assign if there's no existing name
            node.name = f"Internal_{i}"
    return tree

def find_mrca(tree, tips):
    # get MRCA for given list of tips labels
    return tree.common_ancestor(tips)

def get_finetuned_polarization(anc_to_tip_dict, der_to_anc_dict, lineage_dict):
    ## check if a lineage is a subset of another lineage to finetune polarization (ancestry -> derived state)
    anc_to_tip_finetuned_dict = {}
    for anc, der in anc_to_tip_dict.items():
        if len(der) > 1:
            for mut1, mut2 in combinations(der, 2):
                if (len(lineage_dict[mut1]) < len(lineage_dict[mut2])) & all([smpl in lineage_dict[mut2] for smpl in lineage_dict[mut1]]):  
                    anc_to_tip_finetuned_dict[mut2] = [mut1] ## add mutation as new child of ancestor 
                    anc_to_tip_finetuned_dict[anc] = [mut for mut in anc_to_tip_dict[anc] if mut != mut1] ## remove mutation from old ancestor
                    der_to_anc_dict[mut1] = mut2
                elif (len(lineage_dict[mut2]) < len(lineage_dict[mut1])) & all([smpl in lineage_dict[mut1] for smpl in lineage_dict[mut2]]):  
                    anc_to_tip_finetuned_dict[mut1] = [mut2] ## add mutation as new child of ancestor 
                    anc_to_tip_finetuned_dict[anc] = [mut for mut in anc_to_tip_dict[anc] if mut != mut2] ## remove mutation from old ancestor
                    der_to_anc_dict[mut2] = mut1
                else:
                    anc_to_tip_finetuned_dict[anc] = der
        else:
            anc_to_tip_finetuned_dict[anc] = der
    return der_to_anc_dict, anc_to_tip_finetuned_dict

def get_closest_anc_per_genotype_of_interest(tree, lineage_dict, root_name = 'root', finetune_polarization = True):
    # get MRCA nodes for each lineage
    mrca_nodes = {name: find_mrca(tree, tips) for name, tips in lineage_dict.items()}

    der_to_anc_dict = {}  

    for child, child_mrca in mrca_nodes.items():
        closest_parent = None
        closest_depth = -1  # variable for identifying depth of 
        for parent, parent_mrca in mrca_nodes.items():
            if (parent == child) | (parent_mrca.name == child_mrca.name): ## check that mutations are not the same (selfcomparisons) and the nodes are different between 
                continue
            if parent_mrca.is_parent_of(child_mrca):  # check if ancestor
                parent_depth = tree.distance(tree.root, parent_mrca)  # calculate phylogenetic distance between
                if parent_depth > closest_depth:  # keep closest ancestor
                    closest_parent = parent
                    closest_depth = parent_depth
                
        # assign closest node (or 'root' if none found)
        if child != root_name:
            der_to_anc_dict[child] = closest_parent if closest_parent else root_name
        anc_to_tip_dict = {}
        for key, value in der_to_anc_dict.items():
            if value not in anc_to_tip_dict:
                anc_to_tip_dict[value] = []
            anc_to_tip_dict[value].append(key)

    if finetune_polarization:
        ## check if a lineage is a subset of another lineage to finetune polarization (ancestry -> derived state)
        der_to_anc_dict, anc_to_tip_dict = get_finetuned_polarization(anc_to_tip_dict, der_to_anc_dict, lineage_dict)

    return der_to_anc_dict, anc_to_tip_dict

########################################
## Get allele frequencies at timepoint

def get_allele_freq_at_tp(alleles_i, timepoint_subsetting_bool, alt_nti, nt_arr = range(4), ):
    allelesi_at_timepoint = alleles_i[ timepoint_subsetting_bool ] ## extract allele for all samples from timepoint (or location, ...)
    allelei_covered = np.sum(allelesi_at_timepoint != 4)
    if allelei_covered != 0: # test denominator not 0 (all samples == N)
        allelei_is_alt = (allelesi_at_timepoint == alt_nti)
        allele_freq = np.array([np.sum(np.all(((allelesi_at_timepoint == nt_idx), allelei_is_alt ), axis=0)) / allelei_covered for nt_idx in nt_arr], dtype = float)
    else:
        allele_freq = np.array([0 for _ in nt_arr], dtype = float)
    isolate_cnts = sum(timepoint_subsetting_bool)
    return allele_freq, isolate_cnts

########################################
## Smoothing function

def softmax_interpolation(timepoints_array, frequencies, num_interp=100, groups_overall_frequencies = [], log_space = True, np_log_err = 1e-6):
    x_interp = np.linspace(timepoints_array[0], timepoints_array[-1] - 1, num_interp)

    interpolated = []
    ## loop over every genotype individually
    for i in range(frequencies.shape[1]):
        if log_space:
            ## smooth frequencies in log space and ensure 0 values are not causing errors by adding np_log_err
            f = PchipInterpolator(timepoints_array, np.log(frequencies[:, i] + np_log_err))
            ## extract from log space number of expected (=num_interp) values to original space
            interpolated.append(np.exp(f(x_interp)))
        else:
            ## smooth frequencies
            f = PchipInterpolator(timepoints_array, frequencies[:, i])
            ## extract number of expected (=num_interp) values
            interpolated.append(f(x_interp))
    
    interpolated = np.array(interpolated).T
    ## restrict frequencies to 0-1 space
    interpolated = np.clip(interpolated, 0, None)  # Ensure non-negative values
    if len(groups_overall_frequencies) > 0:
        remaining_frequency_of_parent = (groups_overall_frequencies-np.sum(interpolated, axis = 1).reshape(-1, 1))
        remaining_frequency_of_parent = np.clip(remaining_frequency_of_parent, 0, None)  # Ensure non-negative values
        print(remaining_frequency_of_parent.shape)
        interpolated_w_anc = np.concatenate([remaining_frequency_of_parent, interpolated], axis = 1) ## add remaining frequencies of ancestral 
        interpolated_w_anc /= np.sum(interpolated_w_anc, axis=1, keepdims=True) ## normalize to 1 
        interpolated_w_anc *= groups_overall_frequencies ## normalize to maximal frequency of ancestral at time
        return x_interp, interpolated_w_anc
    else:
        interpolated /= np.sum(interpolated, axis=1, keepdims=True)
        return x_interp, interpolated

# =============================================================================
# %% Parameters for interpolation
# =============================================================================

allele_freq_change_min = 0.3 # the minimum shift in observed allele freq between timepoints required for reporting
num_interp = 1000 ## number of steps across timeline to interploate
np_log_err = 1e-3 ## error to add to all values to enable np.log() from 0
root_name = 'basal'
locations_abbreviations = ['Nasal', 'Oral', 'Rectal', 'Skin'] #, 'Lung', 'Blood']
inf_locations = ['Lung', 'Blood', 'Urine']
# clean_background_on_empty_tps = False ## set to true if you want to remove root on timepoints with no isolates collected; --> deprecated
add_missing_timepoints = False ## set to true to add the empty timepoints where no isolate was observed, Note: Empty timepoints need to be listed in a dictionary per species with timepoints as list in `missing_timepoints`
missing_timepoints = {'P07_Ehormaechei-c1': [35] ## days post infection
}
drop_empty_tp = True ## Set to true if you want to drop timepoints with no isolates at location --> interpolation then ignores those as unobserved timepoints
convert_alleles_to_genotypes = True ## Set to true if you want to have genotypes instead of individual mutations in the final files (e.g. a genotype with 2 mutations will be named by both mutations!); NOTE: frequency shifts are still calculated on alleles only!

## NOTE: 
# Frequencies below the allele_freq_change_min are excluded from the plotting and are just shown as part of their basal genotype!

# =============================================================================
# %% Working directory and data
# =============================================================================
res_folder = os.path.expanduser("~/Documents/Data/Drylab/mf_2020_hap/muller_plots")
analysis_path = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/backup')#P0007_P07_Ehormaechei-c1_240825/2025_03_10_P07_Ehormaechei-c1_240825')
analysis_paths = glob.glob(f'{analysis_path}/P0*')

for species_path in sorted(analysis_paths):
    species = species_path.split('/')[-1]
    species = '_'.join(species.split('_')[1:-1]) ## get core species string of path
    tree_oi_file = sorted(glob.glob(f'{species_path}/*/*_ML*.raxml.bestTreeCollapsed')) ## path need to have ML tree
    if len(tree_oi_file) == 0:
        print(f'{species} has no ML tree in path. Skipping...\n')
        continue
    else:
        tree_oi_file = tree_oi_file[-1]
        most_recent_path = '/'.join(tree_oi_file.split('/')[:-1])
    
    annotation_SNV_file = glob.glob(f'{most_recent_path}/*_annotation_SNVs.csv') ## wildcard needed for timestamp
    if len(annotation_SNV_file) == 0:
        print(f'No SNV annotation file found for species {species} in path {most_recent_path}')
        continue
    else:
        annotation_SNV_file = annotation_SNV_file[0]
    snp_table_file = f'{most_recent_path}/snp_table.csv'
    indel_table_file = f'{most_recent_path}/indel_table.csv'

    ## output folder
    os.system(f'mkdir -p {res_folder}/{species}')
    os.chdir(f'{res_folder}/{species}')

    ## write known ancestry table --> get that from tree file!
    tree = load_tree(tree_oi_file) ## load tree
    
    ## extract the frequencies independently from tree
    isolate_labels = [tip.name for tip in tree.get_terminals() if not ( tip.name.endswith('_ref' ) | (tip.name.startswith('OG_') & tip.name.endswith('_TP0')) )]
    annotation_mutations_dedup = pd.read_csv(annotation_SNV_file)
    timepoints_tree = np.unique([label.split('_')[1] for label in isolate_labels])
    num_timepoints_winf = len(timepoints_tree)

    snp_table = pd.read_csv(snp_table_file)
    snp_table_metadata_cols = ['chr', 'pos', 'type', 'muts', 'locustag', 'gene', 'loc1', 'loc2', 'strand', 'product', 'nt_pos', 'aa_pos', 'nt_ref', 'nt_anc', 'nt_alt']
    snp_table_isolate_cols = [col for col in snp_table.columns if col not in snp_table_metadata_cols]
    
    snp_table['muts'] = snp_table['muts'].str.replace("\\[|\\]|'|", '', regex = True).str.replace('.', '')
    snp_table['gene'] = snp_table['gene'].replace('citG,mdcB', 'citG')

    ## explode df to get for evey individual mutation a spearate line 
    for col in ['type', 'nt_alt']:
        snp_table[col] = snp_table[col].apply(lambda x: [nt for nt in x])
    snp_table['muts'] = snp_table['muts'].fillna('').apply(lambda x: [mut if mut == mut else '' for mut in x.split(', ')])
    snp_table = snp_table.explode(['type', 'nt_alt', 'muts'])

    ## correct polarization of mutation in Ehormaechei; note need to also change the polarization for the other allele on the same locus!!! (required for snp table without repolarization on the ancestral sequence (polarized based on outgroup!))
    if species == 'P07_Ehormaechei-c1':
        is_gyrA_wrong_pol = (snp_table['chr'] == 7) & (snp_table['pos'] == 108226)
        snp_table.loc[is_gyrA_wrong_pol, 'muts'] = snp_table.loc[is_gyrA_wrong_pol, 'muts'].str.replace('F83S', 'S83F').str.replace('F83Y', 'S83Y')
        snp_table.loc[(snp_table['nt_anc'] == snp_table['nt_alt']), 'nt_alt'] = snp_table.loc[snp_table['nt_anc'] == snp_table['nt_alt'], 'nt_ref']

    snp_table = rename_locustags(snp_table, locustag_str = 'MDIMIJ')

    ## read in indel table
    if os.path.exists(indel_table_file):
        indel_table = pd.read_csv(indel_table_file)
        indel_table_metadata_cols = ['gene_num', 'gene_num_global', 'chr', 'pos', 'chr_boundary', 'indel_ref', 'indel_anc', 'indel_alt', 'max_indel_GL_diff_to_ref', 'indel_size_gp', 'type', 'product', 'gene', 'protein_id', 'strand', 'loc1', 'loc2', 'sequence', 'note', 'locustag', 'orthologtag', 'translation', 'indel_pos_start', 'indel_pos_end', 'alt_start', 'alt_stop', 'alt_stop_pos', 'frameshift', 'mut_translation', 'mut']
        for col in ['alt_stop', 'indel_alt', 'indel_pos_start', 'indel_pos_end']:
            indel_table[col] = indel_table[col].str.replace("\\[|\\]|'|,|\\(|\\)", '', regex = True).str.replace('np.int64', '').str.replace('.', '')
            if col in ['indel_pos_start', 'indel_pos_end']:
                indel_table[col] = indel_table[col].fillna(-1).astype(str)
            else: 
                indel_table[col] = indel_table[col].fillna('')
            indel_table[col] = indel_table[col].apply(lambda x: [allele for allele in x.split(' ')])
        indel_table = indel_table.explode(['alt_stop', 'indel_alt', 'indel_pos_start', 'indel_pos_end'])
        indel_table['indel_pos_aa_start'] = (indel_table[['indel_pos_start', 'indel_pos_end']].astype(int).min(axis = 1) / 3).astype(int)
        indel_table.loc[indel_table['indel_pos_aa_start'] == 0, 'indel_pos_aa_start'] = ''
        indel_table['muts'] = indel_table['indel_pos_aa_start'].astype(str) + indel_table['indel_alt']
        indel_table_metadata_cols += ['muts', 'indel_pos_aa_start']

        indel_table = rename_locustags(indel_table, locustag_str = 'MDIMIJ')

        ## rename samples in indel_table
        smpl_names_dict = {smpl.split('__')[-1]: smpl for smpl in snp_table_isolate_cols}
        indel_table.rename(smpl_names_dict, axis = 1, inplace = True)


    snp_freq_shifts_anno_lod = [] ## delete that line after troubleshooting!
    af_change_no_singletons_list = {} ## list to append to to plot after a histogram showing the distribution of af changes 
    af_change_singletons_list = {} ## list to append to to plot after a histogram showing the distribution of af changes 
    count = 0
    unique_locs = locations_abbreviations# sorted(np.unique([smpl.split('_')[4] for smpl in snp_table_isolate_cols]))
    
    subj_timepoints = np.array([label.split('_')[1] for label in snp_table_isolate_cols])
    unique_tps = sorted(np.unique(subj_timepoints))
    inf_date = [label.split('_')[2] for label in snp_table_isolate_cols if label.split('_')[3] == 'D']
    inf_tp = np.unique([label.split('_')[1].replace('TP', '') for label in snp_table_isolate_cols if label.split('_')[3] == 'D'])
    if len(inf_date) == 0:
        print(f'{species} is no pathogenic lineage. Using the first observance of that lineage')
        first_inf_date = min([label.split('_')[2] for label in snp_table_isolate_cols])
    else:
        first_inf_date = min(inf_date)
    first_inf_date = datetime.strptime(first_inf_date, "%Y%m%d")
    timepoints_to_date_dict = {label.split('_')[1]: datetime.strptime(label.split('_')[2], "%Y%m%d") for label in snp_table_isolate_cols}
    timepoints_to_inf_days_dict = {tp.replace('TP', ''): (date - first_inf_date).total_seconds() / 86400 for tp, date in timepoints_to_date_dict.items()}
    ## add small error to subsequent timepoints with same date --> Need to read in later the exact timepoints!
    timepoints_to_inf_days_l = sorted(timepoints_to_inf_days_dict.items(), key=lambda x: x[0])
    timepoints_to_inf_days_dict_adj = {}
    prev_val = None
    offset = 0
    margin = 0.1
    for key, val in timepoints_to_inf_days_l:
        if float(val) == prev_val:
            offset += margin
            timepoints_to_inf_days_dict_adj[key] = float(val) + offset
        else:
            offset = 0
            prev_val = float(val)
            timepoints_to_inf_days_dict_adj[key] = float(val)
    inf_days = [day for tp, day in timepoints_to_inf_days_dict_adj.items() if tp in inf_tp]
    subj_location = np.array([label.split('_')[4] for label in snp_table_isolate_cols])

    ## get lookup table 
    smpl_at_tp_dict = {tp: [] for tp in unique_tps}
    smpl_at_loc_dict = {loc: [] for loc in unique_locs}
    for smpl in snp_table_isolate_cols: 
        smpl_split = smpl.split('_')
        smpl_at_tp_dict[smpl_split[1]].append(smpl)
        if smpl_split[4] in smpl_at_loc_dict.keys():
            smpl_at_loc_dict[smpl_split[4]].append(smpl)
        else:
            print(f'Samples in {smpl_split[4]} will not be recorded for per location evaluation')

    if os.path.exists(indel_table_file):
        mutations_to_evaluate = {'SNVs': [snp_table, snp_table_metadata_cols],
                                'Indel': [indel_table, indel_table_metadata_cols]}
    else:
        mutations_to_evaluate = {'SNVs': [snp_table, snp_table_metadata_cols]}

    for muttype, (mut_df, metadata_cols) in mutations_to_evaluate.items():
        af_change_no_singletons_list[muttype] = [] ## list to append to to plot after a histogram showing the distribution of af changes 
        af_change_singletons_list[muttype] = [] ## list to append to to plot after a histogram showing the distribution of af changes 
        for idx, mutation_info_row in mut_df.iterrows():
            freq_timepoint_dict = {tp: {loc: {} for loc in np.concatenate([['collapsed'], unique_locs])} for tp in unique_tps}
            count_spl_timepoint = {tp: {loc: 0 for loc in np.concatenate([['collapsed'], unique_locs])} for tp in unique_tps}

            unique_alleles = mutation_info_row[[col for col in metadata_cols if col.endswith('_alt')]].values

            for tp in unique_tps:
                smpls_at_tp = smpl_at_tp_dict[tp]
                alleles_at_tp = mutation_info_row[smpls_at_tp]
                for allele in unique_alleles:
                    freq_timepoint_dict[tp]['collapsed'][allele] = np.sum(alleles_at_tp == allele) / np.sum(alleles_at_tp != '?')
                count_spl_timepoint[tp]['collapsed'] = len(smpls_at_tp)
                for loc in locations_abbreviations:
                    smpls_at_tp_loc = [smpl for smpl in smpls_at_tp if smpl in smpl_at_loc_dict[loc]]
                    if len(smpls_at_tp_loc) > 0:
                        alleles_at_tp_loc = mutation_info_row[smpls_at_tp_loc]
                        for allele in unique_alleles:
                            freq_timepoint_dict[tp][loc][allele] = np.sum(alleles_at_tp_loc == allele) / np.sum(alleles_at_tp_loc != '?')
                    count_spl_timepoint[tp][loc] = len(smpls_at_tp_loc)
            
            # calc max allele frequency chronological (min >> max from vi to vn; subjects with only one timepoint not considered)
            valid_tps = [tp for tp in count_spl_timepoint.keys() if count_spl_timepoint[tp]['collapsed'] >= 2]
            
            if len(valid_tps) > 1: 
                max_af_change = 0
                for allele in unique_alleles: 
                    allele_freqs = [loc_dict['collapsed'][allele] for tp, loc_dict in freq_timepoint_dict.items() if tp in valid_tps]
                    af_change_allele = np.max(allele_freqs) - np.min(allele_freqs)
                    max_af_change = np.max([af_change_allele, max_af_change])
                
                ## exclude singletons from dataset 
                num_smpls_with_support = np.unique([nt for nt in mutation_info_row[snp_table_isolate_cols] if (nt != '?')], return_counts = True)
                num_smpls_with_mut_minor_alleles = np.sum([cnt for cnt in num_smpls_with_support[1] if (cnt != np.max(num_smpls_with_support[1]))])
                if num_smpls_with_mut_minor_alleles == 1: ## accounts for homoplasies --> if 2 different alleles have each singletons, the site is kept! 
                    if max_af_change >= allele_freq_change_min: ## just print message if samples would exceed allele_freq_change_min
                        print(f"Excluding singleton in gene {mutation_info_row['locustag']} at chr {mutation_info_row['chr']} and pos {mutation_info_row['pos']} with max allele frequency change of {max_af_change} ")
                    af_change_singletons_list[muttype].append(max_af_change)
                    continue
                else:
                    af_change_no_singletons_list[muttype].append(max_af_change)
                
                if max_af_change >= allele_freq_change_min:
                    print(f"Found gene {mutation_info_row['locustag']} at chr {mutation_info_row['chr']} and pos {mutation_info_row['pos']} with max allele frequency change of {max_af_change}")
                    cand_mut_dict = {}
                    # build dict for pd.dataframe
                    cand_mut_dict['_chr'] = mutation_info_row['chr']
                    cand_mut_dict['_pos'] = mutation_info_row['pos']
                    cand_mut_dict['_tag'] = mutation_info_row['locustag']
                    
                    for tp, loc_dict in freq_timepoint_dict.items():
                        for loc, allele_dict in loc_dict.items():
                            for allele, freq in allele_dict.items():
                                cand_mut_dict[f'_timepoint-{tp}_location-{loc}_{allele}_freqs'] = freq
                            cand_mut_dict[f'_timepoint-{tp}_location-{loc}_counts'] = count_spl_timepoint[tp][loc]
                    for col in metadata_cols: # get all data from annotation_mutations_dedup dataframe
                        key_name = re.sub(r'^nt_|^indel_', '', col)
                        cand_mut_dict[key_name] = mutation_info_row[col]
                    snp_freq_shifts_anno_lod.append(cand_mut_dict)

    if len(snp_freq_shifts_anno_lod) == 0:
        print(f"\nNo loci exceeded min AF change of {allele_freq_change_min}. Skipping {species}")
        continue
    else:
        print(f"\n{len(snp_freq_shifts_anno_lod)} locis exceeded min AF change of {allele_freq_change_min} and where no singletons.")

    snp_freq_shifts_anno = pd.DataFrame(snp_freq_shifts_anno_lod) 

    ## make df to long format and split by site 
    snp_freq_shifts_anno_long = snp_freq_shifts_anno.melt([col for col in snp_freq_shifts_anno if not col.startswith('_timepoint')]) #['_chr', '_pos', '_tag', '_nt', '_species'] + snp_table_metadata_cols)
    snp_freq_shifts_anno_long['_timepoint'] = snp_freq_shifts_anno_long['variable'].str.split('_location-').str[0].str.split('-').str[-1]
    snp_freq_shifts_anno_long['_location'] = snp_freq_shifts_anno_long['variable'].str.split('_location-').str[-1].str.split('_').str[0]
    snp_freq_shifts_anno_long['_freq_cnt'] = snp_freq_shifts_anno_long['variable'].str.split('_').str[-1]
    snp_freq_shifts_anno_long['_timepoint_freq_cnt'] = snp_freq_shifts_anno_long['_timepoint'] + '_' + snp_freq_shifts_anno_long['_freq_cnt']
    snp_freq_shifts_anno_long = snp_freq_shifts_anno_long.drop(['_timepoint', '_freq_cnt', 'variable'], axis = 1)
    index_cols = [col for col in snp_freq_shifts_anno_long.columns if col not in ['value', '_timepoint_freq_cnt']]
    snp_freq_shifts_anno_long[index_cols] = snp_freq_shifts_anno_long[index_cols].fillna('')
    ## remove all entries with nan frequencies and counts --> no counts/frequencies as this allele is not present for the given position
    snp_freq_shifts_anno_long = snp_freq_shifts_anno_long.dropna(subset=['value'])
    
    snp_freq_shifts_anno_longwide = snp_freq_shifts_anno_long.pivot(values = 'value', index = index_cols, columns = '_timepoint_freq_cnt').reset_index()
    snp_freq_shifts_anno_longwide['genetag'] = np.where(snp_freq_shifts_anno_longwide['gene'] != '.', 
                                                        snp_freq_shifts_anno_longwide['gene'],
                                                        snp_freq_shifts_anno_longwide['_tag'])
    snp_freq_shifts_anno_longwide['Mutation'] = snp_freq_shifts_anno_longwide[['_chr', '_pos', 'type', 'genetag', 'muts', 'alt']].astype(str).agg('_'.join, axis = 1)
    snp_freq_shifts_anno_longwide.columns.name = None

    freq_cols = [col for col in snp_freq_shifts_anno_longwide.columns if col.endswith('_freqs')]
    cnt_cols = [col for col in snp_freq_shifts_anno_longwide.columns if col.endswith('_counts')]
    snp_freq_shifts_anno_longwide[freq_cols] = snp_freq_shifts_anno_longwide[freq_cols].astype(float)
    snp_freq_shifts_anno_longwide[cnt_cols] = snp_freq_shifts_anno_longwide[cnt_cols].astype(int)

    ## drop all empty locations with no isolates 
    snp_freq_shifts_anno_longwide = snp_freq_shifts_anno_longwide[snp_freq_shifts_anno_longwide[cnt_cols].sum(axis = 1) > 0]
    
    cols_of_interest = {'_location' : 'Location',
                        '_chr':'Chromosome', 
                        '_pos': 'Position', 
                        'alt': 'Alt_allele', 
                        'type': 'Class', 
                        'gene': 'Gene',
                        'Mutation': 'Mutation'}
    
    muller_tbl_freqs = snp_freq_shifts_anno_longwide[list(cols_of_interest.keys()) + freq_cols] 
    muller_tbl_freqs = muller_tbl_freqs.rename(cols_of_interest, axis = 1)
    muller_tbl_freqs = muller_tbl_freqs.rename({col: col.replace('_freqs', '').replace('TP', '') for col in freq_cols}, axis = 1)
    muller_tbl_freqs = muller_tbl_freqs.rename(timepoints_to_inf_days_dict_adj, axis = 1)
    muller_tbl_cnts = snp_freq_shifts_anno_longwide[list(cols_of_interest.keys()) + cnt_cols]
    muller_tbl_cnts = muller_tbl_cnts.rename(cols_of_interest, axis = 1)
    muller_tbl_cnts = muller_tbl_cnts.rename({col: col.replace('_counts', '').replace('TP', '') for col in cnt_cols}, axis = 1)
    muller_tbl_cnts = muller_tbl_cnts.rename(timepoints_to_inf_days_dict_adj, axis = 1)

    if add_missing_timepoints:
        if species in missing_timepoints.keys():
            for tp in missing_timepoints[species]:
                muller_tbl_freqs[tp] = 0
                muller_tbl_cnts[tp] = 0

    ## shorten the mutation names
    muller_tbl_freqs['Mutation'] = muller_tbl_freqs['Mutation'].str.split('_').str[3:5].str.join('[') + ']'

    ## plot number of snvs and indels per frequency change
    bins = np.arange(0, 1.02, 0.02)
    num_cols = len(mutations_to_evaluate.keys())
    fig, axs = plt.subplots(figsize=(3*num_cols+2, 4), ncols = num_cols, sharey = True)
    axs = axs.flatten()
    for ax, muttype in zip(axs, mutations_to_evaluate.keys()):
        af_change_hist_df = pd.DataFrame()
        af_change_hist_df['af_change'] = af_change_no_singletons_list[muttype] + af_change_singletons_list[muttype]
        af_change_hist_df['singletons'] = ['clonal'] * len(af_change_no_singletons_list[muttype]) + ['singleton'] * len(af_change_singletons_list[muttype])
        sns.histplot(data=af_change_hist_df, x="af_change", hue="singletons", edgecolor = 'k', lw = 0.4,
                    multiple="stack", 
                    bins = bins,
                    alpha = 0.8,
                    ax = ax, 
                    zorder = 5)
        ax.axvline(allele_freq_change_min, ls = '--', color = 'grey', zorder = 7)
        ax.text(allele_freq_change_min, ax.get_ylim()[1] * 0.95, 'min AF shift threshold', 
                color='k', ha='right', va='top', rotation=90, size = 10, zorder = 10)
        ax.grid(True, lw = 0.4, zorder = 0)
        ax.set_xticks(np.round(np.arange(0, 1.1, 0.1), 2))
        ax.set_xticklabels(np.round(np.arange(0, 1.1, 0.1), 2))
        ax.set_xlabel('Maximum allele frequency change')
        ax.set_title(muttype)
    plt.suptitle(f'Maximum allele frequency change\nper allele in {species}')
    plt.tight_layout()
    fig.savefig(f'{species}_allele_frequency_shift.pdf')
    
    ## get timepoints
    timepoint_cols_all = np.array(sorted(muller_tbl_freqs.columns[7:]))
    timepoints = timepoint_cols_all.astype(float)
    if len(timepoints) < 2:
        print(f'{species} has too few timepoints to interpolate frequencies')
        continue

    ## Loop over every location in the muller table generated above
    for loc in sorted(muller_tbl_freqs['Location'].unique()):
        muller_tbl = muller_tbl_freqs[muller_tbl_freqs['Location'] == loc]

        ## remove from non-collapsed plot the timepoint(s) of infection
        if (loc != 'collapsed'):
            is_mbiome_only = muller_tbl_cnts['Location'].isin(['Nasal', 'Rectal', 'Oral', 'Skin']) ## get only mbiome sites without inf isolates
            # num_isolates_on_tp = muller_tbl_cnts.loc[is_mbiome_only, timepoint_cols_all].sum(axis = 0)
            timepoint_cols_of_inf = np.array([col for col in timepoint_cols_all if col in inf_days])
            timepoint_cols = np.array([col for col in timepoint_cols_all if col not in inf_days])
            muller_tbl = muller_tbl.drop(inf_days, axis = 1)

            if drop_empty_tp: 
                cnts_on_tp = muller_tbl_cnts.loc[muller_tbl_freqs['Location'] == loc, timepoint_cols]
                tp_has_no_counts = cnts_on_tp.columns[(cnts_on_tp.sum(axis=0) == 0)]
                tp_has_no_counts = [tp for tp in tp_has_no_counts if tp != min(timepoint_cols)] ## leave first timepoint for interpolation in
                timepoint_cols = np.array([col for col in timepoint_cols if col not in tp_has_no_counts])
                muller_tbl = muller_tbl.drop(tp_has_no_counts, axis = 1)

            if len(timepoint_cols) < 2:
                print(f'{species} has too few timepoints at site {loc} to interpolate frequencies')
                continue
        else:
            timepoint_cols = timepoint_cols_all
            if drop_empty_tp: 
                cnts_on_tp = muller_tbl_cnts.loc[muller_tbl_freqs['Location'] == loc, timepoint_cols]
                tp_has_no_counts = cnts_on_tp.columns[(cnts_on_tp.sum(axis=0) == 0)]
                tp_has_no_counts = [tp for tp in tp_has_no_counts if tp != min(timepoint_cols)] ## leave first timepoint for interpolation in
                timepoint_cols = np.array([col for col in timepoint_cols if col not in tp_has_no_counts])
                muller_tbl = muller_tbl.drop(tp_has_no_counts, axis = 1)
        
        ## get only mutations with any frequency
        muller_tbl = muller_tbl[muller_tbl[timepoint_cols].sum(axis = 1) > 0]
        
        if muller_tbl.empty:
            ## no mutation at location identified
            continue
        
        gt_to_isolate_dict = {}
        for muttype, (mut_df, _) in mutations_to_evaluate.items():
            tmp_gt_to_isolate_dict = get_tips_w_mutation(mut_df, snp_table_isolate_cols, muller_tbl) ## get dict with isolates which have mutation
            gt_to_isolate_dict = {**gt_to_isolate_dict, **tmp_gt_to_isolate_dict}

        ## check if any mutation is linked and if so deduplicate
        muller_tbl, gt_to_isolate_dict_dedup = collapse_linked_mutations(muller_tbl, gt_to_isolate_dict, muller_table_columns_to_combine = muller_tbl.columns[1:7])

        ## dictionary with genotype (child) to parental genotype's 
        ## NOTE: Does not handle homoplasies! --> check before if that does not apply for the current dataset!
        gt_ancestral_dict, ancestral_to_gt_dict = get_closest_anc_per_genotype_of_interest(tree, gt_to_isolate_dict_dedup, root_name=root_name, finetune_polarization=True) ## get for each mutaton its closest ancestry amongst the given other mutations 

        if convert_alleles_to_genotypes:
            allele_to_gt = {root_name: root_name}
            for derived_allele, parent_allele in gt_ancestral_dict.items():
                if parent_allele == root_name:
                    allele_to_gt[derived_allele] = derived_allele
                else:
                    der = derived_allele
                    derived_gt = [derived_allele]
                    while gt_ancestral_dict[der] != root_name:
                        derived_gt.append(gt_ancestral_dict[der])
                        der = gt_ancestral_dict[der]
                    allele_to_gt[derived_allele] = '; '.join(derived_gt[::-1])

        ## Parameters for interpolation
        ancestral_gt_to_walk = [root_name]
        tmp_dfs = []
        interpolated_df = pd.DataFrame()
        
        ## add basal genotype
        basal_gt = muller_tbl.iloc[0].copy()
        
        basal_gt['Mutation'] = root_name
        basal_gt[['Chromosome', 'Position', 'Alt_allele', 'Class']] = None 
        basal_gt[timepoint_cols] = 1 ## Note: basal genotype is parent/ancestor of all genotypes in tree!!! 
        basal_gt = basal_gt.to_frame().T
        basal_gt = basal_gt[basal_gt.columns[1:]]
        muller_tbl_mod = pd.concat([basal_gt, muller_tbl]).reset_index(drop = True)
        muller_tbl_mod[timepoint_cols] = muller_tbl_mod[timepoint_cols].astype(float)
        muller_tbl_mod[min(timepoint_cols)] = muller_tbl_mod[min(timepoint_cols)].fillna(0) ## no frequency on timepoint but kept to allow inference from beginning

        ## pad last timepoint by one day (important as otherwise the last timepoint would be cut off)
        muller_tbl_mod_pad = muller_tbl_mod.copy()
        timepoint_cols_pad = list(timepoint_cols) + [str(int(timepoint_cols[-1]) +1)]
        muller_tbl_mod_pad[timepoint_cols_pad[-1]] = muller_tbl_mod_pad[timepoint_cols[-1]] ## pad df at the last timepoint
        timepoints = [float(tp) for tp in timepoint_cols_pad]

        while len(ancestral_gt_to_walk) > 0:
            ancestral_gt = ancestral_gt_to_walk[0] ## get first entry to access
            print(ancestral_gt)
            if ancestral_gt == root_name:
                ## infer root at first
                trajectories_der = [root_name]
                derived_freqs = np.array(muller_tbl_mod_pad.loc[muller_tbl_mod_pad['Mutation'].astype(str).isin(trajectories_der), timepoint_cols_pad].T)
                x_softmax, softmax_data = softmax_interpolation(timepoints, derived_freqs, num_interp, np_log_err = np_log_err)
                tmp_df = pd.DataFrame(data = softmax_data.T, columns = x_softmax, index = trajectories_der)
                interpolated_df = pd.concat([interpolated_df, tmp_df])
            if ancestral_gt in ancestral_to_gt_dict.keys():
                trajectories_der = [traj for traj in ancestral_to_gt_dict[ancestral_gt]]
                trajectort_anc = ancestral_gt
                ## sort mullertable based on trajectories_der
                muller_row_is_der = muller_tbl_mod_pad['Mutation'].astype(str).isin(trajectories_der)
                derived_freqs = np.array(muller_tbl_mod_pad[muller_row_is_der].set_index('Mutation').loc[trajectories_der, timepoint_cols_pad].T)
                ancestral_freqs = np.array(interpolated_df.loc[trajectort_anc]).reshape(-1, 1)
                x_softmax, softmax_data = softmax_interpolation(timepoints, derived_freqs, num_interp, groups_overall_frequencies=ancestral_freqs, np_log_err = np_log_err)
                interpolated_df.loc[trajectort_anc] = softmax_data[:, 0] ## update ancestral frequencies
                tmp_df = pd.DataFrame(data = softmax_data[:, 1:].T, columns = x_softmax, index = trajectories_der)
                interpolated_df = pd.concat([interpolated_df, tmp_df])
                ## add the other genotypes to 
                if ancestral_gt in ancestral_to_gt_dict.keys():    
                    ancestral_gt_to_walk += ancestral_to_gt_dict[ancestral_gt]
            # remove current genotype --> it was accessed
            ancestral_gt_to_walk = ancestral_gt_to_walk[1:]
        
        ## validate the data 
        ## get timepoints most closely related
        interpolated_idx_of_tps = [np.argmin(np.abs(interpolated_df.columns - int(col))) for col in timepoint_cols]
        interpolated_df_of_tps = interpolated_df[interpolated_df.columns[interpolated_idx_of_tps]].reset_index()
        interpolated_df_of_tps.columns = ['Mutation'] + list(timepoint_cols)

        obs_data_frequencies = muller_tbl_mod[['Mutation'] + list(timepoint_cols)].melt('Mutation')
        obs_data_frequencies['Mutation'] = obs_data_frequencies['Mutation'].astype(str)
        est_data_frequencies = interpolated_df_of_tps.melt('Mutation')
        est_data_frequencies['Mutation'] = est_data_frequencies['Mutation'].astype(str)

        comparison_df = pd.merge(obs_data_frequencies,
                                est_data_frequencies, 
                                how = 'outer', on = ['Mutation', 'variable'], suffixes = ['_obs', '_est'])
        comparison_df['diff'] = comparison_df['value_obs'] - comparison_df['value_est']
        comparison_df_major_diff = comparison_df[comparison_df['diff'] >= 1e-3] ## note the major genotypes of which some derived genotypes appeared have higher differences as the derived frequency needed to be removed for the ggmuller!

        ## save the inferred frequencies
        edges_df = pd.DataFrame([(anc, der) for der, anc in gt_ancestral_dict.items()], columns = ['Parent', 'Identity'])
        population_df = interpolated_df.reset_index().melt('index')
        population_df.columns = ['Identity', 'Generation', 'Population']

        ## get counts table 
        cnt_table = muller_tbl_cnts.loc[muller_tbl_freqs['Location'] == loc, timepoint_cols]
        cnt_table = cnt_table.reset_index(drop = True).iloc[0] ## counts are the same for every mutation --> pic first entry
        cnt_table = pd.DataFrame(cnt_table).reset_index()
        cnt_table.columns = ['Days', 'cnts']
        cnt_table['Species'] = species
        cnt_table['Location'] = loc
        cnt_table = cnt_table[['Species', 'Location', 'Days', 'cnts']]
        
        if convert_alleles_to_genotypes:
            edges_df['Parent'] = edges_df['Parent'].map(allele_to_gt)
            edges_df['Identity'] = edges_df['Identity'].map(allele_to_gt)
            population_df['Identity'] = population_df['Identity'].map(allele_to_gt)

        ## check if every identity in edges_df is present in population_df; else combine!
        unique_pop_identities = population_df['Identity'].unique()


        edges_df.to_csv(f'{species}_edges_{loc}_freq{allele_freq_change_min}_windels.tsv', index = False, sep = '\t')
        population_df.to_csv(f'{species}_pop_{loc}_freq{allele_freq_change_min}_windels.tsv', index = False, sep = '\t')
        cnt_table.to_csv(f'{species}_cnt_{loc}_freq{allele_freq_change_min}_windels.tsv', index = False, sep = '\t')
