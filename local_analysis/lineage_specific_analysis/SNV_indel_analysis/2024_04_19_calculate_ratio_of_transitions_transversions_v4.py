
## identify per isolate the number of transitions and traversions and generate a per node plot 

## TODO:
## Need to perform simulation to test enrichment of pattern 
## (given the underlying sequence and test for enrichment)
## plot pval for enrichment into square (or asterisks into normalized enrichment table)

## NOTE:
## it is always comparing to ancestral state but one should compare to node of interest (i would assume to allow accounting for reversions)

## NOTE: 
## for some sequences the ancestral sequence and the anc sequence from the reconstruction is not the same as the outgroup might have a basal branch for itself!

import numpy as np
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import seaborn as sns
from Bio import SeqIO
from Bio import Phylo
import networkx
from scipy.stats import binom, chi2_contingency

plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

sense_antisense_redundancy_dict = {
    'A>T': 'T>A\nA>T',
    'A>C': 'T>G\nA>C',
    'A>G': 'T>C\nA>G',
    'T>A': 'T>A\nA>T',
    'T>C': 'T>C\nA>G',
    'T>G': 'T>G\nA>C',
    'C>A': 'C>A\nG>T',
    'C>T': 'C>T\nG>A',
    'C>G': 'C>G\nG>C',
    'G>A': 'C>T\nG>A',
    'G>T': 'C>A\nG>T',
    'G>C': 'C>G\nG>C',
}

trans_trav_dict = {'T>A\nA>T': 'Transversions',
                   'T>C\nA>G': 'Transitions',
                   'T>G\nA>C': 'Transversions',
                   'C>A\nG>T': 'Transversions',
                   'C>T\nG>A': 'Transitions',
                   'C>G\nG>C': 'Transversions'}

nondirectional_mutsig = {'T>A\nA>T': 'T>A:A>T',
                         'T>C\nA>G': 'T>C:A>G\nC>T:G>A',
                         'T>G\nA>C': 'T>G:A>C\nC>A:G>T',
                         'C>A\nG>T': 'T>G:A>C\nC>A:G>T',
                         'C>T\nG>A': 'T>C:A>G\nC>T:G>A',
                         'C>G\nG>C': 'C>G:G>C'}

nts = ['A', 'T', 'C', 'G']
exclude_leaf_with_str=['_ASM_', '_GCF_', '_GCA_', '_ref'] ## outgroup/ref

trans_trav_df = pd.DataFrame().from_dict(trans_trav_dict, orient = 'index').reset_index()
trans_trav_df.columns = ['trans_trav', 'trans_trav_class']

candidate_hypermutators = {}
threshold_cand_hypermutator = {'min_mut_total_tree': 5}
outpath_candidate_hypermutator = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/candidate_hypermutators_based_on_outgroup_background.csv')

def parse_tree(tree_path):
    """helper function to parse tree file, with error catching. Returns first tree in a file, if multiple"""
    try:
        parsed_tree=Phylo.read(tree_path,'newick')
    except ValueError:
        parsed_tree=None
        #print("Multiple trees in file, trying to parse first (filled) tree in file.")
        trees= Phylo.parse(tree_path,'newick')
        for tree in trees:
            if tree.count_terminals() > 1:
                parsed_tree=tree
                #print("Found filled tree, with length", tree.count_terminals())
                break
    return parsed_tree

def read_in_fasta(ancestral_reconstruction_fasta):
    calls_per_node_leaf_dict={}
    for record in SeqIO.parse(ancestral_reconstruction_fasta, 'fasta'):
        calls_per_node_leaf_dict[record.id] = str(record.seq)
    return calls_per_node_leaf_dict

def create_tree_network(ancestral_reconstruction_labelled_tree):
    # parse tree structure into adjacency graph to traverse
    parsed_tree = parse_tree(ancestral_reconstruction_labelled_tree) 
    
    net = Phylo.to_networkx(parsed_tree)
    tree_as_dict_of_lists=networkx.to_dict_of_lists(net) ## NOTE: creates bidirectional graph from each node! 
    return parsed_tree, tree_as_dict_of_lists

def get_root_of_tree(parsed_tree, exclude_leaf_with_str = []):
    all_leaves = parsed_tree.get_terminals()
    leaves_of_interest = []
    for leaf in all_leaves:
        if any([exclude_str in leaf.name for exclude_str in exclude_leaf_with_str]):
            continue
        leaves_of_interest.append(leaf.name)
    root = parsed_tree.common_ancestor(leaves_of_interest)
    return root

def read_background_signature(background_mut_signature_file, combine_transtravs = False, trans_trav_dict = trans_trav_dict):
    mut_background = pd.read_csv(background_mut_signature_file)
    if combine_transtravs:
        mut_background['mut_type'] = mut_background['mut_sig'].map(trans_trav_dict)
        mut_background = mut_background.groupby(['mut_type'])[['cnts', 'freq']].agg('sum').reset_index()
        mut_background_dict = {muttype: cnt for muttype, cnt in mut_background[['mut_type', 'cnts']].values}
    else:
        mut_background_dict = {muttype: cnt for muttype, cnt in mut_background[['mut_sig', 'cnts']].values}
    return mut_background_dict

def count_mutations_on_each_node(parsed_tree, root, tree_as_dict_of_lists, calls_per_node_leaf_dict, exclude_leaf_with_str = [], skip_preterminals = False):
    # start tree traversal, checking if each internal node has same call as parent, if not iterate the val for that base +1
    transitions_transversion_dict={}
    to_visit=[(x, root) for x in tree_as_dict_of_lists[root]]

    while len(to_visit)>0:
        currently_processing=to_visit[0]
        parent=currently_processing[1]
        current_node=currently_processing[0]
        ## check if parent is really parent of node, if not remove entry! 
        if not parsed_tree.find_any(parent.name).is_parent_of(current_node.name):
            to_visit=to_visit[1:]
            continue
        ## check if descendent has subpattern --> remove 
        elif any([exclude_str in current_node.name for exclude_str in exclude_leaf_with_str]):
            to_visit=to_visit[1:]
            continue
        transitions_transversion_dict[current_node.name] = {mut_type: 0 for mut_type in trans_trav_dict.keys()}
        is_preterminal=len(tree_as_dict_of_lists[current_node])==1
        if skip_preterminals and is_preterminal:
            pass
        else:
            for ancestral_call, derived_call in zip(calls_per_node_leaf_dict[parent.name], calls_per_node_leaf_dict[current_node.name]):
                if ancestral_call != derived_call and ancestral_call != 'N' and derived_call != 'N': ## check if non-ambigous mutation
                    mut_type = sense_antisense_redundancy_dict[f'{ancestral_call}>{derived_call}'] ## translate mutation type
                    transitions_transversion_dict[current_node.name][mut_type] += 1
            for children in tree_as_dict_of_lists[current_node]:
                if children.name not in transitions_transversion_dict.keys():
                    to_visit.append((children,current_node))
        to_visit=to_visit[1:]
    return transitions_transversion_dict

def combine_mutation_type_counts_per_node(parsed_tree, transitions_transversion_dict, root):
    ## collapse counts for every node 
    transitions_transversion_per_node = {}
    for node in parsed_tree.get_nonterminals():
        ## check if node is descendent of root or root itself
        if not root.is_parent_of(node.name):
            continue
        transitions_transversion_per_node[node.name] = {mut_type: 0 for mut_type in trans_trav_dict.keys()}
        descendents_of_node = node.find_clades()
        for descendent in descendents_of_node:
            if descendent.name in transitions_transversion_dict.keys():
                for mut_type, cnt in transitions_transversion_dict[descendent.name].items():
                    transitions_transversion_per_node[node.name][mut_type] += cnt
    return transitions_transversion_per_node

def normalize_mutation_counts(transitions_transversion_per_node):
    norm_transitions_transversion_per_node = {}
    mutations_per_node = {}
    for node in transitions_transversion_per_node.keys():
        mutations_per_node[node] = sum(list(transitions_transversion_per_node[node].values()))
        if sum(list(transitions_transversion_per_node[node].values())) > 0:
            norm_transitions_transversion_per_node[node] = {mut_type: cnt/sum(list(transitions_transversion_per_node[node].values())) for mut_type, cnt in transitions_transversion_per_node[node].items()}
        else:
            print('No normalization possible')
            norm_transitions_transversion_per_node[node] = transitions_transversion_per_node[node]
    return norm_transitions_transversion_per_node, mutations_per_node

def generate_df_from_dict(dict, trans_trav_dict, value_type = 'cnt'):
    trans_tranv_df = pd.DataFrame(dict).T.reset_index()
    trans_tranv_df = trans_tranv_df.rename({'index': 'node_id'}, axis = 1)
    trans_tranv_df = trans_tranv_df.melt('node_id', var_name='mut_sig', value_name=value_type)
    trans_tranv_df = trans_tranv_df.sort_values('node_id')
    trans_tranv_df['mut_sig_class'] = trans_tranv_df['mut_sig'].map(trans_trav_dict)
    trans_tranv_df['mut_sig'] = pd.Categorical(trans_tranv_df['mut_sig'], trans_trav_dict.keys(), ordered = True) ## sort df by dict order 
    return trans_tranv_df

def plt_barplot_of_node_w_background(df, background_df, ylim = '', node_of_interest = '', x = 'mut_sig', y = 'freq', hue = 'mut_sig_class', palette = 'BrBG', title = '', annotation_bg = '', annotation_obs = '', filename = ''):
    if node_of_interest: 
        df = df[df['node_id'] == node_of_interest]
    if not ylim: 
        ylim = (0, max(df[y])*1.1)

    fig, axs = plt.subplots(figsize = (7, 3), ncols = 2, sharex = True, sharey = True)
    sns.barplot(data = background_df, x = x, y = y, hue = hue, palette = palette, dodge = False, edgecolor = 'k', lw = 0.5, ax = axs[0])
    axs[0].annotate(annotation_bg, xy = (0.98, 0.96), xycoords = 'axes fraction', ha = 'right', va = 'top')
    axs[0].set_ylabel('Frequency')
    axs[0].legend().remove()
    axs[0].set_title('Background')
    ## observed frequencies
    sns.barplot(data = df, x = x, y = y, hue = hue, palette = palette, dodge = False, edgecolor = 'k', lw = 0.5, ax = axs[1])
    axs[1].annotate(annotation_obs, xy = (0.98, 0.96), xycoords = 'axes fraction', ha = 'right', va = 'top')
    axs[1].set_ylabel(None)
    axs[1].legend(title = 'Substitution type', loc = 'center left', bbox_to_anchor = (1, 0.5))
    axs[1].set_title('Observed')
    for ax in axs:
        ax.set_ylim(ylim)
        ax.set_xlabel(None)
        ax.yaxis.grid(which = 'both', linewidth = 0.5, linestyle = ':', color = 'grey', zorder = 0)
        ax.set_axisbelow(True)
    fig.suptitle(title, y = 0.92)
    plt.tight_layout()
    if filename:
        fig.savefig(os.path.expanduser(filename))

## load snp tables
analysis_paths = glob.glob(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/backup/P0*'))
mut_signature_background_path = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/refbased_2024_08/analysis/outgroup_mut_signature/')
pathogen_lineage_list = ['P07_Apittii-c1', 
                         'P07_Bthetaiotaomicron-c1', 
                         'P07_Ecoli-ST131-c4', 
                         'P07_Ehormaechei-c1', 
                         'P07_Paeruginosa-c1', 
                         'P10_Kmichiganensis-c1', 
                         'P10_Smaltophilia-c1', 
                         'P21_Ecoli-ST23-c2', 
                         'P21_Paeruginosa-c1', 
                         'P21_Pmirabilis-c2', 
                         'P21_Saureus-ST45-c2', 
                         'P21_Saureus-ST97-c3']
plotting_dict = {} ## saves all info during loop for one final plotting graphic

## USES Tree and internal nodes as ancestral sequence and extracts from there all mutations
## This helps that only on each node the mutational signature will be investigated and this 
## allows then to infer if at the node some genomic changes appeared which altered mutational signature
##
## implemented background signature extracted from outgroups mapped to reference genome and 
## calculated pvalues based on transitions:transversions only as well as on a nondirectional 
## mutation signature (1 transition, 3 transversions) to increase sensitivitiv in low mutational counts

for analysis_path in sorted(analysis_paths): 
    genome = analysis_path.split('/')[-1]
    genome = '_'.join(genome.split('_')[1:-1])
    avail_analysis_path = glob.glob(f'{analysis_path}/*/*.tree')
    if avail_analysis_path == []:
        continue
    avail_analysis_path = ['/'.join(path.split('/')[:-1]) for path in avail_analysis_path]
    mostrecent_analysis_path = sorted(avail_analysis_path, key = lambda x: x.rsplit('/')[-1])[-1]
    
    ancestral_tree_inference_file = f'{mostrecent_analysis_path}/node_sequence_inference/annotated_tree.nexus'
    ancestral_fasta_inference_file = f'{mostrecent_analysis_path}/node_sequence_inference/ancestral_sequences.fasta'
    outpath = f'{mostrecent_analysis_path}/pdf/mutational_signature/'
    ## select correct mutational background path
    patient = genome.split('_')[0]
    patient_long = patient[0] + patient[1:].zfill(4)
    species = genome.split('_')[1].split('-')[0]
    background_mut_signature_file = glob.glob(f'{mut_signature_background_path}{patient_long}_{species}*/*/mutational_spectrum_outgroups_to_consensus_call_core_genome.csv')
    mostrecent_bg_mut_sig_path = sorted(background_mut_signature_file, key = lambda x: x.rsplit('/')[-2])[-1]

    if not os.path.exists(ancestral_tree_inference_file):
        print(f'No ancestral reconstruction found for {genome}')
        continue 

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    calls_per_node_leaf_dict = read_in_fasta(ancestral_fasta_inference_file)
    parsed_tree, tree_as_dict_of_lists = create_tree_network(ancestral_tree_inference_file)
    root = get_root_of_tree(parsed_tree, exclude_leaf_with_str)
    mut_sig_background = read_background_signature(mostrecent_bg_mut_sig_path, True, {key: key for key in trans_trav_dict.keys()})
    mut_sig_background_transtravs_dict = read_background_signature(mostrecent_bg_mut_sig_path, True, trans_trav_dict)
    mut_sig_background_nondirmutsig_dict = read_background_signature(mostrecent_bg_mut_sig_path, True, nondirectional_mutsig)
    transitions_transversion_dict = count_mutations_on_each_node(parsed_tree, root, tree_as_dict_of_lists, calls_per_node_leaf_dict, exclude_leaf_with_str)
    transitions_transversion_per_node = combine_mutation_type_counts_per_node(parsed_tree, transitions_transversion_dict, root)
    norm_transitions_transversion_per_node, mut_per_node = normalize_mutation_counts(transitions_transversion_per_node)

    ## generate df
    norm_trans_tranv_per_node_df = generate_df_from_dict(norm_transitions_transversion_per_node, trans_trav_dict, value_type='freq')
    
    ## generate df from background dict
    bg_df = pd.DataFrame.from_dict(mut_sig_background, orient = 'index').reset_index()
    bg_df.columns = ['mut_sig', 'cnts']
    bg_df['mut_sig'] = pd.Categorical(bg_df['mut_sig'], trans_trav_dict.keys(), ordered = True)
    bg_df['freq'] = bg_df['cnts'] / sum(bg_df['cnts'])
    bg_df['mut_sig_class'] = bg_df['mut_sig'].map(trans_trav_dict)


    ## stats if distribution is skewed from background distribution
    ## combine transitions/traversions
    observed_counts_transtravs = {muttype: 0 for muttype in np.unique(list(trans_trav_dict.values()))}
    total_num_muts = sum(list(transitions_transversion_per_node[root.name].values()))
    for mutsig, cnt in transitions_transversion_per_node[root.name].items():
        observed_counts_transtravs[trans_trav_dict[mutsig]] += cnt
    observed_counts = []
    background_counts = []
    for muttype in observed_counts_transtravs.keys():
        observed_counts.append(observed_counts_transtravs[muttype])
        background_counts.append(mut_sig_background_transtravs_dict[muttype])
    observed_counts = np.array(observed_counts)
    background_counts = np.array(background_counts)
    if (sum(observed_counts) > 0) & (sum(background_counts) > 0):
        chi2_stat_transtravs, p_value_transtravs, _, _ = chi2_contingency([observed_counts, background_counts])
    else:
        p_value_transtravs = np.nan
    print(f'\n\nChi-squared result: {genome}:\t{p_value_transtravs}\t{True if (p_value_transtravs < 0.05) else False}')

    ## stats if distribution is skewed from background distribution
    ## combine nondirectional mutationtypes
    observed_counts_nondirmutsig = {muttype: 0 for muttype in np.unique(list(nondirectional_mutsig.values()))}
    total_num_muts = sum(list(transitions_transversion_per_node[root.name].values()))
    for mutsig, cnt in transitions_transversion_per_node[root.name].items():
        observed_counts_nondirmutsig[nondirectional_mutsig[mutsig]] += cnt
    observed_counts = []
    background_counts = []
    for muttype in observed_counts_nondirmutsig.keys():
        observed_counts.append(observed_counts_nondirmutsig[muttype])
        background_counts.append(mut_sig_background_nondirmutsig_dict[muttype])
    observed_counts = np.array(observed_counts)
    background_counts = np.array(background_counts)
    if (sum(observed_counts) > 0) & (sum(background_counts) > 0):
        chi2_stat_nondirmutsig, p_value_nondirmutsig, _, _ = chi2_contingency([observed_counts, background_counts])
    else:
        p_value_transtravs = np.nan
    print(f'\n\nChi-squared result: {genome}:\t{p_value_nondirmutsig}\t{True if (p_value_nondirmutsig < 0.05) else False}')

    print(genome)
    if threshold_cand_hypermutator['min_mut_total_tree'] <= mut_per_node[root.name]:
        is_cand_hypermutator = p_value_transtravs <= 0.05
        num_mutations = sum(transitions_transversion_per_node[root.name].values())
        candidate_hypermutators[genome] = [p_value_transtravs, p_value_nondirmutsig, is_cand_hypermutator, list(transitions_transversion_per_node[root.name].values())]
    ## Save all info to dict for final plotting
    if genome in pathogen_lineage_list:
        plotting_dict[genome] = {'mut_sig_df': norm_trans_tranv_per_node_df, 
                                    'node_of_interest': root.name,
                                    'num_muts': mut_per_node[root.name],
                                    'pval': p_value_transtravs,
                                    'background_mut_sig': mut_sig_background}
    
    plt_barplot_of_node_w_background(norm_trans_tranv_per_node_df, 
                                        bg_df,
                                        ylim = (0, 1.1),
                                        node_of_interest = root.name, 
                                        y = 'freq',
                                        palette = {'Transitions': 'grey', 'Transversions': 'lightgrey'},
                                        annotation_bg = f'n(mut)={sum(bg_df["cnts"])}', 
                                        annotation_obs = f'n(mut)={mut_per_node[root.name]}\np-val={round(p_value_transtravs, 3)}', 
                                        title = f'Normalized number on root\n on {genome}')#, 
                                        filename = f'{outpath}/2025_08_25_{genome}_mut_signature_wbackground_{root.name}.pdf')

with open(outpath_candidate_hypermutator, 'w') as fo:
    mutational_signatures = ','.join(trans_trav_dict.keys()).replace('\n', '_')
    fo.write(f'genome,pval_chi2_transitiontransversion,pval_chi2_nondirectionalmutationsignature,is_candidate_hypermutator,total_num_mutations,{mutational_signatures}\n')
    for genome, values in candidate_hypermutators.items():
        num_mutations = sum(values[3])
        mutational_counts_str = ','.join([str(cnt) for cnt in values[3]])
        fo.write(f'{genome},{round(values[0], 5)},{round(values[1], 5)},{values[2]},{num_mutations},{mutational_counts_str}\n')
