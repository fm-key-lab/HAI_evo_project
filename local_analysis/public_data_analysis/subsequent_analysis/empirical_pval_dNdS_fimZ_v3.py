
## calculate per gene dnds score
## extract empirical p value for dn ds value of fimbrial genes of public data

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import seaborn as sns
import glob 
import numpy as np

plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'


import sys
SCRIPTS_DIRECTORY = os.getcwd() + '/../../modules/'
sys.path.insert(0, SCRIPTS_DIRECTORY)

import analysispy_module as apy


### VARIABLES
annogenes_file = glob.glob(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2023/public_data/2023_12_analysis/analysis/2025_08_25_P07_Ehormaechei-c1_240825_gubbinsdefault_ambig10/*annotation_genes.csv'))
annogenes_file = sorted(annogenes_file)[-1]
anno_genes = pd.read_csv(annogenes_file)
anno_genes.loc[~anno_genes['codon_start'].isna() & anno_genes['codon_start'].str.isnumeric(), 'codon_start'] = anno_genes.loc[~anno_genes['codon_start'].isna() & anno_genes['codon_start'].str.isnumeric(), 'codon_start'].astype(int)
idx_to_extract = []
anno_genes_list = []
for idx, oldidx in anno_genes['Unnamed: 0'].items():
    if oldidx == 0:
        if idx_to_extract != []:
            anno_genes_list.append(anno_genes.iloc[idx_to_extract].dropna(axis = 1, how = 'all'))
        idx_to_extract = [idx]
    else:
        idx_to_extract.append(idx)
    if idx == max(anno_genes.index):
        anno_genes_list.append(anno_genes.iloc[idx_to_extract].dropna(axis = 1, how = 'all'))
annomut_file = glob.glob(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2023/public_data/2023_12_analysis/analysis/2025_08_25_P07_Ehormaechei-c1_240825_gubbinsdefault_ambig10/*annotation_mutations.csv'))
annomut_file = sorted(annomut_file)[-1]
anno_mutations = pd.read_csv(annomut_file)

ref_genome_folder = '~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/metadata/reference_genomes/'
refgenome = 'P07_Ehormaechei-c1_240825'

parameters = {
    'substitution_spectrum': os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2023/public_data/2023_12_analysis/analysis/2025_08_25_P07_Ehormaechei-c1_240825_gubbinsdefault_ambig10/mutationalspectrum.py.pk1')
    }

## calculate per gene the probability of nonsynonymous mutations
prob_N_per_gene_dict = apy.cacluate_expected_dn_ds_per_gene(parameters, anno_genes_list)

## get locust tag translations
snp_table = pd.read_csv(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2023/public_data/2023_12_analysis/analysis/2025_08_25_P07_Ehormaechei-c1_240825_gubbinsdefault_ambig10/snp_table.csv'))
## mark genes by the frequency of which they are covered at a snp 
non_smpl_cols = ['chr', 'pos', 'type', 'muts', 'locustag', 'gene', 'loc1', 'loc2', 'strand', 'product', 'nt_pos', 'aa_pos', 'nt_ref', 'nt_anc', 'nt_alt']
smpl_cols = [col for col in snp_table.columns if not col in non_smpl_cols]
snp_table['num_Ns'] = (snp_table[smpl_cols] == '?').sum(axis = 1)
snp_table['across_smpl_support'] = 1 - (snp_table['num_Ns'] / len(smpl_cols))
snp_table['Locustag'] = snp_table['chr'].astype(str) + '_' + snp_table['locustag']

gene_spec_N_cnt_dict = {}
gene_spec_S_cnt_dict = {}
for Locustag in snp_table['Locustag'].unique():
    gene_spec_N_cnt_dict[Locustag] = 0;
    gene_spec_S_cnt_dict[Locustag] = 0;
    cand_mut_anno = snp_table.loc[ snp_table['Locustag'] == Locustag ]
    for i,row in cand_mut_anno.iterrows():
        ## ensures only count if in genic regions
        gene_spec_N_cnt_dict[Locustag] += sum([muttype=='N' for muttype in row['type']])
        gene_spec_S_cnt_dict[Locustag] += sum([muttype=='S' for muttype in row['type']])

snp_table = snp_table[['locustag', 'Locustag', 'gene', 'product', 'across_smpl_support']] # .drop_duplicates() ## drop columns and duplictaes
snp_table = snp_table[~snp_table['Locustag'].str.contains(';')].reset_index(drop = True) ## drop columns of snps in intergenic regionms (two concatenated locus tags)
snp_table = snp_table.groupby(['locustag', 'Locustag', 'gene', 'product'])['across_smpl_support'].mean().reset_index()

snp_table['num_N'] = snp_table['Locustag'].map(gene_spec_N_cnt_dict)
snp_table['num_S'] = snp_table['Locustag'].map(gene_spec_S_cnt_dict)
snp_table['num_muts_genic'] = snp_table['num_S'] + snp_table['num_N']

## map probability of nonsyn muts and calculate dnds score

snp_table_genic = snp_table[snp_table['num_muts_genic'] > 0]
snp_table_genic['prob_N'] = snp_table_genic['locustag'].map(prob_N_per_gene_dict).astype(float)
snp_table_genic['dN'] = snp_table_genic['num_N'] / (snp_table_genic['num_muts_genic']*snp_table_genic['prob_N'])
snp_table_genic['dS'] = snp_table_genic['num_S'] / (snp_table_genic['num_muts_genic']* (1-snp_table_genic['prob_N']) )
snp_table_genic['dnds'] = snp_table_genic['dN'] / snp_table_genic['dS']
## remove those which have nonsyn only for the plot
print(f'Identified {sum(np.isinf(snp_table_genic["dnds"]))} with nonsyn mutations only')
snp_table_genic = snp_table_genic[~np.isinf(snp_table_genic['dnds'])]
snp_table_genic = snp_table_genic.sort_values('dnds').reset_index(drop = True)
snp_table_genic['empirical_pval_of_dnds'] = snp_table_genic.index / max(snp_table_genic.index)
snp_table_genic['empirical_pval_of_dnds_perc'] = snp_table_genic['empirical_pval_of_dnds'] * 100 ## convert to percentages
snp_table_genic['dnds_log'] = np.log(snp_table_genic['dnds']).fillna(np.nanmin(np.log(snp_table_genic['dnds']))) ## Note: neutrality is at 0 then!

snp_table_genic['dnds_plt'] = np.where(np.isinf(snp_table_genic['dnds']), 
                                 snp_table_genic.loc[~np.isinf(snp_table_genic['dnds']), 'dnds'].max() + 1,
                                 snp_table_genic['dnds'])

#change_ytick_dict = {3: '', 
#                     3.5: 'Infinite'}

fimZ_locustag = 'MDIMIJ_05850' # fimZ
fim_locustags = ['MDIMIJ_05855', # sfmF
                'MDIMIJ_05860', # fimH
                'MDIMIJ_05865', # fimD
                'MDIMIJ_05870', # fimC
                'MDIMIJ_05875', # fimI
                'MDIMIJ_05880'] # fimA

fig, ax = plt.subplots(figsize = (8, 3))
ax.axhline(y = 1, ls = '--', lw = 0.5, color = 'k', zorder = 5) ## annotate neutrality
sns.scatterplot(data = snp_table_genic, x = 'empirical_pval_of_dnds_perc', y = 'dnds_plt', color = 'lightgrey', edgecolor = 'k', lw = 0.1, alpha = 0.4, zorder = 10, ax = ax)
sns.scatterplot(data = snp_table_genic[snp_table_genic['locustag'].isin(fim_locustags)], x = 'empirical_pval_of_dnds_perc', y = 'dnds_plt', color = 'orange', edgecolor = 'k', lw = 0.1, zorder = 15, ax = ax)
sns.scatterplot(data = snp_table_genic[snp_table_genic['locustag'] == fimZ_locustag], x = 'empirical_pval_of_dnds_perc', y = 'dnds_plt', color = 'red', edgecolor = 'k', lw = 0.1, zorder = 20, ax = ax)
ax.set_xlim((0, 100))
ax.set_ylim((None, None))
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
#ax.set_yticks(np.concatenate([np.arange(0, 3, 0.5), np.array([3.5])]))
#ax.set_yticklabels([l if l not in change_ytick_dict.keys() else change_ytick_dict[l] for l in ax.get_yticks()])
legend_info = [Line2D([0], [0], marker = 'o', color = 'red', markerfacecolor = 'red', markeredgecolor = 'k', markeredgewidth = 0.1, linestyle="None", markersize = 8, label = 'fimZ')] 
legend_info += [Line2D([0], [0], marker = 'o', color = 'orange', markerfacecolor = 'orange', markeredgecolor = 'k', markeredgewidth = 0.1, linestyle="None", markersize = 8, label = 'Gene in\nfim operon')] 
legend_info += [Line2D([0], [0], marker = 'o', color = 'lightgrey', markerfacecolor = 'lightgrey', markeredgecolor = 'k', markeredgewidth = 0.1, linestyle="None", markersize = 8, label = 'other')]
#legend_labels = ['fimZ', 'fim operon', 'other'] 
ax.legend(title = 'Genes', handles = legend_info, loc='center left', bbox_to_anchor=(1, 0.5), fancybox=False)
ax.set_ylabel('dNdS')
ax.set_xlabel('Ordered dNdS score percentile [%]')
ax.set_title(f'dNdS score of {len(snp_table_genic)} genes {len(smpl_cols)}\nEnterobacter hormaechei isolates')
plt.tight_layout()
fig.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2023/public_data/2023_12_analysis/analysis/2025_08_25_P07_Ehormaechei-c1_240825_gubbinsdefault_ambig10/pdf/2025_06_18_empirical_pval_dNdS_sorted.pdf'))

## print genes which are under positive selection: 
print(f"Genes under positive selection are:\n{snp_table_genic[['Locustag', 'gene', 'product', 'num_muts_genic', 'num_N', 'num_S', 'dnds', 'across_smpl_support']].sort_values('dnds', ascending = False).to_string(index = False)}")

print(f"FimZ is in the following dNdS percentiles:\n{snp_table_genic.loc[(snp_table_genic['locustag'] == fimZ_locustag), ['gene', 'num_muts_genic', 'num_N', 'num_S', 'dnds', 'empirical_pval_of_dnds_perc', 'across_smpl_support']].to_string(index = False)}")

print(f"Fim operon genes are in the following dNdS percentiles:\n{snp_table_genic.loc[snp_table_genic['locustag'].isin(fim_locustags), ['gene', 'num_muts_genic', 'num_N', 'num_S', 'dnds', 'empirical_pval_of_dnds_perc', 'across_smpl_support']].to_string(index = False)}")