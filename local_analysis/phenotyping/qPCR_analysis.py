
########################
## Load modules 

import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np 
import scipy.stats as stats
from statannotations.Annotator import Annotator
from statannotations.stats.StatTest import StatTest

plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

########################
## If you want to save plots new, set to true
save_plt = True

########################
## Load all datasets 
## get all files starting with pattern
files = glob.glob(os.path.expanduser(f'~/Desktop/fimZ_phenotype_validation/2025_05_qPCR/2025_06*Ehorm*rpoDfimZfimA_P*.csv'))
files = [f for f in files if not f.endswith('_raw.csv')]

samples_to_ids_dict = {'1': 'P07_2048',
                       '2': 'P07_2055',
                       '4': 'P07_2052',
                       '5': 'P07_2237',
                       '6': 'L00071_01_1',
                       '8': 'L00071_02_1'}
map_muts = {'P07_2048': 'wt',
            'P07_2055': 'wt',
            'P07_2052': 'FimZ[F126L]',
            'P07_2237': 'FimZ[F126L]',
            'L00071_01_1': 'FimZ[F126L]',
            'L00071_02_1': 'FimZ[F126L]'}

primer_ordered = ['rpoD', 'fimZ', 'fimA']

## generate color map
genotype_col_dict = {'wt': '#c8c8c8', 
                    'FimZ[F126L]': '#f9af92'}

remove_high_gDNA_contaminations = True ## if samples with a Ct difference < 5 between sample and RTNC should be removed (True) or not (False)

## read in all files and concatenate them 
all_data_df = pd.DataFrame()
for file in files:
    date = file.split('/')[-1][:10] ## get the date of each run to separate individual replicates
    plate = file.split('/')[-1].split('_')[-1].replace('.csv', '') ## get the date of each run to separate individual replicates
    growth_phase = file.split('/')[-1].split('_')[-3] ## get the date of each run to separate individual replicates
    tmp_df = pd.read_csv(file, sep = ';', decimal = ',')
    tmp_df['date'] = date
    tmp_df['plate'] = plate
    tmp_df['growth_phase'] = growth_phase
    all_data_df = pd.concat([all_data_df, tmp_df])

## map the genes
all_data_df['col'] = all_data_df['Well'].str[1:].astype(int)
all_data_df['row'] = all_data_df['Well'].str[0]
is_fimA_p1 = (all_data_df['plate'] == 'P1') & (all_data_df['row'].isin(['E', 'F']) | all_data_df['Well'].isin(['G1', 'G2', 'G3']))
is_fimZ_p1 = (all_data_df['plate'] == 'P1') & (all_data_df['row'].isin(['C', 'D']) | all_data_df['Well'].isin(['E1', 'E2']))
is_rpoD_p1 = (all_data_df['plate'] == 'P1') & (all_data_df['row'].isin(['A', 'B']) | all_data_df['Well'].isin(['B1']))
is_fimA_p2 = (all_data_df['plate'] == 'P2') & (all_data_df['row'].isin(['C', 'G']) | all_data_df['Well'].isin(['D1', 'D2', 'D3', 'H1', 'H2', 'H3']))
is_fimZ_p2 = (all_data_df['plate'] == 'P2') & (all_data_df['row'].isin(['B', 'F']) | all_data_df['Well'].isin(['C1', 'C2', 'G1', 'G2']))
is_rpoD_p2 = (all_data_df['plate'] == 'P2') & (all_data_df['row'].isin(['A', 'E']) | all_data_df['Well'].isin(['B1', 'F1']))
is_fimA_PE = (all_data_df['plate'].str.contains('Prim')) & (all_data_df['row'].isin(['G', 'H']) | all_data_df['col'].isin([12]) | all_data_df['Well'].isin(['F11', 'G11', 'H11']))
is_fimZ_PE = (all_data_df['plate'].str.contains('Prim')) & ((all_data_df['row'].isin(['D', 'E', 'F']) & all_data_df['col'].isin(np.arange(1, 11))) | all_data_df['Well'].isin(['E11']))
is_rpoD_PE = (all_data_df['plate'].str.contains('Prim')) & ((all_data_df['row'].isin(['A', 'B', 'C']) & all_data_df['col'].isin(np.arange(1, 11))) | all_data_df['Well'].isin(['B11']))
all_data_df.loc[is_fimA_p1 | is_fimA_p2 | is_fimA_PE, 'Gene'] = 'fimA'
all_data_df.loc[is_fimZ_p1 | is_fimZ_p2 | is_fimZ_PE, 'Gene'] = 'fimZ'
all_data_df.loc[is_rpoD_p1 | is_rpoD_p2 | is_rpoD_PE, 'Gene'] = 'rpoD'

## map growth phase
is_early =  (all_data_df['plate'] == 'P2') & (all_data_df['row'].isin(['A', 'B', 'C', 'D']))
is_late =  (all_data_df['plate'] == 'P2') & (all_data_df['row'].isin(['E', 'F', 'G', 'H']))
all_data_df.loc[is_early & (all_data_df['growth_phase'] == 'earlyAndlatelog'), 'growth_phase'] = 'earlylog'
all_data_df.loc[is_late & (all_data_df['growth_phase'] == 'earlyAndlatelog'), 'growth_phase'] = 'latelog'

all_data_df = all_data_df[all_data_df['Exclusion'].isna()]

all_data_df.reset_index(drop = True, inplace = True)

## set undetermined ct values to 40 to allow calculations on them!
all_data_df.loc[all_data_df['Ct'].isna(), 'Ct'] = 40
all_data_df['Ct'] = all_data_df['Ct'].astype(float)

## annotate data
all_data_df['sampleid'] = all_data_df['Sample'].str[0].map(samples_to_ids_dict)
all_data_df['replicate'] = all_data_df['Sample'].str[1]
all_data_df['sample_rep'] = all_data_df['Sample'].str[0] + '.' + all_data_df['replicate']
all_data_df.loc[all_data_df['replicate'] == 'T', 'replicate'] = ''
all_data_df['sample_type'] = np.where(all_data_df['Sample'].str.contains('-'), 'RTNC', 
                                      np.where(all_data_df['Sample'] == 'NTC', 'NTC', all_data_df['sample_rep']))

########################
## Primer efficiency 
########################
## measured rpoD(783-890) & fimZ(466-595): 24/12/15; fimA(246-369): 23/08/21
efficiency_dict = {'fimZ': 103.12860110426207, 'rpoD': 106.90652565606507, 'fimA': 104.12957565472598}
r_val_dict = {'fimZ': 0.992816941116989, 'rpoD': 0.9882236015156362, 'fimA': 0.9970334188582403}

############################################################
############################################################
############################################################
## analyse data 
df_samples = all_data_df[~all_data_df['Sample'].isna()]

## extract different sample types from df
df_smpl = df_samples.loc[~df_samples['sample_type'].isin(['RTNC', 'NTC']), ['sampleid', 'sample_rep', 'Gene', 'Ct', 'growth_phase']]

df_ctrs = df_samples[df_samples['sample_type'].isin(['RTNC', 'NTC'])]
df_ntc = df_ctrs.loc[df_ctrs['sample_type'] == 'NTC', ['Gene', 'Ct', 'growth_phase']]
df_rtnc = df_ctrs.loc[df_ctrs['sample_type'] == 'RTNC', ['sampleid', 'sample_rep', 'Gene', 'Ct', 'growth_phase']]
df_rtnc['sample_rep'] = df_rtnc['sample_rep'].str.replace('-', '')

## compare RTNC against NTC
df_ctr_diff = pd.merge(df_rtnc, df_ntc, how = 'left', on = ['Gene', 'growth_phase'], suffixes=['_rtnc', '_ntc'])
df_ctr_diff['Ct_diff'] = df_ctr_diff['Ct_ntc'] - df_ctr_diff['Ct_rtnc']

## compare Samples against RTNC
df_rtnc_smpl_diff = pd.merge(df_smpl, df_rtnc, how = 'left', on = ['sampleid', 'sample_rep', 'Gene', 'growth_phase'], suffixes=['_smpl', '_rtnc'])
df_rtnc_smpl_diff['Ct_diff'] = df_rtnc_smpl_diff['Ct_rtnc'] - df_rtnc_smpl_diff['Ct_smpl']

## remove all samples if their houskeeping gene ct difference between sample and RT- control is >5 (--> gDNA contamination is <3%!)
## cutoff based on: https://pmc.ncbi.nlm.nih.gov/articles/PMC3326333/
df_rtnc_smpl_diff_low_gDNA = df_rtnc_smpl_diff.loc[(df_rtnc_smpl_diff['Gene'] == 'rpoD') & (df_rtnc_smpl_diff['Ct_diff'] >= 5), ['sampleid', 'sample_rep', 'growth_phase', 'Ct_diff']].drop_duplicates()
df_rtnc_smpl_diff_filt = pd.merge(df_rtnc_smpl_diff, df_rtnc_smpl_diff_low_gDNA[['sampleid', 'sample_rep', 'growth_phase']], how = 'right', on = ['sampleid', 'sample_rep', 'growth_phase'])
df_rtnc_smpl_diff_high_gDNA = df_rtnc_smpl_diff.loc[(df_rtnc_smpl_diff['Gene'] == 'rpoD') & (df_rtnc_smpl_diff['Ct_diff'] < 5), ['sampleid', 'sample_rep', 'growth_phase', 'Ct_diff']].drop_duplicates()
df_rtnc_smpl_diff_rmvd = pd.merge(df_rtnc_smpl_diff, df_rtnc_smpl_diff_high_gDNA[['sampleid', 'sample_rep', 'growth_phase']], how = 'right', on = ['sampleid', 'sample_rep', 'growth_phase'])

print(f"Removed {len(df_rtnc_smpl_diff_high_gDNA)}/{sum(df_rtnc_smpl_diff['Gene'] == 'rpoD')} samples due to too high gDNA load:\n{df_rtnc_smpl_diff_high_gDNA.sort_values(['growth_phase', 'sample_rep', 'Ct_diff']).to_string(index = False)}")

## normalize values by removing gDNA background contamination 
if remove_high_gDNA_contaminations:
    df_rtnc_smpl_diff_corrected = df_rtnc_smpl_diff_filt.copy()
else:
    df_rtnc_smpl_diff_corrected = df_rtnc_smpl_diff.copy()

df_rtnc_smpl_diff_corrected['Ct_smpl'] = df_rtnc_smpl_diff_corrected.apply(lambda x: -np.log2( 2**( -x['Ct_smpl'] ) - 2**( -x['Ct_rtnc'] ) ), axis = 1)

########
## calculate delta ct values
dct_smpls = df_rtnc_smpl_diff_corrected.copy()
dct_smpls['efficiency'] = dct_smpls['Gene'].map(efficiency_dict) / 100 ## map calculated efficiencies and return to ratios

## calculate mean ct value per biological replicate and then per sample id
dct_smpls_mean_ct_smpl = df_rtnc_smpl_diff_corrected.groupby(['sample_rep', 'sampleid', 'Gene', 'growth_phase'])['Ct_smpl'].agg('mean').reset_index()
## get mean ct value per primer and growth stage for wildtype
ref_smpls_mean_ct_smpl = dct_smpls_mean_ct_smpl.loc[dct_smpls_mean_ct_smpl['sampleid'].isin(['P07_2048', 'P07_2055']), ['Gene', 'growth_phase', 'Ct_smpl']].groupby(['Gene', 'growth_phase']).mean().reset_index()

dct_smpls = pd.merge(dct_smpls, ref_smpls_mean_ct_smpl, how = 'left', on = ['Gene', 'growth_phase'], suffixes = ['', '_ref'])

dct_smpls['dct_diff'] = dct_smpls['Ct_smpl_ref'] - dct_smpls['Ct_smpl']

housekeeping_rpoD = dct_smpls[dct_smpls['Gene'] == 'rpoD']
housekeeping_rpoD.drop(['Gene'], axis = 1, inplace = True)
housekeeping_rpoD.drop_duplicates(inplace = True)

## merge data back to df
ddct_smpls = pd.merge(dct_smpls, housekeeping_rpoD, how = 'left', on = ['sample_rep', 'growth_phase'], suffixes=['', '_hk_rpoD'])

## calculate fold-change per isolate and biological replicate
ddct_smpls['fc_rpoD_raw'] = ( (2*ddct_smpls['efficiency']) ** (ddct_smpls['dct_diff']) ) / ( (2*ddct_smpls['efficiency_hk_rpoD']) ** (ddct_smpls['dct_diff_hk_rpoD']) )
ddct_smpls['l2fc_rpoD_raw'] = np.log2(ddct_smpls['fc_rpoD_raw'])

## get mean per isolate and condition
ddct_smpls_grouped = ddct_smpls.groupby(['sample_rep', 'sampleid', 'Gene', 'growth_phase'], as_index = False)[['Ct_smpl', 'Ct_rtnc', 'Ct_diff', 'efficiency', 'Ct_smpl_ref', 'dct_diff', 'Ct_smpl_hk_rpoD', 'Ct_rtnc_hk_rpoD', 'Ct_diff_hk_rpoD', 'efficiency_hk_rpoD', 'Ct_smpl_ref_hk_rpoD', 'dct_diff_hk_rpoD', 'fc_rpoD_raw']].agg('mean')
ddct_smpls_grouped['l2fc_rpoD_raw'] = np.log2(ddct_smpls_grouped['fc_rpoD_raw'])

## map different sample collections together
ddct_smpls['mut_group'] = ddct_smpls['sampleid'].map(map_muts)
ddct_smpls = ddct_smpls.sort_values('sampleid').reset_index(drop = True)

ddct_smpls_grouped['mut_group'] = ddct_smpls_grouped['sampleid'].map(map_muts)
ddct_smpls_grouped = ddct_smpls_grouped.sort_values('sampleid').reset_index(drop = True)

####################################
## group over all growth phases and locations!!
####################################

mannwhitneyu_long_name = 'mannwhitneyu test'
mannwhitneyu_short_name = ''
mannwhitneyu_func = stats.mannwhitneyu
mannwhitneyu_test = StatTest(mannwhitneyu_func, mannwhitneyu_long_name, mannwhitneyu_short_name)

sorted_order_fim = ['rpoD', 'fimZ', 'fimA']
mutations = ['wt', 'FimZ[F126L]']

mutation_pairs = []
for gene in sorted_order_fim:
    mutation_pairs.append(((gene, mutations[0]), (gene, mutations[1])))

boxplt_df = ddct_smpls_grouped.copy()


boxplt_fim_df = ddct_smpls_grouped[ddct_smpls_grouped['Gene'].isin(sorted_order_fim)]
boxplt_fim_df['x_axis'] = pd.Categorical(boxplt_fim_df['Gene'], sorted_order_fim, ordered = True)

for pval_method in ['star', 'text']:
    fig, ax = plt.subplots(figsize = (4.5, 4))
    ax.axhline(y = 0, lw = 1, ls = '--', c = 'grey')
    sns.boxplot(data = boxplt_fim_df, x = 'x_axis', y = 'l2fc_rpoD_raw', hue = 'mut_group', 
                    hue_order = genotype_col_dict.keys(), palette = genotype_col_dict, 
                    linewidth = 1, fliersize = 0, zorder = 10, ax = ax)
    sns.stripplot(data = boxplt_fim_df, x = 'x_axis', y = 'l2fc_rpoD_raw', hue = 'mut_group', 
                    hue_order = genotype_col_dict.keys(), palette = genotype_col_dict, 
                    alpha = 0.8, linewidth = 1, dodge = True, zorder = 20, ax = ax)
    ## add significance vals
    annot = Annotator(ax, mutation_pairs, data=boxplt_fim_df, x='x_axis', y='l2fc_rpoD_raw', order=sorted_order_fim, hue='mut_group', hue_order=genotype_col_dict.keys())
    if pval_method == 'star': 
        annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=True, line_height=0.01, text_offset=-1, pvalue_thresholds=[[1e-4, '****'], [1e-3, '***'], [1e-02, '**'], [0.05, '*'], [1, 'ns']]).apply_test()
    else:
        annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=True, line_height=0.01, text_offset=-1, text_format='full', pvalue_format_string='{:.1e}').apply_test()
    annot.apply_test()
    annot.annotate(line_offset_to_group=0.04, line_offset=-1)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.grid(which='major', axis = 'y', linestyle='-', linewidth=0.4, color='black', zorder = 1)
    ax.grid(which='minor', axis = 'y', linestyle=':', linewidth=0.4, color='black', zorder = 1)
    ax.set_ylabel('Fold change gene expression\n(log2)')
    ax.set_xlabel(None)
    h, l = ax.get_legend_handles_labels() ## get legend handles and labels
    len_legend = int(len(h)/2) ## cut legend in half --> second half are strip plot values which are redundant! --> remove
    ax.legend(title = 'Genotype',loc = 'center left', bbox_to_anchor = (1, 0.5), 
            handles = h[:len_legend], labels = l[:len_legend])
    plt.xticks(fontstyle = 'italic')
    plt.tight_layout()
    if save_plt:
        plt.savefig(os.path.expanduser(f'~/Desktop/fimZ_phenotype_validation/2025_05_qPCR/proc/plots/2025_06_12_RTqPCR_P07_Ehorm_fimZ_pval{pval_method}.svg'))

## save the final df 
if save_plt:
    boxplt_fim_df.to_csv(os.path.expanduser('~/Desktop/fimZ_phenotype_validation/2025_05_qPCR/proc/qPCR_analysis_values.csv'), index = False)