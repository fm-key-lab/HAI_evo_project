
## analyse infection assays 

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import scipy.stats as stats
from statannotations.Annotator import Annotator
from statannotations.stats.StatTest import StatTest
import itertools

plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

## sanity check of column names (and columns)
def sanity_check_df(df, cols_to_melt, date = False):
    ## sanity check and convert to long df shape
    df.columns = [col.strip() for col in df.columns]
    for col in df.columns[:3]:
        df[col] = df[col].astype(str).str.strip()
    ## generate long df 
    df_long = df.melt(cols_to_melt)
    if date:
        df_long['date'] = date
    return df_long

## read in data
datafiles_dict = {'2024_01_17_adherenceassay_cfu_wt_vs_fimZ_mod.xlsx': '240117',
                  '2024_01_18_adherenceassay_cfu_wt_vs_fimZ_mod.xlsx': '240118',
                  '2024_01_19_adherenceassay_cfu_wt_vs_fimZ_mod.xlsx': '240119',
                  '2024_02_13_4_adherenceassay_cfu_wt_vs_fimZ_mod.xlsx': '240213',
                  '2024_02_14_4_adherenceassay_cfu_wt_vs_fimZ_mod.xlsx': '240214_1',
                  '2024_02_14_adherenceassay_cfu_wt_vs_fimZ_mod.xlsx': '240214_2'}

cols_to_melt = ['Well', 'Strain', 'replicate', 'cfu_replicate', 'step', 'cells_bacteria']

translate_sampleids = {'P07_2048': 'wt',
                        'P07_2055': 'wt',
                        'P07_2052': 'FimZ[F126L]',
                        'P07_2237': 'FimZ[F126L]',
                        'L00071_01_1': 'FimZ[F126L]',
                        'L00071_02_1': 'FimZ[F126L]'}
strain_order = ['P07_2048', 'P07_2055', 'P07_2052', 'P07_2237', 'L00071_01_1', 'L00071_02_1']
genotype_col_dict = {'wt': '#c8c8c8', 
                    'FimZ[F126L]': '#f9af92'}

df_l = []
for file, date in datafiles_dict.items():
    df = pd.read_excel(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/epithelial_adherence_assays/cfu/raw/{file}'), skiprows=1)
    if 'SampleID' in df.columns:
        df = df.drop('SampleID', axis = 1)
    if 'washing' in df.columns:
        df = df[~(df['cells_bacteria'].str.contains('calu3') & df['washing'].str.contains('vacuum'))].reset_index(drop = True) ## remove vacuum results as the vacuum was sucking all cells away (no cells post washing) --> random low cfu counts while pipetting showed reporducible results with cells being attached post washing!
        df = df.drop('washing', axis = 1)
    df = sanity_check_df(df, cols_to_melt, date)
    df_l.append(df)

df_long = pd.concat(df_l).reset_index(drop = True)
cols_to_melt += ['date']

## drop any columns with na or strings instead of values, then get the lowest dilution which still was countable 
df_long = df_long.sort_values('variable').reset_index(drop = True)
df_long = df_long[ df_long['value'].astype(str).str.replace('.', '').str.isdigit() & ~df_long['value'].isna() & (df_long['value'] != 0)]
df_long = df_long.fillna('nan')
df_long[['replicate', 'cfu_replicate']] = df_long[['replicate', 'cfu_replicate']].astype(float).astype(int)
df_long = df_long.groupby(cols_to_melt).head(1).reset_index(drop = True)

## translate the different experiments to replicates of experiments 
unique_strain_dates = df_long.sort_values('date')[['Strain', 'date']].drop_duplicates()
unique_strain_dates['exp_replicate'] = unique_strain_dates.groupby(['Strain']).cumcount() + 1
df_long = pd.merge(df_long, unique_strain_dates, how = 'left', on = ['Strain', 'date'])

## get the right dilutions in 
df_long['dilution'] = df_long['variable'].str[-1].astype(int) + 2
df_long['value'] = df_long['value'].astype(int) * 10**df_long['dilution']

## qc of replicated cfus
qc_df = df_long.copy()
qc_df = qc_df.sort_values(['date', 'Strain', 'replicate', 'cells_bacteria']).reset_index(drop = True)
qc_df['combined_x'] = qc_df['Strain'] + '_rep' + qc_df['replicate'].astype(str) + '_' + qc_df['cells_bacteria'] ## get a combined name for x axis 
 ## shorten the x axis names 
shorten_xaxis_dict = {'L00071_01_1': 'fimZ_lung1', 
                      'L00071_02_1': 'fimZ_lung2', 
                      'P07_2052': 'fimZ_rectal', 
                      'P07_2237': 'fimZ_oral', 
                      'P07_2055': 'wt_rectal1', 
                      'P07_2048': 'wt_rectal2', 
                      'bacteria_only': 'bac', 
                      'calu3_bacteria': 'wcalu3', 
                      '-': 'input'}
for key, val in shorten_xaxis_dict.items():
    qc_df['combined_x'] = qc_df['combined_x'].str.replace(key, val)
## calculate the offset of the mean per replicated cfu
qc_df['norm_cfu'] = qc_df.groupby(['date', 'combined_x'])['value'].transform(lambda x: x / x.mean())


## get average of replicated cfu assay
df_long = df_long.groupby(['Well', 'Strain', 'replicate', 'step', 'cells_bacteria', 'date', 'exp_replicate'])['value'].agg('mean').reset_index()


## get input values
input_vals = df_long.loc[df_long['step'] == 'input', ['date', 'Strain', 'value']]
input_vals = input_vals.groupby(['date', 'Strain']).mean().reset_index()
## get mean bacteria only to remove from 
baconly_vals = df_long.loc[df_long['cells_bacteria'] == 'bacteria_only', ['date', 'Strain', 'value']]
baconly_vals = baconly_vals.groupby(['date', 'Strain']).mean().reset_index()
df_long = pd.merge(df_long, input_vals, how = 'left', on = ['date', 'Strain'], suffixes=['', '_input'])
df_long = pd.merge(df_long, baconly_vals, how = 'left', on = ['date', 'Strain'], suffixes=['', '_bac_only'])
df_long['background_removed_value'] = df_long['value'] - df_long['value_bac_only']
df_long['rel_adherence_bg_rm'] = (df_long['background_removed_value'] / df_long['value_input']) * 100
df_long['rel_adherence'] = (df_long['value'] / df_long['value_input']) * 100
df_long['cells_date'] = df_long['cells_bacteria'] + '_' + df_long['date']
df_long['cells_exprep'] = df_long['exp_replicate'].astype(str) + ' (' + df_long['cells_bacteria'] + ')'

df_long_rel = df_long[df_long['step'] != 'input']
df_long_rel['Strain'] = pd.Categorical(df_long_rel['Strain'], strain_order, ordered = True)

## group values together
df_long_rel['genotype_loc'] = df_long_rel['Strain'].map(translate_sampleids)

df_long_rel_calu3_only = df_long_rel[df_long_rel['cells_bacteria'] == 'calu3_bacteria']
df_long_rel_calu3_only['genotype_loc'] = pd.Categorical(df_long_rel_calu3_only['genotype_loc'], ['wt', 'FimZ[F126L]'], ordered = True)

## generate the test function (import from stats.mannwhitneyu_test)
mannwhitneyu_long_name = 'mannwhitneyu test'
mannwhitneyu_short_name = ''
mannwhitneyu_func = stats.mannwhitneyu
mannwhitneyu_test = StatTest(mannwhitneyu_func, mannwhitneyu_long_name, mannwhitneyu_short_name)

isolate_pairs = list(itertools.combinations(set(list(translate_sampleids.values())), 2))

for pval_method in ['star', 'text']:
    fig, ax = plt.subplots(figsize = (3.5, 4), sharey = True)
    sns.boxplot(data = df_long_rel_calu3_only, 
                x = 'genotype_loc', y = 'rel_adherence_bg_rm', hue = 'genotype_loc',
                palette = genotype_col_dict, legend = True,
                fliersize = 0, dodge = False,
                ax = ax)
    sns.stripplot(data = df_long_rel_calu3_only, 
                    x = 'genotype_loc', y = 'rel_adherence_bg_rm', hue = 'genotype_loc',
                    palette = genotype_col_dict,
                    alpha = 0.8, dodge = False,
                    s = 8, linewidth = 0.4, edgecolor ='k', 
                    ax = ax)
    annot = Annotator(ax, isolate_pairs, data=df_long_rel_calu3_only, x='genotype_loc', y='rel_adherence_bg_rm')
    if pval_method == 'star': 
        annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=True, line_height=0.01, text_offset=-1, pvalue_thresholds=[[1e-4, '****'], [1e-3, '***'], [1e-02, '**'], [0.05, '*'], [1, 'ns']]).apply_test()
    else:
        annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=True, line_height=0.01, text_offset=-1, text_format='full', pvalue_format_string='{:.1e}').apply_test()
    annot.apply_test()
    annot.annotate(line_offset_to_group=0.04, line_offset=-1)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax.grid(which='major', axis = 'y', linestyle='-', linewidth=0.4, color='black', zorder = 1)
    ax.grid(which='minor', axis = 'y', linestyle=':', linewidth=0.4, color='black', zorder = 1)
    ax.set_ylim((0, None))
    ax.set_ylabel('Relative adherence [%]')
    ax.set_xlabel(None)
    h, l = ax.get_legend_handles_labels() ## get legend handles and labels
    ax.legend(title = 'Genotype',loc = 'center left', bbox_to_anchor = (1, 0.5), 
            handles = h, labels = l)
    plt.tight_layout()
    fig.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/2024_05_29_P07_Ehorm_adherence_assay_results_pval{pval_method}.svg'))
