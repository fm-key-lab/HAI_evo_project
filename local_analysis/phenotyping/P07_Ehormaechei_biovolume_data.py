

## analyse biovolume data and CV data from carey Nadell 2024/02/20

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import seaborn as sns
import scipy.stats as stats
from statannotations.Annotator import Annotator
from statannotations.stats.StatTest import StatTest
import itertools

plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'


save_plt = False

###############################################################
###############################################################
###############################################################
## load data 
df = pd.read_excel('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/biofilm_formation/2024_02_20_flow_biofilm_CV_carey_nadell.xlsx',
                sheet_name=1)

## extract the strain names 
strain_names = df.loc[6:, 'Time (h)'].str.split(' \\(').values
strain_dict = {strainid: sampleid.strip(')') for strainid, sampleid in strain_names}

data_entry_stop = np.where(df.isna().sum(axis = 1) == len(df.columns))[0][0]
df = df[:data_entry_stop]
## rename columns 
colname = ''
colnames = []
for colname in df.columns:
    if colname.startswith('Strain '):
        strainname = colname
        i = 1
        colname_new = strainname + '.' + str(i)
    elif colname.startswith('Unnamed: '):
        i += 1
        colname_new = strainname + '.' + str(i)
    else:
        colname_new = colname
    colnames.append(colname_new)
df.columns = colnames

df_long = df.melt('Time (h)')
df_long['strain'] = df_long['variable'].str.split('.').str[0]
df_long['sampleid'] = df_long['strain'].map(strain_dict)

###############################################################
###############################################################
###############################################################

df_cv = pd.read_excel('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/biofilm_formation/2024_02_20_flow_biofilm_CV_carey_nadell.xlsx',
                sheet_name=0)

## extract the strain names 
strain_names = df_cv.loc[6:, 'Strain 1'].str.split(' \\(').values
strain_dict = {strainid: sampleid.strip(')') for strainid, sampleid in strain_names}
data_entry_stop = np.where(df_cv.isna().sum(axis = 1) == len(df_cv.columns))[0][0]
df_cv = df_cv[:data_entry_stop]
df_cv.rename(strain_dict, axis =1, inplace = True)
df_cv['replicate'] = df_cv.index+1
df_cv_long = df_cv.melt('replicate')
df_cv_long['value'] = df_cv_long['value'].astype(float)

############################
## Plots of interest 
############################

## group values together
translate_sampleids = { 'P07 2048 S3': 'wt',
                        'P07 2055 S3': 'wt',
                        'P07 2052 S3': 'FimZ[F126L]',
                        'P07 2237 S3': 'FimZ[F126L]',
                        'L00071 01 1 S3': 'FimZ[F126L]',
                        'L00071 02 1 S3': 'FimZ[F126L]'}

genotype_col_dict = {'wt': '#c8c8c8', 
                    'FimZ[F126L]': '#f9af92'}


## get only the isolates of interest used across experiments 
df_long = df_long[df_long['sampleid'].isin(translate_sampleids.keys())]
df_long['genotype'] = df_long['sampleid'].map(translate_sampleids)
df_cv_long = df_cv_long[df_cv_long['variable'].isin(translate_sampleids.keys())]
df_cv_long['genotype'] = df_cv_long['variable'].map(translate_sampleids)

## map the genotype to the data 
df_long['genotype'] = df_long['sampleid'].map(translate_sampleids)
df_long['genotype'] = pd.Categorical(df_long['genotype'], ['wt', 'FimZ[F126L]'], ordered = True)

## generate the test function (import from stats.mannwhitneyus)
mannwhitneyu_long_name = 'mannwhitneyu test'
mannwhitneyu_short_name = ''
mannwhitneyu_func = stats.mannwhitneyu
mannwhitneyu_test = StatTest(mannwhitneyu_func, mannwhitneyu_long_name, mannwhitneyu_short_name, alternative = 'two-sided')


isolate_pairs = list(itertools.combinations(set(list(translate_sampleids.values())), 2))

for pval_method in ['star', 'text']:
    
    np.random.seed(123)
    ####################
    ## plot biovolume
    fig, ax = plt.subplots(figsize = (3.5, 4))
    sns.boxplot(data = df_long[df_long['Time (h)'] == 24], 
                x = 'genotype', y = 'value', hue = 'genotype',
                palette = genotype_col_dict, legend = True,
                fliersize = 0, dodge = False,
                ax = ax)
    sns.stripplot(data = df_long[df_long['Time (h)'] == 24], 
                    x = 'genotype', y = 'value', hue = 'genotype',
                    palette = genotype_col_dict, 
                    alpha = 0.8, dodge = False,
                    s = 8, linewidth = 0.4, edgecolor = 'k', 
                    ax = ax)
    annot = Annotator(ax, isolate_pairs, data=df_long[df_long['Time (h)'] == 24], x='genotype', y='value')
    if pval_method == 'star': 
        annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=True, line_height=0.01, text_offset=-1, pvalue_thresholds=[[1e-4, '****'], [1e-3, '***'], [1e-02, '**'], [0.05, '*'], [1, 'ns']]).apply_test()
    else:
        annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=True, line_height=0.01, text_offset=-1, text_format='full', pvalue_format_string='{:.1e}').apply_test()
    annot.apply_test()
    annot.annotate(line_offset_to_group=0.04, line_offset=-1)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(10000))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(5000))
    ax.grid(which='major', axis = 'y', linestyle='-', linewidth=0.4, color='black', zorder = 1)
    ax.grid(which='minor', axis = 'y', linestyle=':', linewidth=0.4, color='black', zorder = 1)
    ax.set_ylim((-1000, None))
    ax.set_ylabel('Biovolume [µm^3]')
    ax.set_xlabel(None)
    h, l = ax.get_legend_handles_labels() ## get legend handles and labels
    len_legend = int(len(h)/2) ## cut legend in half --> second half are strip plot values which are redundant! --> remove
    ax.legend(title = 'Genotype',loc = 'center left', bbox_to_anchor = (1, 0.5), 
            handles = h[:], labels = l[:])
    plt.tight_layout()
    if save_plt:
        fig.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/biofilm_formation/2024_05_29_flow_biofilm_CN_loc_genotype_24h_pval{pval_method}.svg'))

    ####################
    ## plot crystal violet
    
    fig, ax = plt.subplots(figsize = (3.5, 4), sharey = True)
    sns.boxplot(data = df_cv_long, 
                x = 'genotype', y = 'value', hue = 'genotype',
                palette = genotype_col_dict, legend = True,
                fliersize = 0, dodge = False,
                ax = ax)
    sns.stripplot(data = df_cv_long, 
                    x = 'genotype', y = 'value', hue = 'genotype',
                    palette = genotype_col_dict, 
                    alpha = 0.8, dodge = False,
                    s = 8, linewidth = 0.4, edgecolor = 'k', 
                    ax = ax)
    annot = Annotator(ax, isolate_pairs, data=df_cv_long, x='genotype', y='value')
    if pval_method == 'star': 
        annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=True, line_height=0.01, text_offset=-1, pvalue_thresholds=[[1e-4, '****'], [1e-3, '***'], [1e-02, '**'], [0.05, '*'], [1, 'ns']]).apply_test()
    else:
        annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=True, line_height=0.01, text_offset=-1, text_format='full', pvalue_format_string='{:.1e}').apply_test()
    annot.apply_test()
    annot.annotate(line_offset_to_group=0.04, line_offset=-1)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.02))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
    ax.grid(which='major', axis = 'y', linestyle='-', linewidth=0.4, color='black', zorder = 1)
    ax.grid(which='minor', axis = 'y', linestyle=':', linewidth=0.4, color='black', zorder = 1)
    ax.set_ylim((0, None))
    h, l = ax.get_legend_handles_labels() ## get legend handles and labels
    len_legend = int(len(h)/2) ## cut legend in half --> second half are strip plot values which are redundant! --> remove
    ax.legend(title = 'Genotype',loc = 'center left', bbox_to_anchor = (1, 0.5), 
            handles = h[:], labels = l[:])

    ax.set_ylabel('OD550')
    ax.set_xlabel(None)
    plt.tight_layout()
    if save_plt:
        fig.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/biofilm_formation/2024_05_29_CV_CN_per_loc_genotype_collapsed_pval{pval_method}.svg'))