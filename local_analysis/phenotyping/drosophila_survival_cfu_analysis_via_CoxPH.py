

## Analysis of Drosophila survivals - 2025

############################
## Load modules
############################
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import seaborn as sns
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
import scipy.stats as stats
from statannotations.Annotator import Annotator
from statannotations.stats.StatTest import StatTest
import itertools

############################
## set plotting style
############################
## change plotting style
plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

############################
## general variables
############################

## generate sample id to dictionary
id_to_genotype_dict = {'P07_2048': 'wt',
                       'P07_2055': 'wt',
                       'P07_2052': 'FimZ[F126L]',
                       'L00071_01_1': 'FimZ[F126L]'}
genotype_col_dict = {'wt': '#c8c8c8', 
                    'FimZ[F126L]': '#f9af92'}
isolate_col_dict = {'P07_2048': '#dddddd', 
                   'P07_2055': '#ababab',
                   'P07_2052': '#FAC4AF', ## 30% saturation of #f9af92
                   'L00071_01_1': '#FA8557' ## 65% saturation of #f9af92
                   }
date_to_exp_dict = {'2025_01_23': 'Experiment 1',
                    '2025_02_10': 'Experiment 2',
                    '2025_02_13': 'Experiment 3'}
## state path for the input files
drosophila_phenotype_path = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/Drosophila_survival_assay/')

## state the dilutions used for the CFU assays
cfu_initial_dilutionfactor = 10/500 # one fly was put in 500 µl PBS and 10 µl was plated = 10^0 or diluted 1:10 for the subsequent dilutions
cfu_stock_dilutionfactor = 10*1000 # 10 µl were plated and we want to know the CFU per nl

## set to true if you want to overwrite 
saveme = False
remove_signif = True

## remove significant different starting cultures:
## list od ['start_date', 'isolate', 'od', 'fly']
col_order_to_mask_surv = ['start_date', 'isolate', 'od', 'fly']
col_order_to_mask_cfu = ['date', 'sampleid', 'OD', 'fly']
## list in the order of the columns as above 
pairs_to_remove = [['2025-01-23', 'L00071_01_1', 5, 'RelishE20 iso'],
                    ['2025-02-10', 'P07_2048', 5, 'RelishE20 iso'],
                    ['2025-02-10', 'P07_2055', 5, 'RelishE20 iso'],
                    ['2025-01-23', 'L00071_01_1', 5, 'DrosDel w1118 iso'],
                    ['2025-02-10', 'P07_2048', 5, 'DrosDel w1118 iso']]
pairs_to_remove_tuple = set(map(tuple, pairs_to_remove))

## dict to overwrite max time for plotting of specific graphs
max_time_dict = {('RelishE20 iso', 5): 15,
                 ('DrosDel w1118 iso', 5): 175}

############################
## functions
############################

## General functions for data manipulation
def subset_df_to_group(df, **conditions):
    ## NOte: conditions is a dictionary of the variable name (key) and the variable (value)
    mask = pd.Series(True, index=df.index)  # generate a mask to keep all if no filter is applied
    for col, value in conditions.items(): ## loop over all conditions
        mask &= (df[col] == value)  # update filter given the condition
    return df[mask].copy()

## plotting specific functions
def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "ns"

## cfu plotting
def annotate_stats(plt_df_date, sampleid_col = 'sampleid', pval_method = 'star'):
    isolate_tp_pairs = []
    for tp in plt_df_date['timepoint'].unique():
        pairs = list(itertools.combinations(plt_df_date[sampleid_col].unique(), 2))
        isolate_tp_pairs.extend([((tp, p1), (tp, p2)) for p1, p2 in pairs])
    annot = Annotator(ax, isolate_tp_pairs, data = plt_df_date, x = 'timepoint', y = 'cfu_of_fly', hue = sampleid_col)
    if pval_method == 'star': 
        annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=True, line_height=0.01, text_offset=-1, pvalue_thresholds=[[1e-4, '****'], [1e-3, '***'], [1e-02, '**'], [0.05, '*'], [1, 'ns']]).apply_test()
    else:
        annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=True, line_height=0.01, text_offset=-1, text_format='full', pvalue_format_string='{:.1e}').apply_test()
    annot.annotate(line_offset_to_group=0.04, line_offset=-1)

def adjust_legend(ax, is_most_left_axis, isolate_mapping_dict = {}, title = 'Isolate'):
    if is_most_left_axis:
        h, l = ax.get_legend_handles_labels()
        if isolate_mapping_dict != {}:
            l = [f'{isolate} ({isolate_mapping_dict[isolate]})' for isolate in l]
        legend_length = int(len(h)/2)
        ax.legend(handles=h[:legend_length], labels=l[:legend_length],
            title = title, loc = 'center left', bbox_to_anchor = (1, 0.5))
    else:
        ax.legend().remove()

def annotate_figure(ylabel, xlabel, title):
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.yaxis.grid(which = 'both', linewidth = 0.25, linestyle = '--', color = 'grey', zorder = 0)
    ax.set_axisbelow(True)
    ax.set_title(title)

## survival plotting functions 
def fit_kaplan_meier_curves(ax1, ax2, fun_death_table, kmf_group_vals_l, label_id_col_name, col_dict):
    for group_vars, _ in fun_death_table.groupby(kmf_group_vals_l):
        if len(kmf_group_vals_l) > 1:
            sampleid = np.array(group_vars)[np.array(kmf_group_vals_l) == label_id_col_name][0] ## get label of the group 
        else:
            sampleid = group_vars[0]
            group_vars = [group_vars[0]]
        ## fit kaplan meier curves on both axes
        kmf = KaplanMeierFitter()
        fun_death_table_subset = subset_df_to_group(fun_death_table, **dict(zip(kmf_group_vals_l, group_vars)))
        kmf.fit(fun_death_table_subset["time"], 
                event_observed= fun_death_table_subset["event"])
        ax1 = kmf.plot_survival_function(ax=ax1, label=sampleid, color = col_dict[sampleid], lw = 2, ci_show=False) ## remove CIs --> estimates using Greenwood’s Exponential formula; but we want to show observed data within that plot only!
        ax2 = kmf.plot_survival_function(ax=ax2, label=sampleid, color = col_dict[sampleid], lw = 2, ci_show=False) ## remove CIs --> estimates using Greenwood’s Exponential formula; but we want to show observed data within that plot only!
    return ax1, ax2

def annotate_cfu_tp(cfu_tps_fly_od_date, ax, xdatarange):
    [ax.axvline(x = cfu_tp, color = 'k', linestyle = '--', lw = 0.75, zorder = 0) for cfu_tp in cfu_tps_fly_od_date if (cfu_tp > xdatarange[0]) & (cfu_tp < xdatarange[1])] ## annotate CFU timepoints

def modify_axes(ax1, ax2, xdatarange_low, xdatarange_high, widths_ratio):
    ax1.legend().remove()
    ax1.set_xlim(xdatarange_low)
    ax2.set_xlim((xdatarange_high[0], xdatarange_high[1]))
    ax1.set_ylim((val/100 for val in ydatarange))
    ax2.set_ylim((val/100 for val in ydatarange))
    ax1.set_yticks([tick for tick in ax1.get_yticks() if (tick >= ydatarange[0]/100) & (tick <= ydatarange[1]/100)])
    ax2.set_yticks([tick for tick in ax1.get_yticks() if (tick >= ydatarange[0]/100) & (tick <= ydatarange[1]/100)])
    ax1.set_xticks([0])
    if xdatarange_high[1] >= 72:
        ax2.set_xticks(np.arange(np.ceil(xdatarange_high[0]/24)*24, np.ceil(xdatarange_high[1]), 24))
    elif xdatarange_high[1] <= 20:
        ax2.set_xticks(np.arange(np.ceil(xdatarange_high[0]/2)*2, np.ceil(xdatarange_high[1]), 2))
    else:
        ax2.set_xticks(np.arange(np.ceil(xdatarange_high[0]/8)*8, np.ceil(xdatarange_high[1]), 8))
    ax1.set_yticklabels([int(tick*100) for tick in ax1.get_yticks()])
    ax1.spines.right.set_visible(False)
    ax2.spines.left.set_visible(False)
    ax2.set_yticks([])
    d = .4  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([1, 1], [1, 0], transform=ax1.transAxes, **kwargs)
    ax2.plot([0, 0], [1, 0], transform=ax2.transAxes, **kwargs)
    ax1.set_ylabel('Survival rate [%]')
    ax1.set_xlabel(None)
    ax2.set_xlabel('Time [hours]')
    ax2.xaxis.set_label_coords(0.5*(widths_ratio[1]/sum(widths_ratio)), -0.1)

############################
## set stat test for CFU
############################

## generate the test function (import from stats.ranksums)
mannwhitneyu_long_name = 'mannwhitneyu test'
mannwhitneyu_short_name = ''
mannwhitneyu_func = stats.mannwhitneyu
mannwhitneyu_test = StatTest(mannwhitneyu_func, mannwhitneyu_long_name, mannwhitneyu_short_name)

############################
## read & combine data 
############################
## read in data
drosophila_survival_files = [f'{drosophila_phenotype_path}raw/2025_01_23_survivals_and_cfu.xlsx', f'{drosophila_phenotype_path}raw/2025_02_10_survivals_and_cfu.xlsx', f'{drosophila_phenotype_path}raw/2025_02_13_survivals_and_cfu.xlsx']
survival_curves_l = []
cfus_l = []
for file in drosophila_survival_files:
    exp_dat = file.split('/')[-1].split('_survivals_')[0]
    survival_curve = pd.read_excel(file, sheet_name = f'{exp_dat}_survivals')
    survival_curve['start_date'] = survival_curve['start_date'].astype(str)
    survival_curves_l.append(survival_curve)
    cfu = pd.read_excel(file, sheet_name = f'{exp_dat}_cfus', header = 1)
    cfu['date'] = cfu['date'].astype(str)
    cfus_l.append(cfu)
survival_curves = pd.concat(survival_curves_l).reset_index(drop = True)
cfus = pd.concat(cfus_l).reset_index(drop = True)

## rename rel to relish
survival_curves.loc[survival_curves['fly'] == 'rel', 'fly'] = 'RelishE20 iso'
cfus.loc[cfus['fly'] == 'rel', 'fly'] = 'RelishE20 iso'
survival_curves.loc[survival_curves['fly'] == 'wt', 'fly'] = 'DrosDel w1118 iso'
cfus.loc[cfus['fly'] == 'wt', 'fly'] = 'DrosDel w1118 iso'


############################
## process CFUs
############################
melting_cols = ['date', 'sampleid', 'fly', 'OD', 'fly_rep', 'cfu_rep']
cfu_long = cfus.melt(melting_cols)
cfu_long = cfu_long[~cfu_long['value'].isna() & (cfu_long['value'] != '>')]
cfu_long['timepoint'] = cfu_long['variable'].str.split('h_').str[0].astype(int)
cfu_long['dilution'] = cfu_long['variable'].str.split('h_').str[-1].astype(float)


## select the highest value per timepoint
cfu_long = cfu_long.sort_values('dilution', ascending = False).groupby(melting_cols + ['timepoint']).head(1)
cfu_long = cfu_long[cfu_long['value'] != 'mistake']

## save unique timepoints for survival curve annotations
cfu_tps = cfu_long[['date', 'sampleid', 'fly', 'OD', 'timepoint']].drop_duplicates()

cfu_long['cfu_of_fly'] = cfu_long['value'].astype(float) * 1/(cfu_initial_dilutionfactor * cfu_long['dilution'])
cfu_long = cfu_long.groupby(['date', 'sampleid', 'fly', 'OD', 'fly_rep', 'timepoint'])['cfu_of_fly'].mean().reset_index() ## mean of the 3 replicates per CFU
cfu_long['genotype'] = cfu_long['sampleid'].map(id_to_genotype_dict)

cfu_long['sampleid'] = pd.Categorical(cfu_long['sampleid'], list(isolate_col_dict.keys()), ordered = True)
cfu_long['genotype'] = pd.Categorical(cfu_long['genotype'], ['wt', 'FimZ[F126L]'], ordered = True)

cfu_long_flies = cfu_long[cfu_long['fly'] != 'stock']
cfu_long_stock = cfu_long[cfu_long['fly'] == 'stock'].copy()
cfu_long_stock = cfu_long_stock.reset_index(drop = True)
cfu_long_stock['cfu_of_fly'] = cfu_long_stock['cfu_of_fly']*cfu_initial_dilutionfactor / cfu_stock_dilutionfactor ## correct dilution factor
cfu_long_stock['timepoint'] = 'stock'

############################
## identify significantly enriched datapoints at 0h 
############################

p_vals = {}
for group, grouped_df in cfu_long_flies.groupby(['date', 'fly', 'OD', 'timepoint']):
    #grouped_df_tp0 = grouped_df[grouped_df['timepoint'] == 0]
    p_vals[group] = {}
    # pairwise Mann-Whitney U tests
    for sample1, sample2 in itertools.combinations(grouped_df['sampleid'].unique(), 2):
        cfu_s1 = grouped_df.loc[(grouped_df['sampleid'] == sample1), 'cfu_of_fly']
        cfu_s2 = grouped_df.loc[(grouped_df['sampleid'] == sample2), 'cfu_of_fly']
        stat, p = stats.mannwhitneyu(cfu_s1, cfu_s2, alternative='two-sided')
        p_vals[group][(sample1, sample2)] = p

if saveme:
    with open(f'{drosophila_phenotype_path}/proc/cfus/stats/P07_Ehorm_drosophila_CFU_per_isolate_mannwhitneyus_pvals.csv', 'w') as fid:
        fid.write('date,fly,OD,timepoint,sample1,sample2,p_val\n')
        for (date, fly, OD, timepoint), sub_pval_dict in p_vals.items():
            for (s1, s2), pval in sub_pval_dict.items():
                fid.write(f'{date},{fly},{OD},{timepoint},{s1},{s2},{pval}\n')


############################
## PLOT CFUs
############################
for fly, od in cfu_long_flies[['fly', 'OD']].drop_duplicates().values:
    plt_df = cfu_long_flies[(cfu_long_flies['fly'] == fly) & (cfu_long_flies['OD'] == od)]
    num_days = len(np.unique(plt_df['date']))

    for pval_show in ['text', 'star']:
    ## plot per isolate 
        fig, axs = plt.subplots(figsize = (2+3.25*num_days,4.5), ncols = num_days, sharey = True)
        if num_days == 1:
            axs = [axs]
        for idx, (date, ax) in enumerate(zip(plt_df['date'].unique(), axs)):
            plt_df_date = plt_df[plt_df['date'] == date]

            sns.boxplot(data = plt_df_date, x = 'timepoint', y = 'cfu_of_fly', hue = 'sampleid', palette = isolate_col_dict, fliersize = 0, zorder = 10, ax = ax)
            sns.stripplot(data = plt_df_date, x = 'timepoint', y = 'cfu_of_fly', hue = 'sampleid', palette = isolate_col_dict, alpha = 0.6, size = 6, edgecolor = 'k', linewidth = 0.5, dodge = True, zorder = 20, ax = ax)
            #ax.set_yscale('log') ## need to be done before anootation of stats!
            ax.set_yscale('symlog', linthresh=5)
            ax.set_ylim((0, None))
            annotate_stats(plt_df_date, 'sampleid', pval_show)
            max_val_axis = int(np.ceil(np.log10(ax.get_ylim()[1]))) ## the max log base of the axis to set minor ticks in that range
            ax.yaxis.set_minor_locator(FixedLocator([k * 10**n for n in range(max_val_axis) for k in range(1, 10)])) ## set minor ticks for all integer values (remove minor ticks < 1!)
            adjust_legend(ax, (idx+1 == num_days), id_to_genotype_dict, title = 'Isolate (Genotype)')
            annotate_figure(ylabel = 'CFU/fly', xlabel = 'Time since infection [h]', title = date_to_exp_dict[date.replace('-', '_')])
        plt.suptitle(f'CFU per isolate\nof single {fly} flies\ninfected with OD {int(od)}')
        plt.tight_layout()
        if saveme:
            fig.savefig(f'{drosophila_phenotype_path}/proc/cfus/plots/P07_Ehorm_drosophila_{fly}_survival_CFUs_od{od}_per_isolate_pval{pval_show}.pdf', 
                        bbox_inches = 'tight')


############################
## add metadata to survivals
############################

## drop empty rows
survival_curves = survival_curves[~survival_curves['dead_flies'].isna()].reset_index(drop = True)

## modify date and number of dead flies 
survival_curves['start'] = pd.to_datetime(survival_curves['start_date'].astype(str) + ' ' + survival_curves['start_time'].astype(str))
survival_curves['measurement'] = pd.to_datetime(survival_curves['measure_date'].astype(str) + ' ' + survival_curves['measure_time'].astype(str))
survival_curves['time_hours'] = (survival_curves['measurement'] - survival_curves['start']).dt.total_seconds() / 3600
survival_curves['num_flies_recov'] = survival_curves['num_flies'] - survival_curves['non_recov']
survival_curves['num_flies_died_exp'] = survival_curves['dead_flies'] - survival_curves['non_recov']
survival_curves['genotype'] = survival_curves['isolate'].map(id_to_genotype_dict)

############################
## remove significant enriched data 
############################

if remove_signif:
    ## just remove significant from survivals (final results) but show differences in cfus!!
    suff_label = '_remove_signif_start'
    survival_curves = survival_curves[~survival_curves[col_order_to_mask_surv].apply(tuple, axis=1).isin(pairs_to_remove_tuple)]
    ## drop experiments with just one genotype left 
    exp_to_remove_dict = survival_curves.groupby(['start_date', 'fly', 'od'])['genotype'].nunique()
    exp_to_remove_dict = exp_to_remove_dict[exp_to_remove_dict <= 1].to_dict()
    for date, fly, od in exp_to_remove_dict.keys():
        survival_curves = survival_curves[~((survival_curves['start_date'] == date) & (survival_curves['fly'] == fly) & (survival_curves['od'] == od))]
else:
    suff_label = ''

############################
## process survivals
############################

## filter to the set maximal time of the assays 
survival_curves_l = []
for key, filter_val in max_time_dict.items():
    subset_df = subset_df_to_group(survival_curves, fly = key[0], od = key[1])
    survival_curves_l.append(subset_df[subset_df['time_hours'] <= filter_val])
survival_curves = pd.concat(survival_curves_l).reset_index(drop = True)

## cut survivals after 15 h for relish 
survival_curves = survival_curves[~((survival_curves['fly'] == 'rel') & (survival_curves['time_hours'] > 15))].reset_index(drop = True)

## fit data for kaplan meier extimator 
## net table with number of individuals and their length of survival 
group_cols = ['start_date', 'start_time', 'fly', 'inf_type', 'fly_sex', 'isolate', 'replicate', 'od', 'temp']
survival_curves['num_died_prev_period'] = survival_curves.sort_values('time_hours').groupby(group_cols)['num_flies_died_exp'].transform(lambda x: x - x.shift(1))
survival_curves['num_died_prev_period'] = survival_curves['num_died_prev_period'].fillna(0).astype(int)

kmf_group_vals = ['genotype', 'isolate', 'od', 'start', 'fly']
kmf_group_vals_grouped = ['genotype', 'od', 'start', 'fly']
death_table_list_dict = {}
for group_vars, group_df in survival_curves.groupby(kmf_group_vals): ## do it for every genotype individually to keep track of max time 
    death_table_list_dict[group_vars] = []
    max_time = max(group_df['time_hours'])
    

    for time in group_df['time_hours'].unique():
        row_of_tp = group_df['time_hours'] == time
        num_died_flies = group_df.loc[row_of_tp, 'num_died_prev_period'].sum()
        if num_died_flies > 0: ## check if not nan
            death_table_list_dict[group_vars].extend([[time, True] for i in range(num_died_flies)]) ## generate for every death during previous preiod a single event 
        ## check if some survived entire assay
        if time == max_time:
            num_flies_survived_assay = group_df.loc[row_of_tp, 'num_flies_recov'].sum() - group_df.loc[row_of_tp, 'num_flies_died_exp'].sum() 
            if num_flies_survived_assay > 0:
                death_table_list_dict[group_vars].extend([[max_time, False] for i in range(int(num_flies_survived_assay))]) ## generate for every survived fly a single "False" event

data_tuples = [(*key, *value) for key, values in death_table_list_dict.items() for value in values]
death_table_df = pd.DataFrame(data_tuples, columns= kmf_group_vals + ['time', 'event'])


############################
## PLOT Survivals per experiment
############################
for fly, od, date in survival_curves[['fly', 'od', 'start']].drop_duplicates().values:
    print(f'###\nFly genotype: {fly}; OD {od}; {date}')
    str_date = str(date)[:10].replace('-', '_')

    ## subset data to fly, od and date
    cfu_tps_fly_od_date = subset_df_to_group(cfu_tps, fly=fly, OD=od, date=str(date)[:10])
    cfu_tps_fly_od_date = cfu_tps_fly_od_date['timepoint'].unique()

    date_survival_df = subset_df_to_group(survival_curves, fly=fly, od=od, start=date)
    date_survival_df['isolate_rep'] = date_survival_df['isolate'] + '(rep. ' + date_survival_df['replicate'].astype(str) + ')'
    
    date_death_table = subset_df_to_group(death_table_df, fly=fly, od=od, start=date)
    
    ## calculate stats based on cox proportional hazards model 
    date_death_table["genotype_cat"] = date_death_table["genotype"].astype("category").cat.codes  
    gt_to_gt_cat_dict = {cat: gt for gt, cat in date_death_table[['genotype', 'genotype_cat']].drop_duplicates().values}

    ## calculate median flies which survived to set the annotation in the plot right
    max_time = date_survival_df['time_hours'].max()
    is_last_tp = date_survival_df['time_hours'] == max_time
    median_rate_surv_flies = date_survival_df[is_last_tp].groupby('genotype').apply(lambda x: 1- (x['num_flies_died_exp'] / x['num_flies_recov']).median())
    upper_val = max(median_rate_surv_flies)+0.02
    lower_val = min(median_rate_surv_flies)-0.02
    stat_annotaton = np.mean([upper_val, lower_val])+0.005
    if stat_annotaton == 0:
        stat_annotaton += 0.025
        upper_val += 0.025
        lower_val += 0.025
        
    ## select plotting boundaries 
    xdatarange_low = (-0.2, 0.5)
    max_time_plot = max_time
    if (fly, od) in max_time_dict.keys():
        max_time_plot = max_time_dict[(fly, od)]
    xdatarange_high = (7.5, max_time_plot*1.065)
    wspace = 0.05
    ydatarange = (-3, 103)
    if fly == 'DrosDel w1118 iso':
        xdatarange_low = (-0.2, 8)
        xdatarange_high = (22, max_time_plot*1.065)
    widths_ratio = [max-min for min, max in [xdatarange_low, xdatarange_high]]
    
    ############
    ## survivals grouped per isolate
    date_death_table_iso = date_death_table.copy()
    date_death_table_iso['genotype'] = pd.Categorical(date_death_table_iso['genotype'], ['wt', 'FimZ[F126L]'], ordered = True)
    date_death_table_iso['isolate'] = pd.Categorical(date_death_table_iso['isolate'], list(isolate_col_dict.keys()), ordered = True)

    fig, (ax1, ax2) = plt.subplots(figsize = (3.5, 3.5), ncols = 2, 
                                    gridspec_kw = {'width_ratios': widths_ratio, 'wspace': wspace})
    fit_kaplan_meier_curves(ax1, ax2, date_death_table_iso, kmf_group_vals, 'isolate', isolate_col_dict)
    annotate_cfu_tp(cfu_tps_fly_od_date, ax1, xdatarange_low)
    annotate_cfu_tp(cfu_tps_fly_od_date, ax2, xdatarange_high)

    modify_axes(ax1, ax2, xdatarange_low, xdatarange_high, widths_ratio)
    h, l = ax2.get_legend_handles_labels()
    l = [f'{isolate} ({id_to_genotype_dict[isolate]})' for isolate in l]
    ax2.legend(handles=h, labels=l, title = 'Isolate (Genotype)', loc = 'center left', bbox_to_anchor = (1, 0.5))
    plt.title(f'{fly} Drosophila infected with OD {od}\n{date_to_exp_dict[str_date]}')
    if saveme:
        fig.savefig(f'{drosophila_phenotype_path}/proc/survival_curves/plots/{str_date}_P07_Ehorm_drosophila_{fly}_survival_curves_od{od}_per_isolate{suff_label}.pdf', 
                    bbox_inches = 'tight')
    
########################
## PLOT Survivals across experiments and stratify by experiment
########################
for fly, od in survival_curves[['fly', 'od']].drop_duplicates().values:
    print(f'###\nFly genotype: {fly}; OD {od}')

    ## subset data to fly, od and date
    cfu_tps_fly_od = subset_df_to_group(cfu_tps, fly=fly, OD=od)
    cfu_tps_fly_od = cfu_tps_fly_od['timepoint'].unique()

    sub_survival_df = subset_df_to_group(survival_curves, fly=fly, od=od)
    sub_survival_df['isolate_rep'] = sub_survival_df['isolate'] + '(rep. ' + sub_survival_df['replicate'].astype(str) + ')'
    
    sub_death_table = subset_df_to_group(death_table_df, fly=fly, od=od)
    
    ## calculate stats based on cox proportional hazards model 
    sub_death_table['start'] = sub_death_table['start'].astype(str).str.split(' ').str[0]
    sub_death_table['experiment'] = sub_death_table['fly'].astype(str) + ' OD' + sub_death_table['od'].astype(str) + ' ' + sub_death_table['start'].astype(str)
    sub_death_table['genotype_cat'] = sub_death_table['genotype'].astype('category').cat.codes  
    sub_death_table['isolate_cat'] = sub_death_table['isolate'].astype('category').cat.codes  
    sub_death_table['exp_cat'] = sub_death_table['experiment'].astype('category').cat.codes  
    
    gt_to_gt_cat_dict = {cat: gt for gt, cat in sub_death_table[['genotype', 'genotype_cat']].drop_duplicates().values}
    exp_to_exp_cat_dict = {cat: exp for exp, cat in sub_death_table[['experiment', 'exp_cat']].drop_duplicates().values}

    # fit cox proportional hazards model with stratification over experiments
    cph = CoxPHFitter()
    cph.fit(sub_death_table[['time', 'event', 'genotype_cat', 'exp_cat']], duration_col='time', event_col='event')
    print(f'For final CoxPH fit the following number of flies have been used: {sub_death_table.groupby(['genotype'], as_index = False).size().to_string(index = False)}')
    
    # Display results
    cph.print_summary()
    if saveme:
        cph.summary.to_csv(f'{drosophila_phenotype_path}/proc/survival_curves/stats/P07_Ehorm_drosophila_{fly}_cph_summary_od{od}_per_genotype{suff_label}.csv')
    p_val = cph.summary.loc["genotype_cat", "p"]
    
    ## calculate median flies which survived to set the annotation in the plot right
    num_exp = len(sub_survival_df['start'].unique())
    max_time = sub_survival_df['time_hours'].max()
    is_last_tp = sub_survival_df['time_hours'] == max_time
    median_rate_surv_flies = sub_survival_df[is_last_tp].groupby('genotype').apply(lambda x: 1- (x['num_flies_died_exp'] / x['num_flies_recov']).median())
    upper_val = max(median_rate_surv_flies)+0.02
    lower_val = min(median_rate_surv_flies)-0.02
    stat_annotaton = np.mean([upper_val, lower_val])+0.005
    if stat_annotaton == 0:
        stat_annotaton += 0.025
        upper_val += 0.025
        lower_val += 0.025
        
    ## select plotting boundaries 
    max_time_plot = max_time
    if (fly, od) in max_time_dict.keys():
        max_time_plot = max_time_dict[(fly, od)]
    xdatarange_low = (-0.2, 0.5)
    xdatarange_high = (7.5, max_time_plot*1.065)
    wspace = 0.05
    ydatarange = (-3, 103)
    if fly == 'DrosDel w1118 iso':
        xdatarange_low = (-0.2, 8)
        xdatarange_high = (22, max_time_plot*1.065)
    widths_ratio = [max-min for min, max in [xdatarange_low, xdatarange_high]]

    ############
    ## survivals predicted by CoxPH grouped per genotype incl stat annotation
    ## predict survivals with integrated cox function with the replicates not as separate covariates
    
    ## predict survivals for every condition and plot
    pred_survival_df = sub_death_table[['genotype', 'genotype_cat']].drop_duplicates().copy().reset_index(drop = True)
    ## add the experiment categorical and real values from the central values
    pred_survival_df['exp_cat'] = cph._central_values.reset_index()['exp_cat'].values[0] ## just one value!
    pred_survival_df['experiment'] = pred_survival_df['exp_cat'].map(exp_to_exp_cat_dict)

    # generate survival curve estimation
    pred_survival_curve = cph.predict_survival_function(pred_survival_df[['genotype_cat', 'exp_cat']])
    ## add start of assay & add survival probability to 1
    if 0 not in pred_survival_curve.index:
        zero_time = pd.DataFrame({col: [1] for col in pred_survival_curve.columns}, index=[0]) 
        pred_survival_curve = pd.concat([zero_time, pred_survival_curve])

    # plot the estimated survival function
    fig, (ax1, ax2) = plt.subplots(figsize = (5, 4.5), ncols = 2, 
                                    gridspec_kw = {'width_ratios': widths_ratio, 'wspace': wspace})
    for ax_idx, ax in enumerate([ax1, ax2]):
        if ax_idx == 0:
            idx_to_keep = pred_survival_curve.index.get_loc(pred_survival_curve[pred_survival_curve.index <= xdatarange_low[-1]].index[0]) + 2 ## avoids large clipped data
            plt_curve = pred_survival_curve.iloc[:idx_to_keep]
        else:
            idx_to_keep = pred_survival_curve.index.get_loc(pred_survival_curve[pred_survival_curve.index >= xdatarange_high[0]].index[0]) - 1 ## avoids large clipped data
            plt_curve = pred_survival_curve.iloc[idx_to_keep:]
        for col in pred_survival_curve.columns:
            exp = pred_survival_df.loc[col, 'experiment']
            gt = pred_survival_df.loc[col, 'genotype']
            ax.step(plt_curve.index, plt_curve[col], label = gt, color = genotype_col_dict[gt], where='post', lw = 2, alpha = 1)
    
    ## annotate statistics
    ax2.plot([max_time_plot*1.02, max_time_plot*1.02], 
            [lower_val, upper_val], color = 'k')
    ax2.text(max_time_plot*1.0275, stat_annotaton, convert_pvalue_to_asterisks(p_val),
            fontsize=10, verticalalignment='center', horizontalalignment='left', rotation = 90)
    modify_axes(ax1, ax2, xdatarange_low, xdatarange_high, widths_ratio)
    ax2.legend(title = 'Genotype', loc = 'center left', bbox_to_anchor = (1, 0.5))
    plt.title(f'{fly} Drosophila infected with OD {od}\n({num_exp} experiments combined; p = {round(p_val, 4)})')
    if saveme:
        fig.savefig(f'{drosophila_phenotype_path}/proc/survival_curves/plots/P07_Ehorm_drosophila_{fly}_survival_curves_od{od}_per_genotype_CoxPHpred{suff_label}.pdf', 
                    bbox_inches = 'tight')

