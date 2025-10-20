
## analyse growth rate of different isolates 

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import scipy.stats as stats
#from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from statannotations.Annotator import Annotator
from statannotations.stats.StatTest import StatTest
import itertools

plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

# Assuming your dataframe has columns 'Time' and 'OD' (Optical Density) representing time and growth measurements
# Replace 'your_dataframe.csv' with the actual file path or dataframe variable

##############################
## FUNCTIONS
##############################


# Load your dataframe
def read_and_reshape(tecan_file, skiprows='auto', experiment_id_range = [0, 10]):
    tecan_file = os.path.expanduser(tecan_file)
    
    ## if read in auto --> get index where skip row occurs first 
    if skiprows == 'auto':
        growth_curves = pd.read_excel(tecan_file)
        ## get first col name and identify row where first instance of 'Cycles / Well' appears
        first_colname = growth_curves.columns[0]
        skiprows = growth_curves[growth_curves[first_colname] == 'Cycles / Well'].index[0] + 1
    
    ## read in data
    growth_curves = pd.read_excel(tecan_file, skiprows = skiprows)

    ## create a column which indicates new wells
    growth_curves['well_split'] = growth_curves['Cycles / Well'].isna().cumsum()
    ## remove all nas and cycles well rows --> well id is then in first row per well_spli group
    growth_curves = growth_curves[~((growth_curves['Cycles / Well'] == 'Cycles / Well') | (growth_curves['Cycles / Well'].isna()))]
    ## assign well id
    growth_curves['well_id'] = growth_curves.groupby('well_split')['Cycles / Well'].transform('first')
    growth_curves['well_id'] = growth_curves['well_id'].str[0] + growth_curves['well_id'].str[1:].str.zfill(2)
    growth_curves.drop(['well_split'], axis = 1, inplace = True)
    growth_curves = growth_curves.groupby('well_id').apply(lambda x: x.iloc[1:]).reset_index(drop = True)

    ## reshape to long format and then to wide to have per time, temp, mean, stdev,... a single column
    growth_curves_long = growth_curves.melt(id_vars = ['well_id', 'Cycles / Well'])
    growth_curves_long['variable'] = growth_curves_long['variable'].str.replace('Unnamed: ', '')
    growth_curves_wide = growth_curves_long.pivot(index = ['well_id', 'variable'], columns='Cycles / Well').droplevel(0, axis=1).reset_index()
    if experiment_id_range:
        growth_curves_wide['experiment_id'] = tecan_file.split('/')[-1][experiment_id_range[0]:experiment_id_range[1]] ## experimentid gives range between which the id should be extracted from the filename
    else:
        growth_curves_wide['experiment_id'] = tecan_file.split('/')[-1].replace('.xlsx', '')
    growth_curves_wide = growth_curves_wide.sort_values(['well_id', 'Time [s]'])
    return growth_curves_wide

def exp_setup_wide(setup_file, sheet=0, pos_ctr = 'ATCC', neg_ctr = ['nan', 'neg.', 'negative', '-', 'media']):
    setup_file = os.path.expanduser(setup_file)
    exp_setup = pd.read_excel(setup_file, sheet_name = sheet, header = None)
    ## create a column which indicates a new descriptor
    exp_setup['setup_split'] = exp_setup[0].isna().cumsum()
    exp_setup.columns = ['rowdescription', 'well_ids'] + exp_setup.loc[exp_setup[0].str.lower() == 'columns'].values.tolist()[0][2:-1] + ['setup_split']
    ## drop empty rows 
    exp_setup = exp_setup[~(exp_setup['rowdescription'].isna() & exp_setup['well_ids'].isna())]
    ## assign descriptor to the full plate
    exp_setup['descriptor'] = exp_setup.groupby('setup_split')['well_ids'].transform('first')
    exp_setup.drop(['rowdescription', 'setup_split'], axis = 1, inplace = True)
    ## make df to long format and remove unecessary rows
    exp_setup_long = exp_setup.melt(id_vars = ['well_ids', 'descriptor'], var_name='col_id') 
    exp_setup_long = exp_setup_long[~((exp_setup_long['well_ids'].isna()) | (exp_setup_long['well_ids'] == exp_setup_long['descriptor']))]
    exp_setup_long['descriptor'] = exp_setup_long['descriptor'].str.lower()
    ## generate wide format to have per well one row
    exp_setup_wide = exp_setup_long.pivot(index = ['well_ids', 'col_id'], columns = 'descriptor').droplevel(0, axis = 1).reset_index()
    exp_setup_wide.rename(columns = {'sampleid/species': 'sampleid', 'well_ids': 'row_id'}, inplace = True)
    exp_setup_wide['sampleid'] = exp_setup_wide['sampleid'].astype(str)
    exp_setup_wide.loc[exp_setup_wide['sampleid'].str.contains('|'.join(pos_ctr)), 'control'] = 'positive'
    exp_setup_wide.loc[exp_setup_wide['sampleid'].str.lower().str.contains('|'.join(neg_ctr)), 'control'] = 'negative'
    exp_setup_wide.loc[exp_setup_wide['sampleid'].str.lower().str.contains('|'.join(neg_ctr)), 'sampleid'] = 'negative control'
    exp_setup_wide.rename(columns = {'well_ids' : 'row_id'}, inplace = True)
    exp_setup_wide['col_id'] = exp_setup_wide['col_id'].astype(int)
    exp_setup_wide['well_id'] = exp_setup_wide['row_id'] + exp_setup_wide['col_id'].astype(str).str.zfill(2)
    return exp_setup_wide

def annotate_contaminated_wells(df, threshold_contamination = 0.2, cutoff_contaminated_wells_stop = 3):
    ## QC
    ## identify contaminated control wells
    df = df.sort_values('Time [s]')
    df['time_h'] = (df['Time [s]'] / 3600).astype(float).round(2) ## calculate time in hours to account for small (second differences) in measurement
    df['raw_OD_diff'] = df.groupby('well_id')['Mean'].transform(lambda x: x.iat[-1] - x.iat[0]) ## calculate the difference between first and last measured timepoint
    df.loc[((df['raw_OD_diff'] > threshold_contamination) & (df['control'] == 'negative')), 'contaminated'] = True
    contaminated_neg_controls = df.loc[df['contaminated'] == True, ['well_id', 'raw_OD_diff', 'contaminated']].drop_duplicates()
    if not contaminated_neg_controls.empty:
        print(f'Contamination found in negative controls of experiment {df["experiment_id"].values[0]}!')
    if (len(contaminated_neg_controls) > cutoff_contaminated_wells_stop):
        SystemExit(f'Found too many ({len(contaminated_neg_controls)}) contaminated negative control wells! Script will stop analysis since results might be wrong!')
    return df

def remove_wells_without_growth(df, threshold_growth = 0.4, cutoff_min_samples_remaining = -1):
    df['norm_od_diff'] = df.groupby('well_id')['norm_od'].transform(lambda x: x.max() - x.min()) ## calculate the maximum OD increase
    ##sampleids with no growth 
    smpl_no_growth = df.loc[df['norm_od_diff'] <= threshold_growth, 'sampleid'].unique()
    if len(smpl_no_growth) > 0:
        print(f'The following samples have replicates in dataset which indicate no growth based on {threshold_growth}: {", ".join(smpl_no_growth)}')
    df_cleaned = df[df['norm_od_diff'] >= threshold_growth].copy() ## remove wells with no growth 
    num_replicates_smpl = df_cleaned.groupby('sampleid')['well_id'].transform('nunique')
    if cutoff_min_samples_remaining == -1: ## all replicates should remain
        replicates_per_sample = df.groupby('sampleid')['well_id'].unique().str.len().to_dict() ## get number of replicates per isolate
        df_cleaned = df_cleaned[num_replicates_smpl == df_cleaned['sampleid'].map(replicates_per_sample)]
    else: ## at least cutoff_min_samples_remaining should remain
        df_cleaned = df_cleaned[num_replicates_smpl >= cutoff_min_samples_remaining]
    return df_cleaned

def normalize_ods(df, remove_contamonated = True, remove_stdev = 0.4):
    negative_ctr_bool = (df['control'] == 'negative')
    if remove_contamonated:
        negative_ctr_bool = (negative_ctr_bool & (df['contaminated'] != True))
    if remove_stdev:
        negative_ctr_bool = (negative_ctr_bool & (df['StDev'] <= 0.4))
    neg_controls = df[negative_ctr_bool] ## select negative controls which are not contaminated
    df_samples = df[df['control'].isna()] ## select neither pos nor neg controls
    od_bg_blank = neg_controls.groupby('Time [s]')['Mean'].mean().reset_index() ## calculate balnk value for each timepoint
    od_bg_blank.rename(columns = {'Mean': 'background_od'}, inplace = True) 
    df_samples = pd.merge(df_samples, od_bg_blank, how = 'left', on = 'Time [s]')
    df_samples['norm_od'] = df_samples['Mean'] - df_samples['background_od'] ## remove blanks from each timpoint
    return df_samples

def combine_and_normalize(growth_curve_file, setupfile, stop_if_contaminated = False, stop_if_no_growth = False, threshold_contamination = 0.2, cutoff_contaminated_wells_stop = 3, remove_stdev = 0.4, threshold_growth = 0.4, cutoff_min_samples_remaining = -1):
    df = read_and_reshape(growth_curve_file)
    exp_setup = exp_setup_wide(setupfile)
    df = pd.merge(df, exp_setup, how = 'left', on = 'well_id')
    df = df[~df['sampleid'].isna()] ## check that all wells have a sample id

    ## QC 
    df = annotate_contaminated_wells(df, threshold_contamination, cutoff_contaminated_wells_stop)
    if stop_if_contaminated and any(df['contaminated'] == True):
        print(f'\nNOTE: The calculation was force stopped for {df["experiment_id"].unique()[0]} due to contamination.\n')
        return pd.DataFrame()
    ## normalize data NOTE: normalization is important to account for media which might alter the OD over time!
    df_norm = normalize_ods(df, remove_stdev = remove_stdev)
    
    ## after normalization check if in wells growth was detected!
    df_norm_cleaned = remove_wells_without_growth(df_norm, threshold_growth, cutoff_min_samples_remaining)
    if stop_if_no_growth and len(df_norm_cleaned) < len(df_norm):
        print(f'\nNOTE: The calculation was force stopped for {df["experiment_id"].unique()[0]} due to limited growth in some samples!\n')
        return pd.DataFrame()
    return df_norm_cleaned

def analyse_growth_sliding_window(df, timecol='Time [s]', group_col=False, od_col='Mean', window_size=3):
    # Calculate the time difference between consecutive timepoints
    df[timecol] = df[timecol].astype(float)
    df = df.sort_values(timecol).reset_index(drop = True)
    if group_col:
        df['smoothed_od'] = df.groupby(group_col)[od_col].rolling(window=window_size, center=True).mean().reset_index(0,drop=True)
        df[['time_diff_sliding', 'growth_diff_sliding']] = df.groupby(group_col)[[timecol, 'smoothed_od']].diff()
        # Calculate the growth rate
        df['growth_rate_sliding'] = df['growth_diff_sliding'] / df['time_diff_sliding']
    else:
        df['smoothed_od'] = df[od_col].rolling(window=window_size, center=True).mean().reset_index(0,drop=True)
        df[['time_diff_sliding', 'growth_diff_sliding']] = df[[timecol, 'smoothed_od']].diff()
        # Calculate the growth rate
        df['growth_rate_sliding'] = df['growth_diff_sliding'] / df['time_diff_sliding']
    return df

##############################
## VARIABLES
##############################

## dict contains rawdatapath as keys with dict as values of result file (platereader) and its setup file as key and value respectively 
platereader_setupfile_dict = {
    '~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/growth_curves/raw/': {
    '2024_12_13_Ehorm_wt_fimZ_growthcurve_LB.xlsx': '2024_12_13_Ehorm_wt_fimZ_growthcurve_LB_setup.xlsx',
    '2024_12_14_Ehorm_wt_fimZ_growthcurve_LB.xlsx': '2024_12_14_Ehorm_wt_fimZ_growthcurve_LB_setup.xlsx',
    '2024_12_15_Ehorm_wt_fimZ_growthcurve_LB.xlsx': '2024_12_15_Ehorm_wt_fimZ_growthcurve_LB_setup.xlsx'
    }
}

## set parameters for troubleshooting
growth_rate_threshold = 0.1 
clean_data = True
reduce_to_replicate_num = 3
genotype_col_dict = {'wt': '#f1f0f0', 
                    'FimZ[F126L]': '#f9af92'}

##############################
## MAIN SCRIPT
##############################


growth_rate_df_all = pd.DataFrame()
growth_curve_df_all = pd.DataFrame()
for rawdata_path, rawfilenames in platereader_setupfile_dict.items():
    for growth_curve_file, setupfile in rawfilenames.items():
        growthrate_results = combine_and_normalize(rawdata_path + growth_curve_file, rawdata_path + setupfile, stop_if_contaminated=True, stop_if_no_growth=True)
        if growthrate_results.empty:
            continue
        growthrate_results['sampleid'] = growthrate_results['sampleid'].str.replace(' ', '_')
        df_growthrate_simple_sliding = analyse_growth_sliding_window(growthrate_results, timecol='Time [s]', group_col='well_id', od_col='norm_od', window_size = 3)

        ## extract maximal od 
        max_od = df_growthrate_simple_sliding.groupby('well_id')['norm_od'].agg('max').to_dict()
        max_od_diff = df_growthrate_simple_sliding.groupby('well_id')['norm_od'].agg(lambda x: max(x) - min(x)).to_dict()

        ## extract max growth rates
        growth_rate_df = df_growthrate_simple_sliding.sort_values('growth_rate_sliding', ascending = False).groupby('well_id').head(1)
        growth_rate_df['od_p_h'] = growth_rate_df['growth_rate_sliding'] * 3600 ## normalize to get to od change per hour
        growth_rate_df['method'] = 'sliding'
        growth_rate_df['max_od'] = df_growthrate_simple_sliding['well_id'].map(max_od)
        growth_rate_df['max_od_diff'] = df_growthrate_simple_sliding['well_id'].map(max_od_diff)

        ## normalize the time around the maximal grouth rate
        df_growthrate_simple_sliding = pd.merge(df_growthrate_simple_sliding, growth_rate_df[['well_id', 'Time [s]']], how ='left', on = 'well_id', suffixes=['', '_max_growth'])
        df_growthrate_simple_sliding['norm_time_max_growth'] = df_growthrate_simple_sliding['Time [s]'] - df_growthrate_simple_sliding['Time [s]_max_growth']
        df_growthrate_simple_sliding = df_growthrate_simple_sliding.drop('Time [s]_max_growth', axis = 1)

        growth_rate_df_all = pd.concat([growth_rate_df_all, growth_rate_df])
        growth_curve_df_all = pd.concat([growth_curve_df_all, df_growthrate_simple_sliding[['experiment_id', 'well_id', 'sampleid', 'Time [s]', 'norm_time_max_growth', 'norm_od', 'growth_rate_sliding']]])

growth_rate_df_all = growth_rate_df_all.reset_index(drop = True)
growth_curve_df_all = growth_curve_df_all.reset_index(drop = True)

translate_sampleids = { 'P07_2048': 'wt',
                        'P07_2055': 'wt',
                        'P07_2052': 'FimZ[F126L]',
                        'P07_2237': 'FimZ[F126L]',
                        'L00071_01_1': 'FimZ[F126L]',
                        'L00071_02_1': 'FimZ[F126L]'}
growth_rate_df_all['border'] = np.where(growth_rate_df_all['well_id'].str.contains('01') | growth_rate_df_all['well_id'].str.contains('12') | growth_rate_df_all['well_id'].str.contains('A') | growth_rate_df_all['well_id'].str.contains('H'), 'border', 'inner')
growth_rate_df_all['genotype'] = growth_rate_df_all['sampleid'].map(translate_sampleids)
growth_rate_df_all['sampleid_ext'] = growth_rate_df_all['sampleid'] + ' (' + growth_rate_df_all['genotype'] + ')'
growth_rate_df_all['media'] = 'LB'


growth_rate_df_all['sampleid'] = pd.Categorical(growth_rate_df_all['sampleid'], translate_sampleids.keys(), ordered=True)
growth_rate_df_all['sampleid_ext'] = pd.Categorical(growth_rate_df_all['sampleid_ext'], [f'{key} ({val})' for key, val in translate_sampleids.items()], ordered=True)
growth_rate_df_all['genotype'] = pd.Categorical(growth_rate_df_all['genotype'], ['wt', 'FimZ[F126L]'], ordered=True)
growth_curve_df_all['sampleid_ext'] = growth_curve_df_all['sampleid'] + ' (' + growth_curve_df_all['sampleid'].map(translate_sampleids) + ')'
growth_curve_df_all['sampleid_ext'] = pd.Categorical(growth_curve_df_all['sampleid_ext'], [f'{key} ({val})' for key, val in translate_sampleids.items()], ordered=True)
unique_sampleids = growth_curve_df_all.sort_values('sampleid_ext')['sampleid_ext'].unique()
cmap = plt.get_cmap('tab20')
color_dict = {sampleid: cmap(i % 20) for i, sampleid in enumerate(growth_curve_df_all['sampleid_ext'].unique())}


########################################
# Combine genotypes and add stats to plots
######################################## 
## The plots below are used for the MS

## generate the test function (import from stats.ranksums)
mannwhitneyu_long_name = 'mannwhitneyu test'
mannwhitneyu_short_name = ''
mannwhitneyu_func = stats.mannwhitneyu
mannwhitneyu_test = StatTest(mannwhitneyu_func, mannwhitneyu_long_name, mannwhitneyu_short_name)

genptype_pairs = list(itertools.combinations(growth_rate_df_all['genotype'].unique(), 2))

## plot the max growth rates 
for pval_method in ['star', 'text']:
    fig, ax = plt.subplots(figsize = (3.5, 4))
    sns.boxplot(data = growth_rate_df_all, x = 'genotype', y = 'od_p_h', hue = 'genotype', palette = genotype_col_dict, dodge = False, fliersize = 0, legend = True, ax = ax, zorder = 5)
    sns.stripplot(data = growth_rate_df_all, x = 'genotype', y = 'od_p_h', hue = 'genotype', palette = genotype_col_dict, dodge = False, alpha = 0.8, jitter = 0.25, s = 8, edgecolor = 'k', linewidth = 0.4, legend = False, ax = ax, zorder = 10)
    annot = Annotator(ax, genptype_pairs, data=growth_rate_df_all, x='genotype', y='od_p_h')
    if pval_method == 'star': 
        annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=False, line_height=0.01, text_offset=-0.5, pvalue_thresholds=[[1e-4, '****'], [1e-3, '***'], [1e-02, '**'], [0.05, '*'], [1, 'ns']]).apply_test()
    else:
        annot.configure(test=mannwhitneyu_test, verbose=2, hide_non_significant=False, line_height=0.01, text_offset=-0.5, text_format='full', pvalue_format_string='{:.1e}').apply_test()
    annot.annotate(line_offset_to_group=0.2, line_offset=-2)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.04))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.02))
    ax.grid(which='major', axis = 'y', linestyle='-', linewidth=0.4, color='black', zorder = 1)
    ax.grid(which='minor', axis = 'y', linestyle=':', linewidth=0.4, color='black', zorder = 1)
    ax.set_ylim((0, ax.get_ylim()[1]*1.08))
    h, l = ax.get_legend_handles_labels() ## get legend handles and labels
    ax.legend(title = 'Genotype',loc = 'center left', bbox_to_anchor = (1, 0.5))
    ax.set_ylabel('max(Delta[OD595]/h)')
    plt.xlabel(None)
    plt.tight_layout()
    fig.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/growth_curves/proc/2025_01_08_P07_Ehorm_wt_fimZ_max_growthrate_per_genotype_pval{pval_method}.pdf'))

## growth rate distribution
## plot growth rate vs time (non-smoothed but sliding window)
isolate_col_dict = {'P07_2048 (wt)': '#dddddd', 
                   'P07_2055 (wt)': '#ababab',
                   'P07_2052 (FimZ[F126L])': '#FAC4AF', ## 30% saturation of #f9af92
                   'P07_2237 (FimZ[F126L])': '#FAA07D',             ## 50% saturation of #f9af92
                   'L00071_01_1 (FimZ[F126L])': '#FA8557', ## 65% saturation of #f9af92
                   'L00071_02_1 (FimZ[F126L])': '#FA6225' ## 85% saturation of #f9af92
                   }

result_df_tmp = growth_curve_df_all.copy()
result_df_tmp['norm_time_max_growth'] = np.trunc(result_df_tmp['norm_time_max_growth']/10)*10 ## round the seconds closer to 0 (negative values np.ceil, positive np.floor) to keep the 20 min for proper estimation of mean values
result_df_tmp = result_df_tmp[~result_df_tmp['growth_rate_sliding'].isna()] ## remove nan values from sliding window (first and last 2 datapoints)
result_df_tmp = result_df_tmp.groupby(['sampleid_ext', 'norm_time_max_growth']).filter(lambda x: len(x) > 1) ## require that at least 2 datapoints are present per time to allow plotting (to remove outliers)
result_df_tmp['norm_time_max_growth'] = result_df_tmp['norm_time_max_growth']/3600 ## convert to hours
result_df_tmp['growth_rate_sliding'] = result_df_tmp['growth_rate_sliding']*3600 ## convert to hours
result_df_tmp['genotype'] = result_df_tmp['sampleid'].map(translate_sampleids)

fig, ax = plt.subplots(figsize = (6, 4))
sns.lineplot(data = result_df_tmp, x = 'norm_time_max_growth', y = 'growth_rate_sliding', hue = 'sampleid_ext', palette = isolate_col_dict, estimator = 'mean', errorbar=('ci', 95), alpha = 0.6, linewidth = 2, ax = ax)
ax.xaxis.set_major_locator(ticker.MultipleLocator(4))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.set_ylabel('Delta[OD595]/h')
ax.set_xlabel('Time to maximum growth rate [h]')
ax.legend(title = 'Isolate (Genotype)',loc = 'center left', bbox_to_anchor = (1, 0.5))
plt.tight_layout()
fig.savefig(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/growth_curves/proc/2025_01_08_P07_Ehorm_wt_fimZ_growthrate_sliding_curve.pdf'))

