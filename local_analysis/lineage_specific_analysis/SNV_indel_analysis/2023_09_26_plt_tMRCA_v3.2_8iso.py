
## I have the underlying data of the mutations
## I have the mutation rate with the cis
## I can bootstrap the mutations at each locus with replacement (and one time just mean)
## and can calculate the underlying confidence interval and compare results!

## NOTE: This is an adaption of the version from 2023_09_26!
## With this script confidence intervals in a positive range can be evaluated only!


import os
import glob
import numpy as np
from scipy import stats
from scipy.optimize import minimize
import pandas as pd
import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

sys.path.insert(0, os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/analysis_scripts/modules/'))

import analysispy_module as apy

### functions

def read_mol_clock_data_file(file):
    with open(file, 'r') as fid:
        mol_clock_data = np.loadtxt(fid, skiprows = 1, delimiter = ',')
    if len(np.shape(mol_clock_data)) == 1: ## if 1d array --> reshape to 2d 
        mol_clock_data = np.reshape(mol_clock_data, (-1, 4))
    return mol_clock_data

def calculate_molecular_clock_rate(mol_clock_data, subj_spec, unique_tps):

    num_timepoints = len(unique_tps)

    print(f'\n\nAnalysing {subj_spec}')
    ## calculate patient specific mutation rate
    if num_timepoints > 1:
        num_samples = len(mol_clock_data)
        tps = mol_clock_data[:,1]
        num_mut_scaled = mol_clock_data[:,2]/mol_clock_data[:,3] * basescale

        regress_dict = apy.calc_regression(np.column_stack( (tps, num_mut_scaled) ),subj_spec,num_samples, num_timepoints)

        mut_rate_obs = regress_dict['slope'] * 365 ## mutation rate based on regression over all isolates 
        ts = abs(stats.t.ppf( (1-confidence) / 2, len(tps)-2)) ## (confidence interval and the degrees of freedom minus num of parameters to estimate)
        SE_mut_rate_obs = regress_dict['std_err'] * 365 
        margin_of_error = SE_mut_rate_obs * ts ## moe = std_err * t-score
        
        ci_lower = mut_rate_obs - margin_of_error
        ci_upper = mut_rate_obs + margin_of_error
        mut_rate_obs_ci95 = [ci_lower, ci_upper]
        
    else:
        print(f'{subj_spec} has only 1 tp --> skipping')
        mut_rate_obs = np.nan
        mut_rate_obs_ci95 = [np.nan, np.nan]
        num_samples = np.nan
    return [mut_rate_obs, mut_rate_obs_ci95, num_samples]

def calculate_tmrca_via_MCsimulation(num_muts_per_smpl, mut_rate, ci95_vals, n_sim=10000, seed=123, mut_rate_fit_method = 'gamma', plot_histogram = False, subj_spec = ''):
    ## calculate the observed point estimate for the tmrca as well as 95% CI values given Monte Carlo simulations of the mutation rate and 
    ## bootstraps of the number of observed SNVs
    ## To fit the distribution for the mutation rate, choose either a 'normal' (0-bounded), 'log-normal' or 'gamma' distribution
    ## plot histograms to evaluate model fits

    if mut_rate != mut_rate:
        print('No mutation rate given. Skipping calculation')
        return [np.nan, [np.nan, np.nan]]
    
    ## calculate patient specific tMRCA with oberved mutation rate of species
    mean_mut = np.mean(num_muts_per_smpl)
    tmrca_mutrate = mean_mut / mut_rate ## calculate observed point estimate; should be close to `np.mean(bootstrapped_tmrca)`!

    ## return just tmrca if no CI for mutationrate are given
    if ci95_vals == [0, 0]:
        print('WARNING! No confidence intervals were given')
        return [tmrca_mutrate, [np.nan, np.nan]]

    ci_diff_to_meanpoint = (ci95_vals[1]-ci95_vals[0])/2
    ci_diff_to_mean_lower = 1-(mut_rate / (ci95_vals[0]+ci_diff_to_meanpoint))
    ci_diff_to_mean_upper = 1-(mut_rate / (ci95_vals[1]-ci_diff_to_meanpoint))
    print(f'The CI value is skewed to left(negative) or right(posiive) by {ci_diff_to_mean_lower} and {ci_diff_to_mean_upper} ({ci95_vals[0]}, {mut_rate}, {ci95_vals[1]})')
    
    if (mut_rate_fit_method == 'normal') or plot_histogram:
        np.random.seed(seed)
        # fit truncated (0-bounded) normal distribution to sample from n_sim times
        sd_mut_rate = (ci95_vals[1] - ci95_vals[0]) / (2 * 1.96)
        mut_rate_dist_sampled_normal = np.random.normal(loc=mut_rate, scale=sd_mut_rate, size=n_sim)
        mut_rate_dist_sampled_truncnormal = stats.truncnorm((0 - mut_rate) / sd_mut_rate, (np.inf - mut_rate) / sd_mut_rate, loc=mut_rate, scale=sd_mut_rate).rvs(n_sim)
        if mut_rate_fit_method == 'normal':
            mut_rate_dist_sampled = mut_rate_dist_sampled_truncnormal

    if (mut_rate_fit_method == 'log-normal') or plot_histogram:
        np.random.seed(seed)
        # fit distribution to sample from n_sim times (--> mutation-rates are 0-bounded! to avoid negative values, use log-normal distribution)
        log_mut_rate_ci95_l, log_mut_rate, log_mut_rate_ci95_u = np.log([ci95_vals[0], mut_rate, ci95_vals[1]])
        sd_log_mut_rate_log = (log_mut_rate_ci95_u - log_mut_rate_ci95_l) / (2 * 1.96) ## get standard deviation to sample from distribution
        mut_rate_dist_sampled_lognormal = np.random.lognormal(mean=log_mut_rate, sigma=sd_log_mut_rate_log, size=n_sim) ## sample from mutation rate distribution n_sim times 
        if mut_rate_fit_method == 'log-normal':
            mut_rate_dist_sampled = mut_rate_dist_sampled_lognormal

    if (mut_rate_fit_method == 'gamma') or plot_histogram:
        np.random.seed(seed)
        # fit distribution to sample from n_sim times (--> mutation-rates are 0-bounded! to avoid negative values, use gamma distribution)
        ## NOTE: While this fits the CI values truely, it also performs well on credible intervals given those parameters. HOWEVER: This still does not mean that it is the true shape, but it is still rather the better guess than just using the HDP interval estimates without error propagation!!!
        gamma_shape, gamma_scale = fit_gamma_from_ci(mut_rate, ci95_vals[0], ci95_vals[1], mean_weight=2, ci_lower_weight=1, ci_upper_weight=1) ## fit gamma distribution given the mean and CI values via minimization of the difference between guessed and observed values; apply weights if one value is more important
        mut_rate_dist_sampled_gamma_fit = stats.gamma.rvs(a=gamma_shape, scale=gamma_scale, size=n_sim) ## sample from fitted gamma distribution
        if mut_rate_fit_method == 'gamma':
            mut_rate_dist_sampled = mut_rate_dist_sampled_gamma_fit
        
    if plot_histogram:
        trunc_norm_y, trunc_norm_bin_edges = np.histogram(mut_rate_dist_sampled_normal, bins=100)
        trunc_norm_x = 0.5 * (trunc_norm_bin_edges[1:] + trunc_norm_bin_edges[:-1])
        log_normal_y, log_normal_bin_edges = np.histogram(mut_rate_dist_sampled_lognormal, bins=100)
        log_normal_x = 0.5 * (log_normal_bin_edges[1:] + log_normal_bin_edges[:-1])
        gamma_fit_y, gamma_fit_bin_edges = np.histogram(mut_rate_dist_sampled_gamma_fit, bins=100)
        gamma_fit_x = 0.5 * (gamma_fit_bin_edges[1:] + gamma_fit_bin_edges[:-1])

        fig, ax = plt.subplots(figsize = (9,3.5))
        ## plot actual values
        ax.axvline(x = mut_rate, lw = 3, color = 'k', ls = '-', label = 'Observed')
        [ax.axvline(x = ci, lw = 3, color = 'k', ls = '--') for ci in ci95_vals]
        ## plot truncated normal distribution
        [ax.axvline(x = ci, lw = 2, color = '#1b9e77', ls = '--', zorder = 5) for ci in np.percentile(mut_rate_dist_sampled_truncnormal, [2.5, 97.5])]
        ax.axvline(x = np.mean(mut_rate_dist_sampled_truncnormal), lw = 2, color = '#1b9e77', ls = '-')
        ax.plot(trunc_norm_x, trunc_norm_y, color='#1b9e77', lw = 1.5, label = '0-bounded\nnormal')
        ## plot log-normal distribution
        [ax.axvline(x = ci, lw = 2, color = '#7570b3', ls = '--', zorder = 10) for ci in np.percentile(mut_rate_dist_sampled_lognormal, [2.5, 97.5])]
        ax.axvline(x = np.mean(mut_rate_dist_sampled_lognormal), lw = 2, color = '#7570b3', ls = '-')
        ax.plot(log_normal_x, log_normal_y, color='#7570b3', lw = 1.5, label = 'log-normal')
        ## plot gamma distribution
        [ax.axvline(x = ci, lw = 2, color = '#d95f02', ls = '--', zorder = 15) for ci in np.percentile(mut_rate_dist_sampled_gamma_fit, [2.5, 97.5])]
        ax.axvline(x = np.mean(mut_rate_dist_sampled_gamma_fit), lw = 2, color = '#d95f02', ls = '-')
        ax.plot(gamma_fit_x, gamma_fit_y, color='#d95f02', lw = 1.5, label = 'gamma')
        ax.set_ylim((0, None))
        ax.legend(title='Distributions with\nmean (solid lines) and\nCI95% (dashed lines)', loc = 'center left', bbox_to_anchor = (1, 0.5))
        ax.set_title(subj_spec)
        plt.tight_layout()
        fig.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2025_06_09_{subj_spec}_mol_clock_distribution_fit_for_MCsimulation.svg'))

    ## bootstrap n_sim times the number of mutations and get mean number of mutation 
    bootstrap_idx = np.random.randint(0, len(num_muts_per_smpl), size=(n_sim, len(num_muts_per_smpl)))
    bootstrapped_mean_mut = np.mean(num_muts_per_smpl[bootstrap_idx], axis = 1)
    bootstrapped_tmrca = bootstrapped_mean_mut / mut_rate_dist_sampled 
    tmrca_mutrate_ci95 = np.percentile(bootstrapped_tmrca, [2.5, 97.5])

    ## calculate point estimate 
    print(f'The difference between the observed tMRCA point estimate and the MC simulated point estimate is: {np.mean(bootstrapped_tmrca) - tmrca_mutrate}')
    print(f'The difference between the observed tMRCA point estimate and the median of the MC simulations is: {np.median(bootstrapped_tmrca) - tmrca_mutrate}')

    return [tmrca_mutrate, tmrca_mutrate_ci95]

def objective_function(gamma_args_tuple, mut_rate, ci_lower, ci_upper, mean_weight=1, ci_lower_weight=1, ci_upper_weight=1):
    ## calculate difference between guessed and observed mutation rate (mean and CIs) and 
    ## sum up (weighted) squared errors to minimize upon 

    # validate that agruments are ≥0
    if gamma_args_tuple[0] <= 0 or gamma_args_tuple[1] <= 0:
        return np.inf

    # get mean of gamma distribution
    gamma_mean = gamma_args_tuple[0] * gamma_args_tuple[1]

    # get CI of gamma distribution
    gamma_lower_quantile = stats.gamma.ppf(0.025, a=gamma_args_tuple[0], scale=gamma_args_tuple[1])
    gamma_upper_quantile = stats.gamma.ppf(0.975, a=gamma_args_tuple[0], scale=gamma_args_tuple[1])

    ## get errors and return sum of errors with weights (optional) to apply minimization
    error_mean = (gamma_mean - mut_rate)**2
    error_lower_ci = (gamma_lower_quantile - ci_lower)**2
    error_upper_ci = (gamma_upper_quantile - ci_upper)**2

    return error_mean*mean_weight + error_lower_ci*ci_lower_weight + error_upper_ci*ci_upper_weight

def fit_gamma_from_ci(mut_rate, ci_lower, ci_upper, mean_weight=1, ci_lower_weight=1, ci_upper_weight=1):
    ## fitting a gamma distribution to a mutation rate given the mean and 95% confidence intervals
    if ci_upper - ci_lower == 0:
        return None, None
    else:
        sd_mut_rate = (ci_upper - ci_lower) / (2 * 1.96) ## zsocre of 1.96 corresponds to 95% confidence interval; Note it is a good start but gamma will differ!
        var_mut_rate = sd_mut_rate**2

    ini_scale = var_mut_rate / mut_rate ## https://library.virginia.edu/data/articles/getting-started-with-gamma-regression
    ini_shape = mut_rate**2 / var_mut_rate ## https://library.virginia.edu/data/articles/getting-started-with-gamma-regression
        
    # Bounds for the parameters (shape and scale must be positive)
    bounds = ((1e-16, None), (1e-16, None)) # Small positive lower bound

    
    minimizer_result = minimize(objective_function,
                                (ini_shape, ini_scale),
                                args=(mut_rate, ci_lower, ci_upper, mean_weight, ci_lower_weight, ci_upper_weight),
                                bounds=bounds,
                                method='L-BFGS-B'  # default method using bounds
                                )

    if minimizer_result.success:
        return minimizer_result.x
    else:
        print(f'Optimization failed: {minimizer_result.message}')
        return None, None

def get_tmrca_per_tp(mol_clock_data, unique_tps, mut_rate_obs = np.nan, mut_rate_obs_ci95 = np.nan, mut_rate_literature = np.nan, mut_rate_ci95_literature = np.nan, basescale = 1000000, req_isolate_num = 5, subj_spec = '', plot_histogram=False):
    
    mut_rate_dict = {}
    ## loop over all timepoints 
    ploted_hist = False
    for tp in unique_tps:
        tp_idx = np.where(mol_clock_data[:, 1] == tp)[0]
        ## check if tp has at least num required isolates to calculate tMRCA
        if len(tp_idx) >= req_isolate_num:
            
            ## calculate scaled number of mutations
            num_mut_scaled = mol_clock_data[tp_idx,2]/mol_clock_data[tp_idx,3] * basescale

            ## plot histogram only once!
            if ploted_hist:
                plot_histogram = False
            ploted_hist = True
            [tmrca_pat_spec_mutrate, tmrca_pat_spec_mutrate_ci95] = calculate_tmrca_via_MCsimulation(num_mut_scaled * 365, mut_rate_obs, mut_rate_obs_ci95, mut_rate_fit_method = 'gamma', n_sim=10000, seed=123, plot_histogram=plot_histogram, subj_spec = subj_spec) ## note: *365 is to convert mutation rate from a yearly basis to a daily! By multiplying num_mut_scaled it is easier than implementing it to the list of CIs as well!
            [tmrca_lit_mutrate, tmrca_lit_mutrate_ci95] = calculate_tmrca_via_MCsimulation(num_mut_scaled * 365, mut_rate_literature, mut_rate_ci95_literature, mut_rate_fit_method = 'gamma', n_sim=10000, seed=123, plot_histogram=plot_histogram, subj_spec = subj_spec) ## note: *365 is to convert mutation rate from a yearly basis to a daily! By multiplying num_mut_scaled it is easier than implementing it to the list of CIs as well!
            
            ## get how many days since hosp admission to tp:
            d_since_hospadmission = np.unique(mol_clock_data[tp_idx, 0])[0]

            mut_rate_dict[tp + 1] = {'days_in_hosp': d_since_hospadmission,
                                    'days_to_root': np.unique(mol_clock_data[tp_idx,1])[0],
                                    'num_isolates': len(tp_idx),
                                    'mut_rate_obs': mut_rate_obs,
                                    'mut_rate_obs_ci95': mut_rate_obs_ci95,
                                    'mut_rate_lit': mut_rate_literature,
                                    'mut_rate_lit_ci95': mut_rate_ci95_literature,
                                    'mean_mut': np.mean(num_mut_scaled/basescale), 
                                    'ci95_mut': np.percentile(num_mut_scaled/basescale, [2.5, 97.5]), 
                                    'num_muts': f"[{','.join([str(mut) for mut in num_mut_scaled])}]",
                                    'tmrca_meanmut_pat_spec': tmrca_pat_spec_mutrate,
                                    'ci95_meanmut_tmrca_pat_spec': tmrca_pat_spec_mutrate_ci95,
                                    'tmrca_meanmut_lit': tmrca_lit_mutrate,
                                    'ci95_meanmut_tmrca_lit': tmrca_lit_mutrate_ci95
                                    }
    ## convert data to df 
    tmp_mutation_rate_df = pd.DataFrame(mut_rate_dict).T.reset_index()
    tmp_mutation_rate_df.rename(columns={'index':'timepoint'}, inplace = True)
    tmp_mutation_rate_df['subj_spec'] = subj_spec

    return tmp_mutation_rate_df

def add_metadata_to_df(tmrca_spec, inf_day_dict):

    tmrca_spec['inf_date'] = tmrca_spec['subj_spec'].map(inf_day_dict)
    tmrca_spec['subj_spec'] = tmrca_spec['subj_spec'].str.split('_').str[1:3].str.join('_')

    tmrca_spec['patient'] = tmrca_spec['subj_spec'].str.split('_').str[0]
    tmrca_spec['species'] = tmrca_spec['subj_spec'].str.split('_').str[1].str.split('-').str[0]
    tmrca_spec['species'] = tmrca_spec['species'].str[0] + '. ' + tmrca_spec['species'].str[1:]
    tmrca_spec['clade'] = tmrca_spec['subj_spec'].str.split('_').str[1].str.split('-').str[-1].str.replace('c', '').astype(int)
    tmrca_spec['spec_clade_subj_plt'] = tmrca_spec['species'] + ' clade' + tmrca_spec['clade'].astype(str) + ' (' + tmrca_spec['patient'] + ')'
    tmrca_spec['spec_subj_plt'] = tmrca_spec['species'] + ' (' + tmrca_spec['patient'] + ')'

    tmrca_spec = tmrca_spec.sort_values(['patient', 'species', 'clade'])
    return tmrca_spec

def combine_literature_inferred_clock_for_plot(tmrca_spec_df, clade_to_use_inferred_clocks = [], annotate_clock = {}, color_of_clock = {}):

    ## bring the respective plots to one column for plot
    tmrca_spec_df['tmrca_fin'] = tmrca_spec_df['tmrca_meanmut_lit']
    tmrca_spec_df['ci95_tmrca_fin'] = tmrca_spec_df['ci95_meanmut_tmrca_lit']
    tmrca_spec_df['mut_rate_fin'] = tmrca_spec_df['mut_rate_lit']
    tmrca_spec_df['ci95_mut_rate_fin'] = tmrca_spec_df['mut_rate_lit_ci95']
    is_ehorm_bool = (tmrca_spec_df['spec_clade_subj_plt'].isin(clade_to_use_inferred_clocks))
    tmrca_spec_df.loc[is_ehorm_bool, 'tmrca_fin'] = tmrca_spec_df.loc[is_ehorm_bool, 'tmrca_meanmut_pat_spec']
    tmrca_spec_df.loc[is_ehorm_bool, 'ci95_tmrca_fin'] = tmrca_spec_df.loc[is_ehorm_bool, 'ci95_meanmut_tmrca_pat_spec']
    tmrca_spec_df.loc[is_ehorm_bool, 'mut_rate_fin'] = tmrca_spec_df.loc[is_ehorm_bool, 'mut_rate_obs']
    tmrca_spec_df.loc[is_ehorm_bool, 'ci95_mut_rate_fin'] = tmrca_spec_df.loc[is_ehorm_bool, 'mut_rate_obs_ci95']
    
    ## map which species has which clock gotten // hardcoded
    tmrca_spec_inf_clades['clock_for_inference'] = tmrca_spec_inf_clades['species'].map(annotate_clock)
    tmrca_spec_inf_clades['plt_color'] = tmrca_spec_inf_clades['clock_for_inference'].map(color_of_clock)
    
    return tmrca_spec_df

def extract_listentry_pd_col(pd_col, n):
    return pd_col.apply(lambda x: x[n] if len(x) > n else None)

############################################################
############################################################
############################################################

## get genome Lengths and calculate frequency where genome length is required
genome_length_dict = {}
for ref_genome_folder in glob.glob(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/metadata/reference_genomes/*240825/')): 
    species = ref_genome_folder.split('/')[-2].split('_')[1]
    [_, genomeLength, _] = apy.genomestats(ref_genome_folder);
    genome_length_dict[species] = genomeLength


############################################################
############################################################
############################################################
## Set variables

confidence=0.95
basescale = 1000000 # rescale genome-wide molecular rate to human-readable numbers
mutation_rate_df = pd.DataFrame()
mutation_rate_inf_df = pd.DataFrame()

## ncbi search: ((("Pseudomonas aeruginosa") OR (Pseudomonas)) AND (("mutation rate") OR ("molecular clock") OR ("substitution rate") AND (("/year") OR ("per year")) 

mutationrates = {   'Apittii': [1.286e-6, [1.125e-6, 1.449e-6]],                                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8875366/ # 5/genomeLength, # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5320584/, # https://www.researchsquare.com/article/rs-80066/v1, #https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000050 
                    'Bthetaiotaomicron': [0.9/genome_length_dict['Bthetaiotaomicron-c1']],      # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6749991/
                    'Ecoli': [6.9e-7, [3.14e-7, 1.4e-6]],                                       # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5835743/ # https://link.springer.com/article/10.1007/s00239-019-09912-5#additional-information / # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6015860/
                    'Ehormaechei': [2.7/genome_length_dict['Ehormaechei-c1'], [2.5/genome_length_dict['Ehormaechei-c1'], 3.0/genome_length_dict['Ehormaechei-c1']]], # https://journals.asm.org/doi/10.1128/mbio.00542-18 # 4.53/genomeLength, # https://journals.asm.org/doi/full/10.1128/aac.02091-17?casa_token=txr7x-GjJz0AAAAA%3A-BhhwwF1Xqkt2ASZg5FJcZ5_Sw8N4g7EbKkL8XWTuHzFKZefp-0ATYGMMERHTSZGEhOyti7XaDwbqg # 3*10E-08, # https://pubmed.ncbi.nlm.nih.gov/29566147/#:~:text=Results%3A%20Mutation%20rates%20were%20high,Providencia%20spp.%2C%20Serratia%20spp.
                    'Kmichiganensis': [1.9e-6, [1.1e-6, 2.9e-6]],                               #[4.18e-7, [0.87e-7, 4.89e-7]] # https://academic.oup.com/gbe/article/9/3/574/2977333 ## [1.9e-6, [1.1e-6, 2.9e-6]], # https://journals.asm.org/doi/full/10.1128/aac.04292-14 : 1.9 × 10−6 substitutions/called site/year (95% credibility interval [95% CI], 1.1 × 10−6 to 2.9 × 10−6) # https://www.embopress.org/doi/full/10.15252/emmm.201404767#
                    'Smaltophilia': [np.nan],
                    'Mmorganii': [np.nan],
                    'Pmirabilis': [np.nan],
                    'Paeruginosa': [2.92/genome_length_dict['Paeruginosa-c1'], [2.22/genome_length_dict['Paeruginosa-c1'], 3.62/genome_length_dict['Paeruginosa-c1']]],                  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4556809/ # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6015860/
                    'Saureus': [2.72/genome_length_dict['Saureus-ST97-c3'], [1.64/genome_length_dict['Saureus-ST97-c3'], 4.42/genome_length_dict['Saureus-ST97-c3']]] #[1.22e-6, [6.04e-7, 1.86e-6]]                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1401006111 ## https://www.science.org/doi/full/10.1126/science.1182395?casa_token=e2Ki2OmURhwAAAAA%3A8boxk8Ng3kLy5KBXyrk0jipPU6tXmRGOu-vkVec6HkMjs2pexzNdhrCpW5h8Qz7fDgS3oxIoihWGFw # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4752609/ # https://academic.oup.com/mbe/article/28/5/1593/1267325 # https://www.pnas.org/doi/epdf/10.1073/pnas.1113219109 / 2.05*10E-6: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6015860/
                }

## get a dict of the tp of infections in respect to hospital admission (note: only the clades with infectious isolates will be stated!)
inf_day_dict = {'P0007_P07_Apittii-c1_240825': 50.89,
                'P0007_P07_Bthetaiotaomicron-c1_240825': 50.831,
                'P0007_P07_Ecoli-ST4260-c2_240825': np.nan,
                'P0007_P07_Ecoli-ST58-c5_240825': np.nan, 
                'P0007_P07_Ecoli-ST95-c3_240825': np.nan,
                'P0007_P07_Ehormaechei-c1_240825': 38.688,
                'P0010_P10_Kmichiganensis-c1_240825': 3.127,
                'P0021_P21_Ecoli-ST23-c2_240825': 2.011,
                'P0021_P21_Ecoli-STnfx-c1_240825': np.nan,
                'P0021_P21_Pmirabilis-c2_240825': 2.011,
                'P0021_P21_Saureus-ST97-c3_240825': 2.011}

## plot specific variables
fontsize = 12
markersize = 25
ci_lower_idx = 0
inf_color = '#d73027'
tmrca_col = '#EAECCC' # '#F0E6B2' ## paul tol's Sunset/nightfall midpoint color schemes 
mol_clock_col = {'Inferred': '#484848',   #'#696969', ## grey shades with same contrast ratio apart (8x30): https://dev.to/finnhvman/grayscale-color-palette-with-equal-contrast-ratios-2pgl
                 'Species-specific': '#818181', #'#A9A9A9',
                 'Genus-specific': '#c8c8c8'} #'#DCDCDC'}
species_clock = {'A. pittii': 'Genus-specific', 
                 'B. thetaiotaomicron': 'Genus-specific', 
                 'E. hormaechei': 'Inferred', 
                 'K. michiganensis': 'Genus-specific', 
                 'E. coli': 'Species-specific', 
                 'P. mirabilis': np.nan, 
                 'S. aureus': 'Species-specific'}

############################################################
############################################################
############################################################
## fetch files 
files = glob.glob(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/backup/*/*/*_mol_clock_data_noinf.csv'))
files_inf = glob.glob(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/backup/*/*/*_mol_clock_data_winf.csv'))

## subset data to certain date
files = sorted([file for file in files if '2025-08-24' in file.split('/')[-1]])
files_inf = sorted([file for file in files_inf if '2025-08-24' in file.split('/')[-1]])

############################################################
############################################################
############################################################

mutrate_noinf_obs_dict = {}
for file in files:
    req_isolate_num = 8

    subj_spec = file.split('/')[-3]
    ## read in molecular clock without infectious isolates
    mol_clock_data = read_mol_clock_data_file(file)

    ## read in molecular clock _with_ infectious isolates
    file_inf = [file for file in files_inf if subj_spec in file][0]
    mol_clock_data_inf = read_mol_clock_data_file(file_inf)
    
    ## mol_clock_data[: ,0] = time since hospital admission of patient
    ## mol_clock_data[: ,1] = time in days since first positive culture
    ## mol_clock_data[: ,2] = number of SNVs to MRCA (inferred via treetime)
    ## mol_clock_data[: ,3] = number of bases covered on the genome ≥ 8x

    unique_tps = set(mol_clock_data[:, 1]) ## get all unique dates since first observation of pathogen (without clinical dates!)
    
    ## remove from inf dataset the microbiome timepoints
    inf_entries_bool = ~np.isin(mol_clock_data_inf[:, 1], list(unique_tps))
    mol_clock_data_inf_only = mol_clock_data_inf[inf_entries_bool, :]
    unique_inf_tps = set(mol_clock_data_inf_only[:, 1]) ## note: it is days to the first isolate!
    ## get date if first clinical observation of pathogen
    if len(unique_inf_tps) > 0:
        first_inf_tp = min(unique_inf_tps)
    else:
        first_inf_tp = np.nan
    mutrate_noinf_obs_dict[subj_spec] = calculate_molecular_clock_rate(mol_clock_data, subj_spec, unique_tps)
    mut_rate_obs = mutrate_noinf_obs_dict[subj_spec][0]
    mut_rate_obs_ci95 = mutrate_noinf_obs_dict[subj_spec][1]
    if not subj_spec in ['P0007_P07_Ehormaechei-c1_240825']:
        print(f'Observed mutation rate for {subj_spec}: {mut_rate_obs} CI95 {mut_rate_obs_ci95}')
        mut_rate_obs = np.nan
        mut_rate_obs_ci95 = np.nan 
    
    mut_rate_dict = {}
    species_clade = subj_spec.split('_')[2]
    species = species_clade.split('-')[0]
    mut_rate_literature = mutationrates[species][0] * basescale
    try:
        mut_rate_ci95_literature = [ci * basescale for ci in mutationrates[species][1]]
    except:
        mut_rate_ci95_literature = [0, 0]

    if np.isnan(mut_rate_literature):
        print('No mutation rate for literature found, using the mean of the identified mutation rates across isolates')
        continue
        
    ## check if there is any timepoint within ≤48 h (needed to account for rounding errors due to the use of integers (days!); cutoff need to be checked again manually for ≤24h!) which has more than the number of isolates required
    if first_inf_tp == first_inf_tp:
        if not any([(len(np.where(mol_clock_data[:, 1] == tp)[0]) >= req_isolate_num) & (tp <= first_inf_tp+2) for tp in unique_tps]): 
            print(f'Did not find any timepoint with more than {req_isolate_num} isolates')
            continue
    else:
        ## lineage does not cause an infection
        if not any([(len(np.where(mol_clock_data[:, 1] == tp)[0]) >= req_isolate_num) for tp in unique_tps]): 
            print(f'Did not find any timepoint with more than {req_isolate_num} isolates')
            continue

    ## calculate tmrca
    tmp_mutation_rate_df = get_tmrca_per_tp(mol_clock_data, unique_tps, mut_rate_obs, mut_rate_obs_ci95, mut_rate_literature, mut_rate_ci95_literature, basescale, req_isolate_num, subj_spec, plot_histogram=True)

    ## combine all data to one df
    mutation_rate_df = pd.concat([mutation_rate_df, tmp_mutation_rate_df])

mutation_rate_df = mutation_rate_df.reset_index(drop = True)

## plot per species just the first TP and its tmrca with repect to the infection 
## get first tps
tmrca_spec_inf = mutation_rate_df.sort_values('timepoint').groupby('subj_spec').head(1)

tmrca_spec_inf = add_metadata_to_df(tmrca_spec_inf, inf_day_dict)

####################
## plot with all clades which have enough data 
legend_info_tmrca = [Line2D([0], [0], marker = 'x', color = inf_color, linewidth = 0, markeredgewidth = 2, label = 'HAI\ndiagnosis', markersize = 8)]
# legend_info_tmrca += [plt.Rectangle((0,0),1,1, facecolor = 'w', edgecolor = 'k', linewidth = 1, label = 'Inferred acquisition time\nof lineage in microbiome\n(mean & 95% CI)')]
legend_info_tmrca += [plt.Rectangle((0,0),1,1, facecolor = tmrca_col, edgecolor = 'w', linewidth = 0, label = 'tMRCA of pathogenic\nlineage in microbiome\n(mean & 95% CI)')]

legend_info_mol_clock = [Line2D([0], [0], linewidth = 0, label = 'Molecular clock\n(mean & 95% CI)')]
for label, col in mol_clock_col.items():
    legend_info_mol_clock += [Line2D([0], [0], marker = 'o', color = col, markerfacecolor = col, linewidth = 2, label = label, markersize = 8)]

## plot with clades/lineages which contributed to inf and had enough data only
tmrca_spec_inf_clades = tmrca_spec_inf[~tmrca_spec_inf['inf_date'].isna()].copy()

## literature clock, for Ehormaechei --> inferred clock 
tmrca_spec_inf_clades = combine_literature_inferred_clock_for_plot(tmrca_spec_inf_clades, ['E. hormaechei clade1 (P07)'], species_clock, mol_clock_col)

## Removal from hypermutator candidates 
candidate_hypermutators = pd.read_csv('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/candidate_hypermutators_based_on_outgroup_background.csv')
candidate_hypermutators = candidate_hypermutators.loc[candidate_hypermutators['is_candidate_hypermutator'] == True, 'genome'].values

tmrca_spec_inf_clades_nohypmut = tmrca_spec_inf_clades[~tmrca_spec_inf_clades['subj_spec'].isin(candidate_hypermutators) | ~tmrca_spec_inf_clades['mut_rate_obs'].isna()]

fig, axs = plt.subplots(figsize = (7.5, 6), nrows = 2, sharex = True, gridspec_kw={'height_ratios': [3, 1]})
ax = axs.flatten()
ax[0].axhline(0, ls = '--', lw = 1, c = 'k', alpha = 0.5)
ax[0].annotate(text = 'Hospital admission', xy = (4.4, 0.25), fontsize=fontsize, ha = 'right')
ax[0].scatter(x= tmrca_spec_inf_clades_nohypmut['species'], y = tmrca_spec_inf_clades_nohypmut['inf_date'], linestyle='None', s = 100, marker='x', color = inf_color)
ax[0].bar(x = tmrca_spec_inf_clades_nohypmut['species'], 
            height = extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['ci95_tmrca_fin'], 1) - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['ci95_tmrca_fin'], ci_lower_idx),
            bottom = tmrca_spec_inf_clades_nohypmut['days_in_hosp'] - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['ci95_tmrca_fin'], 1),
            color = tmrca_col, width = 0.65)
ax[0].scatter(x= tmrca_spec_inf_clades_nohypmut['species'], y = tmrca_spec_inf_clades_nohypmut['days_in_hosp'] - tmrca_spec_inf_clades_nohypmut['tmrca_fin'], color = '#555555', linestyle='None', s = 1400, marker='_')
ax[0].set_ylim((np.nanmin(tmrca_spec_inf_clades_nohypmut['days_in_hosp'] - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['ci95_tmrca_fin'], 1)) - 300,
            np.nanmax(tmrca_spec_inf_clades_nohypmut['days_in_hosp'] - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['ci95_tmrca_fin'], ci_lower_idx)) + 150,
            ))
ax[0].set_ylabel('Time since hospital admission\n(in days)')
ax[0].set_yscale('symlog')
ax[0].get_yaxis().set_major_formatter(ScalarFormatter())
ax[0].legend(handles = legend_info_tmrca, loc='center left', bbox_to_anchor=(1, 0.5), fancybox=False, shadow=False, fontsize = fontsize)
ax[0].set_ylabel(ax[0].get_ylabel(), fontsize=fontsize+1)
ax[0].set_title('Acquisition time', loc='right', ha = 'right', fontsize=fontsize-2)

ax[1].scatter(x= tmrca_spec_inf_clades_nohypmut['species'], y = tmrca_spec_inf_clades_nohypmut['mut_rate_fin'], linestyle='None', s = 35, lw = 1, marker='o', color = tmrca_spec_inf_clades_nohypmut['plt_color'], zorder = 10)
ax[1].vlines(x = tmrca_spec_inf_clades_nohypmut['species'], 
           ymin = extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['ci95_mut_rate_fin'], 0),
           ymax = extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['ci95_mut_rate_fin'], 1),
           color = tmrca_spec_inf_clades_nohypmut['plt_color'], lw = 2, zorder = 15)
ax[1].legend(handles = legend_info_mol_clock, loc='center left', bbox_to_anchor=(1, 0.5), fancybox=False, shadow=False, fontsize = fontsize)
ax[1].set_yticks(np.arange( np.ceil(max( 
                                    np.nanmax(extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['ci95_mut_rate_fin'], 1)), 
                                    np.nanmax(tmrca_spec_inf_clades_nohypmut['mut_rate_fin'])) 
                                    + 1))) 
ax[1].set_ylabel('Mutations/\nMb/year')
ax[1].set_ylabel(ax[1].get_ylabel(), fontsize=fontsize+1)
ax[1].set_title('Substitution rates', loc='right', ha = 'right', fontsize=fontsize-2)
[ax[idx].spines['top'].set_visible(False) for idx, a in enumerate(ax)]
[ax[idx].spines['right'].set_visible(False) for idx, a in enumerate(ax)]
plt.xticks(rotation = 90, rotation_mode = 'anchor', ha = 'right', va = 'center', fontstyle='italic')
plt.tick_params(labelsize=fontsize)
plt.tight_layout()
plt.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2025_06_09_tmrca_nonhypermut_patho_clades_min8isolates_lit_clocks_excEhclock_wclocks_10kMCsim.svg'))


## to get the exact dates when in respect to hospitalization what happened: 
tmrca_spec_inf_clades_nohypmut_sub = tmrca_spec_inf_clades_nohypmut[['species', 'inf_date', 'days_in_hosp', 'tmrca_fin', 'ci95_tmrca_fin']].copy()
tmrca_spec_inf_clades_nohypmut_sub['mrca_date_since_hosp'] = tmrca_spec_inf_clades_nohypmut['days_in_hosp'] - tmrca_spec_inf_clades_nohypmut['tmrca_fin']
tmrca_spec_inf_clades_nohypmut_sub['ci_lower_mrca_date_since_hosp'] = tmrca_spec_inf_clades_nohypmut['days_in_hosp'] - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut_sub['ci95_tmrca_fin'], 0)
tmrca_spec_inf_clades_nohypmut_sub['ci_upper_mrca_date_since_hosp'] = tmrca_spec_inf_clades_nohypmut['days_in_hosp'] - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut_sub['ci95_tmrca_fin'], 1)
tmrca_spec_inf_clades_nohypmut_sub.loc[tmrca_spec_inf_clades_nohypmut_sub['ci_lower_mrca_date_since_hosp'].isna(), 'ci_lower_mrca_date_since_hosp'] = tmrca_spec_inf_clades_nohypmut_sub.loc[tmrca_spec_inf_clades_nohypmut_sub['ci_lower_mrca_date_since_hosp'].isna(), 'days_in_hosp']
tmrca_spec_inf_clades_nohypmut_sub.loc[tmrca_spec_inf_clades_nohypmut_sub['ci_upper_mrca_date_since_hosp'].isna(), 'ci_upper_mrca_date_since_hosp'] = tmrca_spec_inf_clades_nohypmut_sub.loc[tmrca_spec_inf_clades_nohypmut_sub['ci_upper_mrca_date_since_hosp'].isna(), 'days_in_hosp']

tmrca_spec_inf_clades_nohypmut_sub['mrca_date_before_inf'] = tmrca_spec_inf_clades_nohypmut_sub['mrca_date_since_hosp'] - tmrca_spec_inf_clades_nohypmut_sub['inf_date']
tmrca_spec_inf_clades_nohypmut_sub['ci_lower_mrca_date_before_inf'] = tmrca_spec_inf_clades_nohypmut_sub['ci_lower_mrca_date_since_hosp'] - tmrca_spec_inf_clades_nohypmut_sub['inf_date']
#tmrca_spec_inf_clades_nohypmut_sub.loc[tmrca_spec_inf_clades_nohypmut_sub['ci_lower_mrca_date_before_inf'].isna(), 'ci_lower_mrca_date_before_inf'] = tmrca_spec_inf_clades_nohypmut_sub.loc[tmrca_spec_inf_clades_nohypmut_sub['ci_lower_mrca_date_before_inf'].isna(), 'days_in_hosp']

tmrca_spec_inf_clades_nohypmut_sub['mrca_date_since_hosp_cilower'] =  tmrca_spec_inf_clades_nohypmut['days_in_hosp'] - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['ci95_tmrca_fin'], 0)
tmrca_spec_inf_clades_nohypmut_sub['mrca_date_since_hosp_cihigher'] =  tmrca_spec_inf_clades_nohypmut['days_in_hosp'] - extract_listentry_pd_col(tmrca_spec_inf_clades_nohypmut['ci95_tmrca_fin'], 1)

for rowidx, row in tmrca_spec_inf_clades_nohypmut_sub.iterrows():
    print(f"{row['species']} was first observed {row['days_in_hosp']} days post hospitalization and estimated to be acquired {round(row['mrca_date_since_hosp'], 1)} days (CI: {round(row['mrca_date_since_hosp_cilower'], 1)}-{round(row['mrca_date_since_hosp_cihigher'], 1)}) post hospitalization. The inferred MRCA predates the infection onset by {round(row['mrca_date_before_inf'], 1)} days.")

## correlate molecular clocks with tmrca values
tmrca_spec_inf_clades_nohypmut[['mut_rate_fin', 'tmrca_fin']] = tmrca_spec_inf_clades_nohypmut[['mut_rate_fin', 'tmrca_fin']].astype(float)
# stats.pearsonr(tmrca_spec_inf_clades_nohypmut['mut_rate_fin'], tmrca_spec_inf_clades_nohypmut['tmrca_fin'])
print(stats.spearmanr(tmrca_spec_inf_clades_nohypmut['mut_rate_fin'], tmrca_spec_inf_clades_nohypmut['tmrca_fin']))
