
## Analyse MIC data

## TODO: 
## Make a separate genotype list to use for input into plot!!!

import re
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

def replace_second_occurrence(text, pattern, replacement, replace_every_nth = 2):
    parts = re.split(f'({pattern})', text)
    count = 0
    for i, part in enumerate(parts):
        if part == pattern:
            count += 1
            if count % 2 == 0:
                parts[i] = replacement 
    return ''.join(parts)

def get_fold_change(df, timepoints_l, ab_l, genotype, exclude_genotype = []):
    time_of_interest = df['colonization_days'].isin(timepoints_l)
    ab_of_interest = df['AB'].isin(ab_l)
    mutation_of_interest = df['genotype_all'].str.contains(genotype)
    if exclude_genotype:
        for idx, gt in enumerate(exclude_genotype):
            if idx == 0:
                exclude_mut = df['genotype_all'].str.contains(gt)
            else:
                exclude_mut += df['genotype_all'].str.contains(gt)
        max_mic_no_mut = df[time_of_interest & ab_of_interest & ~mutation_of_interest & ~exclude_mut].groupby('AB')['MIC'].max()
    else:
        max_mic_no_mut = df[time_of_interest & ab_of_interest & ~mutation_of_interest].groupby('AB')['MIC'].max()
    max_mic_mut = df[time_of_interest & ab_of_interest & mutation_of_interest].groupby('AB')['MIC'].max()
    fold_change_ampD = pd.merge(max_mic_no_mut, max_mic_mut, how = 'outer', left_index=True, right_index=True, suffixes = [f'_no_{genotype}', f'_{genotype}'])
    fold_change_ampD['fold_change'] = 2 ** (fold_change_ampD[f'MIC_{genotype}'] - fold_change_ampD[f'MIC_no_{genotype}'])
    
    print(fold_change_ampD.to_string())

################################################
## set plot params
################################################
plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'


################################################
## load AB treatment file 
################################################
ehorm_treatment = pd.read_csv('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/MIC_Assays/metadata/2024_02_ABtreatment_P07_in_respect_to_Ehorm_colonization.csv')

################################################
## load patient 7 metadata
################################################
p07_enddate = pd.read_csv('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/RedCap_Data/Proc/P0007/endpoint/P0007_redcap_endpoint.csv')
p07_ehorm_inf = pd.read_csv('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/RedCap_Data/Proc/P0007/P0007_diagnostic_smpls.csv')
last_day = p07_enddate.loc[p07_enddate['redcap_event_name'] == 'studienende_te_arm_1', 'endpoint_endpointdate'].astype('datetime64[s]').values[0]
inf_date = p07_ehorm_inf.loc[p07_ehorm_inf['species'] == 'Enterobacter_cloacae', 'date'].astype('datetime64[s]').values[0]
last_day_post_inf = (last_day - inf_date).astype('timedelta64[D]').astype(int)+1

################################################
## load mic data
################################################
ehorm_mic = pd.read_csv('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/MIC_Assays/Raw/2024_07_04_MIC_Results_long_w_image_info_sorted.csv')
ehorm_mic = ehorm_mic[(ehorm_mic['ID'] != 'P07_2244') & (ehorm_mic['Species'] == 'E. hormaechei')] ## non existent lineage ID!!! Should be 2240 but cannot be certainly be recovered!

################################################
## load metadata file
################################################
metadata = pd.read_csv('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/MIC_Assays/Raw/2024_07_metadata_HAIpatientsP07P21_MIC_isolates.csv')
metadata = metadata.drop(['Patient', 'Species'], axis = 1)

## merge mic results with metadata of strains
ehorm_mic = pd.merge(ehorm_mic, metadata, how = 'left', on = 'ID')

## generate long df for plotting and remove nans
ehorm_mic['AB'] = ehorm_mic['AB'].replace('Ceftazidim', 'Ceftazidime').replace('Cefotaxim', 'Cefotaxime').replace('Piperacillin/Tazobactam', 'Piperacillin/\nTazobactam')
ehorm_mic = ehorm_mic[(~ehorm_mic['MIC'].isna()) & (ehorm_mic['MIC'] != '-')]
ehorm_mic['MIC'] = ehorm_mic['MIC'].astype(str).str.replace('>', '').astype(float)

## rename genotypes for lineages where no genotype given in metadata (--> separate lineages)
ehorm_mic.loc[ehorm_mic['genotype'].isna(), 'genotype'] = 'Lineage ' + ehorm_mic.loc[ehorm_mic['genotype'].isna(), 'branch'].astype(str)


AB_color_dict = {'Piperacillin/\nTazobactam': "#44BB99",  #P // 1 saturation of penicillin
                 'Imipenem': '#909000',                 #C // 1 saturation of Carbapenem; v*0.85
                 'Meropenem': '#AAAA00',                #C // 1 saturation of Carbapenem
                 'Cefotaxime': '#BBCC33',               #ce // 1 saturation of Cephalosporin
                 'Ceftazidime': '#9ead2b',              #ce // 1 saturation of Cephalosporin; v*0.85
                 'Vancomycin': '#EEDD88',               #Gly // 1 saturation of Glycopeptide
                 'Ciprofloxacin': '#FFAABB',            #Flu // 1 saturation of Fluorochinolone
                 'Erythromycin': '#77AADD',             #Mac // 1 saturation of Macrolide
}

## see https://www.eucast.org/clinical_breakpoints version 15.0; given the values for resistance (R)
eucast_breakpoints = {'Piperacillin/\nTazobactam': 8,
                      'Imipenem': 4,
                      'Meropenem': 8,
                      'Cefotaxime': 2,
                      'Ceftazidime': 4,
                      'Vancomycin': np.nan,
                      'Ciprofloxacin': 0.25,
                      'Erythromycin': np.nan}

################################################
## Load SNPs
################################################
## update genotypes from metadata file to ensure validity & get frequency of a given genotype
snp_table = pd.read_csv('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/backup/P0007_P07_Ehormaechei-c1_240825/2025_03_10_P07_Ehormaechei-c1_240825/snp_table.csv')
non_isolate_cols = ['chr', 'pos', 'type', 'muts', 'locustag', 'gene', 'loc1', 'loc2', 'strand', 'product', 'nt_pos', 'aa_pos', 'nt_ref', 'nt_anc', 'nt_alt']
isolate_cols = [col for col in snp_table.columns if col not in non_isolate_cols]
isolate_to_tp_dict = {isolate_col.split('__')[-1]: isolate_col.split('_')[1] for isolate_col in isolate_cols}
snp_table['muts'] = snp_table['muts'].str.replace("\\[|\\]|'|", '', regex = True).str.replace('.', '')
snp_table['gene'] = snp_table['gene'].replace('citG,mdcB', 'citG')

## replace locustag prefix with LT
for row_idx, row in snp_table.iterrows():
    new_gene_tag = ''
    for gene_idx, gene in enumerate(row['gene'].split(';')):
        if gene == '.':
            lt = row["locustag"].split(';')[gene_idx].split(':')[-1]
            lt = lt.replace('MDIMIJ', '').replace('_', '')
            new_gene_tag += f';LT{lt}'
        else:
            new_gene_tag += f';{gene}'
    snp_table.loc[row_idx, 'gene'] = new_gene_tag.strip('; ')


genotype_dict = {}
for isolate in isolate_cols:
    isolate_id = isolate.split('__')[-1]
    genotype_dict[isolate_id] = []
    for row_idx, row in snp_table.iterrows():
        if len(np.unique(row[isolate_cols])) <= 1: ## all isolates have same allele --> no mutation amonst the isolates but just between isolates and outgroup --> mutation is not ancestral!!!
            continue
        anc_nt = row['nt_anc']
        if anc_nt == '?':
            anc_nt = row['nt_ref']
        if len(row['type']) > 1:
            muts = row['muts'].split(', ')
            for mut, alt_nt in zip(muts, row['nt_alt']):
                if ( (row[isolate] != anc_nt) and (row[isolate] == alt_nt) and (anc_nt != alt_nt) ) or ( (row[isolate] == row['nt_ref']) and (anc_nt == alt_nt) ): ## ancestral is written as alternative in snp_table --> repolarize
                    mut = row['gene'] + '[' + mut + ']'
                    mut = mut.replace('gyrA[F83S]', 'gyrA[S83F]').replace('gyrA[F83Y]', 'gyrA[S83Y]').replace(';', '|') ## give correct orientation to ancestral not reference!
                    genotype_dict[isolate_id].append(mut)
        if row_idx == 59:
            break
        elif ( (row[isolate] != anc_nt) and (row[isolate] == row['nt_alt']) and (anc_nt != row['nt_alt']) ) or ( (row[isolate] == row['nt_ref']) and (anc_nt == row['nt_alt']) ): ## ancestral is written as alternative in snp_table --> repolarize
            gene = row['gene']
            if row['muts'] == row['muts']:
                mut = gene + '[' + row['muts'] + ']'
            else:
                mut = gene
            mut = mut.replace('gyrA[F83S]', 'gyrA[S83F]').replace('gyrA[F83Y]', 'gyrA[S83Y]').replace(';', '|') ## give correct orientation to ancestral not reference!
            genotype_dict[isolate_id].append(mut)
genotype_dict_string = {gt: '; '.join(muts) for gt, muts in genotype_dict.items()}

ehorm_mic['genotype'] = ehorm_mic['ID'].map(genotype_dict_string)


################################################
## Load Indels
################################################
indel_table = pd.read_csv('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/backup/P0007_P07_Ehormaechei-c1_240825/2024_10_09_P07_Ehormaechei-c1_240825/indel_table.csv')

## add ampD genptypes which have been additionally identified 
## positions are nt positions!
additional_ampD_indels =  {'P07_1149': 'ampD_480',
                            'P07_1049': 'ampD_279',
                            'P07_1553': 'ampD_279'}

## indel table --> change seq to its mut type 
indel_table_metadata = ['gene_num', 'gene_num_global', 'chr', 'pos', 'chr_boundary', 'indel_ref', 'indel_anc', 'indel_alt', 'max_indel_GL_diff_to_ref',
                        'indel_size_gp', 'type', 'product', 'gene', 'protein_id', 'strand', 'loc1', 'loc2', 'sequence', 'note', 'locustag', 'orthologtag',
                        'translation', 'indel_pos_start', 'indel_pos_end', 'alt_start', 'alt_stop', 'alt_stop_pos', 'frameshift', 'mut_translation', 'mut',
                        'num_isolates', 'gene1', 'locustag1', 'orthologtag1', 'product1', 'distance1', 'gene2', 'locustag2', 'orthologtag2', 'product2', 'distance2']
indel_table_metadata_avail = [col for col in indel_table.columns if col in indel_table_metadata]
indel_table_samplecol = [col for col in indel_table.columns if col not in indel_table_metadata]
indel_table.loc[indel_table['gene'].isna(), 'gene'] = indel_table.loc[indel_table['gene'].isna(), 'gene1'] + '|' + indel_table.loc[indel_table['gene'].isna(), 'gene2']
indel_table.loc[indel_table['locustag'].isna(), 'locustag'] = (indel_table.loc[indel_table['locustag'].isna(), 'locustag1'] + ':' + indel_table.loc[indel_table['locustag'].isna(), 'distance1'].astype(str) + '|' + 
                                                        indel_table.loc[indel_table['locustag'].isna(), 'locustag2'] + ':' + indel_table.loc[indel_table['locustag'].isna(), 'distance2'].astype(str))

for row_idx, row in indel_table.iterrows():
    new_gene_tag = ''
    for gene_idx, gene in enumerate(row['gene'].split('|')):
        if gene == '.':
            lt = row["locustag"].split('|')[gene_idx].split(':')[0]
            lt = lt.replace('MDIMIJ', '').replace('_', '')
            new_gene_tag += f'|LT{lt}'
        else:
            new_gene_tag += f'|{gene}'
    indel_table.loc[row_idx, 'gene'] = new_gene_tag.strip('| ')

cols_to_explode_on = ['indel_alt', 'alt_stop', 'indel_pos_start', 'indel_size_gp']
## convert list strings (e.g.: "['1', 'a', 'c']") to real lists in df
indel_table.loc[indel_table['type'] != 'G', 'alt_stop'] = indel_table.loc[indel_table['type'] != 'G', 'indel_alt'].apply(lambda x: ', '.join( ['intergenic' for i in (x.split(' ')) ]) )
indel_table.loc[indel_table['type'] != 'G', 'indel_pos_start'] = indel_table.loc[indel_table['type'] != 'G', 'indel_alt'].apply(lambda x: ', '.join( ['' for i in (x.split(' ')) ]) )
for col in cols_to_explode_on:
    indel_table[col] = indel_table[col].str.replace(r"\[|\]|'|,", '', regex=True).str.replace(r'\s+', ', ', regex = True).str.split(', ')


indel_table_exploded = indel_table.explode(cols_to_explode_on).reset_index(drop = True)
for col in indel_table_samplecol:
    indel_table_exploded[col] = np.where(indel_table_exploded[col] == indel_table_exploded['indel_alt'], indel_table_exploded['alt_stop'].values, np.nan)
indel_table_exploded['genotype'] = indel_table_exploded['gene'] + '_' + indel_table_exploded['indel_pos_start'].astype(str) + '_' + (indel_table_exploded.groupby(['chr', 'pos'])['indel_ref'].cumcount()+1).astype(str)
## generate dict of samples which have the indel genotypes as values
indel_table_exploded_long = indel_table_exploded[['genotype'] + indel_table_samplecol].melt('genotype')
sample_indel_genotype_dict = indel_table_exploded_long.dropna(subset='value').groupby('variable')['genotype'].agg(lambda x: '; '.join(x)).to_dict()

for smpl, indel in additional_ampD_indels.items():
    if smpl in sample_indel_genotype_dict.keys():
        sample_indel_genotype_dict[smpl] = sample_indel_genotype_dict[smpl] + '; ' + indel
    else:
        sample_indel_genotype_dict[smpl] = indel

################################################
## Extract genotype frequency per timepoint for plot
################################################
unique_genotypes = {sid: f'{snv_gt}; {sample_indel_genotype_dict[sid]}' if sid in sample_indel_genotype_dict.keys() else snv_gt for sid, snv_gt in genotype_dict_string.items()}

gt_freq_p_tp = {tp: {} for tp in sorted(np.unique(list(isolate_to_tp_dict.values())))}
for isolate, gt in unique_genotypes.items():
    if gt in gt_freq_p_tp[isolate_to_tp_dict[isolate]].keys():
        gt_freq_p_tp[isolate_to_tp_dict[isolate]][gt] += 1
    else:
        gt_freq_p_tp[isolate_to_tp_dict[isolate]][gt] = 1

for tp in gt_freq_p_tp.keys():
    total_cnt = sum([cnt for cnt in list(gt_freq_p_tp[tp].values())])
    gt_freq_p_tp[tp] = {gt: cnt/total_cnt for gt, cnt in gt_freq_p_tp[tp].items()}


################################################
## process treatment data for plotting
################################################
ehorm_treatment['ab'] = ehorm_treatment['ab'].str.title() ## capitalize first letter
ehorm_treatment['ab'] = ehorm_treatment['ab'].replace('Ceftazidim', 'Ceftazidime').replace('Cefotaxim', 'Cefotaxime').replace('Piperacillin/Tazobactam', 'Piperacillin/\nTazobactam')
ehorm_treatment = ehorm_treatment.sort_values(['ab'])
## generate new df structure for plotting
ehorm_treatment = pd.concat([pd.DataFrame({'days': [row['start_respective_colonization'], row['end_respective_colonization']], 
                                                'AB_period_ID': [i/10, i/10], ## necessary to 
                                                'AB': [row['ab'], row['ab']]}) 
                                for i, (_, row) in enumerate(ehorm_treatment.iterrows())], ignore_index=True)

## set timepoint of infection to tp 0
timepoint_of_infection = ehorm_mic.loc[ehorm_mic['ID'].str.startswith('L000'), 'colonization_days'].unique()[0]
ehorm_mic['colonization_days'] = ehorm_mic['colonization_days']-timepoint_of_infection
ehorm_treatment['days'] = ehorm_treatment['days']-timepoint_of_infection

negative_days = 14
ehorm_mic = ehorm_mic[~ehorm_mic['AB'].isin(['Ceftazidime', 'Erythromycin'])]
ehorm_mic['MIC'] = np.log2(ehorm_mic['MIC'].astype(str).str.replace('>', '').astype(float))

#############################################
#############################################
#############################################
## Plot multiple axes & all ABs at once With heatmap as genotype indicator
#############################################
#############################################
#############################################

snv_indel_marker = {'SNV': 'o', 'Indel': 's'}
snv_indel_marker_sizes = {'SNV': 100, 'Indel': 100}


exclude_mpip = True
unique_days = sorted(ehorm_mic['colonization_days'].unique())

if exclude_mpip:
    unique_abs = ['Ciprofloxacin', 'Cefotaxime', 'Piperacillin/\nTazobactam']
    ab_plotted = 'CP_CTX_PT'
    mutation_of_interest = ['fimZ', 'gyrA', 'LT12530', 'ampD', 'iscR', 'acrR', 'romA|acrR', 'yjfF']
else:
    unique_abs = ['Ciprofloxacin', 'Cefotaxime', 'Piperacillin/\nTazobactam', 'Imipenem', 'Meropenem']
    ab_plotted = 'CP_CTX_PT_IP_MP'
    mutation_of_interest = [] ## empty list to plot all genotypes

## add indels if required 
plt_df = ehorm_mic.copy()
plt_df['genotype_indel'] = plt_df['ID'].map(sample_indel_genotype_dict).fillna('')
plt_df['genotype_all'] = plt_df['genotype'] + '; ' + plt_df['genotype_indel']
plt_df['genotype_all'] = plt_df['genotype_all'].str.strip('; ')
plt_df['genotype_short'] = plt_df['genotype_all'].str.replace(r'\[.{1,5}\]|_..{1,5}_', '', regex = True)
## sort genotypes based on their first occurences
ordered_gt_dict = {gt:idx for idx, gt in enumerate(plt_df.sort_values('colonization_days')['genotype_all'].unique())}
plt_df = plt_df.sort_values('genotype_all', key = lambda x: x.map(ordered_gt_dict)).reset_index(drop = True)

plt_df_gt_to_d = plt_df[['colonization_days', 'genotype_all', 'genotype', 'genotype_indel', 'ID']].copy()

## extract the order of genotypes --> sorted to follow the increase of Ciprofloxaicin per date
ordered_GT_per_day_over_CP = plt_df[plt_df['AB'] == 'Ciprofloxacin'].sort_values(['colonization_days', 'MIC'])['genotype_all'].unique()

## prepare for the heatmap below 
## get unique mutations
unique_snv_gt = plt_df_gt_to_d.sort_values('colonization_days')['genotype'].unique()
snv_mut = list(dict.fromkeys([mut for genotype in unique_snv_gt for mut in genotype.split('; ')]))

unique_indel_gt = plt_df_gt_to_d.sort_values('colonization_days')['genotype_indel'].unique()
indel_mut = list(dict.fromkeys([mut for genotype in unique_indel_gt for mut in genotype.split('; ') if mut != '']))
genotype_to_idx = {gt: idx for idx, gt in enumerate((snv_mut+indel_mut)[::-1])}

## generate long df with  mutation per row
plt_df_gt_to_d['genotype_all_s'] = plt_df_gt_to_d['genotype_all'].str.replace(r'; $', '', regex = True).str.replace(r'; ', ';', regex = True).str.split(';')
plt_df_gt_to_d = plt_df_gt_to_d.explode(['genotype_all_s']).reset_index(drop = True) 
plt_df_gt_to_d = plt_df_gt_to_d.drop_duplicates(subset=['colonization_days', 'genotype_all_s', 'genotype_all'])

plt_df_gt_to_d['mutated_gene'] = plt_df_gt_to_d['genotype_all_s'].str.split('[').str[0].str.split('_').str[0] ## get gene name
plt_df_gt_to_d.loc[plt_df_gt_to_d['genotype_all_s'].str.contains('\\['), 'mut_pos'] = '[' + plt_df_gt_to_d.loc[plt_df_gt_to_d['genotype_all_s'].str.contains('\\['), 'genotype_all_s'].str.split('[').str[-1] ## extract SNV position in protein
plt_df_gt_to_d.loc[plt_df_gt_to_d['genotype_all_s'].str[4] == '_', 'mut_pos'] = (plt_df_gt_to_d.loc[plt_df_gt_to_d['genotype_all_s'].str[4] == '_', 'genotype_all_s'].str.split('_').str[1].astype(int) / 3).astype(int) ## get all indels
unique_snvs = np.unique([gt.strip('\n ') for gtl in plt_df_gt_to_d['genotype'].str.split(';') for gt in gtl])
plt_df_gt_to_d['mutation'] = np.where(plt_df_gt_to_d['genotype_all_s'].isin(unique_snvs), 'SNV', 'Indel') ## map back if mutation is indel or SNV

if mutation_of_interest == []:
    mutation_of_interest = [mut for mut in sorted(plt_df_gt_to_d['mutated_gene'].unique()) if 'basal' not in mut] ## remove 'basal' annotation as distinct mutation
mut_gene_to_idx = {gt: idx for idx, gt in enumerate(mutation_of_interest[::-1]) if gt != ''} ## get dict to convert mutations to unique list entries for mutations of interest
plt_df_gt_to_d['mutated_gene_idx'] = plt_df_gt_to_d['mutated_gene'].map(mut_gene_to_idx)


## calculate the number of genotypes per timepoint 
width_ratio_l = plt_df.sort_values('colonization_days').groupby('colonization_days')['genotype_all'].nunique().to_list()
width_ratio_l = width_ratio_l ## add a width for the first empty ax
## calculate the number of ABs 
max_mic_diff_per_ab = {ab: max(plt_df.loc[plt_df['AB'] == ab, 'MIC'].max(), np.log2(eucast_breakpoints[ab])) - min(plt_df.loc[plt_df['AB'] == ab, 'MIC'].min(), np.log2(eucast_breakpoints[ab])) for ab in plt_df['AB'].unique()} #plt_df.groupby('AB')['MIC'].apply(lambda x: max(x) - min(x)).to_dict()
height_ratio_l = [max_mic_diff_per_ab[ab] for ab in unique_abs]
height_ratio_l = [3] + height_ratio_l + [len(mutation_of_interest)*1.25] ## add a height for the first empty ax and the last

## MICs per timepoint (--> note there are timepoints where there are multiple different genotypes)
fig = plt.figure(figsize=(14, sum(height_ratio_l)/3))
if ab_plotted == '': #'CP_CTX_PT': ## just plot frequencies if we plot all!
    height_ratio_l = height_ratio_l + [5]
    gs = gridspec.GridSpec(len(unique_abs)+3, len(unique_days), height_ratios=height_ratio_l, width_ratios= width_ratio_l, hspace = 0.3, wspace = 0.2)
else:
    gs = gridspec.GridSpec(len(unique_abs)+2, len(unique_days), height_ratios=height_ratio_l, width_ratios= width_ratio_l, hspace = 0.3, wspace = 0.2)
ax1 = plt.subplot(gs[0, :])  # The first row spans all columns
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
sns.lineplot(data = ehorm_treatment[~ehorm_treatment.duplicated(keep = False)], x = 'days', y = 'AB', 
            hue = 'AB', palette = AB_color_dict, 
            units = 'AB_period_ID', lw = 5,
            estimator = None, 
            ax = ax1)
## add single treatments
sns.scatterplot(data = ehorm_treatment[ehorm_treatment.duplicated(keep = 'last')], x = 'days', y = 'AB', 
            hue = 'AB', marker = '|', palette = AB_color_dict, 
            linewidth = 5, s = 50, 
            ax = ax1)
ax1.set_ylim(-2, len(ehorm_treatment['AB'].unique())+1)
ax1.set_xlim((min(ehorm_treatment['days']), max(max(unique_days)+3, last_day_post_inf)))
ax1.xaxis.set_major_locator(ticker.MultipleLocator(7)) ## offset based on decimals in ehorm_treatment_shifted['days']
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax1.grid(which='major', axis = 'x', linestyle='-', linewidth=0.5, color='black', zorder = 1)
ax1.grid(which='minor', axis = 'x', linestyle=':', linewidth=0.5, color='black', zorder = 1)
ax1.xaxis.tick_top()
ax1.tick_params(width=0, pad=0)
ax1.xaxis.set_label_position('top')
ax1.tick_params(axis='both',which='both',bottom=False,left= False, top=True, labelbottom=False, labeltop=True, labelleft=False)
ax1.set_ylabel('Antibiotic\ntreatment', rotation = 90)
ax1.set_xlabel('Time (in days)')

# Second row: n subplots
mic_axes = []
xy_start_figs = []
xy_end_figs = []
for ab_idx, AB in enumerate(unique_abs):
    plt_df_tmp = plt_df[plt_df['AB'] == AB]
    if plt_df_tmp.empty:
        continue
    for day_idx, day in enumerate(unique_days):
        ## extract day of interest and calculate the mean mic
        plt_df_tmp_date = plt_df_tmp[plt_df_tmp['colonization_days'] == day].copy()
        ordered_GT_over_CP_present = [gt for gt in ordered_GT_per_day_over_CP if gt in plt_df_tmp_date['genotype_all'].unique()] ## get just the genotypes present on the date to sort the df by that
        plt_df_tmp_date['genotype_all'] = pd.Categorical(plt_df_tmp_date['genotype_all'], ordered_GT_over_CP_present, ordered = True)
        mean_mics = plt_df_tmp_date.groupby('genotype_all', observed = False)['MIC'].mean().reset_index()
        ## add subplot
        ax = plt.subplot(gs[ab_idx+1, day_idx])  # Each column gets its own plot in the second row
        sns.scatterplot(data = mean_mics, x = 'genotype_all', y = 'MIC', marker = '_', color = AB_color_dict[AB], linewidth = 3.5, s = 250, zorder = 10, ax = ax)
        sns.scatterplot(data = plt_df_tmp_date, x = 'genotype_all', y = 'MIC', color = AB_color_dict[AB], s=40, edgecolor = 'k', linewidth = 0.2, alpha = 0.8, zorder = 15, ax = ax)

        ax.axhline(y=np.log2(eucast_breakpoints[AB]), ls = '--', lw = 2, c = 'darkgrey')
        ax.set_xticks(np.arange(len(plt_df_tmp_date['genotype_all'].cat.categories)))
        ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, ha = 'right', rotation_mode = 'anchor')
        ax.set_xlim(-0.5, len(plt_df_tmp_date['genotype_all'].unique())-0.5) 
        ax.set_xlabel(None)
        ax.set_ylim(min(min(plt_df_tmp['MIC']), np.log2(eucast_breakpoints[AB]))-0.5, max(max(plt_df_tmp['MIC']), np.log2(eucast_breakpoints[AB])) +0.5) 
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.grid(which='major', axis = 'y', linestyle=':', linewidth=0.5, color='black', zorder = 1)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if day_idx == 0:
            ax.set_ylabel(f'{AB}\nMIC [log2(µg/ml)]')
        else:
            ax.tick_params(axis='y',which='both',left= False,labelleft=False)
            ax.set_ylabel(None)
        ax.tick_params(axis='x',which='both',bottom=False,labelbottom=False)
        mic_axes.append(ax)

        if ab_idx == 0:
            ## add a connection from subplot to the antibiotic treatment plot
            xy_start = ax.transAxes.transform((0.5, 1))  ## center of current axis
            xy_end = ax1.transData.transform((day, -1))  # data coordinate of antibiotic treatment plot
            xy_start = xy_start - np.array([65, 120])
            xy_end = xy_end - np.array([65, 120])
            xy_start_figs.append(fig.transFigure.inverted().transform(xy_start))
            xy_end_figs.append(fig.transFigure.inverted().transform(xy_end))

gt_axes = []
for day_idx, day in enumerate(unique_days):
    plt_df_gt = plt_df_gt_to_d[(plt_df_gt_to_d['colonization_days'] == day)].copy()

    ordered_GT_over_CP_present = [gt for gt in ordered_GT_per_day_over_CP if gt in plt_df_gt['genotype_all'].unique()] ## get just the genotypes present on the date to sort the df by that
    plt_df_gt['genotype_all'] = pd.Categorical(plt_df_gt['genotype_all'], ordered_GT_over_CP_present, ordered = True)
    plt_df_gt = plt_df_gt[~plt_df_gt_to_d['mutated_gene_idx'].isna()] ## just after creating the categorical for the plt, remove unmapped 

    if ab_plotted == '': #'CP_CTX_PT':                              
        axhm = plt.subplot(gs[-2, day_idx])  # The last row spans all columns
    else:
        axhm = plt.subplot(gs[-1, day_idx])  # The last row spans all columns
    axhm.spines['right'].set_visible(False)
    axhm.spines['bottom'].set_visible(False)
    sns.scatterplot(data = plt_df_gt, x = 'genotype_all', y = 'mutated_gene_idx', style = 'mutation', markers = snv_indel_marker, size = 'mutation', sizes = snv_indel_marker_sizes, color = 'k', linewidth = 0, zorder = 10, ax = axhm) # ,
    for i, row in plt_df_gt.iterrows():
        if row['mut_pos'] != row['mut_pos']: ## skip positions which are intergenic
            continue
        axhm.annotate(row['mut_pos'], 
                    (row['genotype_all'], row['mutated_gene_idx']),
                    textcoords="offset points", 
                    xytext=(0, -11), size = 6.5,  
                    ha = 'center') 
    axhm.set_ylim((-0.5, max(mut_gene_to_idx.values())+0.5)) 
    axhm.set_yticks(list(mut_gene_to_idx.values()))
    axhm.set_yticklabels(list(mut_gene_to_idx.keys()))
    axhm.grid(which='major', axis = 'y', linestyle=':', linewidth=0.5, color='black', zorder = 1)
    axhm.set_xticks(np.arange(len(plt_df_gt['genotype_all'].cat.categories)))
    axhm.set_xticklabels(plt_df_gt['genotype_all'].cat.categories, rotation = 90, ha = 'right', rotation_mode = 'anchor')
    axhm.set_xlim(-0.5, len(plt_df_gt['genotype_all'].cat.categories)-0.5) 
    axhm.set_xlabel(None)
    axhm.spines['top'].set_visible(False)
    axhm.spines['right'].set_visible(False)
    axhm.tick_params(axis='x',which='both',bottom=False,labelbottom=False)
    if day_idx == 0:
        axhm.set_ylabel(f'Major mutation')
    else:
        axhm.tick_params(axis='y',which='both',left= False,labelleft=False)
        axhm.set_ylabel(None)
    if day_idx < len(unique_days):
        try:
            axhm.get_legend().remove()
        except:
            print()
    gt_axes.append(axhm)

## plot frequency of genotypes
if ab_plotted == '': #'CP_CTX_PT': ## just plot frequencies if we plot all!
    gt_freq_axes = []
    for day_idx, (day, tp) in enumerate(zip(unique_days, gt_freq_p_tp.keys())):
        plt_df_gt = plt_df_gt_to_d[(plt_df_gt_to_d['colonization_days'] == day)].copy()
        ordered_GT_over_CP_present = [gt for gt in ordered_GT_per_day_over_CP if gt in plt_df_gt['genotype_all'].unique()] ## get just the genotypes present on the date which were tested to sort the df by that
        
        plt_df_gt_freq = pd.DataFrame(gt_freq_p_tp[tp], index = ['freq']).T.reset_index()
        plt_df_gt_freq.columns = ['genotype_all', 'freq']
        plt_df_gt_freq['genotype_all'] = pd.Categorical(plt_df_gt_freq['genotype_all'], ordered_GT_over_CP_present, ordered = True)
        plt_df_gt_freq = plt_df_gt_freq[~plt_df_gt_freq['genotype_all'].isna()] ## remove unmapped values

        axbar = plt.subplot(gs[-1, day_idx])  # The last row spans all columns
        axbar.spines['right'].set_visible(False)
        axbar.spines['bottom'].set_visible(False)
        sns.barplot(data = plt_df_gt_freq, x = 'genotype_all', y = 'freq', color = 'k', linewidth = 0, zorder = 10, ax = axbar) # ,
        
        axbar.set_ylim((0, 1)) 
        axbar.grid(which='major', axis = 'y', linestyle=':', linewidth=0.5, color='black', zorder = 1)
        axbar.set_xticks(np.arange(len(plt_df_gt_freq['genotype_all'].cat.categories)))
        axbar.set_xticklabels(plt_df_gt_freq['genotype_all'].cat.categories, rotation = 90, ha = 'right', rotation_mode = 'anchor')
        axbar.set_xlim(-0.5, len(plt_df_gt_freq['genotype_all'].cat.categories)-0.5) 
        axbar.set_xlabel(None)
        axbar.spines['top'].set_visible(False)
        axbar.spines['right'].set_visible(False)
        axbar.tick_params(axis='x',which='both',bottom=False,labelbottom=False)
        if day_idx == 0:
            axbar.set_ylabel(f'Genotype frequency')
        else:
            axbar.tick_params(axis='y',which='both',left= False,labelleft=False)
            axbar.set_ylabel(None)
        if day_idx < len(unique_days):
            try:
                axbar.get_legend().remove()
            except:
                print()
        gt_freq_axes.append(axbar)

## generate legend
legend_info = [Line2D([], [], marker='None', linestyle='None')] ## header
legend_labels = ['Antibiotic']
for ab in sorted(AB_color_dict.keys()):
    if ab in ehorm_treatment['AB'].unique():
        legend_info += [Line2D([0], [0], marker = '_', color = AB_color_dict[ab], markerfacecolor = AB_color_dict[ab], label = '', linewidth=5, markersize = 10)] ## blood
        legend_labels += [ab]
legend_info += [Line2D([], [], marker='None', linestyle='None')] ## empty line
legend_info += [Line2D([], [], marker='None', linestyle='None')] ## header
legend_labels += ['', 'Mutation type\n(with AA position)']
for muttype in sorted(snv_indel_marker.keys()):
    if muttype in plt_df_gt_to_d['mutation'].unique():
        legend_info += [Line2D([0], [0], marker = snv_indel_marker[muttype], color = 'k', markerfacecolor = None, label = '', linewidth=0, markersize = 10)] ## blood
        legend_labels += [muttype]

ax1.legend(legend_info, legend_labels, loc = 'center left', bbox_to_anchor = (1, -1))
[plt.annotate('',xy=xy_end_fig, xycoords='figure fraction',
            xytext=xy_start_fig, textcoords='figure fraction',
            arrowprops=dict(arrowstyle="->", color='k', lw=1.5)) for xy_start_fig, xy_end_fig in zip(xy_start_figs, xy_end_figs)]
fig.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2025_06_09_P07_Ehormaechei_MIC_{ab_plotted}_genotypes.pdf'),bbox_inches = 'tight', transparent=True)
fig.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2025_06_09_P07_Ehormaechei_MIC_{ab_plotted}_genotypes.svg'),bbox_inches = 'tight', transparent=True)


## print mic differences of interest

## 12d post treatment
get_fold_change(plt_df, [4.75, 11.75], ['Cefotaxime', 'Piperacillin/\nTazobactam'], 'ampD')
get_fold_change(plt_df, [4.75, 11.75], ['Ciprofloxacin'], 'gyrA')

## post 2nd infection
get_fold_change(plt_df, [25.75, 39.75, 46.75], ['Ciprofloxacin'], 'romA\\|acrR', exclude_genotype = ['acrR_283']) ## intergenic
get_fold_change(plt_df, [25.75, 39.75, 46.75], ['Ciprofloxacin'], 'acrR_283', exclude_genotype = ['romA|acrR']) ## genic AA 94
