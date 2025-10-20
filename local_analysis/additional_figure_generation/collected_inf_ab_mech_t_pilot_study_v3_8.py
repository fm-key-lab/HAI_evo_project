## generate figure for sampled patients

## load modules
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from matplotlib.ticker import ScalarFormatter, FixedLocator
from datetime import timedelta
import seaborn as sns

plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

## infos and paths
saveme = True ## set to true if you want to save the plot
patients = ['7', '10', '21']
outpath = "~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/RedCap_Data/Proc/"

## set what to plot
figure_suffix = ''
plt_mbiome = True ## set to true if Mbiome sampling should be plotted
if not plt_mbiome:
    figure_suffix += '_noMbiome'
plt_treatment = True ## set to true if AB treatment should be plotted
if not plt_treatment:
    figure_suffix += '_noAB'
plt_legend = True ## set to true if legend should be plotted
if not plt_legend:
    figure_suffix += '_nolegend'


######################################
## FUNCTIONS
######################################

def to_end_of_day(date_str):
    ## read a date and set the time to the last minute of the day --> for discharge needed!
    date_obj = pd.to_datetime(date_str)
    if len(date_str) <= 11: ## if only date is provided --> get to end of the day
        date_obj += timedelta(hours=23, minutes=59, seconds=59)
    return date_obj  # return time
    
def get_first_last_day(outpath, patient):

    ## First day at hospital
    ## import hospital admission date 
    redcap_instance_timepoints = pd.read_csv(outpath + patient + "/" + patient + "_redcap_instance_timepoints.csv")
    if redcap_instance_timepoints['entry_type'].str.contains('basedata_hospadm').any():
        first_day = redcap_instance_timepoints.loc[redcap_instance_timepoints['entry_type'].str.contains('basedata_hospadm'), 'date']
    elif redcap_instance_timepoints['entry_type'].str.contains('enrollment_inclusiondate').any(): 
        first_day = redcap_instance_timepoints.loc[redcap_instance_timepoints['entry_type'].str.contains('enrollment_inclusiondate'), 'date']
        print(f'Using enrollment date instead of hospital admission date for patient {patient}')
    else:
        first_day = redcap_instance_timepoints.sort_values('date')['date'].head(1)
        print(f'Hospital admission date and enrollment date are missing. Using first date which was entered in database for patient {patient} instead.')
    first_day = first_day.reset_index(drop = True)
    first_day = first_day.astype('datetime64[ns]')[0]

    ## Last day at hospital
    ## end of study, including death and hospital discharge
    endpoint_df = pd.read_csv(outpath + patient + "/endpoint/" + patient + "_redcap_endpoint.csv")
    ## subset df to data of relevance 
    if 'endpoint_deathdate' in endpoint_df.columns:
        death = endpoint_df[['redcap_event_name', 'redcap_repeat_instance', 'endpoint_deathdate']]
        death = death[~death['endpoint_deathdate'].isna()]
        if not death.empty:
            death['date'] = death['endpoint_deathdate'].apply(to_end_of_day)
            death['condition'] = 'death'
            last_day_reason = death[['date', 'condition']].reset_index(drop = True)
    elif 'endpoint_pretermenddate' in endpoint_df.columns:
        study_pretermend = endpoint_df[['redcap_event_name', 'redcap_repeat_instance', 'endpoint_pretermenddate', 'endpoint_pretermexclusion']]
        study_pretermend = study_pretermend[~study_pretermend['endpoint_pretermenddate'].isna()]
        if not study_pretermend.empty:
            study_pretermend['date'] = study_pretermend['endpoint_pretermenddate'].apply(to_end_of_day)
            study_pretermend['condition'] = study_pretermend['endpoint_pretermexclusion'] 
            study_pretermend['condition'] = study_pretermend['condition'].str.replace('discharge', ' discharge') ## add whitespace 
            last_day_reason = study_pretermend[['date', 'condition']].reset_index(drop = True)
    elif 'endpoint_endpointdate' in endpoint_df.columns:
        study_end = endpoint_df[['redcap_event_name', 'redcap_repeat_instance', 'endpoint_endpointdate']]
        study_end = study_end[~study_end['endpoint_endpointdate'].isna()]
        if not study_end.empty:
            study_end['date'] = study_end['endpoint_endpointdate'].apply(to_end_of_day)
            study_end['condition'] = 'study end'
            last_day_reason = study_end[['date', 'condition']].reset_index(drop = True)
    else:
        last_day_reason = pd.DataFrame()
    return first_day, last_day_reason

## load isolated species
def load_isolated_species(isolate_metadata):

    ## check if any species for that patient have been reported 
    if not os.path.exists(os.path.expanduser(isolate_metadata)):
        return pd.DataFrame()
    
    isolate_df = pd.read_excel(isolate_metadata, skiprows=1)
    isolate_df = isolate_df[~isolate_df['is_resequenced_isolate']]
    isolate_df.loc[(isolate_df['hospitalization_days'] == 59.71) & (isolate_df['patient'] == 'P0007') & (isolate_df['focal_species'] == 'Bacteroides_thetaiotaomicron'), 'isolation_source'] = 'Pleura'
    ## group and summarize df
    isolate_df = isolate_df.groupby(['patient', 'focal_species', 'hospitalization_days', 'isolation_source'])['sampleID'].count().reset_index()
    isolate_df_wide = isolate_df.pivot_table(index = ['patient', 'focal_species', 'hospitalization_days'], columns = 'isolation_source', fill_value = 0).droplevel(0, axis=1).reset_index()
    
    return isolate_df_wide


def load_physical_cond(outpath, patient, first_day, no_species):

    ## read in ventilation and sedation/reanimation df
    try: 
        vent_df = pd.read_csv(outpath + patient + "/" + patient + "_redcap_ventilation.csv")
    except: 
        vent_df = pd.DataFrame()
    try: 
        sedation_df = pd.read_csv(outpath + patient + "/" + patient + "_redcap_patient_condition.csv")
    except:
        sedation_df = pd.DataFrame()

    ## merche mechanical treatment df together
    mech_treatment_df = pd.concat([vent_df,  sedation_df])
    if not mech_treatment_df.empty:
        ## ensure to have all columns which are needed 
        for col_needed in ['begin', 'stop', 'date', 'type', 'condition', 'measurement']:
            if col_needed not in mech_treatment_df.columns:
                mech_treatment_df[col_needed] = float('nan')
        mech_treatment_df.loc[mech_treatment_df['type'].isna(), 'type'] = mech_treatment_df.loc[mech_treatment_df['type'].isna(), 'condition']
        mech_treatment_df.dropna(subset = 'type', inplace = True)
        ## not efficient, but allows to avoid errors
        for col_to_drop in ['redcap_event_name', 'redcap_repeat_instance', 'condition', 'dur']: ## note duration of ventilation is dropped --> might be usable in later plots as it contains the duration of ventilation in hours on the latest timepoint
            if col_to_drop in mech_treatment_df.columns: 
                mech_treatment_df.drop(col_to_drop, axis = 1, inplace = True) 

        mech_treatment_df[['begin', 'stop', 'date']] = mech_treatment_df[['begin', 'stop', 'date']].astype('datetime64[ns]')
        mech_treatment_df['stop'] = (mech_treatment_df['stop'] - mech_treatment_df['begin']) / np.timedelta64(1, 'D')
        mech_treatment_df['begin'] = (mech_treatment_df['begin'] - first_day) / np.timedelta64(1, 'D')
        mech_treatment_df['date'] = (mech_treatment_df['date'] - first_day) / np.timedelta64(1, 'D')
        mech_treatment_df[['begin', 'stop', 'date']] = mech_treatment_df[['begin', 'stop', 'date']].round(3)
        
        mech_treatment_df.drop_duplicates(inplace = True)

        print("NOTE: ANALYSE all different types for dict")
        mech_treatment_df['marker'] = np.where(mech_treatment_df['measurement'] == 'trach', 'd', 
                                              np.where(mech_treatment_df['type'].str.contains('discharge|end'), 'o', 'X'))
        ## set colors
        ## Note this is still implemented to be quickly tweakable for other plots even though we just want intubation/extubation and heartfailure
        ## Note 'heart_failure' and death are set individually to colors --> Not stated in the list below!; study end gets same color as hospital discharge!
        unique_mechanical_treatments = ['covid19', 'sedation', 'PDT', 'CST', 'mechanical_ventilation', 'non_invasive_ventilation', 'nasal_high_flow', 'oxygenation', 'hospital discharge']
        cmap_decorations = plt.get_cmap('tab10')(np.linspace(0, 1, 10))[:len(unique_mechanical_treatments)] ## get all colors from colormap, use just first n for unique entries --> therefore not the entire cmap is used, but the first n colors instead!
        color_dict_decorations = dict(zip(unique_mechanical_treatments, cmap_decorations))
        color_dict_decorations['cardiac_arrest'] = [0.85098, 0.11764, 0.09411, 1] ## set to color red
        color_dict_decorations['death'] = [0, 0, 0, 1] ## set to color black
        mech_treatment_df['colors'] = mech_treatment_df['type'].map(color_dict_decorations).astype('object')
        if mech_treatment_df['type'].str.contains('study end').any():
            tmp_df = pd.DataFrame({'color_np_array': [color_dict_decorations['hospital discharge']]}) ## need to convert dict to pd.df to allow replacement!
            mech_treatment_df.loc[mech_treatment_df['type'] == 'study end', 'colors'] = tmp_df['color_np_array'] ## set study end to same color as hospital discharge 
        mech_treatment_df.loc[mech_treatment_df['type'] == 'PDT', 'label'] = 'PDT (tracheostomy)'
        mech_treatment_df.loc[mech_treatment_df['type'] == 'CST', 'label'] = 'CST (tracheostomy)'
        mech_treatment_df.loc[mech_treatment_df['type'] == 'begin mechanical_ventilation', 'label'] = 'intubation'
        mech_treatment_df.loc[mech_treatment_df['type'] == 'end mechanical_ventilation', 'label'] = 'extubation'
        mech_treatment_df.loc[mech_treatment_df['type'] == 'begin oxygenation', 'label'] = 'begin oxygenation (face mask)'
        mech_treatment_df['label'] = mech_treatment_df['type'].str.replace('_', ' ')
        mech_treatment_df['label'] = mech_treatment_df['type'].str[0].str.upper() + mech_treatment_df['label'].str[1:]

        mech_treatment_df = mech_treatment_df.sort_values('type')
        mech_treatment_df = mech_treatment_df[mech_treatment_df['label'] != 'study end']
        
    return mech_treatment_df


def get_ab_treatment_data(outpath, patient, hosp_admission_tp, unique_ab_class_color):
    ## process antiinfectives df 
    ab_treatment_df = pd.read_csv(outpath + patient + "/antiinfectives/" + patient + "_redcap_antiinfectives.csv")
    
    #####
    ## Drop empty rows and empty columns
    ab_treatment_df = ab_treatment_df[~ab_treatment_df.iloc[:, 3:].isin([0, np.nan]).all(axis = 1)]
    ab_treatment_df = ab_treatment_df.loc[:, (ab_treatment_df != 0).any(axis = 0)]

    ## stop function if emtpy df
    if ab_treatment_df.empty:
        return
    
    ## change df to long format, tweaek col names revert to wide format
    ab_treatment_df_long = ab_treatment_df.melt(id_vars = ['Unnamed: 0', 'redcap_event_name', 
                                                           'redcap_repeat_instance'], 
                                                var_name = "measurement")
    
    ab_treatment_df_long["measurement"] = ab_treatment_df_long["measurement"].str.replace("antiinfectives_", "")
    
    str_repl = ["start", "end", "dose", "tdm"]

    for str_match in str_repl:
        ab_treatment_df_long["measurement"] = ab_treatment_df_long["measurement"].str.replace(str_match, "_" + str_match)
    ab_treatment_df_long[["var1", "var2"]] = ab_treatment_df_long["measurement"].str.split("_", n=1, expand = True)
    ab_treatment_df_long["var2"] = ab_treatment_df_long["var2"].str.replace("__", "")
    ## makrolids have wrong pattern (makrop) for end points instead of makrol
    ab_treatment_df_long['var1'] = np.where(ab_treatment_df_long['var1'].str.contains('makrop'), 'makrol', ab_treatment_df_long['var1'])
    ab_treatment_df_long['var1'] = np.where(ab_treatment_df_long['var1'].str.contains('linco'), 'lincosamide', ab_treatment_df_long['var1'])

    ab_treatment_df_long['var1'],ab_treatment_df_long['var2'] = np.where(ab_treatment_df_long['var2'].isin(str_repl),
                                                                        (ab_treatment_df_long['var2'], ab_treatment_df_long['var1']),
                                                                        (ab_treatment_df_long['var1'], ab_treatment_df_long['var2']))
    ab_treatment_df_long = ab_treatment_df_long.drop(['Unnamed: 0', "measurement"], axis = 1)

    ## generate dict for abs which have been used at which timepoint
    ab_treatment_tps = ab_treatment_df_long[ab_treatment_df_long['var1'].isin(["antibiotic"]) & (ab_treatment_df_long['value'] == 1)] ## if antimycotica should be shown to, use also phrase "antimycotic"

    ## get names of antibiotic classes 
    ab_unique = ab_treatment_tps["var2"].unique()
    if len(ab_treatment_df_long[ab_treatment_df_long['var1'].isin(["eradication"])]['var2'].unique()) > 0:
        ab_unique = np.append(ab_unique, 'oralab')

    ab_df = pd.DataFrame(ab_unique, columns=(["ab_grp"]))
    
    ## stop function if no AB
    if ab_df.empty: 
        return

    ## oral AB treatment has pattern eradication for AB type (but might not be consistent!)
    ab_treatment_df_long['var1'] = np.where(ab_treatment_df_long['var1'] == 'eradication', 'oralabtype', ab_treatment_df_long['var1'])

    ## for some antibiotics, full name of type is provided in value but not in var2
    ab_treatment_df_long['var2'] = np.where((ab_treatment_df_long['var1'].str.contains('type') & ab_treatment_df_long['var2'].isnull()), 
                                            ab_treatment_df_long['value'], ab_treatment_df_long['var2']) ## if pattern 'type' in var1 and var2 is nan, store value in var2
    ab_treatment_df_long['value'] = np.where((ab_treatment_df_long['var1'].str.contains('type') & ab_treatment_df_long['var2'].isnull()), 
                                            0, np.where( (ab_treatment_df_long['var1'].str.contains('type') & (ab_treatment_df_long['value'] == ab_treatment_df_long['var2']) ), 
                                                        1,  ab_treatment_df_long['value'])) 
    
    ab_treatment_df_long = ab_treatment_df_long[~(ab_treatment_df_long['var1'].str.contains('type') & ab_treatment_df_long['var2'].isnull())]
    
    ab_types = ab_treatment_df_long[ab_treatment_df_long['var1'].str.contains("type")][["var1", "var2"]].drop_duplicates()
    ab_types['var2'] = ab_types['var2'].str.lower()
    ab_types["var1"] = ab_types["var1"].str.replace("type", "")
    ab_pat = '|'.join(ab_types['var1'])
    
    ## combine AB that zou have the real AB and not just the AB group as value 
    ab_df['subpat'] = ab_df['ab_grp'].str.extract('('+ ab_pat + ')', expand=False)
    if ab_df['ab_grp'].str.contains('cephalosporin').any():
        ab_df.loc[ab_df['ab_grp'].str.contains('cephalosporin'), ['subpat', 'ab_grp']] = ab_df.loc[ab_df['ab_grp'].str.contains('cephalosporin'), 'ab_grp'].str.replace('cephalosporin', 'ceph') ## keep number in subpattern
    ab_df = pd.merge(ab_df, ab_types, left_on= 'subpat', right_on='var1', how = 'left').drop(['subpat', 'var1'], axis = 1)
    ab_df['var2'] = ab_df['var2'].fillna(ab_df['ab_grp'])
    
    ## Subselect interesting columns only
    ab_treatment_df_long = ab_treatment_df_long.loc[ab_treatment_df_long['var1'].isin(str_repl)].reset_index(drop = True)
    
    ## Insert the full AB name 
    ab_treatment_subpat = '|'.join(ab_treatment_df_long['var2'].unique())
    ab_df['subpat'] = ab_df['var2'].str.extract('('+ ab_treatment_subpat + ')', expand=False)
    ab_df['subpat1'] = ab_df['ab_grp'].str.extract('('+ ab_treatment_subpat + ')', expand=False)
    ab_df['subpat'] = ab_df['subpat'].fillna(ab_df['subpat1'])
    ab_df['subpat'] = np.where(ab_df['ab_grp'] == 'oralab', 'oralab', ab_df['subpat']) ## change subpattern ofor oral ABs
    ab_df['ab_class'] = ab_df['ab_grp'].str.replace(r'^ceph', 'cephalosporin', regex = True)
    ab_df.loc[ab_df['ab_grp'] == 'clindamycin', 'ab_grp'] = 'lincosamide'

    ab_treatment_df_long[["ab", "ab_class"]] = np.nan

    for row in range(len(ab_df)):
        ab_treatment_df_long["ab"] = np.where(ab_treatment_df_long['var2'].isin(ab_df.iloc[row]), 
                                              ab_df.iloc[row]['var2'], 
                                              ab_treatment_df_long["ab"])
        ab_treatment_df_long["ab_class"] = np.where(ab_treatment_df_long['var2'].isin(ab_df.iloc[row]), 
                                              ab_df.iloc[row]['ab_class'], 
                                              ab_treatment_df_long["ab_class"])
    
    ab_treatment_df_long['ab_cleaned_class'] = ab_treatment_df_long['ab'].map(ab_class_dict)
    ## generate color map
    
    ab_treatment_df_long['colors'] = ab_treatment_df_long['ab_cleaned_class'].map(unique_ab_class_color).astype(str) ## conversion to string to allow pivot table on np.array
    
    ab_treatment_df_long = ab_treatment_df_long.drop('var2', axis = 1)
        
    ab_treatment_df_wide = ab_treatment_df_long.pivot_table(index = ['redcap_event_name', 
                                                               'redcap_repeat_instance', 'ab', 'ab_class', 'ab_cleaned_class', 'colors'],
                                                      columns = 'var1', aggfunc = 'first').droplevel(0, axis=1).reset_index()
    
    
    ## change names of antimycotics; note: AB's should be changed just before plotting or also in the diagnostic samples csv, otherwise mapping of resistancies to treated AB's does not work!
    antibiotic_name_change = {'lamb': 'amphotericin B\n(liposomal)', 
                                'amb': 'amphotericin B\n(non-liposomal)',
                                'cotrim': 'cotrimoxazol',
                                'unacid': 'ampicillin/sulbactam',
                                'tazobac': 'piperacillin/tazobactam'}
                                
    for key, value in antibiotic_name_change.items():
        ab_treatment_df_wide.loc[ab_treatment_df_wide['ab'] == key, 'ab'] = value
    
    ## Add missing columns
    ab_treatment_df_wide = ab_treatment_df_wide.reindex(ab_treatment_df_wide.columns.union(['start', 'end', 'dose', 'tdm'], 
                                                                                           sort=False), 
                                                        axis=1, fill_value= np.nan)
    ab_treatment_df_wide.dropna(subset = ['start', 'end', 'dose', 'tdm'], how = 'all', inplace = True)
    ab_treatment_df_wide[['start', 'end']] = (ab_treatment_df_wide[['start', 'end']].astype('datetime64[ns]') - hosp_admission_tp) / np.timedelta64(1, 'D')
    ## concat dose and tdm for instances where we don't have exact time resolution
    ab_treatment_df_wide = ab_treatment_df_wide.sort_values('redcap_repeat_instance').groupby(['start', 'ab', 'ab_class', 'ab_cleaned_class', 'colors']).agg({'end': 'last',
                                                                                                                                                              'dose': lambda x: ', '.join(x.dropna().astype(str)),
                                                                                                                                                              'tdm': lambda x: ', '.join(x.dropna().astype(str))}).reset_index()
    
    if ab_treatment_df_wide.size > 0:
        ab_treatment_df_wide['info'] = 'dose: ' + ab_treatment_df_wide['dose'] + '\n tdm: ' + ab_treatment_df_wide['tdm']
    else: 
        return
    ab_treatment_df_wide = ab_treatment_df_wide.drop(['dose', 'tdm'], axis = 1).sort_values('ab', ascending = False).reset_index(drop = True)
    color_is_list = ab_treatment_df_wide['colors'].str.contains('\[') & ab_treatment_df_wide['colors'].str.contains('\]')
    ab_treatment_df_wide.loc[color_is_list, 'colors'] = ab_treatment_df_wide.loc[color_is_list, 'colors'].apply(lambda x: np.fromstring(x.replace('\n','').replace('[','').replace(']','').replace('  ',' '), sep=' '))
    ## replace oral AB with class if possible 
    ab_treatment_df_wide.loc[ab_treatment_df_wide['ab_class'] == 'oralab', 'oral_hatch'] = '/'
    ab_treatment_df_wide.loc[ab_treatment_df_wide['ab_class'] != 'oralab', 'oral_hatch'] = ''
    oraltreated_abs = ab_treatment_df_wide.loc[ab_treatment_df_wide['ab_class'] == 'oralab', 'ab'].tolist()
    ab_class_oralabs = ab_treatment_df_wide.loc[(ab_treatment_df_wide['ab'].isin(oraltreated_abs)) & (ab_treatment_df_wide['ab_class'] != 'oralab'), ['ab', 'ab_class', 'colors']] ## extract ab_class of oral given ab which is NOT marked as oralab --> that AB is then given non-oral and therefore can be remapped to class

    for idx, oralab in ab_class_oralabs.iterrows():
        ab_treatment_df_wide.loc[ab_treatment_df_wide['ab'] == oralab['ab'], 'ab_class'] = oralab['ab_class']
            
    return ab_treatment_df_wide

def read_lineage_association(isolate_metadata, pairwise_snp_diff_file, lineage_snp_cutoff = 30, color_dict = False):
    isolate_df = pd.read_excel(isolate_metadata, skiprows=1)
    isolate_df = isolate_df[~isolate_df['is_resequenced_isolate']]
    
    pairwise_snp_diff = pd.read_csv(os.path.expanduser(pairwise_snp_diff_file))
    ## extract pathogen vs mbiome samples
    pathogen_species_list = isolate_df.loc[~isolate_df['isolation_source'].isin(['Rectal', 'Nasal', 'Oral', 'Skin', 'Gastric']), 'sampleID'].values ## note: exclude gastric joice samples!
    mbiome_species_list = isolate_df.loc[isolate_df['isolation_source'].isin(['Rectal', 'Nasal', 'Oral', 'Skin']), 'sampleID'].values ## note: exclude gastric joice samples!
    smpl1_p_smpl2_m = (pairwise_snp_diff['sample_1'].isin(pathogen_species_list) & pairwise_snp_diff['sample_2'].isin(mbiome_species_list) ) 
    smpl1_m_smpl2_p = (pairwise_snp_diff['sample_1'].isin(mbiome_species_list) & pairwise_snp_diff['sample_2'].isin(pathogen_species_list) ) 
    pairwise_snp_diff_pm = pairwise_snp_diff[smpl1_p_smpl2_m | smpl1_m_smpl2_p].reset_index(drop = True)
    pairwise_snp_diff_pm['pathogen'] = np.where(pairwise_snp_diff_pm['sample_1'].isin(pathogen_species_list), 
                                                pairwise_snp_diff_pm['sample_1'].values,
                                                pairwise_snp_diff_pm['sample_2'].values)
    pairwise_snp_diff_pm['mbiome'] = np.where(pairwise_snp_diff_pm['sample_1'].isin(mbiome_species_list), 
                                                pairwise_snp_diff_pm['sample_1'].values,
                                                pairwise_snp_diff_pm['sample_2'].values)
    pairwise_snp_diff_pm['patient'] = pairwise_snp_diff_pm['subject_species'].str.split('_').str[0]
    pairwise_snp_diff_pm['species_short'] = pairwise_snp_diff_pm['subject_species'].str.split('_').str[1].str.replace('YT3', '')
    pairwise_snp_diff_pm = pairwise_snp_diff_pm[['patient', 'species_short', 'pathogen', 'mbiome', 'snp_diff']]
    ################################################################
    ## NOTE: SHOULD JUST INCLUDE SAMPELS WHERE MBIOME SAMPLE IS BEFORE OR AT TP OF PATHOGEN SAMPLE!!!
    ## get closest pathogen sample to mbiome sample
    pairwise_snp_diff_pm = pairwise_snp_diff_pm.sort_values('snp_diff').groupby('mbiome').first().reset_index()
    snp_diff_to_pathogen = pd.merge(isolate_df, pairwise_snp_diff_pm, how = 'left', left_on = ['sampleID', 'patient'], right_on = ['mbiome', 'patient'])
    ## subset to mbiome samples of inetrest only 
    snp_diff_to_pathogen = snp_diff_to_pathogen[~snp_diff_to_pathogen['mbiome'].isna()].reset_index(drop = True)
    snp_diff_to_pathogen['pathogen_lineage'] = np.where(snp_diff_to_pathogen['snp_diff'] <= lineage_snp_cutoff, 
                                                        True, False)
    ## group by kit and select if pathogen lineage only, other lineage or both are at timepoint and sample location
    lineage_at_tp_loc = snp_diff_to_pathogen.groupby(['patient', 'species_short', 'hospitalization_days', 'isolation_source'])['pathogen_lineage'].agg(['sum', lambda x: len(x) - x.sum()]) 
    lineage_at_tp_loc = lineage_at_tp_loc.rename({'sum': 'true', '<lambda_0>': 'false'}, axis = 1).reset_index() ## true == num_patho_lin_isolates ; false == num_nonpatho_lin_isolates
    lineage_at_tp_loc['pathogen_lineage'] = np.where((lineage_at_tp_loc['true'] > 0), 'true', 'false')
    if color_dict:
        lineage_at_tp_loc['pathogen_lineage_color'] = lineage_at_tp_loc['pathogen_lineage'].map(color_dict)
    
    ## split in 2nd df to return --> lineage counts 
    lineage_patho_count = lineage_at_tp_loc.groupby(['patient', 'species_short'])[['true', 'false']].sum().reset_index() ## aggregate to single species and patient clusters
    lineage_patho_count = lineage_patho_count.melt(['patient', 'species_short'], var_name = 'pathogen_lineage', value_name = 'count') ## make long df and map color dict 
    if color_dict:
        lineage_patho_count['pathogen_lineage_color'] = lineage_patho_count['pathogen_lineage'].map(color_dict)
    lineage_at_tp_loc = lineage_at_tp_loc.drop(['true', 'false'], axis = 1)
    return lineage_at_tp_loc, lineage_patho_count

def clinical_sample_legend(legend_info, legend_labels):
    legend_info += [Line2D([], [], marker='None', linestyle='None')] ## header
    legend_info += [Line2D([0], [0], marker = 'X', color = 'k', markerfacecolor = clinic_cnt_cmap['Blood'][0], linestyle="None", label = 'Blood culture', markersize = 10)] ## blood
    legend_info += [Line2D([0], [0], marker = 'X', color = 'k', markerfacecolor = clinic_cnt_cmap['Pleura'][0], linestyle="None", label = 'Pleural fluid', markersize = 10)] ## pleural fluid
    legend_info += [Line2D([0], [0], marker = 'X', color = 'k', markerfacecolor = clinic_cnt_cmap['Lung'][0], linestyle="None", label = 'Respiratory tract\nmaterial', markersize = 10)] ## lung

    legend_info += [Line2D([], [], marker='None', linestyle='None')] ## space to fill to same length of entries per column
    legend_labels += ['Culture-positive\nclinical sample', 'Blood culture', 'Pleural fluid', 'Respiratory tract\nmaterial', '']
    return legend_info, legend_labels

######################################
## Control input
######################################


## Control for input 
if outpath[-1] != '/':
    outpath = outpath + '/'
if outpath[0] == '~':
    outpath = os.path.expanduser(outpath)

## generate dict to rename mislabbeled species per patient
rename_species_d = {'P0007': [['Enterobacter cloacae', 'Enterobacter hormaechei'],
                             ['Acinetobacter pitii', 'Acinetobacter pittii']],
                    'P0010': [['Klebsiella oxytoca', 'Klebsiella michiganensis']]}

lineage_snp_cutoff = 30 ## num snpds to be used for lineage cutoff

## Ab classes cleaned
ab_class_dict = {'cefotaxim': 'Cephalosporin',# '3rd gen. Cephalosporin', 
                'ceftazidim': 'Cephalosporin',# '3rd gen. Cephalosporin', 
                'tazobac': 'Penicillin', 
                'ciprofloxacin': 'Fluorochinolone', 
                'imipenem': 'Carbapenem', 
                'meropenem': 'Carbapenem', 
                'vancomycin': 'Glycopeptide', 
                'erythromycin': 'Macrolide', 
                'flucloxacillin': 'Penicillin', 
                'unacid': 'Penicillin', 
                #'linezolid': 'Oxazolidinone', 
                'cefuroxim': 'Cephalosporin', # '1st gen. Cephalosporin', 
                'levofloxacin': 'Fluorochinolone', 
                'cotrim': 'Antifolate', 
                'clarithromycin': 'Macrolide', 
                'azithromycin': 'Macrolide'}
ab_groups = {'Beta-lactam': ['Penicillin', 'Carbapenem', 'Cephalosporin'],# ['Penicillin', 'Carbapenem', '1st gen. Cephalosporin', '3rd gen. Cephalosporin'],
'Glycopeptide': ['Glycopeptide'], 
'Antifolate': ['Antifolate'],
'Fluorochinolone': ['Fluorochinolone'],
'Macrolide': ['Macrolide'],
'Oxazolidinone': ['Oxazolidinone']}
ab_targets = {'Cell wall synthesis': ['Penicillin', 'Carbapenem', 'Cephalosporin', 'Glycopeptide'], # ['Penicillin', 'Carbapenem', '1st gen. Cephalosporin', '3rd gen. Cephalosporin', 'Glycopeptide'],
'Protein synthesis': ['Macrolide'], #, 'Oxazolidinone'],
'Nucleic acid synthesis': ['Antifolate', 'Fluorochinolone']}

## sorted for 4x cell wall ABS, 2x protein targeting ABs, 2x nucelic acid ABs
tol_muted_palette = [
                    '#44BB99', 
                    '#AAAA00',
                    '#BBCC33',
                    '#EEDD88',
                    '#77AADD',
                    #'#99DDFF',
                    '#EE8866',
                    '#FFAABB',
                    '#DDDDDD']

unique_ab_class_color = {ab: color for color, ab in zip(tol_muted_palette, sum([ab_targets[ab_grp] for ab_grp in ab_targets.keys()], []))}

clinic_cnt_cmap = {'B': ["#ca0020"],
                   'L': ["#0571b0"],
                   'G': ["green"]}

clinic_cnt_cmap = {'Blood': ["#ca0020"],
                   'Pleura': ["#c97985"],
                   'Lung': ["#0571b0"],
                   'Gastric': ["green"]}


## get lineage associations per lineage
## generate grey scale 
color_list = colors.LinearSegmentedColormap.from_list("", ["black","white"])(np.linspace(0, 1, 4)) #plt.get_cmap('Greys')(np.linspace(1, 0, 5)) ## 
pathogen_color = '#eb4824'              ## color for diagnostic samples of pathogen species (== pathogen)
pathogen_lineage_color = color_list[0]      ## color for mbiome samples of pathogen lineage (== true)
nonpathogen_lineage_color = color_list[2]   ## color for mbiome samples of non-pathogen lineage (== false)
lineage_color_dict = {'true': pathogen_lineage_color, 'false': nonpathogen_lineage_color, 'pathogen': pathogen_color}
isolate_df = load_isolated_species(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/tables/2025_09_08_hai_isolates.xlsx')
lineage_association, lineage_assoc_count = read_lineage_association(isolate_metadata=f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/tables/2025_09_08_hai_isolates.xlsx',
                                                                    pairwise_snp_diff_file='~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/refbased_2024_08/analysis/2025_06_09_pairwise_SNP_diff.csv',
                                                                    lineage_snp_cutoff = lineage_snp_cutoff,
                                                                    color_dict = lineage_color_dict)

# Create figure and plot a stem plot with the date
fontsize = 9
linewidth = 0.25 

## alignment of the squares
xalignment = 1.05 ## for square size 50 and figsize (5,7): 1.25
yalignment = 0.17 ## for square size 50 and figsize (5,7): 0.23
routine_scatter_align = {'Nasal': [-xalignment, yalignment], ## upper left 
                         'Oral': [xalignment, yalignment], ## upper right
                         'Rectal': [-xalignment, -yalignment], ## lower left
                         'Skin': [xalignment, -yalignment]} ## lower right

adjust_species = 0.85 ## instead of 1 x distance between species, use adjust_species value 

patients_plt_order = ['7', '21', '10'] ## order patients by the length of hospital stay


#######################################
####################################### 
######################################
## PLOTTING
######################################

plt_offset = 0 ## offset to account for subplot offset of next patient to be plotted
max_xval = 0 ## will be automatically updated during plotting
ab_treatment_offset = 1 ## ofset for AB treatment annotation from timeline
ytick_label_dict = {} ## keep track of all species for ylabel setting
num_samples_in_patient_l = [] ## list of df for figure 1C
lineage_appearance_vs_inf_l = [] ## list of df for figure 1D
figheight = 5.75
if not plt_treatment:
    figheight *= 20/25
fig, ax = plt.subplots(figsize=(5.5, figheight))
                                       
plt.rc('font', family='Helvetica')

## plot AB treatment below other stuff

plt.setp(ax.get_xticklabels(), fontsize=fontsize)
plt.setp(ax.get_yticklabels(), fontsize=fontsize, verticalalignment = "center")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)


## mark the 48h mark for HAI definition

for plt_id, patient in enumerate(patients_plt_order):
    patient_id = 'P' + patient.zfill(4)

    first_day, last_day_reason_df = get_first_last_day(outpath, patient_id)
    
    ## load routine swabs info 
    routine_swabs_df = pd.read_csv(outpath + patient_id + "/" + patient_id + "_kits.csv", 
                                    sep = ",")
    
    ## remove second kit as it was not processed as it should have been (>> 4 h after it was taken put into -80°C)
    if patient_id == 'P0007':
        routine_swabs_df = routine_swabs_df[routine_swabs_df['sampleid'] != 'K00009-L']

    ## if empty df and no admission swab was taken, go to next patient
    if routine_swabs_df.empty:
        SystemError(f"Routine sample df is empty for patient {patient_id}")
        continue 
    routine_swabs_df = routine_swabs_df[['sampleid', 'date', 'date_d']]
    routine_swabs_df['sampleid'] = routine_swabs_df['sampleid'].str.replace('-L', '')
    
    ## subset isolates to patient of interest 
    isolate_df_pat = isolate_df[isolate_df['patient'] == patient_id]
    
    ## get pseudo last day
    if last_day_reason_df.empty:
        last_day_reason_df.loc[0, 'date'] = max(routine_swabs_df['date'].astype('datetime64[ns]'))
        print('Note: No actual last day was insterted. We use the last day of any sample taken from the patient as current last day!')

    hosp_stay_length = (last_day_reason_df['date'] - first_day)[0].days
    if max(isolate_df_pat['hospitalization_days']) > hosp_stay_length:
        hosp_stay_length = max(isolate_df_pat['hospitalization_days'])

    ## translate species to numeric values for plotting
    unique_species = isolate_df_pat['focal_species'].dropna().sort_values(ascending = False).unique()
    spec_to_num = {spec: idx*adjust_species for idx, spec in enumerate(unique_species)}
    spec_short_to_num = {spec[0] + spec.split('_')[1]: idx*adjust_species for spec, idx in spec_to_num.items()} ## make a short species name for lineage association (we just have the short name there!)
    spec_short_label_to_num = {f'{spec[0]}. {spec.split("_")[1]}' if (spec[-4:]!=' sp.') else spec: idx*adjust_species for idx, spec in enumerate(unique_species)} ## we need that for putting y label later on plot

    plt_offset += (len(unique_species)*adjust_species+1.25)
    if plt_treatment:
        plt_offset += ab_treatment_offset

    isolate_df_pat['spec_id'] = isolate_df_pat['focal_species'].map(spec_to_num)
    num_samples_in_patient = lineage_assoc_count[lineage_assoc_count['patient'] == patient_id]
    num_diagnostic_samples_in_patient = isolate_df_pat.groupby(['focal_species'])[['Blood', 'Pleura', 'Lung', 'Gastric']].sum().sum(axis = 1).reset_index() ## NOTE: We count for figure S3 gastric samples as clinical sample!
    num_diagnostic_samples_in_patient = num_diagnostic_samples_in_patient.rename({0: 'count'}, axis = 1)
    num_diagnostic_samples_in_patient['species_short'] = num_diagnostic_samples_in_patient['focal_species'].str[0] + num_diagnostic_samples_in_patient['focal_species'].str.split('_').str[1]
    num_diagnostic_samples_in_patient[['patient', 'pathogen_lineage']] = [patient_id, 'pathogen']
    num_diagnostic_samples_in_patient['pathogen_lineage_color'] = num_diagnostic_samples_in_patient['pathogen_lineage'].map(lineage_color_dict)
    num_diagnostic_samples_in_patient = num_diagnostic_samples_in_patient[num_diagnostic_samples_in_patient.columns[1:]] ## drop spec_clust column
    if num_samples_in_patient.empty:
        num_samples_in_patient = num_diagnostic_samples_in_patient.reset_index(drop = True)
    else:
        num_samples_in_patient = pd.concat([num_samples_in_patient, num_diagnostic_samples_in_patient]).reset_index(drop = True)
    num_samples_in_patient = num_samples_in_patient.fillna(value = 0)
    num_samples_in_patient['spec_id'] = num_samples_in_patient['species_short'].map(spec_short_to_num) - plt_offset
    if any(num_samples_in_patient['spec_id'].isna()):
        print('NOTE: Some species associations for the barplot have failed --> recheck!')
    num_samples_in_patient_l.append(num_samples_in_patient)
    lineage_appearance_vs_inf_l.append(num_samples_in_patient)


    ## set locations of plots:
    timeline_ypos = len(unique_species)*adjust_species - plt_offset
    if plt_treatment:
        timeline_ypos += 0.75
    subplot_locs = {"Timeline": timeline_ypos}
    
    ## get clinical_course data
    mechanical_treatment_df = load_physical_cond(outpath, patient_id, first_day, no_species = subplot_locs["Timeline"])

    ## load treatment pattern of patient
    ab_treatment_df = get_ab_treatment_data(outpath, patient_id, first_day, unique_ab_class_color)


    ####################################### 
    ## plot timeline of Pat and the clinical course
    if patient_id == 'P0021': ## account for unsampled stay of Patient 21 
        ax.plot([0, hosp_stay_length],[subplot_locs["Timeline"]] *2, ls='-',lw=1.5,color='k')
    else:
        ax.plot([0, hosp_stay_length],[subplot_locs["Timeline"]] *2,'o',ls='-',lw=1.5,color='k',ms=5, markevery=[-1])
    
    ## set the vertical lines for sampling timepoints
    [ax.plot([date, date], 
             [subplot_locs["Timeline"] - 0.25, subplot_locs["Timeline"] + 0.25], 
             color = 'k', 
             lw = 1.5) 
        for date in routine_swabs_df['date_d'].values]
    
    
    ## add for ventilation an tickened line above the timeline
    if not mechanical_treatment_df.empty:
        print('NOTE: need to check if criteria are ok for selection of timeperiods!!!')
        intubation_periods = mechanical_treatment_df.loc[~mechanical_treatment_df['begin'].isna() & 
                                                         ~mechanical_treatment_df['stop'].isna() &
                                                         (mechanical_treatment_df['label'].str.contains('ventilation') | mechanical_treatment_df['label'].str.contains('oxygenation') ), 
                                                         ['begin', 'stop']].values
        if patient_id == 'P0021':
            intubation_periods = [[start, end] if end < hosp_stay_length else [start, hosp_stay_length]for start, end in intubation_periods]
        for period in intubation_periods:
            ax.plot([period[0], period[1]],[subplot_locs["Timeline"]] *2,ls='-',lw=3,color='k')

    if plt_treatment:
        ## add AB treatment below timeline
        for ab_idx, ab in enumerate(sorted(ab_treatment_df['ab_cleaned_class'].unique())):
            ab_treatment_curr_ab = ab_treatment_df[ab_treatment_df['ab_cleaned_class'] == ab]
            intubation_periods = ab_treatment_curr_ab.loc[~ab_treatment_curr_ab['start'].isna() & 
                                                            ~ab_treatment_curr_ab['end'].isna(), 
                                                            ['start', 'end', 'colors']].values
            ab_offset_plt = 0.125*ab_idx
            if patient_id == 'P0021':
                intubation_periods = [[start, end, color] if end < hosp_stay_length else [start, hosp_stay_length, color] for start, end, color in intubation_periods]
            for period in intubation_periods:
                ax.plot([period[0], period[1]],[subplot_locs["Timeline"]-ab_treatment_offset + ab_offset_plt] *2,ls='-',lw=2,color=period[2])

    ####################################### 
    ## plot isolates collected from routine samples 
    if not isolate_df_pat.empty and plt_mbiome: ## if any sample has been collected for patient

        for sample_type, xy_adj in routine_scatter_align.items():
            ## get colors based on lineage association
            tmp_plt_df = isolate_df_pat[isolate_df_pat[routine_scatter_align.keys()].sum(axis = 1) > 0].copy() ## select only time stamps with collected mbiome samples
            tmp_plt_df['species_short'] = tmp_plt_df['focal_species'].str[0] + tmp_plt_df['focal_species'].str.split('_').str[1]
            tmp_lineage_association = lineage_association[(lineage_association['isolation_source'] == sample_type) & (lineage_association['patient'] == patient_id)]
            tmp_plt_df = pd.merge(tmp_plt_df, tmp_lineage_association, how = 'left', on = ['species_short', 'hospitalization_days'])
            tmp_plt_df.loc[tmp_plt_df['pathogen_lineage_color'].isna(), 'pathogen_lineage_color'] = '#ffffff' ## set all to white if no isolate was identified!

            ax.scatter(tmp_plt_df['hospitalization_days'] + xy_adj[0], 
                            tmp_plt_df['spec_id'] + (xy_adj[1]) - plt_offset,
                            color = tmp_plt_df['pathogen_lineage_color'],
                            edgecolors = 'k',
                            lw = linewidth,
                            zorder = 10,
                            s = 45, 
                            marker = 's')
    ####################################### 
    ## plot isolates collected from clinical samples 
    for sample_type, sample_color in clinic_cnt_cmap.items():
        if sample_type not in isolate_df_pat.columns:
            continue 
        tmp_clin_subset_df = isolate_df_pat[isolate_df_pat[sample_type] > 0]
        
        if sample_type != 'Gastric':
            markersize = 40    
            ax.scatter(tmp_clin_subset_df['hospitalization_days'], 
                            tmp_clin_subset_df['spec_id'] - plt_offset,
                            color = 'w',
                            edgecolors = 'w',
                            lw = 0.75,
                            zorder = 20,
                            alpha = 0.9,
                            s = markersize+10, 
                            marker = 'X')
            ax.scatter(tmp_clin_subset_df['hospitalization_days'], 
                            tmp_clin_subset_df['spec_id'] - plt_offset,
                            color = sample_color,
                            edgecolors = 'k',
                            lw = linewidth,
                            zorder = 20,
                            alpha = 0.9,
                            s = markersize, 
                            marker = 'X')
    
    ## annotate subparts of plot
    for x, y in zip([-35, -35, -35], subplot_locs.values()): ## if full species names are used use [-60, -60, -60] instead of [-35, -35, -35]
        ax.annotate(f'Patient {patient}', 
                    xy=(x, y), xytext=(0,0),
                    textcoords="offset points",
                    verticalalignment='bottom', 
                    horizontalalignment = 'left',
                    annotation_clip=False,
                    fontsize = fontsize + 2)

    
    ## store information for next loop rounds
    if max_xval < max(routine_swabs_df['date_d']): ## get max value of xlabel for tick setting
        max_xval = max(routine_swabs_df['date_d'])
    if max_xval < mechanical_treatment_df[['begin', 'date', 'stop']].max().max():
        max_xval = mechanical_treatment_df[['begin', 'date', 'stop']].max().max()
    ##NOTE implement mechanical treatment as well!
    if plt_treatment:
        ytick_label_dict[subplot_locs['Timeline'] - ab_treatment_offset + 0.2] = 'Antibiotic treatment'
    for species, locid in spec_short_label_to_num.items():
        ytick_label_dict[locid - plt_offset] = species
    

## set legend
if plt_legend:
    
    legend_labels = []
    legend_columns = 1
    legend_info = [Line2D([], [], marker='None', linestyle='None')] ## header timeline
    legend_info += [Line2D([0], [0], marker = '_', color = 'k', markerfacecolor = 'k', linestyle="None", markeredgewidth = 1.5, markersize = 20)] ## ICU stay
    legend_info += [Line2D([0], [0], marker = '_', color = 'k', markerfacecolor = 'k', linestyle="None", markeredgewidth = 3, markersize = 20)] ## intubation
    legend_info += [Line2D([0], [0], marker = '|', color = 'k', markerfacecolor = 'k', linestyle="None", markeredgewidth = 1.5, markersize = 10)] ## Microbiome sampling timepoint
    legend_info += [Line2D([0], [0], marker = 'o', color = 'k', markerfacecolor = 'k', linestyle="None", markersize = 5)] ## discharge label
    legend_info += [Line2D([], [], marker='None', linestyle='None')] ## space to legend group before
    legend_labels += ['Timeline', 'ICU stay', 'Mechanical\nventilation', 'Microbiome\nsampling', 'End of study', ''] 
    if plt_treatment:
        legend_info += [Line2D([], [], marker='None', linestyle='None')] ## header
        legend_labels += ['Antibiotic classes']
        for ab, ab_color in unique_ab_class_color.items():
            legend_info += [Line2D([0], [0], marker = '_', color = ab_color, markerfacecolor = ab_color, linestyle="None", markeredgewidth = 3, markersize = 20)] ## 
            if 'Cephalosporin' in ab:
                ab = ab.replace('gen. ', 'gen.\n')
            legend_labels += [ab]
    else:
        legend_info, legend_labels = clinical_sample_legend(legend_info, legend_labels)
    
    

    legend_columns += 1
    legend_info += [Line2D([], [], marker='None', linestyle='None')] ## header
    legend_info += [Line2D([0], [0], marker = 's', color = 'k', markerfacecolor = 'w', linestyle="None", markersize = 20)] ## routine sample, no label yet
    legend_info += [Line2D([], [], marker='None', linestyle='None')] ## space to legend group before
    legend_info += [Line2D([], [], marker='None', linestyle='None')] ## header
    marker_shape = 's'
    legend_info += [Line2D([0], [0], marker = marker_shape, color = 'k', markerfacecolor = lineage_color_dict['true'], linestyle="None", label = 'Pathogenic', markersize = 10)] ## lineage type
    legend_info += [Line2D([0], [0], marker = marker_shape, color = 'k', markerfacecolor = lineage_color_dict['false'], linestyle="None", label = 'Non-pathogenic', markersize = 10)] ## lineage type
    legend_info += [Line2D([0], [0], marker = marker_shape, color = 'k', markerfacecolor = '#ffffff', linestyle="None", label = 'None', markersize = 10)] ## lineage type
    legend_info += [Line2D([], [], marker='None', linestyle='None')] ## space to legend group before
    legend_labels += ['Microbiome niche', 'Nasal\tOral\nRectal\tSkin', '', 'Distance to closest\nclinical isolate', f'< {lineage_snp_cutoff} SNV', f'≥ {lineage_snp_cutoff} SNV', 'None', '']
    
    if plt_treatment:
        legend_info, legend_labels = clinical_sample_legend(legend_info, legend_labels)

ax.xaxis.tick_top()
ax.tick_params(width=0, pad=0)
ax.xaxis.set_label_position('top')
xtick_width = 7
ax.set_xticks(np.arange(max_xval + xtick_width-1, step = xtick_width)) ## xtick_width-1 to account for max value which would be one below next xtick step 
ax.set_xticklabels(np.arange(max_xval + xtick_width-1, step = xtick_width, dtype = 'int8'))
ax.set_yticks(list(ytick_label_dict.keys()))
ax.set_yticklabels(list(ytick_label_dict.values()))
ax.xaxis.grid(True, lw = 0.2, zorder = 0)
ax.yaxis.grid(True, ls = ':', lw = 0.2, zorder = 0)
xlabel_text = 'Days since hospital admission'
ax.set_xlabel(xlabel_text, size = fontsize)

if plt_legend:
    plt.legend(legend_info, legend_labels, 
                ncol=legend_columns, loc='center left', 
                bbox_to_anchor=(1, 0.5), fancybox=False, 
                shadow=False, columnspacing = 0.2, labelspacing = 0.3, 
                fontsize = fontsize) # Two columns, vertical group labels


# plt.tight_layout()
filename = f'Panel1B_pilot_study_patients_timelines_250609_v5{figure_suffix}'

if saveme:
    plt.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/{filename}.pdf'),
                bbox_inches = 'tight', transparent=True)
    plt.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/{filename}.svg'),
                bbox_inches = 'tight', transparent=True, format='svg')


########################################
########################################
########################################
## FIGURE S4
########################################
########################################
########################################

def split_label(label):
    parts = label.split("__")
    if len(parts) == 2:
        return parts[1], parts[0]  # Return (primary_label, secondary_label)
    return label, "" 

def nested_xlabels(ax):
    prim_sec_labels = {}
    ax_labels = ax.get_xticklabels()
    for label in ax_labels:
        primary, secondary = split_label(label.get_text())
        if secondary in prim_sec_labels.keys():
            prim_sec_labels[secondary].append(primary)
        else:
            prim_sec_labels[secondary] = [primary]
    ## set primary label
    ax.set_xticks(np.arange(len(ax_labels)))
    prim_labels = [label for labels in prim_sec_labels.values() for label in labels]
    ax.set_xticklabels(prim_labels, rotation = 90)
    ## set secondary label 
    sec_pos = []
    prev_count = [-0.5]
    for prim_label in prim_sec_labels.values():
        sec_pos.append(prev_count[-1] + 0.5*len(prim_label))
        sec_lab = [label[0] + label[3:] for label in prim_sec_labels.keys()]
        prev_count.append(prev_count[-1] + len(prim_label))
    sec = ax.secondary_xaxis(location=-0.65)
    sec.set_xticks(sec_pos, labels=sec_lab)
    sec.tick_params('x', length=0, width = 0)
    sec.spines[['bottom']].set_visible(False)
    sec2 = ax.secondary_xaxis(location=0)
    sec2.set_xticks(prev_count, labels=[])
    sec2.tick_params('x', length=115, width=0.5)


num_samples_in_patient_all = pd.concat(num_samples_in_patient_l).reset_index(drop = True)
num_samples_in_patient_all['patient_spec'] = num_samples_in_patient_all['patient'] + '__' + num_samples_in_patient_all['species_short'].str[0] + '. ' + num_samples_in_patient_all['species_short'].str[1:]
num_samples_in_patient_all = num_samples_in_patient_all.sort_values('patient_spec').reset_index(drop = True)
num_samples_in_patient_all['pathogen_lineage'] = pd.Categorical(num_samples_in_patient_all['pathogen_lineage'], ['pathogen', 'true', 'false'], ordered = True)


fig, ax = plt.subplots(figsize = (7.5, 5))
## plot on second axis the bar plots for isolate counts
sns.barplot(data = num_samples_in_patient_all, x = 'patient_spec', y = 'count', hue = 'pathogen_lineage', palette = lineage_color_dict, edgecolor = 'k', lw = 0.5, width = 0.65, ax = ax, zorder = 10)
ax.set_yscale('symlog')
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.get_major_formatter().set_scientific(False)
ax.yaxis.get_major_formatter().set_useOffset(False)
ax.set_ylim((0, ax.yaxis.get_view_interval()[1]*1.5))
max_val_axis = int(np.ceil(np.log10(ax.yaxis.get_view_interval()[1]))) ## the max log base of the axis to set minor ticks in that range
ax.yaxis.set_minor_locator(FixedLocator([k * 10**n for n in range(max_val_axis) for k in range(1, 10)])) ## set minor ticks for all integer values (remove minor ticks < 1!)
ax.grid(axis = 'y', linewidth = 0.25, c = 'k')
ax.grid(axis = 'y', which = 'minor', linewidth = 0.05, c = 'k')
legend_info = [mpatches.Patch(facecolor = lineage_color_dict['pathogen'], edgecolor = 'k', lw = 0.5)] ## lineage type 
legend_info += [mpatches.Patch(facecolor = lineage_color_dict['true'], edgecolor = 'k', lw = 0.5)] ## lineage type
legend_info += [mpatches.Patch(facecolor = lineage_color_dict['false'], edgecolor = 'k', lw = 0.5)] ## lineage type
legend_labels = ['Clinical isolate', f'< {lineage_snp_cutoff} SNV', f'≥ {lineage_snp_cutoff} SNV']
ax.legend(title = 'Distance to closest\nclinical isolates', handles = legend_info, labels = legend_labels, 
          loc='center left', bbox_to_anchor=(1, 0.5), fancybox=False, shadow=False, fontsize = fontsize)
ax.set_xlabel(None)
ax.set_ylabel('Number of isolates')
nested_xlabels(ax)
plt.tight_layout()

filename = f'Panel1C_pilot_study_num_isolates_barplot'


if saveme:
    plt.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/{filename}.pdf'),
                bbox_inches = 'tight', transparent=True)
    plt.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/{filename}.svg'),
                bbox_inches = 'tight', transparent=True, format='svg')
