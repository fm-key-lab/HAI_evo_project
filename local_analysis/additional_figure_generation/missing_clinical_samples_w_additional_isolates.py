
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
import seaborn as sns
import numpy as np

plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

redcap_path = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/RedCap_Data/Proc/')
## get the fixed list 
corrected_medmibi_df = pd.read_csv('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/RedCap_Data/2024_12_18_HAIpilot_infections_with_missing_entries_wadditional_medmibi_isolates.csv', skiprows = 1)
corrected_medmibi_df = corrected_medmibi_df[[ col for col in corrected_medmibi_df.columns[1:] if not col.startswith('Unnamed')]]
corrected_medmibi_df = corrected_medmibi_df[~corrected_medmibi_df['patient'].isna() & ~corrected_medmibi_df['date processing'].isna()]
for col in ['date sampling', 'date processing']:
    corrected_medmibi_df[col] = corrected_medmibi_df[col].astype(str).str.split(':').str[:2].str.join(':') ## cut seconds away for proper transformation
    corrected_medmibi_df[col] = pd.to_datetime(corrected_medmibi_df[col], format='%d.%m.%Y %H:%M')
corrected_medmibi_df['species'] = corrected_medmibi_df['species'].str.replace(' ', '_')
corrected_medmibi_df = corrected_medmibi_df.rename({'date sampling': 'date', 'Redcap sampleid': 'sampleid'}, axis = 1)
corrected_medmibi_df.loc[corrected_medmibi_df['Coments'].isna(), 'Coments'] = ''
is_contamination = (corrected_medmibi_df['Coments'].str.lower().str.contains('ontamination') | corrected_medmibi_df['Source of infection'].str.lower().str.contains('ontamination'))
corrected_medmibi_df = corrected_medmibi_df[~is_contamination] ## remove all contaminations
processing_has_date = (corrected_medmibi_df['date'].isna() & ~corrected_medmibi_df['date processing'].isna())
corrected_medmibi_df.loc[processing_has_date, 'date'] = corrected_medmibi_df.loc[processing_has_date, 'date processing']
corrected_medmibi_df['species'] = corrected_medmibi_df['species'].str.replace('ESBL_E.coli', 'Escherichia_coli_ESBL')

patients = ['P0007', 'P0010', 'P0021']
patients_nonhai_species = {'P0007': ['Staphylococcus_epidermidis', 'Serratia_marcescens'], 
                           'P0010': [], 
                           'P0021': ['Staphylococcus_epidermidis']}
patients_rename_species = {'P0007': {'Enterobacter_cloacae': 'Enterobacter_hormaechei'}, 
                           'P0010': {'Klebsiella_oxytoca': 'Klebsiella_michiganensis'}, 
                           'P0021': {}}
target_hai_material_of_study = ['Urinary tract infection', 'UTI',
                                'Pneumonia', 'HAP', 'VAP',
                                'blood stream infection', 'BSI']
non_hai_material = ['Gastric aspirate', ## self checked in lab and not considered as infection itself
                    'Wound smear', ## not in our IRB as focus
                    'Routine swab']
diagnostic_smpls_d = {}
plt_heights = []

for patient in patients:
    patient_path = f'{redcap_path}{patient}'

    diagnostic_file = f'{patient_path}/{patient}_diagnostic_smpls.csv'
    redcap_timepoint_map_file = f'{patient_path}/{patient}_redcap_instance_timepoints.csv'
    
    redcap_tps = pd.read_csv(redcap_timepoint_map_file)
    for date in ['date']:
        redcap_tps[date] = redcap_tps[date].astype('datetime64[ns]')
    if not redcap_tps.empty:
        firstday = min(redcap_tps['date'])
    else:
        firstday = 0

    ## read diagnostic sample file in 
    if os.path.exists(diagnostic_file):
        diagnostic_smpls = pd.read_csv(diagnostic_file) 
        diagnostic_smpls = diagnostic_smpls[[col for col in diagnostic_smpls if not col.startswith('resistence')]]
        diagnostic_smpls['date'] = diagnostic_smpls['date'].astype('datetime64[ns]')
        ## subset to identified bacterial HAI species only!
        diagnostic_smpls = diagnostic_smpls[~diagnostic_smpls['species'].isna()]
        diagnostic_smpls['species'] = diagnostic_smpls['species'].str.replace(' ', '_')
        if patient == 'P0007':
            diagnostic_smpls.loc[diagnostic_smpls['sampleid'] == 'B00096', 'species'] = 'Bacteroides_thetaiotaomicron' ## correct Pleuraempyem (contextual sample for BSI caused by B theta)
        
        ## add the additional isolates not entered in redcap, but from where we got isolates
        if patient == 'P0021':
            add_pathos = diagnostic_smpls.loc[diagnostic_smpls['sampleid'].isin(['L00020', 'L00026'])]
            add_pathos['cfu'] = np.nan
            add_pathos['spec_id'] += 1
            add_pathos.loc[add_pathos['sampleid'] == 'L00020', 'species'] = 'Proteus_mirabilis'
            add_pathos.loc[add_pathos['sampleid'] == 'L00026', 'species'] = 'Escherichia_coli'
            diagnostic_smpls = pd.concat([diagnostic_smpls, add_pathos]).sort_values('date').reset_index(drop = True)


        diagnostic_smpls = pd.merge(diagnostic_smpls, corrected_medmibi_df[corrected_medmibi_df['patient'] == patient], how = 'outer', on = ['date', 'species', 'sampleid'])
        diagnostic_smpls['missing'] = np.where((diagnostic_smpls['sampleid']=='-') | (diagnostic_smpls['Coments']=='nicht am MPI angekommen'), 'missing', 'available')
        
        ## deduplicate samples
        diagnostic_smpls = diagnostic_smpls.sort_values(['missing', 'date']).groupby(['date', 'species']).head(1)

        ## remove nonhai species and material
        for ms_spec_name, wgs_spec_name in patients_rename_species[patient].items():
            diagnostic_smpls['species'] = diagnostic_smpls['species'].str.replace(ms_spec_name, wgs_spec_name) ## rename to properly identified species
        
        ## get only target samples (BSI/HAP/VAP/UTI/Gastritis) and get all associated samples!
        source_is_target = diagnostic_smpls['Source of infection'].fillna('').str.lower().str.contains('|'.join(target_hai_material_of_study).lower(), regex = True)
        diagnostic_smpls = diagnostic_smpls[source_is_target | (diagnostic_smpls['sampleid'].str.len() == 6)] ## get only species identified in target locus and its associated (contextual) samples

        diagnostic_smpls = diagnostic_smpls[~diagnostic_smpls['species'].str.contains('Standortflora') & ~diagnostic_smpls['species'].isin(patients_nonhai_species[patient])] 
        diagnostic_smpls = diagnostic_smpls[~diagnostic_smpls['Sample material'].isin(non_hai_material)]
        diagnostic_smpls.loc[diagnostic_smpls['Sample material'].str.lower().str.contains('pleural'), 'Sample material'] = 'Pleural fluid'

        uniq_spec = sorted(diagnostic_smpls['species'].unique())
        diagnostic_smpls['days'] = (diagnostic_smpls['date'] - firstday).dt.total_seconds() / (24*60*60) ## convert to days  
        ## restrict to max days of study (90d!)
        diagnostic_smpls = diagnostic_smpls[diagnostic_smpls['days'] <= 90]

        diagnostic_smpls['species'] = diagnostic_smpls['species'].str.replace('_', ' ')
        diagnostic_smpls['freq_avail'] = diagnostic_smpls.groupby('species')['missing'].transform(lambda x: (x == 'available').sum()/len(x))
        general_freq = round((diagnostic_smpls['missing'] == 'available').sum() / len(diagnostic_smpls), 1)*100
        diagnostic_smpls['species_short'] = diagnostic_smpls['species']
        has_species_name = diagnostic_smpls['species_short'].str.split(' ').str.len() > 1
        diagnostic_smpls.loc[has_species_name, 'species_short'] = diagnostic_smpls.loc[has_species_name, 'species_short'].str[0] + '. ' + diagnostic_smpls.loc[has_species_name, 'species_short'].str.split(' ').str[1]
        diagnostic_smpls['label'] = diagnostic_smpls['species_short'] + ' (' + (diagnostic_smpls['freq_avail']*100).round(0).astype(int).astype(str)+'%)'
        

    ## get number of samples (not isolates!) which have been collected/missed
    uniq_smpls = diagnostic_smpls.sort_values('missing').groupby('Sample ID MedMIBI').head(1)
    patients_overall_freq = int( round( (sum(uniq_smpls['missing'] == 'available') / len(uniq_smpls))*100 , 0))
    diagnostic_smpls_d[patient] = [diagnostic_smpls, patients_overall_freq]
    plt_heights.append(len(uniq_spec))

total_samples = sum([len(diagnostic_smpls_d[patient][0].groupby('Sample ID MedMIBI').head(1)) for patient in diagnostic_smpls_d.keys()])
total_samples_avail = sum([sum(diagnostic_smpls_d[patient][0].groupby('Sample ID MedMIBI').head(1)['missing'] == 'available') for patient in diagnostic_smpls_d.keys()])

############################
############################
############################
max_xlim = 4

material_to_color_d = {}

color_to_material_d = {'#ca0020': ['Blood culture', 'Pleural fluid'], ## pleura just included as contextual sample
                       '#0571b0': ['Bronchial lavage', 'Tracheal aspirate']
}


for base_hex, items in color_to_material_d.items():
    num_items = len(items)
    base_rgb = mcolors.to_rgb(base_hex)  # RGB in 0–1 range
    h, s, v = mcolors.rgb_to_hsv(base_rgb)
    if base_hex == '#0571b0':
        saturations = np.linspace(s/4, s, num_items)[::-1]  # Prevent oversaturation
    else:
        saturations = np.linspace(s/2.5, s, num_items)[::-1]  # Prevent oversaturation
    if num_items == 1:
        saturations = np.array([s])
    for item, new_s in zip(items, saturations):
        material_to_color_d[item] = mcolors.to_hex(mcolors.hsv_to_rgb([h, new_s, v]))


## Plot
fig, axs = plt.subplots(figsize=(7, 3.75), nrows = 3, sharex=True, gridspec_kw = {'hspace': 0.3, 'height_ratios': plt_heights}) 
axs = axs.flatten()
for ax, patient in zip(axs, patients):
    plt_df = diagnostic_smpls_d[patient][0].sort_values('date').reset_index(drop = True)
    plt_df['label'] = pd.Categorical(plt_df['label'], sorted(plt_df['label'].unique(), reverse = True), ordered = True)
    
    sns.scatterplot(data = plt_df[plt_df['missing']=='available'], x = 'days', y = 'label', hue = 'Sample material', palette = material_to_color_d, marker = 'X', s=50, edgecolor = 'k', lw = 0.5, alpha = 0.8, legend = False, zorder = 10, ax = ax)
    sns.scatterplot(data = plt_df[plt_df['missing']=='missing'], x = 'days', y = 'label', hue = 'Sample material', palette = material_to_color_d, marker = 'o', s=25, edgecolor = 'k', lw = 0.5, alpha = 0.8, legend = False, zorder = 15, ax = ax)
    ax.axvline(2, color='grey', ls = '--', lw = 0.5, zorder = 0) ## mark HAI cutoff
    ax.set_ylim(-0.5, len(np.unique(plt_df['species']))-0.5)
    ax.margins(y = 0.1)
    ax.set_ylabel(None)
    ax.xaxis.set_major_locator(MultipleLocator(7))
    ax.xaxis.grid(True, ls = ':', lw = 0.4, zorder = 0)
    min_xlim = min(0, plt_df['days'].min()-1)
    max_xlim = max(max_xlim, plt_df['days'].max())+1
    ax.set_xlim(min_xlim, max_xlim)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title(f'{patient} ({diagnostic_smpls_d[patient][1]}% of clinical samples available)', loc = 'right', va = 'top', fontsize = 10)

h = [Line2D([], [], marker='None', linestyle='None', label = 'Patient material')] ## 
for l, c in material_to_color_d.items():
    h += [Line2D([], [], marker='s', linestyle='None', markersize=10, color=c, markeredgecolor='black', label=l)]
h += [Line2D([], [], marker='None', linestyle='None', label = 'Sample availability\nas study material')] ## 
h += [Line2D([], [], marker='X', linestyle='None', markersize=10, color='w', lw = 0.5, markeredgecolor='black', label='Available')]
h += [Line2D([], [], marker='o', linestyle='None', markersize=6, color='w', lw = 0.5, markeredgecolor='black', label='Unavailable')]
axs[1].legend(handles = h, loc = 'center left', bbox_to_anchor = (1, 0.5))
axs[-1].set_xlabel('Stay at ICU [days]')
fig.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2025_07_09_HAI_pilot_missing_clinical_samples.pdf'),
                bbox_inches = 'tight', transparent=True)
fig.savefig(os.path.expanduser(f'~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2025_07_09_HAI_pilot_missing_clinical_samples.svg'),
                bbox_inches = 'tight', transparent=True)

## print the total number and frequency
print(f'{total_samples_avail}/{total_samples} ({total_samples_avail/total_samples*100}%) samples are available across patients')