
## make frequency chart of species infecting before, during or undetected a patient

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import seaborn as sns
import pandas as pd
import numpy as np

plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

## plotting functions 
def strain_before_or_at_symptome_onset(grouped_df, sites = ['nasal', 'oral', 'rectal', 'skin', 'anysite']):
    my_data = []
    for groupname, group in grouped_df:
        patient = groupname[0]
        species = groupname[1]
        for site in sites:
            if group.loc[(group['variable'] == f'Mbiome_detected_{site}'), 'value'].values[0] & group.loc[(group['variable'] == f'before_or_during_infection_{site}'), 'value'].values[0]:
                if group.loc[group['variable'] == f'same_strain_as_infection_{site}', 'value'].values[0]:
                    status = f'< {snp_cutoff} SNVs'
                else:
                    status = f'≥ {snp_cutoff} SNVs'
            else:
                status = 'None'
            my_data.append([patient, species, site, status])
    my_df = pd.DataFrame(my_data, columns = ['patient', 'species', 'site', 'status'])
    return my_df

def split_label(label):
    parts = label.split("__")
    if len(parts) == 2:
        return parts[1], parts[0]  # Return (primary_label, secondary_label)
    return label, "" 

def nested_labels(ax, label_loc = 1):
    prim_sec_labels = {}
    ax_labels = ax.get_yticklabels()
    for label in ax_labels:
        primary, secondary = split_label(label.get_text())
        if secondary in prim_sec_labels.keys():
            prim_sec_labels[secondary].append(primary)
        else:
            prim_sec_labels[secondary] = [primary]
    ## set primary label
    ax.set_yticks(np.arange(len(ax_labels))+0.5)
    prim_labels = [label for labels in prim_sec_labels.values() for label in labels]
    ax.set_yticklabels(prim_labels)
    ## set secondary label 
    sec_pos = []
    prev_count = [0]
    for prim_label in prim_sec_labels.values():
        sec_pos.append(prev_count[-1] + 0.5*len(prim_label))
        sec_lab = [label for label in prim_sec_labels.keys()]
        prev_count.append(prev_count[-1] + len(prim_label))
    sec = ax.secondary_yaxis(location=-1.3*label_loc)
    sec.set_yticks(sec_pos, labels=sec_lab)
    sec.tick_params('y', length=0, width = 0)
    sec.spines[['left']].set_visible(False)
    sec2 = ax.secondary_yaxis(location=0)
    sec2.set_yticks(prev_count, labels=[])
    sec2.tick_params('y', length=120, width=0.5)
    
## Load data
strains_detected = pd.read_excel(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/metadata/2024_08_18_pilot_study_species_detected_site_strain.xlsx'))
## remove any species not identified by clinic
strains_detected_clin = strains_detected[strains_detected['identified_by_clinic'] == True]

## get the three different classes: 1. Same strain in Mbiome detected, 2. Different strain in Mbiome detected, 3. Not in Mbiome detected
plt_dict = {}
## plot data in a combined stacked and dodged barplot
site = 'anysite'
mbiome_col = f'Mbiome_detected_{site}'
before_inf_col = f'before_or_during_infection_{site}'
same_strain_col = f'same_strain_as_infection_{site}'
snp_cutoff = 30

color_list = mcolors.LinearSegmentedColormap.from_list("", ["black","white"])(np.linspace(0, 1, 4))
pathogen_lineage_color = color_list[0]      ## color for mbiome samples of pathogen lineage (== true)
nonpathogen_lineage_color = color_list[2]   ## color for mbiome samples of non-pathogen lineage (== false)
none_identified = color_list[3]
lineage_color_dict = {f'< {snp_cutoff} SNVs': pathogen_lineage_color, f'≥ {snp_cutoff} SNVs': nonpathogen_lineage_color, 'None': none_identified}
presence_absence_color_dict = {f'Observed': pathogen_lineage_color, 'Unobserved': none_identified}

strains_detected_clin_l = strains_detected_clin.melt(['Patient', 'Species'])
strains_detected_clin_l = strain_before_or_at_symptome_onset(strains_detected_clin_l.groupby(['Patient', 'Species']))
strains_detected_clin_l['patient_spec'] = strains_detected_clin_l['patient'] + '__' + strains_detected_clin_l['species']
strains_detected_clin_l_freq = (strains_detected_clin_l.groupby(['site', 'status'])['patient_spec'].count() / len(strains_detected_clin_l['patient_spec'].unique())).reset_index()
strains_detected_clin_l_freq_hm = strains_detected_clin_l_freq.pivot(index = 'status', columns = 'site', values = 'patient_spec').fillna(0)

strains_detected_clin_w = strains_detected_clin_l.pivot(index = ['patient_spec'], columns = 'site', values = 'status')
heatmap_numeric = strains_detected_clin_w.applymap(lambda x: list(lineage_color_dict.keys()).index(x))
heatmap_numeric.index = [f'{idx[:6]}. {idx.split(" ")[1]}' for idx in heatmap_numeric.index]

## sort heatmap and df for barplot the same way
mbiome_order_dict = {mbiome: idx for idx, mbiome in enumerate(['Nasal', 'Oral', 'Rectal', 'Skin', 'Microbiome'])}

heatmap_numeric = heatmap_numeric[['nasal', 'oral', 'rectal', 'skin', 'anysite']]
strains_detected_clin_l_freq_hm = strains_detected_clin_l_freq_hm[['nasal', 'oral', 'rectal', 'skin', 'anysite']]
heatmap_numeric.columns = mbiome_order_dict.values()
strains_detected_clin_l_freq_hm.columns = mbiome_order_dict.values()


snv_to_pathogen_presence = {'< 30 SNVs': 'Observed', 
                            'None': 'Unobserved', 
                            '≥ 30 SNVs': 'Unobserved'}


species_observed_mbiome = strains_detected.copy()
species_observed_mbiome.index = species_observed_mbiome['Patient'] + '__' + species_observed_mbiome['Species']
species_observed_mbiome = list(species_observed_mbiome[species_observed_mbiome['Mbiome_detected_anysite'] & species_observed_mbiome['identified_by_clinic']].index)
heatmap_numeric_obs = heatmap_numeric.copy()
cmap_culture_pos = mcolors.ListedColormap(presence_absence_color_dict.values())
heatmap_numeric = strains_detected_clin_w.applymap(lambda x: list(lineage_color_dict.keys()).index(x))

strains_detected_clin_l_obs = strains_detected_clin_l.copy()
strains_detected_clin_l_obs['patient_spec_short'] = strains_detected_clin_l_obs['patient'] + '__' + strains_detected_clin_l_obs['species']
strains_detected_clin_l_obs['patient_spec_short'] = [f'{idx[:6]}. {idx.split(" ")[1]}' for idx in strains_detected_clin_l_obs['patient_spec_short']]
strains_detected_clin_l_obs['status'] = strains_detected_clin_l_obs['status'].map(snv_to_pathogen_presence)
strains_detected_clin_l_obs = (strains_detected_clin_l_obs.groupby(['site', 'status'])['patient_spec_short'].count() / len(strains_detected_clin_l_obs['patient_spec_short'].unique())).reset_index()
strains_detected_clin_l_obs_hm = strains_detected_clin_l_obs.pivot(index = 'status', columns = 'site', values = 'patient_spec_short').fillna(0)

strains_detected_clin_l_obs_hm = strains_detected_clin_l_obs_hm[['nasal', 'oral', 'rectal', 'skin', 'anysite']]
heatmap_numeric_obs.columns = mbiome_order_dict.values()
strains_detected_clin_l_obs_hm.columns = mbiome_order_dict.values()


fig, (ax1, ax2) = plt.subplots(figsize = (2, 5.75), nrows = 2, sharex = True, gridspec_kw = {'height_ratios': [5, len(species_observed_mbiome)], 'hspace':0.05})
bottom_tracker = np.zeros((5), dtype = 'float')
for dist in presence_absence_color_dict.keys():
    tmp_df = strains_detected_clin_l_obs_hm[strains_detected_clin_l_obs_hm.index == dist]
    ax1.bar(x = tmp_df.columns+0.5, height = tmp_df.values[0], bottom=bottom_tracker, color = presence_absence_color_dict[dist], edgecolor = 'k', linewidth = 0.5)
    bottom_tracker += tmp_df.values[0]
ax1.set_ylim(0, 1)
ax1.tick_params(axis='x',which='major',bottom=False,labelbottom=False)
ax1.set_yticks(np.arange(0, 1.05, 0.2))
ax1.set_yticklabels(np.arange(0, 105, 20))
ax1.set_ylabel('Frequency\n of observaion [%]')
sns.heatmap(heatmap_numeric_obs, cmap=cmap_culture_pos, cbar=False, linecolor = 'k', linewidths=0.5, ax = ax2)
ax2.set_ylabel(None)
ax2.set_xlabel(None)
ax2.axvline(x = 4, lw = 0.5, c = 'k', clip_on=False)
ax2.set_xticks([val+0.5 for val in mbiome_order_dict.values()])
ax2.set_xticklabels(list(mbiome_order_dict.keys()), rotation = 90)
plt.tight_layout()
sec2 = ax2.secondary_xaxis(location=0)
sec2.set_xticks([4], labels=[])
sec2.tick_params('x', length=70, width=0.5)
for _, spine in ax2.spines.items():
    spine.set_visible(True)
nested_labels(ax2, label_loc = 0.6)

## generate legend
legend_info = [Line2D([], [], marker='None', linestyle='None')] ## header
legend_labels = ['Pathogen lineage presence\nbefore or within 6 h\nof HAI onset']
for dist, color in presence_absence_color_dict.items():
    legend_info += [Line2D([0], [0], marker = 's', color = 'k', markerfacecolor = color, label = '', linewidth=0, markersize = 12.5)] ## blood
    legend_labels += [dist]
ax2.legend(legend_info, legend_labels, loc = 'upper left', bbox_to_anchor = (1, 1))
fig.savefig(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2025_06_09_HAI_patients_patho_strain_presence_in_Mbiome_plt.pdf'),bbox_inches = 'tight', transparent=True)
fig.savefig(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/2025_06_09_HAI_patients_patho_strain_presence_in_Mbiome_plt.svg'),bbox_inches = 'tight', transparent=True)

