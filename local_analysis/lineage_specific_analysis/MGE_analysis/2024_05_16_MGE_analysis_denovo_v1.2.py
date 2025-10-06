
## analysis of MGE 

## NOTE:
## MGEs which are spread across non-consectutetive contigs are stated as separate MGEs! --> implementation of clustering (based on median, stdev and coverage correlation would be nice!)

## NOTE:
## MGE start and end positions as well as coordinates in final csv file are 0 based!

## NOTE:
# conda env 'anapy_312' is required

################
#%% IMPORT PACKAGES
################
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.sparse # inmporting sparse matrix
import gzip, pickle, sys, subprocess
import pandas as pd 
from Bio import Phylo
import os
import glob

# %% Import APY/MGE modules
SCRIPTS_DIRECTORY = os.getcwd() + "/modules/"
sys.path.insert(0, SCRIPTS_DIRECTORY)

import analysispy_module as apy
import mobileelements_module as mge

################
# %% PARAMETERs 
################

general_params_dict = {'scaling_factor': 1000, ## scaling factor used in snakemake to store floats as integer
                        'min_length': 750, # 750 ## min length of region to be identified
                        'overlap_bp': 10, ## max distance between regions individual mge regions to be called as one mge
                        'avg_cov2be_present': 3, ## mean coverage (absolute) across region to be identified as present
                        'group_nts_in_plot': 25, ## number of nucleotides to group in coverage plot
                        'use_ml_tree': True ## set true to use latest raxml tree for sorting isolates
                        }

## thresholds deletions 
deletion_params_dict = {'loose_threshold_del': -0.5 * general_params_dict['scaling_factor'],
                        'stringent_threshold_del': -1 * general_params_dict['scaling_factor'],
                        'threshold_raw_cov': 5 ## at least one sample's cov must be above this threshold in candidate deletion region
                        }

amplification_params_dict = {'loose_threshold_amp': 0.5 * general_params_dict['scaling_factor'],
                             'stringent_threshold_amp': 1 * general_params_dict['scaling_factor'],
                             'threshold_amp_mulitplier': 2.0, ## NOTE: amplifications will be determined by multiplier * mean_coverage(covered position) to remove non covered position which increase false positiove hits
                             'threshold_amp_multiplier_std': 0, ## multiplier for including the standard deviation to remove false positive hits in high variable coverage samples; set to 0 if this cutoff should not be used
                             'threshold_cov_step_at_region_boundaries': 1.25, ## coverage increase/drop at start/end of region (1/100th of region will be considered to evaluate this); value should be smaller or equal to 'threshold_amp_mulitplier'; set to 0 if this cutoff should not be used
                             'length_around_region_for_cov_step': 1000, ## length before/after region to calculate coverage in and compare it with coverage at very start/end of region --> step expected
                             'min_length_in_region_for_cov_step': 10, ## min length at start/end of region to calculate coverage in and compare it with coverage before/after region --> step expected
                             'fraction_of_region_for_cov_step': 0.10 ## fraction of region to calculate mean cov in region and compare it to cov argound region; if below min_length_in_region_for_cov_step, then min_length_in_region_for_cov_step will be taken
                             }

save_plot_params_dict = {'save_plots_of_all_positions': True, ## if set to true, plots of all merged mobile elements will be generated and saved
                         'save_plots_of_amp_only': False, ## if True, only plots of amplifications will be saved
                         'save_plots_of_del_only': True ## if True, only plots of deletions will be saved
                         }

dataset_for_interaction = [] ## list of reference genomes/datasets set by user to stop loop for interactive plotting of candidate MGEs until user presses enter to resume loop; set to 'all' to stop at every dataset 

################
# %% PATHS
################

analysis_run = 'denovo_2024_08'

main_workingdir = os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/' + analysis_run + '/analysis/backup/')
os.chdir(os.path.expanduser(main_workingdir))

cand_mut_table_dirs = '/../../case/2-candidate_mutation_table_2024_10_18'

## get species and their reference genomes (stored in subdir name) + 
## remove hidden files from list
refgenomes_ls = [dir for dir in os.listdir(os.getcwd() + cand_mut_table_dirs) if not dir.startswith('.')]
refgenomes_ls = list(sorted(refgenomes_ls)) ## sort refgenomes 
reference_path = '~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/metadata/reference_genomes/'

if 'all' in dataset_for_interaction:
    dataset_for_interaction = refgenomes_ls ## update if script should step at every ref_genome

for refgenome in refgenomes_ls:
    print(f'\n\n\n################################################\nAnalysing {refgenome}...')
    
    refgenome_path_full = os.path.expanduser(reference_path + refgenome)

    analysis_path = ''
    if general_params_dict['use_ml_tree']:
        analysis_paths = glob.glob(f'{main_workingdir}*{refgenome}/*/*_ML*.raxml.bestTreeCollapsed')
        if len(analysis_paths) > 0:
            analysis_path = sorted(analysis_paths, key = lambda x: x.rsplit('/')[-2])[-1] ## extract latest analysis
            analysis_path = '/'.join(analysis_path.split('/')[:-1]) + '/'
    if not general_params_dict['use_ml_tree'] or (analysis_path == ''):
        analysis_paths = glob.glob(f'{main_workingdir}*{refgenome}/*/')
        analysis_paths = [path for path in analysis_paths if not 'fimZ_' in path.split('/')[-2]] ## exlcude unwanted paths 
        analysis_path = sorted(analysis_paths, key = lambda x: x.rsplit('/')[-2])[-1] ## extract latest analysis

    if general_params_dict['use_ml_tree']:
        outpath = f'{analysis_path}mge_analysis_ml_tree_based'
    else:
        outpath = f'{analysis_path}mge_analysis'
    subprocess.run(f'mkdir -p {outpath}', shell = True)
    print(f'Using {outpath} for analysis.')

    ################
    # %% Save set parameters
    mge.save_set_params(f'{outpath}/{refgenome}', general_params_dict, deletion_params_dict, amplification_params_dict, save_plot_params_dict, is_interactive_analysis=(refgenome in dataset_for_interaction))

    ################
    # %% Get genome stats 
    chrstarts,genome_length,scaf_names = apy.genomestats(refgenome_path_full)

    ################
    # %% read sample Names from candidate_mutation_table (same order as cov_matrix) & get good samples from analysis
    with gzip.open(f'{main_workingdir}{cand_mut_table_dirs}/{refgenome}/candidate_mutation_table.pickle.gz', 'rb') as f:
        cmt = pickle.load(f)
    sampleNames = np.array(cmt[0])
    del cmt

    ## get order of samples in tree (== good nsample list as well)

    # read high-qual goodSamples from latest tree (analysis.py)
    # read in tip_label order of desired tree
    tree_paths = []
    if general_params_dict['use_ml_tree']:
        tree_paths = glob.glob(f'{analysis_path}*_ML.raxml.bestTreeCollapsed')
    if len(tree_paths) == 0:
        tree_paths = glob.glob(f'{analysis_path}{refgenome}_latest_rerooted_*_TP0.nwk.tree')
    if len(tree_paths) > 0: ## should just be one instance !
        if len(tree_paths) > 1:
            print('Mulitple trees have been identified. Please check!')
        mytree=sorted(tree_paths, key = lambda x: x.rsplit('/')[-1])[-1] ## extract latest tree
        tree = Phylo.read(mytree, "newick")
        tip_labels = np.array([i.name for i in tree.get_terminals()]) # tip labels ordered as tree
        
        ## get timepoints of samples (second instance in tip labels separated by '_')   
        tip_labels_ingroup = np.array([label for label in tip_labels if not ('_ref' in label) and not label.startswith('OG_')]) ## remove ref and outgroup samples
        tip_labels_tp_dates = np.array([f"{tp}_{date[6:8]}/{date[4:6]}/{date[2:4]}" for tp, date in (item.split('_')[1:3] for item in tip_labels_ingroup)])
        unique_tps_sorted = np.sort(np.unique(tip_labels_tp_dates))
        unique_tp_dict = {tp_date.split('_')[0]: tp_date.split('_')[1] for tp_date in unique_tps_sorted}
    else:
        print('No tree file was found. proceeding with the sampleNames as tiplabels')
        tip_labels = np.sort(np.array(sampleNames))
        unique_tp_dict = {}

    ## extract sample ids from tip labels (from tree)
    goodsamples_analysis = np.array([mge.extract_and_rename_sampleIDs_from_tree(sample) for sample in tip_labels])

    sampleNames = np.array([mge.rename_outgroups_sampleNames(sample) for sample in sampleNames])


    # %% LOAD cov data & filter samplesNames-goodpos-noOut & recalculate dbl-norm
    ################
    ## Data/Sim cov in format: rows: samples;  columns: positions
    cov_mat_csr = scipy.sparse.load_npz(f'{main_workingdir}{cand_mut_table_dirs}/{refgenome}/cov_raw_sparsecsr_mat.npz')
    cov_matrix = np.array(cov_mat_csr.todense())
    del cov_mat_csr

    cov_norm_csr = scipy.sparse.load_npz(f'{main_workingdir}{cand_mut_table_dirs}/{refgenome}/cov_norm_sparsecsr_mat.npz')
    cov_norm = np.array(cov_norm_csr.todense())
    del cov_norm_csr

    print('Coverage matrix read')

    ################
    # %% Filter samples
    ################
    # filter sampleNames for desired samples, and filter cov_matrix
    # NOTE: always keep sampleNames and cov_matrix row 'labels' in sync!
    bool_samples_filt = mge.names_of_a_found_in_b(sampleNames,goodsamples_analysis)
    sampleNames = sampleNames[ bool_samples_filt ]
    cov_matrix = cov_matrix[bool_samples_filt,]
    cov_norm = cov_norm[bool_samples_filt,]

    print('Filtering done')

    ################
    # %% Order matrices based on phylogenetic order 
    ################

    # read in tip_label order of desired tree
    [tip_label_indexCovMat,tip_label_absent] = mge.get_tip_label_order_for_cov_matrix(goodsamples_analysis,sampleNames)
    tip_label_indexCovMat_present = tip_label_indexCovMat[~tip_label_absent] ## geneate idx of tip labels present
    ## extract for later where to insert empty rows to insert empty data
    idx_to_insert_data = np.sort(np.where(tip_label_absent)[0])
    idx_to_insert_data = np.array([idx-i for i, idx in enumerate(idx_to_insert_data)]) ## remove on every idx one, as the np.insert used for those idx will not update idx before inserting data! (e.g. idx to insert: [2, 4] should [old_data, old_data, inserted_data, old_data, inserted_data])

    ## reorder matrices based on pyhlogenetic order
    cov_matrix_ordered = cov_matrix[tip_label_indexCovMat_present, :]
    del cov_matrix
    cov_norm_ordered = cov_norm[tip_label_indexCovMat_present, :]
    del cov_norm
    sampleNames_ordered = sampleNames[tip_label_indexCovMat_present]

    print('Ordered matrices by phylogeny')

    ################
    # %% Calculate mean of samples on coverer positions
    ################
    ## calculate mean and standard deviation of samples excluding all positions which are not covered
    ## This allows to remove effectively sites (and therefore in the end false positive amplifications) from samples with many deletions
    mean_coverage_of_covered_pos, std_coverage_of_covered_pos = mge.mean_std_cov_above_coverage(cov_matrix_ordered, cov_threshold = 0)
    mean_std_covered_pos_samples = np.concatenate((mean_coverage_of_covered_pos.reshape(-1, 1), std_coverage_of_covered_pos.reshape(-1, 1)),axis=1)

    ################
    # %% Detect Candidates
    ################
    obs_del = mge.candidate_deletions(cov_norm_ordered, deletion_params_dict['loose_threshold_del'] , deletion_params_dict['stringent_threshold_del'] , general_params_dict['min_length'] , cov_matrix_ordered , deletion_params_dict['threshold_raw_cov'])
    obs_amp = mge.candidate_amplifications(cov_norm_ordered, amplification_params_dict['loose_threshold_amp'] , amplification_params_dict['stringent_threshold_amp'] , general_params_dict['min_length'] , cov_matrix_ordered , 
                                            amplification_params_dict['threshold_amp_mulitplier'], amplification_params_dict['threshold_amp_multiplier_std'], 
                                            mean_std_covered_pos_samples, 
                                            amplification_params_dict['threshold_cov_step_at_region_boundaries'], amplification_params_dict['fraction_of_region_for_cov_step'], 
                                            amplification_params_dict['length_around_region_for_cov_step'], amplification_params_dict['min_length_in_region_for_cov_step']) ## use mean_coverage_of_covered_pos to remove false positive callings from samples which map poorly

    cand_amp_del = np.concatenate((obs_del,obs_amp))

    print('Identified candidate deletions and amplifications')

    ################
    # %% Explore candidate amplifications/deletions. Interactive plotting
    ################
    # Report includes:
    # - clickable overview of detected candidate mobile elements
    # - simple click can provide following info about individual candidates:
    #   - coverage (absolute/norm) for target (+-bumper zone) region 
    #   - report gene content of target region
    #   - plot coverage heatmap in order of tree (NOTE: reshaping tree often not poperly stored!) 


    ## generate reversed samplenames for heatmap plotting and correct extraction of sample ids
    heatmap_samplenames = [f'{samplename}_cov{int(mean_coverage_of_covered_pos[sampleNames == samplename][0])}' for samplename in sampleNames_ordered]

    # %% Evaluate candidates - Prerequisites
    # Generate variables necessary for interactive plotting
    line_bumper = 2000; # to make lines wide enough to be clickable
    [lines, lines_colors, n_genome, n_samples,list_samples,list_starts,list_ends,list_type] = mge.interactive_plotting_prerequisite(cov_matrix_ordered,cand_amp_del,line_bumper)

    if lines == []:
        print('No candidate MGEs identified.')
        with open(f'{outpath}/{refgenome}_no_cand_MGE_identified.txt', 'w') as fo:
            fo.write('')
        continue 

    # wrap variables and arguments for clickable graphic
    my_data_dict = {}
    my_data_dict['list_samples']=list_samples
    my_data_dict['list_starts']=list_starts
    my_data_dict['list_ends']=list_ends
    my_data_dict['lines_colors']=lines_colors
    my_data_dict['list_type']=list_type
    # Coverage matrices
    my_data_dict['cov_norm']=cov_norm_ordered
    my_data_dict['cov_matrix']=cov_matrix_ordered
    my_data_dict['mean_std_cov_samples']=mean_std_covered_pos_samples
    # Info on samples and genome
    my_data_dict['n_samples']=n_samples
    my_data_dict['n_genome']=n_genome
    my_data_dict['SampleNames']=sampleNames_ordered
    my_data_dict['reffolder']=refgenome_path_full # required for gene content report.
    # Info on samples in tree vs samples in coverage matrix
    my_data_dict['tip_label_indexCovMat'] = tip_label_indexCovMat_present # sample order as in tree 
    my_data_dict['tip_label_absent'] = tip_label_absent # samples in tree but not in data
    my_data_dict['idx_to_insert_data'] = idx_to_insert_data # indices of samples in tree which are not in cov matrix (or have been excluded) to insert data based on pyhlogeny 
    # Clickable plot options
    my_data_dict['report_genes']=True # Report genes present in mobile_element candidate region. Requires reffolder!
    my_data_dict['tree_aligned_cov_plot'] = True # Plot coverage distribution per base for each sample ordered as *latest.nwk.tree.
    my_data_dict['bool_plot_relative_cov_tree_aligned'] = True # plot relative coverage increase/decrease in amp/del respectively. Tree aligned.
    my_data_dict['plot2file'] = None # flag to print cov to plot.png 


    # How to make a plot of candidate amp/del segments
    # enable plotting in separate windows for the interactive plotting 
    if refgenome in dataset_for_interaction:
        %matplotlib qt 
    fig = mge.fig_candidate_segments( lines, lines_colors, n_genome, n_samples,my_data_dict, sampleNames=heatmap_samplenames, chrstarts=chrstarts)
    ## This will stop the loop at specific user defined datasets and just proceeds the script if all plotting windows were closed
    # or "enter" was clicked on the main plotting window created by fig_candidate_segments()
    if refgenome in dataset_for_interaction:
        fig.canvas.mpl_connect('key_press_event', mge.on_key)
        fig.savefig(f'{outpath}/{refgenome}_candidate_amplifications_deletions.pdf')
        plt.show(block = True)
        print('Script proceeds... ')
    else:
        fig.savefig(f'{outpath}/{refgenome}_candidate_amplifications_deletions.pdf')
        plt.close()


    #############plt.close('all')
    # %% Merge Candidate Regions
    #############
    # merge candidate regions based on overlap
    # use csv for taking notes about visual exploration
    annotation_genes = apy.parse_gff(refgenome_path_full,scaf_names,forceReDo=False)
    contig_edges = np.append(chrstarts, np.array([n_genome]))

    me_merged = mge.cand_mobile_element_region_merger(cand_amp_del,general_params_dict['overlap_bp'])
    me_merged = np.unique(me_merged[:,3:7],axis=0)

    ## read snp_table to identify snps falling into mge region
    if os.path.exists(f'{analysis_path}snp_table.csv'):
        snp_table = mge.read_snp_table(f'{analysis_path}snp_table.csv')    
    else:
        snp_table = pd.DataFrame()

    ## read indel_table to identify snps falling into mge region
    if os.path.exists(f'{analysis_path}indel_table.csv'):
        indel_table = mge.read_indel_table(f'{analysis_path}indel_table.csv')    
    else:
        indel_table = pd.DataFrame()

    ## save mobile elements detected with metadata 
    me_element_dict_list = []
    me_mean_cov_dict_list = []
    me_breadth_cov_dict_list = []
    matplotlib.use('agg')
    for i in range(me_merged.shape[0]):
        me_element_dict, isolates_w_ampdel_region_bool, me_mean_cov_dict, me_breadth_cov_dict = mge.extract_metadata_of_mge(me_merged[i, :], cov_matrix_ordered, cand_amp_del, mean_coverage_of_covered_pos, annotation_genes, sampleNames_ordered, chrstarts, snp_table, indel_table)
        me_element_dict_list.append( me_element_dict )
        me_mean_cov_dict_list.append( me_mean_cov_dict )
        me_breadth_cov_dict_list.append( me_breadth_cov_dict )
        
        ## if true, for all positions within the merged mobile elements, a plot will be generated
        if save_plot_params_dict['save_plots_of_all_positions'] or save_plot_params_dict['save_plots_of_amp_only'] or save_plot_params_dict['save_plots_of_del_only']:
            if not save_plot_params_dict['save_plots_of_all_positions'] and save_plot_params_dict['save_plots_of_amp_only'] and (me_element_dict['TYPE'] == 0):
                continue 
            if not save_plot_params_dict['save_plots_of_all_positions'] and save_plot_params_dict['save_plots_of_del_only'] and (me_element_dict['TYPE'] == 1):
                continue 
            if i % 10 == 0:
                print('.', end = '')
            cand_amp_del_in_region = cand_amp_del[isolates_w_ampdel_region_bool, 0]
            segment_ids_of_mge = np.where(isolates_w_ampdel_region_bool)[0]
            mge_name = f"{me_element_dict['TYPE']}_{str(me_element_dict['START'])}_{str(me_element_dict['END'])}"
            ## set for each grouped MGE a new subfolder 
            outpath_of_mge = f'{outpath}/{mge_name}'
            subprocess.run(f'mkdir -p {outpath_of_mge}', shell = True)
            
            plot_region_start, plot_region_end = mge.generate_plot_region(me_element_dict['START'], me_element_dict['END'], me_element_dict['LENGTH'], contig_edges[-1])
            cov_of_region_matrix = cov_matrix_ordered[:,me_element_dict['START']:me_element_dict['END']]
            mean_cov_of_regio_mat = np.mean(cov_of_region_matrix,axis=1)

            ## downsample the matrices to reduce plotting overhead 
            cov_of_region_matrix_ext = mge.downsample_mean_axis_1(cov_matrix_ordered[:,plot_region_start:plot_region_end], factor = general_params_dict['group_nts_in_plot'])
            cov_of_region_matrix_norm_ext = mge.downsample_mean_axis_1(cov_norm_ordered[:,plot_region_start:plot_region_end], factor = general_params_dict['group_nts_in_plot'])
            
            for sample_idx, segment_id in zip(cand_amp_del_in_region, segment_ids_of_mge):
                cov_plot_title = f'Sample ID: {sampleNames_ordered[sample_idx]},\nGenomic position: {me_element_dict["START"]}:{me_element_dict["END"]}, Segment length: {me_element_dict["LENGTH"]}'
                color_dict = {'cand_mge_edges_color': (0, 1, 0, 1), 
                            'highlight_sample_color': lines_colors[segment_id],
                            'samples_color': (0, 0, 0, .25), 
                            'contig_edges_color': (0.6, 0.2, 0.8, .5)}
                mge.create_norm_abs_cov_plot(sample_idx, sampleNames_ordered, 
                                            cov_of_region_matrix_ext, cov_of_region_matrix_norm_ext, 
                                            mean_std_covered_pos_samples[:,0], ## mean coverage of covered positions 
                                            me_element_dict['START'], me_element_dict['END'], 
                                            contig_edges, color_dict, cov_plot_title, scaling_factor = general_params_dict['group_nts_in_plot'], 
                                            outpath = f'{outpath_of_mge}/{mge_name}_{sampleNames_ordered[sample_idx]}_cov.png', dpi = 300)
            
            mge.plot_enrichment_of_region_heatmap(mean_cov_of_regio_mat, mean_std_covered_pos_samples[:,0], idx_to_insert_data, me_element_dict['TYPE'], 
                                                yticklabels=[tip_label.split('__')[-1] for tip_label in tip_labels], 
                                                xticklabels=[f"{str(me_element_dict['START'])}:{str(me_element_dict['END'])}"], 
                                                title = mge_name, 
                                                outpath = f'{outpath_of_mge}/{mge_name}_cov_enrichment_of_region.svg')
            
            mge.plot_cov_across_sorted_samples(cov_of_region_matrix,idx_to_insert_data,ymax=1,title=mge_name, 
                                                outpath = f'{outpath_of_mge}/{mge_name}_cov_across_phylo.png', dpi = 300)     
            ## save per mge the coverage which can be used for later plotting with ggtree 
            n = 25
            rows, cols = cov_of_region_matrix.shape
            remainder = cols % n
            if remainder != 0:
                pad_width = n - remainder
                padded = np.pad(cov_of_region_matrix.astype(float), ((0, 0), (0, pad_width)), mode='constant', constant_values=np.nan)
            else:
                padded = cov_of_region_matrix
            num_bins = padded.shape[1] // n
            binned = padded.reshape(rows, num_bins, n)
            binned_means = np.nanmean(binned, axis=2)
            cov_df = pd.DataFrame(binned_means / mean_std_covered_pos_samples[:,0].reshape(-1, 1), index=sampleNames_ordered)
            cov_df.to_csv(f'{outpath_of_mge}/{mge_name}_cov_norm_against_mean_genomewide_cov_{n}binwidth.csv.gz', compression='gzip')
            
            ## group samples per timepoint and plot enrichment of coverage as well as frequency of mge on timepoint
            samples_of_tp_df, samples_w_mge = mge.group_samples_by_timepoint(mean_cov_of_regio_mat, mean_std_covered_pos_samples[:,0], tip_labels[~tip_label_absent], unique_tp_dict, general_params_dict['avg_cov2be_present'])
            mge.plot_region_enrichment_per_timepoint(samples_of_tp_df, mge_name, outpath = f'{outpath_of_mge}/{mge_name}_cov_enrichment_of_region_per_tp.pdf')
            mge.plot_freq_of_mge_presence(samples_w_mge, mge_name, outpath = f'{outpath_of_mge}/{mge_name}_mge_presence_freq_tp.pdf')

    # save unique identified elements for annotation
    mobile_element = pd.DataFrame(me_element_dict_list)
    mobile_element.to_csv(f'{outpath}/mobile_genetic_elements_summary.csv',index=False)
    me_mean_cov = pd.DataFrame(me_mean_cov_dict_list)
    me_mean_cov.to_csv(f'{outpath}/mobile_genetic_elements_mean_cov.csv',index=False) ## normalized for sample specific coverage 
    me_breadth_cov = pd.DataFrame(me_breadth_cov_dict_list)
    me_breadth_cov.to_csv(f'{outpath}/mobile_genetic_elements_cov_breadth.csv',index=False) ## breadth per sample (cov >= 1)
    
    ## get a dict of isolates representing which isolate has a deletion: 
    isolate_w_deletion = {samplename: [] for samplename in np.sort(sampleNames_ordered)}
    for entry in me_element_dict_list:
        if entry['TYPE'] == 0:
            for sample in entry['ISOLATES'].split(';'):
                isolate_w_deletion[sample].append(entry['START'])
