#%% IMPORT PACKAGES

import numpy as np
from matplotlib import collections as mc
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import sys
import pandas as pd 
import os

SCRIPTS_DIRECTORY = os.getcwd() + '/'
sys.path.insert(0, SCRIPTS_DIRECTORY)

import analysispy_module as apy


# %% FUNCTIONS

############ 
## Identify candidate amplifications and deletions
############

def candidate_amplifications( array_cov_norm, threshold_loose, threshold_stringent, min_length, array_cov_raw, threshold_amp_multiplier, threshold_amp_multiplier_std, mean_std_cov_samples, threshold_cov_step_at_region_boundaries, fraction_of_region_for_cov_step, length_around_region_for_cov_step, min_length_in_region_for_cov_step):
    ''' get amp candidates loose threshold (sequence of min_length bp need to be above threshold) '''
    ''' output matrix header: sample index, start, end, amplification identifier (ie. 1) '''
    ''' candidate apmplifications in sample need to have at least x times (threshold_amp_multiplier) coverage in target region to be considered TP '''
    base_candidates = (array_cov_norm > threshold_loose).astype(int) # everything above threshold is 1 (True) and everything below threshold is 0 (False)
    genome_buffer = np.zeros((array_cov_norm.shape[0],1),dtype=np.int8)
    isone = np.concatenate((genome_buffer, base_candidates, genome_buffer),axis = 1) # Create an array that is 1 where array is 1 (ie. below threshold) and otherwise 0; add to each end an extra 0. 
    del base_candidates ## clean large matrix from memory 
    absdiff = np.abs(np.diff(isone,axis=1)) # mark border pos as 1 for all sequences of 1 (incl. sequence of single occurence) 
    del isone ## clean large matrix from memory 
    tpl_sample_ranges = np.where(absdiff == 1) # Runs start and end where absdiff is 1. reports tuple with sequence: [0] samples; [1] coord
    del absdiff ## clean large matrix from memory 
    samples = tpl_sample_ranges[0].reshape(-1, 2)[:,0] # extract sample indices
    ranges = tpl_sample_ranges[1].reshape(-1, 2) # get ranges or regions
    del tpl_sample_ranges ## clean large matrix from memory 
    boolMinLength = ranges[:,1]-ranges[:,0] > min_length # which ranges fullfill min_length
    samples = samples[boolMinLength] # filter samples
    ranges = ranges[boolMinLength,:] # filter ranges
    # filter stringent: avg norm2 has to exceed threshold_stringent AND minimum 
    boolStringentPass = np.zeros( len(samples) , dtype=bool) # all False
    min_coverage_for_amp = (threshold_amp_multiplier * mean_std_cov_samples[:, 0]) + (threshold_amp_multiplier_std * mean_std_cov_samples[:, 1])
    for i,s in enumerate(samples):
        if np.mean(array_cov_norm[s,ranges[i,0]:ranges[i,1]]) > threshold_stringent:
            region_length_10th = int(abs(ranges[i,0] - ranges[i,1]) * fraction_of_region_for_cov_step) ## calculate one tenth of the region length
            if region_length_10th < min_length_in_region_for_cov_step:
                region_length_10th = min_length_in_region_for_cov_step
            mean_cov_region = np.mean(array_cov_raw[s,ranges[i,0]:ranges[i,1]])
            mean_cov_before_region = np.mean( array_cov_raw[s, (ranges[i,0]-length_around_region_for_cov_step):ranges[i,0]] ) ## mean coverage of n bases before amp region
            mean_cov_first_nth_region = np.mean( array_cov_raw[s, ranges[i,0]:(ranges[i,0]+region_length_10th)] ) ## mean coverage of 1/nth of region length at first 1/nth of region
            mean_cov_last_nth_region = np.mean( array_cov_raw[s, (ranges[i,1]-region_length_10th):ranges[i,1]] ) ## mean coverage of 1/nth of region length at last 1/nth of region
            mean_cov_after_region = np.mean( array_cov_raw[s, ranges[i,1]:(ranges[i,1]+length_around_region_for_cov_step)] ) ## mean coverage of n bases after amp region
            
            mean_cov_above_exp_cov = (mean_cov_region > min_coverage_for_amp[s]) 
            cov_step_at_region_start = ((mean_cov_before_region * threshold_cov_step_at_region_boundaries) <= mean_cov_first_nth_region)
            cov_drop_at_region_end = ((mean_cov_after_region * threshold_cov_step_at_region_boundaries) <= mean_cov_last_nth_region)
            if (mean_cov_above_exp_cov & cov_step_at_region_start & cov_drop_at_region_end): # filter amps to have mean cov at least x times larger than mean cov(covered positions) + standard deviation (for high variable samples) [remove FP amps with marginal cov spikes]
                boolStringentPass[i] = True
    # return matrix for filtered amps
    res = np.array((samples[boolStringentPass],ranges[boolStringentPass,0],ranges[boolStringentPass,1])).T
    return np.append(res, np.ones([res.shape[0],1]), axis=1).astype(int)

def candidate_deletions( array_cov_norm, threshold_loose, threshold_stringent, min_length, array_cov_raw, threshold_raw_cov ):
    ''' get del candidates loose threshold (sequence of min_length bp need to be below threshold) '''
    ''' output matrix header: sample index, start, end, deletion identifier (ie. 0) '''
    ''' candidate deletion need to have at least one sample with a certain coverage (array_cov_raw , threshold_raw_cov) to filter false positives due to ''deletions'' in all samples '''
    ''' candidate deletion has to have median cov of 0 in target sample. Otherwise FPs occur due to overall low cov'''
    base_candidates = (array_cov_norm < threshold_loose).astype(int) # everything below threshold is 1 (True) and everything above threshold is 0 (False)
    genome_buffer = np.zeros((array_cov_norm.shape[0],1),dtype=np.int8)
    isone = np.concatenate((genome_buffer, base_candidates, genome_buffer),axis = 1) # Create an array that is 1 where array is 1 (ie. below threshold) and otherwise 0; add to each end an extra 0. 
    del base_candidates ## clean large matrix from memory 
    absdiff = np.abs(np.diff(isone,axis=1)) # mark border pos as 1 for all sequences of 1 (incl. sequence of single occurence) 
    del isone ## clean large matrix from memory 
    tpl_sample_ranges = np.where(absdiff == 1) # Runs start and end where absdiff is 1. reports tuple with sequence: [0] samples; [1] coord
    del absdiff ## clean large matrix from memory 
    samples = tpl_sample_ranges[0].reshape(-1, 2)[:,0] # extract sample indices
    ranges = tpl_sample_ranges[1].reshape(-1, 2) # get ranges or regions
    del tpl_sample_ranges ## clean large matrix from memory 
    boolMinLength = ranges[:,1]-ranges[:,0] > min_length # which ranges fullfill min_length
    samples = samples[boolMinLength] # filter samples
    ranges = ranges[boolMinLength,:] # filter ranges
    # filter stringent: avg norm2 has to exceed threshold_stringent 
    boolStringentPass = np.zeros( len(samples) , dtype=bool) # all False
    for i,s in enumerate(samples):
        if np.mean(array_cov_norm[s,ranges[i,0]:ranges[i,1]]) < threshold_stringent:
            if (np.max(np.median(array_cov_raw[:,ranges[i,0]:ranges[i,1]],axis=1)) > threshold_raw_cov): # remove cand if all samples no cov (median robust to artifactual cov spikes observed (thus no mean!))
                if (np.median(array_cov_raw[s,ranges[i,0]:ranges[i,1]]) == 0): # median in target region has to be 0 (remove FPs due to low cov in region, which can result otherwise in a hit)
                    boolStringentPass[i] = True
    # return matrix for filtered amps
    res = np.array((samples[boolStringentPass],ranges[boolStringentPass,0],ranges[boolStringentPass,1])).T
    return np.append(res, np.zeros([res.shape[0],1]), axis=1).astype(int)

def cand_mobile_element_region_merger(cand_amp_del, overlap_bp):
    ''' returns sorted amp and del candidates array with 3 additional columns: [merged_region_num,merged_region_start,merged_region_end] '''
    ''' merges amp /del candidates [format: sample, start, end, type] into continues regions if overlap by x bp '''
    ''' sorts by start pos to find overlapping region. requires start is continues numeric across genome (no separation by contigs)! '''
    ''' merges amps and dels separately '''
    merged_amp_dels = []
    for amp_del_type in range(2):
        amp_del = cand_amp_del[cand_amp_del[:,3] == amp_del_type] ## extract either del (0) or amps (1)
        cand_sort = amp_del[amp_del[:,1].argsort()] # sort by start position. requires continues pos...no separate contigs!
        mrg_ctr = 0
        mrg_roi_lol = [] # list of list with merged region, iteratively updated
        mrg_cand_row_lol = [] # list of list with rows of cand_amp_del part of merged_roi
        for i in range(cand_sort.shape[0]):
            start = cand_sort[i,1]
            end = cand_sort[i,2]
            if i == 0:
                mrg_roi_lol.append([start,end])
                mrg_cand_row_lol.append([i])
            else:
                if (mrg_roi_lol[mrg_ctr][1]-overlap_bp) > start: # overlap w/ last region (mrg_roi_lol[mrg_ctr][1] == end of last region); NOTE: Regions are sorted by start positions! 
                    mrg_cand_row_lol[mrg_ctr].append(i)
                    if mrg_roi_lol[mrg_ctr][1] < end: # iteration regions extends prev region
                        mrg_roi_lol[mrg_ctr][1] = end
                else: # no overlap. new region
                    mrg_ctr += 1 # new region identifier in lol's
                    mrg_roi_lol.append([start,end])
                    mrg_cand_row_lol.append([i])
        # add merged region data to cand_amp_del array
        ls_res = [] # build list with all data, turn to array, reshape
        for idx_mrg_roi,ls_rows in enumerate(mrg_cand_row_lol): # idx_mrg_roi == merged-region identifier
            for idx_row in ls_rows: # this double loop has length of num rows cand_amp_del and is ordered as such
                ls_res.extend([idx_mrg_roi, mrg_roi_lol[idx_mrg_roi][0],mrg_roi_lol[idx_mrg_roi][1] ])
        merged_amp_dels.append( np.concatenate( (cand_sort, np.array(ls_res).reshape( cand_sort.shape[0] , 3 ) ) ,axis = 1) )
    
    final_merged_amp_dels = np.concatenate(merged_amp_dels).astype(int)
    return final_merged_amp_dels

def split_regions_on_contig_edges(mge_regions_array, chrstarts, min_region_length):
    split_mge_regions_array = []
    
    for region in mge_regions_array:
        current_start_pos = region[2] 
        end_pos = region[3]

        for chrstart in chrstarts:
            if current_start_pos < chrstart <= end_pos:
                ## add entry only if it is longer than minimal region length to report!
                if ((chrstart-1) - current_start_pos) > min_region_length:
                    split_mge_regions_array.append([region[0], region[1], current_start_pos, chrstart - 1])
                current_start_pos = chrstart

        ## add entry only if it is longer than minimal region length to report!
        if (end_pos - current_start_pos) > min_region_length:
            split_mge_regions_array.append([region[0], region[1], current_start_pos, end_pos])

    return np.array(split_mge_regions_array)

############ 
## Genomic region parsing
############

def mean_std_cov_above_coverage(cov_matrix, cov_threshold = 0):
    covered_pos_bool = cov_matrix > cov_threshold
    filtered_matrix = cov_matrix * covered_pos_bool ## set values < thershold to 0
    sums = np.sum(filtered_matrix, axis=1)  # Sum values greater than threshold for each row
    counts = np.sum(covered_pos_bool, axis=1)  # Count of values greater than threshold for each row
    mean_of_mat = np.divide(sums, counts, where= counts!=0) ## calculate mean, if no position > threshold --> mean = 0
    
    # Calculate the standard deviation
    squared_diffs = (filtered_matrix - mean_of_mat[:, None])**2
    sum_squared_diffs = np.sum(squared_diffs * covered_pos_bool, axis = 1)
    std_of_mat = np.sqrt(np.divide(sum_squared_diffs, counts, where=counts != 0))
    return mean_of_mat, std_of_mat

def annotations_in_genomic_region(chr_pos_array,annotation_genes,print2screen=True,print2csv=None):
    ''' prints genes, distance to contig/chr start/end for genomic region (based on continues genome between first entry and last)'''
    ''' returns subset of annotation genes for genes in roi '''
    mychr = np.unique(chr_pos_array[:,0]) 
    pd.set_option('display.expand_frame_repr', False) # print all data to screen. prevents cut via "..."
    if len(mychr) > 1: # region spans multiple contigs
        gene_num_ls = []
        length_ls = []
        contig_pos_data_ls = []
        if print2screen:
            print('REGIONs GENE CONTENT')
        for i,chrom in enumerate(mychr):
            chr_data = chr_pos_array[ chr_pos_array[:,0]==chrom, : ]
            start_pos = chr_data[0,1] # first pos.
            end_pos = chr_data[-1,1] #last pos.
            contig_pos_data_ls.append([chrom, start_pos, end_pos])
            chr_anno = annotation_genes[ int(chrom-1) ] # 0-based
            if chr_anno.empty: # no genes on chr! skip!
                if i == 0: 
                    table_genes = pd.DataFrame() ## make empty dataframe to append if mutliple chromosomes
                continue
            #print(chr_data,start_pos,end_pos)
            chr_anno.rename({'indices': 'pos_start_end'}, axis = 1, inplace = True)
            chr_anno['contig'] = int(chrom) ## insert contig number in df
            chr_anno_roi = chr_anno.loc[(chr_anno['loc1'] >= start_pos) & (chr_anno['loc2'] <= end_pos)]
            if i == 0: # build continous dataframe with gene annotation across different contigs
                table_genes = chr_anno_roi
            else:
                table_genes = pd.concat([table_genes, chr_anno_roi], ignore_index=True, sort=True) 
        num_genes = table_genes.shape[0]
        total_length = end_pos - start_pos
        gene_num_ls.append(num_genes)
        length_ls.append(total_length)
        target_region_str = '; '.join([f'{str(chr)}:{str(start)}-{str(end)}' for chr,start,end in contig_pos_data_ls])
        if print2screen:
            print("Target region contig: "+ target_region_str)
            if not table_genes.empty:
                print(table_genes.loc[:,['gene','contig','pos_start_end','product','locustag']].to_string(index=False))
            else:
                print('No genes in that region on contig')
        if print2screen:
            print("\nNumber of contigs: "+str(len(mychr)))
            print('Number of genes: '+str(sum(gene_num_ls))+'. Length region (in bp): '+str(sum(length_ls)))
        if print2csv:
            print('NOTE: THE SEGMENT IS NOT SAVED TO CSV --> NO CURRENT IMPLEMENTATION FOR THAT!')
            with open(print2csv,'w') as file: # open and overwrite
                file.write("##Target region of interest: " +target_region_str+'\n')
                file.write('##Number of genes: '+str(num_genes)+'. Length region (in bp): '+str(total_length)+'\n')
                if not table_genes.empty:
                    file.write('#gene,contig,start,end,product,locustag\n')
                    np.savetxt(file, table_genes.loc[:,['gene','contig','loc1','loc2','product','locustag']],fmt='%s',delimiter=',') 
    else: # region spans single contigs
        start_pos = chr_pos_array[0,1] # first pos.
        end_pos = chr_pos_array[-1,1] #last pos.
        chr_anno = annotation_genes[ int(mychr-1) ] # 0-based
        if chr_anno.empty:
            print('No Genes in Region!') # no genes on chr! skip!
            return pd.DataFrame()
        chr_anno.rename({'indices': 'pos_start_end'}, axis = 1, inplace = True)
        chr_anno['contig'] = int(mychr[0]) ## insert contig number in df
        chr_anno_roi = chr_anno.loc[(chr_anno['loc1'] >= start_pos) & (chr_anno['loc2'] <= end_pos)]
        num_genes = chr_anno_roi.shape[0]
        total_length = end_pos - start_pos
        table_genes = chr_anno_roi
        if print2screen:
            print('REGIONs GENE CONTENT')
            print("Target region of interest: "+str(mychr[0])+":"+str(start_pos)+"-"+str(end_pos))
            if not chr_anno_roi.empty:
                print(chr_anno_roi.loc[:,['gene','contig','pos_start_end','product','locustag']].to_string(index=False))
            else:
                print('No genes in that region on contig')
            print('\nNumber of genes: '+str(num_genes)+'.\nLength region (in bp): '+str(total_length))
        if print2csv:
            with open(print2csv,'w') as file: # open and overwrite
                file.write("##Target region of interest:"+str(mychr[0])+":"+str(start_pos)+"-"+str(end_pos)+'\n')
                file.write('##Number of genes: '+str(num_genes)+'. Length region (in bp): '+str(total_length)+'\n')
                if not chr_anno_roi.empty:
                    file.write('#gene,contig,start,end,product,locustag\n')
                    np.savetxt(file, chr_anno_roi.loc[:,['gene','contig','loc1','loc2','product','locustag']],fmt='%s',delimiter=',')    
    return table_genes

def get_relative_coverage_across_mobile_elements(mobile_elements,cov_matrix,mean_std_samples,tip_label_indexCovMat,tip_label_absent,amp_cap=4):
    ''' returns matrix of relative coverge observed in ME's ordered as in phylogeny; returns tip labels absent in cov as 0, amplifications capped at amp-cap'''
    ''' cap of amplification necessary to avoid extreme amplifications bias too much color scale (high values only extreme cases)'''
    me_cov = np.array([])
    for me_id in np.unique(mobile_elements[:,4]): # loop over mobile element identifies
        start=mobile_elements[ mobile_elements[:,4]==me_id , :][0,5] # get from first row merged_start and merged_stop
        end=mobile_elements[ mobile_elements[:,4]==me_id , :][0,6]
        me_cov = np.append(me_cov, np.median(cov_matrix[:, start:end ],axis=1)/mean_std_samples[:,0])
    me_cov = me_cov.reshape( cov_matrix.shape[0] , len(np.unique(mobile_elements[:,4])) , order='F' ) # order='F' reshape column after column
    me_cov = np.round(me_cov,0) # turn to integer for heatmap.     
    # cap at amp_cap. 
    me_cov[ me_cov > amp_cap ] = amp_cap
    # me_cov = me_cov/10
    me_cov[ me_cov > amp_cap ] = amp_cap    
    # order by tree and put samples not in cov to 0
    me_cov_tree = me_cov[tip_label_indexCovMat,:]
    me_cov_tree[tip_label_absent,:] = 0
    return me_cov_tree

def extract_metadata_of_mge(me_coord, cov_raw, cand_amp_del, mean_coverage_isolate, annotation_genes, sampleNames_ordered, chrstarts, snp_table = pd.DataFrame(), indel_table = pd.DataFrame()):
    mobile_element_dict = {'TYPE':me_coord[0],'NUM':me_coord[1],'START':int(me_coord[2]),'END':int(me_coord[3]),'LENGTH':(me_coord[3]-me_coord[2]),
                            'COORDINATES':'','GENES':'','PRODUCTS':'',
                            'SNP_POS':'','SNP_TYPES':'','SNP_LOCUSTAGS':'','SNPs':'',
                            'INDEL_START':'','INDEL_END':'','INDEL_TYPES':'','INDEL_LOCUSTAGS':'','INDELs':'',
                            'NUM_ISOLATES':0,'ISOLATES':''}
    posonchr = apy.p2chrpos(np.arange(mobile_element_dict['START'],mobile_element_dict['END'],1),chrstarts)
    genes_in_me = annotations_in_genomic_region(posonchr,annotation_genes,print2screen=False)
    if not genes_in_me.empty:
        # unique_genes, index_genes, gene_counts = np.unique(genes_in_me['product'], return_index=True, return_counts=True) # get unique genes and count
        # sorted_index = np.argsort(index_genes)
        # unique_genes_sorted = unique_genes[sorted_index] ## sort genes and counts based on how they appear on MGE
        # gene_counts_sorted = gene_counts[sorted_index]
        # genes_uq_ct = np.unique(genes_in_me['product'],return_counts=True) 
        # mobile_element_dict['GENES'] = '; '.join([f'{str(cnt)}_{gene}' for cnt,gene in zip(gene_counts_sorted,unique_genes_sorted)])
        mobile_element_dict['GENES'] = '; '.join(genes_in_me['gene'].fillna(''))
        mobile_element_dict['PRODUCTS'] = '; '.join(genes_in_me['product'].fillna(''))
    contigs_target = np.unique(posonchr[:,0])
    string_w_target_roi = ''
    
    for c in contigs_target:
        bool_contig = (posonchr[:,0]==c)
        contig_id = int(c)
        posonchr_of_interest = posonchr[bool_contig,:]
        posonchr_of_interest_chr = posonchr_of_interest[0, 0]
        posonchr_of_interest_start = int(posonchr_of_interest[0,1])
        posonchr_of_interest_end = int(posonchr_of_interest[-1,1])
        string_w_target_roi = f'{string_w_target_roi}{str(contig_id)}:{str(posonchr_of_interest_start)}-{str(posonchr_of_interest_end)};'
        ## identify snps falling into mge region
        if not snp_table.empty:
            snp_covered_by_mge_bool = (posonchr_of_interest_chr == snp_table['chr'].astype(int)) & (snp_table['pos'].astype(int).between(posonchr_of_interest_start+1, posonchr_of_interest_end+1)) ## between: (start <= pos <= end) // +1 as in snp_table the pos are 1-based, on mge are 0-based!
            if sum(snp_covered_by_mge_bool) > 0:
                mobile_element_dict['SNP_POS'] += (';'.join(snp_table.loc[snp_covered_by_mge_bool, 'snp_pos'].values) + ';')
                mobile_element_dict['SNP_TYPES'] += (';'.join(snp_table.loc[snp_covered_by_mge_bool, 'type'].values) + ';')
                mobile_element_dict['SNP_LOCUSTAGS'] += (';'.join(snp_table.loc[snp_covered_by_mge_bool, 'locustag'].values) + ';')
                mobile_element_dict['SNPs'] += (';'.join(snp_table.loc[snp_covered_by_mge_bool, 'mutation'].values) + ';')
        if not indel_table.empty:
            chr_is_indel_chr = (posonchr_of_interest_chr == indel_table['chr'].astype(int))
            indel_within_mge = (indel_table['loc1'].astype(int).between(posonchr_of_interest_start+1, posonchr_of_interest_end+1)) | (indel_table['loc2'].astype(int).between(posonchr_of_interest_start+1, posonchr_of_interest_end+1))
            indel_covered_by_mge_bool = (chr_is_indel_chr & indel_within_mge) ## between: (start <= pos <= end) // +1 as in snp_table the pos are 1-based, on mge are 0-based!
            if sum(snp_covered_by_mge_bool) > 0:
                mobile_element_dict['INDEL_START'] += (';'.join(indel_table.loc[indel_covered_by_mge_bool, 'indel_start'].values) + ';')
                mobile_element_dict['INDEL_END'] += (';'.join(indel_table.loc[indel_covered_by_mge_bool, 'indel_end'].values) + ';')
                mobile_element_dict['INDEL_TYPES'] += (';'.join(indel_table.loc[indel_covered_by_mge_bool, 'type'].values) + ';')
                mobile_element_dict['INDEL_LOCUSTAGS'] += (';'.join(indel_table.loc[indel_covered_by_mge_bool, 'locustag'].values) + ';')
                mobile_element_dict['INDELs'] += (';'.join(indel_table.loc[indel_covered_by_mge_bool, 'mutation'].values) + ';')
    for key in ['SNP_POS', 'SNP_TYPES', 'SNP_LOCUSTAGS', 'SNPs', 'INDEL_START', 'INDEL_END', 'INDEL_LOCUSTAGS', 'INDELs']:
        mobile_element_dict[key] = mobile_element_dict[key].strip(';')
    mobile_element_dict['COORDINATES'] = string_w_target_roi.strip(';')
    isolates_w_ampdel_region = (mobile_element_dict['START'] <= cand_amp_del[:, 1]) & (mobile_element_dict['END'] >= cand_amp_del[:, 2]) & (mobile_element_dict['TYPE'] == cand_amp_del[:, 3])## identify isolates with amplification/deletion region
    cand_amp_del_in_region = cand_amp_del[isolates_w_ampdel_region, :]
    samples_with_region = np.unique(cand_amp_del_in_region[:, 0])
    mobile_element_dict['NUM_ISOLATES'] = int(len(samples_with_region))
    mobile_element_dict['ISOLATES'] = ';'.join(np.sort(sampleNames_ordered[samples_with_region]))
    if np.std(cand_amp_del[isolates_w_ampdel_region, 3]) != 0:
        # print(f'Amplification and deletion was found on genomic position {int(mobile_element_dict["START"])}:{int(mobile_element_dict["END"])}.', end = '\t')
        # print('Overwriting type to deletion')
        mobile_element_dict['TYPE'] = 0 ## type for deletions = 0, amplifications = 1
    ## extract the mean coverage within region per sample
    mean_cov_per_isolate = np.mean(cov_raw[:,int(me_coord[2]):int(me_coord[3])], axis = 1) / mean_coverage_isolate
    mobile_element_mean_cov_dict = {isolateID: mean_cov for isolateID, mean_cov in zip(sampleNames_ordered, mean_cov_per_isolate)}
    mobile_element_mean_cov_dict = {**{'TYPE':me_coord[0],'START':me_coord[2],'END':me_coord[3],'LENGTH':(me_coord[3]-me_coord[2])}, **mobile_element_mean_cov_dict}
    ## extract breadth of reagion covered
    breadth_cov_per_isolate = np.mean(cov_raw[:,int(me_coord[2]):int(me_coord[3])] >= 1, axis = 1)
    mobile_element_breadth_cov_dict = {isolateID: breadth_cov for isolateID, breadth_cov in zip(sampleNames_ordered, breadth_cov_per_isolate)}
    mobile_element_breadth_cov_dict = {**{'TYPE':me_coord[0],'START':me_coord[2],'END':me_coord[3],'LENGTH':(me_coord[3]-me_coord[2])}, **mobile_element_breadth_cov_dict}
    return mobile_element_dict, isolates_w_ampdel_region, mobile_element_mean_cov_dict, mobile_element_breadth_cov_dict

def read_snp_table(snp_table_path):
    if not os.path.exists(snp_table_path):
        return pd.DataFrame()
    snp_table = pd.read_csv(snp_table_path)    
    snp_table['snp_pos'] = snp_table[['chr', 'pos']].astype(str).apply(lambda x:'__'.join(x), axis = 1)
    snp_table['locustag'] = snp_table['locustag'].str.replace(';', ',')
    snp_table['mutation'] = snp_table['gene'].str.replace(';', ',')
    snp_table.loc[snp_table['mutation'] == '.', 'mutation'] = snp_table.loc[snp_table['mutation'] == '.', 'locustag']
    snp_table.loc[snp_table['type'].str.contains('N'), 'mutation'] += snp_table.loc[snp_table['type'].str.contains('N'), 'muts']
    snp_table = snp_table.fillna('')
    return snp_table[['chr', 'pos', 'snp_pos', 'type', 'locustag', 'mutation']]

def read_indel_table(indel_table_path):
    if not os.path.exists(indel_table_path):
        return pd.DataFrame()
    indel_table = pd.read_csv(indel_table_path)    
    indel_table['indel_start'] = indel_table[['chr', 'loc1']].astype(str).apply(lambda x:'__'.join(x), axis = 1)
    indel_table['indel_end'] = indel_table[['chr', 'loc2']].astype(str).apply(lambda x:'__'.join(x), axis = 1)
    indel_table['locustag'] = indel_table['locustag'].str.replace(';', ',')
    indel_table['mutation'] = indel_table['gene'].str.replace(';', ',')
    indel_table['alt_stop'] = indel_table['alt_stop'].str.replace(r"\[|\]| |'", '', regex = True)
    is_intergenic_indel = indel_table['loc1'].isna()
    indel_table.loc[is_intergenic_indel, 'loc2'] = indel_table.loc[is_intergenic_indel, 'indel_size_gp'].str.replace(r"\[|\]| |'", '', regex = True).astype(float).max() 
    indel_table.loc[is_intergenic_indel, 'loc1'] = indel_table.loc[is_intergenic_indel, 'pos']
    indel_table.loc[is_intergenic_indel, 'loc2'] = indel_table.loc[is_intergenic_indel, 'loc2'] + indel_table.loc[is_intergenic_indel, 'loc1']
    indel_table.loc[indel_table['mutation'] == '.', 'mutation'] = indel_table.loc[indel_table['mutation'] == '.', 'locustag']
    indel_table.loc[indel_table['type'].str.contains('G'), 'mutation'] += '[' + indel_table.loc[indel_table['type'].str.contains('G'), 'alt_stop'] + ']'
    indel_table = indel_table.fillna('')
    return indel_table[['chr', 'pos', 'indel_start', 'indel_end', 'type', 'locustag', 'mutation', 'loc1', 'loc2']]

def downsample_mean_axis_1(arr, factor=10):
    ## see of array size on axis 1 is a multiple of factor 
    remainder = arr.shape[1] % factor
    # if no remained, directly calculate means of size factor and return
    if remainder == 0:
        return arr.reshape(arr.shape[0], -1, factor).mean(axis=2)
    else:
        # Compute mean for the complete chunks and remainder separately
        complete_mean = arr[:, :-remainder].reshape(arr.shape[0], -1, factor).mean(axis=2)
        remainder_mean = arr[:, -remainder:].mean(axis=1, keepdims=True)
        return np.hstack((complete_mean, remainder_mean))


def group_samples_by_timepoint(mean_cov_of_regio_mat, mean_cov_samples, tip_labels_present, unique_tp_dict, avg_cov2be_present):
    ## group samples per timepoint together
    samples_of_tp_dict = {}
    samples_w_mge = {}
    
    if unique_tp_dict != {}:    
        for tp, date in unique_tp_dict.items():
            samples_of_tp_bool = [(f'_{tp}_' in tip_label) for tip_label in tip_labels_present]
            samples_of_tp_dict[f'{tp}\n{date}'] = mean_cov_of_regio_mat[samples_of_tp_bool] / mean_cov_samples[samples_of_tp_bool]
            num_samples_w_mge = np.sum(mean_cov_of_regio_mat[samples_of_tp_bool] > avg_cov2be_present)
            num_samples_at_tp = np.sum(samples_of_tp_bool)
            samples_w_mge[f'{tp}\n{date}'] = num_samples_w_mge / num_samples_at_tp
    else:
        samples_of_tp_dict['all'] = mean_cov_of_regio_mat
        samples_w_mge['all'] = np.sum(np.sum(mean_cov_of_regio_mat > avg_cov2be_present)) / mean_cov_of_regio_mat.shape[0]
    samples_of_tp_df = pd.DataFrame([(key, value) for key, values in samples_of_tp_dict.items() for value in values], columns = ['Timepoint', 'values'])
    return samples_of_tp_df, samples_w_mge



############ 
## Tweak sample names
############

def get_tip_label_order_for_cov_matrix(tip_labels,sampleNames):
    ''' return index of cov matrix rows in order of tree tips, return cov matrix row indices that are NOT in tree'''
    # get amp/del matrices index present in tree and in leaf order of tree
    tip_label_indexCovMat = []
    tip_label_absent = np.zeros(( tip_labels.shape[0] ), dtype=bool)
    for idx,label in enumerate(tip_labels):
        if label in sampleNames:
            index = np.where(sampleNames == label)[0][0]
            tip_label_indexCovMat.append(index)
        else: # tip_label in tree but not present in sampleNames (ie. amp-del inference)
            tip_label_indexCovMat.append(-1) # put index -1 as placeholder for thse cases
            tip_label_absent[idx] = True # bool of tip_label index that are not in cov matrix. needed to plot grey.
    tip_label_indexCovMat = np.array(tip_label_indexCovMat) # array of cov indices that are present in tree, and in order of tree leafes
    return [tip_label_indexCovMat,tip_label_absent]


def names_of_a_found_in_b(a,b):
    ''' Returns bool of length(a) with True when string found in b'''
    bool_a = np.zeros(( a.shape[0] ), dtype=bool)
    for i,spl in enumerate(a):
        if spl in b:
            bool_a[i]=True
    return bool_a


def extract_and_rename_sampleIDs_from_tree(sample):
    '''
    Specific function to extract sample names (separated by '[some_label]__[SampleID]') and rename outgroup names (e.g. OG_GCF_003408555_P07-c1_TP0 --> GCF_003408555_P07-c1)
    '''
    sample_name = sample.split('__')[-1]
    if sample_name.startswith('OG_'):
        return # sample_name = sample_name.replace('OG_', '').replace('_TP0', '')  # Extract from outgroups only sample ID
    elif sample_name.endswith('_ref'): ## remove ref from tips
        return 
    return sample_name

def rename_outgroups_sampleNames(sample):
    '''
    Specific function to extract sample names from outgroup names (e.g. Ehormaechei_GCF_010319625_P07-c1 --> GCF_010319625_P07-c1)
    '''
    if 'GCF_' in sample:
        return 'GCF_' + sample.split('GCF_')[-1]
    elif 'GCA_' in sample:
        return 'GCA_' + sample.split('GCA_')[-1]
    elif 'ASM_' in sample:
        return 'ASM_' + sample.split('ASM_')[-1]
    return sample

############ 
## PLOTTING
############

def generate_plot_region(mge_start, mge_end, bumper_region, genome_length):
    # Calculate bumper zone of cov plot around del/amp region
    bumper_region = int(bumper_region)
    # Start minus bumper
    mge_start_plot = mge_start-bumper_region
    if mge_start_plot < 0:
        mge_start_plot = 0
    # End plus bumper
    mge_end_plot = mge_end+bumper_region
    if mge_end_plot > genome_length:
        mge_end_plot = genome_length
    return mge_start_plot, mge_end_plot

def plot_coverage( ax, segment_sample_idx, sampleNames, region_of_interest_mat, mean_coverage,
                  segment_pos_start, segment_pos_end, 
                  segment_pos_start_bumper, segment_pos_end_bumper,
                  str_cov_type, colordict , contig_edges, scaling_factor):
    ''' makes a plot showing coverage of all samples over the region of interest, with the sample of interest highlighted (red for candidate deletion; blue for candidate amplification) '''
    # Basic setup
    ax.margins(.01)
    
    xcoord_start_of_segment = (segment_pos_start-segment_pos_start_bumper) / scaling_factor
    xcoord_end_of_segment = (segment_pos_end-segment_pos_start_bumper) / scaling_factor
    xcoord_max_x = (segment_pos_end_bumper-segment_pos_start_bumper) / scaling_factor

    # Plot all samples in semi-transparent black 
    for sid, samplename in enumerate(sampleNames):
        if samplename != sampleNames[segment_sample_idx]: ## do not plot sample which have been clicked
            ax.plot(region_of_interest_mat[sid,:],color=colordict['samples_color'],linewidth=.25)
    # Plot this sample in color
    ax.plot(region_of_interest_mat[segment_sample_idx,:],color=colordict['highlight_sample_color'],linewidth=1)
    # Plot vertical lines at edge of candidate region
    ax.axvline(xcoord_start_of_segment, color=colordict['cand_mge_edges_color'],linewidth=1,linestyle='--')
    ax.axvline(xcoord_end_of_segment, color=colordict['cand_mge_edges_color'],linewidth=1,linestyle='--')
    # Put lines in background indicating contig edges
    for pos in contig_edges:
        if (pos > segment_pos_start_bumper) & (pos < segment_pos_end_bumper):
            xcoord_contig_edge = (pos-segment_pos_start_bumper) / scaling_factor
            ax.axvline(xcoord_contig_edge, color=colordict['contig_edges_color'],linewidth=1)            
    # abs cov plots: add vertical for mean cov and 2x mean cov
    if str_cov_type == 'abs':
        y_mean_val = mean_coverage[segment_sample_idx]
        ax.axhline(y_mean_val, color='slategray', linestyle='-.')
        ax.axhline(2*y_mean_val, color='slategray', linestyle='dotted')

    # Label plot
    ax.set_xticks([0, xcoord_start_of_segment, xcoord_end_of_segment, xcoord_max_x])
    ax.set_xticklabels([segment_pos_start_bumper, segment_pos_start, segment_pos_end, segment_pos_end_bumper])
    ax.set_xlabel('Genome position (bp)')
    ax.set_ylabel('Coverage ('+str_cov_type+')')


def create_norm_abs_cov_plot(sample_idx, sampleNames, region_of_interest_matrix, region_of_interest_norm, mean_coverage, mge_start, mge_end, contig_edges, color_dict, plot_title = '', scaling_factor = 1, outpath = '', dpi = 300):
    plot_region_start, plot_region_end = generate_plot_region(mge_start, mge_end, mge_end-mge_start, contig_edges[-1])

    fig, (ax1, ax2) = plt.subplots(figsize = (6, 5), nrows = 2)
    # Create subplot for normalized coverage
    plot_coverage( ax1, sample_idx, sampleNames, region_of_interest_norm, mean_coverage, 
                        mge_start, mge_end, plot_region_start, plot_region_end,
                        'norm', color_dict, contig_edges, scaling_factor )

    # Create subplot for actual coverage
    plot_coverage( ax2, sample_idx, sampleNames, region_of_interest_matrix, mean_coverage,
                        mge_start, mge_end, plot_region_start, plot_region_end,
                        'abs', color_dict, contig_edges, scaling_factor)
    ymax_cov_plot = 2*np.max(region_of_interest_matrix[sample_idx,:])
    if ymax_cov_plot == 0:
        ymax_cov_plot = 20 
    ax2.set_ylim([0,ymax_cov_plot])
    fig.suptitle(plot_title, fontsize=12)
    fig.subplots_adjust( hspace=.4 )
    plt.tight_layout()
    if outpath: 
        fig.savefig(outpath, dpi = dpi)
        plt.close()
    else:
        fig.show()


def annotate_heatmap(heatmap, data=None, valfmt="{x:.2f}", textcolors=["black", "white"], threshold=None, **textkw):

    if not isinstance(data, (list, np.ndarray)):
        data = heatmap.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = heatmap.norm(threshold)
    else:
        threshold = heatmap.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(heatmap.norm(data[i, j]) > threshold)])
            text = heatmap.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

def plot_enrichment_of_region_heatmap(mean_cov_region_mat, mean_cov_samples, idx_to_insert_data, mge_type, 
                                      yticklabels='', xticklabels='', xlab=None, ylab=None, title = '', 
                                      put_vals_into_heatmap = False, outpath = '', dpi = 300, **kwargs):
    
    rel_cov_roi = (mean_cov_region_mat / mean_cov_samples).reshape(mean_cov_samples.shape[0], 1)
    # np.mean(cov_matrix_ordered[:,me_element_dict['START']:me_element_dict['END']],axis=1)/mean_std_cov_samples[:,0]).reshape(n_samples,1)
    if int(mge_type) == 1: ## amplification
        my_colorbar = 'YlGnBu'
        rel_cov_roi = np.insert(rel_cov_roi, idx_to_insert_data, [0], axis = 0) ## insert missing data
    else: ## deletion
        my_colorbar = 'YlOrRd'
        rel_cov_roi = np.insert(rel_cov_roi, idx_to_insert_data, [1], axis = 0) ## insert missing data
    cbarlabel = "Enrichment of region above mean sample coverage"
    plt_height = rel_cov_roi.shape[0]/6

    fig, ax = plt.subplots(figsize=(6,plt_height))
    # Plot the heatmap
    heatmap = ax.imshow(rel_cov_roi, cmap = my_colorbar, interpolation='nearest', aspect='auto', **kwargs)
    # Create colorbar
    cbar = fig.colorbar(heatmap) 
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    
    if put_vals_into_heatmap: 
        annotate_heatmap(heatmap)

    ax.set_xticks(np.arange(rel_cov_roi.shape[1]))
    ax.set_yticks(np.arange(rel_cov_roi.shape[0]))
    ax.set_xticklabels(xticklabels)
    ax.set_yticklabels(yticklabels)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if title:
        ax.set_title(title)
    
    # # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    plt.tight_layout()
    if outpath:
        fig.savefig(outpath, bbox_inches = 'tight', dpi = dpi)
        plt.close()
    else:
        fig.show()


def plot_cov_across_sorted_samples(cov_of_region_matrix,idx_to_insert_data,ymax,title='',outpath='', dpi = 300):
    ''' return plot coverage distribution across all samples (ordered as tip labels in tree), tips not in matrix cov == ymax (uniform,grey)'''
    ''' all subplots y-axis from 0 to ymax'''
    ''' write plot.png if desired '''    
    me_cov_tree_cand = np.insert(cov_of_region_matrix, idx_to_insert_data, [ymax], axis = 0)
    (num_samples, length_of_region) = np.shape(me_cov_tree_cand)

    fig = plt.figure(figsize=(2, num_samples/6))
    fig.suptitle(title)
    for i in range(num_samples):
        ax = plt.subplot(num_samples,1,1+i,ylim=(0,ymax))
        if i in idx_to_insert_data:
            ax.fill_between(np.arange(0,length_of_region, 1), 0, me_cov_tree_cand[i,],color='grey') # plot NAs grey (samples not in cov matrix (at least no name match))
        else:
            ax.fill_between(np.arange(0,length_of_region, 1), 0, me_cov_tree_cand[i,])
        plt.subplots_adjust(hspace = .001)
        plt.tick_params(
            axis='both',       # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False, # labels along the bottom edge are off   
            left=False,        # y ticks off
            labelleft=False)   # y labels off
    if outpath:
        fig.savefig(outpath, dpi = dpi)
        plt.close()
    else:
        fig.show()

# How to make a plot of candidate amp/del segments
def interactive_plotting_prerequisite(cov_matrix,cand_amp_del,line_bumper):
    ''' creates line collection consisting of candidate segments; prerequisite calc in function for candidate segment figure'''
    
    n_samples = cov_matrix.shape[0]
    print('Number of samples: '+str(n_samples))
    n_genome = cov_matrix.shape[1]
    print('Size of genome: '+str(n_genome))
    
    # Make lists used below
    list_samples = cand_amp_del[:,0]
    list_starts = cand_amp_del[:,1]
    list_ends = cand_amp_del[:,2]
    list_is_amp = cand_amp_del[:,3]
    
    # Generate lines for line collection
    # Note: this could be over/under the genome length but it doesn't really matter
    lines = []
    for i,sample in enumerate(list_samples):
        lines.append([(list_starts[i]-line_bumper,sample),(list_ends[i]+line_bumper,sample)])
        # TODO: figure out how to avoid appending things

    # Generate colors depending on deletion or amplification
    red = (1, 0, 0, 1)
    green = (0, 1, 0, 1)
    blue = (0, 0, 1, 1)
    black = (0, 0, 0, .25)
    lines_colors = []
    for i,sample in enumerate(list_samples):
        if list_is_amp[i]: # amplification
            lines_colors.append( blue )
        else: # deletion
            lines_colors.append( red )
        # TODO: figure out how to avoid appending things
    return [lines, lines_colors, n_genome, n_samples,list_samples,list_starts,list_ends,list_is_amp]


# Set up the click callback for line segments
def onclick_plot_report(event,args_dict):
    ''' used to make the lines on the candidate amplification and deletion plot clickable 
        added optional report: 
        - report annotated genes in region
        - plot tree-sorted heatmap of coverage of candidate region for each sample
    '''
    
    # Note: Add more checks here if you are putting more keys in args_dict

    # Unwrap variables
    # Info on candidate amplifications and deletions
    list_samples = args_dict['list_samples']
    list_starts = args_dict['list_starts']
    list_ends = args_dict['list_ends']
    lines_colors = args_dict['lines_colors']
    list_type = args_dict['list_type']
    # Coverage matrices
    array_cov_norm = args_dict['cov_norm']
    array_cov = args_dict['cov_matrix']
    mean_std_cov_samples=args_dict['mean_std_cov_samples']
    # Info on samples and genome
    n_samples = args_dict['n_samples']
    n_genome = args_dict['n_genome']
    sampleNames = args_dict['SampleNames']
    reffolder = args_dict['reffolder']
    # Info on samples in tree vs samples in coverage matrix
    tip_label_indexCovMat = args_dict['tip_label_indexCovMat'] # sample order as in tree 
    tip_label_absent = args_dict['tip_label_absent'] # samples in tree but not in data
    idx_to_insert_data = args_dict['idx_to_insert_data'] # indices of samples in tree which are not in cov matrix (or have been excluded) to insert data based on pyhlogeny 
    # Clickable plot options
    bool_report_annotation = args_dict['report_genes']
    bool_tree_aligned_cov_plot = args_dict['tree_aligned_cov_plot'] # Plot coverage for each sample ordered as *latest.nwk.tree
    bool_plot_relative_cov_tree_aligned = args_dict['bool_plot_relative_cov_tree_aligned']
    plot2file_flag = args_dict['plot2file'] # flag to print cov to plot.png 
    
    ## get genome stats 
    chrstarts,_,scaf_names = apy.genomestats(reffolder)
    contigedges = np.append(chrstarts, np.array([n_genome]))

    # Get index (or indices) of the line(s) clicked
    print('click event(s):', event.ind)

    for dataind in event.ind:
        # Get candidate segment and sample indices
        clicked_segment_id = int(dataind)
        clicked_segment_sample_idx = int(list_samples[clicked_segment_id])
        clicked_segment_sample_name = sampleNames[clicked_segment_sample_idx] 
        clicked_segment_type = int(list_type[clicked_segment_id]) # 0: del , 1: amp
        
        # Get relevant positions
        clicked_segment_pos_start = int(list_starts[clicked_segment_id])
        clicked_segment_pos_end = int(list_ends[clicked_segment_id])
        clicked_segment_length = clicked_segment_pos_end - clicked_segment_pos_start
        
        # print some info about ROI
        print('\n'+'Current segment index: '+str(dataind)+',  current sample: '+str(clicked_segment_sample_idx))
        print('\n'+'Current sample name: '+clicked_segment_sample_name)
        print('Genomic region: ', str(clicked_segment_pos_start)," - ", str(clicked_segment_pos_end))
        # print genes annotated in region 
        if bool_report_annotation:
            annotation_genes = apy.parse_gff(reffolder,scaf_names)
            posonchr = apy.p2chrpos(np.arange(clicked_segment_pos_start,clicked_segment_pos_end,1),chrstarts)
            annotations_in_genomic_region(posonchr,annotation_genes)

        cov_of_region_matrix = array_cov[:,clicked_segment_pos_start:clicked_segment_pos_end]
        plot_region_start, plot_region_end = generate_plot_region(clicked_segment_pos_start, clicked_segment_pos_end, clicked_segment_pos_end-clicked_segment_pos_start, contigedges[-1])
        ## Plot coverage of selected sample at selected region
        plot_title = f'Sample index: {str(clicked_segment_sample_idx)}, Sample ID: {str(clicked_segment_sample_name)},\nSegment index: {str(clicked_segment_id)}, Segment length: {str(clicked_segment_length)}'
        color_dict = {'cand_mge_edges_color': (0, 1, 0, 1), 
                          'highlight_sample_color': lines_colors[clicked_segment_id],
                          'samples_color': (0, 0, 0, .25), 
                          'contig_edges_color': (0.6, 0.2, 0.8, .5)}
        create_norm_abs_cov_plot(clicked_segment_sample_idx, sampleNames, 
                                array_cov[:,plot_region_start:plot_region_end], 
                                array_cov_norm[:,plot_region_start:plot_region_end], 
                                mean_std_cov_samples[:, 0], 
                                clicked_segment_pos_start, clicked_segment_pos_end, 
                                contigedges, color_dict, plot_title)
        
        # plot heatmap cov for each sample (ordered as tree -- precalculated)
        # plots relative coverage target region / mean genome-wide
        if bool_plot_relative_cov_tree_aligned:
            mean_cov_region_mat = np.mean(cov_of_region_matrix,axis=1)
            plot_enrichment_of_region_heatmap(mean_cov_region_mat, mean_std_cov_samples[:,0], idx_to_insert_data, clicked_segment_type, 
                                      xticklabels=[f'{str(clicked_segment_pos_start)}:{str(clicked_segment_pos_end)}'], 
                                      put_vals_into_heatmap = True)

        ## Plot coverage distribution for target region aligned to tree
        # ymax set to 1; > provides info on presence/absence only; ymax used to be ymax_cov_plot/2 but that is not usable for amp's
        if bool_tree_aligned_cov_plot:
            ymax = 1
            plot_cov_across_sorted_samples(cov_of_region_matrix,idx_to_insert_data,ymax=ymax,
                                           title=f'{str(clicked_segment_pos_start)}:{str(clicked_segment_pos_end)}\nymax: {str(np.round(ymax,1))}')   
    return True


def fig_candidate_segments( lines, lines_colors, n_genome, n_samples,my_data_dict, sampleNames=None, chrstarts=None):
    ''' makes a plot showing candidate amplification and deletion segments over the genome for each sample 
        NOTE: sampleNames need to be in reverted order to be in the correct shape for plotting! (e.g. list(reversed([1, 2, 3]))) '''
    plt.close(1) # Necessary so that you don't get multiple instances of the clicker running at the same time
    plot_width = (n_genome/5e5) + 1
    if plot_width < 4:
        plot_width = 4
    plot_height = n_samples / 8
    if plot_height < 4:
        plot_height = 4
    fig = plt.figure(1,figsize=(plot_width,plot_height)) 
    ax = fig.add_subplot(111)
    ax.set_xlim([0, n_genome])
    ax.set_ylim([-1, n_samples])
    ax.set_ylim([n_samples, -1]) ## reverse order of plot 
    # Put lines in background indicating contig edges
    if chrstarts is None:
        print('Warning! chrstarts not provided.')
    else:
        lines_contigs = []
        lines_contigs_colors = []
        gray = (0, 0, 0, .1)
        for i,pos in enumerate(chrstarts):
            lines_contigs.append([(pos,-1),(pos,n_samples)])
            lines_contigs_colors.append( gray )
        lc_contigs = mc.LineCollection( lines_contigs, colors=lines_contigs_colors, linewidths=1 )
        ax.add_collection(lc_contigs)
    # Put line segments for each candidate region on plot
    lc = mc.LineCollection(lines, colors=lines_colors, linewidths=5, picker=True)
    ax.add_collection(lc)
    # Make the plot look pretty
    ax.set_title('Candidate amplification and deletion segments')
    ax.set_xlabel('Genome position (bp)')
    ax.set_ylabel('Sample')
    ax.set_yticks(range(n_samples))
    if sampleNames is None:
        ax.set_yticklabels(range(0,n_samples,1)) # old labels for sample number
    else:
        ax.set_yticklabels(sampleNames)
    ax.margins(0.1)
    fig.tight_layout()
    ## Initialize plot picker & explore candidates
    fig.canvas.mpl_connect('pick_event', lambda event: onclick_plot_report( event, args_dict=my_data_dict ))
    print('\n\n\nNOTE:')
    print('You can click on a segment, further plots and reports will be generated with more detailed information about this region.')
    print('Depending on the amount of data this might take a couple of seconds.\n\n\n')
    return fig # need to return fig for use with callback function


def on_key(event):
    '''
    press "enter" to close main plotting window and proceed 
    '''
    # Check if the pressed key is 'enter'
    if event.key == 'enter':
        print('NOTE: Enter clicked, but you need to close all windows to resume loop!')
        plt.close('all') 


def plot_region_enrichment_per_timepoint(samples_of_tp_df, title, outpath = '', dpi = 300):
    plt_width = len(samples_of_tp_df['Timepoint'].unique()) * 0.6
    if plt_width < 4:
        plt_width = 4
    fig, ax = plt.subplots(figsize = (plt_width, 4))
    sns.boxplot(data = samples_of_tp_df, x = 'Timepoint', y = 'values', 
                boxprops={'edgecolor':'k'}, whiskerprops={'color':'k'}, capprops={'color':'k'}, medianprops={'color':'k'}, 
                linewidth = 0.75, color='white', showfliers=False, ax = ax)
    sns.stripplot(data = samples_of_tp_df, x = 'Timepoint', y = 'values', jitter=True, color='grey', alpha=0.5, ax = ax)
    ax.set_ylim(ymin = 0)
    ax.set_ylabel("Enrichment of region\nabove mean sample coverage", fontsize=12)
    ax.set_xlabel(None)
    plt.yticks(fontsize=10) 
    plt.xticks(fontsize=10, rotation = 45, ha = 'right', rotation_mode = 'anchor') 
    ax.set_title(title)
    plt.tight_layout()
    if outpath: 
        fig.savefig(outpath, dpi = dpi)
        plt.close()
    else:
        fig.show()

def plot_freq_of_mge_presence(samples_w_mge, title, outpath = '', dpi = 300):
    plt_width = len(samples_w_mge.keys()) * 0.6
    if plt_width < 4:
        plt_width = 4
    fig, ax = plt.subplots(figsize = (plt_width, 4))
    ax.bar(samples_w_mge.keys(), samples_w_mge.values(),color='grey',width=0.6)
    plt.ylabel('Percent isolates with MGE')
    plt.yticks(fontsize=10) 
    plt.xticks(fontsize=10, rotation = 45, ha = 'right', rotation_mode = 'anchor') 
    ax.set_title(title)
    plt.tight_layout()
    if outpath: 
        fig.savefig(outpath, dpi = dpi)
        plt.close()
    else:
        fig.show()


############ 
## MISSCELLANEOUS
############

def save_set_params(outpath, general_params_dict, deletion_params_dict, amplification_params_dict, save_plot_params_dict, is_interactive_analysis):
    linebreak = '\n'
    with open(f'{outpath}_parameters.txt', 'w') as fo:
        fo.write(linebreak.join([f'{param_name}\t{param}' for param_name, param in general_params_dict.items()]) + linebreak)
        fo.write(linebreak.join([f'{param_name}\t{param}' for param_name, param in deletion_params_dict.items()]) + linebreak)
        fo.write(linebreak.join([f'{param_name}\t{param}' for param_name, param in amplification_params_dict.items()]) + linebreak)
        fo.write(f'interactive\t{is_interactive_analysis}{linebreak}')
        fo.write(linebreak.join([f'{param_name}\t{param}' for param_name, param in save_plot_params_dict.items()]))