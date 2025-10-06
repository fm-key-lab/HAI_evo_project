'''
Input: Path to log samples (for mapping and picard processing)

outputs: alignment_stat.csv and deduplication_stat.csv
'''

import os
import re
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

## argument parser
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='''\
                            Generates QC stats for mapping with bowtie2 and deduplication with picard
                               ''',
                               epilog="Questions or comments? --> fenk@mpiib-berlin.mpg.de")
parser.add_argument("-s", dest="snakemake", help="Script is called from snakemake, default: True", required=False, default=True)
parser.add_argument("-p", dest="path", help="Current working directory or absolut path to subdirectories with log files. default: current working dir", required=False, default=os.getcwd())
parser.add_argument("-m", dest="mapping_path", help="Relative path to log files from bowtie2 including its mapping rates", required=False)
parser.add_argument("-d", dest="dedup_path", help="Relative path to log files from picard including deduplication metrics", required=False)
args = parser.parse_args()

## path quality check
def check_path(path_to_stats):
    if path_to_stats[-1] != '/':
        path_to_stats = path_to_stats + '/'
    return path_to_stats

## get mapping infos
def get_mapping_stats(path_to_stats):
    ## check path 
    path_to_stats = check_path(path_to_stats)
    
    ## get all files in directory
    stat_files = os.popen('ls ' + path_to_stats + '*.bowtie_stats.txt').read().split('\n')

    ## remove empty strings in file list
    stat_files = list(filter(None, stat_files))

    ## generate pattern for re.subt() to grep the info of interest in line  
    grep_pattern = '|'.join(['.*\\(', '%.*']) ## info is between '(' and '%'. Note: parenthesis need escape character '\\'
    
    ## get stats from each file 
    mapping_stat_ls = []
    for file in stat_files: 
        filename = file.split('/')[-1]
        sample_id = re.sub('_ref_.*', '', filename) ## extract sample ids (before reference genome strings), Note: Split does not work, as sample ids might have '_' in names
        ref_genome = re.sub('.*_ref_', '', re.sub('_aligned.*', '', filename)) ## extract reg genome strings (after sample id and before '_aligned*' string)
        with open(file, 'r') as fid:
            ## generate tracker to have same shapes in all entries. If one might be corrupt and one pattern is missing get an emtpy string instead
            hit_once = 0
            hit_overall = 0
            for line in fid:
                lst = line.strip()
                ## check for pattern in line 
                print(line)
                print(lst)
                if 'exactly' in lst:
                    print('in exactly')
                    read_mapping_rate_once = re.sub(grep_pattern, '', lst)
                    read_mapping_rate_once = float(read_mapping_rate_once)
                    hit_once = 1
                elif 'overall' in lst:
                    print('in overall')
                    read_mapping_rate_overall = re.sub(grep_pattern, '', lst)
                    read_mapping_rate_overall = float(read_mapping_rate_overall)
                    hit_overall = 1
            if hit_once == 0:
                print('in not once')
                read_mapping_rate_once = float('nan')
                print(f'Warning: Line with exactely one time map is missing in {file}')
            if hit_overall == 0:
                print('in not overall')
                read_mapping_rate_overall = float('nan')
                print(f'Warning: Line with overall mapping rate is missing in {file}')
        
        mapping_stat_ls.append([sample_id, ref_genome, read_mapping_rate_overall, read_mapping_rate_once])

    mapping_stat_df = pd.DataFrame(mapping_stat_ls, columns = ['sample_id', 'ref_genome', 'overall_map_rate', 'once_map_rate'])
    ## save dataframe 
    mapping_stat_df.to_csv(path_to_stats + 'alignment_stat.csv', index = False)
    
    return mapping_stat_df

def get_dedup_stats(path_to_stats):
    
    ## get all files in directory
    stat_files = os.popen('ls ' + path_to_stats + '*.pic_metrics.txt').read().split('\n')

    ## remove empty strings in file list
    stat_files = list(filter(None, stat_files))

    ## generate pattern for re.subt() to grep the info of interest in line  
    grep_pattern = '|'.join(['.*\\(', '%.*']) ## info is between '(' and '%'. Note: parenthesis need escape character '\\'

    ## get stats from each file 
    dedup_stat_ls = []

    for file in stat_files:
        filename = file.split('/')[-1]
        sample_id = re.sub('_ref_.*', '', filename) ## extract sample ids (before reference genome strings), Note: Split does not work, as sample ids might have '_' in names
        ref_genome = re.sub('.*_ref_', '', re.sub('_aligned.*', '', filename)) ## extract reg genome strings (after sample id and before '_aligned*' string)
        with open(file, 'r') as fid:
            line_grepper = 0
            for line in fid: 
                ## get line after match 
                if line_grepper == 1:
                    lst = line.strip().split('\t')
                    dedup_stat_ls.append([sample_id, ref_genome, int(lst[2]), int(lst[6]), int(lst[7])])
                    line_grepper = 0 ## set to 0 to not grep any other line in file
                ## save info if match occured 
                if 'UNPAIRED_READS_EXAMINED' in line: 
                    line_grepper = 1

    ## generate df 
    dedup_stat_df = pd.DataFrame(dedup_stat_ls, columns = ['sample_id', 'ref_genome', 'no_read_pairs', 'no_red_pair_dups', 'no_red_pair_optical_dups'])
    ## calculate percentage of duplicate reads
    dedup_stat_df['duplicate_rate'] = dedup_stat_df['no_red_pair_dups'] / dedup_stat_df['no_read_pairs'] *100
    dedup_stat_df['opt_duplicate_rate'] = dedup_stat_df['no_red_pair_optical_dups'] / dedup_stat_df['no_read_pairs'] *100
    ## save df 
    dedup_stat_df.to_csv(path_to_stats + 'deduplication_stat.csv', index = False)
    
    return dedup_stat_df


def plot_map_stats(map_stat_df, y_col, title, pdf_out):
    ## set size of figure 
    fig_width = len(set(map_stat_df['ref_genome']))/2
    if fig_width < 5:
        fig_width = 5

    sns.set_style("ticks",{'axes.grid' : True})
    fig = plt.figure(figsize=(fig_width, 8))
    
    ## plt
    bplot = sns.boxplot(data = map_stat_df, 
                        x = 'ref_genome', 
                        y = y_col,
                        showfliers = False)
    ## add points to boxplot
    bplot = sns.stripplot(data = map_stat_df, 
                        x = 'ref_genome', 
                        y = y_col, 
                        jitter=True, 
                        marker='o', 
                        alpha=0.3,
                        color='black')

    ## modfiy and label axes
    bplot.set_ylim(0, 100)
    bplot.set_ylabel("Read alignment rate to genome [%]", fontsize=12)
    bplot.set_xlabel("Reference Genomes", fontsize=12)
    bplot.set_xticklabels(bplot.get_xticklabels(), rotation=45, horizontalalignment='right')
    bplot.axes.set_title(title, fontsize=16)
    bplot.grid(True)
    fig.tight_layout()

    ## save figue
    bplot.figure.savefig(pdf_out, format='pdf')

def plot_dedup_stats(dedup_stat_df, ycol, title, pdf_out, y_label = 'Duplicates [%]', reg_line = True):
    
    sns.set_style("ticks",{'axes.grid' : True})
    ## set size of figure 
    fig = plt.figure(figsize=(8, 8))
    ## plt 
    ax = sns.regplot(data = dedup_stat_df, 
                    x = 'no_read_pairs',
                    y = ycol, 
                    ci = None,
                    scatter_kws={'color': "none",
                                'edgecolor': 'black'},
                    line_kws = {'ls': '--',
                                'color': 'red'},
                    fit_reg = reg_line)
    
    ## modify labels
    ax.set_ylabel(y_label)
    ax.set_xlabel('Number of paired reads')
    ax.set_title(title)
    ax.grid(True)
    fig.tight_layout()

    ## save figure
    fig.savefig(pdf_out, format = 'pdf')


if __name__ == '__main__':
    ## check path 
    path_to_stats = check_path(args.path)

    ## if mapping path was given or script is run from snakemake, extract data and plot
    if (args.mapping_path is not None) or (args.snakemake == True):
        ## set path if no path was stated but script was run in snakemake mode
        if args.mapping_path is None:
            mapping_path = "3-bowtie2/mapping_stats/"
        else:
            ## check path
            mapping_path = check_path(args.mapping_path)
        
        ## get mapping stats and save them
        mapping_stats = get_mapping_stats(path_to_stats + mapping_path)
        ## check if outdir exists, else create
        if not os.path.exists('pdf/'):
            os.makedirs('pdf/')

        ## plot mapping stats
        plot_map_stats(mapping_stats, 
                        'overall_map_rate', 
                        'Overall Mapping Rate', 
                        path_to_stats + 'pdf/overall_map_rate.pdf')
        plot_map_stats(mapping_stats, 
                        'once_map_rate', 
                        'Rate of all reads mapping once', 
                        path_to_stats + 'pdf/once_map_reads_rate.pdf')

    ## if dedup path was provided or script started from within snakemake, extract data and plot 
    if (args.dedup_path is not None) or (args.snakemake == True):
        ## set path if no path was stated but script was run in snakemake mode
        if args.dedup_path is None:
            dedup_path = "3-bowtie2/dedup_stats/"
        else:
            ## check path
            dedup_path = check_path(args.dedup_path)
        
        ## get mapping stats and save them
        dedup_stats = get_dedup_stats(path_to_stats + dedup_path)
        ## check if outdir exists, else create
        if not os.path.exists('pdf/'):
            os.makedirs('pdf/')

        ## plot dedup stats
        plot_dedup_stats(dedup_stats, 
                        'duplicate_rate', 
                        'Number of reads vs. rate of duplicates', 
                        path_to_stats + 'pdf/read_duplication_rate.pdf')
        plot_dedup_stats(dedup_stats, 
                        'opt_duplicate_rate', 
                        'Number of reads vs. rate of optical duplicates', 
                        path_to_stats + 'pdf/read_optical_duplication_rate.pdf', 
                        y_label = 'Optical Duplicates [%]', 
                        reg_line = False)
    
    print('Script finished!')


