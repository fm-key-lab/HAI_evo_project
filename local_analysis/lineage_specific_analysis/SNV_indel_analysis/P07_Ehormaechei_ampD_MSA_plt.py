
## extract ampD gene sequences from assemblies of E. hormaechei genomes

import glob
import gzip
import Bio.SeqIO as SeqIO
from Bio import Seq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def read_gff(gff_file):
    gff = []
    if gff_file[-2:] == 'gz':
        with gzip.open(gff_file, 'rt') as fid:
            for line in fid:
                if line.startswith('##FASTA'):
                    break
                elif line.startswith('#') or (line.strip() == ''):
                    continue
                else:
                    gff.append(line.strip().split('\t'))
    else:
        with open(gff_file, 'rt') as fid:
            for line in fid:
                if line.startswith('##FASTA'):
                    break
                elif line.startswith('#') or (line.strip() == ''):
                    continue
                else:
                    gff.append(line.strip().split('\t'))
    gff_df = pd.DataFrame(gff)   
    gff_df.columns = ['chr', 'tool', 'type', 'start', 'end', 'pval', 'strand', 'val', 'attr']
    gff_df[['start', 'end']] = gff_df[['start', 'end']].astype(int)
    return gff_df 

def find_gene_in_gff(gff, goi, interactive = False, write_annotation = False, annotation_outfile = ''):
    gois_in_gff = gff[gff['attr'].str.lower().str.contains(goi.lower())].reset_index(drop = True)
    try:
        gois_in_gff['gene'] = gois_in_gff['attr'].str.split('gene=').str[1].str.split(';').str[0]
    except:
        gois_in_gff['gene'] = ''
    gois_in_gff['ID'] = gois_in_gff['attr'].str.split('ID=').str[1].str.split(';').str[0]
    gois_in_gff['product'] = gois_in_gff['attr'].str.split('Name=').str[1].str.split(';').str[0]
    if gois_in_gff.empty:
        print(f'\nNo genes with the given identifier "{goi}" observed')
        return pd.DataFrame()
    gois_in_gff['is_chr_start'] = (gois_in_gff['start'] == 1)
    gois_in_gff['chr_end'] = gois_in_gff.apply(lambda x: gff.loc[(gff['chr'] == x['chr']) & (gff['tool'] == 'Bakta') & (gff['type'] == 'region'), 'end'].values[0], axis = 1)
    gois_in_gff['is_chr_end'] = (gois_in_gff['end'] == gois_in_gff['chr_end'])
    if len(gois_in_gff) > 1:
        print(f'\nMultiple genes with the given identifier "{goi}" observed:')
        print(gois_in_gff.loc[:,['gene','ID','chr','start','end', 'strand', 'product']].to_string(index=False))
        if interactive:
            locustag_to_extract = input('Which gene should be selected? Please specify the unique locustag (2nd column)!\n')
            gois_in_gff = gois_in_gff[gois_in_gff['ID'].str.lower() == locustag_to_extract.lower()].reset_index(drop = True)
        if write_annotation:
            gois_in_gff.loc[:,['gene','ID','chr','start','end', 'strand', 'product']].to_csv(os.path.expanduser(annotation_outfile.replace('.fasta', '.csv')), index = False)
        return gois_in_gff
    else:
        print(f'\nOne gene with the given identifier "{goi}" observed:')
        print(gois_in_gff.loc[:,['gene','ID','chr','start','end', 'strand', 'product']].to_string(index=False))
        if write_annotation:
            gois_in_gff.loc[:,['gene','ID','chr','start','end', 'strand', 'product']].to_csv(os.path.expanduser(annotation_outfile.replace('.fasta', '.csv')), index = False)
        return gois_in_gff

def extract_fasta(fasta_file, gff_of_goi, translate):
    seq_dict = {}
    if translate:
        print('\nNOTE: During translation, non-CDS will be skipped!\n')
    with gzip.open(fasta_file, 'rt') as fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            for _, entry in gff_of_goi.iterrows():
                if translate and (entry['type'] != 'CDS'):
                    continue
                if entry['chr'] == record.id: 
                    locus = f"{entry['chr']}_{entry['start']}_{entry['end']}_{entry['ID']}_{entry['gene']}"
                    if entry['strand'] == '-':
                        if entry['end'] > len(record.seq):
                            entry['end'] = len(record.seq)
                        seq_dict[locus] = str(Seq.reverse_complement(record.seq[entry['start']-1:entry['end']]))
                    else:
                        if entry['start'] < 0:
                            entry['start'] = 0
                        seq_dict[locus] = str(record.seq[entry['start']-1:entry['end']])
    return seq_dict

samplecsv = pd.read_csv('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/mapping/2024_08_25_samples.csv')
ehorm_sampleids = samplecsv.loc[samplecsv['ReferenceGenome'].str.contains('P07_Ehormaechei'), 'Sample'].values
contaminated_smpls = ['L00071_02_1', 'L00071_07_1']

## MIC results
mic_results = pd.read_csv('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/2021_09_HAP_1st_batch_UKL/Data/phenotyping/MIC_Assays/Raw/2024_07_04_MIC_Results_long_w_image_info_sorted.csv')
mic_results = mic_results[(mic_results['Species'] == 'E. hormaechei') & (mic_results['ID'] != 'P07_2244')] ## ID 2244 is non existent --> type --> remove!
mic_results = mic_results[['ID', 'AB', 'MIC']]
isolate_l = np.unique(mic_results['ID'])
mic_results['MIC'] = np.log2(mic_results['MIC'].astype(str).str.replace('>', '').astype(float))
pt_isolates = mic_results[mic_results['AB'].str.contains("Piper")]
pt_isolates['MIC_freq'] = pt_isolates['MIC'] / max(pt_isolates['MIC'])
cef_isolates = mic_results[mic_results['AB'].str.contains("Cefo")]
cef_isolates['MIC_freq'] = cef_isolates['MIC'] / max(cef_isolates['MIC'])

## get snp table and indel table to get ampD variants
snp_table = pd.read_csv('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/backup/P0007_P07_Ehormaechei-c1_240825/2025_08_09_P07_Ehormaechei-c1_240825/snp_table.csv')
indel_table = pd.read_csv('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/backup/P0007_P07_Ehormaechei-c1_240825/2025_08_09_P07_Ehormaechei-c1_240825/indel_table.csv')
snp_table.columns = [col.replace('_S2', '').replace('_S3', '') for col in snp_table.columns]
indel_table.columns = [col.replace('_S2', '').replace('_S3', '') for col in indel_table.columns]
indel_metadata_cols = ['gene_num', 'gene_num_global', 'chr', 'pos', 'chr_boundary', 'indel_ref', 'indel_anc', 'indel_alt', 'max_indel_GL_diff_to_ref', 'indel_size_gp', 'type', 'product', 'gene', 'protein_id', 'strand', 'loc1', 'loc2', 'sequence', 'note', 'locustag', 'orthologtag', 'translation', 'indel_pos_start', 'indel_pos_end', 'alt_start', 'alt_stop', 'alt_stop_pos', 'frameshift', 'mut_translation', 'mut', 'num_isolates']
smpl_cols = [col for col in indel_table.columns if (col not in indel_metadata_cols) & ('_GCA_' not in col) & ('_GCF_' not in col)]
has_ampD_variant = {id:[np.nan, np.nan] for id in smpl_cols}
for isolate in has_ampD_variant.keys():
    ## SNV 
    isolate_col = snp_table.columns[[isolate in col for col in snp_table.columns]][0]
    snv_is_ampD = (snp_table['gene'] == 'ampD')
    smpl_has_snv_at_pos = (snp_table.loc[snv_is_ampD, 'nt_anc'] != snp_table.loc[snv_is_ampD, isolate_col]) & (snp_table.loc[snv_is_ampD, isolate_col] != '?')
    if sum(snv_is_ampD & smpl_has_snv_at_pos) > 0:
        has_ampD_variant[isolate][0] = snp_table.loc[(snv_is_ampD & smpl_has_snv_at_pos), 'nt_pos'].values
    ## Indel
    indel_is_ampD = (indel_table['gene'] == 'ampD')
    smpl_has_indel_at_pos = (indel_table.loc[indel_is_ampD, 'indel_ref'] != indel_table.loc[indel_is_ampD, isolate])
    if sum(indel_is_ampD & smpl_has_indel_at_pos) > 0:
        has_ampD_variant[isolate][1] = indel_table.loc[(indel_is_ampD & smpl_has_indel_at_pos), 'indel_pos_start'].str.replace(r'\[|np.int64\(|\)|\]|', '', regex = True).str.split(', ').str[0].values

## assembled genomes 
genome_dirs = glob.glob(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/metadata/reference_genomes/isolates_250507/*'))

outfasta_ampE = '~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2023/small_analyses/P07_Ehormaechei_ampDE_region/fasta/P07_Ehormaechei_denovoassembly_ampEminus750bp.fasta'
outfasta_ampD = '~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2023/small_analyses/P07_Ehormaechei_ampDE_region/fasta/P07_Ehormaechei_denovoassembly_ampD.fasta'

ampE_seqs = {}
ampD_seqs = {}
ampD_at_chr_boundary = {}

## sample ids with new identified indels
smpls_to_check = ['P07_1047', 'P07_1540', 'P07_1049', 'P07_1545', 'P07_1553', 'P07_1046', 'P07_1201', 'P07_1149']
gff_of_new_indel_smpls = {}
for genome_dir in genome_dirs:
    sampleid = genome_dir.split('/')[-1][:-7]
    if sampleid in ehorm_sampleids:
        gff = read_gff(f'{genome_dir}/bakta_annot.gff.gz')
        ## find ampE in gff as anchor of ampD identificatiuon 
        gff_ampE = find_gene_in_gff(gff, 'ampE')
        ## extend region to include ampD gene
        gff_ampE.loc[gff_ampE['strand'] == '+', 'start'] = gff_ampE.loc[gff_ampE['strand'] == '+', 'start'] - 750
        gff_ampE.loc[gff_ampE['strand'] == '-', 'end'] = gff_ampE.loc[gff_ampE['strand'] == '-', 'end'] + 750

        fasta_seq = extract_fasta(f'{genome_dir}/genome.fasta.gz', gff_ampE, translate = False)
        ampE_seqs[sampleid] = list(fasta_seq.values())
        ## search for ampD specifically
        gff_ampD = find_gene_in_gff(gff, 'ampD')
        if sampleid in smpls_to_check:
            gff_of_new_indel_smpls[sampleid] = gff_ampD[['chr', 'start', 'end', 'strand', 'chr_end', 'is_chr_start', 'is_chr_end']]
        fasta_seq = extract_fasta(f'{genome_dir}/genome.fasta.gz', gff_ampD, translate = False)
        ampD_seqs[sampleid] = list(fasta_seq.values())
        if not gff_ampD.empty:
            ampD_at_chr_boundary[sampleid] = list((gff_ampD['is_chr_start'] | gff_ampD['is_chr_end']).values)
        else:
            ampD_at_chr_boundary[sampleid] = np.nan

with open(os.path.expanduser(outfasta_ampE), 'w') as fo:
    for id, seq in ampE_seqs.items():
        if id not in contaminated_smpls:
            for seq_id, seq in enumerate(seq):
                fo.write(f'>{id}_seqid{seq_id+1}\n{seq}\n')
with open(os.path.expanduser(outfasta_ampD), 'w') as fo:
    for id, seq in ampD_seqs.items():
        if id not in contaminated_smpls:
            for seq_id, seq in enumerate(seq):
                fo.write(f'>{id}_seqid{seq_id+1}\n{seq}\n')

## sequences were aligned via AliView via mafft_conda `--out TEMP_OUT_FILE CURRENT_ALIGNMENT_FASTA`

## read in aligned sequences and extract from genomes with mutiple sequences the ones which are clustering the best 
sample_seq_dict = {}
with open(os.path.expanduser(outfasta_ampD.replace('.fasta', '.aligned.fasta')), 'r') as fasta:
    for record in SeqIO.parse(fasta, 'fasta'):
        sampleid = '_'.join(record.id.split('_')[:-1])
        if sampleid in sample_seq_dict.keys():
            sample_seq_dict[sampleid].append(str(record.seq))
        else:
            sample_seq_dict[sampleid] = [str(record.seq)]

## get all samples with a unique sequence and get major representative sequence from that 
major_seq_l = []
for seq_of_sample in sample_seq_dict.values():
    if len(seq_of_sample) == 1:
        major_seq_l.append(np.array([nt for nt in seq_of_sample[0]]))

major_alleles = pd.DataFrame(major_seq_l).mode(axis = 0).to_numpy()[0]

## add samples with empty sequences back in
for sampleid, seq in ampD_seqs.items():
    if (seq == []) & (sampleid not in sample_seq_dict.keys()):
        sample_seq_dict[sampleid] = [''.join(['-' for _ in major_alleles])]

major_seq_ampD = {}
major_seq_ampD_chr_boundaries = {}
for sample, seq_of_sample in sample_seq_dict.items():
    if len(seq_of_sample) == 1:
        major_seq_ampD[sample] = np.array(['.' if (nt == mnt) & (mnt != '-') else nt for mnt, nt in zip(major_alleles, seq_of_sample[0])])
        if np.isnan(ampD_at_chr_boundary[sampleid]):
            major_seq_ampD_chr_boundaries[sample] = np.nan
        else:
            major_seq_ampD_chr_boundaries[sampleid] = ampD_at_chr_boundary[sampleid][0]
    else:
        min_dist_to_major_allele = len(major_alleles)
        for seq_id, seq in enumerate(seq_of_sample):
                
            aligned_seq = np.array(['.' if nt == mnt else nt for mnt, nt in zip(major_alleles, seq)])
            dist_to_major_allele = np.sum(aligned_seq != '.')
            if dist_to_major_allele < min_dist_to_major_allele:
                selected_seq = aligned_seq
                min_dist_to_major_allele = dist_to_major_allele
                major_seq_ampD_chr_boundaries[sample] = ampD_at_chr_boundary[sample][seq_id]
        major_seq_ampD[sample] = np.array(['.' if (nt == mnt) & (mnt != '-') else nt for mnt, nt in zip(major_alleles, seq_of_sample[0])])


heatmap_df = pd.DataFrame(major_seq_ampD, index = np.arange(len(major_alleles))).T
heatmap_df = heatmap_df.loc[[idx for idx in heatmap_df.index if (idx[-3:-1] != '_S') | idx.startswith('L00071_02_1') | idx.startswith('L00071_07_1')]]
## remove contaminated samples (have been resequenced later!)
heatmap_df = heatmap_df[~heatmap_df.index.isin(contaminated_smpls)]
for col in heatmap_df.columns:
    if (len(np.unique(heatmap_df[col])) == 1) & (np.unique(heatmap_df[col])[0] == '-'):
        heatmap_df = heatmap_df.drop(col, axis = 1)

df_long = heatmap_df.reset_index().melt('index')
df_long.loc[df_long['value'] == '.', 'value'] = float('nan')
df_long = df_long.sort_values('index').reset_index(drop = True)
color_dict = {idx:color for idx, color in zip('-atcg', ['lightgrey', 'green', 'red', 'blue', 'yellow'])}

## get some statistics:
major_allele_len = sum([nt!='-' for nt in major_alleles])
ampd_stat = {}
maf_non_gap_pos = np.where(major_alleles != '-')
for sid, seq in heatmap_df.iterrows():
    length_major_allele = sum([nt=='.' for nt in seq])
    length_diff_allele = sum([nt not in ['.', '-'] for nt in seq])
    alignment_start = np.where(seq == '.')
    if np.shape(alignment_start)[1] > 0:
        alignment_start = np.min(alignment_start)
        first_difference = np.where(seq[alignment_start:] != '.') ## account for contig breaks and therefore alignment issues
        if np.shape(first_difference)[1] > 1: ## require at least 2 different nts, otherwise it is an SNV
            seq_pos_variant_first = first_difference[0]+alignment_start
            idx_pos_variant_first = np.where(~((major_alleles[seq_pos_variant_first] == '-') & (seq[seq_pos_variant_first] == '-'))) ## check if either the sample has an insertion or the major allele is covered --> this allows to "jump" over insertions generated by other samples and stick to the true nt count!
            if np.shape(idx_pos_variant_first)[1] > 1: ## require at least 2 different nts, otherwise it is an SNV
                first_difference = seq_pos_variant_first[np.min(idx_pos_variant_first)]
            else:
                first_difference = len(seq)
        else:
            first_difference = len(seq)-1
        first_difference = first_difference - sum(major_alleles[:first_difference] == '-') ## remove any gap which was introduced to the alignment to get true position
        
        ## get also the last difference before the gene ends
        last_non_difference = np.where(seq[alignment_start:] == '.') ## ≥ 20 is for contig break issues at start 
        if np.shape(last_non_difference)[1] > 1: ## require at least 2 different nts, otherwise it is an SNV
            seq_pos_variant_last = last_non_difference[0]+alignment_start
            idx_pos_variant_last = np.where(~((major_alleles[seq_pos_variant_last] == '-') & (seq[seq_pos_variant_last] == '-'))) ## check if either the sample has an insertion or the major allele is covered --> this allows to "jump" over insertions generated by other samples and stick to the true nt count!
            if np.shape(idx_pos_variant_last)[1] > 1: ## require at least 2 different nts, otherwise it is an SNV
                last_non_difference = seq_pos_variant_last[np.max(idx_pos_variant_last)]
            else:
                last_non_difference = len(seq)
        else:
            last_non_difference = len(seq)-1
        last_non_difference = last_non_difference - sum(major_alleles[:last_non_difference] == '-') ## remove any gap which was introduced to the alignment to get true position
    else:
        alignment_start = np.nan
        first_difference = np.nan
        last_non_difference = np.nan

    ampd_stat[sid] = [length_major_allele, length_diff_allele, length_major_allele/major_allele_len, length_diff_allele/major_allele_len, first_difference, last_non_difference]


ampd_stat_df = pd.DataFrame(ampd_stat, index = ['num_major_allele', 'num_alt_allele', 'freq_major_allele', 'freq_alt_allele', 'first_difference', 'last_non_difference']).T
ampd_stat_df = ampd_stat_df.sort_values(['freq_major_allele', 'freq_alt_allele'], ascending = False)

## add positions of snvs and indels 
for smpl, (snv_pos, indel_pos) in has_ampD_variant.items():
    if np.shape(snv_pos) != ():
        ampd_stat_df.loc[smpl, 'snv_pos'] = '; '.join(snv_pos.astype(int).astype(str))
    if np.shape(indel_pos) != ():
        ampd_stat_df.loc[smpl, 'indel_pos'] = '; '.join(indel_pos.astype(int).astype(str))

ampd_stat_df.to_csv(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2023/small_analyses/P07_Ehormaechei_ampDE_region/P07_Ehormaechei_ampD_additional_indels.csv'))

fig, ax = plt.subplots(figsize = (30, 7.5), nrows = 2, sharex = True, gridspec_kw={'height_ratios':[2, 1], 'hspace': 0.1})
bottom_tracker = np.zeros(len(ampd_stat_df), dtype = 'float')
for col, mdata in {'freq_major_allele': ['Major', 'darkgray'], 'freq_alt_allele': ['Alternative', 'white']}.items():
    l, c = mdata
    ax[0].bar(x = ampd_stat_df.index, height = ampd_stat_df[col], bottom=bottom_tracker, width = 1, color = c, label = l, edgecolor = 'k', lw = 0.25)
    bottom_tracker += ampd_stat_df[col]
ax[0].set_ylabel('Frequency of covered alleles')
ax[0].legend(title = 'Allele', bbox_to_anchor = (1, 0.5), loc = 'center left')
ax[1].scatter(x = pt_isolates['ID'], y = pt_isolates['MIC'], c = 'r', edgecolor = 'k', marker = 'X', label = 'Piperacillin/\nTazobactam')
ax[1].scatter(x = cef_isolates['ID'], y = cef_isolates['MIC'], c = 'purple', edgecolor = 'k', marker = 'X', label = 'Cefotaxime')
ax[1].legend(title = 'Antibiotic', bbox_to_anchor = (1, 0.5), loc = 'center left')
ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation = 90, ha = 'center')
ax[1].set_xlim(-1, len(ampd_stat_df)+1)
ax[1].set_ylabel('log2(MIC)')
fig.suptitle('AmpD length in de novo assembled genomes')
fig.savefig(os.path.expanduser('~/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/MIC_analysis/2025_05_02_ampD_genelength_vs_MIC.pdf'), bbox_inches = 'tight')

ampd_length = {sid: sum([True for nt in seq if nt == '.']) for sid, seq in heatmap_df.iterrows()}
print(ampd_length)


## print some statistics
unique_indels = ampd_stat_df.groupby('first_difference').agg({'num_major_allele':'size',
                                                              'indel_pos': lambda x: '; '.join(np.unique([val for val in x if val == val])),
                                                              'snv_pos': lambda x: '; '.join(np.unique([val for val in x if val == val]))
                                                              }).reset_index()
unique_indels = unique_indels[unique_indels['first_difference'] < major_allele_len-1]
unique_indels['first_difference'] = unique_indels['first_difference'].astype(int) + 1 ## convert to 1-based positions (note, indel_pos is 1-based!)
unique_indels.rename({'num_major_allele': 'smpl_cnts'}, axis = 1, inplace = True)
num_samples_w_indel = sum(unique_indels["smpl_cnts"])
freebayes_identified_indels = unique_indels[unique_indels['indel_pos'] != '']
num_samples_annotated_indels = num_samples_w_indel - len(ampd_stat_df.loc[~ampd_stat_df['indel_pos'].isna() & ampd_stat_df['first_difference'].isna()])
msa_indels = unique_indels[unique_indels['indel_pos'] == '']
num_smpls_msa = sum(msa_indels['smpl_cnts'])
smpls_w_indels = ampd_stat_df[ampd_stat_df['first_difference'] < major_allele_len-1]
smpls_w_indels['first_difference'] = smpls_w_indels['first_difference'].astype(int) + 1 ## convert to 1-based positions (note, indel_pos is 1-based!)
indel_to_smpl_dict = smpls_w_indels.groupby('first_difference').apply(lambda x: x.index.tolist()).to_dict()


print(f'Identified {len(unique_indels)} indels within {num_samples_w_indel} samples. This adds to {sum([len(snp_table[snv_is_ampD]), len(unique_indels)])} mutations ({len(snp_table[snv_is_ampD])} incl. SNVs) within ampD')
print(f'{len(freebayes_identified_indels)} of those indels have been identified by freebayes with {num_samples_annotated_indels} have been found by the MSA as well. (NOTE: Not necessarily the same number of indels due to alignment issues!)')

print(f'{len(msa_indels)} of those indels have been identified by MSA only within {num_smpls_msa} samples.')

print(f'Unique indel positions are (1-based positions):\n{unique_indels.to_string(index = False)}')
print(f'\nSmpls with indels (1-based positions):')
for indel, smpls in indel_to_smpl_dict.items():
    print(f'{int(indel)}: {", ".join(smpls)}')

###
## OUTPUT:
"""
Identified 8 indels within 43 samples. This adds to 11 mutations (3 incl. SNVs) within ampD
5 of those indels have been identified by freebayes with 43 have been found by the MSA as well. (NOTE: Not necessarily the same number of indels due to alignment issues!)
3 of those indels have been identified by MSA only within 8 samples.
Unique indel positions are (1-based positions):
 first_difference  smpl_cnts indel_pos snv_pos
               48          1        47        
               51          4        47        
              207         28       206        
              272          1       265        
              279          6                  
              348          1       345        
              391          1                  
              480          1                  

Smpls with indels (1-based positions):
48: P07_1548
51: P07_1583, P07_1155, P07_1585, P07_1160
207: P07_1084, P07_1546, P07_1086, P07_1094, P07_1564, P07_1550, P07_1534, P07_1561, P07_1087, P07_1541, P07_1093, P07_1088, P07_1562, P07_1090, P07_1082, P07_1542, P07_1080, P07_1092, P07_1052, P07_1089, P07_1563, P07_1085, P07_1557, P07_1095, P07_1050, P07_1091, P07_1565, P07_1083
272: P07_1544
279: P07_1047, P07_1540, P07_1049, P07_1545, P07_1553, P07_1046
348: P07_1555
391: P07_1201
480: P07_1149
"""