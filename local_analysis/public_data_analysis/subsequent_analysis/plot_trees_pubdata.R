## Load libraries
library(ape) ## to modify tree file
library(tidyverse)
library(Biostrings)
library(stringr)
library(ggplot2)
library(ggtree)
library(ggstar) ## to plot stars as shape
library(ggnewscale)  # for multiple fill scales in ggplot
library(RColorBrewer)
# library(treeio)

## Functions
get_descendant_tips <- function(tree, node) {
  # Find all descendants
  get_desc <- function(n) {
    children <- tree$edge[tree$edge[, 1] == n, 2]
    if (length(children) == 0) return(integer(0))
    c(children, unlist(lapply(children, get_desc)))
  }
  
  descendants <- get_desc(node)
  
  # Return tip labels only
  tip_ids <- descendants[descendants <= length(tree$tip.label)]
  tree$tip.label[tip_ids]
}

## Translate study and source abbreviations 
source_l = list('urine' = 'Urinary tract', ## Includes: urine, urinary, urethral
                'rectum' = 'Gastrointestinal tract', ## Includes: rectum, rectal, stool, anal, feces
                'hos_ass' = 'Hospital-associated', ## Includes: patient or hospital, clinical
                'blood' = 'Bloodstream', ## Includes: blood
                'wound' = 'Wound', ## Includes: wound
                'lung' = 'Respiratory tract', ## Includes: sputum, respiratory, tracheal, bronchoalveolar, lavage, bronchial, aspirate, bal, trach asp
                'other_body' = 'Other body sources', ## Includes: body, bodily, non-wound, tissue, abdomen, bone, foot, skin, endocervical, leg, ulcer, breast, abscess, abdominal, gall, bile, peritoneal, throat, groin
                'env' = 'Environment', ## Includes: environmental, wastewater, sink, drain, shower, bathroom
                'other' = 'Other/Unclassified' ## Includes: all remaining and missing data
) 
study_l = list('ford22' = 'Forde et al. (2022)',
               'van20' = 'Van Duin et al. (2020)',
               'chav16' = 'Chavda et al. (2016)',
               'cdc' = 'CDC HAI-Seq',
               'stoe21' = 'Stoesser et al. (2024)' ## note--> now published
)

## Okabe ito extended by 9th color and with grey
source_colors = list('Respiratory tract' = '#56B4E9',
                     'Urinary tract' = '#F0E442',
                     'Bloodstream' = '#AA4499',
                     'Wound' = '#CC79A7',
                     'Gastrointestinal tract' = '#E69F00',
                     'Other body sources' = '#009E73',
                     'Hospital-associated' = '#D55E00',
                     'Environment' = '#0072B2',
                     'Other/Unclassified' = '#DDDDDD'
) 

study_colors = list('CDC HAI-Seq' = '#00767B',
                    'Chavda et al. (2016)' = '#60BCE9',
                    'Forde et al. (2022)' = '#DEE6E7',
                    'Stoesser et al. (2024)' = '#F9D576',
                    'Van Duin et al. (2020)' = '#E94C1F'
)
nonsyn_syn_color = list('Nonsynonymous' = '#454545',
                        'Synonymous' = '#cfcfcf')


## NOTE:
# One need to let at first run refbased plottings through to get all the first observational dates of that species and then the denovo dates!
outpath <- '/Users/ad_lemur/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/phylogenetic_trees'

analysis_path <- "/Users/ad_lemur/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2023/public_data/2023_12_analysis/analysis/2025_08_25_P07_Ehormaechei-c1_240825_gubbinsdefault_ambig10"
species_observation_range <- c()

patient_species <- gsub('.*/', '', analysis_path)
patient <- 'Public data'
species <- 'Ehormaechei-c1_240825'

## Load tree
tree_files <- list.files(path = paste0(analysis_path, '/raxml_bootstrapped_tree1000'), pattern = ".*1000bs_ML_s123_with_bootstraps.raxml.support", full.names = TRUE)
tree <- read.tree(tree_files[length(tree_files)])

## ensure tree is rooted
outgroup_tip <- grep("^Equasihorm", tree$tip.label, value = TRUE)
rerooted_tree <- root(tree, outgroup = outgroup_tip, resolve.root = TRUE)
rerooted_tree <- rotateConstr(rerooted_tree, sort(rerooted_tree$tip.label, decreasing = TRUE)) ## sort internal branches by sampleid

## remove reference from tree
tips_to_remove <- grep("_ref$", rerooted_tree$tip.label, value = TRUE)
tips_to_remove <- append(tips_to_remove, outgroup_tip)
# rerooted_tree_wo_ref <- rerooted_tree
rerooted_tree_wo_ref <- drop.tip(rerooted_tree, tips_to_remove)

## Extract metadata from tree
split_tip_labels <- strsplit(rerooted_tree_wo_ref$tip.label, "_")
study <- sapply(split_tip_labels, function(x) x[1])
accession_id <- sapply(split_tip_labels, function(x) x[length(x)])
source <- sapply(split_tip_labels, function(x) gsub('^eq_|^20_', '', paste(x[3:length(x)-1], collapse = '_')))

tree_annotation_heatmap_df <- tibble(
  tip_label = rerooted_tree_wo_ref$tip.label,
  study = gsub('^van', 'van20', study),
  source = source,
  accession_id = accession_id
) %>%
  as.data.frame()

## map the data to correct names
map_source_vec <- unlist(source_l)
tree_annotation_heatmap_df$source <- map_source_vec[tree_annotation_heatmap_df$source]
tree_annotation_heatmap_df$source <- factor(tree_annotation_heatmap_df$source, levels = names(source_colors))
map_study_vec <- unlist(study_l)
tree_annotation_heatmap_df$study <- map_study_vec[tree_annotation_heatmap_df$study]


## set everything for outgroup to NA
rows_to_na <- grepl("^Equasihorm", tree_annotation_heatmap_df$tip_label)

## get nodes to highlight 
tips_to_highlight_branch <- c()
if (length(tips_to_highlight_branch) > 0){
  highlight_clade <- MRCA(rerooted_tree_wo_ref, tips_to_highlight_branch)
}

## get nodes to annotate
snp_table_file <- paste0(analysis_path, "/snp_table.csv")
snp_table <- read.csv(snp_table_file, header=TRUE, row.names = NULL, sep = ',')

annotation_mutation_file <- list.files(path = analysis_path, pattern = ".*annotation_mutations.csv", full.names = TRUE) 
annotation_mutation <- read.csv(annotation_mutation_file[length(annotation_mutation_file)], header=TRUE, row.names = NULL, sep = ',')
annotation_mutation <- annotation_mutation[, c('chr', 'pos', 'nt_pos', 'AA_gt_ref', 'AA_gt_anc', 'AA_gt_alt')]
snp_table <- left_join(snp_table, annotation_mutation, by = c('chr', 'pos', 'nt_pos'))

## get ancestral seq from treetime reconstruction for proper rooting
anc_reconstruction_fasta <- paste0(analysis_path, "/node_sequence_inference/ancestral_sequences.fasta")
anc_seq <- readDNAStringSet(anc_reconstruction_fasta)$NODE_0000001
anc_seq <- unlist(strsplit(as.character(anc_seq), split = ""))
snp_table$nt_anc <- anc_seq

snp_table_fimZ <- snp_table[(snp_table$gene == 'fimZ'), ]
rownames(snp_table_fimZ) <- NULL
metadata_colnames <- c(colnames(snp_table_fimZ)[1:15], c('AA_gt_ref', 'AA_gt_anc', 'AA_gt_alt'))
isolate_colnames <- colnames(snp_table_fimZ)[!(colnames(snp_table_fimZ) %in% metadata_colnames)]

print(snp_table_fimZ[,c('AA_gt_ref', 'AA_gt_anc', 'AA_gt_alt')])

tips_to_mark_mut <- list()
for (i in 1:nrow(snp_table_fimZ)){
  obs_nts <- snp_table_fimZ[i, c('nt_anc', 'nt_ref', isolate_colnames)]
  unique_nts <- unique(na.omit(unlist(obs_nts)))
  ## get samples which contain mutation
  smpl_with_mut <- names(obs_nts)[!(obs_nts == obs_nts$nt_anc | obs_nts == '?')] 
  ## check if mutation is polarized in snp table 
  if (obs_nts$nt_anc != obs_nts$nt_ref){
    print('Mutation is not polarized well. Please fix!')
    break
  }
  ## clean mutations ([, ] and ') as well as empty dots
  mutation <- gsub("\\[|\\]|'", "", snp_table_fimZ[i, 'muts'])       # remove brackets and quotes
  mutation <- gsub("\\b\\.\\b", "", mutation)                   # remove dots
  mut_name <- ifelse(mutation == '', paste0(snp_table_fimZ[i, 'AA_gt_alt'], snp_table_fimZ[i, 'aa_pos']), mutation)
  tips_to_mark_mut[as.character(mut_name)] <- list(sort(smpl_with_mut))
}

####################
## plot initial tree 
####################
max_tree_depth <- max(node.depth.edgelength(rerooted_tree_wo_ref))
number_tips <- length(rerooted_tree_wo_ref$tip.label)
## heatmap specific variables
colnames_offset <- number_tips/250
width = 0.1
width_offset = max_tree_depth*(width*1.1)
plot_offset = 0
open_angle = 20
tree_rotation = 90+open_angle ## degree how much the tree is rotated
tree_scale_loc = round((360-tree_rotation)*number_tips/(360-open_angle), 0)

## plot tree
tree_plt <- ggtree(rerooted_tree_wo_ref, size = 0.8, layout="fan", open.angle = open_angle) +
  geom_tiplab(align = TRUE, linetype = "dotted", size = 0, linesize = 0.2) +
  geom_rootedge(max_tree_depth/50, size = 0.8) + 
  geom_treescale(x = max_tree_depth*1.4, y = tree_scale_loc) + ## add tree scale
  geom_point2(aes(subset = !isTip, fill = as.numeric(label)), shape = 21, color = 'black', size = 2) +
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap value") + ## add bootstrap values
  ggtitle(paste(patient, species))

  
####################
## annotate tree 
####################
tree_plt_data <- tree_plt$data


## get nodes to mark via a star
if (length(tips_to_mark_mut) > 0){
  ## check how many nodes need to be annotated multiple times and how often
  node_with_multi_annotation <- list()
  mut_to_node_list <- list()
  for (mut in names(tips_to_mark_mut)){
    ## subset list of isolates just to those which are present in tree (required to avoid errors); Note: If wrong tree loaded, then some mutations might not be annotated!!
    isolates_on_node = tips_to_mark_mut[[mut]][tips_to_mark_mut[[mut]] %in% rerooted_tree_wo_ref$tip.label]
    mark_mut_node <- MRCA(rerooted_tree_wo_ref, isolates_on_node)
    mark_mut_node_str <- as.character(mark_mut_node)
    
    ## verify the number of tips on node (in case of homoplaies to resolve manually)
    tips_on_node <- get_descendant_tips(rerooted_tree_wo_ref, mark_mut_node)
    if (length(tips_on_node) > length(isolates_on_node)){
      print(paste0('Identified for mutational position ', mut, ' ', length(tips_on_node), ' tips, but had ', length(isolates_on_node), ' mutated samples only'))
    }
    ## manual correction of homoplasy
    if (mut %in% c('S38')){
      isolates_aa38 <- list(c("cdc_eq_lung_24749090", "cdc_eq_other_25176724", "cdc_eq_urine_16893795", "cdc_eq_urine_16893810", "cdc_eq_urine_24066631", "cdc_eq_urine_26371314", "stoe21_hos_ass_18031499", "stoe21_hos_ass_18031642", "stoe21_hos_ass_18032051"),
                          c("cdc_eq_urine_19034561", "cdc_eq_urine_19572613")
      )
      mut_to_node_list[[mut]] <- c()
      for (isolates in isolates_aa38){
        mark_mut_node <- MRCA(rerooted_tree_wo_ref, isolates)
        mark_mut_node_str <- as.character(mark_mut_node)  
        if (!(mark_mut_node %in% names(node_with_multi_annotation))){
          node_with_multi_annotation[[mark_mut_node_str]] <- c(mut)
        }else{
          node_with_multi_annotation[[mark_mut_node_str]] <- c(node_with_multi_annotation[[mark_mut_node_str]], mut)
        }
        mut_to_node_list[[mut]] <- c(mut_to_node_list[[mut]], mark_mut_node)
      }
    }else{
      ## get a list of nodes with multiple annotations
      if (!(mark_mut_node %in% names(node_with_multi_annotation))){
        node_with_multi_annotation[[mark_mut_node_str]] <- c(mut)
      }else{
        node_with_multi_annotation[[mark_mut_node_str]] <- c(node_with_multi_annotation[[mark_mut_node_str]], mut)
      }
      mut_to_node_list[[mut]] <- mark_mut_node  
    }
    
    
  }
  for (mut in names(tips_to_mark_mut)){
    for (mut_node in mut_to_node_list[[mut]]){
      # Find midpoint of branch given the node
      derived <- tree_plt_data[tree_plt_data$node == mut_node, ]
      ancestral <- tree_plt_data[tree_plt_data$node == derived$parent, ]
      ## get the correct x to annotate freely the mutation
      num_muts_on_branch <- length(node_with_multi_annotation[[as.character(mut_node)]])
      id_in_list <- which(mut == node_with_multi_annotation[[as.character(mut_node)]])
      x_fraction <- (id_in_list/(1+num_muts_on_branch)) ## get the fraction where to place the star on the branch 
      mid_x <- sum(c(ancestral$x, derived$x)) * x_fraction
      ## mark mutations
      mut_is_syn <- ifelse(grepl("[0-9]$", mut), 'Synonymous', 'Nonsynonymous') ## check if last character is digit --> syn
      tree_plt <- tree_plt +
        geom_star(x = mid_x, y = min(derived$y), fill = nonsyn_syn_color[[mut_is_syn]], size = 7.5, alpha = 0.8) +
        annotate('text', x = mid_x, y = min(derived$y)+2, label = mut, size = 3)
    }
  }
  tree_plt <- tree_plt + 
    scale_color_identity(guide = "legend", labels = names(nonsyn_syn_color), breaks = unname(unlist(nonsyn_syn_color)),name = "Mutation type") +
    theme(legend.position = "right")
}

## generate heatmap dataframes 
source_df <- data.frame('Source' = tree_annotation_heatmap_df$source) %>%
  `rownames<-`(tree_annotation_heatmap_df$tip_label)

study_df <- data.frame('Study' = tree_annotation_heatmap_df$study) %>%
  `rownames<-`(tree_annotation_heatmap_df$tip_label)


tree_plt <- rotate_tree(tree_plt, angle = tree_rotation)
## add source heatmap 
tree_plt <- tree_plt + new_scale_fill()
tree_plt <- gheatmap(tree_plt, study_df, 
                     offset = plot_offset, width = width,
                     color = 'black',
                     colnames_angle=0, colnames_position = "top", 
                     colnames_offset_y = colnames_offset, hjust = 1) +
  scale_fill_manual(values = scales::alpha(study_colors, 0.6), name="Study", drop = FALSE)
  
## add new scale and days since observation heatmap 
tree_plt <- tree_plt + new_scale_fill()
tree_plt <- gheatmap(tree_plt, source_df, 
                     offset=plot_offset+width_offset, width=width,
                     color = 'black', 
                     colnames_angle=0, colnames_position = "top", 
                     colnames_offset_y = colnames_offset, hjust = 1) +
  scale_fill_manual(values = scales::alpha(source_colors, 0.8), name="Source", drop = FALSE) 

print(tree_plt)

## save plot
pdf(paste0(outpath, "/", gsub(' ', '_', patient), '-', species, '.pdf'), 
    width = 15, height = 15, paper = "special", onefile = F)
print(tree_plt)
dev.off()