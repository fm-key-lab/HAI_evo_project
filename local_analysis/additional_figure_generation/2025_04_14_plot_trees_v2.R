## Load libraries
library(ape) ## to modify tree file
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggtree)
library(ggstar) ## to plot stars as shape
library(ggnewscale)  # for multiple fill scales in ggplot
library(RColorBrewer)

## NOTE:
# One need to let at first run refbased plottings through to get all the first observational dates of that species and then the denovo dates!
use_ml_tree = FALSE
analysis_types <- c("refbased", "denovo")
outpath <- '/Users/ad_lemur/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/phylogenetic_trees'

for (analysis_type in analysis_types){
  
  if (analysis_type == "refbased"){
    analysis_path <- "/Users/ad_lemur/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/refbased_2024_08/analysis"
    species_observation_range <- c()
    all_STs <- c("23", "45", "58", "95", "97", "117", "131", "4260", "NF") ## manually defined upfront after running it once before so fixed colors per ST are used
    palette_colors <- brewer.pal(n = length(all_STs)+1, name = "Set3")
    palette_colors <- c(palette_colors[length(palette_colors)], palette_colors[3:length(palette_colors)-1]) ## remove first color which is similar to location color!
    color_map <- setNames(palette_colors, all_STs)
    color_map <- c(color_map, 'NA' = 'white') ## add NA value
    
  }else if(analysis_type == "denovo"){
    analysis_path <- "/Users/ad_lemur/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/backup"
    species_observation_range <- list(
      "P07 A. pittii" = c("2021-07-27", "2021-08-17"),
      "P07 B. thetaiotaomicron" = c("2021-08-10", "2021-08-26"),
      "P07 E. coli" = c("2021-06-29", "2021-07-20"),
      "P07 E. hormaechei" = c("2021-07-27", "2021-09-21"),
      "P07 P. aeruginosa" = c("2021-07-16", "2021-08-17"),
      "P10 K. michiganensis" = c("2021-07-26", "2021-08-02"),
      "P21 E. coli" = c("2021-09-02", "2021-10-19"),
      "P21 P. mirabilis" = c("2021-09-02", "2021-11-02"),
      "P21 S. aureus" = c("2021-09-02", "2021-09-08"))
  }
  species_paths <- list.files(path = analysis_path, pattern = "^P00.*", full.names = TRUE)
  
  ## NOTE: FROM S aureus the blood culture positive strain will be removed due to contamination!
  remove_contamination_isolate <- c("P0021_Saureus_ASM1342v1"=c("P0021_TP12_210915_D_Blood_ST5__P21_MIBI03"))
  
  for (species_path in sort(species_paths)){
    patient_species <- gsub('.*/', '', species_path)
    patient <- strsplit(patient_species, '_')[[1]][1]
    patient <- paste0(substring(patient, 1,1), substring(patient, 4))
    if (analysis_type == "refbased"){
      species_label <- strsplit(patient_species, '_')[[1]][2]
      species <- paste0(substring(species_label, 1,1), '. ', substring(species_label, 2))
      species <- gsub('YT3$', '', species) ## applies to Ehorm only 
    }else if(analysis_type == "denovo"){
      species_label <- strsplit(patient_species, '_')[[1]][3]
      species <- paste0(substring(species_label, 1,1), '. ', substring(species_label, 2))
    }
    
    most_recent_analysis_path <- sort(list.files(path = species_path, pattern = "", full.names = TRUE))
    if (length(most_recent_analysis_path) == 0) {
      next
    }else{most_recent_analysis_path <- most_recent_analysis_path[[length(most_recent_analysis_path)]]}
    ## Load tree
    if (use_ml_tree == TRUE){
      #tree_files <- list.files(path = most_recent_analysis_path, pattern = ".*.bestTree", full.names = TRUE)
      tree_files <- list.files(path = most_recent_analysis_path, pattern = ".*.bestTreeCollapsed", full.names = TRUE)
      tree_type = 'MLtree'
    }else{
      tree_files <- list.files(path = most_recent_analysis_path, pattern = ".*_TP0.nwk.tree", full.names = TRUE)
      tree_type = 'PStree'
    }
    if (length(tree_files) == 0) {
      next
    }else{tree <- read.tree(tree_files[length(tree_files)])}
      
    ## ensure tree is rooted
    outgroup_tip <- grep("^OG", tree$tip.label, value = TRUE)
    rerooted_tree <- root(tree, outgroup = outgroup_tip, resolve.root = TRUE)
    rerooted_tree <- rotateConstr(rerooted_tree, sort(rerooted_tree$tip.label, decreasing = TRUE)) ## sort internal branches by sampleid
    
    ## remove reference from tree
    tips_to_remove <- grep("_ref$", rerooted_tree$tip.label, value = TRUE)
    tips_to_remove <- append(tips_to_remove, outgroup_tip)
    if (patient_species %in% names(remove_contamination_isolate)){
      tips_to_remove <- append(tips_to_remove, unname(remove_contamination_isolate[patient_species]))
    }
    rerooted_tree_wo_ref <- drop.tip(rerooted_tree, tips_to_remove)
    
    ## Extract metadata from tree
    tip_metadata <- strsplit(rerooted_tree_wo_ref$tip.label, "_")
    tip_sampleid <- strsplit(rerooted_tree_wo_ref$tip.label, "__")
    
    if (analysis_type == "refbased"){
      sampleid <- sapply(tip_sampleid, function(x) x[2])
      date <- sapply(tip_metadata, function(x) x[3])
      date <- as.Date(date, format = "%y%m%d")
      min_species_date <- min(date, na.rm = TRUE)
      max_species_date <- max(date, na.rm = TRUE)
      species_observation_range[[paste(patient, species, collapse = '_')]] <- c(min_species_date, max_species_date)
      days_since_obs = as.numeric(date - min_species_date)
      STs <- sapply(tip_metadata, function(str_list) {
        str_contains_ST <- grepl('ST', str_list)
        if (any(str_contains_ST)) {
          st = str_list[str_contains_ST]
          return(sub('\\?$', '', sub('^ST', '', st))) ## remove low confidence and ST strings
          } else {return(NA)}
        })
      if (grepl("Ehormaechei", patient_species)){
        STs <- rep(NA, length(tip_metadata)) ## has not STs@
        }
      ST_variant <- grepl('x$', STs, perl = TRUE) & !grepl('NFx$', STs, perl = TRUE) ## grep true STs (not NF) and check if they are an allele of the ST
      STs <- sub('x$', '', STs) ## remove the x at the end which is saved in ST_variant
    }else if(analysis_type == "denovo"){
      sampleid <- sapply(tip_sampleid, function(x) x[2])
      date <- sapply(tip_metadata, function(x) x[3])
      date <- as.Date(date, format = "%Y%m%d")
      refbased_pat_spec <- paste(patient, strsplit(species, '-')[[1]][1], collapse = '_')
      min_species_date <- as.Date(species_observation_range[[refbased_pat_spec]][1], format = "%Y-%m-%d")
      max_species_date <- as.Date(species_observation_range[[refbased_pat_spec]][2], format = "%Y-%m-%d")
      days_since_obs = as.numeric(date - min_species_date)
      STs <- rep(NA, length(tip_metadata))
      ST_variant <- rep(NA, length(tip_metadata))
    }else(sampleid <- NA)
    
    tree_annotation_heatmap_df <- tibble(
      tip_label = rerooted_tree_wo_ref$tip.label,
      patient = sapply(tip_metadata, function(x) x[1]),
      timepoint = sapply(tip_metadata, function(x) x[2]),
      date = date,
      sampletype = sapply(tip_metadata, function(x) x[4]),
      location = sapply(tip_metadata, function(x) x[5]),
      sampleid = sampleid,
      days_since_obs = days_since_obs,
      ST = STs,
      ST_variant = ST_variant
    ) %>%
      mutate(sampleid = if_else(str_detect(tip_label, "^OG_"), 
                                gsub("_P.*_TP0$|^OG_|seq", "", tip_label), 
                                sampleid),
             ) %>%
      as.data.frame()
    
    ## set everything for outgroup to NA
    rows_to_na <- grepl("^OG_", tree_annotation_heatmap_df$tip_label)
    cols_to_na <- setdiff(names(tree_annotation_heatmap_df), c("sampleid", "tip_label"))
    tree_annotation_heatmap_df[rows_to_na, cols_to_na] <- NA
    
    rownames(tree_annotation_heatmap_df) <- tree_annotation_heatmap_df$sampleid
    
    ## replace tree tip labels with sampelid only 
    rerooted_tree_wo_ref$tip.label <- tree_annotation_heatmap_df$sampleid[match(rerooted_tree_wo_ref$tip.label, tree_annotation_heatmap_df$tip_label)] 
    
    ## get nodes to highlight 
    tips_to_highlight_branch <- c()
    if (length(tips_to_highlight_branch) > 0){
      highlight_clade <- MRCA(rerooted_tree_wo_ref, tips_to_highlight_branch)
    }
    
    ## get nodes to annotate
    if (grepl("Ehormaechei", patient_species)){
      gt_to_isolate_file <- "/Users/ad_lemur/Documents/Data/Drylab/mf_2020_hap/muller_plots/P07_Ehormaechei-c1/P07_Ehormaechei_genotype_to_isolates.csv"
      gt_to_isolates <- read.csv(gt_to_isolate_file, header=FALSE, row.names = NULL, sep = ',')
      
      muts_to_annotate <- c('rpiR[T138P]', 'iscR[G60D]', 'gyrA[S83Y]', 'gyrA[S83F]', 'LT12530[G275W]', 'fimZ[F126L]', 'citG[L68Q]', 'gyrA[D87G]', 
                            "dacA[T243I]", "gyrA[G81C]", "narJ[E78*]", "narL[E49*]", "romA|acrR")
      
      tips_to_mark_mut <- list()
      for (mut in muts_to_annotate){
        isolates_with_mut <- gt_to_isolates[gt_to_isolates$V1 == mut, 'V2']
        tips_to_mark_mut[[mut]] <- strsplit(isolates_with_mut, ";")[[1]]
      }
      
    }else{tips_to_mark_mut <- c()}
    
    ## generate for heatmaps a color scale
    location_colors <- c("Nasal" = "#72e0ba",
                         "Oral" = "#41c49e",
                         "Skin" = "#1b9e77",
                         "Rectal" = "#146252",
                         "Gastric" = "#0f3f33",
                         "Lung" = "#eb4824", 
                         "Blood" = "#eb4824")
    mandatory_sites_for_legend <- c("Nasal", "Oral", "Skin", "Rectal")
    
    ####################
    ## plot initial tree 
    ####################
    max_tree_depth <- max(node.depth.edgelength(rerooted_tree_wo_ref))
    number_tips <- length(rerooted_tree_wo_ref$tip.label)
    ## heatmap specific variables
    colnames_offset <- number_tips/250
    width = 0.1
    width_offset = max_tree_depth*(width*1.1)
    plot_offset = max_tree_depth*0.23
    
    print(species)
    print(max_tree_depth)
    
    ## plot tree
    tree_plt <- ggtree(rerooted_tree_wo_ref, size = 0.8) +
      geom_tiplab(align = TRUE, linetype = "dotted", size = 3, linesize = 0.2) +
      geom_rootedge(max_tree_depth/50, size = 0.8) + 
      geom_treescale(x = max_tree_depth*1.4, y = -2) + ## add tree scale
      xlim(-max_tree_depth/25, max_tree_depth*1.5) + ## expand on right of plot
      vexpand(1/number_tips*8, 1) + ## expand on top of plot
      vexpand(1/number_tips/2.5, -1) + ## expand on bottom of plot
      coord_cartesian(clip = "off") + ## do not clip stuff outside of plotting area
      ggtitle(paste(patient, species))
      
    ####################
    ## annotate tree 
    ####################
    tree_plt_data <- tree_plt$data
    
    ## get nodes to mark via a star
    if (length(tips_to_mark_mut) > 0){
      for (mut in names(tips_to_mark_mut)){
        ## subset list of isolates just to those which are present in tree (required to avoid errors); Note: If wrong tree loaded, then some mutations might not be annotated!!
        isolates_on_node = tips_to_mark_mut[[mut]][tips_to_mark_mut[[mut]] %in% rerooted_tree_wo_ref$tip.label]
        mark_mut_node <- MRCA(rerooted_tree_wo_ref, isolates_on_node)
        # Find midpoint of branch given the node
        derived <- tree_plt_data[tree_plt_data$node == mark_mut_node, ]
        ancestral <- tree_plt_data[tree_plt_data$node == derived$parent, ]
        mid_x <- mean(c(ancestral$x, derived$x))
        ## mark mutations
        tree_plt <- tree_plt +
          geom_star(x = mid_x, y = min(derived$y), size = 5, fill = "grey") +
          annotate('text', x = mid_x, y = min(derived$y)+2, label = mut)
      }
    }
    
    ## generate heatmap dataframes 
    level_order <- union(mandatory_sites_for_legend, unique(na.omit(tree_annotation_heatmap_df$location)))
    loc_df <- data.frame(Location = tree_annotation_heatmap_df$location) %>%
      `rownames<-`(tree_annotation_heatmap_df$sampleid)
    
    date_df <- data.frame(Days = tree_annotation_heatmap_df$days_since_obs) %>%
      `rownames<-`(tree_annotation_heatmap_df$sampleid)
    
    if (!all(is.na(tree_annotation_heatmap_df$ST))){
      st_df <- data.frame(ST = tree_annotation_heatmap_df$ST) %>%
        `rownames<-`(tree_annotation_heatmap_df$sampleid) %>%
        mutate(ST = ifelse(is.na(ST), 'NA', ST),
               ST = factor(ST, levels = names(color_map)))
      st_variant_df <- tree_annotation_heatmap_df[(tree_annotation_heatmap_df$ST_variant == TRUE) & !is.na(tree_annotation_heatmap_df$ST_variant) , c('sampleid', 'ST_variant')]
      st_variant_df <- left_join(st_variant_df, select(tree_plt$data, c('label', 'x', 'y')), by = c('sampleid' = 'label'))
    }else(st_df <- data.frame())
    
    ## add location heatmap 
    tree_plt <- gheatmap(tree_plt, loc_df, 
                         offset = plot_offset, width = width,
                         color = 'black',
                         colnames_angle=90, colnames_position = "top", 
                         colnames_offset_y = colnames_offset, hjust = 0) +
      scale_fill_manual(values = location_colors, breaks = level_order, na.value="white", name="Location", drop = FALSE) 
      
    ## add new scale and days since observation heatmap 
    tree_plt <- tree_plt + new_scale_fill()
    tree_plt <- gheatmap(tree_plt, date_df, 
                         offset=plot_offset+width_offset, width=width,
                         color = 'black',
                         colnames_angle=90, colnames_position = "top", 
                         colnames_offset_y = colnames_offset, hjust = 0) +
      scale_fill_gradient(low = "white", high = "black", na.value="white",
                          limits = c(0, max_species_date-min_species_date), 
                          name="Days since\nobservation")
    if (any(!is.na(st_df$ST))) {
      ## add new scale and ST heatmap 
      tree_plt <- tree_plt + new_scale_fill() +
        xlim(0, max_tree_depth + plot_offset+width_offset*2 + 0.25)
      tree_plt <- gheatmap(tree_plt, st_df, 
                           offset=plot_offset+width_offset*2, width=width,
                           color = 'black',
                           colnames_angle=90, colnames_position = "top", 
                           colnames_offset_y = colnames_offset, hjust = 0) +
        scale_fill_manual(values = color_map, name="Sequence type") +
        geom_point(data = st_variant_df, aes(y = y), x = max_tree_depth+plot_offset+width_offset*2+width*0.94, shape = 21, fill = "black", size = 2)
    }
    
    ## save plot
    pdf(paste0(outpath, "/", patient, '_', species_label, '_', analysis_type, '_' , tree_type, '.pdf'), 
        width = 8, height = number_tips/10+2.5, paper = "special", onefile = F)
    print(tree_plt)
    dev.off()
  }
}
