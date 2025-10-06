
library(ggplot2)
library(ggmuller)
library(cowplot)

build_flexible_grid <- function(plot_list, global_title = "") {
  n <- length(plot_list)
  if (is.null(names(plot_list))) {
    names(plot_list) <- seq_along(plot_list)
  }
  
  # add individual titles and remove legends after first
  plot_list <- lapply(seq_along(plot_list), function(i) {
    plot_list[[i]] +
      ggtitle(names(plot_list)[i]) +
      if (i > 1) theme(legend.position = "none") else theme()
  })
  
  ## build grid based on number of plots
  # helper to plot 2-column rows 
  row_2col <- function(a, b = NULL) plot_grid(a, b, ncol = 2)
  layout <- switch(
    as.character(n),
    "1" = plot_list[[1]],
    "2" = plot_grid(plotlist = plot_list, ncol = 1, rel_heights = c(3, 2)),
    "3" = plot_grid(plot_list[[1]], row_2col(plot_list[[2]], plot_list[[3]]), ncol = 1, rel_heights = c(3, 2)),
    "4" = plot_grid(plot_list[[1]],
                    row_2col(plot_list[[2]], plot_list[[3]]),
                    row_2col(plot_list[[4]], NULL),
                    ncol = 1,
                    rel_heights = c(3, 2, 2)),
    plot_grid(plot_list[[1]],
              row_2col(plot_list[[2]], plot_list[[3]]),
              row_2col(plot_list[[4]], plot_list[[5]]),
              ncol = 1,
              rel_heights = c(3, 2, 2))
  )
  
  layout <- plot_grid(
    ggdraw() + draw_label(global_title, fontface = 'bold', size = 16, x = 0, hjust = 0),
    layout,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )
  
  return(layout)
}

background_color <- c("background" = 'white')

species_path <- "/Users/ad_lemur/Documents/Data/Drylab/mf_2020_hap/muller_plots/P07_Ehormaechei-c1"

species <- gsub('.*/', '', species_path)
pop_files <- list.files(path = species_path, pattern = ".*_pop_.*0.2_windels.tsv", full.names = TRUE)
if (length(pop_files) == 0) {next}
plot_list <- list()
for (pop_file in pop_files){
  location = gsub(".*pop_", "", gsub("_freq.*", "", pop_file))
  edge_file <- gsub("_pop_", "_edges_", pop_file)
  iso_cnts_file <- gsub("_pop_", "_cnt_", pop_file)
  
  population <- read.table(pop_file, header=TRUE, sep = '\t')
  edges <- read.table(edge_file, header=TRUE, sep = '\t')
  iso_cnts <- read.table(iso_cnts_file, header=TRUE, sep = '\t')
  
  Muller_df <- get_Muller_df(edges, population)
  num_muts <- length(unique(Muller_df$Identity))
  ## get min/max x values to share them across plots
  if (location == 'collapsed'){
    xmin <- min(population$Generation)
    xmax <- max(population$Generation)
  }
  
  ## SNVs with Indels
  final_colors <- c(
                  "basal" = "#F0F0F0",
                  "rpiR[T138P]" = "#c8d5e1", 
                  "ampD[68AG]; iscR[G60D]" = "#a1bad3",
                  "gyrA[S83Y]" = "#7a9fc5",
                  "gyrA[S83F]" = "#5384b6",
                  "gyrA[S83F]; acrR[94TGGGGCCA]; LT12530[G275W]" = "#2c69a8",
                  "gyrA[S83F]; acrR[94TGGGGCCA]; LT12530[G275W]; eptA[273TC]" = "#054f9a", 
                  "fimZ[F126L]" = "#F9AE92",
                  "fimZ[F126L]; citG[L68Q]" = "#F26B4E",
                  "fimZ[F126L]; gyrA[D87G]" = "#DF3027"
  )
  ## remove indel sequence (squishes plot and does not add anything)
  new_labels <- c(
                  "basal" = "basal",
                  "rpiR[T138P]" = "rpiR[T138P]",
                  "ampD[68AG]; iscR[G60D]" = "ampD[68]; iscR[G60D]",
                  "gyrA[S83Y]" = "gyrA[S83Y]",
                  "gyrA[S83F]" = "gyrA[S83F]",
                  "gyrA[S83F]; acrR[94TGGGGCCA]; LT12530[G275W]" = "gyrA[S83F]; acrR[94]; LT12530[G275W]",
                  "gyrA[S83F]; acrR[94TGGGGCCA]; LT12530[G275W]; eptA[273TC]" = "gyrA[S83F]; acrR[94];\nLT12530[G275W]; eptA[273]",
                  "fimZ[F126L]" = "fimZ[F126L]",
                  "fimZ[F126L]; citG[L68Q]" = "fimZ[F126L]; citG[L68Q]",
                  "fimZ[F126L]; gyrA[D87G]" = "fimZ[F126L]; gyrA[D87G]"
  )
  
  inf_days <- c(0)
  inf_color <- c("#ED4724")
  ## sort colors in their given way
  print(unique(Muller_df$Identity))
  Muller_df$Identity <- factor(Muller_df$Identity, levels = names(final_colors))
  
  ## select timepoints to mark as sampled via dashed lines
  x_loc <- iso_cnts$Days
  x_annotation <- paste0(iso_cnts$Days, '\nn=', iso_cnts$cnts)
  
  exp_plt <- ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) +
    geom_area() +
    geom_vline(xintercept=inf_days, color = inf_color, linewidth = 0.5) + 
    geom_vline(xintercept = x_loc, color = "black", linetype = "dotted", linewidth = 0.2) +
    guides(linetype = FALSE, color = FALSE, fill = guide_legend(byrow = TRUE)) +
    scale_y_continuous(name = "Frequency\nof major genotypes", 
                       labels = 25 * (0:4), 
                       expand=c(0,0)) +
    scale_x_continuous(name = "Days since observation",
                       breaks = x_loc, 
                       labels = x_annotation, 
                       expand = c(0, 0),
                       limits = c(xmin, xmax)) +
    scale_fill_manual(name = "Genotype", values = final_colors, labels = new_labels) +
    scale_color_manual(values = final_colors, labels = new_labels) + 
    theme_minimal() +
    theme(legend.position = "right", panel.grid = element_blank()) +
    ggtitle(paste(species, location))
  
  plot_list[[location]] <- exp_plt
}

pdf(paste0(species_path, "/Muller_plot_combined_", species, "_windels.pdf"), 
    width = 10, height = 7.5, paper = "special", onefile = F)
print(build_flexible_grid(plot_list, global_title = species))
dev.off()