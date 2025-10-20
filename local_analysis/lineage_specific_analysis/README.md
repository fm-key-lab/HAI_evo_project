# Lineage-specific analysis

The analysis is spit into 3 parts: the SNV and indel-based analysis, the MGE identification, and the subsequent analyses. All analyses are designed to run on a local machine interactively. 

## Prerequisites
Prerequisites: 
- Lineage-specific raw data processing (required output files: `candidate_mutation_table.pickle.gz`, `cov_raw_sparsecsr_mat.npz`, `cov_norm_sparsecsr_mat.npz`)
- Installation of the analysis conda environment (`envs/analysis_conda_env.yaml`)
    - Note: If the conda environment cause problems to install, a minimal environment without versions is given within `envs/analysis_conda_env_minimal.yaml`
- Access to the analysis-specific module `../modules/analysispy_module.py` & `../modules/mobileelements_module.py`

## Analysis scripts

1. `SNV_indel_analysis/hap2020_analysispy_denovo_v2.py`: Filters samples and variants based on relevant statistics, reports SNVs and indels observed within each sample, reconstructs the phylogeny (SNV-based), fits the molecular clock, and identifies signatures of parallel evolution. As input, it requires the output files from the lineage-specific raw data processing (`candidate_mutation_table.pickle.gz`, `cov_raw_sparsecsr_mat.npz`).
2. `SNV_indel_analysis/P07_Ehormaechei_ampD_MSA_plt.py`: Uses individual sample-specific *de novo* assemblies (generated within `../../raw_data_processing/kraken2`) to generate multiple-sequence alignments of the annotated *ampD* genes to validate and identify additional indels within *E. hormaechei* samples.  
3. `MGE_analysis/MGE_analysis_denovo_v1.2.py`: Identifies gain/loss of mobile genetic elements using the raw coverage matrix (`cov_raw_sparsecsr_mat.npz`) containing absolute coverage and double-normalized coverage matrix (`cov_norm_sparsecsr_mat.npz`) across all isolates within each lineage. 
4. `subsequent_analysis/calculate_ratio_of_transitions_transversions_v4.py`: Calculates the deviation of each lineage's mutational spectrum from the observed background mutational spectrum observed across the used outgroups. 
5. `subsequent_analysis/tMRCA_v3.2_8iso.py`: Inferes and plots the time to the most common recent ancestor (TMRCA) using observed and literature based molecular clocks. 
6. `subsequent_analysis/hap_muller_softmax_glv_smoothing_v5_1_windels.py`: Uses reconstructed phylogenies to identify variants of interest (≥ 30% frequency shift), polarize and smooth the frequencies over the time course.
7. `subsequent_analysis/plot_ggmuller_smoothed_v2_Ehorm.r`: Smoothed longitudinal variant frequencies are used to plot the Muller plots
