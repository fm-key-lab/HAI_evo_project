# Species-specific analysis

The analysis is designed to run on a local machine interactively. It performs additional quality filtering on the samples and identified variants to reconstruct a phylogeny. This phylogeny allow the identification of lineage clusters required to *de novo* assemble lineage-specific genomes. 

## Prerequisites
Prerequisites: 
- Species-specific raw data processing (required output file: `candidate_mutation_table.pickle.gz`, `cov_raw_sparsecsr_mat.npz`)
- Installation of the analysis conda environment (`envs/analysis_conda_env.yaml`)
    - Note: If the conda environment cause problems to install, a minimal environment without versions is given within `envs/analysis_conda_env_minimal.yaml`
- Access to the analysis-specific module `../modules/analysispy_module.py`

## Analysis scripts

1. `SNV_analysis/hap2020_analysispy_refbased_v3.py`: Filters samples and variants based on relevant statistics and reconstructs the phylogeny. As input, it requires both output files from the species-specific raw data processing (`candidate_mutation_table.pickle.gz`, `cov_raw_sparsecsr_mat.npz`).
2. `subsequent_analysis/pairwise_SNV_diff_refbased_all_specs_v2.1.py`: Given the reconstructed tree, the distance to the most recent common ancestor is used to calculated the pairwise SNV difference within each pathogen-patient pair. 


