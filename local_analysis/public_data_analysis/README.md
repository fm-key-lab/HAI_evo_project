# Public-data analysis

The analysis is designed to run on a local machine with ≥64 GB RAM. It performs quality filtering on the samples and identified variants, removes recombinant sites using `Gubbins`, reconstructs a maximum-likelihood phylogeny and identifies the dN/dS values for each gene. 

## Prerequisites
Prerequisites: 
- Public data raw data processing (required output files: `candidate_mutation_table.pickle.gz`)
- Installation of the analysis conda environment (`envs/analysis_conda_env.yaml`)
    - Note: If the conda environment cause problems to install, a minimal environment without versions is given within `envs/analysis_conda_env_minimal.yaml`
- Access to the analysis-specific module `../modules/analysispy_module.py`
- Install [`Gubbins`](https://github.com/nickjcroucher/gubbins) using `conda create -n "gubbins" -c bioconda -c conda-forge -c r gubbins`

## Analysis scripts

1. `SNV_analysis/trigger_anapy_w_BSMLtree.sh`: This file triggers the quality filtering of samples and SNVs (`SNV_analysis/hap2020_analysispy_public_data_logan_v1_2_gubbins.py`) and reconstructs a Felsenstein-bootstrapped maximum-likelihood phylogeny. As input, it requires the multi-dimensional matrix from the raw data processing (`candidate_mutation_table.pickle.gz`).
    - Note: `Gubbins` can be run on multiple cores and depending on the availability of threads need to be adapted within the `gubbins_params` of `SNV_analysis/hap2020_analysispy_public_data_logan_v1_2_gubbins.py`
2. `subsequent_analysis/plot_trees_pubdata.R`: Used to visualize reconstructed phylogeny (Figure S20a)
3. `subsequent_analysis/empirical_pval_dNdS_fimZ_v3.py`: Used to evaluate the empirical percentile of the dN/dS scores of genes within the *fim* operon. 