# Online Code/Metadata to analzye the evolution of pathogens causing hospital-associated infections within critically ill patients

Here we present the computational code part of the analysis presented in: 

Fenk, M., Hrdina, A., Winans, J. B., Soerensen, M., Ostertag, L., Coquery, E., Sow, F., Petros, S., Stingu, C.-S., Lippmann, N., Nadell, C. D., Iatsenko, I., Pasieka, B., & Key, F. M. (2025). Colonization, translocation, and evolution of opportunistic pathogens during hospital-associated infections. bioRxiv. https://doi.org/10.1101/2025.10.23.683921

We hope you find the code useful. In case you recycle it for your own analyses please cite our study.

## Overview

The analysis is grouped into three parts

- [Species-specific analyses](#species-specific-analyses)
    - [Raw data processing](#species-specific-raw-data-processing)
    - [SNV-based analysis](#species-specific-snv-based-analysis)
- [Lineage-specific analyses](#lineage-specific-analyses)
    - [Raw data processing](#lineage-specific-raw-data-processing)
    - [SNV- and indel-based analysis](#lineage-specific-snv--and-indel-based-analysis)
    - [MGE-based analysis](#lineage-specific-mge-based-analysis)
    - [Phenotype analysis](#phenotype-analysis)
- [Public data analyses](#public-data-analyses)
    - [Raw data processing](#public-data-raw-data-processing)
    - [SNV-based analysis](#public-data-snv-based-analysis)

<br/><br/>

## Prerequisites

The analyses are built upon [`Snakemake`](https://snakemake.readthedocs.io/en/stable/) pipelines which are used to process raw sequencing reads on a computing cluster with a `slurm` workload manager. `Snakemake` can be installed using [`miniconda`](https://www.anaconda.com/docs/getting-started/miniconda/main) given the provided conda environment (`envs/snakemake_conda_env.yaml`). Subsequent analyses are mainly written using `Python` (version 3.12; see provided conda environment `envs/analysis_conda_env.yaml`) with some parts done in `R` (version 4.4.0; see required library versions `envs/R_library_versions.txt`). Prior to execution of any script, make sure the correct paths are provided and the input is in the correct format.

The raw sequencing reads generated within this study are available on the ENA BioProject [PRJEB98220](https://www.ebi.ac.uk/ena/browser/view/PRJEB98220) and the Biosample IDs of utilized public data are available within `metadata/public_data/Pub_data_biosampleid_list.txt`. All sequencing data should be downloaded using e.g. `nf-core/fetchngs` or `fasterq-dump` prior to analysis.

<br/><br/>

## Species-specific analyses

The first part of the analysis processes raw sequencing reads on a computing cluster using `Snakemake` pipelines with published bioinformatic tools, classifies the taxonomy of every sample, aligns the quality-controlled reads against a reference genome, calls single-nucleotide variants (SNVs) and proccesses those with related summary metrics into a multi-dimensional matrix. The second part can be performed on a local machine. It filters samples and SNVs based on which the fine-grained evolutionary analysis can be performed.

### Species-specific raw data processing

Prerequisites: 
- Download of sequencing reads [PRJEB98220](https://www.ebi.ac.uk/ena/browser/view/PRJEB98220)
- Installation of the `Snakemake` conda environment (`envs/snakemake_conda_env.yaml`)
- Access to a computing cluster with `slurm` workload manager

Raw data processing is implemented to be conducted on a computing cluster (see [Prerequisites](#prerequisites)) using multiple `Snakemake` pipelines. 

1. `raw_data_processing/kraken2`: Taxonomic classification to estimate the abundance of reads of the focal species within each sample, sequence typing & *de novo* genome assemblies of each individual sample.
    - prior to run the pipeline, a database for kraken and bracken need to be built
    - to run the samples on the data, use `raw_data_processing/kraken2/HAI_samples.csv`
2. `raw_data_processing/species_specific/mapping`: Alignment of quality-filtered reads against a reference genome
3. `raw_data_processing/species_specific/case`: Generation of a multi-dimensional matrix (`candidate_mutation_table.pickle.gz`) for each reference genome containing candidate SNVs and respective summary metrics

Further information can be found [here](raw_data_processing/). 

### Species-specific SNV-based analysis

Prerequisites: 
- Species-specific raw data processing (required output file: `candidate_mutation_table.pickle.gz`, `cov_raw_sparsecsr_mat.npz`)
- Installation of the analysis conda environment (`envs/analysis_conda_env.yaml`)

Following the raw data processing, the species-specific evolutionary analysis (`local_analysis/species_specific_analysis/SNV_analysis/hap2020_analysispy_refbased_v3.py`) utilizes the output from the raw data processing (`candidate_mutation_table.pickle.gz`) to generate pathogen-patient pair specific maximum likelihood phylogenies to identify lineage clusters.

Further information can be found [here](local_analysis/species_specific_analysis/). 

<br/><br/>

## Lineage-specific analyses 

After identification of lineage clusters, quality-filtered sequencing reads are utilized on a computing cluster using `Snakemake` pipelines to assemble *de novo* lineage-specific genomes against which the reads are aligned and obtained single-nucleotide variants (SNVs) and indels are processed with related summary metrics into a multi-dimensional matrix, while coverage matrices are generated for mobile genetic element (MGE) identification. The subsequent analyses are designed to run on a local machine and filter samples and candidate SNVs and indels and identifies MGEs.

### Lineage-specific raw data processing

Prerequisites: 
- Download of sequencing reads [PRJEB98220](https://www.ebi.ac.uk/ena/browser/view/PRJEB98220)
- Installation of the `Snakemake` conda environment (`envs/snakemake_conda_env.yaml`)
- Access to a computing cluster with `slurm` workload manager
- Samples clustered to distinct lineages based on the identified lineage cutoff and the phylogenies from the species-specific analysis (`raw_data_processing/lineage_specific/assembly/sample_lineage_files`)

Raw data processing is implemented to be conducted on a computing cluster (see [Prerequisites](#prerequisites)) using multiple `Snakemake` pipelines. 

0. `raw_data_processing/kraken2`: (required from [Species-specific analysis](#species-specific-analyses))
1. `raw_data_processing/lineage_specific/assembly`: Assembly of lineage-specific *de novo* genomes
    - based on the reconstructed phylogenies and the taxonomic classification from the species-specific analysis, high-quality samples of each lineage need to be provided (see `raw_data_processing/lineage_specific/assembly/sample_lineage_files`)
2. `raw_data_processing/lineage_specific/mapping`: Alignment of quality-filtered reads against the respective *de novo* assembled lineage-specific genome
3. `raw_data_processing/lineage_specific/case`: Generation of a multi-dimensional matrix for each *de novo* assembled lineage-specific genome containing candidate SNVs, indels and respective summary metrics for each variant and two coverage matrices

Further information can be found [here](raw_data_processing/). 

### Lineage-specific SNV- and indel-based analysis

Prerequisites: 
- Lineage-specific raw data processing (required output files: `candidate_mutation_table.pickle.gz`, `cov_raw_sparsecsr_mat.npz`)
- Installation of the analysis conda environment (`envs/analysis_conda_env.yaml`)

The lineage-specific SNV and indel analysis (`local_analysis/lineage_specific_analysis/SNV_indel_analysis/hap2020_analysispy_denovo_v2.py`) utilizes the output from the raw data processing (`candidate_mutation_table.pickle.gz`) to perform evolutionary analyses (phylogenetic reconstruction, ancestral allele inference, and molecular clock and parallel evolution analysis).

Further information can be found [here](local_analysis/lineage_specific_analysis). 

### Lineage-specific MGE-based analysis

Prerequisites: 
- Lineage-specific raw data processing (required output files: `cov_raw_sparsecsr_mat.npz`, `cov_norm_sparsecsr_mat.npz`)
- Installation of the analysis conda environment (`envs/analysis_conda_env.yaml`)

The MGE analysis (`local_analysis/lineage_specific_analysis/MGE_analysis/MGE_analysis_denovo_v1.2.py`) utilizes two coverage matrices generated at the end of the raw data processing (`cov_raw_sparsecsr_mat.npz`, `cov_norm_sparsecsr_mat.npz`) to identify gained/lost mobile genetic elements within each lineage.

Further information can be found [here](local_analysis/lineage_specific_analysis). 

### Phenotype analysis

Individual analysis of measured phenotype data (`local_analysis/lineage_specific_analysis/subsequent_analysis/phenotyping`).

<br/><br/>

## Public data analyses

Re-analysis of public data to gain insights about the frequency of the identified infection-associated mutation (*fimZ*[F126L]).

### Public data raw data processing

Prerequisites: 
- Download of sequencing data `metadata/public_data/Pub_data_biosampleid_list.txt`
- Installation of the `Snakemake` conda environment (`envs/snakemake_conda_env.yaml`)
- Access to a computing cluster with `slurm` workload manager
- Annotated *de novo* assembled lineage-specific *E. hormaechei* genome (`metadata/lineage_specific_assembled_genomes/P07_Ehormaechei-c1`)

Raw data processing is implemented to be conducted on a computing cluster (see [Prerequisites](#prerequisites)) using multiple `Snakemake` pipelines. 

1. `raw_data_processing/kraken2`: Taxonomic classification to estimate the abundance of reads of the focal species within each sample, sequence typing & *de novo* genome assemblies of each individual sample.
    - to run the samples on the data, use `raw_data_processing/kraken2/Pub_data_samples.csv`
2. `raw_data_processing/species_specific/mapping`: Alignment of quality-filtered reads (≥ 80% reads are assigned as *E. hormaechei*) against the *de novo* assembled lineage-specific *E. hormaechei* genome (see [Lineage-specific analysis](#lineage-specific-analyses))
3. `raw_data_processing/species_specific/case`: Generation of a multi-dimensional matrix for each reference genome containing candidate SNVs and respective summary metrics

Further information can be found [here](raw_data_processing). 

### Public data SNV-based analysis

Prerequisites: 
- Public data raw data processing (required output files: `candidate_mutation_table.pickle.gz`)
- Installation of the analysis conda environment (`envs/analysis_conda_env.yaml`) & [`Gubbins`](https://github.com/nickjcroucher/gubbins) (`conda create -n "gubbins" -c bioconda -c conda-forge -c r gubbins`)
- ≥64 GB RAM

The public data SNV analysis (`local_analysis/public_data_analysis/SNV_analysis/trigger_anapy_w_BSMLtree.sh`) utilizes the output from the raw data processing (`candidate_mutation_table.pickle.gz`) to perform filter isolates and candidate SNVs and generates the final phylogenetic tree.

Further information can be found [here](local_analysis/public_data_analysis). 

<br/><br/>

## Additional figures

Code for additional figures generated for the publication can be found in `local_analysis/additional_figure_generation`.

Further information can be found [here](local_analysis/additional_figure_generation). 
