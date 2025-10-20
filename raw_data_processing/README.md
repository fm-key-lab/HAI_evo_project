# Raw data processing

The first step of the three analysis parts start with the processing of raw sequencing reads before follwing the analysis via local executable scripts. 

## `Snakemake` pipelines 

- `kraken2`: Taxonomic classification, sequence typing and *de novo* assembly for every individual sample
- `species_specific/mapping`: Alignment of quality-filtered reads against a reference genome & identification of SNVs
- `species_specific/case`: Generation of a multi-dimensional matrix for each reference genome containing candidate SNVs and respective summary metrics for each variant, and two coverage matrices
- `lineage_specific/assembly`: Assembly of lineage-specific *de novo* genomes
- `lineage_specific/mapping`: Alignment of quality-filtered reads against the respective *de novo* assembled lineage-specific genome & identification of SNVs and indels
- `lineage_specific/case`: Generation of a multi-dimensional matrix for each *de novo* assembled lineage-specific genome containing candidate SNVs, indels and respective summary metrics for each variant and two coverage matrices

## Installation of `Snakemake`

To install `Snakemake` you can use the `Snakemake` conda environment (`../envs/snakemake_conda_env.yaml`). Please note, that `Snakemake` requires a [`mamba`](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) installation. 

## `Snakemake` pipeline organization

Each `Snakemake` pipeline follows the same structure:

- `Snakefile`: `Snakemake` workflow with individual steps (called rules)
- `snakemakeslurm.sh`: trigger file to start the `Snakemake` pipeline
- `samples.csv`: comma-separated input file for the pipeline with one sample per row
    - Depending on the pipeline, the file contains the following columns: 
        - `Path`: absolute path of the directory where the raw sequencing reads of the sample is located
        - `Sample`: sample ID 
        - `ReferenceGenome`: directory name containing the reference genome (must contain `genome.fasta`)
        - `ProviderName`: unique raw data file name within `Path` without `_R[1-2].fastq.gz` suffix
        - `Subject`: subject/patient ID
        - `Outgroup`: boolean for outgroup (0 = ingroup sample; 1 = outgroup sample)
- `cluster.slurm.json`: cluster parameters and required resources by individual `Snakemake` rules
- `scripts/`: scripts and functions used within the pipeline
- `envs/`: conda environments required for individual rules

## Mandatory modifications

Before running the `Snakemake` pipeline, minor modifications need to be conducted. 

- `Snakefile`: Update the paths pointing to the reference genomes & databases 
- `snakemakeslurm.sh`: Update the `tmp/` directories and the `conda-prefix` path 
- `samples.csv`: Update paths, sample names etc for individual usage. 
    - The samples csvs used for the various analyses are stored within the respective directories
    - Note: The samples csvs for the `case` pipeline will be generated at the end of the `mapping` pipeline 
    - Note: If 2 fastqs should be combined for subsequent analysis, one need to provide the respective paths and provider names whitespace separated on one row. 
- `cluster.slurm.json`: update the email address and change the parameters to make the resources fit to your computing cluster

## Run `Snakemake`

1. Generate a log directory within the pipeline's directory `mkdir log`. 
2. Run the snakemake pipeline by calling `bash snakemakeslurm.sh 2>&1 | tee logs/0_primary_${date}.log` from within the pipeline's directory. 



