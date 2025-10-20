#!/bin/bash
timestamp=$(date +%Y_%m_%d)

WORKING_DIR="${HOME}/data_analysis/2023/public_data_analysis"
SCRIPTS_DIR="${WORKING_DIR}/public_data/scripts"
RAXML_WORKING_DIR="${WORKING_DIR}/public_data/analysis/case_Pub_data_P07_Ehormaechei-c1_240825/${timestamp}_P07_Ehormaechei-c1_240825_gubbinsdefault_ambig10/raxml_bootstrapped_tree1000/"
RAXML_ALIGNMENT_FILE="${WORKING_DIR}/public_data/analysis/case_Pub_data_P07_Ehormaechei-c1_240825/${timestamp}_P07_Ehormaechei-c1_240825_gubbinsdefault_ambig10/snp_msa.fa"
RAXML_TREE_PREFIX="Ehormaechei_pubdata"
RAXML_OUTGROUP="Equasihormaechei_GCF_004331385_1"


## activate conda environment 
activate_conda_path=$(echo $CONDA_EXE | sed 's#bin/conda#bin/activate#g')
source ${activate_conda_path} anapy_312

cd ${SCRIPTS_DIR}

echo "Starting analysis py"
date

## run python unbuffered to get immediate stdout prints
python -u hap2020_analysispy_public_data_logan_v1_2_gubbins.py

date

echo "Starting RaxML"

source ${activate_conda_path} raxml

mkdir -p ${RAXML_WORKING_DIR}
cd ${RAXML_WORKING_DIR}

echo "full tree, 1000 bootstraps"
raxml-ng --all --threads 70 --seed 123 --msa ${RAXML_ALIGNMENT_FILE} --outgroup ${RAXML_OUTGROUP} --model GTR+G --prefix ${RAXML_TREE_PREFIX}_1000bs_ML_all_s123 --bs-trees 1000

raxml-ng --support --seed 123 --tree ${RAXML_TREE_PREFIX}_1000bs_ML_all_s123.raxml.bestTree --bs-trees ${RAXML_TREE_PREFIX}_1000bs_ML_all_s123.raxml.bootstraps --prefix ${RAXML_TREE_PREFIX}_1000bs_ML_s123_with_bootstraps
