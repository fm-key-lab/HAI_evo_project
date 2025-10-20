#!/bin/bash

#################
##  SETUP RUN  ##
#################

## GET VARIABLES
## general 
WORKING_DIR=$1 ## Working directory where to save all outputs
## AF3 
AF3_JSON_IN=$2 ## full path of json file (NOTE: basename will be the project directory!)
DEFAULT_AF_PARAMETERS=${HOME}/alphafold_3_0_0/parameters.inc
PROJECT_NAME=$( basename ${AF3_JSON_IN} .json )
AF3_WORKING_DIR=${WORKING_DIR}/${PROJECT_NAME}
## APBS
## APBS_EXECUTABLE=${HOME}/library/apbs/APBS-3.4.1/bin/apbs
## DEFAULT_APBS_PARAMETERS=${HOME}/library/apbs/apbs_input_template.in
## AF3_APBS_WORKING_DIR=${AF3_MSA_OUTPUT_DIR}/apbs
## CIF_FILE=${AF3_MSA_OUTPUT_DIR}/${PROJECT_NAME_LOWER}/${PROJECT_NAME_LOWER}_model.cif

## LOAD MODULES
module load alphafold/3.0.0

## HOW TO RUN:
## bash run_AF3_MPCDF.sh /ptmp/mfenk/data/mf_2020_hap/batch1/prot_prediction/alphafold3 fimZ/fimZ.json

#################
##  START RUN  ##
#################

echo -e "\n\n\nSTARTING\n" 
echo "   ######   ##      >=====>=====>       ####  "
echo " ##     ## ##       ||    ||   ||      #   ## "
echo "##      ###         <=====<=====<         ##  "
echo " ##     ## ##       ||    ||   ||      #   ## "
echo "  ######    ##      >=====>=====>       ####  "
echo -e "\n\n\n"

echo "Project directory was set to ${AF3_WORKING_DIR} accordingly to the json file name"
mkdir -p ${AF3_WORKING_DIR}

## SANITY CHECK
## check if path of json file is full path
if [[ ${AF3_JSON_IN} != /* ]]; then
  AF3_JSON_IN=${PWD}/${AF3_JSON_IN}
fi
if [[ $( dirname ${AF3_JSON_IN} ) != ${AF3_WORKING_DIR} ]]; then
  cp ${AF3_JSON_IN} ${AF3_WORKING_DIR}/${PROJECT_NAME}.json
  AF3_JSON_IN=${AF3_WORKING_DIR}/${PROJECT_NAME}.json
fi

## change to working directory
cd ${AF3_WORKING_DIR}

##################
## RUN AF3 Server (https://alphafoldserver.com/welcome) or via a local installation using the AF3_input/fimZ.json file

##################
## START ABPS

## ELECTROSTATIC ESTIMATION
sbatch ~/alphafold_3_0_0/jobscript-alphafold3-APBS-step_3-electrostatic.sh


##################
## INFO:
##################

## APBS INSTALLATION
##   conda create -n pqr python=3.12
##   conda activate pqr
##   cd ~/library/apbs
##  ## install propka (dependency of pdb2pqr)
##   wget https://github.com/jensengroup/propka/archive/refs/tags/v3.5.1.tar.gz
##   untar v3.5.1.tar.gz
##   cd propka-3.5.1/ ; pip install . ; cd ../
##   pip install pdb2pqr ## requires python >= 3.8
##   conda install -c conda-forge gemmi
## 
##  ## install apbs
##   wget https://github.com/Electrostatics/apbs/releases/download/v3.4.1/APBS-3.4.1.Linux.zip
##   unzip APBS-3.4.1.Linux.zip


## APBS INPUT FILE
## Example APBS in file (= ${DEFAULT_APBS_PARAMETERS}); taken from APBS website -- calculates electrostatic potential with defaults:

##  read
##      mol pqr PQRINFILE
##  end
##  elec
##      mg-auto
##      dime 161 129 161
##      cglen 95.3445 70.9172 101.9609
##      fglen 76.085 61.716 79.977
##      cgcent mol 1
##      fgcent mol 1
##      mol 1
##      lpbe
##      bcfl sdh
##      pdie 2.0
##      sdie 78.54
##      srfm smol
##      chgm spl2
##      sdens 10.0
##      srad 1.4
##      swin 0.3
##      temp 298.15
##      calcenergy total
##      calcforce no
##      write pot dx DXOUTFILE
##  end
##  quit

