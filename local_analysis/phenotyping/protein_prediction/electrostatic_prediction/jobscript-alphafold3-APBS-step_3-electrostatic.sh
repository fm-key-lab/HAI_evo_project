#!/bin/bash -l
#SBATCH --job-name=AF3_APBS
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=UPDATE_YOUR_EMAIL_ADDRESS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --time=01:00:00


#################
##  SETUP RUN  ##
#################

## LOAD CONDA ENV
source /u/mfenk/library/miniconda3/bin/activate pqr

source parameters.inc

## GENERATE SUBDIR 
if [[ ${AF3_APBS_WORKING_DIR} = ${HOME}/ ]] || [[ ${HOME}/ == ${AF3_APBS_WORKING_DIR}* ]]; then
    echo "Given path is empty, parent of HOME or HOME. Exiting script"
    exit 1
else
    rm -rf ${AF3_APBS_WORKING_DIR}/
    mkdir -p ${AF3_APBS_WORKING_DIR}/
fi

##################
## RUN APBS ON EVERY CHAIN SEPARATELY

##################
## START CIF TO PQR TRANSLATION

## convert CIF file to PDB (required for pdb2pqr which is in turn required as input for ABPS)
gemmi convert ${CIF_FILE}  ${AF3_APBS_WORKING_DIR}/${PROJECT_NAME_LOWER}_model_tmp.pdb
sed '/^CRYST1/d' ${AF3_APBS_WORKING_DIR}/${PROJECT_NAME_LOWER}_model_tmp.pdb > ${AF3_APBS_WORKING_DIR}/${PROJECT_NAME_LOWER}_model.pdb ## remove crystallography header added by gemmi
rm ${AF3_APBS_WORKING_DIR}/${PROJECT_NAME_LOWER}_model_tmp.pdb

# Extract all unique chain IDs
chain_ids=$(awk '$1 == "ATOM" || $1 == "HETATM" {print substr($0, 22, 1)}' ${AF3_APBS_WORKING_DIR}/${PROJECT_NAME_LOWER}_model.pdb | sort -u)

# Split into one file per chain
for chain in $chain_ids; 
do
    awk -v c="$chain" '($1 == "ATOM" || $1 == "HETATM") && substr($0,22,1)==c {print} END {print "END"}' ${AF3_APBS_WORKING_DIR}/${PROJECT_NAME_LOWER}_model.pdb > ${AF3_APBS_WORKING_DIR}/${PROJECT_NAME_LOWER}_model_chain${chain}.pdb
    
    ## translate into pqr files
    ## --ff=PARSE specifies the force field (PARSE commonly used other options: CHARMM, AMBER).
    ## --with-ph=7.0 specifies the pH.
    ## --chain keeps the chain information in the output file.
    pdb2pqr \
        --ff=PARSE \
        --with-ph=7.0 \
        ${AF3_APBS_WORKING_DIR}/${PROJECT_NAME_LOWER}_model_chain${chain}.pdb \
        ${AF3_APBS_WORKING_DIR}/${PROJECT_NAME_LOWER}_model_chain${chain}.pqr


    ##################
    ## START APBS

    # Run psize.py once and store output to get the gridsize parameter required for calculation
    ## NOTE: change line 174 from `if nsmall <= 0:` to `if nsmall[i] <= 0:`
    PSIZE_OUTPUT=$( python ${HOME}/library/apbs/APBS-3.4.1/share/apbs/tools/manip/psize.py ${AF3_APBS_WORKING_DIR}/${PROJECT_NAME_LOWER}_model_chain${chain}.pqr ) 

    # Extract needed values from stored output
    read CGX CGY CGZ <<< $(echo "$PSIZE_OUTPUT" | awk '/Coarse grid dims/ {print $5, $7, $9}')
    read FGX FGY FGZ <<< $(echo "$PSIZE_OUTPUT" | awk '/Fine grid dims/ {print $5, $7, $9}')
    read DX DY DZ   <<< $(echo "$PSIZE_OUTPUT" | awk '/Num. fine grid pts/ {print $6, $8, $10}')

    # Replace placeholders in template and write to output
    ## modify apbs input file for calculation
    sed \
    -e "s/cglen.*/cglen $CGX $CGY $CGZ/" \
    -e "s/fglen.*/fglen $FGX $FGY $FGZ/" \
    -e "s/dime.*/dime $DX $DY $DZ/" \
    -e "s|PQRINFILE|${AF3_APBS_WORKING_DIR}/${PROJECT_NAME_LOWER}_model_chain${chain}.pqr|g" \
    -e "s|DXOUTFILE|${AF3_APBS_WORKING_DIR}/${PROJECT_NAME_LOWER}_model_chain${chain}_apbs|g" \
    ${DEFAULT_APBS_PARAMETERS} > ${AF3_APBS_WORKING_DIR}/apbs_${PROJECT_NAME_LOWER}_model_chain${chain}.in

    ## run apbs 
    ${APBS_EXECUTABLE} ${AF3_APBS_WORKING_DIR}/apbs_${PROJECT_NAME_LOWER}_model_chain${chain}.in

done

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
##  ## Modify psize.py script:
##  ## open file psize.py ( APBS-3.4.1/share/apbs/tools/manip/psize.py )
##  ## change line 174 from `if nsmall <= 0:` to `if nsmall[i] <= 0:`


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
