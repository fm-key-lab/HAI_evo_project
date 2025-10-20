AF3_MSA_JSON_PATH=/ptmp/mfenk/data/mf_2020_hap/batch1/prot_prediction/alphafold3/fimZ/fimZ.json                    ## full path of json file including .json extension
AF3_MSA_OUTPUT_DIR=/ptmp/mfenk/data/mf_2020_hap/batch1/prot_prediction/alphafold3/fimZ           ## absolute path for output of entire AF3 prediction
PROJECT_NAME_LOWER=$( basename ${AF3_MSA_JSON_PATH} .json | tr '[:upper:]' '[:lower:]' ) ## is in lower case!

APBS_EXECUTABLE=${HOME}/library/apbs/APBS-3.4.1/bin/apbs                                ## APBS executable
DEFAULT_APBS_PARAMETERS=${HOME}/library/apbs/apbs_input_template.in                     ## APBS input file
AF3_APBS_WORKING_DIR=${AF3_MSA_OUTPUT_DIR}/apbs                                         ## WORKING DIRECTORY (default in AF3 output dir)
CIF_FILE=${AF3_MSA_OUTPUT_DIR}/${PROJECT_NAME_LOWER}/${PROJECT_NAME_LOWER}_model.cif    ## Predicted protein by AF3
