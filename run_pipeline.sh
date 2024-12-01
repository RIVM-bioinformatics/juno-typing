#!/bin/bash

# Wrapper for juno typing pipeline

set -euo pipefail

#----------------------------------------------#
# User parameters
input_dir="${1%/}"
output_dir="${2%/}"

#----------------------------------------------#
## make sure conda works

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/mnt/miniconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/mnt/miniconda/etc/profile.d/conda.sh" ]; then
        . "/mnt/miniconda/etc/profile.d/conda.sh"
    else
        export PATH="/mnt/miniconda/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<export -f conda
export -f __conda_activate
export -f __conda_reactivate
export -f __conda_hashr


#----------------------------------------------#
# Create the environment

# we can use the base installation of mamba to create the environment. 
# Swapping to a parent env is not necessary anymore.
mamba env create -f envs/juno_typing.yaml --name pipeline_env
conda activate pipeline_env

#----------------------------------------------#
# Run the pipeline

echo -e "\nRun pipeline..."

if [ ! -z ${irods_runsheet_sys__runsheet__lsf_queue} ]; then
    QUEUE="${irods_runsheet_sys__runsheet__lsf_queue}"
else
    QUEUE="bio"
fi

set -euo pipefail

python juno_typing.py --queue "${QUEUE}" -i "${input_dir}" -o "${output_dir}"

result=$?

python cr_meta.py "${output_dir}"

# Propagate metadata

set +euo pipefail 

SEQ_KEYS=
SEQ_ENV=`env | grep irods_input_sequencing`
for SEQ_AVU in ${SEQ_ENV}
do
    SEQ_KEYS="${SEQ_KEYS} ${SEQ_AVU%%=*}"
done

for key in $SEQ_KEYS irods_input_illumina__Flowcell irods_input_illumina__Instrument \
    irods_input_illumina__Date irods_input_illumina__Run_number irods_input_illumina__Run_Id
do
    if [ ! -z ${!key} ] ; then
        attrname=${key:12}
        attrname=${attrname/__/::}
        echo "${attrname}: '${!key}'" >> ${output_dir}/metadata.yml
    fi
done

set -euo pipefail 

exit ${result}
# Produce svg with rules
# snakemake --config sample_sheet=config/sample_sheet.yaml \
#             --configfiles config/pipeline_parameters.yaml config/user_parameters.yaml \
#             -j 1 --use-conda \
#             --rulegraph | dot -Tsvg > files/DAG.svg
