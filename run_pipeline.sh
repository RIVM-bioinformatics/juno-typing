#!/bin/bash

# Wrapper for juno typing pipeline

#----------------------------------------------#
# Create/update necessary environments
PATH_MAMBA_YAML="envs/mamba.yaml"
PATH_MASTER_YAML="envs/master_env.yaml"
MAMBA_NAME=$(head -n 1 ${PATH_MAMBA_YAML} | cut -f2 -d ' ')
MASTER_NAME=$(head -n 1 ${PATH_MASTER_YAML} | cut -f2 -d ' ')

envs_list=$(conda env list)

if ! $(echo $envs_list | grep -q mamba)
then
    conda env update -f "${PATH_MAMBA_YAML}"
fi

source activate "${MAMBA_NAME}"

if ! $(echo $envs_list | grep -q "${MASTER_NAME}")
then
    mamba env update -f "${PATH_MASTER_YAML}"
fi

source activate "${MASTER_NAME}"

#----------------------------------------------#
# Run the pipeline

if [ ! -z ${irods_runsheet_sys__runsheet__lsf_queue} ]; then
    QUEUE="${irods_runsheet_sys__runsheet__lsf_queue}"
else
    QUEUE="bio"
fi

python juno_typing.py --queue "${QUEUE}" ${@}

# Produce svg with rules
# snakemake --config sample_sheet=config/sample_sheet.yaml \
#             --configfiles config/pipeline_parameters.yaml config/user_parameters.yaml \
#             -j 1 --use-conda \
#             --rulegraph | dot -Tsvg > files/DAG.svg