set -euo pipefail

export irods_runsheet_sys__runsheet__lsf_queue='bio'
export irods_input_illumina__Flowcell='XXXX'
export irods_input_illumina__Instrument='nextseq'
export irods_input_illumina__Date='20220214' 
# export irods_input_illumina__Run_number='0001' 
# export irods_input_illumina__Run_Id='XXXX'
export irods_input_sequencing__fake_arg='aaaa'

input_dir="/data/BioGrid/hernanda/test_data_per_pipeline/Enteric/Juno-typing/"
output_dir="test_output"

mkdir -p "${output_dir}"

set +euo pipefail 

./run_pipeline.sh "${input_dir}" "${output_dir}"
result="$?"

if [ "$result" != 0 ]
then
    echo "Something went wrong with the pipeline! (Exit ${result})"
    # rm -rf "${output_dir}"
    exit 1
fi

if [ ! -f "${output_dir}/metadata.yml" ]
then
    echo "Metadata not produced!"
    # rm -rf "${output_dir}"
    exit 1
fi

echo "Run sucessful!"
rm -rf "${output_dir}"