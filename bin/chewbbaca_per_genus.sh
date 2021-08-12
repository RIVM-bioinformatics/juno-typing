set -euo pipefail

# Input user
input_files=$(realpath "$1")
threads="$2"
output_dir="$3"
db_dir="$4"
genus="$5"

# Make new variables
downloaded_scheme="${db_dir}/downloaded_schemes/${genus}"
prepared_scheme="${db_dir}/prepared_schemes/${genus}/scheme"
script_path="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

echo "Deleting any previous results from old ChewBBACA runs if existing in ${output_dir}...\n"
rm -rf "results_*"  

if [ ! -d "${prepared_scheme}" ]; then
    echo "\nPreparing scheme for running it with ChewBBACA...\n"
    chewBBACA.py PrepExternalSchema -i "${downloaded_scheme}" \
        --output-directory "${prepared_scheme}" \
        --cpu ${threads}
    echo "Prepared scheme can be found at ${prepared_scheme}.\n"
fi

echo "Making output directory ${output_dir}...\n"
mkdir -p ${output_dir}
cd "${output_dir}"

echo "Running ChewBBACA for ${genus} scheme...\n"
chewBBACA.py AlleleCall --cpu ${threads} \
                -i "${input_files}" \
                -g "${prepared_scheme}" \
                -o "." \
                --fr

find "." -type f -name "results_alleles.tsv" -exec mv {} "." \;