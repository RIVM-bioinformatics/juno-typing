#############################################################################
#####                               MLST7                               #####
#############################################################################

rule download_cgmlst_scheme:
    output:
        scheme = CGMLST_DB + '/downloaded_schemes/downloaded_schemes.yaml'
    log:
        OUT + '/log/cgmlst/download_schemes.log'
    threads: config['threads']['chewbbaca_preparation']
    resources: mem_gb=config['mem_gb']['chewbbaca_preparation']
    shell: 
        """
set -euo pipefail
python bin/download_cgmlst_scheme.py -t {threads} \
    -o $(dirname {output}) > {log}
touch {output}
        """


rule prepare_cgmlst_scheme:
    input:
        species = CGMLST_DB + '/downloaded_schemes/downloaded_schemes.yaml'
    output:
        scheme = CGMLST_DB + '/prepared_schemes/prepared_schemes.yaml'
    conda:
        '../../envs/chewbbaca.yaml'
    log:
        OUT + '/log/cgmlst/prepare_cgmlst_schemes.log'
    threads: config['threads']['chewbbaca_preparation']
    resources: mem_gb=config['mem_gb']['chewbbaca_preparation']
    shell: 
        """
set -euo pipefail
schemes_dir=$(dirname "{input}")
output_dir=$(dirname "{output}")
for folder in $(find "${{schemes_dir}}" -mindepth 1 -maxdepth 1 -type d); do

    species_name=$(basename "${{folder}}")
    output_species="${{output_dir}}/${{species_name}}"
    mkdir -p "${{output_species}}"

    chewBBACA.py PrepExternalSchema -i "${{folder}}" \
        -o "${{output_species}}" --cpu {threads} &> {log}

done

touch {output}
        """

#---------------- Choose cgMLST scheme according to genus/species ------------#

# The cgMLST for all samples needs to be run in one single step instead of in 
# one step per sample. This is because Chewbbaca creates a tmp folder in the 
# database folder and it does not accept that there are more than one running
# simultaneously. 
# TODO: If two genera are present in the same run, it is likely that they can 
# run simultaneously. I didn't take the time to look into it because in our use
# case this is hardly ever the case (at least not intentionally)
rule cgmlst:
    input:
        sample_sheet = sample_sheet,
        besthits = expand(OUT + '/identify_species/{sample}/best_species_hit.txt', sample=SAMPLES)
    output:
        result = OUT + '/cgmlst/cgmlst_finished.txt'
    conda:
        '../../envs/chewbbaca.yaml'
    log:
        OUT + '/log/cgmlst/chewbbaca.log'
    threads: config['threads']['chewbbaca']
    resources: mem_gb=config['mem_gb']['chewbbaca']
    params:
        output_dir = OUT + '/cgmlst',
        cgmlst_db_dir = CGMLST_DB + '/prepared_schemes'
    shell:
        """
python bin/chewbbaca_per_genus.py --best-hit-kmerfinder-files {input.besthits} \
                                --sample-sheet {input.sample_sheet} \
                                --output-dir {params.output_dir} \
                                --cgmlst-db-dir {params.cgmlst_db_dir} \
                                --threads {threads} &> {log}

chewbbaca_success=$?

if [ "$chewbbaca_success" == 0 ];then
    touch {output.result}
fi
        """


# rule cgmlst:
#     input:
#         assembly = lambda wildcards: SAMPLES[wildcards.sample]['assembly'],
#         scheme = CGMLST_DB + '/prepared_schemes/prepared_schemes.yaml',
#         species = OUT + '/identify_species/{sample}/best_species_hit.txt'
#     output:
#         result = OUT + '/cgmlst/{sample}_chewbbaca.txt',
#         input_file = OUT + '/cgmlst/{sample}/input_file.txt'
#     conda:
#         '../../envs/chewbbaca.yaml'
#     log:
#         OUT + '/log/cgmlst/chewbbaca_{sample}.log'
#     threads: config['threads']['chewbbaca']
#     resources: mem_gb=config['mem_gb']['chewbbaca']
#     params:
#         output_dir = OUT + '/cgmlst'
#     shell:
#         """
# genus=$(cat {input.species} | cut -f1 -d' ')
# scheme_dir="$(dirname "{input.scheme}")/${{genus}}"
# echo "{input.assembly}" > "{output.input_file}"

# if [[ -d "${{scheme_dir}}" ]]; then
#     chewBBACA.py AlleleCall --cpu 8 \
#                 -i "{output.input_file}" \
#                 -g "${{scheme_dir}}" \
#                 -o "{params.output_dir}" \
#                 --fr &> {log}
#     chewbbaca_success=$?
# else
#     echo -e "There is no cgMLST scheme downloaded and prepared for the genus ${{genus}}. This analysis was skipped." > {log}
#     chewbbaca_success=0
# fi

# # Since the output directory produced by chewbbaca has a timestamp 
# # and it is therefore difficult to find out which one belongs to which
# # sample, then a 'fake' sample needs to be produced and in a later step
# # are the output directories per sample renamed
# if [ "$chewbbaca_success" == 0 ];then
#     touch {output.result}
# fi
#         """

# rule cgmlst_multireport:
#     input:
#         expand(OUT + '/cgmlst/{sample}_chewbbaca.txt', sample=SAMPLES)
#     output:
#         multireport = OUT + '/cgmlst/cgmlst_multireport.csv',
#         # genus = expand(OUT + '/cgmlst/{sample}/results_alleles.tsv')
#     log:
#         OUT + '/log/cgmlst/multireport.log'
#     threads: config['threads']['other']
#     resources: mem_gb=config['mem_gb']['other']
#     params:
#         cgmlst_dir = OUT + '/cgmlst'
#     shell:
#         """
# python bin/chewbbaca_multireport.py --cgmlst_dir {params.cgmlst_dir} \
#                                     --output_multireport_file {output.multireport}
#         """