#############################################################################
#####                               MLST7                               #####
#############################################################################

checkpoint enlist_samples_for_cgmlst_scheme:
    input:
        species = expand(OUT + '/identify_species/{sample}/best_species_hit.txt', sample=SAMPLES),
        sample_sheet = sample_sheet
    output:
        directory(OUT + '/cgmlst')
    log:
        OUT + '/log/cgmlst/list_samples_per_cgmlst_scheme.log'
    threads: config['threads']['other']
    resources: mem_gb=config['mem_gb']['other']
    shell: 
        """
mkdir -p {output}
python bin/chewbbaca_input_files.py --sample-sheet {input.sample_sheet} \
                                --output-dir {output} \
                                --best-hit-kmerfinder-files {input.species} &> {log}
if [ ! $? == 0 ]; then
    rm -rf {output}
    exit 1
fi
        """

#----------------------- Choose cgMLST scheme per genus ----------------------#

rule cgmlst_per_genus:
    input:
        input_files = OUT + '/cgmlst/{genus}_samples.txt'
    output:
        chewbbaca_result = OUT + '/cgmlst/{genus}/results_alleles.tsv'
    conda:
        '../../envs/chewbbaca.yaml'
    log:
        OUT + '/log/cgmlst/chewbbaca_{genus}.log'
    threads: config['threads']['chewbbaca']
    resources: mem_gb=config['mem_gb']['chewbbaca']
    params:
        output_dir = abspath(OUT + '/cgmlst/{genus}'),
        db_dir = CGMLST_DB,
        genus = '{genus}'
    shell:
        """
bash bin/chewbbaca_per_genus.sh {input} \
                                {threads} \
                                {params.output_dir} \
                                {params.db_dir} \
                                {params.genus} &> {log}
        """


rule hash_cgmlst:
    input:
        OUT + '/cgmlst/{genus}/results_alleles.tsv'
    output:
        OUT + '/cgmlst/{genus}/hashed_results_alleles.csv'
    log:
        OUT + '/log/cgmlst/hashed_chewbbaca_{genus}.log'
    threads: config['threads']['other']
    resources: mem_gb=config['mem_gb']['other']
    params:
        db_dir = CGMLST_DB + '/prepared_schemes/',
        genus = '{genus}'
    shell:
        """
python bin/get_allele_hashes.py --scheme-name {params.genus} \
                                --output {output} \
                                --db-dir {params.db_dir} \
                                --threads {threads} \
                                {input} &> {log}
        """

#----------------------- Finish cgMLST rule ----------------------#
# This rule is just for the checkpoint to work and choose the 
# correct cgMLST scheme. Does not produce real outupt.

def cgmlst_output(wildcards):
    checkpoint_output = checkpoints.enlist_samples_for_cgmlst_scheme.get(**wildcards).output[0]
    return expand(OUT + '/cgmlst/{genus}/hashed_results_alleles.csv',
            genus=glob_wildcards(os.path.join(checkpoint_output, "{genus}_samples.txt")).genus)


rule cgmlst_multireport:
    input:
        cgmlst_output
    output:
        temp(OUT + '/cgmlst_multireport.csv')
    threads: config['threads']['other']
    resources: mem_gb=config['mem_gb']['other']
    shell:
        """
cat {input} > {output}
        """