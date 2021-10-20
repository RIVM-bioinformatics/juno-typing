#############################################################################
#####                               MLST7                               #####
#############################################################################

checkpoint enlist_samples_for_cgmlst_scheme:
    input:
        sample_sheet
    output:
        directory(OUT + '/cgmlst')
    log:
        OUT + '/log/cgmlst/list_samples_per_cgmlst_scheme.log'
    threads: config['threads']['other']
    resources: mem_gb=config['mem_gb']['other']
    shell: 
        """
mkdir -p {output}
python bin/chewbbaca_input_files.py --sample-sheet {input} \
                                --output-dir {output} &> {log}
# If no file was produced (namely because none of the genera are supported)
# then a file is produced saying that no genera were supported.
# This is necessary to not break next steps.
if [ -z "$(ls -A {output})" ]; then
    echo "None of the genera of the samples in this run are supported for cgMLST" > "{output}/unsupported_genus_samples.txt"
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
        db_dir = CGMLST_DB
    shell:
        """
if [ {wildcards.genus} == "unsupported_genus" ]; then
    touch {output}
else
    bash bin/chewbbaca_per_genus.sh {input} \
                                {threads} \
                                {params.output_dir} \
                                {params.db_dir} \
                                {wildcards.genus} &> {log}
fi
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