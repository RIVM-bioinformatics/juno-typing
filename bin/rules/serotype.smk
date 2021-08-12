# ## SEROTYPE ACCORDING TO GENUS ##

# --------------- Choose serotyper according to genus/species ----------------#

def choose_serotyper(wildcards):
    with checkpoints.which_species.get(sample=wildcards.sample).output[0].open() as f:
        species_res = f.read().strip()
        is_salmonella = species_res.find("salmonella") != -1
        is_ecoli = species_res.find("escherichia") != -1
        is_strepto = species_res.find("streptococcus") != -1
        is_shigella = species_res.find("shigella") != -1
        if is_salmonella:
            return [OUT+'/serotype/{sample}/SeqSero_result.tsv',
                    OUT + "/serotype/{sample}/final_salmonella_serotype.tsv"]
        elif is_ecoli:
            return [OUT + '/serotype/{sample}/data.json',
                    OUT + '/serotype/{sample}/result_serotype.csv']
        elif is_strepto:
            return [OUT + '/serotype/{sample}/pred.tsv']
        elif is_shigella:
            #TODO fill in what to return here for shigella, depends on the tool
            return [OUT + 'serotype/command_line/{sample}.txt']
        else:
            return OUT + "/serotype/{sample}/no_serotype_necessary.txt"

#-----------------------------------------------------------------------------#
# This is just a mock rule to make the multiserotypers work
# Similar to the aggregate samples of https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
rule aggregate_serotypes:
    input:
        choose_serotyper
    output:
        temp(OUT+'/serotype/{sample}_done.txt')
    threads: 1
    resources: mem_gb=config["mem_gb"]["other"]
    shell:
        'touch {output}'

#-----------------------------------------------------------------------------#
### Salmonella serotyper ###

rule salmonella_serotyper:
    input:
        r1 = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        species = OUT + "/identify_species/{sample}/best_species_hit.txt"
    output:
        seqsero = OUT+'/serotype/{sample}/SeqSero_result.tsv',
        seqsero_tmp1 = temp(OUT+'/serotype/{sample}/SeqSero_result.txt'),
        seqsero_tmp2 = temp(OUT+'/serotype/{sample}/blasted_output.xml'),
        seqsero_tmp3 = temp(OUT+'/serotype/{sample}/data_log.txt'),
        final_serotype = OUT + "/serotype/{sample}/final_salmonella_serotype.tsv"
    benchmark:
        OUT+'/log/benchmark/serotype_salmonella/{sample}.txt'
    log:
        OUT+'/log/serotype_salmonella/{sample}.log'
    params:
        output_dir = OUT + '/serotype/{sample}/',
        min_cov = config['salmonellamonophasic']['min_cov']
    threads: 
        config["threads"]["seqsero2"]
    resources: 
        mem_gb=config["mem_gb"]["seqsero2"]
    conda:
        '../../envs/seqsero.yaml'
    shell:
        """
# Run seqsero2 
# -m 'a' means microassembly mode and -t '2' refers to separated fastq files (no interleaved)
SeqSero2_package.py -m 'a' -t '2' -i {input.r1} {input.r2} -d {params.output_dir} -p {threads} &> {log}
python bin/check_salmmonophasic.py -n {wildcards.sample} \
                                    -i {output.seqsero} \
                                    -f {input.r1} \
                                    -r {input.r2} \
                                    -o {output.final_serotype} \
                                    -s {params.min_cov} \
                                    -b {params.min_cov} \
                                    -t {threads} &>> {log}
        """

#-----------------------------------------------------------------------------#
### E. coli serotyper ###

rule ecoli_serotyper:
    input: 
        assembly = lambda wildcards: SAMPLES[wildcards.sample]['assembly'],
        species = OUT + "/identify_species/{sample}/best_species_hit.txt"
    output: 
        json = OUT + '/serotype/{sample}/data.json',
        csv = OUT + '/serotype/{sample}/result_serotype.csv'
    log:
        OUT+'/log/serotype_ecoli/{sample}.log'
    # benchmark:
    #     OUT+'/log/benchmark/serotype_ecoli/{sample}.txt'
    conda: 
        '../../envs/serotypefinder.yaml'
    threads: config["threads"]["serotypefinder"]
    resources: mem_gb=config["mem_gb"]["serotypefinder"]
    params: 
        ecoli_db = config['serotypefinder_db'],
        min_cov = config['serotypefinder']['min_cov'],
        identity_thresh = config['serotypefinder']['identity_thresh'],
        output_dir = OUT + '/serotype/{sample}/'
    shell:
        """
python bin/serotypefinder/serotypefinder.py -i {input.assembly} \
    -o {params.output_dir} \
    -p {params.ecoli_db} \
    -l {params.min_cov} \
    -t {params.identity_thresh} &> {log}

python bin/serotypefinder/extract_alleles_serotypefinder.py {output.json} {output.csv} &>> {log}
        """

#-----------------------------------------------------------------------------#
### Streptococcus pneumoniae serotyper ###

rule seroba:
    input:
        r1 = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        species = OUT + "/identify_species/{sample}/best_species_hit.txt"
    output:
        OUT + "/serotype/{sample}/pred.tsv"
    log:
        OUT+'/log/serotype_spneumoniae/{sample}.log'
    # benchmark:
    #     OUT+'/log/benchmark/serotype_spneumoniae/{sample}.txt'
    conda:
        "../../envs/seroba.yaml"
    threads: config["threads"]["seroba"]
    resources: mem_gb=config["mem_gb"]["seroba"]
    params:
        min_cov = config["seroba"]["min_cov"],
        seroba_db = config["seroba_db"]
    shell:
        """
rm -rf {wildcards.sample} 
OUTPUT_DIR=$(dirname {output})
mkdir -p $OUTPUT_DIR

seroba runSerotyping --coverage {params.min_cov} {params.seroba_db}/database {input.r1} {input.r2} {wildcards.sample} &> {log}

mv {wildcards.sample}/* $OUTPUT_DIR
        """

#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
### Shigella serotyper ###

#TODO change name to tool being used
rule shigatyper:
    input:
        r1 = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        species = OUT + "/identify_species/{sample}/best_species_hit.txt"
    output:
        #TODO give correct output, same as at the start of this file
        #OUT + 'serotype/{sample}/test.csv'
        OUT + 'serotype/command_line/{sample}.txt'
    log:
        OUT+'/log/serotype_shigella/{sample}.log'
    conda:
        #TODO make shigella yaml file of the env, depends on tool
        "../../envs/shigatyper.yaml"
    #TODO threads depending on the tool being used
    #threads: config["threads"]["seroba"]
    #TODO where are these coming from?
    resources: mem_gb=config["mem_gb"]["seroba"]

    shell:
    #TODO add shigella pipeline command
        """
shigatyper {input.r1} {input.r2} > {output}
        """

#-----------------------------------------------------------------------------#

## No serotyper necessary

rule no_serotyper:
    input: 
        assembly = lambda wildcards: SAMPLES[wildcards.sample]['assembly'],
        species = OUT + "/identify_species/{sample}/best_species_hit.txt"
    output: 
        OUT + "/serotype/{sample}/no_serotype_necessary.txt"
    threads: 1
    resources: mem_gb=config["mem_gb"]["other"]
    shell:
        """
touch {output}
        """