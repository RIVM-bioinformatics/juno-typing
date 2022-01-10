# ## SEROTYPE ACCORDING TO GENUS ##

# --------------- Choose serotyper according to genus/species ----------------#

def choose_serotyper(wildcards):
    if SAMPLES[wildcards.sample]['genus'] == 'salmonella':
        return [OUT+'/serotype/{sample}/SeqSero_result.tsv']
    elif SAMPLES[wildcards.sample]['genus'] == 'escherichia' or SAMPLES[wildcards.sample]['genus'] == 'shigella':
        return [OUT + '/serotype/{sample}/data.json',
                OUT + '/serotype/{sample}/result_serotype.csv',
                OUT + '/serotype/{sample}/shigatyper.csv',
                OUT + '/serotype/{sample}/command.txt']
    elif SAMPLES[wildcards.sample]['genus'] == 'streptococcus':
        return [OUT + '/serotype/{sample}/pred.tsv']
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
    message: "Checking correct serotyper ran properly for {wildcards.sample}"
    threads: 1
    resources: mem_gb=config["mem_gb"]["other"]
    shell:
        'touch {output}'

#-----------------------------------------------------------------------------#
### Salmonella serotyper ###

rule salmonella_serotyper:
    input:
        r1 = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output:
        seqsero = OUT+'/serotype/{sample}/SeqSero_result.tsv',
        seqsero_tmp1 = temp(OUT+'/serotype/{sample}/SeqSero_result.txt'),
        seqsero_tmp2 = temp(OUT+'/serotype/{sample}/blasted_output.xml'),
        seqsero_tmp3 = temp(OUT+'/serotype/{sample}/data_log.txt')
    message: "Running Salmonella serotyper for {wildcards.sample}."
    log:
        OUT+'/log/serotype/{sample}_salmonella.log'
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
        """

#-----------------------------------------------------------------------------#
### E. coli serotyper ###

rule ecoli_serotyper:
    input: 
        assembly = lambda wildcards: SAMPLES[wildcards.sample]['assembly']
    output: 
        json = OUT + '/serotype/{sample}/data.json',
        csv = OUT + '/serotype/{sample}/result_serotype.csv'
    message: "Running E. coli serotyper for {wildcards.sample}."
    log:
        OUT+'/log/serotype/{sample}_ecoli.log'
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
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output:
        OUT + "/serotype/{sample}/pred.tsv"
    message: "Running S. pneumoniae serotyper for {wildcards.sample}."
    log:
        OUT+'/log/serotype/{sample}_spneumoniae.log'
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


rule shigatyper:
    input:
        r1 = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output:
        sample_out = OUT + '/serotype/{sample}/shigatyper.csv',
        command_out = OUT + '/serotype/{sample}/command.txt'
    message: "Running Shigella serotyper for {wildcards.sample}."
    log:
        OUT+'/log/serotype/{sample}_shigella.log'
    conda:
        "../../envs/shigatyper.yaml"
    resources: mem_gb=config["mem_gb"]["shigatyper"]
    params:
        output_dir = OUT + "/serotype/{sample}"
    shell:
        """
CURRENT_DIR=$(pwd)
cd "{params.output_dir}"

shigatyper {input.r1} {input.r2} > "$(basename {output.command_out})" 2> {log}

for file in *.csv
do 
    file_stem=$(echo ${{file%\.*}})
    if [[ {wildcards.sample} = *"$file_stem"* ]]; then
        mv "${{file}}" "${{file/*/shigatyper.csv}}"
    fi
done
        """

#-----------------------------------------------------------------------------#

## No serotyper necessary

rule no_serotyper:
    input: 
        assembly = lambda wildcards: SAMPLES[wildcards.sample]['assembly']
    output: 
        temp(OUT + "/serotype/{sample}/no_serotype_necessary.txt")
    message: "Skipping serotyper step for {wildcards.sample}."
    threads: 1
    resources: mem_gb=config["mem_gb"]["other"]
    shell:
        """
touch {output}
        """