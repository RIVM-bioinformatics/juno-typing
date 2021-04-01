# ## SEROTYPE ACCORDING TO GENUS ##

### Salmonella serotyper ###

rule salmonella_serotyper:
    input:
        r1 = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        species = OUT + "/identify_species/{sample}/best_species_hit.txt"
    output:
        OUT+'/serotype/{sample}/SeqSero_result.tsv',
        temp(OUT+'/serotype/{sample}/SeqSero_result.txt'),
        temp(OUT+'/serotype/{sample}/blasted_output.xml'),
        temp(OUT+'/serotype/{sample}/data_log.txt')
    benchmark:
        OUT+'/log/benchmark/serotype_salmonella/{sample}.log'
    log:
        OUT+'/log/serotype_salmonella/{sample}.log'
    params:
        output_dir = OUT + '/serotype/{sample}/'
    threads: 
        config["threads"]["seqsero2"]
    resources: 
        mem_mb=config["mem_mb"]["seqsero2"]
    conda:
        '../../envs/seqsero.yaml'
    shell:
        """
# Run seqsero2 
# -m 'a' means microassembly mode and -t '2' refers to separated fastq files (no interleaved)
SeqSero2_package.py -m 'a' -t '2' -i {input.r1} {input.r2} -d {params.output_dir} -p {threads}
        """



### E. coli serotyper ###

rule ecoli_serotyper:
    input:
        r1 = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        species = OUT + "/identify_species/{sample}/best_species_hit.txt"
    output:
        OUT + "/serotype/{sample}/serotypefinder_result.txt"
    benchmark:
        OUT+'/log/benchmark/serotype_ecoli/{sample}.txt'
    log:
        OUT+'/log/serotype_ecoli/{sample}.log'
    params:
        output_dir = OUT + '/serotype/{sample}/'
    threads: 
        1
    resources: 
        mem_mb=2000
    # conda:
    #     '../../envs/seqsero.yaml'
    shell:
        """
touch {output}
        """
