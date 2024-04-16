# ## SEROTYPE ACCORDING TO GENUS ##

# --------------- Choose serotyper according to genus/species ----------------#


def choose_serotyper(wildcards):
    if SAMPLES[wildcards.sample]["genus"] == "salmonella":
        return [
            OUT + "/serotype/{sample}/SeqSero_result_with_context.tsv",
            OUT + "/serotype/{sample}/SeqSero_extra_hits.csv",
        ]
    elif (
        SAMPLES[wildcards.sample]["genus"] == "escherichia"
        or SAMPLES[wildcards.sample]["genus"] == "shigella"
    ):
        return [
            OUT + "/serotype/{sample}/data.json",
            OUT + "/serotype/{sample}/result_serotype.csv",
            OUT + "/serotype/{sample}/shigatyper.csv",
            OUT + "/serotype/{sample}/command.txt",
        ]
    elif SAMPLES[wildcards.sample]["genus"] == "streptococcus":
        return [OUT + "/serotype/{sample}/pred.tsv"]
    elif SAMPLES[wildcards.sample]["genus"] == "neisseria":
        # TODO is the output folder enough
        return [OUT + "/serotype/{sample}"]
    elif SAMPLES[wildcards.sample]["genus"] == "bordetella":
        return [OUT + "/vaccine_antigen_mlst/{sample}.tsv"]
    else:
        return OUT + "/serotype/{sample}/no_serotype_necessary.txt"


# -----------------------------------------------------------------------------#
# This is just a mock rule to make the multiserotypers work
# Similar to the aggregate samples of https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
rule aggregate_serotypes:
    input:
        choose_serotyper,
    output:
        temp(OUT + "/serotype/{sample}_done.txt"),
    message:
        "Checking correct serotyper ran properly for {wildcards.sample}"
    threads: 1
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        "touch {output}"


# -----------------------------------------------------------------------------#
### Salmonella serotyper ###


rule salmonella_serotyper:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
    output:
        seqsero=OUT + "/serotype/{sample}/SeqSero_result.tsv",
        seqsero_tmp1=temp(OUT + "/serotype/{sample}/SeqSero_result.txt"),
        seqsero_tmp2=temp(OUT + "/serotype/{sample}/data_log.txt"),
    message:
        "Running Salmonella serotyper for {wildcards.sample}."
    log:
        OUT + "/log/serotype/{sample}_salmonella.log",
    params:
        output_dir=OUT + "/serotype/{sample}/",
    threads: config["threads"]["seqsero2"]
    resources:
        mem_gb=config["mem_gb"]["seqsero2"],
    conda:
        "../../envs/seqsero.yaml"
    shell:
        """
        # Run seqsero2 
        # -m 'a' means microassembly mode and -t '2' refers to separated fastq files (no interleaved)
        SeqSero2_package.py -m 'a' -t '2' -i {input.r1} {input.r2} -d {params.output_dir} -p {threads} &> {log}
        """


rule add_context_salmonella_serotyper:
    input:
        seqsero=OUT + "/serotype/{sample}/SeqSero_result.tsv",
    output:
        seqsero=OUT + "/serotype/{sample}/SeqSero_result_with_context.tsv",
    message:
        "Adding context to salmonella serotype report for {wildcards.sample}"
    log:
        OUT + "/log/add_context_salmonella_serotyper/{sample}.log",
    params:
        seqsero_context=config["seqsero_context"],
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    conda:
        "../../envs/python.yaml"
    shell:
        """
        python bin/add_context_seqsero.py \
            --input {input.seqsero} \
            --output {output.seqsero} \
            --context {params.seqsero_context} \
            --verbose 2>&1>{log}
        """


rule convert_blastxml_to_csv:
    input:
        seqsero=OUT + "/serotype/{sample}/SeqSero_result.tsv",
    output:
        seqsero=OUT + "/serotype/{sample}/SeqSero_extra_hits.csv",
    message:
        "Converting and filtering blasted_output.xml for {wildcards.sample}"
    log:
        OUT + "/log/convert_blastxml_to_csv/{sample}.log",
    params:
        mincov=0.6,
        minid=0.8,
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    conda:
        "../../envs/python.yaml"
    shell:
        """
# missing blasted_output.xml file does not mean the analysis failed,
# so we need to check manually if the file exists

BLAST_XML=$(dirname {input.seqsero})/blasted_output.xml
if [ -f $BLAST_XML ]
then
    python bin/convert_blastxml_to_csv.py \
    $BLAST_XML \
    {output.seqsero} \
    --minid {params.minid} \
    --mincov {params.mincov} \
    --verbose 2>&1>{log}
else
    touch {output.seqsero}
fi
        """


# -----------------------------------------------------------------------------#
### E. coli serotyper ###


rule ecoli_serotyper:
    input:
        assembly=lambda wildcards: SAMPLES[wildcards.sample]["assembly"],
    output:
        json=OUT + "/serotype/{sample}/data.json",
        csv=OUT + "/serotype/{sample}/result_serotype.csv",
    message:
        "Running E. coli serotyper for {wildcards.sample}."
    log:
        OUT + "/log/serotype/{sample}_ecoli.log",
    conda:
        "../../envs/serotypefinder.yaml"
    threads: config["threads"]["serotypefinder"]
    resources:
        mem_gb=config["mem_gb"]["serotypefinder"],
    params:
        ecoli_db=config["serotypefinder_db"],
        min_cov=config["serotypefinder"]["min_cov"],
        identity_thresh=config["serotypefinder"]["identity_thresh"],
        output_dir=OUT + "/serotype/{sample}/",
    shell:
        """
        python bin/serotypefinder/serotypefinder.py -i {input.assembly} \
            -o {params.output_dir} \
            -p {params.ecoli_db} \
            -l {params.min_cov} \
            -t {params.identity_thresh} &> {log}

        python bin/serotypefinder/extract_alleles_serotypefinder.py {output.json} {output.csv} &>> {log}
        """


# -----------------------------------------------------------------------------#
### Streptococcus pneumoniae serotyper ###


rule build_seroba_db:
    output:
        config["seroba_db"] + "/database/kmer_size.txt",
    conda:
        "../../envs/seroba.yaml"
    params:
        seroba_db=config["seroba_db"],
        kmer_size=config["seroba"]["kmer_size"],
    shell:
        """
        cd {params.seroba_db}
        seroba createDBs database {params.kmer_size}
        """


rule seroba:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        check_db=config["seroba_db"] + "/database/kmer_size.txt",
    output:
        OUT + "/serotype/{sample}/pred.tsv",
    message:
        "Running S. pneumoniae serotyper for {wildcards.sample}."
    log:
        OUT + "/log/serotype/{sample}_spneumoniae.log",
    conda:
        "../../envs/seroba.yaml"
    threads: config["threads"]["seroba"]
    resources:
        mem_gb=config["mem_gb"]["seroba"],
    params:
        min_cov=config["seroba"]["min_cov"],
        seroba_db=config["seroba_db"],
    shell:
        """
        rm -rf {wildcards.sample} 
        OUTPUT_DIR=$(dirname {output})
        mkdir -p $OUTPUT_DIR

        seroba runSerotyping --coverage {params.min_cov} {params.seroba_db}/database {input.r1} {input.r2} {wildcards.sample} &> {log}

        mv {wildcards.sample}/* $OUTPUT_DIR
        """


# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
### Shigella serotyper ###


rule shigatyper:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
    output:
        sample_out=OUT + "/serotype/{sample}/shigatyper.csv",
        command_out=OUT + "/serotype/{sample}/command.txt",
    message:
        "Running Shigella serotyper for {wildcards.sample}."
    log:
        OUT + "/log/serotype/{sample}_shigella.log",
    conda:
        "../../envs/shigatyper.yaml"
    resources:
        mem_gb=config["mem_gb"]["shigatyper"],
    params:
        output_dir=OUT + "/serotype/{sample}",
    shell:
        """
        CURRENT_DIR=$(pwd)
        cd "{params.output_dir}"

        shigatyper --R1 {input.r1} --R2 {input.r2} | grep -E -A 1 '^sample' > "$(basename {output.command_out})" 2> {log}

        if [ -f {wildcards.sample}.csv ]
        then
            mv {wildcards.sample}.csv shigatyper.csv
        else
            # save header without data in expected output file
            echo ",Hit,Number of reads,Length Covered,reference length,% covered,Number of variants,% accuracy" > shigatyper.csv
        fi
        """


# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
### Neisseria serotyper ###


rule characterize_neisseria_capsule:
    input:
        assembly=lambda wildcards: SAMPLES[wildcards.sample]["assembly"],
    output:
        output_dir=directory(OUT + "/serotype/{sample}"),
    message:
        "Running characterize neisseria capsule for {wildcards.sample}."
    log:
        OUT + "/log/serotype/{sample}_neisseria.log",
    conda:
        "../../envs/characterize_neisseria_capsule.yaml"
    resources:
        mem_gb=config["mem_gb"]["characterize_neisseria_capsule"],
    threads: config["threads"]["characterize_neisseria_capsule"]
    params:
        fasta_dir=OUT + "/de_novo_assembly_filtered/",
        output_dir=OUT + "/serotype/{sample}/serogroup",
        copy_to=OUT + "/serotype/{sample}",
    # For this tool we need an input directory and not files, so I copied the filtered assemblies to a new directory, with a directory per sample and use that sample directory as the input
    shell:
        """
        sample=$(awk -F/ '{{print $NF}}' <<< {input.assembly})
        dir_name=$(awk -F. '{{print $1}}' <<< $sample)
        final_name="{params.fasta_dir}$dir_name" 

        mkdir -p $final_name
        cp {input.assembly} "$final_name/"

        python3 bin/characterize_neisseria_capsule/characterize_neisseria_capsule.py -d $final_name -o {output.output_dir}

        cd  {params.output_dir}
        for file in *.tab
        do
            cp "${{file}}" "{params.copy_to}/${{file/*/neisseriatyper.tab}}"
        done
        """


# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
### Bordetella vaccine antigen MLST ###


rule vaccine_antigen_mlst_bordetella:
    input:
        assembly=lambda wildcards: SAMPLES[wildcards.sample]["assembly"],
    output:
        OUT + "/vaccine_antigen_mlst/{sample}.tsv",
    message:
        "Running vaccine antigen mlst for {wildcards.sample}."
    log:
        OUT + "/log/vaccine_antigen_mlst/{sample}_bordetella.log",
    conda:
        "../../envs/tseemann_mlst.yaml"
    resources:
        mem_gb=config["mem_gb"]["tseemann_mlst"],
    threads: config["threads"]["tseemann_mlst"]
    params:
        datadir=config["db_dir"],
        scheme=config["bordetella_vaccine_antigen_scheme"],
        blastdb=config["bordetella_vaccine_antigen_blastdb"],
    shell:
        """
        mlst --nopath --datadir {params.datadir} --blastdb {params.blastdb} --scheme {params.scheme} {input.assembly} > {output} 2> {log}
        """


# -----------------------------------------------------------------------------#

## No serotyper necessary


rule no_serotyper:
    input:
        assembly=lambda wildcards: SAMPLES[wildcards.sample]["assembly"],
    output:
        temp(OUT + "/serotype/{sample}/no_serotype_necessary.txt"),
    message:
        "Skipping serotyper step for {wildcards.sample}."
    threads: 1
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
        touch {output}
        """
