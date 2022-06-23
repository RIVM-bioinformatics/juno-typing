#############################################################################
#####                               MLST7                               #####
#############################################################################


rule mlst7:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        db=config["mlst7_db"] + "/senterica/senterica.length.b",
    output:
        json=temp(OUT + "/mlst7/{sample}/data.json"),
        txt=OUT + "/mlst7/{sample}/results.txt",
        fasta=OUT + "/mlst7/{sample}/MLST_allele_seq.fsa",
        hits=temp(OUT + "/mlst7/{sample}/Hit_in_genome_seq.fsa"),
        tab=temp(OUT + "/mlst7/{sample}/results_tab.tsv"),
    message:
        "Calculating the 7 locus-MLST for {wildcards.sample}"
    conda:
        "../../envs/mlst7.yaml"
    log:
        OUT + "/log/mlst7/{sample}.log",
    threads: config["threads"]["cgemlst"]
    resources:
        mem_gb=config["mem_gb"]["cgemlst"],
    params:
        species=lambda wildcards: SAMPLES[wildcards.sample]["species-mlst7"],
        mlst7_db=config["mlst7_db"],
    shell:
        """
        if [ {params.species} == 'None' ]
        then
            echo -e "The species of this sample is not supported by the MLST7 tool." > {log}
            touch {output}
            cp files/no_mlst7.json {output.json}
        else
            python bin/cge-mlst/mlst.py -i {input.r1} {input.r2} \
            -o $(dirname {output.json}) \
            -s {params.species} \
            --database {params.mlst7_db} \
            -mp kma \
            -x &>> {log}
        fi
        """
