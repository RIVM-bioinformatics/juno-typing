#############################################################################
#####                               MLST7                               #####
#############################################################################

rule cgemlst:
    input:
        r1 = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        db_dir = OUT + "/checker_cgemlst_builddb.txt"
    output:
        OUT + "/{sample}/results_tab.tsv"
    conda:
        "../../envs/cgemlst.yaml"
    log:
        OUT + "/log/cgemlst/{sample}.log"
    benchmark:
        OUT + "/log/benchmark/cgemlst_{sample}.txt"
    threads: 1
    params:
        species = lambda wildcards: SAMPLES[wildcards.sample]["species"]
    shell:
        """
python bin/cge-mlst/mlst.py -i {input.r1} {input.r2} \
-o $(dirname {output}) \
-s {params.species} \
--database $(dirname {input.db_dir})/mlst_db \
-mp bin/kma/kma \
-x 2> {log}
        """




