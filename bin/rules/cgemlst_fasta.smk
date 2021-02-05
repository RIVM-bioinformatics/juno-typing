#############################################################################
#####                               MLST7                               #####
#############################################################################

rule cgemlst:
    input:
        assembly = lambda wildcards: SAMPLES[wildcards.sample]["assembly"],
        db_dir = OUT + "/checker_cgemlst_builddb.txt"
    output:
        OUT + "/{sample}/results_tab.tsv"
    conda:
        "../../envs/cgemlst.yaml"
    log:
        OUT + "/log/cgemlst/{sample}.log"
    benchmark:
        OUT + "/log/benchmark/cgemlst_{sample}.txt"
    threads: config["threads"]["cgemlst"]
    resources: mem_mb=config["mem_mb"]["cgemlst"]
    params:
        species = lambda wildcards: SAMPLES[wildcards.sample]["species"]
    shell:
        """
python bin/cge-mlst/mlst.py -i {input.r1} {input.r2} \
-o $(dirname {output}) \
-s senterica \
--database $(dirname {input.db_dir})/mlst_db \
-x 2> {log}
        """

