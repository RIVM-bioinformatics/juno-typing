#############################################################################
#####                               MLST7                               #####
#############################################################################

rule mlst7:
    input:
        assembly = lambda wildcards: SAMPLES[wildcards.sample]["assembly"],
        db = config["mlst7_db"] + "/senterica/senterica.length.b",
        species = OUT + "/identify_species/{sample}/best_species_hit.txt"
    output:
        jjson = temp(OUT + "/mlst7/{sample}/data.json"),
        txt = OUT + "/mlst7/{sample}/results.txt",
        fasta = OUT + "/mlst7/{sample}/MLST_allele_seq.fsa",
        hits = temp(OUT + "/mlst7/{sample}/Hit_in_genome_seq.fsa"),
        tab = temp(OUT + "/mlst7/{sample}/results_tab.tsv")
    conda:
        "../../envs/mlst7.yaml"
    log:
        OUT + "/log/mlst7/{sample}.log"
    threads: config["threads"]["cgemlst"]
    resources: mem_gb=config["mem_gb"]["cgemlst"]
    params:
        species = lambda wildcards: SAMPLES[wildcards.sample]["species-mlst7"],
        mlst7_db = config["mlst7_db"]
    shell:
        """
if [ {params.species} == "NotProvided" ]; then
    SPECIES=$(python bin/get_mlst7_db_name.py {input.species}) &> {log}
else
    SPECIES={params.species}
fi

python bin/cge-mlst/mlst.py -i {input.r1} {input.r2} \
-o $(dirname {output}) \
-s $SPECIES \
--database {params.mlst7_db} \
-x &>> {log}
        """

