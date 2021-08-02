#############################################################################
#####                          Identify species                         #####
#############################################################################

rule identify_species:
    input:
        r1 = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        db = config["kmerfinder_db"] + "/bacteria/bacteria.ATG.length.b"
    output:
        kmerfinder = OUT + "/identify_species/{sample}/data.json"
    log:
        OUT + "/log/identify_species/{sample}_identify_species.log"
    benchmark:
        OUT + "/log/benchmark/identify_species_{sample}.txt"
    threads: config["threads"]["kmerfinder"]
    resources: mem_gb=config["mem_gb"]["kmerfinder"]
    shell:
        """
DB_DIR=$(dirname {input.db})

python bin/kmerfinder/kmerfinder.py -i "{input.r1}" "{input.r2}" \
-o "$(dirname {output.kmerfinder})" \
-db "${{DB_DIR}}/bacteria.ATG" \
-tax "${{DB_DIR}}/bacteria.tax" \
-x 2> {log}

if `cat {output} | grep -q "species_hits': {{}}"`; then
    echo "No species were detected." 2> {log}
    rm -f {output}
fi
        """

# KmerFinder often makes errors in which it finds no hits. 
# These errors cause a problem in the following steps of the pipeline.
# However, often just re-running KmerFinder fixes it.
# TODO: Find and fix the reason why it fails instead of just doing this quick fix.


checkpoint which_species:
    input:
        OUT + "/identify_species/{sample}/data.json"
    output:
        OUT + "/identify_species/{sample}/best_species_hit.txt"
    log:
        OUT + "/log/identify_species/{sample}_which_species.log"
    threads: 
        1
    resources: 
        mem_gb=config["mem_gb"]["other"]
    shell:
        """
python bin/get_species.py {input} > {output} 2> {log}
        """