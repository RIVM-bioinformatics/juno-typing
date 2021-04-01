# ---------Rscript to summarize results for all samples ----------------------#

# ---------- Choose serotyper and make the multireport accordingly -----------#
def choose_serotyper(wildcards):
    with checkpoints.which_species.get(sample=wildcards.sample).output[0].open() as f:
        species_res = f.read().strip()
        is_salmonella = species_res.find("salmonella") != -1
        is_ecoli = species_res.find("escherichia") != -1
        if is_salmonella:
            return OUT+'/serotype/{sample}/SeqSero_result.tsv'
        elif is_ecoli:
            return OUT + "/serotype/{sample}/serotypefinder_result.txt"
        else:
            return OUT + "/serotype/{sample}/serotypefinder_result.txt"

rule aggregate_serotypes:
    input:
        choose_serotyper
    output:
        temp(OUT+'/serotype/{sample}_done.txt')
    threads: 1
    resources: mem_mb=2000
    shell:
        'touch {output}'

# ----------------------- Salmonella multireport -----------------------------#

rule salmonella_serotype_multireport:
    input:
        expand(OUT+'/serotype/{sample}_done.txt', sample = SAMPLES)
    output:
        OUT+'/serotype/salmonella_serotype_multireport.csv'
    benchmark:
        OUT+'/log/benchmark/serotype_salmonella/salmonella_serotype_multireport.txt'
    log:
        OUT+'/log/serotype_salmonella/salmonella_serotype_multireport.log'
    shell:
        """
    INPUT='{input}'
    declare -a StringArray=(${{INPUT//_done\.txt/\/SeqSero_result.tsv}})
    # Iterate the string array using for loop
    SEROTYPE_INPUT=""
    for val in ${{StringArray[@]}}; do
        if [ -f $val ]; then
            SEROTYPE_INPUT="${{SEROTYPE_INPUT}} ${{val}}"
        fi
    done 2> {log}

    if [[ $test == "" ]]; then
        touch {output}
    else
        python bin/seqsero2_multireport.py -i $SEROTYPE_INPUT -o {output} 2> {log}
    fi
        """



