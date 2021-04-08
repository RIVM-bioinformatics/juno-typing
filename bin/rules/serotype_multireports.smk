# ----------------------- Serotypers multireport -----------------------------#

rule serotype_multireports:
    input:
        expand(OUT+'/serotype/{sample}_done.txt', sample = SAMPLES)
    output:
        salmonella = OUT + '/serotype/salmonella_serotype_multireport.csv',
        ecoli = OUT + '/serotype/ecoli_serotype_multireport.csv'
    benchmark:
        OUT+'/log/benchmark/serotype_salmonella/salmonella_serotype_multireport.txt'
    log:
        OUT+'/log/serotype_salmonella/salmonella_serotype_multireport.log'
    shell:
        """
INPUT=$(tr ' ' $'\n' <<< '{input}')


## Salmonella serotype multireport
declare -a StringArray=(${{INPUT//_done\.txt/\/SeqSero_result.tsv}})
SEROTYPE_INPUT=""
for val in ${{StringArray[@]}}; do
    if [ -f $val ]; then
        SEROTYPE_INPUT="${{SEROTYPE_INPUT}} ${{val}}"
    fi
done 2> {log}

if [[ ${{#SEROTYPE_INPUT}} -eq 0 ]]; then
    touch {output}
else
    python bin/seqsero2_multireport.py -i ${{SEROTYPE_INPUT}} -o {output.salmonella} 2> {log}
fi

# E. coli serotype multireport
declare -a StringArray=(${{INPUT//_done\.txt/\/result_serotype.csv}})
SEROTYPE_INPUT=""
for val in ${{StringArray[@]}}; do
    if [ -f $val ]; then
        SEROTYPE_INPUT="${{SEROTYPE_INPUT}} ${{val}}"
    fi
done 2> {log}

if [[ ${{#SEROTYPE_INPUT}} -eq 0 ]]; then
    touch {output}
else
    python bin/serotypefinder_multireport.py ${{SEROTYPE_INPUT}} {output.ecoli} 2> {log}
fi
        """


