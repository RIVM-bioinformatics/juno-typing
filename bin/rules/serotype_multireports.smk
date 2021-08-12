# ----------------------- Serotypers multireport -----------------------------#

rule serotype_multireports:
    input:
        expand(OUT+'/serotype/{sample}_done.txt', sample = SAMPLES)
    output:
<<<<<<< HEAD
        salmonella = OUT + '/serotype/salmonella_serotype_multireport.csv',
        ecoli = OUT + '/serotype/ecoli_serotype_multireport.csv',
        spneumoniae = OUT + '/serotype/spneumoniae_serotype_multireport.csv'        
    # benchmark:
    #     OUT+'/log/benchmark/serotype/serotype_multireport.txt'
    threads: 1
=======
        OUT + '/serotype/serotyper_multireport.csv'
    threads: config["threads"]["other"]
>>>>>>> 9f0c35587ec974eca5e4d64feacc9ca05f72253c
    resources: mem_gb=config["mem_gb"]["other"]
    log:
        OUT+'/log/serotype/serotype_multireport.log'
    params:
        output_dir = OUT + '/serotype'
    shell:
        """
input_serotype=""
for subfolder in $(ls {params.output_dir})
do
    sample_subfolder="{params.output_dir}/${{subfolder}}"
    result_sample=$(find "${{sample_subfolder}}" \
            -type f \
            -name "final_salmonella_serotype.tsv" \
            -o -name "result_serotype.csv" \
            -o -name "pred.tsv")
    echo $result_sample &> {log}
    input_serotype="${{input_serotype}} ${{result_sample}}"
done

echo "Input serotype:"
echo $input_serotype &>> {log}

if [ "$input_serotype" == "" ]; then
    touch {output}
else
    python bin/serotyper_multireport.py -i ${{input_serotype}} -o {params.output_dir} &>> {log}
fi
        """
