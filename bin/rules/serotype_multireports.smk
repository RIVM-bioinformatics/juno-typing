# ----------------------- Serotypers multireport -----------------------------#


rule serotype_multireports:
    input:
        expand(OUT + "/serotype/{sample}_done.txt", sample=SAMPLES),
    output:
        OUT + "/serotype/serotyper_multireport.csv",
    message:
        "Making multireport(s) for the serotyping results (if any)."
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    log:
        OUT + "/log/serotype/serotype_multireport.log",
    params:
        output_dir=OUT + "/serotype",
    shell:
        """
        input_serotype=""
        for subfolder in $(ls {params.output_dir})
        do
            sample_subfolder="{params.output_dir}/${{subfolder}}"
            result_sample=$(find "${{sample_subfolder}}" \
                    -type f \
                    -name "SeqSero_result_with_context.tsv" \
                    -o -name "result_serotype.csv" \
                    -o -name "command.txt" \
                    -o -name "shigatyper.csv" \
                    -o -name "neisseriatyper.tab" \
                    -o -name "pred.tsv")
            echo $result_sample &> {log}
            input_serotype="${{input_serotype}} ${{result_sample}}"
        done

        echo "Input serotype:"
        echo $input_serotype &>> {log}

        if [[ -z "${{input_serotype// }}" ]]; then
            touch {output}
        else
            python bin/serotyper_multireport.py -i ${{input_serotype}} -o {params.output_dir} &>> {log}
        fi
        """
