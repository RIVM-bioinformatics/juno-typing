# -------------------------- MLST7 multireport -------------------------------#

rule mlst7_multireport:
    input:
        expand(OUT + "/mlst7/{sample}/data.json", sample = SAMPLES)
    output:
        OUT + "/mlst7/mlst7_multireport.csv"
    benchmark:
        OUT+'/log/benchmark/mlst7/mlst7_multireport.txt'
    log:
        OUT+'/log/mlst7/mlst7_multireport.log'
    threads: 1
    resources: mem_gb=config["mem_gb"]["other"]
    shell:
        """
python bin/mslt7_multireport.py  -i {input} -o {output} &> {log}
        """
