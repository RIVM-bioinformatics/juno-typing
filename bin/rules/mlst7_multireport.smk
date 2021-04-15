# -------------------------- MLST7 multireport -------------------------------#

rule mlst7_multireport:
    input:
        expand(OUT + "/mlst7/{sample}/results_tab.tsv", sample = SAMPLES)
    output:
        OUT + "/mlst7/mlst7_multireport.csv"
    benchmark:
        OUT+'/log/benchmark/mlst7/mlst7_multireport.txt'
    log:
        OUT+'/log/mlst7/mlst7_multireport.log'
    shell:
        """
python bin/mslt7_multireport.py  -i {input} -o {output}
        """
