# ---------------------------- 16S extraction --------------------------------#


rule barrnap:
    input:
        assembly=lambda wildcards: SAMPLES[wildcards.sample]["assembly"],
    output:
        OUT + "/16s/{sample}/barrnap_result.fasta",
        temp(OUT + "/16s/{sample}/temp_barrnap_input.fasta"),
        temp(OUT + "/16s/{sample}/temp_barrnap_input.fasta.fai"),
    message:
        "Running Barrnap for {wildcards.sample}."
    log:
        OUT + "/log/16s/{sample}_barrnap.log",
    conda:
        "../../envs/16s.yaml"
    threads: config["threads"]["barrnap"]
    resources:
        mem_gb=config["mem_gb"]["barrnap"],
    shell:
        """
        cp {input.assembly:q} {output[1]:q}
        barrnap {output[1]:q} --outseq {output[0]:q} 2>&1 {log:q}
        """


from Bio import SeqIO


rule extract_16s_from_barrnap:
    input:
        OUT + "/16s/{sample}/barrnap_result.fasta",
    output:
        OUT + "/16s/{sample}/16S_seq.fasta",
    resources:
        mem_gb=1,
    run:
        records = SeqIO.parse(input[0], "fasta")
        SeqIO.write(
            (r for r in records if r.id.startswith("16S_rRNA")), output[0], "fasta"
        )
