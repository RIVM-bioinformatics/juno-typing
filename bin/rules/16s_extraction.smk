# ---------------------------- 16S extraction --------------------------------#
from Bio import SeqIO

rule barrnap:
    input: 
        assembly = lambda wildcards: SAMPLES[wildcards.sample]['assembly']
    output:
        OUT + '/16s/{sample}/barrnap_result.fasta'
    message: "Running Barrnap for {wildcards.sample}."
    log:
        OUT+'/log/16s/{sample}_barrnap.log'
    conda: 
        '../../envs/16s.yaml'
    threads: config["threads"]["barrnap"]
    resources: mem_gb=config["mem_gb"]["barrnap"]
    shell:
        """
barrnap {input.assembly:q} --outseq {output:q} &> {log:q}
        """

rule extract_16s_from_barrnap:
    input:
        OUT + '/16s/{sample}/barrnap_result.fasta'
    output:
        OUT + '/16s/{sample}/16S_seq.fasta'
    resources: mem_gb=1
    run:
        records = Bio.SeqIO.parse(input[0], 'fasta')
        Bio.SeqIO.write((r for r in records if r.id.startswith('16S_rRNA')), output[0], 'fasta')