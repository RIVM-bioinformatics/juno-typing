#############################################################################
#####                               MLST7                               #####
#############################################################################

rule mlst7:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        db=config["mlst7_db"] + "/senterica/senterica.length.b",
    output:
        json=temp(OUT + "/mlst7/{sample}/data.json"),
        txt=OUT + "/mlst7/{sample}/results.txt",
        fasta=OUT + "/mlst7/{sample}/MLST_allele_seq.fsa",
        hits=temp(OUT + "/mlst7/{sample}/Hit_in_genome_seq.fsa"),
        tab=temp(OUT + "/mlst7/{sample}/results_tab.tsv"),
    message:
        "Calculating the 7 locus-MLST for {wildcards.sample}"
    conda:
        "../../envs/mlst7.yaml"
    log:
        OUT + "/log/mlst7/{sample}.log",
    threads: config["threads"]["cgemlst"]
    resources:
        mem_gb=config["mem_gb"]["cgemlst"],
    params:
        species=lambda wildcards: SAMPLES[wildcards.sample]["species-mlst7"],
        mlst7_db=config["mlst7_db"],
    shell:
        """
        if [ {params.species} == 'None' ]
        then
            echo -e "The species of this sample is not supported by the MLST7 tool." > {log}
            touch {output}
            cp files/no_mlst7.json {output.json}
        else
            python bin/cge-mlst/mlst.py -i {input.r1} {input.r2} \
            -o $(dirname {output.json}) \
            -s {params.species} \
            --database {params.mlst7_db} \
            -mp kma \
            -x &>> {log}
        fi
        """

rule mlst7_nanopore:
    input:
        # fastq=lambda wildcards: SAMPLES[wildcards.sample]["nanopore_input"],
        assembly=lambda wildcards: SAMPLES[wildcards.sample]["assembly"],
        db=config["mlst7_db"] + "/senterica/senterica.length.b",
    output:
        json=temp(OUT + "/mlst7/{sample}/data.json"),
        txt=(OUT + "/mlst7/{sample}/results.txt"),
        fasta=(OUT + "/mlst7/{sample}/MLST_allele_seq.fsa"),
        hits=temp(OUT + "/mlst7/{sample}/Hit_in_genome_seq.fsa"),
        tab=temp(OUT + "/mlst7/{sample}/results_tab.tsv"),
        tmp_mlst=directory(OUT + "/mlst7/{sample}/temp")
    message:
        "Calculating the 7 locus-MLST for {wildcards.sample}"
    conda:
        "../../envs/mlst7.yaml"
    log:
        OUT + "/log/mlst7/{sample}.log",
    threads: config["threads"]["cgemlst"]
    resources:
        mem_gb=config["mem_gb"]["cgemlst"],
    params:
        species=lambda wildcards: SAMPLES[wildcards.sample]["species-mlst7"],
        mlst7_db=config["mlst7_db"],
    shell:
        """
        if [ {params.species} == 'None' ]
        then
            echo -e "The species of this sample is not supported by the MLST7 tool." > {log}
            touch {output.json}
            cp files/no_mlst7.json {output.json}
        else
            python bin/cge-mlst/mlst.py -i {input.assembly} \
            -o $(dirname {output.json}) \
            -s {params.species} \
            --database {params.mlst7_db} \
            -mp blastn \
            -t {output.tmp_mlst}\
            -x &>> {log}
        fi
        """

# rule mlst7_nanopore:
#     input:
#         # fastq=lambda wildcards: SAMPLES[wildcards.sample]["nanopore_input"],
#         assembly=lambda wildcards: SAMPLES[wildcards.sample]["assembly"],
#     output:
#         mlst=OUT + "/mlst7/{sample}/mlst_tseemann.tsv",
        
#     message:
#         "Calculating the 7 locus-MLST for {wildcards.sample}"
#     conda:
#         "../../envs/tseemann_mlst.yaml"
#     log:
#         OUT + "/log/mlst7/{sample}.log",
#     threads: config["threads"]["tseemann_mlst"]
#     resources:
#         mem_gb=config["mem_gb"]["tseemann_mlst"],
#     params:
#         species=lambda wildcards: SAMPLES[wildcards.sample]["species-mlst7"],
#         # mlst7_db=config["mlst7_db"],
#     shell:
#         """
#         mlst {input.assembly}  > {output.mlst} 2>> {log}
#         """

# rule parse_mlst7_tseemann:
#     input:
#         OUT + "/mlst7/{sample}/mlst_tseemann.tsv",
#     output:
#         OUT + "/mlst7/{sample}/mlst_report.tsv",
#     threads: config["threads"]["tseemann_mlst"]
#     resources:
#         mem_gb=config["mem_gb"]["tseemann_mlst"],
#     run:
#         import pandas as pd
#         import re
        
#         rows = []
#         genes = set()

#         with open(input[0]) as f:
#             for line in f:
#                 fields = line.strip().split("\t")

#                 row = {
#                     "sample": fields[0].split("/")[-1].split(".")[0],
#                     "scheme": fields[1],
#                     "ST": fields[2],
#                 }

#                 for allele in fields[3:]:
#                     m = re.match(r"(.+)\((.+)\)", allele)
#                     if m:
#                         gene, allele_number = m.groups()
#                         row[gene] = allele_number
#                         genes.add(gene)
#                 rows.append(row)
        
#         cols = ["sample", "scheme", "ST"] + sorted(genes)
#         pd.DataFrame(rows).reindex(columns=cols).to_csv(output[0], sep="\t", index=False)
