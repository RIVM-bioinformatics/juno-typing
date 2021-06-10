"""
Juno-typing
Authors: Alejandra Hernandez-Segura, Maaike van der Beld
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium Surveillance (IDS), Bacteriologie (BPD)
Date: 12-01-2021
Documentation: 
Snakemake rules (in order of execution):
    1. CGE-MLST
    2. Bacterial serotyping
        - SeqSero2 for Salmonella
        - SerotypeFinder for E. coli
        - Seroba for S. pneumoniae
"""
#################################################################################
##### Import config file, sample_sheet and set output folder names          #####
#################################################################################

import os
import yaml


#################################################################################
#####     Load samplesheet, load genus dict and define output directory     #####
#################################################################################

# Loading sample sheet as dictionary 
# ("R1" and "R2" keys for fastq, and "assembly" for fasta)
sample_sheet = config["sample_sheet"]
SAMPLES = {}
with open(sample_sheet) as sample_sheet_file:
    SAMPLES = yaml.safe_load(sample_sheet_file) 

# OUT defines output directory for most rules.
OUT = config["out"]


#@################################################################################
#@####                              Processes                                #####
#@################################################################################

include: "bin/rules/identify_species.smk"
include: "bin/rules/mlst7_fastq.smk"
include: "bin/rules/mlst7_multireport.smk"
include: "bin/rules/serotype.smk"
include: "bin/rules/serotype_multireports.smk"

#@################################################################################
#@####              Finalize pipeline (error/success)                        #####
#@################################################################################

onerror:
    shell("""
find -maxdepth 1 -type d -empty -exec rm -rf {{}} \;
echo -e "Something went wrong with Juno-typing pipeline. Please check the logging files in {OUT}/log/"
    """)


onsuccess:
    shell("""
        find -maxdepth 1 -type d -empty -exec rm -rf {{}} \;
        find {OUT}/serotype -type f -empty -exec rm {{}} \;
        echo -e "\tGenerating Snakemake report..."
        snakemake --config sample_sheet={sample_sheet} \
                    --configfile config/pipeline_parameters.yaml config/user_parameters.yaml \
                    --cores 1 --unlock
        snakemake --config sample_sheet={sample_sheet} \
                    --configfile config/pipeline_parameters.yaml config/user_parameters.yaml \
                    --cores 1 --report '{OUT}/audit_trail/snakemake_report.html'
        echo -e "Juno-typing finished successfully!"
        """)


#################################################################################
#####                       Specify final output                            #####
#################################################################################

localrules:
    all,
    no_serotyper

rule all:
    input:
        expand(OUT + "/mlst7/{sample}/results.txt", sample = SAMPLES),
        OUT+'/serotype/salmonella_serotype_multireport.csv',
        OUT + '/serotype/ecoli_serotype_multireport.csv',
        OUT + '/serotype/spneumoniae_serotype_multireport.csv',
        OUT + "/mlst7/mlst7_multireport.csv"


