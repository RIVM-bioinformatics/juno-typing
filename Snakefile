"""
Juno-typing
Authors: Alejandra Hernandez-Segura, Robert Verhagen, Maaike van der Beld
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium Surveillance (IDS), Bacteriologie (BPD)
Date: 12-01-2021
Documentation: 
Snakemake rules (in order of execution):
    1. CGE-MLST
    2. Bacterial serotyping
        - SeqSero2 for Salmonella
        - SerotypeFinder for E. coli
        - Seroba for Shigella
"""
#################################################################################
##### Import config file, sample_sheet and set output folder names          #####
#################################################################################

configfile: "config/pipeline_parameters.yaml"
configfile: "config/user_parameters.yaml"

from pandas import *
import pathlib
import pprint
import os
import yaml
import json


#################################################################################
#####     Load samplesheet, load genus dict and define output directory     #####
#################################################################################

# Loading sample sheet as dictionary 
# ("R1" and "R2" keys for fastq, and "assembly" for fasta)
SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.safe_load(sample_sheet_file) 

# OUT defines output directory for most rules.
OUT = config["out"]
KMERFINDER_DB=config["kmerfinder_db"]
MLST7_DB = config["mlst7_db"]


#@################################################################################
#@####                              Processes                                #####
#@################################################################################

include: "bin/rules/identify_species.smk"
include: "bin/rules/mlst7_fastq.smk"
include: "bin/rules/mlst7_multireport.smk"
include: "bin/rules/serotype.smk"
include: "bin/rules/serotype_multireports.smk"

#@################################################################################
#@####                    Onstart checker codeblock                          #####
#@################################################################################

onstart:
    try:
        print("Checking if all specified files are accessible...")
        important_files = [ config["sample_sheet"], MLST7_DB, KMERFINDER_DB ]
        for filename in important_files:
            if not os.path.exists(filename):
                raise FileNotFoundError(filename)
    except FileNotFoundError as e:
        print("This file is not available or accessible: %s" % e)
        sys.exit(1)
    else:
        print("\tAll specified files are present!")
    shell("""
        mkdir -p {OUT}
        mkdir -p {OUT}/audit_trail
        mkdir -p {OUT}/log/drmaa
        LOG_CONDA="{OUT}/audit_trail/log_conda.txt"
        LOG_CONFIG="{OUT}/audit_trail/log_config.txt"
        LOG_GIT="{OUT}/audit_trail/log_git.txt"
        echo -e "\nLogging pipeline settings..."
        echo -e "This is the link to the code used for this analysis:"  > "${{LOG_GIT}}"
        echo -e "\thttps://github.com/AleSR13/Juno-typing/tree/$(git log -n 1 --pretty=format:"%H")" >> "${{LOG_GIT}}"
        echo -e "\tGenerating full software list of current Conda environment ..."
        echo -e "Master environment:\n\n" > "${{LOG_CONDA}}"
        conda list >> "${{LOG_CONDA}}"
        echo -e "\tGenerating config file log..."
        rm -f "${{LOG_CONFIG}}"
        echo -e "Date run: $(date)" > "${{LOG_CONFIG}}"
        touch "${{LOG_CONFIG}}"
        for file in config/*.yaml
        do
            echo -e "\n==> Contents of file \"${{file}}\": <==" >> "${{LOG_CONFIG}}"
            cat ${{file}} >> "${{LOG_CONFIG}}"
            echo -e "\n\n" >> "${{LOG_CONFIG}}"
        done
        echo -e "\n\nSample sheet:\n" >> "${{LOG_CONFIG}}"
        cat sample_sheet.yaml >> "${{LOG_CONFIG}}"
    """)

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
        snakemake --profile config --cores 1 --unlock
        snakemake --profile config --cores 1 --report '{OUT}/audit_trail/snakemake_report.html'
        echo -e "Juno-typing finished successfully!"
         """)


#################################################################################
#####                       Specify final output                            #####
#################################################################################

localrules:
    all,
    no_serotyper,
    serotype_multireports,
    mlst7_multireport

rule all:
    input:
        expand(OUT + "/mlst7/{sample}/results.txt", sample = SAMPLES),
        OUT+'/serotype/salmonella_serotype_multireport.csv',
        OUT + '/serotype/ecoli_serotype_multireport.csv',
        OUT + '/serotype/spneumoniae_serotype_multireport.csv',
        OUT + "/mlst7/mlst7_multireport.csv"


