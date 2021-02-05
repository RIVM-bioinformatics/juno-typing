"""
Juno-typing
Authors: Alejandra Hernandez-Segura, Robert Verhagen, Maaike van der Beld
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium Surveillance (IDS), Bacteriologie (BPD)
Date: 12-01-2021
Documentation: 
Snakemake rules (in order of execution):
    1. CGE-MLST
"""
#################################################################################
##### Import config file, sample_sheet and set output folder names          #####
#################################################################################

configfile: "config/pipeline_parameters.yaml"
#configfile: "config/variables.yaml"

from pandas import *
import pathlib
import pprint
import os
import yaml
import json


#################################################################################
#####     Load samplesheet, load genus dict and define output directory     #####
#################################################################################

# SAMPLES is a dict with sample in the form sample > read number > file. E.g.: SAMPLES["sample_1"]["R1"] = "x_R1.gz"
SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.safe_load(sample_sheet_file) 

# Add species name
SPECIES_ALL = config["species"].split()
for sample in SAMPLES:
    try:
        SAMPLES[sample]["species"]
    except:
        if SPECIES_ALL != "NotProvided" : 
            SAMPLES[sample]["species"] = (SPECIES_ALL[0][0] + SPECIES_ALL[1]).lower()
        else:
            print("""
            You did not provide a species for one or more of the samples. 
            Please use the --species or --metadata arguments when running the pipeline and try again. 
            """)
            sys.exit(1)


# OUT defines output directory for most rules.
OUT = config["out"]

#@################################################################################
#@####                              Processes                                #####
#@################################################################################

include: "bin/rules/cgemlst_setup.smk"
include: "bin/rules/cgemlst_builddb.smk"
include: "bin/rules/cgemlst_fastq.smk"

#@################################################################################
#@####              The `onstart` checker codeblock                          #####
#@################################################################################

onstart:
    try:
        print("Checking if all specified files are accessible...")
        important_files = [ config["sample_sheet"] ]
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
        mkdir -p {OUT}/results
        echo -e "\nLogging pipeline settings..."
        echo -e "\tGenerating methodological hash (fingerprint)..."
        echo -e "This is the link to the code used for this analysis:\thttps://github.com/AleSR13/Juno-typing/tree/$(git log -n 1 --pretty=format:"%H")" > '{OUT}/results/log_git.txt'
        echo -e "This code with unique fingerprint $(git log -n1 --pretty=format:"%H") was committed by $(git log -n1 --pretty=format:"%an <%ae>") at $(git log -n1 --pretty=format:"%ad")" >> '{OUT}/results/log_git.txt'
        echo -e "\tGenerating full software list of current Conda environment (\"juno_mmaster\")..."
        conda list > '{OUT}/results/log_conda.txt'
        echo -e "\tGenerating config file log..."
        rm -f '{OUT}/results/log_config.txt'
        for file in config/*.yaml
        do
            echo -e "\n==> Contents of file \"${{file}}\": <==" >> '{OUT}/results/log_config.txt'
            cat ${{file}} >> '{OUT}/results/log_config.txt'
            echo -e "\n\n" >> '{OUT}/results/log_config.txt'
        done
    """)

#@################################################################################
#@#### These are the conditional cleanup rules                               #####
#@################################################################################

onerror:
    shell("""
rm -rf "{OUT}/mlst_db"
echo -e "Something went wrong with Juno-typing pipeline. Please check the logging files in {OUT}/log/"
    """)


onsuccess:
    shell("""
        echo -e "\tGenerating Snakemake report..."
        rm -rf "{OUT}/mlst_db"
        snakemake --config out={OUT} species="{SPECIES_ALL}" --profile config --unlock
        snakemake --config out={OUT} species="{SPECIES_ALL}" --profile config --report '{OUT}/results/snakemake_report.html'
        echo -e "Juno-typing finished successfully!"
         """)


#################################################################################
##### Specify final output:                                                 #####
#################################################################################

localrules:
    all,
    cgemlst_setup,
    cgemlst_builddb

rule all:
    input:
        expand(OUT + "/{sample}/results_tab.tsv", sample = SAMPLES)


