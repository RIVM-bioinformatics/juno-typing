"""
Juno-typing
Author(s): Alejandra Hernandez-Segura, Kaitlin Weber, Edwin van der Kind and Maaike van den Beld
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
    3. Multireports for serotyper and CGE-MLST
"""
#################################################################################
##### Import config file, sample_sheet and set output folder names          #####
#################################################################################

from os.path import getsize, exists, abspath
from yaml import safe_load

#################################################################################
#####     Load samplesheet, load genus dict and define output directory     #####
#################################################################################

# Loading sample sheet as dictionary 
# ("R1" and "R2" keys for fastq, and "assembly" for fasta)
sample_sheet = config["sample_sheet"]
SAMPLES = {}
with open(sample_sheet) as sample_sheet_file:
    SAMPLES = safe_load(sample_sheet_file) 

# OUT defines output directory for most rules.
OUT = config["out"]


#@################################################################################
#@####                              Processes                                #####
#@################################################################################

include: "bin/rules/mlst7_fastq.smk"
include: "bin/rules/mlst7_multireport.smk"
include: "bin/rules/serotype.smk"
include: "bin/rules/serotype_multireports.smk"

#@################################################################################
#@####              Finalize pipeline (error/success)                        #####
#@################################################################################

# TODO: eventually these files should be stored somewhere else and included in the pipeline as tmp files
onerror:
    shell("""
find -maxdepth 1 -type d -empty -exec rm -rf {{}} \;
find -maxdepth 1 -type f -name "*.depth.txt*" -exec rm -rf {{}} \;
echo -e "Something went wrong with Juno-typing pipeline. Please check the logging files in {OUT}/log/"
    """)


#################################################################################
#####                       Specify final output                            #####
#################################################################################

localrules:
    all,
    aggregate_serotypes,
    no_serotyper

rule all:
    input:
        expand(OUT + "/mlst7/{sample}/results.txt", sample = SAMPLES),
        OUT+'/serotype/serotyper_multireport.csv',
        OUT + "/mlst7/mlst7_multireport.csv"


