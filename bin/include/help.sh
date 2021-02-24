#!/bin/bash

PIPELINE_NAME="$1"
VERSION="$2"

echo """
${PIPELINE_NAME}
Version used: 'https://github.com/AleSR13/Juno-typing/tree/${VERSION}'
Built with Snakemake

  Usage: bash $(echo "${PIPELINE_NAME}" | awk '{print tolower($0)}') -i <INPUT_DIR> --species <GENUS SPECIES> <parameters>
  N.B. it is designed for Illumina paired-end data only


Input:
  -i, --input [DIR]                 This is the folder containing your input 
                                    fastq or fasta files. Only one type (either
                                    fastq or fasta) should be present. Default 
                                    is raw_data/
  -o, --output [DIR]                This is the folder containing your output 
                                    fastq files. Default is out/ 

  -db [DIR]                         Path to directory containing the databases
                                    for KmerFinder and CGE-MLST. If it does not
                                    exist or the databases are not present, they
                                    will be automatically downloaded.

  --update-db                       If this flag is present, the databases will 
                                    be updated, even if they are present. Note
                                    that due to a bug in the FTP server from the 
                                    CGE, the KmerFinder database cannot be 
                                    updated. Instead, the same version is always 
                                    used (version '20201202'). 

  --metadata [PATH_TO_FILE]         Excel file (.xlsx) containing the information
                                    of the samples and their corresponding genus. 
                                    It should contain at least a column called 
                                    'Sample' containing the sample ID (same name 
                                    than fastq files but removing the suffix 
                                    (_S##_R1|2.fastq.gz or .fasta), another one 
                                    called 'Genus' containing the name of the genus 
                                    and a last one called 'Species' containing the 
                                    name of the species (without the genus). Mind 
                                    the capital at the beginning of each word 
                                    in the column samples. 

  --queue, -q [STR]                 If using a cluster, this is the name of the queue
                                    to which the jobs should be sent. Default
                                    is 'bio'.  

  --cores [INT]                     Number of cores to use to run the pipeline. If
                                    running in a cluster, the default is 300. If 
                                    running locally, the default is 4.
                                  
  --local                           If this flag is present, the pipeline is run 
                                    locally instead of in a cluster. The default 
                                    is to run in a cluster ('bio' queue)


Output (automatically generated):
  <output_dir>                      Contains dir contains the results of every 
                                    step of the pipeline.

  <output_dir> /log/                Contains the log files for every step of 
                                    the pipeline

  <output_dir> /log/drmaa           Contains the .out and .err files of every 
                                    job sent to the grid/cluster.

  <output_dir> /log/results         Contains the log files and parameters that 
                                    the pipeline used for the current run


Parameters:
-h, --help                          Print this help document.

-sh, --snakemake-help               Print help document for snakemake.

-y, --skip-confirmation             Skip confirmation (-y forces 'Yes' on all 
                                    prompts).
  
Other snakemake parameters          Any other parameters will be passed to 
                                    snakemake. Read snakemake help (-sh) to 
                                    see the options.
"""
        