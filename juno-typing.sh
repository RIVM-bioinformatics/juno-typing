#!/bin/bash
###############################################################################################################################################
### juno-typing pipeline                                                                                                                   ### 
### Authors: Alejandra Hernandez-Segura, Maaike van der Beld                                                                                ### 
### Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)                                                                      ### 
### Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium Surveillance (IDS), Bacteriologie (BPD)                                ### 
### Date: 03-02-2021                                                                                                                        ### 
###                                                                                                                                         ### 
### Documentation: https://github.com/AleSR13/Juno-typing                                                                               ### 
###                                                                                                                                         ### 
###                                                                                                                                         ### 
### Snakemake rules (in order of execution):                                                                                                ### 
###     1 cgemlst_builddb   # Download and build database for MLST using the CGE-MLST software:                                             ###
###                         https://bitbucket.org/genomicepidemiology/mlst/src/master/.                                                     ### 
###     2 cgemlst_fastq     # Do the 7-locus MLST taking filtered fastq files as input. Either the cgemlst_fastq or the cgemlst_fasta       ###
###                         are run, according to the user's input.                                                                         ### 
###     2 cgemlst_fasta     # Do the 7-locus MLST taking filtered fasta files (assemblies) as input. Either the cgemlst_fastq or the        ###
###                         cgemlst_fastq are run, according to the user's input.                                                           ###                                                     ###
###############################################################################################################################################


### conda environment
PATH_MASTER_YAML="envs/master_env.yaml"
MASTER_NAME=$(head -n 1 ${PATH_MASTER_YAML} | cut -f2 -d ' ') # Extract Conda environment name as specified in yaml file

## Generate ID and get host info
source bin/include/functions.sh
UNIQUE_ID=$(bin/include/generate_id.sh)
SET_HOSTNAME=$(bin/include/gethostname.sh)

### Default values for CLI parameters
INPUT_DIR="raw_data"
OUTPUT_DIR="out"
SPECIES_ALL="NotProvided"
METADATA="NotProvided"
DB_DIR="db"
SNAKEMAKE_HELP="FALSE"
HELP="FALSE"
MAKE_SAMPLE_SHEET="FALSE"
SHEET_SUCCESS="FALSE"

### Parse the commandline arguments, if they are not part of the pipeline, they get send to Snakemake
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -i|--input)
        INPUT_DIR="${2%/}"
        shift # Next
        shift # Next
        ;;
        -o|--output)
        OUTPUT_DIR="${2%/}"
        shift # Next
        shift # Next
        ;;
        --species)
        SPECIES_ALL="${2} ${3}"
        shift
        shift
        shift
        ;;
        --metadata)
        METADATA="$2"
        shift
        shift
        ;;
        -h|--help)
        HELP="TRUE"
        shift # Next
        ;;
        -sh|--snakemake-help)
        SNAKEMAKE_HELP="TRUE"
        shift # Next
        ;;
        --make-sample-sheet)
        MAKE_SAMPLE_SHEET="TRUE"
        shift # Next
        ;;
        *) # Any other option
        POSITIONAL+=("$1") # save in array
        shift # Next
        ;;
    esac
done
set -- "${POSITIONAL[@]:-}" # Restores the positional arguments (i.e. without the case arguments above) which then can be called via `$@` or `$[0-9]` etc. These parameters are send to Snakemake.


### Print Juno help message
if [ "${HELP:-}" == "TRUE" ]; then
    line
    cat bin/include/help.txt
    exit 0
fi

###############################################################################################################
##### Installation block                                                                                  #####
###############################################################################################################

### Pre-flight check: Assess availability of required files, conda and master environment
if [ ! -e "${PATH_MASTER_YAML}" ]; then # If this yaml file does not exist, give error.
    line
    spacer
    echo -e "ERROR: Missing file \"${PATH_MASTER_YAML}\""
    exit 1
fi

if [[ $PATH != *${MASTER_NAME}* ]]; then # If the master environment is not in your path (i.e. it is not currently active), do...
    line
    spacer
    conda env update -f envs/mamba.yaml -q -v
    source activate mamba
    source activate ${MASTER_NAME}
    if [ ! $? -eq 0 ]; then
    	set +ue # Turn bash strict mode off because that breaks conda
    	if [ "${SKIP_CONFIRMATION}" = "TRUE" ]; then
       		echo -e "\tInstalling master environment..." 
       		mamba env update -f ${PATH_MASTER_YAML} 
       		echo -e "DONE"
    	else
       		while read -r -p "The master environment hasn't been installed yet, do you want to install this environment now? [y/n] " envanswer
       		do
            		envanswer=${envanswer,,}
            		if [[ "${envanswer}" =~ ^(yes|y)$ ]]; then
                		echo -e "\tInstalling master environment..." 
				mamba env update -f ${PATH_MASTER_YAML}
                		echo -e "DONE"
                		break
            		elif [[ "${envanswer}" =~ ^(no|n)$ ]]; then
                		echo -e "The master environment is a requirement. Exiting because Juno cannot continue without this environment"
                		exit 1
            		else
                		echo -e "Please answer with 'yes' or 'no'"
            		fi
        	done
    	fi
    fi
    source activate "${MASTER_NAME}"
    set -ue # Turn bash strict mode on again
fi


### Print Snakemake help
if [ "${SNAKEMAKE_HELP:-}" == "TRUE" ]; then
    line
    snakemake --help
    exit 0
fi


### Check argument validity
if [ ! -d "${INPUT_DIR}" ]; then
    minispacer
    echo -e "The input directory specified (${INPUT_DIR}) does not exist"
    echo -e "Please specify an existing input directory"
    minispacer
    exit 1
fi

if [ ${METADATA} == "NotProvided" ] && [ ${SPECIES_ALL} == "NotProvided" ]; then
    minispacer
    echo "ERROR! You need to provide either a --metadata file or the --species argument."
    minispacer
    exit 1
fi

if [ ${METADATA} != "NotProvided" ]; then
    if [ ! -f ${METADATA} ] || [[ ! ${METADATA} =~ ".csv"$ ]]; then
        minispacer
        echo "ERROR! The provided metadata file (${METADATA}) does not exist or does not have the .csv extension. Please provide an existing metadata .csv file."
        minispacer
        exit 1
    fi
fi

### Generate sample sheet
if [  `ls -A "${INPUT_DIR}" | grep 'R[0-9]\{1\}.*\.f[ast]\{0,3\}q\.\?[gz]\{0,2\}$' | wc -l` -gt 0 ]; then
    minispacer
    echo -e "Files in input directory (${INPUT_DIR}) are present"
    echo -e "Generating sample sheet..."
    if [ ${METADATA} == "NotProvided" ];then 
        python bin/generate_sample_sheet.py "${INPUT_DIR}" > sample_sheet.yaml
    else
        python bin/generate_sample_sheet.py "${INPUT_DIR}" --metadata ${METADATA} > sample_sheet.yaml
    fi
    if [ $(wc -l sample_sheet.yaml | awk '{ print $1 }') -gt 2 ]; then
        SHEET_SUCCESS="TRUE"
    fi
else
    minispacer
    echo -e "The input directory you specified (${INPUT_DIR}) exists but is empty or does not contain the expected input files...\nPlease specify a directory with input-data."
    exit 0
fi

### Checker for succesfull creation of sample_sheet
if [ "${SHEET_SUCCESS}" == "TRUE" ]; then
    echo -e "Succesfully generated the sample sheet"
    echo -e "ready_for_start"
else
    echo -e "Couldn't find files in the input directory that ended up being in a .FASTQ, .FQ or .GZ format"
    echo -e "Please inspect the input directory (${INPUT_DIR}) and make sure the files are in one of the formats listed below"
    echo -e ".fastq.gz (Zipped Fastq)"
    echo -e ".fq.gz (Zipped Fq)"
    echo -e ".fastq (Unzipped Fastq)"
    echo -e ".fq (unzipped Fq)"
    exit 1
fi


### Actual snakemake command with checkers for required files. N.B. here the UNIQUE_ID and SET_HOSTNAME variables are set!
if [ -e sample_sheet.yaml ]; then
    echo -e "Starting snakemake"
    set +ue #turn off bash strict mode because snakemake and conda can't work with it properly
    echo -e "pipeline_run:\n    identifier: ${UNIQUE_ID}" > config/variables.yaml
    echo -e "Server_host:\n    hostname: http://${SET_HOSTNAME}" >> config/variables.yaml
    eval $(parse_yaml config/variables.yaml "config_")
    snakemake --config out="$OUTPUT_DIR" species="${SPECIES_ALL}" \
        --profile config --drmaa " -q bio -n {threads} -R \"span[hosts=1]\"" --drmaa-log-dir ${OUTPUT_DIR}/log/drmaa ${@}
    set -ue #turn bash strict mode back on
else
    echo -e "Sample_sheet.yaml could not be found"
    echo -e "This also means that the pipeline was unable to generate a new sample sheet for you"
    echo -e "Please inspect the input directory (${INPUT_DIR}) and make sure the right files are present"
    exit 1
fi

exit 0 