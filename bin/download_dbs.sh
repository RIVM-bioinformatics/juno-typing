#!/bin/bash

CURRENT_DIR=$(pwd)
DB_DIR=$(realpath $1)
UPDATE_DBS="$2"
LOG_SOFTWARE_VERSIONS="$3"

set -eu

###############################################################################
#####      Download and installing databases and necessary software       #####
###############################################################################
rm -f "$LOG_SOFTWARE_VERSIONS"
mkdir -p $DB_DIR

########################### List software/db to download ######################
KMERFINDER_DB="${DB_DIR%/}/kmerfinder_db"
KMERFINDER_DB_DOWNLOAD="${CURRENT_DIR%/}/bin/kmerfinder_db"
KMERFINDER="${CURRENT_DIR%/}/bin/kmerfinder"
CGEMLST="${CURRENT_DIR%/}/bin/cge-mlst"
MLST7_DB="${DB_DIR%/}/mlst7_db"
KMERFINDER_DB_V="20210228" 

#################### KmerFinder and the database download #####################
# Check available versions here: http://www.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/
# For some reason the 'latest' version does not work
if [ ! -f "${KMERFINDER_DB%/}/config" ] || [ "$UPDATE_DBS" == "TRUE" ]; then
    echo -e "Downloading KmerFinder database"
    rm -rf "${KMERFINDER_DB_DOWNLOAD}"
    git clone -b 'master' --single-branch --depth=1 https://bitbucket.org/genomicepidemiology/kmerfinder_db.git "$KMERFINDER_DB_DOWNLOAD"
    mkdir -p "$KMERFINDER_DB"
    bash "${KMERFINDER_DB_DOWNLOAD%/}/INSTALL.sh" "$KMERFINDER_DB/" bacteria $KMERFINDER_DB_V 
fi

if [ ! -f "${KMERFINDER%}/kmerfinder.py" ] ; then
    echo -e "Downloading KmerFinder software"
    rm -rf "${KMERFINDER}"
    git clone -b '3.0.2' --single-branch --depth=1 https://bitbucket.org/genomicepidemiology/kmerfinder.git "$KMERFINDER"
fi

####################### CGE-MLST and database download ########################

if [ ! -f "${CGEMLST%/}/mlst.py" ]; then
    echo -e "Downloading CGE-MLST software"
    rm -rf "${CGEMLST}" 
    git clone -b '2.0.4' --single-branch --depth 1 https://bitbucket.org/genomicepidemiology/mlst.git "$CGEMLST"
    echo -e "\nCGE-MLST software has been installed\n"
fi

if [ ! -f "${MLST7_DB%/}/senterica/senterica.length.b" ] || [ "$UPDATE_DBS" == "TRUE" ]; then
    echo -e "Downloading CGE-MLST database"
    rm -rf "${MLST7_DB}"
    git clone -b 'master' --single-branch --depth=1 https://bitbucket.org/genomicepidemiology/mlst_db.git "$MLST7_DB"
    CURRENT_DIR=$(pwd)
    cd $MLST7_DB
    python INSTALL.py kma_index
    cd $CURRENT_DIR
    echo -e "\nDatabase for 7-locus MLST has been downloaded\n"
fi

################# Log versions for each software/db downloaded ################
cd $KMERFINDER
KMERFINDER_V="$(git log -n 1 --pretty=format:"%H")"
cd $CGEMLST
CGEMLST_V="$(git log -n 1 --pretty=format:"%H")"
cd $MLST7_DB
MLST7_DB_V="$(git log -n 1 --pretty=format:"%H")"

cd $CURRENT_DIR

echo "
Extra software versions:

- KmerFinder: ${KMERFINDER_V}
- CGE-MLST: ${CGEMLST_V}

Database version:

- KmerFinder: ${KMERFINDER_DB_V}
- CGE-MLST: ${MLST7_DB}
" > "$LOG_SOFTWARE_VERSIONS"