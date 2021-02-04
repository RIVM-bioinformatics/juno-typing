#############################################################################
#####                               MLST7                               #####
#############################################################################

rule cgemlst_builddb:
    input:
        sample_sheet = config["sample_sheet"],
        setup = OUT + "/cgemlst_setup_done.txt"
    output:
        OUT + "/checker_cgemlst_builddb.txt"
    conda:
        "../../envs/cgemlst.yaml"
    benchmark:
        OUT + "/log/benchmark/cgemlst_builddb.txt"
    threads: 1
    log:
        OUT + "/log/cgemlst_builddb.log"
    shell:
        """
CURRENT_DIR=$(pwd)
OUT_DIR=`dirname {output}`
OUT_DIR=`realpath $OUT_DIR`

cd ${{CURRENT_DIR}}/bin/

cd $OUT_DIR

if [ ! -d mlst_db ]; then
    git clone -b 'master' --single-branch --depth=1 https://bitbucket.org/genomicepidemiology/mlst_db.git 2> ${{CURRENT_DIR}}/{log}
    rm -rf mlst_db/.git 2> ${{CURRENT_DIR}}/{log}
fi

cd mlst_db
MLST_DB=$(pwd)

python INSTALL.py ${{CURRENT_DIR}}/bin/kma/kma_index 2> ${{CURRENT_DIR}}/{log}

echo "Downloaded CGE-mlst database" > ${{CURRENT_DIR}}/{output}
        """
