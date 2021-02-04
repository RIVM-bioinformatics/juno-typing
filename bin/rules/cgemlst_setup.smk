#############################################################################
#####                               MLST7                               #####
#############################################################################

rule cgemlst_setup:
    input:
        config["sample_sheet"]
    output:
        OUT + "/cgemlst_setup_done.txt"
    conda:
        "../../envs/cgemlst.yaml"
    benchmark:
        OUT + "/log/benchmark/cgemlst_setup.txt"
    threads: 1
    log:
        OUT + "/log/cgemlst_setup.log"
    shell:
        """
CURRENT_DIR=$(pwd)
LOG_FILE=$(realpath {log})

cd ${{CURRENT_DIR}}/bin/

if [ ! -d kma ]; then
    git clone -b '1.3.9' --single-branch --depth 1 https://bitbucket.org/genomicepidemiology/kma.git 2> ${{LOG_FILE}}
    cd kma && make 2> ${{CURRENT_DIR}}/{log}
    rm -rf "kma/.git" 2> ${{CURRENT_DIR}}/{log}
    echo "\nkma has been installed\n"
fi

cd ${{CURRENT_DIR}}/bin

if [ ! -d cge-mlst ]; then
    git clone -b '2.0.4' --single-branch --depth 1 https://bitbucket.org/genomicepidemiology/mlst.git cge-mlst 2> ${{LOG_FILE}}
    rm -rf "mlst_db/.git" 2> ${{CURRENT_DIR}}/{log}
    echo "\ncge-mlst has been installed\n"
fi

cd ${{CURRENT_DIR}} 

touch {output}
        """