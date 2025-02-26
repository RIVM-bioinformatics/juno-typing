import yaml
import pandas as pd
from pathlib import Path


def return_ssonnei_samples(mock_wildcard):
    """
    Read multireport, select samples on criteria, and return corresponding files.

    Parameters
    ----------
    mock_wildcard : dict
        Wildcard dictionary. If a multireport is the input, this is a mock wildcard.

    Returns
    -------
    list
        List of expected input files.

    Notes
    -----
    This function can serve as an example of how to use an input function after a snakemake
    checkpoint based on a multireport.

    """
    empty_placeholder = "files/empty_lineage_result.csv"
    output_file = OUT + "/lineage_typing/mykrobe_ssonnei/{sample}.json"

    # Step 1. Find and read multireport
    base_path = Path(checkpoints.serotype_multireports.get(**mock_wildcard).output[0])

    # start with empty df
    df = pd.DataFrame()

    # check potential multireport files (multireport.csv, multireport1.csv, multireport2.csv, etc)
    for i in range(10):
        file_path = (
            base_path if i == 0 else base_path.with_name(f"{base_path.stem}{i}.csv")
        )

        if file_path.exists():
            df = pd.read_csv(file_path)

            # Check if multireport is based on Shigatyper
            if ("prediction" in df.columns) and ("ipaB" in df.columns):
                break
            else:
                # reset df to empty if the file is not shigatyper format
                df = pd.DataFrame()

    # Step 2: Extract relevant samples based on condition(s)
    if not df.empty:
        ssonnei_samples = df.loc[
            df["prediction"].str.contains("Shigella sonnei", na=False), "Samplename"
        ].tolist()
    else:
        ssonnei_samples = []

    # Step 3: Return the relevant samples
    if len(ssonnei_samples) == 0:
        return [empty_placeholder]
    else:
        return expand(
            output_file,
            sample=ssonnei_samples,
        )


rule mykrobe_ssonnei:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
    output:
        OUT + "/lineage_typing/mykrobe_ssonnei/{sample}.json",
    conda:
        "../../envs/mykrobe.yaml"
    container:
        "docker://quay.io/biocontainers/mykrobe:0.13.0--py311h4f141b6_1"
    threads: config["threads"]["mykrobe"]
    resources:
        mem_gb=config["mem_gb"]["mykrobe"],
    params:
        species="sonnei",
        sample="{sample}",
    shell:
        """
DIR=$(dirname {output})
mkdir -p $DIR/tmp_{params.sample}
mykrobe predict \
    --format json \
    --sample {params.sample} \
    --species {params.species} \
    --output {output} \
    --threads {threads} \
    --skeleton_dir $DIR/tmp_{params.sample} \
    -i {input.r1} {input.r2}
rm -rf $DIR/tmp_{params.sample}
        """


rule sonneityping:
    input:
        return_ssonnei_samples,
    output:
        OUT + "/lineage_typing/sonneityping/multireport.tsv",
    conda:
        "../../envs/mykrobe.yaml"
    container:
        ""
    threads: config["threads"]["mykrobe"]
    resources:
        mem_gb=config["mem_gb"]["mykrobe"],
    shell:
        """
if [ {input} == files/empty_lineage_result.csv ]
then
    cp {input} {output}
else
    PREFIX=$(dirname {output})
    parse_mykrobe_predict --jsons {input} --prefix $PREFIX/multireport
    mv $PREFIX/multireport_predictResults.tsv {output}
fi
        """
