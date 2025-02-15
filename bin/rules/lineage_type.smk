import yaml
import pandas as pd


def return_ssonnei_samples(mock_wildcard):
    serotype_multireport_filepath = checkpoints.serotype_multireports.get(
        **mock_wildcard
    ).output[0]

    _iteration = 1

    df = pd.read_csv(serotype_multireport_filepath)

    # check if prediction and ipaB are in the columns (shigatyper multireport)
    # if not, read in the file with .csv stripped and 1.csv added
    # if still not, continue incrementing until a file is not found. This should only happen if different species are typed in a single analysis.
    while ("prediction" not in df.columns) and ("ipaB" not in df.columns):
        # check if other serotype multireport exists
        if not os.path.exists(
            serotype_multireport_filepath.replace(".csv", f"{_iteration}.csv")
        ):
            # this means all serotype multireports have been checked and none are in shigatyper format
            # in this case, just return empty placeholder file
            return ["files/empty_lineage_result.csv"]
        else:
            df = pd.read_csv(
                serotype_multireport_filepath.replace(".csv", f"{_iteration}.csv")
            )
            _iteration += 1

    # get a list of Samplename for the samples that have "Shigella sonnei" in the prediction
    SSONNEI_SAMPLES = df.loc[
        df["prediction"].str.contains("Shigella sonnei"), "Samplename"
    ].tolist()
    if len(SSONNEI_SAMPLES) == 0:
        # if the multireport is shigatyper format, but no sonnei detected: also return empty placeholder file
        return ["files/empty_lineage_result.csv"]
    else:
        return expand(
            OUT + "/lineage_typing/mykrobe_ssonnei/{sample}.json",
            sample=SSONNEI_SAMPLES,
        )

    # only keep entries that have "Shigella sonnei" as the species
    SSONNEI_SAMPLES = {
        sample: species
        for sample, species in SPECIES_PER_SAMPLE.items()
        if "Shigella sonnei" in species
    }
    if len(SSONNEI_SAMPLES) == 0:
        output_ssonei_mykrobe = ["files/empty_lineage_result.csv"]
    else:
        output_ssonei_mykrobe = expand(
            OUT + "/lineage_typing/mykrobe_ssonnei/{sample}.json",
            sample=SSONNEI_SAMPLES,
        )
    return output_ssonei_mykrobe


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
