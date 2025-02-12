import yaml


def return_ssonnei_samples(mock_variable):
    with open(
        checkpoints.parse_serotype_multireport.get(**mock_variable).output[0]
    ) as file:
        SPECIES_PER_SAMPLE = yaml.safe_load(file)
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


checkpoint parse_serotype_multireport:
    input:
        OUT + "/serotype/serotyper_multireport.csv",
    output:
        OUT + "/serotype/serotyper_multireport_outputs.yaml",
    run:
        import pandas as pd
        import yaml

        df = pd.read_csv(input[0])

        # check if prediction and ipaB are in the columns
        # if not, read in the file with .csv stripped and 1.csv added
        if ("prediction" not in df.columns) and ("ipaB" not in df.columns):
            df = pd.read_csv(input[0].replace(".csv", "1.csv"))

            # get the species/serotype prediction per sample
        df_species = df[["Samplename", "prediction"]].copy()

        # write to yaml
        with open(output[0], "w") as file:
            yaml.dump(
                df_species.set_index("Samplename").to_dict()["prediction"],
                file,
                default_flow_style=False,
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
mykrobe predict \
    --format json \
    --sample {params.sample} \
    --species {params.species} \
    --output {output} \
    --threads {threads} \
    -i {input.r1} {input.r2}
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
if [[ {input} == files/empty_lineage_result.csv ]]
then
    cp {input} {output}
else
    PREFIX=$(dirname {output})
    parse_mykrobe_predict --jsons {input} --prefix $PREFIX/multireport
    mv $PREFIX/multireport_predictResults.tsv {output}
fi
        """
