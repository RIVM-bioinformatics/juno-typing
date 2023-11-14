#! /usr/bin/env python3
"""
Juno-typing pipeline
Authors: Alejandra Hernandez-Segura
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium
            Surveillance (IDS), Bacteriologie (BPD)     
Date: 18-05-2021   

Documentation: https://rivm-bioinformatics.github.io/ids_bacteriology_man/juno-typing.html 

"""

# Dependencies
import argparse
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import yaml
from juno_library import Pipeline

# Own scripts
import bin.download_dbs
from version import __package_name__, __version__


def main() -> None:
    juno_typing = JunoTyping()
    juno_typing.run()


@dataclass
class JunoTyping(Pipeline):
    pipeline_name: str = __package_name__
    pipeline_version: str = __version__
    input_type: str = "both"

    def _add_args_to_parser(self) -> None:
        super()._add_args_to_parser()

        self.parser.description = "Juno-typing pipeline. Automated pipeline for bacterial subtyping (7-locus MLST and serotyping)."

        self.add_argument(
            "-m",
            "--metadata",
            type=Path,
            default=None,
            required=False,
            metavar="FILE",
            help="Relative or absolute path to the metadata csv file. If "
            "provided, it must contain at least one column named 'sample' "
            "with the name of the sample (same than file name but removing "
            "the suffix _R1.fastq.gz), a column called "
            "'genus' and a column called 'species'. The genus and species "
            "provided will be used to choose the serotyper and the MLST schema(s)."
            "If a metadata file is provided, it will overwrite the --species "
            "argument for the samples present in the metadata file.",
        )
        self.add_argument(
            "-s",
            "--species",
            type=lambda s: s.strip().lower(),
            nargs=2,
            default=[None, None],
            required=False,
            metavar=("GENUS", "SPECIES"),
            help="Species name (any species in the metadata file will overwrite"
            " this argument). It should be given as two words (e.g. --species "
            "Salmonella enterica)",
        )
        self.add_argument(
            "-d",
            "--db_dir",
            type=Path,
            required=False,
            metavar="DIR",
            default="/mnt/db/juno/typing_db",
            help="Relative or absolute path to the directory that contains the databases for all the tools used in this pipeline or where they should be downloaded. Default is: /mnt/db/juno/typing_db",
        )
        self.add_argument(
            "--serotypefinder_mincov",
            type=float,
            metavar="NUM",
            default=0.6,
            help="Minimum coverage to be used for the SerotypeFinder (E. coli serotyping) tool. It accepts values from 0-1. Default is 0.6.",
        )
        self.add_argument(
            "--serotypefinder_identity",
            type=float,
            metavar="NUM",
            default=0.85,
            help="Identity threshold to be used for the SerotypeFinder (E. coli serotyping) tool. It accepts values from 0-1. Default is 0.85",
        )
        self.add_argument(
            "--seroba_mincov",
            type=int,
            metavar="INT",
            default=20,
            help="Minimum coverage to be used for the Seroba (S. pneumoniae serotyping) tool. It accepts values from 0-100. Default is 20",
        )
        self.add_argument(
            "--seroba_kmersize",
            type=int,
            metavar="INT",
            default=71,
            help="Minimum coverage to be used for the SerotypeFinder (E. coli serotyping) tool. It accepts values from 0-100. Default is 20",
        )
        self.add_argument(
            "--bordetella_vaccine_antigen_scheme_name",
            type=str,
            metavar="DIR",
            default="bordetella",
            help="Name for the directory containing the Bordetella vaccine antigen MLST scheme in --db_dir. Should contain a BLAST db with base name bordetella.fa",
        )
        self.add_argument(
            "--update",
            action="store_true",
            help="Force database update even if they are present.",
        )

    def _parse_args(self) -> argparse.Namespace:
        # Remove this if containers can be used with juno-typing
        if "--no-containers" not in self.argv:
            self.argv.append("--no-containers")

        args = super()._parse_args()
        self.db_dir: Path = args.db_dir.resolve()

        self.genus: Optional[str]
        self.species: Optional[str]
        self.genus, self.species = args.species
        self.metadata_file: Path = args.metadata
        self.serotypefinder_mincov: float = args.serotypefinder_mincov
        self.serotypefinder_identity: float = args.serotypefinder_identity
        self.seroba_mincov: int = args.seroba_mincov
        self.seroba_kmersize: int = args.seroba_kmersize
        self.bordetella_vaccine_antigen_scheme: str = (
            args.bordetella_vaccine_antigen_scheme_name
        )
        self.update_dbs: bool = args.update
        return args

    def setup(self) -> None:
        super().setup()
        self.update_sample_dict_with_metadata()

        if self.snakemake_args["use_singularity"]:
            self.snakemake_args["singularity_args"] = " ".join(
                [
                    self.snakemake_args["singularity_args"],
                    f"--bind {self.db_dir}:{self.db_dir}",
                ]
            )
        self.user_parameters = {
            "input_dir": str(self.input_dir),
            "out": str(self.output_dir),
            "db_dir": str(self.db_dir),
            "mlst7_db": str(self.db_dir.joinpath("mlst7_db")),
            "seroba_db": str(self.db_dir.joinpath("seroba_db")),
            "serotypefinder_db": str(self.db_dir.joinpath("serotypefinder_db")),
            "serotypefinder": {
                "min_cov": self.serotypefinder_mincov,
                "identity_thresh": self.serotypefinder_identity,
            },
            "seroba": {
                "min_cov": self.seroba_mincov,
                "kmer_size": self.seroba_kmersize,
            },
            "bordetella_vaccine_antigen_scheme": str(
                self.bordetella_vaccine_antigen_scheme
            ),
            "bordetella_vaccine_antigen_blastdb": str(
                self.db_dir.joinpath(
                    self.bordetella_vaccine_antigen_scheme,
                    "bordetella.fa",
                )
            ),
        }

        with open(
            Path(__file__).parent.joinpath("config/pipeline_parameters.yaml")
        ) as f:
            parameters_dict = yaml.safe_load(f)
        self.snakemake_config.update(parameters_dict)

    def update_sample_dict_with_metadata(self) -> None:
        self.get_metadata_from_csv_file(
            filepath=self.metadata_file,
            expected_colnames=["sample", "genus", "species"],
        )
        # Add metadata
        for sample in self.sample_dict:
            if self.genus is not None and self.species is not None:
                self.sample_dict[sample]["genus"] = self.genus
                self.sample_dict[sample]["species"] = self.species
            else:
                try:
                    self.sample_dict[sample].update(self.juno_metadata[sample])
                except (KeyError, TypeError, AttributeError):
                    raise ValueError(
                        f"One of your samples is not in the metadata file "
                        f"({self.metadata_file}). Please ensure that all "
                        "samples are present in the metadata file or provide "
                        "a --species argument."
                    )
                self.sample_dict[sample]["genus"] = (
                    self.sample_dict[sample]["genus"].strip().lower()
                )
                self.sample_dict[sample]["species"] = (
                    self.sample_dict[sample]["species"].strip().lower()
                )
        # Update self.sample_dict
        with open("files/dictionary_correct_species.yaml") as translation_yaml:
            self.mlst7_species_translation_tbl = yaml.safe_load(translation_yaml)
            for sample in self.sample_dict:
                genus = self.sample_dict[sample]["genus"].lower()
                species = self.sample_dict[sample]["species"].lower()
                if key := self.mlst7_species_translation_tbl.get(f"{genus}_{species}"):
                    self.sample_dict[sample]["species-mlst7"] = key
                else:
                    self.sample_dict[sample][
                        "species-mlst7"
                    ] = self.mlst7_species_translation_tbl.get(genus)

    def run(self) -> None:
        self.setup()
        if not self.dryrun or self.unlock:
            self.path_to_audit.mkdir(parents=True, exist_ok=True)
            downloads_juno_typing = bin.download_dbs.DownloadsJunoTyping(
                self.db_dir,
                update_dbs=self.update_dbs,
                cge_mlst_asked_version="2.0.4",
                mlst7_db_asked_version="master",
                serotypefinder_db_asked_version="master",
                seroba_db_asked_version="master",
            )
            self.downloads_versions = downloads_juno_typing.downloaded_versions
            with open(
                self.path_to_audit.joinpath("database_versions.yaml"), "w"
            ) as file_:
                yaml.dump(self.downloads_versions, file_, default_flow_style=False)

        if not self.dryrun or self.unlock:
            subprocess.run(
                [
                    "find",
                    self.output_dir,
                    "-type",
                    "f",
                    "-empty",
                    "-exec",
                    "rm",
                    "{}",
                    ";",
                ]
            )
            subprocess.run(
                [
                    "find",
                    self.output_dir,
                    "-type",
                    "d",
                    "-empty",
                    "-exec",
                    "rm",
                    "-rf",
                    "{}",
                    ";",
                ]
            )
        super().run()


if __name__ == "__main__":
    main()
