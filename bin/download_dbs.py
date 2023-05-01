import argparse
import pathlib
import subprocess
from base_juno_pipeline import helper_functions


class DownloadsJunoTyping(helper_functions.GitHelpers):
    """Class that performs all necessary software and database downloads for
    the Juno typing pipeline if necessary"""

    def __init__(
        self,
        db_dir,
        update_dbs=False,
        cge_mlst_asked_version="2.0.4",
        characterize_neisseria_capsule_asked_version="master",
        mlst7_db_asked_version="master",
        serotypefinder_db_asked_version="master",
        seroba_db_asked_version="master",
        seroba_kmersize=71,
    ):
        self.db_dir = pathlib.Path(db_dir)
        self.bin_dir = pathlib.Path(__file__).parent.absolute()
        self.update_dbs = update_dbs
        self.seroba_kmersize = seroba_kmersize
        self.downloaded_versions = self.get_downloads_juno_typing(
            cge_mlst_asked_version=cge_mlst_asked_version,
            characterize_neisseria_capsule_asked_version=characterize_neisseria_capsule_asked_version,
            mlst7_db_asked_version=mlst7_db_asked_version,
            serotypefinder_db_asked_version=serotypefinder_db_asked_version,
            seroba_db_asked_version=seroba_db_asked_version,
        )

    def download_software_kmerfinder(self, version):
        """Function to download kmerfinder if it is not present"""
        kmerfinder_software_dir = self.bin_dir.joinpath("kmerfinder")
        if not kmerfinder_software_dir.joinpath("kmerfinder.py").is_file():
            print("\x1b[0;33m Downloading kmerfinder software...\n\033[0;0m")
            self.download_git_repo(
                version,
                "https://bitbucket.org/genomicepidemiology/kmerfinder.git",
                kmerfinder_software_dir,
            )
        return version

    def download_software_characterize_neisseria_capsule(self, version):
        """Function to download characterize_neisseria_capsule if it is not present"""
        characterize_neisseria_capsule_software_dir = self.bin_dir.joinpath(
            "characterize_neisseria_capsule"
        )
        if not characterize_neisseria_capsule_software_dir.joinpath(
            "characterize_neisseria_capsule.py"
        ).is_file():
            print(
                "\x1b[0;33m Downloading characterize_neisseria_capsule software...\n\033[0;0m"
            )
            self.download_git_repo(
                version,
                "https://github.com/ntopaz/characterize_neisseria_capsule.git",
                characterize_neisseria_capsule_software_dir,
            )
        return version  

    def download_software_mlst7(self, version):
        """Function to download MLST (CGE) if it is not present"""
        mlst7_software_dir = self.bin_dir.joinpath("cge-mlst")
        if not mlst7_software_dir.joinpath("mlst.py").is_file():
            print("\x1b[0;33m Downloading MLST (CGE) software...\n\033[0;0m")
            self.download_git_repo(
                version,
                "https://bitbucket.org/genomicepidemiology/mlst.git",
                mlst7_software_dir,
            )
        return version

    def download_db_kmerfinder(self, version):
        """Function to download kmerfinder database if it is not present"""
        kmerfinder_db_dir = self.db_dir.joinpath("kmerfinder_db")
        if not kmerfinder_db_dir.joinpath("config").exists():
            print("\x1b[0;33m Downloading KmerFinder database...\n\033[0;0m")
            self.download_git_repo(
                "master",
                "https://bitbucket.org/genomicepidemiology/kmerfinder_db.git",
                kmerfinder_db_dir,
            )
            try:
                build = subprocess.run(
                    [
                        "bash",
                        "INSTALL.sh",
                        str(kmerfinder_db_dir.absolute()),
                        "bacteria",
                        version,
                    ],
                    cwd=str(kmerfinder_db_dir),
                    check=True,
                    timeout=3000,
                )
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as err:
                build.kill()
                raise
        return version

    def copy_neisseria_db(self):
        """Function to copy the neisseria db from mnt/db/juno to the bin folder of the neisseria tool.
        It is necessary for the tool to run the database from this location."""
        characterize_neisseria_capsule_db_dir = self.bin_dir.joinpath(
            "characterize_neisseria_capsule"
        )
        if not characterize_neisseria_capsule_db_dir.joinpath(
            "neisseria_capsule_DB"
        ).exists():
            try:
                print(f"Copying neisseria db from mnt/db/juno to: {characterize_neisseria_capsule_db_dir}")
                build = subprocess.run(
                    ["cp", "-R", "/mnt/db/juno/neisseria_capsule_DB", f"{characterize_neisseria_capsule_db_dir}"],
                    cwd=str(characterize_neisseria_capsule_db_dir),
                    check=True,
                    timeout=3000,
                )
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as err:
                build.kill()
                raise
        else:
            return print("Neisseria db is available, continue analysis")

    def download_db_mlst7(self, version):
        """Function to download the MLST (CGE) database if it is not present"""
        mlst7_db_dir = self.db_dir.joinpath("mlst7_db")
        if not mlst7_db_dir.joinpath("senterica", "senterica.length.b").is_file():
            print("\x1b[0;33m Downloading 7-locus MLST (CGE) database...\n\033[0;0m")
            self.download_git_repo(
                version,
                "https://bitbucket.org/genomicepidemiology/mlst_db.git",
                mlst7_db_dir,
            )
            try:
                build = subprocess.run(
                    ["python", "INSTALL.py", "kma_index"],
                    check=True,
                    cwd=str(mlst7_db_dir),
                    timeout=800,
                )
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as err:
                build.kill()
                raise
        version = self.get_commit_git(mlst7_db_dir)
        return version

    def download_db_serotypefinder(self, version):
        """Function to download the SerotypeFinder database if it is not present"""
        serotypefinder_db_dir = self.db_dir.joinpath("serotypefinder_db")
        if not serotypefinder_db_dir.joinpath("H_type.seq.b").is_file():
            print("\x1b[0;33m Downloading SerotypeFinder database...\n\033[0;0m")
            self.download_git_repo(
                version,
                "https://bitbucket.org/genomicepidemiology/serotypefinder_db.git",
                serotypefinder_db_dir,
            )
            try:
                build = subprocess.run(
                    ["python", str("INSTALL.py"), "kma_index"],
                    cwd=str(serotypefinder_db_dir),
                    check=True,
                    timeout=800,
                )
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as err:
                build.kill()
                raise
        version = self.get_commit_git(serotypefinder_db_dir)
        return version

    def download_db_seroba(self, version, kmersize=71):
        """Function to download the Seroba database if it is not present"""
        seroba_db_dir = self.db_dir.joinpath("seroba_db")
        if not seroba_db_dir.joinpath("database", "cdhit_cluster").is_file():
            print("\x1b[0;33m Downloading Seroba database...\n\033[0;0m")
            self.download_git_repo(
                version, "https://github.com/sanger-pathogens/seroba.git", seroba_db_dir
            )
            try:
                kmersize = str(kmersize)
                rm_dir = subprocess.run(
                    [
                        "rm",
                        "-rf",
                        str(seroba_db_dir.joinpath("scripts")),
                        "&&",
                        "rm",
                        "-rf",
                        str(seroba_db_dir.joinpath("seroba")),
                    ],
                    check=True,
                    timeout=60,
                )
                build = subprocess.run(
                    ["seroba", "createDBs", "database", kmersize],
                    cwd=str(seroba_db_dir),
                    check=True,
                    timeout=800,
                )
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as err:
                build.kill()
                raise
        version = self.get_commit_git(seroba_db_dir)
        return version

    def get_downloads_juno_typing(
        self,
        cge_mlst_asked_version,
        characterize_neisseria_capsule_asked_version,
        mlst7_db_asked_version,
        serotypefinder_db_asked_version,
        seroba_db_asked_version,
    ):
        if self.update_dbs:
            try:
                rm_dir = subprocess.run(
                    ["rm", "-rf", str(self.db_dir)], check=True, timeout=60
                )
            except:
                rm_dir.kill()
                raise

        software_version = {
            "mlst7": self.download_software_mlst7(version=cge_mlst_asked_version),
            "characterize_neisseria_capsule": self.download_software_characterize_neisseria_capsule(
                version=characterize_neisseria_capsule_asked_version
            ),
            "mlst7_db": self.download_db_mlst7(version=mlst7_db_asked_version),
            "serotypefinder_db": self.download_db_serotypefinder(
                version=serotypefinder_db_asked_version
            ),
            "seroba_db": self.download_db_seroba(
                version=seroba_db_asked_version, kmersize=self.seroba_kmersize
            ),
            "characterize_neisseria_capsule_db": self.copy_neisseria_db(),
        }

        return software_version


if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser(
        description="Download databases for the typing pipeline."
    )
    argument_parser.add_argument(
        "-d",
        "--db-dir",
        type=pathlib.Path,
        default="db",
        help="Database directory where the databases will be stored.",
    )
    argument_parser.add_argument(
        "-mv",
        "--cgemlst-version",
        type=str,
        default="2.0.4",
        help="Version to download for the MLST software from CGE.",
    )
    argument_parser.add_argument(
        "-neis",
        "--neisseria-version",
        type=str,
        default="master",
        help="Version to download characterize neisseria capsule.",
    )
    argument_parser.add_argument(
        "-mdv",
        "--cgemlst-db-version",
        type=str,
        default="master",
        help="Version to download for the MLST database (CGE).",
    )
    argument_parser.add_argument(
        "-sdv",
        "--serotypefinder-db-version",
        type=str,
        default="master",
        help="Version to download for the SerotypeFinder database.",
    )
    argument_parser.add_argument(
        "-serv",
        "--seroba-db-version",
        type=str,
        default="master",
        help="Version to download for the Seroba database.",
    )
    argument_parser.add_argument(
        "--seroba-kmer-size",
        type=int,
        default=71,
        help="Kmer size to be used to build Seroba's database.",
    )
    argument_parser.add_argument("--update", dest="update_dbs", action="store_true")
    args = argument_parser.parse_args()
    downloads = DownloadsJunoTyping(
        db_dir=args.db_dir,
        update_dbs=args.update_dbs,
        cge_mlst_asked_version=args.cgemlst_version,
        characterize_neisseria_capsule_asked_version=args.neisseria_version,
        mlst7_db_asked_version=args.cgemlst_db_version,
        serotypefinder_db_asked_version=args.serotypefinder_db_version,
        seroba_db_asked_version=args.seroba_db_version,
        seroba_kmersize=args.seroba_kmer_size,
    )
    print(downloads.downloaded_versions)
