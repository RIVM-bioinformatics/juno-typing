import csv
import os
import shutil
import tempfile
from pathlib import Path
from sys import path
import unittest

main_script_path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
downloads_db_path = str(Path(__file__).parent.parent.absolute().joinpath("bin"))
path.insert(0, main_script_path)
path.insert(0, downloads_db_path)

from juno_typing import JunoTyping
from download_dbs import DownloadsJunoTyping

# from ..bin.download_dbs import DownloadsJunoTyping


@unittest.skipIf(
    not Path(
        "/data/BioGrid/hernanda/test_data_per_pipeline/Enteric/Juno-typing/"
    ).exists(),
    "Skipped in non-RIVM environments (for the sake of time)",
)
class TestDownloadDbs(unittest.TestCase):
    """Testing the downloading of databases and software used by the pipeline"""

    @classmethod
    def setUpClass(cls) -> None:
        Path("fake_db").mkdir(exist_ok=True)

    @classmethod
    def tearDownClass(cls) -> None:
        os.system("rm -rf fake_db")

    def test_download_all_dbs(self) -> None:
        """Testing all databases/software are downloaded and set up"""

        path_to_db = Path("fake_db")
        path_to_bin = Path("bin")
        downloads = DownloadsJunoTyping(
            path_to_db,
            update_dbs=True,
            cge_mlst_asked_version="2.0.4",
            mlst7_db_asked_version="master",
            serotypefinder_db_asked_version="master",
            seroba_db_asked_version="master",
            seroba_kmersize=50,
        )
        self.assertEqual(downloads.seroba_kmersize, 50)
        self.assertEqual(downloads.downloaded_versions["mlst7"], "2.0.4")
        self.assertTrue(path_to_bin.joinpath("cge-mlst", "mlst.py").exists())
        self.assertTrue(
            path_to_db.joinpath("mlst7_db", "senterica", "senterica.length.b").exists()
        )
        self.assertTrue(
            path_to_db.joinpath("serotypefinder_db", "H_type.seq.b").exists()
        )
        self.assertTrue(
            path_to_db.joinpath("seroba_db", "database", "cdhit_cluster").exists()
        )


class TestCopyNeisseriaDb(unittest.TestCase):

    def setUp(self) -> None:
        self.tmp = Path(tempfile.mkdtemp())
        self.fake_source = self.tmp / "neisseria_capsule_DB"
        self.fake_source.mkdir()
        (self.fake_source / "marker.txt").write_text("hello")
        self.fake_bin = self.tmp / "bin"
        (self.fake_bin / "characterize_neisseria_capsule").mkdir(parents=True)

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp)

    def _make_instance(self, source: Path) -> DownloadsJunoTyping:
        instance = DownloadsJunoTyping.__new__(DownloadsJunoTyping)
        instance.bin_dir = self.fake_bin
        instance.neisseria_db_source = source
        return instance

    def test_copies_from_overridden_source(self) -> None:
        instance = self._make_instance(self.fake_source)
        instance.copy_neisseria_db()
        copied_marker = (
            self.fake_bin
            / "characterize_neisseria_capsule"
            / "neisseria_capsule_DB"
            / "marker.txt"
        )
        self.assertTrue(copied_marker.exists())
        self.assertEqual(copied_marker.read_text(), "hello")

    def test_copies_from_arbitrarily_named_source(self) -> None:
        renamed = self.tmp / "my_custom_db"
        renamed.mkdir()
        (renamed / "marker.txt").write_text("hello")
        instance = self._make_instance(renamed)
        instance.copy_neisseria_db()
        copied_marker = (
            self.fake_bin
            / "characterize_neisseria_capsule"
            / "neisseria_capsule_DB"
            / "marker.txt"
        )
        self.assertTrue(copied_marker.exists())

    def test_missing_source_raises_clear_error(self) -> None:
        instance = self._make_instance(self.tmp / "does_not_exist")
        with self.assertRaises(FileNotFoundError) as ctx:
            instance.copy_neisseria_db()
        self.assertIn("does_not_exist", str(ctx.exception))

    def test_default_preserves_rivm_hpc_path(self) -> None:
        import inspect
        default = inspect.signature(
            DownloadsJunoTyping.__init__
        ).parameters["neisseria_db_source"].default
        self.assertEqual(str(default), "/mnt/db/juno/neisseria_capsule_DB")


class TestJunoTypingDryRun(unittest.TestCase):
    """Testing the JunoTyping class (code specific for this pipeline)"""

    @classmethod
    def setUpClass(cls) -> None:
        fake_dirs = [
            "fake_dir_wsamples",
            "fake_dir_juno",
            "fake_dir_juno/clean_fastq",
            "fake_dir_juno/de_novo_assembly_filtered",
            "fake_db",
            "fake_db/mlst7_db/",
            "fake_db/mlst7_db/senterica/",
        ]

        fake_files = [
            "fake_dir_wsamples/sample1_R1.fastq",
            "fake_dir_wsamples/sample1_R2.fastq.gz",
            "fake_dir_wsamples/sample2_R1_filt.fq",
            "fake_dir_wsamples/sample2_R2_filt.fq.gz",
            "fake_dir_wsamples/sample1.fasta",
            "fake_dir_wsamples/sample2.fasta",
            "fake_dir_juno/clean_fastq/1234_R1.fastq.gz",
            "fake_dir_juno/clean_fastq/1234_R2.fastq.gz",
            "fake_dir_juno/de_novo_assembly_filtered/1234.fasta",
            "fake_db/mlst7_db/senterica/senterica.length.b",
        ]

        for folder in fake_dirs:
            Path(folder).mkdir(exist_ok=True)
        for file_ in fake_files:
            with open(file_, mode="w") as metadata_file:
                metadata_writer = csv.writer(
                    metadata_file,
                    delimiter=",",
                    quotechar='"',
                    quoting=csv.QUOTE_MINIMAL,
                )
                metadata_writer.writerow(["content1"])
                metadata_writer.writerow(["content2"])

        with open("fake_dir_wsamples/fake_metadata.csv", mode="w") as metadata_file:
            metadata_writer = csv.writer(
                metadata_file, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
            )
            metadata_writer.writerow(["sample", "genus", "species"])
            metadata_writer.writerow(["sample1", "Salmonella", "enterica"])
            metadata_writer.writerow(["sample2", "Escherichia", "coli"])

    @classmethod
    def tearDownClass(cls) -> None:
        fake_dirs = ["fake_dir_wsamples", "fake_dir_juno", "fake_db", "test_output"]
        for folder in fake_dirs:
            os.system("rm -rf {}".format(str(folder)))

    def test_junotyping_dryrun(self) -> None:
        """Testing the pipeline runs properly as a dry run"""
        argv = [
            "-i",
            "fake_dir_wsamples",
            "-o",
            "test_output",
            "-n",
            "--species",
            "Salmonella",
            "enterica",
            "--db_dir",
            "fake_db",
        ]
        juno_typing = JunoTyping(argv=argv)
        juno_typing.run()

    def test_junotyping_dryrun_wMetadata(self) -> None:
        """Testing the pipeline runs properly as a dry run when providing a metadata file"""
        argv = [
            "-i",
            "fake_dir_wsamples",
            "-o",
            "test_output",
            "-n",
            "--species",
            "Salmonella",
            "enterica",
            "--db_dir",
            "fake_db",
            "--metadata",
            "fake_dir_wsamples/fake_metadata.csv",
        ]
        juno_typing = JunoTyping(argv=argv)
        juno_typing.run()


@unittest.skipIf(
    not Path(
        "/data/BioGrid/hernanda/test_data_per_pipeline/Enteric/Juno-typing/"
    ).exists(),
    "Skipped in non-RIVM environments (because test data is needed)",
)
class TestJunoTypingPipeline(unittest.TestCase):
    """Testing the JunoTyping class (code specific for this pipeline)"""

    @classmethod
    def tearDownClass(cls) -> None:
        os.system("rm -rf test_output")
        os.system("rm -rf test_output2")
        os.system("rm -rf test_new_conda_envs")

    def test_junotyping_run(self) -> None:
        """Testing the pipeline runs properly with real samples and new conda envs"""

        output_dir = Path("test_output")
        argv = [
            "-i",
            "/data/BioGrid/hernanda/test_data_per_pipeline/Enteric/Juno-typing/",
            "-o",
            str(output_dir),
            "--db_dir",
            "/mnt/db/juno/typing_db",
            "--prefix",
            "test_new_conda_envs",
        ]
        juno_typing = JunoTyping(argv=argv)
        juno_typing.run()

        self.assertTrue(
            output_dir.joinpath("serotype", "serotyper_multireport.csv").exists()
        )
        self.assertTrue(
            output_dir.joinpath("serotype", "serotyper_multireport1.csv").exists()
        )
        self.assertTrue(
            output_dir.joinpath("serotype", "serotyper_multireport2.csv").exists()
        )
        self.assertTrue(
            output_dir.joinpath("serotype", "serotyper_multireport3.csv").exists()
        )
        self.assertTrue(output_dir.joinpath("mlst7", "mlst7_multireport.csv").exists())
        self.assertTrue(output_dir.joinpath("audit_trail", "log_git.yaml").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "log_pipeline.yaml").exists()
        )
        self.assertTrue(output_dir.joinpath("audit_trail", "log_conda.txt").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "snakemake_report.html").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "sample_sheet.yaml").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "user_parameters.yaml").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "database_versions.yaml").exists()
        )

    def test_junotyping_run_wMetadata(self) -> None:
        """Testing the pipeline runs properly with real samples when providing a metadata file"""
        output_dir = Path("test_output2")
        argv = [
            "-i",
            "/data/BioGrid/hernanda/test_data_per_pipeline/Enteric/Juno-typing/",
            "-o",
            str(output_dir),
            "--db_dir",
            "/mnt/db/juno/typing_db",
            "--prefix",
            "test_new_conda_envs",
            "--metadata",
            "/data/BioGrid/hernanda/test_data_per_pipeline/Enteric/Juno-typing/metadata.csv",
        ]
        juno_typing = JunoTyping(argv=argv)
        juno_typing.run()

        self.assertTrue(
            output_dir.joinpath("serotype", "serotyper_multireport.csv").exists()
        )
        self.assertTrue(
            output_dir.joinpath("serotype", "serotyper_multireport1.csv").exists()
        )
        self.assertTrue(
            output_dir.joinpath("serotype", "serotyper_multireport2.csv").exists()
        )
        self.assertTrue(
            output_dir.joinpath("serotype", "serotyper_multireport3.csv").exists()
        )
        self.assertTrue(output_dir.joinpath("mlst7", "mlst7_multireport.csv").exists())
        self.assertTrue(output_dir.joinpath("audit_trail", "log_git.yaml").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "log_pipeline.yaml").exists()
        )
        self.assertTrue(output_dir.joinpath("audit_trail", "log_conda.txt").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "snakemake_report.html").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "sample_sheet.yaml").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "user_parameters.yaml").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "database_versions.yaml").exists()
        )


if __name__ == "__main__":
    unittest.main()
