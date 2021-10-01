"""
Juno-typing pipeline
Authors: Alejandra Hernandez-Segura
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium
            Surveillance (IDS), Bacteriologie (BPD)     
Date: 18-05-2021   

Documentation: https://rivm-bioinformatics.github.io/ids_bacteriology_man/juno-typing.html 

Snakemake rules (in order of execution): 
    1 identify_species  Using kmerFinder 
    2 mlst7             7-locus MLST taking either filtered fastq files
                        or fasta files (assemblies) as input. Both cases
                        using the MLST software from CGE.
    3 Serotyper         Either SerotypeFinder if the sample is E. coli,
                        SeqSero2 if the sample is Salmonella and 
                        Seroba if the sample is S. pneumoniae.
    4 Multiserotyper    Reports collecting the results of each type 
        report          of serotyper.
"""

# Dependencies
from base_juno_pipeline import *
import argparse
import pandas as pd
import pathlib
import subprocess
import yaml

# Own scripts
import bin.download_dbs

class JunoTypingRun(base_juno_pipeline.PipelineStartup,
                    base_juno_pipeline.RunSnakemake):
    """Class with the arguments and specifications that are only for the Juno-typing pipeline but inherit from PipelineStartup and RunSnakemake"""
    
    def __init__(self, 
                input_dir, 
                output_dir, 
                metadata = None, 
                db_dir = "/mnt/db/juno/typing_db", 
                update_dbs=False,
                serotypefinder_mincov=0.6, 
                serotypefinder_identity=0.85,
                seroba_mincov=20, 
                seroba_kmersize=71,
                cores=300,
                local=False,
                queue='bio',
                unlock=False,
                rerunincomplete=False,
                dryrun=False,
                run_in_container=False,
                singularity_prefix=None,
                conda_prefix=None,
                **kwargs):
        """Initiating Juno-typing pipeline"""
        
        # Get proper file paths
        output_dir = pathlib.Path(output_dir).resolve()
        workdir = pathlib.Path(__file__).parent.resolve()
        self.db_dir = pathlib.Path(db_dir).resolve()
        self.path_to_audit = output_dir.joinpath('audit_trail')
        base_juno_pipeline.PipelineStartup.__init__(self,
            input_dir=pathlib.Path(input_dir).resolve(), 
            input_type='both',
            min_num_lines=2) # Min for viable fasta
        base_juno_pipeline.RunSnakemake.__init__(self,
            pipeline_name='Juno_typing',
            pipeline_version='v0.2',
            output_dir=output_dir,
            workdir=workdir,
            cores=cores,
            local=local,
            queue=queue,
            unlock=unlock,
            rerunincomplete=rerunincomplete,
            dryrun=dryrun,
            useconda=not run_in_container,
            conda_prefix=conda_prefix,
            usesingularity=run_in_container,
            singularityargs=f"--bind {self.input_dir}:{self.input_dir} --bind {output_dir}:{output_dir} --bind {db_dir}:{db_dir}",
            singularity_prefix=singularity_prefix,
            restarttimes=0,
            latency_wait=60,
            name_snakemake_report=str(self.path_to_audit.joinpath('juno_typing_report.html')),
            **kwargs)
        # self.pipeline_name = 'Juno_typing'
        # self.pipeline_version = '0.2'
        # self.output_dir = pathlib.Path(output_dir).resolve()
        # self.workdir = pathlib.Path(__file__).parent.resolve()
        # self.sample_sheet = "config/sample_sheet.yaml"
        # self.user_parameters = pathlib.Path('config/user_parameters.yaml')
        # self.fixed_parameters = pathlib.Path('config/pipeline_parameters.yaml')
        # self.snakefile = 'Snakefile'
        # self.cores = cores
        # self.local = local
        # self.snakemake_report = str(self.path_to_audit.joinpath('juno_typing_report.html'))
        # self.queue = queue
        # self.unlock = unlock
        # self.dryrun = dryrun
        # self.rerunincomplete = rerunincomplete
        # if run_in_container:
        #     self.useconda = False
        #     self.usesingularity = True
        # else:
        #     self.useconda = True
        #     self.usesingularity = False
        # self.conda_frontend = 'mamba'
        # self.conda_prefix = conda_prefix
        # self.singularityargs = f"--bind {self.input_dir}:{self.input_dir} --bind {self.output_dir}:{self.output_dir} --bind {self.db_dir}:{self.db_dir}"
        # self.singularity_prefix = singularity_prefix
        # self.restarttimes = 1 # Set to 1 because kmerfinder sometime fails. Once it is replaced, it can be 0
        # self.latency = 60
        # self.kwargs = kwargs

        # Pipeline attributes
        self.metadata_file = metadata
        self.serotypefinder_mincov=serotypefinder_mincov
        self.serotypefinder_identity=serotypefinder_identity
        self.seroba_mincov=seroba_mincov
        self.seroba_kmersize = seroba_kmersize
        self.update_dbs = update_dbs

        # Start pipeline
        self.start_juno_typing_pipeline()
        self.user_params = self.write_userparameters()
        self.get_run_info()
        if not self.dryrun or self.unlock:
            self.path_to_audit.mkdir(parents=True, exist_ok=True)
            self.audit_trail = self.generate_audit_trail()
            downloads_juno_typing = bin.download_dbs.DownloadsJunoTyping(self.db_dir,
                                                                        update_dbs=self.update_dbs,
                                                                        kmerfinder_asked_version='3.0.2',
                                                                        cge_mlst_asked_version='2.0.4',
                                                                        kmerfinder_db_asked_version='20210425',
                                                                        mlst7_db_asked_version='master',
                                                                        serotypefinder_db_asked_version='master',
                                                                        seroba_db_asked_version='master')
            self.downloads_versions = downloads_juno_typing.downloaded_versions
            with open(self.path_to_audit.joinpath('database_versions.yaml'), 'w') as file_:
                yaml.dump(self.downloads_versions, file_, default_flow_style=False)

        self.successful_run = self.run_snakemake()
        assert self.successful_run, f'Please check the log files'
        if not self.dryrun or self.unlock:
            subprocess.run(['find', self.output_dir, '-type', 'f', '-empty', '-exec', 'rm', '{}', ';'])
            subprocess.run(['find', self.output_dir, '-type', 'd', '-empty', '-exec', 'rm', '-rf', '{}', ';'])
            self.make_snakemake_report()

    def get_mlst7_scheme_name(self, sample_dict):
        with open("files/dictionary_correct_species.yaml") as translation_yaml:
            self.mlst7_species_translation_tbl = yaml.safe_load(translation_yaml)
        for sample in sample_dict:
            if sample_dict[sample]['genus'] != 'NotProvided' and sample_dict[sample]['species'] != 'NotProvided':
                mlst7_species = sample_dict[sample]['genus'][0] + sample_dict[sample]['species']
                try:
                    new_species = self.mlst7_species_translation_tbl[mlst7_species]
                    yield new_species
                except KeyError:
                    yield mlst7_species

    def update_sample_dict_with_metadata(self):
        self.get_metadata_from_juno_assembly(filepath=self.metadata_file)
        # Add metadata
        for sample in self.sample_dict:
            try:
                self.sample_dict[sample].update(self.juno_metadata[sample])
            except (KeyError, TypeError):
                self.sample_dict[sample]['genus'] = 'NotProvided'
                self.sample_dict[sample]['species'] = 'NotProvided'
        # The list does not return anything meaningful. It is just to activate
        # the generator
        list(self.get_mlst7_scheme_name(sample_dict = self.sample_dict))

    def start_juno_typing_pipeline(self):
        """
        Function to start the pipeline (some steps from PipelineStartup need 
        to be modified for the Juno-typing pipeline to accept fastq and fasta 
        input
        """
        self.start_juno_pipeline()
        self.update_sample_dict_with_metadata()
        # Write sample_sheet
        with open(self.sample_sheet, 'w') as file_:
            yaml.dump(self.sample_dict, file_, default_flow_style=False)
    
    def write_userparameters(self):

        config_params = {'input_dir': str(self.input_dir),
                        'out': str(self.output_dir),
                        'kmerfinder_db': str(self.db_dir.joinpath('kmerfinder_db')),
                        'mlst7_db': str(self.db_dir.joinpath('mlst7_db')),
                        'seroba_db': str(self.db_dir.joinpath('seroba_db')),
                        'serotypefinder_db': str(self.db_dir.joinpath('serotypefinder_db')),
                        'serotypefinder': {'min_cov': self.serotypefinder_mincov,
                                            'identity_thresh': self.serotypefinder_identity},
                        'seroba': {'min_cov': self.seroba_mincov,
                                    'kmer_size': self.seroba_kmersize},
                        'cgmlst_db': str(self.db_dir.joinpath('cgmlst'))}
        
        with open(self.user_parameters, 'w') as file_:
            yaml.dump(config_params, file_, default_flow_style=False)

        return config_params     
        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Juno-typing pipeline. Automated pipeline for bacterial subtyping (7-locus MLST and serotyping)."
    )
    parser.add_argument(
        "-i",
        "--input",
        type = pathlib.Path,
        required = True,
        metavar = "DIR",
        help = "Relative or absolute path to the input directory. It must either be the output directory of the Juno-assembly pipeline or it must contain all the raw reads (fastq) and assemblies (fasta) files for all samples to be processed."
    )
    parser.add_argument(
        "-m",
        "--metadata",
        type = pathlib.Path,
        default = None,
        metavar = "FILE",
        help = "Relative or absolute path to the metadata csv file. If provided, it must contain at least one column with the 'Sample' name (name of the file but removing _R1.fastq.gz), a column called 'Genus' and a column called 'Species'. The genus and species provided will be used to choose the serotyper and the MLST schema."
    )
    parser.add_argument(
        "-o",
        "--output",
        type = pathlib.Path,
        metavar = "DIR",
        default = "output",
        help = "Relative or absolute path to the output directory. If non is given, an 'output' directory will be created in the current directory."
    )
    parser.add_argument(
        "-d",
        "--db_dir",
        type = pathlib.Path,
        required = False,
        metavar = "DIR",
        default = "/mnt/db/juno/typing_db",
        help = "Relative or absolute path to the directory that contains the databases for all the tools used in this pipeline or where they should be downloaded. Default is: /mnt/db/juno/typing_db",
    )
    parser.add_argument(
        "--serotypefinder_mincov",
        type = float,
        metavar = "NUM",
        default = 0.6,
        help = "Minimum coverage to be used for the SerotypeFinder (E. coli serotyping) tool. It accepts values from 0-1. Default is 0.6.",
    )
    parser.add_argument(
        "--serotypefinder_identity",
        type = float,
        metavar = "NUM",
        default = 0.85,
        help = "Identity threshold to be used for the SerotypeFinder (E. coli serotyping) tool. It accepts values from 0-1. Default is 0.85",
    )
    parser.add_argument(
        "--seroba_mincov",
        type = int,
        metavar = "INT",
        default = 20,
        help = "Minimum coverage to be used for the Seroba (S. pneumoniae serotyping) tool. It accepts values from 0-100. Default is 20",
    )
    parser.add_argument(
        "--seroba_kmersize",
        type = int,
        metavar = "INT",
        default = 71,
        help = "Minimum coverage to be used for the SerotypeFinder (E. coli serotyping) tool. It accepts values from 0-100. Default is 20",
    )
    parser.add_argument(
        "--update",
        action='store_true',
        help="Force database update even if they are present."
    )
    parser.add_argument(
        "-c",
        "--cores",
        type = int,
        metavar = "INT",
        default = 300,
        help="Number of cores to use. Default is 300"
    )
    #TODO: Get from ${irods_runsheet_sys__runsheet__lsf_queue} if it exists
    parser.add_argument(
        "-q",
        "--queue",
        type = str,
        metavar = "STR",
        default = 'bio',
        help = 'Name of the queue that the job will be submitted to if working on a cluster.'
    )
    parser.add_argument(
        "-l",
        "--local",
        action='store_true',
        help="Running pipeline locally (instead of in a computer cluster). Default is running it in a cluster."
    )
    # Snakemake arguments
    parser.add_argument(
        "-u",
        "--unlock",
        action = 'store_true',
        help = "Unlock output directory (passed to snakemake)."
    )
    parser.add_argument(
        "-n",
        "--dryrun",
        action='store_true',
        help="Dry run printing steps to be taken in the pipeline without actually running it (passed to snakemake)."
    )
    parser.add_argument(
        "--rerunincomplete",
        action='store_true',
        help="Re-run jobs if they are marked as incomplete (passed to snakemake)."
    )
    parser.add_argument(
        "--snakemake-args",
        nargs='*',
        default={},
        action=helper_functions.SnakemakeKwargsAction,
        help="Extra arguments to be passed to snakemake API (https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html)."
    )
    args = parser.parse_args()
    JunoTypingRun(input_dir = args.input, 
                    output_dir = args.output, 
                    metadata = args.metadata,
                    db_dir = args.db_dir,
                    update_dbs = args.update,
                    serotypefinder_mincov = args.serotypefinder_mincov,
                    serotypefinder_identity = args.serotypefinder_identity,
                    seroba_mincov = args.seroba_mincov,
                    seroba_kmersize = args.seroba_kmersize,
                    cores = args.cores,
                    local = args.local,
                    queue = args.queue,
                    unlock = args.unlock,
                    rerunincomplete = args.rerunincomplete,
                    dryrun = args.dryrun,
                    run_in_container=False,
                    singularity_prefix=None,
                    conda_prefix=None,
                    **args.snakemake_args)
