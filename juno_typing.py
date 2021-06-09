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
import argparse
import pathlib
import subprocess
from sys import path
import yaml

# Own scripts
path.insert(0, 'bin/')
import general_juno_pipeline
import download_dbs

class JunoTypingRun:
    """Class with the arguments and specifications that are only for the Juno-typing pipeline but inherit from PipelineStartup and RunSnakemake"""

    def __init__(self, 
                input_dir, 
                output_dir, 
                db_dir = "db", 
                serotypefinder_mincov=0.6, 
                serotypefinder_identity=0.85,
                seroba_mincov=20, 
                seroba_kmersize=71,
                cores=8,
                local=True,
                queue='bio',
                unlock=False,
                rerunincomplete=False,
                dryrun=False,
                update=False):
        """Initiating Juno-typing pipeline"""

        # Pipeline attributes
        self.pipeline_info = {'pipeline_name': "Juno-typing",
                                'pipeline_version': "0.1"}
        self.snakefile = "Snakefile"
        self.workdir = pathlib.Path(__file__).parent.absolute()
        self.useconda = True
        self.usesingularity = False
        self.singularityargs = ""
        self.user_parameters = pathlib.Path("config/user_parameters.yaml")
        self.extra_software_versions = pathlib.Path('config/extra_software_versions.yaml')
        self.output_dir = output_dir
        self.restarttimes = 2       
        # Checking if the input_dir comes from the Juno-assembly pipeline 
        self.startup = general_juno_pipeline.PipelineStartup(input_dir, 'both')
        
        # Parse arguments specific from the user
        self.user_params = self.write_userparameters(input_dir,
                                                    self.output_dir,
                                                    db_dir,
                                                    serotypefinder_mincov,
                                                    serotypefinder_identity,
                                                    seroba_mincov, 
                                                    seroba_kmersize,
                                                    self.user_parameters)
        
        # Download databases if necessary
        if not unlock and not dryrun:
            self.download_databases(db_dir, 
                                    update,
                                    self.extra_software_versions,
                                    seroba_kmersize)

        # Run snakemake
        general_juno_pipeline.RunSnakemake(pipeline_name = self.pipeline_info['pipeline_name'],
                                            pipeline_version = self.pipeline_info['pipeline_version'],
                                            sample_dict = self.startup.sample_dict,
                                            output_dir = self.output_dir,
                                            workdir = self.workdir,
                                            snakefile = self.snakefile,
                                            cores = cores,
                                            local = local,
                                            queue = queue,
                                            unlock = unlock,
                                            rerunincomplete = rerunincomplete,
                                            dryrun = dryrun,
                                            useconda = self.useconda,
                                            usesingularity = self.usesingularity,
                                            singularityargs = self.singularityargs,
                                            restarttimes = self.restarttimes)
        
    
    
    def download_databases(self, 
                            db_dir, 
                            updatedbs,
                            log_software_versions,
                            seroba_kmersize):
        """Function to download software and databases necessary for running the Juno-typing pipeline"""
        # if updatedbs:
        #     update = "TRUE"
        # else:
        #     update = "FALSE"

        db_dir.mkdir(parents = True, exist_ok = True)
        download_dbs.get_downloads_juno_typing(db_dir, 'bin', updatedbs, seroba_kmersize)
        # download_dbs = subprocess.run(['bash', 'bin/download_dbs.sh', str(db_dir), update, log_software_versions, str(seroba_kmersize)])
        # download_dbs.wait()


    def start_pipeline(self, 
                        input_sudir, 
                        sample_sheet):
        """Function to start the pipeline (some steps from PipelineStartup need to be modified for the Juno-typing pipeline to accept fastq and fasta input"""
        # Taking fastq input as the Startup just to inherit all the same attributes
        # from parent class (PipelineStartup). The fasta sample sheet is created 
        # separately and appended to the original one
        startup = general_juno_pipeline.PipelineStartup(input_dir, input_type = 'both')

        with open(sample_sheet, 'w') as file:
            yaml.dump(startup.sample_dict, file, default_flow_style=False)
        
        return startup.sample_dict
    
    def write_userparameters(self,
                            input_dir,
                            output_dir,
                            db_dir,
                            serotypefinder_mincov,
                            serotypefinder_identity,
                            seroba_mincov, 
                            seroba_kmersize,
                            user_parameters):

        config_params = {'input_dir': str(input_dir),
                        'out': str(output_dir),
                        'kmerfinder_db': str(db_dir.joinpath('kmerfinder_db')),
                        'mlst7_db': str(db_dir.joinpath('mlst7_db')),
                        'seroba_db': str(db_dir.joinpath('seroba_db')),
                        'serotypefinder_db': str(db_dir.joinpath('serotypefinder_db')),
                        'serotypefinder': {'min_cov': serotypefinder_mincov,
                                            'identity_thresh': serotypefinder_identity},
                        'seroba': {'min_cov': seroba_mincov,
                                    'kmer_size': seroba_kmersize}}
        
        with open(user_parameters, 'w') as file:
            yaml.dump(config_params, file, default_flow_style=False)

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
        "-sh",
        "--snakemake_help",
        action = 'store_true',
        help = "Print Snakemake help (passed to snakemake)."
    )
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
    args = parser.parse_args()
    JunoTypingRun(args.input, 
                    args.output, 
                    args.db_dir,
                    args.serotypefinder_mincov,
                    args.serotypefinder_identity,
                    args.seroba_mincov,
                    args.seroba_kmersize,
                    args.cores,
                    args.local,
                    args.queue,
                    args.unlock,
                    args.rerunincomplete,
                    args.dryrun,
                    args.update)
