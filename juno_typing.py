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
import pandas as pd
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
                metadata, 
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
        self.sample_sheet = "config/sample_sheet.yaml"
        self.input_dir = pathlib.Path(input_dir)
        self.output_dir = pathlib.Path(output_dir)
        self.metadata = metadata
        self.db_dir = pathlib.Path(db_dir)
        self.serotypefinder_mincov=serotypefinder_mincov
        self.serotypefinder_identity=serotypefinder_identity
        self.seroba_mincov=seroba_mincov
        self.seroba_kmersize = seroba_kmersize
        self.update = update
        self.workdir = pathlib.Path(__file__).parent.absolute()
        self.useconda = True
        self.usesingularity = False
        self.singularityargs = ""
        self.user_parameters = pathlib.Path("config/user_parameters.yaml")
        self.extra_software_versions = pathlib.Path('config/extra_software_versions.yaml')
        self.output_dir = output_dir
        self.restarttimes = 2       
        # Checking if the input_dir comes from the Juno-assembly pipeline 
        self.startup = self.start_pipeline()

        # Parse arguments specific from the user
        self.user_params = self.write_userparameters()
        
        # Download databases if necessary
        if not unlock and not dryrun:
            self.download_databases()

        # Run snakemake
        general_juno_pipeline.RunSnakemake(pipeline_name = self.pipeline_info['pipeline_name'],
                                            pipeline_version = self.pipeline_info['pipeline_version'],
                                            sample_sheet = self.sample_sheet,
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
        
    def download_databases(self):
        """Function to download software and databases necessary for running the Juno-typing pipeline"""
        # if updatedbs:
        #     update = "TRUE"
        # else:
        #     update = "FALSE"

        self.db_dir.mkdir(parents = True, exist_ok = True)
        download_dbs.get_downloads_juno_typing(self.db_dir, 'bin', self.update, self.seroba_kmersize)
        # download_dbs = subprocess.run(['bash', 'bin/download_dbs.sh', str(self.db_dir), self.update, self.log_software_versions, str(self.seroba_kmersize)])
        # download_dbs.wait()


    def add_metadata(self, samples_dic):
        assert self.metadata.is_file(), "Provided metadata file ({}) does not exist".format(self.metadata)
        # Load species file
        species_dic = pd.read_csv(self.metadata, index_col = 0, dtype={'Sample':str})
        species_dic.index = species_dic.index.map(str)
        species_dic['genus-abbr'] = species_dic['Genus'].apply(lambda x: x[0])
        species_dic['species-mlst7'] = species_dic['genus-abbr'] + species_dic["Species"]
        species_dic['species-mlst7'] = species_dic['species-mlst7'].apply(lambda x: x.lower())
        species_dic = species_dic.transpose().to_dict()
        # Update dictionary with species
        for sample_name in samples_dic :
            try: 
                samples_dic[sample_name].update(species_dic[sample_name])
            except:
                pass
        return species_dic

    def start_pipeline(self):
        """Function to start the pipeline (some steps from PipelineStartup need to be modified for the Juno-typing pipeline to accept fastq and fasta input"""
        # Taking fastq input as the Startup just to inherit all the same attributes
        # from parent class (PipelineStartup). The fasta sample sheet is created 
        # separately and appended to the original one
        startup = general_juno_pipeline.PipelineStartup(self.input_dir, input_type = 'both')
        # Add species-mlst7 data if a metadata file was provided or indicate so if not provided
        for sample in startup.sample_dict:
            startup.sample_dict[sample]['species-mlst7'] = "NotProvided"
        if self.metadata is not None:
            startup.sample_dict = self.add_metadata(startup.sample_dict)
        # Write sample_sheet
        with open(self.sample_sheet, 'w') as file:
            yaml.dump(startup.sample_dict, file, default_flow_style=False)
        return startup
    
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
                                    'kmer_size': self.seroba_kmersize}}
        
        with open(self.user_parameters, 'w') as file:
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
        "-m",
        "--metadata",
        type = pathlib.Path,
        default = None,
        metavar = "FILE",
        help = "Relative or absolute path to the metadata csv file. If provided, it must contain at least one column with the 'Sample' name (name of the file but removing _R1.fastq.gz), a column called 'Genus' and a column called 'Species'. If a genus + species is provided for a sample, the MLST7 is not run."
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
    args = parser.parse_args()
    JunoTypingRun(args.input, 
                    args.output, 
                    args.metadata,
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
