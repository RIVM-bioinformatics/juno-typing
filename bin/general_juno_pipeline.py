from datetime import datetime
import pathlib
import re
import snakemake
import subprocess
from uuid import uuid4
import yaml

def print_message(message):
    message = '\033[0;33m' + message + '\n\033[0;0m'
    print(message)

def print_error(message):
    message = '\033[0;31m' + message + '\n\033[0;0m'
    print(message)

class PipelineStartup:
    """Class with checkings for the input directory and to generate a dictionary 
        with sample names and their corresponding input files (sample_dict). 
        Works for pipelines accepting fastq or fasta files as input"""
    
    def __init__(self,
                input_dir, 
                input_type = 'fastq'):

        if not isinstance(input_dir, pathlib.PosixPath):
            input_dir = pathlib.Path(input_dir)
        assert input_dir.is_dir(), "The provided input directory ({}) does not exist. Please provide an existing directory".format(input_dir)
        assert input_type in ['fastq', 'fasta', 'both'], "input_type to be checked can only be 'fastq' or 'fasta'"
        
        self.input_dir = input_dir
        self.input_type = input_type
        self.unique_id = uuid4()
        self.__subdirs_ = self.find_input_subdirs()
        self.check_input_dir()
        self.sample_dict = self.make_sample_dict()

    def find_input_subdirs(self):
        """Function to check whether the input is from the Juno annotation pipeline or just a simple input directory"""
        # Only when the input_dir comes from the Juno-annotation pipeline the fastq and fasta files
        # do not need to be in the same folder (fastq are then expected in a subfolder called <input_dir>/clean_fastq
        # and the )
        if self.input_dir.joinpath('clean_fastq').exists() and self.input_dir.joinpath('de_novo_assembly_filtered').exists():
            return {'fastq': self.input_dir.joinpath('clean_fastq'),
                    'fasta': self.input_dir.joinpath('de_novo_assembly_filtered')}
        else:
            return {'fastq': self.input_dir,
                    'fasta': self.input_dir}

    def check_contains_file_wextension(self, input_subdir, extension = ('.fastq', '.fastq.gz', '.fq', '.fq.gz', '.fasta')):
        for item in input_subdir.iterdir():
            if item.is_file():
                if str(item).endswith(extension):
                    return True
                    break
        
        raise ValueError('Input directory ({}) does not contain files that end with one of the expected extensions {}.'.format(input_subdir, extension))

    def check_input_dir(self):
        """Function to check that input directory is indeed an existing directory that contains files with the expected extension (fastq or fasta)"""
    
        if self.input_type == 'fastq':
            return self.check_contains_file_wextension(self.__subdirs_['fastq'], ('.fastq', '.fastq.gz', '.fq', '.fq.gz'))
        elif self.input_type == 'fasta':
            return self.check_contains_file_wextension(self.__subdirs_['fasta'], ('.fasta'))
        else:
            return self.check_contains_file_wextension(self.__subdirs_['fastq'], 
                                                        ('.fastq', '.fastq.gz', '.fq', '.fq.gz')) and self.check_contains_file_wextension(self.__subdirs_['fasta'], 
                                                                                                                                            ('.fasta'))

    def enlist_fastq_samples(self):
        """Function to enlist the fastq files found in the input directory. Returns a dictionary with the form {sample: fastq_file}"""

        pattern = re.compile("(.*?)(?:_S\d+_|_S\d+.|_|\.)(?:p)?R?(1|2)(?:_.*\.|\..*\.|\.)f(ast)?q(\.gz)?")
        samples = {}
        for file_ in self.__subdirs_['fastq'].iterdir():
            if file_.is_dir():
                continue
            match = pattern.fullmatch(file_.name)
            if match:
                sample = samples.setdefault(match.group(1), {})
                sample["R{}".format(match.group(2))] = str(file_)               
        return samples

    def enlist_fasta_samples(self):
        """Function to enlist the fasta files found in the input directory. Returns a dictionary with the form {sample: fasta_file}"""

        pattern = re.compile("(.*?).fasta")
        samples = {}
        for file_ in self.__subdirs_['fasta'].iterdir():
            if file_.is_dir():
                continue
            match = pattern.fullmatch(file_.name)
            if match:
                sample = samples.setdefault(match.group(1), {})
                sample["assembly"] = str(file_)
        return samples


    def make_sample_dict(self):
        """Function to make a sample sheet from the input directory (expecting either fastq or fasta files as input)"""

        if self.input_type == 'fastq':
            samples = self.enlist_fastq_samples()
        elif self.input_type == 'fasta':
            samples = self.enlist_fasta_samples()
        else:
            samples = self.enlist_fastq_samples()
            samples_fasta = self.enlist_fasta_samples()
            for k in samples.keys():
                samples[k]['assembly'] = samples_fasta[k]['assembly']
        return samples



class RunSnakemake:
    """Class with necessary input to run Snakemake"""

    def __init__(self,
                pipeline_name,
                pipeline_version,
                output_dir,
                workdir,
                sample_sheet=pathlib.Path('config/sample_sheet.yaml'), 
                user_parameters=pathlib.Path('config/user_parameters.yaml'), 
                fixed_parameters=pathlib.Path('config/pipeline_parameters.yaml'),
                snakefile='Snakefile',
                cores=300,
                local=False,
                queue='bio',
                unlock=False,
                rerunincomplete=True,
                dryrun=False,
                useconda=True,
                conda_frontend='mamba',
                usesingularity=True,
                singularityargs='',
                restarttimes=0,
                latency_wait=60):
        self.pipeline_name = pipeline_name
        self.pipeline_version = pipeline_version
        self.output_dir = output_dir
        self.workdir = workdir
        self.sample_sheet = sample_sheet
        self.user_parameters = user_parameters
        self.fixed_parameters = fixed_parameters
        self.snakefile = snakefile
        self.path_to_audit = self.output_dir.joinpath('audit_trail')
        self.cores = cores
        self.local = local
        self.queue = queue
        self.unlock = unlock
        self.dryrun = dryrun
        self.rerunincomplete = rerunincomplete
        self.useconda = useconda
        self.conda_frontend = 'mamba'
        self.usesingularity = usesingularity
        self.singularityargs = singularityargs
        self.restarttimes = restarttimes
        self.latency = latency_wait
        # Generate pipeline audit trail
        self.path_to_audit.mkdir(parents = True, exist_ok = True)
        self.audit_trail = self.generate_audit_trail()
        # Run pipeline
        self.run_snakemake()


    def is_not_repo(self):
        git_remote = subprocess.check_output(["git","remote", "-v"])
        return git_remote == b''
        
    def get_git_audit(self, git_file):
        print_message("Collecting information about the Git repository of this pipeline (see {})".format(git_file))
        if not self.is_not_repo():
            git_audit = {"repo": subprocess.check_output(["git","config", "--get", "remote.origin.url"]).strip(),
                        "commit": subprocess.check_output(["git","log", "-n", "1" "--pretty=format:'%H'"]).strip()}
        else:
            git_audit = {"repo": "NA (folder is not a repo or was not cloned through git command)",
                        "commit": "NA (folder is not a repo or was not cloned through git command)"}

        with open(git_file, 'w') as file:
            yaml.dump(git_audit, file, default_flow_style=False)

    def get_pipeline_audit(self, pipeline_file):
        print_message("Collecting information about the pipeline (see {})".format(pipeline_file))
        pipeline_info = {'pipeline_name': self.pipeline_name,
                        'pipeline_version': self.pipeline_version,
                        'timestamp': datetime.now().strftime('%d-%m-%Y %H:%M:%S'),
                        'hostname': str(subprocess.check_output(['hostname']).strip())}
        with open(pipeline_file, 'w') as file:
            yaml.dump(pipeline_info, file, default_flow_style=False)

    def get_conda_audit(self, conda_file):
        print_message("Getting information of the master environment used for this pipeline.")
        conda_audit = subprocess.check_output(["conda","list"]).strip()
        with open(conda_file, 'w') as file:
            file.writelines("Master environment list:\n\n")
            file.write(str(conda_audit))

    def generate_audit_trail(self):
        assert pathlib.Path(self.sample_sheet).exists(), "The sample sheet ({}) does not exist. This sample sheet is generated before starting your pipeline. Something must have gone wrong while creating it.".format(str(sample_sheet))
        assert pathlib.Path(self.user_parameters).exists(), "The provided user_parameters ({}) was not created properly or was deleted before starting the pipeline".format(','.join(configfiles))

        git_file = self.path_to_audit.joinpath('log_git.yaml')
        self.get_git_audit(git_file)
        conda_file = self.path_to_audit.joinpath('log_conda.txt')
        self.get_conda_audit(conda_file)
        pipeline_file = self.path_to_audit.joinpath('log_pipeline.yaml')
        self.get_pipeline_audit(pipeline_file)
        config_file = self.path_to_audit.joinpath('log_parameters.yaml')
        samples_file = self.path_to_audit.joinpath('sample_sheet.yaml')
        audit_sample_sheet = subprocess.Popen(['cp', self.sample_sheet, samples_file])
        audit_sample_sheet.wait(10)
        audit_userparams = subprocess.Popen(['cp', self.user_parameters, config_file])
        audit_userparams.wait(10)
        return [git_file, conda_file, pipeline_file, config_file, samples_file]

    def run_snakemake(self):

        print_message("Running {} pipeline.".format(self.pipeline_name))

        if self.local:
            print_message("Jobs will run locally")
            drmaa = None
        else:
            print_message("Jobs will be sent to the cluster")
            pathlib.Path(str(self.output_dir)).joinpath('log', 'drmaa').mkdir(parents = True, exist_ok = True)
            drmaa = " -q %s \
                    -n {threads} \
                    -o %s/log/drmaa/{name}_{wildcards}_{jobid}.out \
                    -e %s/log/drmaa/{name}_{wildcards}_{jobid}.err \
                    -R \"span[hosts=1]\" \
                    -R \"rusage[mem={resources.mem_mb}]\" " % (str(self.queue), str(self.output_dir), str(self.output_dir))
        
        try:
            snakemake.snakemake(self.snakefile,
                            workdir=self.workdir,
                            configfiles=[self.user_parameters, self.fixed_parameters],
                            config={"sample_sheet": str(self.sample_sheet)},
                            cores=self.cores,
                            nodes=self.cores,
                            drmaa = drmaa,
                            jobname=self.pipeline_name + "_{name}.jobid{jobid}",
                            use_conda=self.useconda,
                            conda_frontend=self.conda_frontend,
                            use_singularity = self.usesingularity,
                            singularity_args = self.singularityargs,
                            keepgoing=True,
                            printshellcmds=True,
                            force_incomplete=self.rerunincomplete,
                            restart_times=self.restarttimes, 
                            latency_wait=self.latency,
                            unlock=self.unlock,
                            dryrun=self.dryrun)
        except:
            print_error("An error occured while running the {} pipeline.".format(self.pipeline_name))
            raise

        print_message("Finished running {} pipeline!".format(self.pipeline_name))
            