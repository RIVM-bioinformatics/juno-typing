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
        
        self.unique_id = uuid4()
        self.__subdirs_ = self.find_input_subdirs(input_dir)
        self.check_input_dir(self.__subdirs_, input_type)
        self.sample_dict = self.make_sample_dict(self.__subdirs_, 
                                                input_type)

    def find_input_subdirs(self, input_dir):
        """Function to check whether the input is from the Juno annotation pipeline or just a simple input directory"""
        # Only when the input_dir comes from the Juno-annotation pipeline the fastq and fasta files
        # do not need to be in the same folder (fastq are then expected in a subfolder called <input_dir>/clean_fastq
        # and the )
        if input_dir.joinpath('clean_fastq').exists() and input_dir.joinpath('de_novo_assembly_filtered').exists():
            return {'fastq': input_dir.joinpath('clean_fastq'),
                    'fasta': input_dir.joinpath('de_novo_assembly_filtered')}
        else:
            return {'fastq': input_dir,
                    'fasta': input_dir}

    def check_contains_file_wextension(self, 
                                        input_dir, 
                                        extension = ('.fastq', '.fastq.gz', '.fq', '.fq.gz', '.fasta')):
        for item in input_dir.iterdir():
            if item.is_file():
                if str(item).endswith(extension):
                    return True
                    break
        
        raise ValueError('Input directory ({}) does not contain files that end with one of the expected extensions {}.'.format(input_dir, extension))

    def check_input_dir(self, 
                        input_subdirs, 
                        input_type):
        """Function to check that input directory is indeed an existing directory that contains files with the expected extension (fastq or fasta)"""
    
        if input_type == 'fastq':
            return self.check_contains_file_wextension(input_subdirs['fastq'], 
                                                        ('.fastq', '.fastq.gz', '.fq', '.fq.gz'))
        elif input_type == 'fasta':
            return self.check_contains_file_wextension(input_subdirs['fasta'], 
                                                        ('.fasta'))
        else:
            return self.check_contains_file_wextension(input_subdirs['fastq'], 
                                                        ('.fastq', '.fastq.gz', '.fq', '.fq.gz')) and  self.check_contains_file_wextension(input_subdirs['fasta'], 
                                                                                                                                            ('.fasta'))

    def enlist_fastq_samples(self, 
                            input_dir):
        """Function to enlist the fastq files found in the input directory. Returns a dictionary with the form {sample: fastq_file}"""

        pattern = re.compile("(.*?)(?:_S\d+_|_S\d+.|_|\.)(?:p)?R?(1|2)(?:_.*\.|\..*\.|\.)f(ast)?q(\.gz)?")
        samples = {}
        for file_ in input_dir.iterdir():
            if file_.is_dir():
                continue
            match = pattern.fullmatch(file_.name)
            if match:
                sample = samples.setdefault(match.group(1), {})
                sample["R{}".format(match.group(2))] = str(file_)               
        return samples

    def enlist_fasta_samples(self, 
                            input_dir):
        """Function to enlist the fasta files found in the input directory. Returns a dictionary with the form {sample: fasta_file}"""

        pattern = re.compile("(.*?).fasta")
        samples = {}
        for file_ in input_dir.iterdir():
            if file_.is_dir():
                continue
            match = pattern.fullmatch(file_.name)
            if match:
                sample = samples.setdefault(match.group(1), {})
                sample["assembly"] = str(file_)
        return samples


    def make_sample_dict(self, 
                        input_subdirs, 
                        input_type):
        """Function to make a sample sheet from the input directory (expecting either fastq or fasta files as input)"""

        if input_type == 'fastq':
            samples = self.enlist_fastq_samples(input_subdirs['fastq'])
        elif input_type == 'fasta':
            samples = self.enlist_fasta_samples(input_subdirs['fasta'])
        else:
            samples = self.enlist_fastq_samples(input_subdirs['fastq'])
            samples_fasta = self.enlist_fasta_samples(input_subdirs['fasta'])
            for k in samples.keys():
                samples[k]['assembly'] = samples_fasta[k]['assembly']
        return samples



class RunSnakemake:
    """Class with necessary input to run Snakemake"""

    def __init__(self,
                pipeline_name,
                pipeline_version,
                sample_dict,
                output_dir,
                workdir,
                sample_sheet=pathlib.Path('sample_sheet.yaml'), 
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
        self.output_dir = output_dir
        self.snakefile = snakefile
        self.workdir = workdir
        self.user_parameters = user_parameters
        self.fixed_parameters = fixed_parameters
        self.path_to_audit = self.output_dir.joinpath('audit_trail')
        self.path_to_audit.mkdir(parents = True, exist_ok = True)
        self.generate_audit_trail(self.path_to_audit,
                                    sample_dict,
                                    sample_sheet,
                                    pipeline_name,
                                    pipeline_version,
                                    user_parameters)
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
        self.run_snakemake(self.output_dir,
                            pipeline_name,
                            self.snakefile,
                            self.workdir,
                            self.user_parameters,
                            self.fixed_parameters,
                            self.cores,
                            self.local,
                            self.queue,
                            self.unlock,
                            self.rerunincomplete,
                            self.dryrun,
                            self.useconda,
                            self.conda_frontend,
                            self.usesingularity,
                            self.singularityargs,
                            self.restarttimes,
                            self.latency)


    def get_git_audit(self, git_file, sample_dict):
        print_message("Collecting information about the Git repository of this pipeline (see {})".format(git_file))
        git_audit = {"repo": subprocess.check_output(["git","config", "--get", "remote.origin.url"]).strip(),
                    "commit": subprocess.check_output(["git","log", "-n", "1" "--pretty=format:'%H'"]).strip()}
        with open(git_file, 'w') as file:
            yaml.dump(sample_dict, file, default_flow_style=False)

    def get_pipeline_audit(self, pipeline_file, pipeline_name, pipeline_version):
        print_message("Collecting information about the pipeline (see {})".format(pipeline_file))
        pipeline_info = {'pipeline_name': pipeline_name,
                        'pipeline_version': pipeline_version,
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

    def generate_audit_trail(self, 
                            path_to_audit,
                            sample_dict,
                            sample_sheet,
                            pipeline_name,
                            pipeline_version,
                            user_parameters):
        assert sample_sheet.exists(), "The sample sheet ({}) does not exist. This sample sheet is generated before starting your pipeline. Something must have gone wrong while creating it.".format(str(sample_sheet))
        assert user_parameters.exists(), "The provided user_parameters ({}) was not created properly or was deleted before starting the pipeline".format(','.join(configfiles))

        git_file = path_to_audit.joinpath('log_git.yaml')
        self.get_git_audit(git_file, sample_dict)
        conda_file = path_to_audit.joinpath('log_conda.txt')
        self.get_conda_audit(conda_file)
        pipeline_file = path_to_audit.joinpath('log_pipeline.yaml')
        self.get_pipeline_audit(pipeline_file, pipeline_name, pipeline_version)
        config_file = path_to_audit.joinpath('log_parameters.yaml')
        samples_file = path_to_audit.joinpath('sample_sheet.yaml')
        audit_sample_sheet = subprocess.Popen(['cp', sample_sheet, samples_file])
        audit_sample_sheet.wait(10)
        audit_userparams = subprocess.Popen(['cp', user_parameters, config_file])
        audit_userparams.wait(10)

    def run_snakemake(self,
                    output_dir,
                    pipeline_name,
                    snakefile,
                    workdir,
                    user_parameters,
                    fixed_parameters,
                    cores,
                    local,
                    queue,
                    unlock,
                    rerunincomplete,
                    dryrun,
                    useconda,
                    conda_frontend,
                    usesingularity,
                    singularityargs,
                    restarttimes,
                    latency):

        print_message("Running {} pipeline.".format(pipeline_name))

        if local:
            print_message("Jobs will run locally")
            drmaa = None
        else:
            print_message("Jobs will be sent to the cluster")
            drmaa = " -q %s \
                    -n {threads} \
                    -o %s/log/drmaa/{name}_{wildcards}_{jobid}.out \
                    -e %s/log/drmaa/{name}_{wildcards}_{jobid}.err \
                    -R \"span[hosts=1]\" \
                    -R \"rusage[mem={resources.mem_mb}]\" " % (str(queue), str(output_dir), str(output_dir))
        
        try:
            snakemake.snakemake(snakefile,
                            workdir=workdir,
                            configfiles=[user_parameters, fixed_parameters],
                            cores=cores,
                            nodes=cores,
                            drmaa = drmaa,
                            jobname=pipeline_name + "_{name}.jobid{jobid}",
                            use_conda=useconda,
                            conda_frontend=conda_frontend,
                            use_singularity = usesingularity,
                            singularity_args = singularityargs,
                            keepgoing=True,
                            printshellcmds=True,
                            force_incomplete=rerunincomplete,
                            restart_times=restarttimes, 
                            latency_wait=latency,
                            unlock=unlock,
                            dryrun=dryrun)
        except:
            print_error("An error occured while running the {} pipeline.".format(pipeline_name))
            raise
        print_message("Finished running pipeline!")
            