import argparse
import os
import pathlib
import pandas as pd
import re
import subprocess


class SalmonellaMonoOrBiphasic():

    def __init__(self,
                seqsero2output,
                forwardreads,
                reversereads,
                sample_name,
                out_report,
                threads=4,
                threshbiphasic=4000,
                threshstyphi=4000):
        self.seqsero2output = pathlib.Path(seqsero2output)
        assert self.seqsero2output.exists(), f'The provided file {seqsero2output} does not exist.'
        self.forwardreads = pathlib.Path(forwardreads)
        assert self.forwardreads.exists(), f'The provided file {forwardreads} does not exist.'
        self.reversereads = pathlib.Path(reversereads)
        assert self.reversereads.exists(), f'The provided file {reversereads} does not exist.'
        self.sample_name = sample_name
        self.out_report = pathlib.Path(out_report)
        self.out_dir = self.out_report.parent
        assert self.out_dir.exists(), f'The provided output directory {self.out_dir} does not exist. Please provide an existing directory.'
        self.threads = int(threads)

        # Serotypes to be considered typhimurium
        self.serotypes_typhimurium = ["4:i:1,2","I 4,[5],12:i:-"]
        self.pipeline_dir = pathlib.Path(__file__).parent.parent.absolute()
        self.amplicons = {"FFLIB_FFLIA": {'file': self.pipeline_dir.joinpath("files/salmonellamonophasic_amplicon/FFLIB_FFLIA.fasta"),
                                        'first_position': 75,
                                        'last_position': 780,
                                        'threshold': int(threshstyphi)},
                        "sense_59_antisense_83": {'file': self.pipeline_dir.joinpath("files/salmonellamonophasic_amplicon/sense_59_antisense_83.fasta"),
                                                    'first_position': 425,
                                                    'last_position': 1130,
                                                    'threshold': int(threshbiphasic)}}

    def suspected_typhimurium(self):
        with open(self.seqsero2output) as results:
            serotypedata = results.read()
            if any(value in serotypedata for value in self.serotypes_typhimurium):
                print(f'\nThe sample is suspected to be a Salmonella typhimurium. In silico PCR will be performed to confirm.')
                return True
            else:
                print(f'\nThe sample is not suspected to be a Salmonella typhimurium. In silico PCR will not be performed.')
                return False

    def __index_amplicon_files(self, amplicon):
        amplicon_file = self.amplicons[amplicon]['file']
        amplicon_index_files = [str(amplicon_file) + extension for extension in ['.amb', '.ann', '.bwt', '.pac', '.sa']]
        index_files_exist = [pathlib.Path(file_).exists() for file_ in amplicon_index_files]
        if not all(index_files_exist):
            subprocess.run(['bwa', 'index', amplicon_file], 
                            check = True, timeout=500)

    def __generate_depth_file(self, amplicon):
        print(f'Performing in-silico PCR...')
        threads = str(self.threads)
        self.__index_amplicon_files(amplicon=amplicon)
        sam_file = self.out_dir.joinpath(self.sample_name + '.sam')
        bam_file = str(self.out_dir.joinpath(self.sample_name + '.bam'))
        sorted_bam_file = str(self.out_dir.joinpath(self.sample_name + '_sorted.bam'))
        depth_file = str(self.out_dir.joinpath(self.sample_name + '_' + amplicon + '_depth.txt'))
        with open(sam_file, 'wb', 0) as outputfile:
            subprocess.check_call(['bwa', 'mem', '-k', '17', 
                                    '-t', threads, 
                                    str(self.amplicons[amplicon]['file']), 
                                    str(self.forwardreads), str(self.reversereads)],
                                    stdout=outputfile)
        sort_sam_by_name = subprocess.Popen(['samtools', 'sort', '-@', threads, 
                            '-n', '-O', 'sam', sam_file],
                            stdout=subprocess.PIPE)
        fixmate = subprocess.run(['samtools', 'fixmate', '-m', 
                            '-m', '-O', 'bam', '-', bam_file],
                            stdin=sort_sam_by_name.stdout)
        sort_sam_by_name.wait()
        subprocess.run(['samtools', 'sort', '-@', threads, 
                            '-O', 'bam', '-o', sorted_bam_file, 
                            bam_file],
                            check = True, timeout=500)
        with open(depth_file, 'wb') as outputfile:
            subprocess.check_call(['samtools', 'depth', sorted_bam_file],
                                    stdout=outputfile)
        pathlib.Path(sam_file).unlink()
        pathlib.Path(bam_file).unlink()
        pathlib.Path(sorted_bam_file).unlink()
        return depth_file


    def determine_allele_presence(self,
                                    depth_file,
                                    amplicon):
        first_position = self.amplicons[amplicon]['first_position']
        last_position = self.amplicons[amplicon]['last_position']
        threshold = self.amplicons[amplicon]['threshold']
        count = 0
        with open(depth_file) as file:
            for read in file:
                splitread = read.split('\t')
                splitread[2].replace('\n', '')
                if int(splitread[1]) > first_position and int(splitread[1]) < last_position:
                    count += int(splitread[2])
        if count > threshold:
            print(f'Positive for {amplicon}.\nRead count: {count}.')
            return True, count
        else:
            print(f'Negative for {amplicon}.\nRead count: {count}.')
            return False, count


        
    # Create and save multi_report
    def monophasic_or_biphasic(self):
        if self.suspected_typhimurium():
            depth_file_sense_59 = self.__generate_depth_file(amplicon='sense_59_antisense_83')
            depth_file_FFLIB = self.__generate_depth_file(amplicon='FFLIB_FFLIA')
            sense_59_antisense_83output, countsense_59_antisense_83 = self.determine_allele_presence(depth_file=depth_file_sense_59,
                                                                                                    amplicon='sense_59_antisense_83')
            FFLIB_FFLIAoutput, countFFLIB_FFLIA = self.determine_allele_presence(depth_file=depth_file_FFLIB,
                                                                                amplicon='FFLIB_FFLIA')
            if FFLIB_FFLIAoutput == True and sense_59_antisense_83output == True:
                variant = "Biphasic"
            elif FFLIB_FFLIAoutput == True and sense_59_antisense_83output == False:
                variant = "Monophasic"
            else:
                variant = "Not confirmed as Typhimurium"
            print(f'\n Sample determined as: {variant}\n')
        else:
            variant = None
            countFFLIB_FFLIA = None
            countsense_59_antisense_83 = None
        return variant, countFFLIB_FFLIA, countsense_59_antisense_83
    
    def write_typhimurium_plus_seqsero2_report(self):
        variant, countFFLIB_FFLIA, countsense_59_antisense_83 = self.monophasic_or_biphasic()                        
        #add the output into a dataframe using pandas
        seqsero_report = pd.read_csv(self.seqsero2output, sep='\t')
        seqsero_report = seqsero_report.iloc[:,[0,3,4,5,6,7,8,9,10]]
        seqsero_report.insert(7,"Monophasic/Biphasic", variant)
        seqsero_report.insert(10,"FFLIB_FFLIAcount", countFFLIB_FFLIA)
        seqsero_report.insert(11,"sense_59_antisense_83count", countsense_59_antisense_83)
        seqsero_report.to_csv(self.out_report, sep='\t', index=False, header=True)
        print(f'Serotype report can be found at: {self.out_report}')
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--sample_name',
                        type = str,
                        required = True,
                        metavar = "STR",
                        help = "Name of the sample (without extension) to be used for naming the sample within the report and intermediate files")
    parser.add_argument('-i','--seqsero2output', 
                        type = str,
                        required = True,
                        metavar = "FILE",
                        help="Csv file containing the results of SeqSero2 for the sample being processed.")
    parser.add_argument('-f','--forwardreads', 
                        type = str,
                        required = True,
                        metavar = "FILE",
                        help="Fastq file containing the forward (R1) reads for one sample.")
    parser.add_argument('-r','--reversereads', 
                        type = str,
                        required = True,
                        metavar = "FILE",
                        help="Fastq file containing the reverse (R2) reads for one sample.")
    parser.add_argument('-o', '--out_report', 
                        type = str,
                        required = True,
                        metavar = "FILE",
                        help="Path (relative or absolute) and name of the output file that contains the multireport")
    parser.add_argument('-s', '--threshstyphi', 
                        type = int,
                        default = 4000,
                        metavar = "INT", 
                        help="Threshold (number of reads) for considering a sample as S. typhimurium")
    parser.add_argument('-b', '--threshbiphasic', 
                        type = int,
                        default = 4000,
                        metavar = "INT",
                        help="Threshold (number of reads) for considering a sample as biphasic")
    parser.add_argument('-t', '--threads', 
                        default = 2,
                        type = int,
                        metavar = "INT",
                        help="Number of threads to be used for blast")
    args = parser.parse_args()
    mono_or_biphasic = SalmonellaMonoOrBiphasic(seqsero2output=args.seqsero2output,
                                                forwardreads=args.forwardreads,
                                                reversereads=args.reversereads,
                                                sample_name=args.sample_name,
                                                out_report=args.out_report,
                                                threads=args.threads,
                                                threshbiphasic=args.threshbiphasic,
                                                threshstyphi=args.threshstyphi)
    mono_or_biphasic.write_typhimurium_plus_seqsero2_report()
