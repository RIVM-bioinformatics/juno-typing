import argparse
import os
import pathlib
import pandas as pd
import re



def suspected_typhimurium(seqsero2output):
    # Serotypes to be considered typhimurium
    serotypes_typhimurium = ["4:i:1,2","I 4,[5],12:i:-"]
    with open(seqsero2output) as results:
        serotypedata = results.read()
        if any(x in serotypedata for x in serotypes_typhimurium):
            return True
        else:
            return False

def generate_depth_file(amplicon_file, 
                        forward_read, 
                        reverse_read, 
                        sample_name, 
                        threads = 1):
    threads = str(threads)
    ampliconname = str(pathlib.Path(amplicon_file).stem)
    os.system("bwa index %s"%(amplicon_file))
    os.system("bwa mem -k 17 -t %s %s %s %s > %s.sam"%(threads,amplicon_file,forward,reverse,sample_name))
    os.system("samtools sort -@ %s -n -O sam %s.sam | samtools fixmate -m -O bam - %s.bam"%(threads,sample_name,sample_name))
    os.system("rm %s.sam"%(sample_name))
    os.system("samtools sort -@ %s -O bam -o %s.sorted.bam %s.bam"%(threads,sample_name,sample_name))
    os.system("rm %s.bam"%(sample_name))
    os.system("samtools depth %s.sorted.bam | gzip > %s.%s.depth.txt.gz"%(sample_name,sample_name,ampliconname))
    os.system("rm %s.sorted.bam"%(sample_name))
    os.system("gunzip %s.%s.depth.txt.gz"%(sample_name,ampliconname))


def determine_allele_presence(depth_file,
                                ampliconname,
                                first_position,
                                last_position):
    with open(depth_file) as file:
        count = 0
        # loop over every line in the depth file
        for read in file:
            # split the line based on tabs
            splitread = read.split('\t')
            # remove enter (\n)
            splitread[2].replace('\n', '')
            # check the positions where the amplicon is supposed to be
            if int(splitread[1]) > first_position and int(splitread[1]) < last_position:
                count += int(splitread[2])
    # cutoff threshold for amount of reads between these positions is set to default 5000
    if count > threshold:
        return True, count
    else:
        return False, count


    
# Create and save multi_report
def monophasic_or_biphasic(seqsero2output,
                            forwardreads,
                            reversereads,
                            sample_name,
                            threads,
                            threshstyphi,
                            threshbiphasic):
    current_dir = pathlib.Path(__file__).parent.absolute()
    amplicons = {"FFLIB_FFLIA": {"file": current_dir.joinpath("amplicon/FFLIB_FFLIA.fasta"),
                            "threshold": threshstyphi, 
                            "first_position": 425,
                            "last_position": 1130},
                "sense_59_antisense_83": {"file": current_dir.joinpath("amplicon/sense_59_antisense_83.fasta"),
                                            "threshold": threshbiphasic, 
                                            "first_position": 425,
                                            "last_position": 1130}}	
    if suspected_typhimurium(seqsero2output):
        sense_59_antisense_83output, countsense_59_antisense_83 = generate_depth_file(amplicons['FFLIB_FFLIA']['file'], 
                                                                                        forwardreads, 
                                                                                        reversereads, 
                                                                                        sample_name, 
                                                                                        threads)
        FFLIB_FFLIAoutput, countFFLIB_FFLIA = generate_depth_file(amplicons['sense_59_antisense_83']['file'], 
                                                                    forwardreads, 
                                                                    reversereads, 
                                                                    sample_name, 
                                                                    threads)
        if FFLIB_FFLIAoutput == True and sense_59_antisense_83output == True:
            variant = "Biphasic"
        elif FFLIB_FFLIAoutput == True and sense_59_antisense_83output == False:
            variant = "Monophasic"
        else:
            variant = "N/A"
    else:
        variant = "N/A"
        countFFLIB_FFLIA = "N/A"
        countsense_59_antisense_83 = "N/A"
    return variant, countFFLIB_FFLIA, countsense_59_antisense_83



def main(args):
    variant, countFFLIB_FFLIA, countsense_59_antisense_83 = monophasic_or_biphasic(args.seqsero2output,
                                                                                    args.forwardreads,
                                                                                    args.reversereads,
                                                                                    args.sample_name,
                                                                                    args.threads,
                                                                                    args.threshbiphasic,
                                                                                    args.threshstyphi)
    #add the output into a dataframe using pandas
    seqsero_report = pd.read_csv(args.seqsero2output, sep='\t')
    seqsero_report = seqsero_report.iloc[:,[0,3,4,5,6,7,8,9,10]]
    seqsero_report.insert(7,"Monophasic/Biphasic", variant)
    seqsero_report.insert(10,"FFLIB_FFLIAcount", countFFLIB_FFLIA)
    seqsero_report.insert(11,"sense_59_antisense_83count", countsense_59_antisense_83)
    seqsero_report.to_csv(args.out_report,sep='\t', index=False, header=True)
    


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
    main(parser.parse_args())
