import argparse
import os
import pandas as pd


class cgMLSTMulitreport():
    """Multireport combining the results of individual samples analyzed with 
    Chewbbaca (for cgMLST)
    """
    
    def __init__(self,
                cgmlst_dir='output/cgmlst',
                multireport_file=None):
        self.cgmlst_dir = cgmlst_dir
        if multireport_file is None:
            self.multireport_file = os.path.join(self.cgmlst_dir, 'cgmlst_multireport.csv')
        else:
            self.multireport_file = multireport_file

        self.chewbbaca_results_per_sample = list(self.get_results_chewbbaca())
        print(self.chewbbaca_results_per_sample)
        self.bind_results_if_same_colnames()

    def get_results_chewbbaca(self):
        for root, dirs, files in os.walk(self.cgmlst_dir, topdown=False):
            for file_name in files:
                if file_name == "results_alleles.tsv":
                    yield os.path.join(root, file_name)

    def read_results_chewbbaca(self):
        for file_ in self.chewbbaca_results_per_sample:
            raw_results = pd.read_csv(file_, sep='\t')
            sample_name = raw_results.loc[0, ['FILE']].iat[0]
            sample_name = sample_name.replace('.fasta', '')
            raw_results = raw_results.set_index('FILE')
            yield sample_name, raw_results

    def bind_results_if_same_colnames(self):
        chewbbaca_results = {sample_name:chewbbaca_results for (sample_name,chewbbaca_results) in self.read_results_chewbbaca()}
        multireport = pd.concat(chewbbaca_results)
        multireport = multireport.reset_index(level=1, drop=True)
        multireport.to_csv(self.multireport_file)
        print(multireport)


if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser(description='Multireport cgMLST with Chewbbaca.')
    argument_parser.add_argument('-d', '--cgmlst_dir', type=str, 
                        default='output/cgmlst',
                        help='Path to directory containing result files for cgMLST (calculated with chewbbaca).')
    argument_parser.add_argument('-o', '--output_multireport_file', type=str, 
                        default=None,
                        help='Path to multireport file that will be produced by this script. It should have the extension .csv')
    args = argument_parser.parse_args()
    cgMLSTMulitreport(cgmlst_dir=args.cgmlst_dir, multireport_file=args.output_multireport_file)