import argparse
from Bio import SeqIO
import dask.bag as db
import hashlib
import pandas as pd
import pathlib
from sys import path

script_path = str(pathlib.Path(pathlib.Path(__file__).parent.absolute()))
path.insert(0, script_path)
from download_cgmlst_scheme import cgmlst_schemes

class Allele():
    '''
    Allele Class. An allele has a sequence (seq), id number (id_number), an id 
    (id), a hash (hash) and the full SeqIO record from which all the other 
    attributes can be extracted (fasta_record)
    '''
    def __init__(self,
                id = None,
                id_number = None,
                seq = None,
                hash = None,
                locus_name = None,
                fasta_record = None):
        self.seq = seq
        self.id_number = id_number
        self.hash = hash
        self.locus_name = locus_name
        self.fasta_record = fasta_record

    def update_seq(self, seq):
        '''Provide new sequence for an Allele instance'''
        self.seq = str(seq)

    def update_hash(self):
        '''
        Calculate the hash for the seq attribute of the current Allele instance
        '''
        if not isinstance(self.seq, str):
            f"The allele sequence should be a string."
        byte_seq = bytes(self.seq, 'utf-8')
        self.hash = hashlib.sha1(byte_seq).hexdigest()

    def extract_attr_from_fasta_record(self):
        '''
        Get the ID from the header of a fasta record. The fasta record should
        be a SeqRecord made by a SeqIO.parse generator
        '''
        if not isinstance(self.fasta_record, SeqIO.SeqRecord):
            f"The provided fasta_record is not a SeqRecord instance."
        self.seq = str(self.fasta_record.seq)
        self.id = self.fasta_record.name
        # The fasta files should have the allele ID in the form:
        # <locus>_<id_number> .
        record_name_split = self.id.split('_')
        self.id_number = int(record_name_split[-1])
        record_name_split.pop(-1)
        # In case the locus name has an _ , I have to split, remove id_number 
        # and rejoin again
        self.locus_name = '_'.join(record_name_split)
        self.update_hash()



class Locus():
    '''
    Class for a locus. The locus has a locus name (name), a path to the fasta
    file that has all the sequences for the different alleles for this tool and
    a dictionary of alleles. This alleles attribute is expected to have objects
    from the class Allele()
    '''
    def __init__(self,
                fasta_file = None,
                name = "",
                alleles = {}):
        self.name = str(name)
        if fasta_file is not None:
            self.validate_fasta(fasta_file)
            self.fasta_file = pathlib.Path(fasta_file)
        else:
            self.fasta_file = None
        self.alleles = alleles

    def update_fasta_file(self, fasta_file):
        self.validate_fasta(fasta_file)
        self.fasta_file = pathlib.Path(fasta_file)

    def validate_fasta(self, fasta_file):
        if fasta_file is not None:
            if not pathlib.Path(fasta_file).is_file():
                f"The provided fasta file {self.fasta_file} does not exist."
            with open(fasta_file, 'r') as file_:
                fasta_record = SeqIO.parse(file_, 'fasta')
                if not any(fasta_record):
                    raise ValueError(f'The provided fasta file {self.fasta_file} does not have the expected fasta format.')    

    def get_locus_name_from_fasta_file(self):
        with open(self.fasta_file, 'r') as file_:
            fasta_gen = SeqIO.parse(file_, 'fasta')
            fasta_record = next(fasta_gen)
            # The fasta files should have the allele ID in the form:
            # <locus>_<id_number> . In case the locus name has an _ , I have to 
            # split, remove id_number and rejoin again
            self.name = '_'.join(fasta_record.name.split('_')[:-1])

    def __get_id_obj_pairs_from_fasta_record(self, record):
        '''
        Get a tupple with two values: the ID number from a fasta record and 
        the allele_obj itself.
        '''
        allele_obj = Allele(fasta_record = record)
        allele_obj.extract_attr_from_fasta_record()
        return (allele_obj.id, allele_obj)

    def get_all_alleles_from_fasta(self):
        with open(self.fasta_file, 'r') as file_:
            parsed_fasta = SeqIO.parse(file_, 'fasta')
            list_alleles = [self.__get_id_obj_pairs_from_fasta_record(record) for record in parsed_fasta]
            self.alleles = dict(list_alleles)
        print(f'{len(self.alleles)} added for locus {self.name}!')

    def get_subset_alleles_from_fasta(self, id_numbers = []):
        '''
        Get only a subset of the alleles contained in the fasta file for
        this locus. The list of alleles to extract should be given as a list
        of integers (id_numbers).
        '''
        if len(self.alleles) == 0:
            self.get_all_alleles_from_fasta()
        self.get_locus_name_from_fasta_file()
        self.allele_subset = {}
        for id_num in id_numbers:
            allele_id = self.name + '_' + str(id_num)
            self.allele_subset[allele_id] = self.alleles[allele_id]

    def make_table_with_hashes(self, subset = False):
        '''
        Creates a file in the same path and with the same name than the 
        fasta_file (but with .csv extension) that contains the list of 
        alleles (by id_number) and the corresponding hash. If subset = True
        then the file will be produced for subset_alleles only
        '''
        csv_file = self.fasta_file.with_suffix('.csv')
        if subset:
            hash_tuple_list = [(x.id_number, x.hash) for x in self.allele_subset.values()]
        else:
            hash_tuple_list = [(x.id_number, x.hash) for x in self.alleles.values()]
        hash_table = pd.DataFrame(hash_tuple_list, columns = ['id_number', 'hash'])
        hash_table.to_csv(csv_file, index = False)
        print(f'Allele table can be found at {csv_file}')
        return hash_table



class CgmlstResult():
    '''
    Class holding the result of a cgMLST call using ChewBBACA. The attributes
    for this class include a file path to a results_alleles.tsv file created
    by Chewbbaca (chewbbaca_result_file), a db_dir with the database files.
    This means, a fasta file per allele and a csv file with the translation
    of allele id numbers to hashes.
    '''
    def __init__(self,
                chewbbaca_result_file,
                scheme_name,
                db_dir = "/mnt/db/juno/typing_db/cgmlst/prepared_schemes",
                threads = 2): 
        self.chewbbaca_result_file = pathlib.Path(chewbbaca_result_file)
        if not self.chewbbaca_result_file.exists(): 
            f"The provided file {chewbbaca_result_file} does not exist."
        self.read_chewbbaca_result_file()
        self.db_dir = self.get_db_dir(db_dir = db_dir, scheme_name = scheme_name)
        self.threads = int(threads)

    def update_chewbbaca_result_file(self, chewbbaca_result_file):
        self.chewbbaca_result_file = pathlib.Path(chewbbaca_result_file)
        self.read_chewbbaca_result_file()

    def read_chewbbaca_result_file(self):
        if not self.chewbbaca_result_file.exists(): 
            f'The provided file (chewbbaca_result_file) does not exist.'
        self.cgmlst_result = pd.read_csv(
            self.chewbbaca_result_file, 
            sep = "\t", index_col = "FILE",
            dtype=str
        )

    def __validate_db_dir(self, db_dir):   
        db_dir = pathlib.Path(db_dir) 
        if not db_dir.exists():
            f'The directory {db_dir} does not exist.'
        db_dir_content = [fasta_file for fasta_file in db_dir.glob('*.fasta')]
        if len(db_dir_content) > 0:
            return True
        else:
            return False

    def get_db_dir(self, db_dir, scheme_name):
        print(f'Finding directory with scheme...')
        db_dir = pathlib.Path(db_dir)
        db_dir_genus = db_dir.joinpath(scheme_name)
        db_dir_scheme_folder = db_dir.joinpath(scheme_name, 'scheme')
        if self.__validate_db_dir(db_dir):
            print(f'Using scheme in directory {db_dir}.')
            return db_dir
        elif self.__validate_db_dir(db_dir_genus):
            print(f'Using scheme in directory {db_dir_genus}.')
            return db_dir_genus
        elif self.__validate_db_dir(db_dir_scheme_folder):
            print(f'Using scheme in directory {db_dir_scheme_folder}.')
            return db_dir_scheme_folder
        else:
            raise ValueError(f'No fasta files were found in the provided db_directory ({db_dir}) or any of its expected subfolders ({db_dir_genus} or {db_dir_scheme_folder})')
    
    def __get_locus_fasta_file_path(self, locus_name):
        locus_name = locus_name.replace('.fasta', '')
        basename_locus_fasta_file = locus_name + '.fasta'
        fasta_file_path = self.db_dir.joinpath(basename_locus_fasta_file)
        if not fasta_file_path.exists():
            f"The fasta file for the locus {locus_name} does not exist or is not in the expected location ({fasta_file_path})"
        return locus_name, fasta_file_path

    def __remake_csv_hash_id_tbl(self, fasta_file):
        locus_obj = Locus(fasta_file = fasta_file)
        locus_obj.get_all_alleles_from_fasta()
        hash_id_tbl = locus_obj.make_table_with_hashes()
        return hash_id_tbl

    def __get_hash_from_hash_id_tbl(self, hash_id_tbl, id_number):
        allele_hash = hash_id_tbl.loc[hash_id_tbl['id_number'] == id_number, 'hash'].tolist()
        allele_hash = list(set(allele_hash))
        if len(allele_hash) > 1:
            f'The id number {id_number} did not deliver a unique hash: {allele_hash}'
        return allele_hash[0]

    def __gen_hash_per_id_num(self, hash_id_tbl, id_numbers, fasta_file):
        for id_num in id_numbers:
            try:
                id_num_int = int(id_num)
                try:
                    allele_hash = self.__get_hash_from_hash_id_tbl(hash_id_tbl=hash_id_tbl, id_number=id_num_int)
                except IndexError:
                    hash_id_tbl = self.__remake_csv_hash_id_tbl(fasta_file=fasta_file)
                    allele_hash = self.__get_hash_from_hash_id_tbl(hash_id_tbl=hash_id_tbl, id_number=id_num_int)
                yield allele_hash
            except ValueError:
                yield 'NA'

    def get_hashes_for_locus(self, locus, id_numbers = []):
        print(f'Getting hashes for {locus}...')
        fasta_file = self.__get_locus_fasta_file_path(locus_name = locus)[1]
        csv_translation = fasta_file.with_suffix('.csv')
        try:
            hash_id_tbl = pd.read_csv(csv_translation, dtype = {'id_number': int, 'hash': str})
        except FileNotFoundError:
            locus_obj = Locus(fasta_file=fasta_file)
            locus_obj.get_all_alleles_from_fasta()
            hash_id_tbl = locus_obj.make_table_with_hashes()
        subset_allele_hashes = [x for x in self.__gen_hash_per_id_num(hash_id_tbl, id_numbers, fasta_file)]
        return subset_allele_hashes

    def __gen_hash_per_locus(self, loci_list):
        for locus in loci_list:
            id_numbers = list(self.cgmlst_result[locus])
            allele_hashes = self.get_hashes_for_locus(locus=locus, id_numbers=id_numbers)
            yield (locus, allele_hashes)
            
    def get_hashed_cgmlst_report(self):
        print(f'Getting dataframe with hashes...')
        loci_list = list(self.cgmlst_result.columns)
        allele_hashes = dict(self.__gen_hash_per_locus(loci_list))
        allele_hashes_tbl = pd.DataFrame.from_dict(allele_hashes)
        allele_hashes_tbl.index = self.cgmlst_result.index
        return allele_hashes_tbl


def main(input,
        scheme_name,
        db_dir,
        output,
        threads):
    print(f'Initializing...')
    chewbbaca_res = CgmlstResult(chewbbaca_result_file = input,
                                scheme_name = scheme_name,
                                db_dir = db_dir,
                                threads = threads)
    hashed_report = chewbbaca_res.get_hashed_cgmlst_report()
    print(f'Writing report to csv file {output}.')
    hashed_report.to_csv(path_or_buf=output)


if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("input",
                        type=pathlib.Path, 
                        help="Csv file containing the results of ChewBBACA.")
    argument_parser.add_argument('-s','--scheme-name', metavar='N', 
                        type=str, 
                        choices=[scheme for scheme in cgmlst_schemes.keys() if not scheme.startswith('test_')],
                        required = True,
                        help='Genus corresponding to the cgMLST scheme used. The supported schemes are: %(choices)s.')
    argument_parser.add_argument("-d", "--db-dir", 
                        type=pathlib.Path, 
                        help="Path to the database directory containing the prepared schemes for ChewBBACA.",
                        default = "/mnt/db/juno/typing_db/cgmlst/prepared_schemes")
    argument_parser.add_argument("-o", "--output",
                        type=pathlib.Path, 
                        help="File name (including path and with .csv extension) for the resulting report containing hashes instead of allele numbers.")
    argument_parser.add_argument('-t', '--threads', type=int, 
                        default=4,
                        help='Number of threads to be used.')
    args = argument_parser.parse_args()
    main(input = args.input,
        scheme_name = args.scheme_name,
        db_dir = args.db_dir,
        output = args.output,
        threads = args.threads)
