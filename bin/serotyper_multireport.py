import argparse
from itertools import chain
import pathlib
from typing_extensions import final
import pandas as pd
import re
from warnings import warn


class SerotyperMultireport:
    """Class that will choose which serotyper multireport to make according
    to the input data
    """

    def __init__(self, input_files, output_file):
        self.input_files = [pathlib.Path(file_) for file_ in input_files]
        assert all(
            [file_.is_file() for file_ in self.input_files]
        ), "One or more of the specified input files do not exist!"
        self.output_file = pathlib.Path(output_file)
        self.sample_names = [file_.parent.name for file_ in self.input_files]

    def write_multireport(self):
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        self.multireport.to_csv(self.output_file, index=False)
        print(f"Serotyper multireport can be found at {self.output_file}.")
        return True


class SeqSero2Multireport(SerotyperMultireport):
    """Class combining results from multiple SeqSero2 serotype results (used
    for Salmonella)
    """

    def extract_from_seqsero2_result(self, input):
        seqsero_report = pd.read_csv(input, sep="\t")
        del seqsero_report["Output directory"]
        del seqsero_report["Input files"]
        return seqsero_report

    def make_multireport(self):
        multireport = [
            self.extract_from_seqsero2_result(seqsero2_result)
            for seqsero2_result in self.input_files
        ]
        multireport = pd.concat(multireport)
        # In file names, remove everything after the first underscore
        multireport["Sample name"] = self.sample_names
        self.multireport = multireport


class SerotypeFinderMultireport(SerotyperMultireport):
    """Class combining results from multiple SerotypeFinder results (used for
    E. coli)
    """

    def find_allele(self, serotype_df, allele_string):
        allele = []
        for column_name in serotype_df.columns.tolist():
            result_allele = serotype_df.loc["serotype", column_name]
            if allele_string in column_name:
                if result_allele not in allele:
                    allele.append(result_allele)
        return "/".join(allele)

    def get_sample_serotype(self, serotype_csv):
        result_csv = pd.read_csv(serotype_csv, index_col=0)
        genes = ["wzx", "wzy", "wzm", "wzt", "fl"]
        result_summary = []
        for gene in genes:
            result_summary.append(self.find_allele(result_csv, gene))
        return result_summary

    def report_o_type(self, row_df):
        row_df = row_df[["wzx", "wzy", "wzm", "wzt"]]
        # Split string in case any locus already consists of two alleles
        reported_alleles = [str(item).split("/") for item in row_df if item != ""]
        if len(reported_alleles) == 0:
            return "Error! No O type found"
        else:
            # Keep only unique alleles and report them all together (separated by a /)
            final_o_type = list(chain(*reported_alleles))
            # By request of Kim van der Zwaluw, only the numbers of the
            # serotype are reported not the letters after it (e.g. O128ac
            # reported as O128)
            final_o_type = [re.findall("\d+", sample) for sample in final_o_type]
            final_o_type = list(set(chain(*final_o_type)))
            # If there is more than one serotype, sort them from smaller to
            # larger (requested by Kim)
            if len(final_o_type) > 1:
                ranking_o_types = sorted(
                    range(len(final_o_type)), key=lambda k: int(final_o_type[k])
                )
                final_o_type = [final_o_type[ind] for ind in ranking_o_types]
            # Give back the serotype name
            final_o_type = [
                "O{}".format(final_o_type[i]) for i in range(0, len(final_o_type))
            ]
            return "/".join(final_o_type)

    def make_multireport(self):
        results = [self.get_sample_serotype(file_) for file_ in self.input_files]
        results_df = pd.DataFrame(results)
        results_df.columns = ["wzx", "wzy", "wzm", "wzt", "fli"]
        results_df.index = [self.sample_names]
        results_df["O type"] = [
            self.report_o_type(row) for index, row in results_df.iterrows()
        ]
        results_df["H type"] = results_df["fli"]
        results_df.loc[results_df["H type"] == "", "H type"] = "Error! No H type found"
        results_df.reset_index(inplace=True)
        self.multireport = results_df


class SerobaMultireport(SerotyperMultireport):
    """Class combining results from multiple Seroba results (used for
    S. pneumoniae)
    """

    def make_multireport(self):
        names = ["Sample", "Serotype", "Contamination"]
        results = [
            pd.read_csv(file, sep="\t", header=None, names=names)
            for file in self.input_files
        ]
        results_df = results[0].append(results[1:])
        results_df.to_csv(self.output_file, index=False)
        self.multireport = results_df


class ShigatyperMultireport(SerotyperMultireport):
    def make_multireport(self):
        dflist_shigatyper = []
        dflist_command = []
        # We combined ecoli and shigatyper, so now this code wants to run for every output file
        # it only needs to run for the shigella output files
        # this code needs to be changed before the pipeline can function again
        for outfile in self.input_files:
            dirname_splitted = str(outfile).split("/")

            if "command" in str(outfile):
                df = pd.read_csv(outfile, delimiter="\t")
                df.drop(columns=["sample"], inplace=True)
                df.insert(0, "Samplename", dirname_splitted[-2])
                dflist_command.append(df)

            if "shigatyper" in str(outfile):
                df = pd.read_csv(outfile, delimiter=",")
                if df.shape[0] == 0:
                    df.loc[len(df)] = "-"
                df.insert(0, "Samplename", dirname_splitted[-2])
                dflist_shigatyper.append(df)

        final_df_shigatyper = pd.concat(dflist_shigatyper, axis=0, ignore_index=True)
        final_df_command = pd.concat(dflist_command, axis=0, ignore_index=True)
        results_df = pd.merge(final_df_shigatyper, final_df_command, on="Samplename")
        results_df.to_csv(self.output_file, mode="a", index=False)
        self.multireport = results_df


class NeisseriaMultireport(SerotyperMultireport):
    def make_multireport(self):
        df_list_neisseria = []
        for outfile in self.input_files:
            # Query column is read as string, this way leading zeros will not be removed from the samplenames
            df = pd.read_csv(outfile, delimiter="\t", dtype={"Query": "string"})
            df_list_neisseria.append(df)

        final_df_neisseria = pd.concat(df_list_neisseria, axis=0, ignore_index=True)
        self.multireport = final_df_neisseria


class ChooseMultireport:
    """Depending on the input file(s) one or more serotype multireport(s) are
    created
    """

    def __init__(self, serotyper_result_files, out_dir):
        self.serotyper_result_files = serotyper_result_files
        self.out_dir = pathlib.Path(out_dir)
        self.make_serotyper_multireports()

    def __classify_serotyper_result_files(self):
        input_files = {
            "seqsero2": [],
            "serotypefinder": [],
            "seroba": [],
            "shigatyper": [],
            "neisseriatyper": [],
        }
        for file_ in self.serotyper_result_files:
            if file_.endswith("SeqSero_result.tsv"):
                input_files["seqsero2"].append(file_)
            elif file_.endswith("result_serotype.csv"):
                input_files["serotypefinder"].append(file_)
            elif file_.endswith("pred.tsv"):
                input_files["seroba"].append(file_)
            elif file_.endswith("command.txt") or file_.endswith("shigatyper.csv"):
                input_files["shigatyper"].append(file_)
            elif file_.endswith("neisseriatyper.tab"):
                input_files["neisseriatyper"].append(file_)
            else:
                raise ValueError(
                    f"The file {file_} is not recognized as a result from a supported serotyper."
                )
        input_files = {k: v for k, v in input_files.items() if len(v) > 0}
        number_of_expected_multireports = len(input_files)
        if number_of_expected_multireports > 1:
            warn(
                "More than one serotyper was run in the input samples and therefore multiple serotypers will be created."
            )
            multireport_files = [
                self.out_dir.joinpath(f"serotyper_multireport{n}.csv")
                for n in range(0, number_of_expected_multireports)
            ]
            multireport_files[0] = self.out_dir.joinpath("serotyper_multireport.csv")
        else:
            multireport_files = [self.out_dir.joinpath("serotyper_multireport.csv")]
        self.input_files = input_files
        self.multireport_files = multireport_files

    def make_serotyper_multireports(self):
        self.__classify_serotyper_result_files()
        # The name of all the multireports to make are stored in a list.
        # Every time one is made, it is 'popped' from the list and the next one
        # is produced. That is why the output_file argument is always output_file[0]
        output_file = self.multireport_files
        for serotyper_tool in self.input_files:
            print(f"Making serotyper multireport for {serotyper_tool}...")
            if serotyper_tool == "seqsero2":
                multireport = SeqSero2Multireport(
                    input_files=self.input_files[serotyper_tool],
                    output_file=output_file[0],
                )
            elif serotyper_tool == "serotypefinder":
                multireport = SerotypeFinderMultireport(
                    input_files=self.input_files[serotyper_tool],
                    output_file=output_file[0],
                )
            elif serotyper_tool == "shigatyper":
                multireport = ShigatyperMultireport(
                    input_files=self.input_files[serotyper_tool],
                    output_file=output_file[0],
                )
            elif serotyper_tool == "neisseriatyper":
                multireport = NeisseriaMultireport(
                    input_files=self.input_files[serotyper_tool],
                    output_file=output_file[0],
                )
            else:
                multireport = SerobaMultireport(
                    input_files=self.input_files[serotyper_tool],
                    output_file=output_file[0],
                )
            multireport.make_multireport()
            multireport.write_multireport()
            output_file.pop(0)
        assert all([file_.exists() for file_ in self.multireport_files])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        help="Input files (tsv) produced by the serotyper that need to be combined into one multireport",
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        type=str,
        default="output/serotype",
        help="Output directory where the serotyper multireport will be saved.",
    )
    args = parser.parse_args()
    ChooseMultireport(args.input, args.out_dir)
