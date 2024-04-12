#!/usr/bin/env python3

import argparse
import logging
from xml.etree import ElementTree as ET

import pandas as pd


def parse_xml(input_file):
    """
    Parse blast XML output and return a list of dictionaries with the following keys:
    gene: query sequence name
    reference_length: length of the query sequence
    contig_name: name of the hit sequence
    hit_start: start position of the hit sequence
    hit_end: end position of the hit sequence
    matches: number of matches
    aln_length: length of the alignment
    gaps_in_hit: number of gaps in the hit sequence
    ref_start: start position of the query sequence
    ref_end: end position of the query sequence
    hit_sequence: hit sequence

    Parameters
    ----------
    input_file : str
        Path to the XML file

    Returns
    -------
    list
        List of dictionaries with the parsed data

    """
    output = []
    tree = ET.parse(input_file)
    root = tree.getroot()
    for child_1 in root.iter("Iteration"):
        logging.info(f"{child_1.tag=}, {child_1.attrib=}")
        iteration_query_def = child_1.find("Iteration_query-def").text
        iteration_query_len = child_1.find("Iteration_query-len").text
        for child_2 in child_1.iter("Hit"):
            logging.info(f"{child_2.tag=}, {child_2.attrib=}")
            hit_def = child_2.find("Hit_def").text
            for child_3 in child_2.iter("Hit_hsps"):
                logging.info(f"{child_3.tag=}, {child_3.attrib=}")
                for child_4 in child_3.iter("Hsp"):
                    logging.info(f"{child_4.tag=}, {child_4.attrib=}")
                    hsp_hit_from = child_4.find("Hsp_hit-from").text
                    hsp_hit_to = child_4.find("Hsp_hit-to").text
                    hsp_identity = child_4.find("Hsp_identity").text
                    hsp_query_from = child_4.find("Hsp_query-from").text
                    hsp_query_to = child_4.find("Hsp_query-to").text
                    hsp_gaps = child_4.find("Hsp_gaps").text
                    hsp_hseq = child_4.find("Hsp_hseq").text
                    hsp_align_len = child_4.find("Hsp_align-len").text
                    output.append(
                        {
                            "gene": iteration_query_def,
                            "reference_length": iteration_query_len,
                            "contig_name": hit_def,
                            "hit_start": hsp_hit_from,
                            "hit_end": hsp_hit_to,
                            "matches": hsp_identity,
                            "aln_length": hsp_align_len,
                            "gaps_in_hit": hsp_gaps,
                            "ref_start": hsp_query_from,
                            "ref_end": hsp_query_to,
                            "hit_sequence": hsp_hseq,
                        }
                    )

    return output


def convert_to_df(list_of_dicts):
    """
    Convert a list of dictionaries to a pandas DataFrame and calculate the percentage of coverage and identity

    Parameters
    ----------
    list_of_dicts : list
        List of dictionaries with the parsed data

    Returns
    -------
    pandas.DataFrame
        DataFrame with the parsed data

    """
    df = pd.DataFrame(list_of_dicts)
    types_dict = {
        "gene": str,
        "reference_length": int,
        "contig_name": str,
        "hit_start": int,
        "hit_end": int,
        "matches": int,
        "aln_length": int,
        "gaps_in_hit": int,
        "ref_start": int,
        "ref_end": int,
        "hit_sequence": str,
    }
    df = df.astype(types_dict)
    df["pct_coverage"] = df["aln_length"] / df["reference_length"]
    df["pct_identity"] = df["matches"] / df["aln_length"]
    df.to_csv("test.csv", index=False)
    return df


def filter(df, mincov, minid):
    """
    Filter the DataFrame based on the minimum coverage and identity

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with the parsed data
    mincov : float
        Minimum coverage, between 0 and 1
    minid : float
        Minimum identity, between 0 and 1

    Returns
    -------
    pandas.DataFrame
        Filtered DataFrame
    """
    mincov_filter = df["pct_coverage"] >= mincov
    minid_filter = df["pct_identity"] >= minid
    return df[mincov_filter & minid_filter]


def get_limits(start, end):
    """
    Get the left and right limits of an interval.

    Blast can output antisense hits, where start > end. This function returns the limits of the interval

    Parameters
    ----------
    start : int
        Start position of the interval
    end : int
        End position of the interval

    Returns
    -------
    tuple
        Tuple with the left and right limits of the interval
    """
    if start < end:
        return start, end
    else:
        return end, start


def check_overlap(start1, end1, start2, end2, overlap_threshold):
    """
    Check if two intervals overlap

    Parameters
    ----------
    start1 : int
        Start position of the first interval
    end1 : int
        End position of the first interval
    start2 : int
        Start position of the second interval
    end2 : int
        End position of the second interval
    overlap_threshold : float
        Minimum overlap threshold, between 0 and 1

    Returns
    -------
    bool
        True if the intervals overlap, False otherwise

    Notes
    -----
    Returns True if the overlap is greater than the overlap_threshold times the minimum length of the two intervals

    Might need some extra checks if one interval is very small.

    """
    # instead of dealing with frames and signs, just get left and right limits of intervals
    left_limit1, right_limit1 = get_limits(start1, end1)
    left_limit2, right_limit2 = get_limits(start2, end2)
    overlap = max(0, min(right_limit1, right_limit2) - max(left_limit1, left_limit2))
    length_1 = right_limit1 - left_limit1
    length_2 = right_limit2 - left_limit2
    return overlap >= overlap_threshold * min(length_1, length_2)


def get_overlapping_groups(df, overlap_threshold):
    """
    Get groups of overlapping hits

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with the parsed data
    overlap_threshold : float
        Minimum overlap threshold, between 0 and 1. Passed to check_overlap function

    Returns
    -------
    list
        List of lists with the overlapping groups

    """
    groups = []
    for i, row1 in df.iterrows():
        group = [row1["gene"]]
        # All vs all comparison on same contig, needs filtering of redundant comparisons
        tmp_df = df[df["contig_name"] == row1["contig_name"]]
        for j, row2 in tmp_df.iterrows():
            logging.info(f"checking {row1['gene']} and {row2['gene']}")
            logging.info(
                f"{row1['hit_start']}, {row1['hit_end']}, {row2['hit_start']}, {row2['hit_end']}"
            )
            overlap_bool = check_overlap(
                row1["hit_start"],
                row1["hit_end"],
                row2["hit_start"],
                row2["hit_end"],
                overlap_threshold,
            )
            if overlap_bool:
                logging.info("overlap found")
                group.append(row2["gene"])
        # Convert to set and then sorted list to avoid duplicates within group (hits to self)
        sorted_group = sorted(set(group))
        # Avoid duplicates between groups (same genes in different groups)
        if sorted_group not in groups:
            groups.append(sorted_group)
    return groups


def select_best_entries(df, groups):
    """
    Select the best entry for each group

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with the parsed data
    groups : list
        List of lists with the overlapping groups

    Returns
    -------
    dict
        Dictionary with the selected entries, where the key is the gene name

    Notes
    -----
    The best entry is the one with the highest percentage of identity. If this is the same, the one with the percentage of coverage closest to 1 is selected.

    """
    selected_entries = {}
    for group in groups:
        highest_score = 0
        selected_entry = None
        for gene in group:
            entry = df[df["gene"] == gene].iloc[0]
            # keep track of highest score and select the one with highest identity
            if entry["pct_identity"] > highest_score:
                highest_score = entry["pct_identity"]
                selected_entry = entry
            # if identity is the same, select the one with highest coverage
            # highest_score remains the same
            elif entry["pct_identity"] == highest_score:
                # compare absolute difference to 1 for coverage (higher cov is not always better)
                if abs(1 - entry["pct_coverage"]) < abs(
                    1 - selected_entry["pct_coverage"]
                ):
                    selected_entry = entry
        selected_entries[selected_entry["gene"]] = selected_entry
    return selected_entries


def convert_sort_and_save(selected_entries, output_path):
    """
    Convert the selected entries to a DataFrame and sort it

    Parameters
    ----------
    selected_entries : dict
        Dictionary with the selected entries, where the key is the gene name
    output_path : str
        Path to the output file

    Returns
    -------
    pandas.DataFrame
        DataFrame with the selected entries, sorted by gene name

    """
    df = pd.DataFrame(selected_entries.values())
    df["pct_identity"] = df["pct_identity"] * 100
    df["pct_coverage"] = df["pct_coverage"] * 100
    df = df[
        [
            "gene",
            "pct_identity",
            "pct_coverage",
            "reference_length",
            "contig_name",
            "hit_start",
            "hit_end",
            "ref_start",
            "ref_end",
            "matches",
            "aln_length",
            "gaps_in_hit",
            "hit_sequence",
        ]
    ]
    df.to_csv(output_path, index=False, float_format="%.2f")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="input file")
    parser.add_argument("output", help="output file")
    parser.add_argument("--mincov", type=float, default=0.6, help="minimum coverage")
    parser.add_argument("--minid", type=float, default=0.8, help="minimum identity")
    parser.add_argument(
        "--overlap_threshold", type=float, default=0.5, help="overlap threshold"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="increase output verbosity"
    )
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    parsed_output = parse_xml(args.input)
    df = convert_to_df(parsed_output)
    filtered_df = filter(df, args.mincov, args.minid)
    overlapping_groups = get_overlapping_groups(filtered_df, args.overlap_threshold)
    selected_entries = select_best_entries(filtered_df, overlapping_groups)
    convert_sort_and_save(selected_entries, args.output)


if __name__ == "__main__":
    main()
