#!/usr/bin/env python3

import logging
from pathlib import Path
import pandas as pd


def df_to_dict(df, column_name):
    """
    Convert a dataframe to a dictionary

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe
    column_name : str
        Column name to use as key
    
    Returns
    -------
    dict
        Dictionary with column values as keys and the rest of the row as values
    """
    df = df[df["Column"] == column_name]
    dict_ = df.set_index("Value").to_dict(orient="index")
    return dict_


def add_context(df_context, value, col_name):
    """
    Add context to a value

    Parameters
    ----------
    df_context : pd.DataFrame
        Dataframe with context for e.g. specific serotypes
    value : str
        Value to check, e.g. serotype name
    col_name : str
        Column name to check, e.g. "Predicted serotype"

    Returns
    -------
    str
        Context for the value
    """
    logging.info(f"Checking context for {col_name}={value}")
    context = None
    dict_context = df_to_dict(df_context, col_name)
    if value in dict_context:
        logging.info(f"Found context for {col_name}={value}")
        context_partial = dict_context[value]["Context"]
        context = f"{col_name}={value}: {context_partial}"
    return context


def main(args):
    logging.info(f"Reading {args.input} and {args.context}")
    df = pd.read_csv(args.input, sep="\t")
    df_context = pd.read_csv(args.context, sep="\t")
    notes = []

    logging.info(f"Check if this is a single sample report")
    if df.shape[0] > 1:
        raise ValueError("This script only works for single sample reports")

    # Add context to O antigen
    O_gene = df["O antigen prediction"].values[0]
    notes.append(add_context(df_context, O_gene, "O antigen prediction"))

    # Add context to serotype
    serotype = df["Predicted serotype"].values[0]
    notes.append(add_context(df_context, serotype, "Predicted serotype"))

    # Combine all notes
    note_str = "|".join([note for note in notes if note is not None])
    df["Notes"] = note_str

    # Write to output
    logging.info(f"Writing to {args.output}")
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Add context to SeqSero report")

    parser.add_argument("-i", "--input", required=True, type=Path)
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument("-c", "--context", required=True, type=Path)
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")

    main(args)
