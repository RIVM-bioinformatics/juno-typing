import os
import pathlib
import argparse
import json
import sys
import re
import csv
import pandas as pd

def changetoavus(metas):
    avus=[]
    for meta in metas:
        avus.append({ "units": "", "attr": meta, "value": metas[meta] })
    avus.append({ "units": "", "attr": "pipeline", "value": "juno-typing" })
    avus.append({ "units": "", "attr": "project", "value": "salm" })
    return avus

def main(args):
    metadata=[]; values={}
    # walk through all files in output dir
    for filename in pathlib.Path(args.output_dir).glob("*/*/*"):
        f=str(filename).split("/")
        if f[len(f)-3] != "log" and f[len(f)-3] != "audit_trail":

           # for now we only analize SeqSero_result.tsv, add all values to the standard dict (analysis, samplename, filename)
           if f[len(f)-1] == "SeqSero_result.tsv":
              df = pd.read_csv(filename, sep='\t')
              values = {'analysis': f[len(f)-3], 'sample_name': f[len(f)-2], 'filename': f[len(f)-1]}
              for n in df.columns:
                  values.update({n: str(df[n][0])})
              # make irods avus from dict (only use attribute and value)
              avus = changetoavus(values)
           else:
              avus = changetoavus({ 'analysis': f[len(f)-3], 'sample_name': f[len(f)-2], 'filename': f[len(f)-1] })
           metadata.append( { 'path': str(os.path.relpath(filename, args.output_dir)), "type": "dataobject", 'metadata': avus })
    with open( os.path.join( args.output_dir,'metadata.json'), 'w') as outfile:
        json.dump( {"objects": metadata}, outfile, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Sample name to metadata script')
    parser.add_argument("output_dir")
    args = parser.parse_args()
    main(args)
