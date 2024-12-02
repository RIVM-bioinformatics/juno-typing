import os
import pathlib
import argparse
import json
import sys
import re
import csv


def changetoavus(metas):
    avus=[]
    for meta in metas:
        #print(meta, metas[meta])
        avus.append({ "units": "", "attr": meta, "value": metas[meta] })
    #print(avus)
    avus.append({ "units": "", "attr": "pipeline", "value": "juno-typing" })
    avus.append({ "units": "", "attr": "project", "value": "salm" })
    return avus

def main(args):
    metadata=[]
    for filename in pathlib.Path(args.output_dir).glob("*/*/*"):
        f=str(filename).split("/")
        if f[len(f)-3] != "log" and f[len(f)-3] != "audit_trail":
           #print(f[len(f)-3], f)
           if f[len(f)-1] == "SeqSero_result.tsv":
              with open(filename) as file:
                  tsv_file = csv.reader(file, delimiter="\t")
                  for a in tsv_file:
                      if a[8] != "Predicted serotype":
                          avus = changetoavus({ 'analysis': f[len(f)-3], 'sample_name': f[len(f)-2], 'filename': f[len(f)-1], 'Predicted serotype': a[8] })
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
