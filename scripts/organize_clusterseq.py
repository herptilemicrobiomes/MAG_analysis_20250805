#!/usr/bin/env python3

import sys
import re
import os
import argparse
import csv

def parse_arguments():
    parser = argparse.ArgumentParser(description="Gather region protein files by cluster-type.")
    parser.add_argument("input_dir", help="Directory containing the AntiSMASH protein processed fasta files.")
    parser.add_argument("output_dir", help="Directory to save the renamed region files.")
    parser.add_argument("-a", "--annotations",
                        default="record_annotations.tsv",
                        help="annotations file renamed region files.")
    return parser.parse_args()

def gather_cluster_files(inputdir, annotations, outdir):
    """
    Process region files fasta and clusters.
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    annot_table = dict()
    with open(annotations, "rt") as infh:
        recordparser = csv.DictReader(infh,delimiter="\t")
        for row in recordparser:
            recordname = row["Record"]
            gbk = row["GBK"]
            clusterclass = row["Class"]
            category = row["Category"]
            if category not in annot_table:
                annot_table[category] = set()
            annot_table[category].add(gbk)
    
    for classtype in annot_table:
        print(f"There are {len(annot_table[classtype])} records for {classtype}",file=sys.stderr)
        outfile = os.path.join(outdir,f"{classtype}.clusters.fas")
        with open(outfile,"w") as outfh:
            n = 0
            for record in annot_table[classtype]:
                m = re.match(r'^(\S+)_R\.bin',record)
                MAG = ""
                if m:
                    MAG = m.group(1)
                else:
                    m = re.match(r'^(\S+)\.bin',record)
                    if m:
                        MAG = m.group(1)
                    else:
                        print(f"cannot determine MAG parent from {record}",file=sys.stderr)
                        continue
                fafile = os.path.join(inputdir,MAG,f"{record}.fasta")
                if not os.path.exists(fafile):
                    print(f"cannot find file {fafile} in type {classtype}",file=sys.stderr)
                with open(fafile,"rt") as fain:
                    for line in fain:
                        if line.startswith(">"):
                            n += 1
                        outfh.write(line)
            print(f"Wrote {n} sequences to Class {classtype} {outfile}",file=sys.stderr)
if __name__ == "__main__":
    args = parse_arguments()
    gather_cluster_files(args.input_dir, args.annotations, args.output_dir)