#!/usr/bin/env python3

import argparse
import os
import re
import csv
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description="Add taxonomy to Bigscape record annotations.")
    parser.add_argument("-i","--input_record", 
                        help="Bigscape record annotation to be updated.",
                        default="record_annotations.tsv")
    parser.add_argument("--gtdbtk",nargs="+",
                        help="GTKDB taxonomy summary.",
                        default="gtdbtk.bac120.summary.tsv")

    parser.add_argument("-o","--output",
                        type=argparse.FileType('w'), 
                        default=sys.stdout,
                        help="Output the new record_annotations file. (default stdout)")

    parser.add_argument("-m","--metadata",nargs="+",
                        default=["lib/animal_metadata.csv","lib/wood_frog.csv"],
                        help="animal metadata for the UHM/Sample ID to host and animal info")
    
    return parser.parse_args()

def parse_gtdbfiles(taxonomy_files):
    table = dict()
    for taxfile in taxonomy_files:
        with open(taxfile, "rt") as intax:
            gtdbparser = csv.reader(intax,delimiter="\t")
            header = next(gtdbparser)
            for row in gtdbparser:
                MAG = row[0]
                classification = row[1]
                organism = ""
                if classification.startswith("Unclassified"):
                    organism = classification
                else:
                    classification_parsed = classification.split(";")
                    for node in reversed(classification_parsed):
                        (rank,taxon_name) = node.split("__")
                        if len(taxon_name):
                            organism = taxon_name
                            break
                table[MAG] = [organism,classification]
    return table

def parse_metadata(metadatafiles):
    table = dict()
    for metadatafile in metadatafiles:
        # print(f'processing file {metadatafile}', file=sys.stderr)
        with open(metadatafile, "rt") as infh:
            metacsv = csv.DictReader(infh,delimiter=",")
            for row in metacsv:
                sample_id = row['sample_id']
                table[sample_id] = row

    return table
def update_record_annotation(recordfile,taxondb,metadata,outfh):
    header = []
    metadata_cols = [ "host_taxon", "host_genus", "animal_ecomode", "Clade_Order", "Family", "Diet", "Habitat", "ecoregion_III"]
    with open(recordfile, "rt") as rf:
        recordparser = csv.reader(rf,delimiter="\t")
        header = next(recordparser)
    
    header.extend(metadata_cols)
    with open(recordfile, "rt") as rf:
        recordparser = csv.DictReader(rf,delimiter="\t")

        recordwriter = csv.DictWriter(outfh,delimiter="\t",
                                    lineterminator='\n', 
                                    fieldnames=header)
        recordwriter.writeheader()       
        for row in recordparser:
            name = row["Record"]
            m = re.match(r'^(\S+\.bin.\d+)\.',name)
            if m:
                name = m.group(1)
                if name in taxondb:
                    row["Organism"] = taxondb[name][0]
                    row["Taxonomy"] = taxondb[name][1]
            else:
                if not name.startswith("BGC"):
                    print(f'Cannot find MAG base name in {name}',file=sys.stderr)
            sampid = name.split('.')[0]
            if sampid in metadata:
                for meta in metadata_cols:
                    row[meta] = metadata[sampid][meta]
            else:
                print(f"sample {sampid} does not have metadata",file=sys.stderr)

            recordwriter.writerow(row)
    
if __name__ == "__main__":
    args = parse_arguments()
    gtdbtable = parse_gtdbfiles(args.gtdbtk)
    metatable = parse_metadata(args.metadata)
    update_record_annotation(args.input_record, gtdbtable, metatable, args.output)
    print(f"Parsed record {args.input_record} with {args.gtdbtk}.")
    
