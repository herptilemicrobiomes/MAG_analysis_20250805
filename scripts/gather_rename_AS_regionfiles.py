#!/usr/bin/env python3

"""
Gather and rename AntiSMASH region files to facilitate clustering and 
linking back to dataset.
"""

import os
import re
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def parse_arguments():
    parser = argparse.ArgumentParser(description="Gather and rename AntiSMASH region files.")
    parser.add_argument("input_dir", help="Directory containing the AntiSMASH region files.")
    parser.add_argument("output_dir", help="Directory to save the renamed region files.")
    return parser.parse_args()

def gather_and_rename_region_files(input_dir, output_dir):
    """
    Gather and rename AntiSMASH region files to facilitate clustering and linking back to dataset.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    gbkoutdir = os.path.join(output_dir, "gbk")
    if not os.path.exists(gbkoutdir):
        os.makedirs(gbkoutdir)
    faoutdir = os.path.join(output_dir, "fasta")
    if not os.path.exists(faoutdir):
        os.makedirs(faoutdir)

    fnamematch = re.compile(r"^(\S+)\.(region\d+)\.gbk$")
    for UHM in os.listdir(input_dir):
        outUHMdir = os.path.join(gbkoutdir, UHM)
        if not os.path.exists(outUHMdir):
            os.makedirs(outUHMdir)
        outFAUHMdir = os.path.join(faoutdir, UHM)
        if not os.path.exists(outFAUHMdir):
            os.makedirs(outFAUHMdir)

        for MAG in os.listdir(os.path.join(input_dir, UHM)):
            MAG_path = os.path.join(input_dir, UHM, MAG)
            for filename in os.listdir(MAG_path):
                m = fnamematch.match(filename)
                if m:
                    # Extract the sample name from the filename
                    contig = m.group(1)
                    region = m.group(2)                    
                    new_filename = f"{MAG}.{contig}.{region}.gbk"
                    input_file_path = os.path.join(MAG_path, filename)
                    output_file_path = os.path.join(outUHMdir, new_filename)
                    output_file_path_fa = os.path.join(outFAUHMdir, new_filename.replace(".gbk", ".fasta"))
                    with open(input_file_path, "r") as input_handle:
                        inseqs = []
                        inseqs_fa = []
                        for record in SeqIO.parse(input_handle, "genbank"):
                            # Extract the region sequence
                            record.id = f"{MAG}.{contig}.{region}"
                            record.locus = record.id
                            #print(record.annotations)
                            record.annotations['accessions'] = [record.id]
                            #print(record.annotations)
                            #print(record.locus)
                            record.annotations['locus'] = record.id
                            record.description = record.id
                            newrecord = SeqIO.SeqRecord(
                                Seq(str(record.seq)),
                                id=record.id,
                                description=record.description,
                            )
                            newrecord.annotations = record.annotations.copy()
                            newrecord.dbxrefs = record.dbxrefs[:]
                            newrecord.features = record.features[:]
                            inseqs.append(newrecord)

                            for feature in record.features:
                                if feature.type == "CDS":
                                    # Extract the protein sequence
                                    protein_seq = feature.qualifiers.get("translation", [""])[0]
                                    locustag = feature.qualifiers.get("locus_tag", [""])[0]
                                    translationtable = feature.qualifiers.get("transl_table", [""])[0]
                                    if protein_seq:
                                        new_record = SeqIO.SeqRecord(
                                            Seq(protein_seq),
                                            id=f"{MAG}.{contig}.{region}.{locustag}",
                                            description=f"transl_table={translationtable}",
                                        )
                                        inseqs_fa.append(new_record)
                            # Write the region sequence to the output file
                            
                        SeqIO.write(inseqs, output_file_path, "genbank")
                        SeqIO.write(inseqs_fa, output_file_path_fa, "fasta")
                        


if __name__ == "__main__":
    args = parse_arguments()
    gather_and_rename_region_files(args.input_dir, args.output_dir)
    print(f"Gathered and renamed region files from {args.input_dir} to {args.output_dir}.")