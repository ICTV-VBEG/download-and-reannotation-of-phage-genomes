#!/usr/bin/env python3
import os

from Bio import SeqIO


def split(fastafile="test_fasta.fasta", outfastadir="splitoutput"):
    """Extract multiple sequence fasta file and write each sequence in separate file"""
    os.system("mkdir -p %s" % (outfastadir))
    with open(fastafile) as FH:
        record = SeqIO.parse(FH, "fasta")
        file_count = 0
        for seq_rec in record:
            file_count = file_count + 1
            header = seq_rec.id
            with open("%s/%s.fasta" % (outfastadir, str(header)), "w") as FHO:
                SeqIO.write(seq_rec, FHO, "fasta")
    if file_count == 0:
        raise Exception("No valid sequence in fasta file")
    return "Done"


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f",
        "--fastafile",
        action="store",
        default="test_fasta.fasta",
        help="Fasta File for parsing",
    )
    parser.add_argument(
        "-d",
        "--outfastadir",
        action="store",
        default="splitoutput",
        help="Fasta File output directory",
    )

    args = parser.parse_args()
    split(fastafile=args.fastafile, outfastadir=args.outfastadir)
