from Bio.Blast.Applications import NcbiblastnCommandline
from os import path
from Bio.SeqIO import parse
from Bio.SeqIO import write
import os
import argparse

EVALUE_THRESHOLD=10e-10

parser = argparse.ArgumentParser()
parser.add_argument("--slug", type=str, help="Name of the fasta file without the .fasta extension")
args = parser.parse_args()

if __name__ == '__main__':

    # Step 3 blast
    cline = NcbiblastnCommandline(query=f"./data/sequences/LTR.fasta", db=f"./data/3_integrase-free/{args.slug}-Step1.fasta",
                                  evalue=EVALUE_THRESHOLD, outfmt=6, out=f"./data/4_blast_out/{args.slug}.tsv")
    a = os.system(str(cline))

    blast_results = open(f"./data/4_blast_out/{args.slug}.tsv").read().strip().split("\n")
    blast_results = [a.split("\t") for a in blast_results]

    # Step 4 read the blast target (integrase free), and replace the blasted sequence by "NNN..."
    records = parse(f"./data/3_integrase-free/{args.slug}-Step1.fasta", "fasta")
    output_records = []
    for record in records:
        related_blasts = [a for a in blast_results if a[1] == record.id]
        masked = record.seq
        for row in related_blasts:
            start = int(row[8])
            end = int(row[9])
            is_on_positive_strand = start < end
            if is_on_positive_strand:
                masked = masked[:start] + "N" * (end - start) + masked[end:]
            else:
                masked = masked[:end] + "N" * (start - end) + masked[start:]
        record.seq = masked
        output_records.append(record)
    write(output_records, f"./data/5_int-ltr-free/{args.slug}-Step4.fasta", "fasta")
    # build database
    os.system(f"makeblastdb -in ./data/5_int-ltr-free/{args.slug}-Step4.fasta -dbtype nucl -hash_index")
    print("Done")
    print()

    print()