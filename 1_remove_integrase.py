from Bio.Blast.Applications import NcbiblastnCommandline
from os import path
from Bio.SeqIO import parse
from Bio.SeqIO import write
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--slug", type=str, help="Name of the fasta file without the .fasta extension")
args = parser.parse_args()

EVALUE_THRESHOLD=10e-10

if __name__ == '__main__':
    # build blast db
    os.system(f"makeblastdb -in ./data/1_assemblies/{args.slug}.fasta -dbtype nucl -hash_index")

    # Step 1 blast
    cline = NcbiblastnCommandline(query=f"./data/sequences/integrase.fasta", db=f"./data/1_assemblies/{args.slug}.fasta", evalue=EVALUE_THRESHOLD, outfmt=6, out=f"./data/2_blast_out/{args.slug}.tsv")
    a = os.system(str(cline))

    blast_results = open(f"./data/2_blast_out/{args.slug}.tsv").read().strip().split("\n")
    blast_results = [a.split("\t") for a in blast_results]

    # Step 2 read the blast target, and replace the blasted sequence by "NNN..."
    # open file using bio python
    records = parse(f"./data/1_assemblies/{args.slug}.fasta", "fasta")
    output_records = []
    for record in records:
        related_blasts = [a for a in blast_results if a[1] == record.id]
        for row in related_blasts:
            start = int(row[8])
            end = int(row[9])
            is_on_positive_strand = start < end
            if is_on_positive_strand:
                record.seq = record.seq[:start] + "N" * (end - start) + record.seq[end:]
            else:
                record.seq = record.seq[:end] + "N" * (start - end) + record.seq[start:]
        output_records.append(record)
    # Write output file using SeqIO write
    write(output_records, f"./data/3_integrase-free/{args.slug}-Step1.fasta", "fasta")
    # build blast db
    os.system(f"makeblastdb -in ./data/3_integrase-free/{args.slug}-Step1.fasta -dbtype nucl -hash_index")
    print("Done")
    print()


