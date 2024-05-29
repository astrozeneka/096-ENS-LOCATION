from Bio.Blast.Applications import NcbiblastnCommandline
from os import path
from Bio.SeqIO import parse
from Bio.SeqIO import write
# Import record
from Bio.SeqRecord import SeqRecord
import os
import argparse

N_THRESHOLD=0.6

parser = argparse.ArgumentParser()
parser.add_argument("--slug", type=str, help="Name of the fasta file without the .fasta extension")
args = parser.parse_args()


if __name__ == '__main__':

    # open the int+ltr free file
    records = parse(f"./data/5_int-ltr-free/{args.slug}-Step4.fasta", "fasta")
    output_records = []
    for record in records:
        # The percentage of masked sequence should be more than N_THRESHOLD
        if record.seq.count("N") / len(record.seq) < N_THRESHOLD:
            continue
        sub_reads = record.seq.split("N")
        sub_reads = [sub_reads[0]] + [sub_reads[-1]]
        sub_reads = [a for a in sub_reads if len(a) > 80]
        for i, sub_read in enumerate(sub_reads):
            record = SeqRecord(sub_read, id=f"{record.id}_{i+1}", description="")
            output_records.append(record)
    # Write output file using SeqIO write
    write(output_records, f"./data/6_cleaned/{args.slug}-cleaned.fasta", "fasta")
    print("Done")