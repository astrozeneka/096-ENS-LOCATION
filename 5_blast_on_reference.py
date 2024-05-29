from Bio.Blast.Applications import NcbiblastnCommandline
import os
import argparse

PERC_IDENTITY_THRESHOLD=98

parser = argparse.ArgumentParser()
parser.add_argument("--slug", type=str, help="Name of the fasta file without the .fasta extension")
args = parser.parse_args()

if __name__ == '__main__':
    # Blast on reference
    cline = NcbiblastnCommandline(query=f"./data/8_cleaned_and_filtered/{args.slug}-Step5.fasta", db=f"./data/genome/GCA_024206055.2_GGswu_genomic.fna", perc_identity=PERC_IDENTITY_THRESHOLD, outfmt=6, out=f"./data/9_blast_on_ref/{args.slug}-blast.tsv")
    a = os.system(str(cline))
    print("Done")