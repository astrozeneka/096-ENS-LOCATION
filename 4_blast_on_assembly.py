from Bio.Blast.Applications import NcbiblastnCommandline
import argparse
import os
from Bio.SeqIO import parse
from Bio.SeqIO import write

HIT_NUMBER_THRESHOLD = 2

# Blast on assembly in order to exclude the potential repeat
parser = argparse.ArgumentParser()
parser.add_argument("--slug", type=str, help="Name of the fasta file without the .fasta extension")
args = parser.parse_args()

if __name__ == '__main__':
    # Open the cleaned sequences
    cline = NcbiblastnCommandline(query=f"./data/6_cleaned/{args.slug}-cleaned.fasta", db=f"./data/1_assemblies/{args.slug}.fasta", outfmt=6, out=f"./data/7_blast_on_assembly/{args.slug}-blast.tsv")
    os.system(str(cline))

    # Open the blast results
    blast_result = open(f"./data/7_blast_on_assembly/{args.slug}-blast.tsv").read().strip().split("\n")
    blast_result = [a.split("\t") for a in blast_result]
    query_list = list(set([a[0] for a in blast_result]))
    hit_number = {query: 0 for query in query_list}

    for query in query_list:
        hit_number[query] = len([a for a in blast_result if a[0] == query])

    filtered_query = [query for query in query_list if hit_number[query] <= HIT_NUMBER_THRESHOLD]

    # Open the cleaned sequences
    records = parse(f"./data/6_cleaned/{args.slug}-cleaned.fasta", "fasta")
    output_records = []
    for record in records:
        if record.id in filtered_query:
            output_records.append(record)
    # Write output file using SeqIO write
    write(output_records, f"./data/8_cleaned_and_filtered/{args.slug}-Step5.fasta", "fasta")

    print()

