
from glob import glob
import os
import pandas as pd

an_to_chrid_dict = """CP100555.1	1
CP100556.1	2
CP100557.1	3
CP100558.1	4
CP100559.1	5
CP100560.1	6
CP100561.1	7
CP100562.1	8
CP100563.1	9
CP100564.1	10
CP100565.1	11
CP100566.1	12
CP100567.1	13
CP100568.1	14
CP100569.1	15
CP100570.2	16
CP100571.1	17
CP100572.1	18
CP100573.1	19
CP100574.1	20
CP100575.1	21
CP100576.1	22
CP100577.1	23
CP100578.1	24
CP100579.1	25
CP100580.1	26
CP100581.1	27
CP100582.1	28
CP100583.2	29
CP100584.2	30
CP100585.2	31
CP100586.2	32
CP100587.1	33
CP100588.2	34
CP100589.2	35
CP100590.2	36
CP100591.2	37
CP100592.2	38
CP100593.1	W
CP100594.1	Z
CP115610.1	MT"""
# convert to dictionary
an_to_chrid_dict = {a.split()[0]:a.split()[1] for a in an_to_chrid_dict.split("\n")}

if __name__ == '__main__':
    # Use the 9_blast_on_ref_data to build the absence/presence binary matrix
    file_list = glob("./data/9_blast_on_ref/*.tsv")
    genome_list = [os.path.basename(a).split("-")[0] for a in file_list]

    output = {}
    for file, genome in zip(file_list, genome_list):
        blast_output = open(file).read().strip().split("\n")
        blast_output = [a.split("\t") for a in blast_output]
        for blast_row in blast_output:
            loci = f"{an_to_chrid_dict[blast_row[1]]}_{blast_row[8]}"
            if loci not in output:
                output_row = {a:0 for a in genome_list}
                output_row[genome] = 1
                output[loci] = output_row
            else:
                output[loci][genome] = 1

    # Convert output to dataframe
    df = pd.DataFrame(output).T
    # export as a csv
    df.to_csv("./data/binary_matrix.csv")
    print()
    print()