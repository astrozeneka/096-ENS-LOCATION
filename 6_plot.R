# install
install.packages("RIdeogram")
require(RIdeogram)

data(human_karyotype, package="RIdeogram")
data(gene_density, package="RIdeogram")
data(Random_RNAs_500, package="RIdeogram")

an_to_chrid_dict = "CP100555.1	1
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
CP115610.1	MT"
# Convert the string above to dictionary
an_to_chrid_dict = strsplit(an_to_chrid_dict, "\n")[[1]]
an_to_chrid_dict = strsplit(an_to_chrid_dict, "\t")
an_to_chrid_dict = setNames(as.numeric(sapply(an_to_chrid_dict, "[[", 2)), sapply(an_to_chrid_dict, "[[", 1))


karyotype = "1	0	200044509
2	0	152127819
3	0	112377238
4	0	91364047
5	0	59473047
6	0	36158469
7	0	36439816
8	0	29594148
9	0	23966833
10	0	21606253
11	0	19822308
12	0	20591337
13	0	19253982
14	0	16259740
15	0	13529386
16	0	4932503
17	0	11322113
18	0	12199395
19	0	11295715
20	0	15361164
21	0	7106554
22	0	5772503
23	0	6970121
24	0	8591846
25	0	4095597
26	0	6534586
27	0	6322964
28	0	5942585
29	0	6839231
30	0	3643476
31	0	5478923
32	0	3130476
33	0	4180360
34	0	4011625
35	0	3056416
36	0	3424470
37	0	2501204
38	0	3700872
W	0	14151827
Z	0	87735852"
# conver to data table
karyotype = read.table(text=karyotype, sep="\t", header=FALSE, stringsAsFactors=FALSE)
# set column names as Chr, Start, End
colnames(karyotype) = c("Chr", "Start", "End")
# plot empty karyotype using Rideogram
ideogram(karyotype = karyotype)
convertSVG("chromosome.svg", device = "png")

# load the slug from the parameter argument and set "__" as default
slug = "SRR24608405"
if (length(commandArgs(trailingOnly=TRUE)) > 0) {
    slug = commandArgs(trailingOnly=TRUE)[1]
}


# Open the "data/7_blast_on_ref/SRR24608420-ENS-1-blast.tsv" file and load it into data table
# the headers are the blast format 6 headers : seqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
#blast_data = read.table("data/7_blast_on_ref/SRR24608420-ENS-1-blast.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)
blast_data = read.table(paste("data/9_blast_on_ref/", slug, "-ENS-1-blast.tsv", sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(blast_data) = c("seqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
# TODO: remove this later
# only rows having e-values < 10e-10
#blast_data = blast_data[blast_data$evalue < 10e-60,]
nrow(blast_data)

# Create a new table having the following column
# type: ENS
# shape: triangle
# chr: sseqid column
# start: sstart column
# end: send column
# color: <whatever>
# (for each row in the blast_data)
loci_data = data.frame(Type="ENS-1", Shape="triangle", Chr=blast_data$sseqid, Start=blast_data$sstart, End=blast_data$send, color="33a02c")

# use the an_to_chrid_dict to conver the "Chr" column to chromosome id
loci_data$Chr = an_to_chrid_dict[loci_data$Chr]
# drop the rows with NA
loci_data = loci_data[!is.na(loci_data$Chr),]
# Convert the Chr column to integer
loci_data$Chr = as.integer(loci_data$Chr)
# Then convert to string again
loci_data$Chr = as.character(loci_data$Chr)
# As well as Start and End
loci_data$Start = as.integer(loci_data$Start)
loci_data$End = as.integer(loci_data$End)

# Row row where End > Start, switch the two values
loci_data[loci_data$End < loci_data$Start, c("Start", "End")] = loci_data[loci_data$End < loci_data$Start, c("End", "Start")]
# Only 100 rows


# Plot on the karyotype and specify outfile
ideogram(karyotype = karyotype, label = loci_data, label_type = "marker", output = paste("data/10_karyotype/", slug, "-ENS-1-blast.svg", sep=""))
# Conver to png and save it in the same directory as the svg
#convertSVG("data/8_karyotype/SRR24608420-ENS-1-blast.svg", device = "png")
convertSVG(paste("data/10_karyotype/",slug,"-ENS-1-blast.svg", sep=""), device = "png")
# Copy the "chromosome.png" to "data/8_karyotype/SRR24608420-ENS-1-blast.png" even if the destination exists
file.copy("chromosome.png", paste("data/10_karyotype/", slug, "-ENS-1-blast.png", sep=""), overwrite=TRUE)

