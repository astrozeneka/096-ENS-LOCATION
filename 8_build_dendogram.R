# Read the csv matrix
data <- read.csv("data/binary_matrix.csv", row.names=1)
# loci are in rows
# samples are in columns

# Delete the row where sum of the row is 1
data <- data[rowSums(data) != 1,]

# Load the breed names tsv file
breed_names = read.table("data/breed-names.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)
# Convert to dictionary (column 1 will map to column 2)
breed_names_dict = setNames(breed_names$V2, breed_names$V1)

# Update the column name to be (Accession number + (breed name))
colnames(data) <- sapply(colnames(data), function(x) paste(x,
  paste('(', breed_names_dict[x], ')', sep=""), sep=" "))

# Calculate the distance matrix between the samples (columns)
distance_matrix_sample <- dist(t(data), method = "euclidean")

# Cluster the sample using hclust (Rotate the figure 90 degree)
dendogram_sample = hclust(distance_matrix_sample, method = "ward.D2")
#dendogram_sample = hclust(distance_matrix_sample, method = "ward.D2")
pdf("data/sample_clustering_v2.pdf", height=12, width=6)
plot(dendogram_sample)
dev.off()