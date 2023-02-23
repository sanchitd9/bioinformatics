# Load the packages
require(seqinr)
require(Biostrings)
require(ORFhunteR)

# Connect to the bank
choosebank("genbank")

# Accession number for SARS Coronavirus MA15 ExoN1
accession_number <- "FJ882953"

# Query
q <- query("q", paste("AC=", accession_number, sep = ""))

# Query attributes
attributes(q)

# Print the details
print(paste("Retrieved", q$nelem, "sequence(s) for Accession Number:", accession_number))
print(paste("Length of the sequence:", getLength(q)))
seq_vector <- getSequence(q$req)[[1]]
seq_str <- DNAString(c2s(seq_vector))
print(seq_str)

# Reverse complement of the sequence
reverse_comp_seq <- reverseComplement(seq_str)
print(reverse_comp_seq)

# Print all the potential ORFs in the reverse complement of the sequence
orfs_rev_comp <- findORFs(as.character(reverse_comp_seq))
print(orfs_rev_comp[, c("start", "end", "length")], max = 10000)

# Plot potential ORFs in the last 1000 bases
plotORFsinSeq(tail(reverse_comp_seq, 1000))

# Extract, translate and print the longest potential gene
max_length <- max(as.numeric(orfs_rev_comp[, "length"]))
longest_gene <- orfs_rev_comp[as.numeric(orfs_rev_comp[, "length"]) == max_length, ]

print(paste("Length of the maximum potential gene:",max_length))
print(paste("There are", dim(longest_gene)[1], "sequences with the maximum length of", max_length))

for (i in 1:dim(longest_gene)[1]) {
  print(paste("Sequence", i, ": Start =", longest_gene[i, "start"], ", End =", longest_gene[i, "end"]))
  l_seq <- DNAString(longest_gene[i, "orf.sequence"])
  print(l_seq)
  protein <- Biostrings::translate(l_seq)
  print("The resulting protein sequence:")
  print(protein)
}

# Identify the significant ORFs using the 95th percentile as the threshold value.
random_seqs <- generateSeqsWithMultinomialModel(c2s(reverse_comp_seq), 30)
random_orfs <- lapply(random_seqs, findORFs)

length_vector <- sapply(random_orfs, get_max_seq_vector)
length_vector <- unlist(length_vector)

threshold <- quantile(length_vector, 0.95)
cat(threshold)

significant_orfs <- orfs_rev_comp[as.numeric(orfs_rev_comp[, "length"]) > threshold, ]
print("Significant ORFs:")

for (i in 1:dim(significant_orfs)[1]) {
  print(paste("Sequence", i, ": Start =", significant_orfs[i, "start"], ", End =", significant_orfs[i, "end"], ", Length =", significant_orfs[i, "length"]))
  print(DNAString(significant_orfs[i, "orf.sequence"]))
}

