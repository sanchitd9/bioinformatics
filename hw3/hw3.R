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
print(orfs_rev_comp[, c("start", "end", "length")])

# Plot potential ORFs in the last 1000 bases
plotORFsinSeq(c2s(tail(seq_vector, 1000)))

# Extract, translate and print the longest potential gene
potential_orfs <- findORFs(as.character(seq_str))
longest_gene <- potential_orfs[as.numeric(potential_orfs[, "length"]) == max(as.numeric(potential_orfs[, "length"]))]
print(paste("The longest gene starts at index", longest_gene[1], "and ends at index", longest_gene[2]))
print(paste("Length of the longest gene:", longest_gene[3]))
print("The longest gene:")
print(DNAString(longest_gene[4]))

protein <- translate(DNAString(longest_gene[4]))
print(paste("The resulting protein sequence of length:", length(protein)))
print(protein)

# Identify the significant ORFs using the 95th percentile as the threshold value.
random_seqs <- generateSeqsWithMultinomialModel(c2s(seq_vector), 30)
random_orfs <- lapply(random_seqs, findORFs)

length_vector <- sapply(random_orfs, get_max_seq_vector)

threshold <- quantile(length_vector, 0.95)

significant_orfs <- potential_orfs[as.numeric(potential_orfs[, "length"]) > threshold, ]
print("Significant ORFs:")
print(significant_orfs[, c("start", "end", "length")])
