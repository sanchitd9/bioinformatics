require(seqinr)
require(ape)
require(Biostrings)
require(dplyr)


# Task A
seqs <- read.fasta(file = "hw6/usflu.fasta")
lapply(unlist(getSequence(seqs, as.string = TRUE)), DNAString)

# Task B
annotations <- read.csv("hw6/usflu.annot.csv")
strains <- annotations %>% group_by(year) %>% slice_head(n = 1)
print(strains)

filtered_seqs <- seqs[strains$accession]
write.fasta(sequences = filtered_seqs, names = getName(filtered_seqs), file.out = "hw6/filtered_seqs.fasta")


# Task C is on clustal


# Task D
aln <- read.alignment(file = "hw6/filtered_seqs.phy", format = "phylip")
aln



