require(seqinr)
require(ape)
require(Biostrings)
require(dplyr)
require(stringi)


# Task A
seqs <- read.fasta(file = "hw6/usflu.fasta")

# Task B
annotations <- read.csv("hw6/usflu.annot.csv")
strains <- annotations %>% group_by(year) %>% slice_head(n = 1)
print(strains)

filtered_seqs <- seqs[strains$accession]
write.fasta(sequences = filtered_seqs, names = getName(filtered_seqs), file.out = "hw6/filtered_seqs.fasta", open = "w")

# Task D
aln <- read.alignment(file = "hw6/filtered_seqs.phy", format = "phylip")
print(aln)

clean_aln <- cleanAlignment(aln, 75, 30)
print(clean_aln)

dist_matrix <- dist.dna(as.DNAbin(clean_aln))
print(dist_matrix)

max(dist_matrix)  # FJ549055, CY012480
min(dist_matrix)  # CY000737, CY001453

# Task E
unrootedNJtree(clean_aln, "DNA", "phylogram")

# Task G
choosebank("genbank")
q <- query("q", "AC=AB257344")
sars_seq <- getSequence(q$req)[[1]]
print(DNAString(c2s(sars_seq)))
write.fasta(sars_seq, getName(q), file.out = "hw6/filtered_seqs.fasta", open = "a")
new_aln <- read.alignment(file = "hw6/filtered_seqs_1.phy", format = "phylip")
print(new_aln)

dist_matrix_1 <- dist.dna(as.DNAbin(new_aln))
print(dist_matrix_1)

rootedNJtree(new_aln, "AB257344", "DNA")

closebank()

# Task I
long_substr <- function (data) {
  is_substr <- function (find, data) {
    if(length(data) < 1 && nchar(find) < 1) {
      return (FALSE)
    }
    
    for (i in 1:length(data)) {
      if (!stri_detect_fixed(data[[i]], find)) {
        return (FALSE)
      }
    }
    
    return (TRUE)
  }
  
  data <- data$seq
  substr <- " "
  if (length(data) > 1 && nchar(data[[1]]) > 0) {
    for (i in 1:nchar(data[[1]])) {
      for (j in 1:(nchar(data[[1]]) - i + 1)) {
        if (j > nchar(substr) && is_substr(substring(data[[1]], i, i + j), data)) {
          substr <- substring(data[[1]], i, i + j)
        }
      }
    }
  }
  
  return (substr)
}

result <- long_substr(aln)
print(result)

aln

