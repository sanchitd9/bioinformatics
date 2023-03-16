require(seqinr)
require(Biostrings)
require(ORFhunteR)
require(ggplot2)
require(gridExtra)

# Setup
choosebank("genbank")

# Accession numbers
ac_1 <- "AY884001"
ac_2 <- "MH940245"

# Task A
q1 <- query("q1", paste("AC=", ac_1))
q2 <- query("q2", paste("AC=", ac_2))

seq1 <- getSequence(q1$req[[1]])
seq2 <- getSequence(q2$req[[1]])

print(DNAString(c2s(seq1)))
print(DNAString(c2s(seq2)))

#-----------------------------------------#
# Use this if connection is broken

# q1 <- read.fasta(file = paste("hw4/", ac_1, ".fasta", sep = ""))
# q2 <- read.fasta(file = paste("hw4/", ac_2, ".fasta", sep = ""))
# 
# seq1 <- getSequence(q1)[[1]]
# seq2 <- getSequence(q2)[[1]]
#-----------------------------------------#

# Task B
p1 <- seqinr::translate(seq1)
p2 <- seqinr::translate(seq2)

seq1_count <-table(p1)
seq2_count <- table(p2)
print(seq1_count)
print(seq2_count)

seq1_prop <- proportions(seq1_count)
seq2_prop <- proportions(seq2_count)
print(seq1_prop)
print(seq2_prop)

df1 <- data.frame(slices = seq1_prop, "amino_acids" = paste(names(seq1_prop), "=", round(seq1_prop * 100, 2), "%", sep = ""))
df2 <- data.frame(slices = seq2_prop, amino_acids = paste(names(seq2_prop), "=", round(seq2_prop * 100, 2), "%", sep = ""))
plot1 <- ggplot(df1, aes(x = "", y = slices.Freq, fill = amino_acids)) + geom_bar(stat = "identity", width = 1) + coord_polar("y", start = 0) + ggtitle(ac_1) + xlab("") + ylab("")
plot2 <- ggplot(df2, aes(x = "", y = slices.Freq, fill = amino_acids)) + geom_bar(stat = "identity", width = 1) + coord_polar("y", start = 0) + ggtitle(ac_2) + xlab("") + ylab("")
grid.arrange(plot1, plot2, ncol = 2)

# Task C
dotPlot(s2c(seq1_orfs[1, "orf.sequence"]), s2c(seq2_orfs[1, "orf.sequence"]), xlab = ac_1, ylab = ac_2)

# Task D
opt_ga <- pairwiseAlignment(DNAString(c2s(seq1)), DNAString(c2s(seq2)), type = "global", substitutionMatrix = nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE), gapOpening = 0, gapExtension = -2)
print(opt_ga)

cat(paste(substr(pattern(opt_ga), 1, 20), substr(subject(opt_ga), 1, 20), sep = '\n'))

# Task E
seq1_random <- generateSeqsWithMultinomialModel(c2s(seq1), 20)
seq2_random <- generateSeqsWithMultinomialModel(c2s(seq2), 20)

scores <- c()
for (i in 1:length(random_seq)) {
  score <- pairwiseAlignment(DNAString(c2s(seq1_random[[i]])), DNAString(c2s(seq2_random[[i]])), type = "global", substitutionMatrix = nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE), gapOpening = 0, gapExtension = -2, scoreOnly = TRUE)
  print(score)
  scores <- append(scores, score)
}

print(scores)
hist(scores)

p_value <- sum(scores > score(opt_ga)) / length(scores)
print(p_value)

# Task F
opt_la <- pairwiseAlignment(DNAString(c2s(seq1)), DNAString(c2s(seq2)), type = "local", substitutionMatrix = nucleotideSubstitutionMatrix(match = 3, mismatch = -2, baseOnly = TRUE), gapOpening = -4, gapExtension = -2)
print(opt_la)

closebank()
