require(seqinr)
require(Biostrings)

choosebank("swissprot")

q1 <- query("q1", "AC=Q9CD83")
q2 <- query("q2", "AC=A0PQ23")

seq1 <- getSequence(q1$req[[1]])
seq2 <- getSequence(q2$req[[1]])

seq1
seq2

dotPlot(seq1, seq2, wsize = 1, wstep = 1)


seq1 <- DNAString("GAATTC")
seq2 <- DNAString("GATTA")

scoring_matrix <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)

global_alignment <- pairwiseAlignment(seq1, seq2, type = "global", substitutionMatrix = scoring_matrix, gapOpening = -8, gapExtension = -2)


pairwiseAlignment(DNAString("ACCAGGTACA"), DNAString("CAG"), type = "local", substitutionMatrix = nucleotideSubstitutionMatrix(baseOnly = T), gapOpening = -1, gapExtension = 0)

pairwiseAlignment(DNAString("CTTCGG"), DNAString("CGTCG"), type = "global", substitutionMatrix = nucleotideSubstitutionMatrix(baseOnly = T), gapOpening = -1, gapExtension = 0)


den1 <- read.fasta(file = "den1.fasta")
den2 <- read.fasta(file = "den2.fasta")

den1 <- DNAString(c2s(den1$NC_001477.1))
den2 <- DNAString(c2s(den2$NC_001474.2))


pairwiseAlignment(den1, den2, type = "global", substitutionMatrix = nucleotideSubstitutionMatrix(baseOnly = TRUE), gapOpening = 0, gapExtension = -1)


choosebank("refseqViruses")
q1 <- query("q1", "AC=NC_001477")
q2 <- query("q2", "AC=NC_001474")

seq1 <- getSequence(q1$req)[[1]]
seq2 <- getSequence(q2$req)[[1]]
seq1
seq2

alignment_dengue <- pairwiseAlignment(toupper(c2s(seq1)), toupper(c2s(seq2)), type = "global", substitutionMatrix = nucleotideSubstitutionMatrix(baseOnly = TRUE), gapOpening = 0, gapExtension = -1, scoreOnly = TRUE)

random_d1 <- generateSeqsWithMultinomialModel(c2s(seq1), 10)
random_d2 <- generateSeqsWithMultinomialModel(c2s(seq2), 10)

score1_v <- c()
score2_v <- c()
for (i in 1:10) {
  score1 <- pairwiseAlignment(toupper(c2s(seq1)), toupper(random_d1[[i]]), type = "global", substitutionMatrix = nucleotideSubstitutionMatrix(baseOnly = TRUE), gapOpening = 0, gapExtension = -1, scoreOnly = TRUE)
  score2 <- pairwiseAlignment(toupper(c2s(seq2)), toupper(random_d1[[i]]), type = "global", substitutionMatrix = nucleotideSubstitutionMatrix(baseOnly = TRUE), gapOpening = 0, gapExtension = -1, scoreOnly = TRUE)
  score1_v <- append(score1_v, score1)
  score2_v <- append(score1_v, score2)
  
}


hist(score1_v)
hist(score2_v) 

p1 <- sum(score1_v > alignment_dengue) / length(score1_v)
p2 <- sum(score2_v > alignment_dengue) / length(score2_v)
p1
p2
