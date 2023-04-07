require(seqinr)
require(Biostrings)

ac <- "q9cd83"
ac1 <- "a0pq23"
data <- read.fasta(file = paste(ac, ".fasta", sep = ""))
seq <- getSequence(data[[1]])
seq <- toupper(seq)
seq
length(seq)

table(seq)

choosebank("swissprot")
q <- query("q", query = paste("AC=", ac, sep = ""))
seq1 <- getSequence(q$req[[1]])
seq1
length(seq1)


q1 <- query("q1", query = paste("AC=", ac1, sep = ""))
seq2 <- getSequence(q1$req[[1]])
seq2
length(seq2)

dotPlot(seq1, seq2)


data(package = "Biostrings")
data("PAM250")
data("BLOSUM62")
data("BLOSUM50")

s1 <- "PAWHEAE"
s2 <- "HEAGAWGHEE"

pairwiseAlignment(pattern = AAString(s1), subject = AAString(s2), substitutionMatrix = PAM250, gapOpening = -2, gapExtension = -8)

pairwiseAlignment(pattern = AAString(s1), subject = AAString(s2), substitutionMatrix = BLOSUM62, gapOpening = -2, gapExtension = -8)

res <- pairwiseAlignment(pattern = AAString(c2s(seq1)), subject = AAString(c2s(seq2)), type = "local", substitutionMatrix = BLOSUM50, gapOpening = -2, gapExtension = -8)
res

printPairwiseAlignment(res, chunksize = 40) 


th <- pairwiseAlignment(pattern = AAString(c2s(seq1)), subject = AAString(c2s(seq2)), substitutionMatrix = BLOSUM50, gapOpening = -2, gapExtension = -8, type = "global", scoreOnly = T)
th

rand1 <- generateSeqsWithMultinomialModel(c2s(seq1), 10)
rand2 <- generateSeqsWithMultinomialModel(c2s(seq2), 10)

scores <- c()

for (i in 1:length(rand1)) {
  temp <- pairwiseAlignment(AAString(c2s(rand1[[i]])), AAString(c2s(rand2[[i]])), substitutionMatrix = BLOSUM50, gapOpening = -2, gapExtension = -8, type = "global", scoreOnly = T)
  scores <- append(scores, temp) 
}

scores
hist(scores)
p <- sum(scores > th) / length(scores)
p
