require(seqinr)
require(Biostrings)
require(ggplot2)
require(gridExtra)

choosebank("swissprot")
data("PAM250")
data("BLOSUM62")

# Task A
ac_1 <- "A0A1I9G1P5"
ac_2 <- "A0A1S0U2K5"

q1 <- query("q1", query = paste("AC=", ac_1, sep = ""))
q2 <- query("q2", query = paste("AC=", ac_2, sep = ""))

seq1 <- getSequence(q1$req[[1]])
seq2 <- getSequence(q2$req[[1]])

print(AAString(c2s(seq1)))
print(AAString(c2s(seq2)))

# Task B
count_seq1 <- table(seq1)
count_seq2 <- table(seq2)

print(count_seq1)
print(count_seq2)

prop_seq1 <- proportions(count_seq1)
prop_seq2 <- proportions(count_seq2)

print(prop_seq1)
print(prop_seq2)

theme_update(legend.key.size = unit(4, "mm"), legend.text = element_text(size = 6), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid  = element_blank())
df1 <- data.frame(slices = prop_seq1, amino_acids = paste(names(prop_seq1), "=", round(prop_seq1 * 100, 2), "%", sep = ""))
df2 <- data.frame(slices = prop_seq2, amino_acids = paste(names(prop_seq2), "=", round(prop_seq2 * 100, 2), "%", sep = ""))
plot1 <- ggplot(df1, aes(x = "", y = slices.Freq, fill = amino_acids)) + geom_bar(stat = "identity", width = 1) + coord_polar("y", start = 0) + ggtitle(ac_1) + xlab("") + ylab("") + guides(fill = guide_legend(title = NULL))
plot2 <- ggplot(df2, aes(x = "", y = slices.Freq, fill = amino_acids)) + geom_bar(stat = "identity", width = 1) + coord_polar("y", start = 0) + ggtitle(ac_2) + xlab("") + ylab("") + guides(fill = guide_legend(title = NULL))
grid.arrange(plot1, plot2, ncol = 2)


# Task C
dotPlot(tail(seq1, 50), tail(seq2, 50), xlab = ac_1, ylab = ac_2)

# Task D
ga <- pairwiseAlignment(AAString(c2s(seq1)), AAString(c2s(seq2)), substitutionMatrix = PAM250, type = "global", gapOpening = -10, gapExtension = -0.5)
print(ga)

printPairwiseAlignment(ga, chunksize = 30)

# Task E
la <- pairwiseAlignment(AAString(c2s(seq1)), AAString(c2s(seq2)), substitutionMatrix = BLOSUM62, type = "local", gapOpening = -10, gapExtension = -0.5)
print(la)

printPairwiseAlignment(la, chunksize = 30)

# Task F
rand_seq1 <- generateSeqsWithMultinomialModel(c2s(seq1), 10)
rand_seq2 <- generateSeqsWithMultinomialModel(c2s(seq2), 10)

scores_global <- c()
scores_local <- c()

for (i in 1:length(rand_seq1)) {
  score_global <- pairwiseAlignment(pattern = AAString(c2s(rand_seq1[[i]])), subject = AAString(c2s(rand_seq2[[i]])), substitutionMatrix = PAM250, gapOpening = -10, gapExtension = -0.5, type = "global", scoreOnly = TRUE)
  score_local <- pairwiseAlignment(pattern = AAString(c2s(rand_seq1[[i]])), subject = AAString(c2s(rand_seq2[[i]])), substitutionMatrix = BLOSUM62, gapOpening = -10, gapExtension = -0.5, type = "local", scoreOnly = TRUE)
  scores_global <- append(scores_global, score_global)
  scores_local <- append(scores_local, score_local)
}

hist(scores_global)
hist(scores_local)

p_global <- sum(scores_global > score(ga)) / length(scores_global)
p_local <- sum(scores_local > score(la)) / length(scores_local)

print(p_global)
print(p_local)
