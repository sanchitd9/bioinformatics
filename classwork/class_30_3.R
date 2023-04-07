require(seqinr)
require(ape)

choosebank("swissprot")
acs <- c("P06747", "P0C569", "O56773", "Q5VKP1")
seqs <- c()

for (x in acs) {
  q <- query("q", paste("AC=", x, sep = ""))
  seqs <- append(seqs, q$req)    
}

for (x in seqs) {
  temp <- getSequence(x)
  write.fasta(temp, getName(x), file.out = paste(getName(x), ".fasta", sep = ""))
}

x <- read.alignment(file = "C:/Users/sanch/Documents/PFW/Bioinformatics/homework/srutest.aln", format = "clustal")
printMultipleAlignment(x)

y <- cleanAlignment(x, 70, 33)
y

dist.alignment(x)
dist.alignment(y)

scores <- read.table(file = "srutestscore.qscores")
scores

mean(scores$V5)

closebank()
###############################################################

choosebank("genbank")

acs <- c("AF049118", "AF049114", "AF049119", "AF049115")
seqs <- c()
names < c()
for (x in acs) {
  q <- query("q", paste("AC=", x, sep = ""))
  seqs <- append(seqs, q$req)
  names <- append(names, getName(q$req))
}
names
x <- lapply(seqs, getSequence)

write.fasta(x, names = names, file.out = "exercise.fasta")

alignment <- read.alignment(file = "exercise.aln", format = "clustal")
alignment

dist.alignment(alignment)

dist.dna(as.DNAbin(alignment))


