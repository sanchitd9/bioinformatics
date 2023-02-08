require(seqinr)

# Set the bank
choosebank("genbank")

# ----------------------------
# Task A
# ----------------------------

# Search for DNA sequences of "Chlamydia trachomatis"
result <- query("query1", query = "SP=Chlamydia trachomatis AND M=DNA")

# Fetch the attributes
attributes(result)

# Print the number of elements in "result"
print(result$nelem)



# ----------------------------
# Task B
# ----------------------------

# Get all the sequences
all_sequences <- result$req
# Get the length vector
length_vector <- getLength(all_sequences)

# Get the shortest sequence(s)
shortest_seq <- all_sequences[length_vector == min(length_vector)]

# Print the length of the shortest sequence(s)
print(min(length_vector))

# Print the accession number(s) of the shortest sequence(s)
print(paste("Accession Numbers:", paste(getName(shortest_seq), collapse = ", ")))

# Export the data to a FASTA file
write.fasta(getSequence(shortest_seq), names = getName(shortest_seq), file.out = "my_output.fasta")

# ----------------------------
# Task C
# ----------------------------

