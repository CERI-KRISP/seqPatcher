library(Biostrings)
library(stringr)

## load fasta file into R 
#inFasta <- readAAStringSet("aminoAcid.fasta") ## for amino acid fasta
inFasta <- readDNAStringSet("run57/50_documents_from_Run_57__73__74.fasta")  ## for dna fasta

## get seq names from fasta 
fa_given_names <- names(inFasta)

## prepare data frame, 
df <- data.frame(seq_name = names(inFasta) , new_name = sapply(str_split(fa_given_names, "_", n = 2, simplify = FALSE), `[`, 1))

## assign new seq names  by mapping fasta seq name to data frame names
names(inFasta) <- df[match(fa_given_names , df$seq_name) , "new_name"]

## write data to fasta file with updated names
writeXStringSet(inFasta , "run57/50_documents_from_Run_57_73_74_clean.fasta")
