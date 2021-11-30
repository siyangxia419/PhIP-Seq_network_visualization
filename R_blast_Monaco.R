library(data.table)
library(rBLAST)
library(seqinr)
library(Biostrings)
library(doParallel)
library(foreach)
library(tidyverse)
library(foreach)
library(doParallel)
setwd("/home-net/home-1/dmonaco1@jhu.edu/Blast/")

# makeblastdb(file = "uniprot.fasta",dbtype = "prot") # point to fasta reference file (what you blast against) - set as nucl or prot
dbp <- blast("uniprot.fasta",type="blastp") # after above line is run use this to read in generated local database
seqp <- readAAStringSet("prot.fasta") #points to query fasta file (either readAAStringset or readDNAstringset)

registerDoParallel(detectCores())

zeta = foreach(R = 1:10,.combine=rbind) %dopar%{ # i used a foreach parallelized to make this faster - this runs x proteins agaisnt the database at the saem time
  
  # the below predict line is the acutal blast search - the "args" are stardard blast seetting you can find online - here we required an evalue of .1 minimum
predp1 <- predict(dbp,seqp[R],BLAST_args = "-evalue .1 -max_hsps 1 -soft_masking false -word_size 7 -max_target_seqs 100000",
                  custom_format = "qseqid sseqid pident length evalue bitscore positive gaps ppos")

# the custom format line just indicates various metrics you would like reported

fwrite(predp1 %>% as.data.frame(),paste("./output/",R,"blast_cv_rep2.csv"))

}

