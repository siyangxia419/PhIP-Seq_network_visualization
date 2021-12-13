blastp_dm <- function(pep_dt, fasta_dir = "C:/Users/siyang_xia/R/"){
  
  # install and load the required packages
  if(!("rBLAST" %in% installed.packages())) install.packages("rBLAST")
  if(!("seqinr" %in% installed.packages())) install.packages("seqinr")
  if(!("tidyverse" %in% installed.packages())) install.packages("tidyverse")
  
  library(rBLAST)
  library(seqinr)
  library(tidyverse)
  
  # write the sequences to FASTA files
  write.fasta(sequences = as.list(pep_dt$pep_aa),
              names = pep_dt$u_pep_id, 
              file.out = paste0(fasta_dir, "db.fasta"))
  
  write.fasta(sequences = as.list(pep_dt$pep_aa),
              names = pep_dt$u_pep_id, 
              file.out = paste0(fasta_dir, "query.fasta"))
  
  # make the database
  makeblastdb(file = paste0(fasta_dir, "db.fasta"),
              dbtype = "prot")
  dbp <- blast(db = paste0(fasta_dir, "db.fasta"), 
               type="blastp")
  
  # query sequences
  seqp <- readAAStringSet(paste0(fasta_dir, "query.fasta"))
  
  # blast calculation
  predp <- predict(dbp, seqp,
                   BLAST_args = "-evalue .1 -max_hsps 1 -soft_masking false -word_size 7 -max_target_seqs 100000",
                   custom_format = "qseqid sseqid pident length evalue bitscore positive gaps ppos")
  
  
  # arrange the peptide ids
  predp <- predp %>% 
    mutate(qseqid = ordered(qseqid, levels = pep_dt$u_pep_id),
           sseqid = ordered(sseqid, levels = pep_dt$u_pep_id)) %>% 
    rename(subject_id = qseqid, pattern_id = sseqid)
  
  predp[predp$subject_id > predp$pattern_id, c("subject_id", "pattern_id")] <- 
    predp[predp$subject_id > predp$pattern_id, c("pattern_id", "subject_id")]
  
  
  # add information of the peptides
  predp <- predp %>% 
    arrange(subject_id, pattern_id) %>% 
    as_tibble() %>% 
    left_join({pep_dt %>% rename_with(~paste0("subject_", .))}, 
              by = c("subject_id" = "subject_u_pep_id")) %>% 
    left_join({pep_dt %>% rename_with(~paste0("pattern_", .))}, 
              by = c("pattern_id" = "pattern_u_pep_id"))
  
  return(predp)
}
