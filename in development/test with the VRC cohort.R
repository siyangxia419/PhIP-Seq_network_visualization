library(virlink)
library(tidyverse)

VRC_hits <- read.table(file = "VRC_VirscanLar_000_Hits_foldchange_annotated.tsv", 
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
VRC_hits <- tibble::as_tibble(VRC_hits)

VRC_example <- VRC_hits %>% slice_sample(n = 1000)

VRC_example_pep <- VRC_example[, 1:10]
VRC_example_ab <- VRC_example[, c(1, 11:ncol(VRC_example))]

VRC_example_pep <- VRC_example_pep %>% 
  rename(id = u_pep_id)

system.time({
VRC_pairwise_alignment <- peptide_pairwise_alignment(peptides = VRC_example_pep, 
                                                     sub_matrix = "BLOSUM62",
                                                     output_str = "tibble")
})
