###
# Prepare the epitope pairwise analysis data for R shiny network visualization 
#
# Siyang Xia (collaborate with Jenn L. Remmel)
# 2021-09-14
###


library(tidyverse)


# load the pairwise analysis dataset
data_dir <- "C:/Users/Siyang Xia/Dropbox/My laptop/Mina lab/Projects/VirScan/VirScan peptide cross-reactivity/Dataset 3 data exploration/"
load(paste0(data_dir, "VirScan_MikaelKnip_test_dataset_pairwise.RData"))



# data cleaning --------------------------------------------------------

# alignment dataset
epitope_align <- epitope_align %>% dplyr::arrange(subject_id, pattern_id)


# cooccurrence dataset
epitope_cooccur$id1 <- as.numeric(epitope_cooccur$id1)
epitope_cooccur$id2 <- as.numeric(epitope_cooccur$id2)
epitope_cooccur[epitope_cooccur$id1 > epitope_cooccur$id2, c("id1", "id2")] <- 
  epitope_cooccur[epitope_cooccur$id1 > epitope_cooccur$id2, c("id2", "id1")]
epitope_cooccur <- epitope_cooccur %>% dplyr::arrange(id1, id2)


# correlation dataset (no temporal structure)
epitope_cor$id1 <- as.numeric(epitope_cor$id1)
epitope_cor$id2 <- as.numeric(epitope_cor$id2)
epitope_cor[epitope_cor$id1 > epitope_cor$id2, c("id1", "id2")] <- 
  epitope_cor[epitope_cor$id1 > epitope_cor$id2, c("id2", "id1")]
epitope_cor <- epitope_cor %>% dplyr::arrange(id1, id2)


# correlation dataset (with temporal structure)
epitope_cor_ts$id1 <- as.numeric(epitope_cor_ts$id1)
epitope_cor_ts$id2 <- as.numeric(epitope_cor_ts$id2)
epitope_cor_ts[epitope_cor_ts$id1 > epitope_cor_ts$id2, c("id1", "id2")] <- 
  epitope_cor_ts[epitope_cor_ts$id1 > epitope_cor_ts$id2, c("id2", "id1")]
epitope_cor_ts <- epitope_cor_ts %>% dplyr::arrange(id1, id2)




# combine multiple datasets -----------------------------------------------

# add an additional measure in pairwise sequence alignment:
# longest strike of consecutive matching amino acid residuals
epitope_align <- epitope_align %>% 
  dplyr::mutate(
    match_seq_length = sapply(X = string_compare,
                              FUN = function(x){
                                max(nchar(unlist(stringr::str_extract_all(string = x, pattern = "[A-Z]+"))))
                              })) %>% 
  dplyr::mutate(match_seq_length = ifelse(match_seq_length == -Inf, 0, match_seq_length))


# epitope information: add protein information and enrichment frequency
epitope_info <- epitope_info %>% 
  dplyr::left_join({subset_epitope %>% 
      dplyr::select(id, UniProt_acc, Protein_name, start_position = start)},
      by = "id") %>% 
  dplyr::left_join({subset_freq %>% 
      dplyr::mutate(id = as.integer(id)) %>% 
      dplyr::arrange(id) %>% 
      dplyr::group_by(id) %>% 
      dplyr::summarise(freq_mean  = mean(freq_hit),
                       freq_total = sum(freq_hit * n) / sum(n),
                       freq_max   = max(freq_hit)) %>% 
      ungroup()},
      by = "id")

# epitope pairs: combine sequence similarity, z-score correlation, and cooccurrence
epitope_pair <- epitope_align %>% 
  dplyr::select(-subject_score, -pattern_score) %>% 
  dplyr::left_join({epitope_cooccur %>% 
      dplyr::select(subject_id = id1, 
                    pattern_id = id2, 
                    # TT, TF, FT, FF, oddsr, p, p_adj,
                    phi, mean_prop)}, 
      by = c("subject_id", "pattern_id")) %>% 
  dplyr::left_join({epitope_cor %>% 
      dplyr::select(subject_id = id1, 
                    pattern_id = id2, 
                    # t, p_adjust,
                    cor)}, 
      by = c("subject_id", "pattern_id"))%>% 
  dplyr::left_join({epitope_cor_ts %>% 
      dplyr::select(subject_id = id1, 
                    pattern_id = id2,  
                    # t, p_adjust,
                    cor_ts = cor)}, 
      by = c("subject_id", "pattern_id"))






# select a small subset for shiny development -----------------------------

# number of epitopes in different genus-species
genus_species <- as.data.frame(table(epitope_info$genus, epitope_info$species)) %>% 
  dplyr::rename(genus = Var1, species = Var2, n = Freq) %>% 
  dplyr::filter(n > 0) %>% 
  dplyr::arrange(genus)

# randomly select 10 epitopes from nine selected viral genus
set.seed(419)
epitope_info <- epitope_info %>% 
  dplyr::filter(genus %in% c("Alphacoronavirus",  # coronavirus
                             "Betacoronavirus",   # coronavirus
                             "Cytomegalovirus",   # CMV
                             "Enterovirus",       # enterovirus
                             "Flavivirus",        # dengue etc.
                             "Morbillivirus",     # measles
                             "Lymphocryptovirus", # Epstein-Barr virus
                             "Rubivirus",         # rubella
                             "Rubulavirus")) %>%  # mumps
  dplyr::group_by(genus) %>% 
  dplyr::slice_sample(n = 10) %>% 
  dplyr::arrange(id)

id_to_select <- epitope_info$id


# filter the pairwise analysis data to include only the selected ids
epitope_pair <- epitope_pair %>% 
  dplyr::filter(subject_id %in% id_to_select,
                pattern_id %in% id_to_select)



# save the objects for shiny ----------------------------------------------

save(epitope_info, epitope_pair,
     file = "shiny_network_visualization_data.RData")
