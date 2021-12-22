###
# Interactive network visualization of viral peptides in PhIP-Seq libraries
# show peptides' pairwise sequence similarity, correlation, and cooccurrence
#
# This is a R shiny app. To launch the app, click "Run App". 
#
# Contributor:
# Siyang Xia
# Jennifer L. Remmel
# Daniel Monaco
# H. Benjamin Larman
#
# Version: 2021-12-22
###



# 1. Preparation ----------------------------------------------------------


# a) packages -------------------------------------------------------------

# install the packages below if not already installed
load_lib <- c("shiny", "shinythemes", "shinyWidgets", 
              "tidyverse", "here", "devtools", "tools", 
              "rBLAST", "seqinr", 
              "igraph", "ggnetwork", "intergraph", 
              "cowplot", "plotly", "htmlwidgets")
install_lib <- load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)

# install virlink (a customized package for pairwise peptide analysis)
if(!("virlink" %in% installed.packages())){
  devtools::install_github(repo = "siyangxia419/virlink", 
                           ref = "main",
                           upgrade = "never",
                           auth_token = "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")
}
load_lib <- c(load_lib, "virlink")

# load all packages
for(lib in load_lib) library(lib, character.only = TRUE)
rm(install_lib, lib)



# b) data and global variables --------------------------------------------

# load the pairwise analysis dataset
load("VRC_full_data.RData")

# manual color palette
manualPalette = c("#0072B2", "#56B4E9", "#009E73", "#F0E442", "#D55E00")

# virus taxa and protein UniProt accession numbers
taxa_protein <- VRC_peptide_info %>% 
  dplyr::select(starts_with("taxon_"), UniProt_acc) %>% 
  dplyr::distinct()



# c) function for blastp ----------------------------------------------------------------------

blastp_dm <- function(pep_dt, 
                      fasta_dir = "C:/Users/siyang_xia/R/",
                      blastp_arg = "", 
                      other_info = TRUE){
  
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
                   BLAST_args = blastp_arg,
                   custom_format = "qseqid sseqid pident length evalue bitscore positive gaps ppos")
  
  
  # arrange the peptide ids
  predp <- predp %>% 
    mutate(qseqid = ordered(qseqid, levels = pep_dt$u_pep_id),
           sseqid = ordered(sseqid, levels = pep_dt$u_pep_id)) %>% 
    rename(id1 = qseqid, id2 = sseqid)
  
  predp[predp$id1 > predp$id2, c("id1", "id2")] <- 
    predp[predp$id1 > predp$id2, c("id2", "id1")]
  
  
  # add information of the peptides
  if(other_info){
    predp <- predp %>% 
      arrange(id1, id2) %>% 
      as_tibble() %>% 
      left_join({pep_dt %>% rename_with(~paste0("id1_", .))}, 
                by = c("id1" = "id1_u_pep_id")) %>% 
      left_join({pep_dt %>% rename_with(~paste0("id2_", .))}, 
                by = c("id2" = "id2_u_pep_id"))
  }
  
  return(predp)
}



# d) function to plot the network -----------------------------------------

epitope_network_visualization <- function(net_df, 
                                          color_var = "genus",
                                          color_title = "virus",
                                          color_pal = c("#0072B2", "#56B4E9", "#009E73", "#F0E442", "#D55E00"),
                                          edge_var = "sim_score",
                                          edge_title = "similarity score",
                                          edge_presentation = "color", 
                                          vertex_size_var = "freq",
                                          fig_title = "peptide network",
                                          interactive_plot = TRUE){
  
  # check if any edge should be plotted
  any_edge <- any(!(is.na(net_df$string_compare)))
  
  # vertex color
  ncolor <- length(unique(net_df[, color_var]))
  vertex_color_pal <- colorRampPalette(color_pal)(ncolor)  # color palette
  
  # range of the edge metric
  if(any_edge){
    edge_range <- range(pretty(range(net_df[, edge_var], na.rm = TRUE)))
  }
  
  # prepare for plotly
  if(interactive_plot){
    
    # add the peptide id of the end point
    net_df <- net_df %>% 
      dplyr::right_join({net_df %>% 
          select(xend = x, yend = y, nameend = name) %>% 
          distinct()},
          by = c("xend", "yend"))
    
    # add a marker point for each edge at the center to anker the hover text
    net_df <- net_df %>% 
      dplyr::mutate(x1 = ifelse(test = is.na(string_compare), 
                                yes = NA, 
                                no = (x + xend) / 2),
                    y1 = ifelse(test = is.na(string_compare), 
                                yes = NA, 
                                no = (y + yend) / 2),
                    text = ifelse(test = is.na(string_compare),
                                  yes  = paste(paste0("peptide id: ", name), 
                                               paste0("family    : ", family),
                                               paste0("genus     : ", genus),
                                               paste0("species   : ", species),
                                               paste0("strain    : ", organism),
                                               paste0("protein   : ", UniProt_acc),
                                               paste0("sequence  : ", pep_aa),
                                               paste0("enrichment frequency: ", round(freq, 3)),
                                               sep = "<br>"), 
                                  no   = paste(paste0("peptide 1: ", name),
                                               paste0("virus    : ", subject_organism),
                                               paste0("peptide 2: ", nameend),
                                               paste0("virus    : ", pattern_organism),
                                               paste0("alignment: ", string_compare), 
                                               paste0("sequence similarity: ", round(sim_score, 3)),
                                               paste0("Ab reactivity correlation: ", round(cor, 3)),
                                               paste0("Ab hit jaccard index: ", round(jaccard, 3)),
                                               sep = "<br>")))
  }
  
  
  # generate the static ggplot
  net_fig <- ggplot(net_df, aes(x = x, y = y, xend = xend, yend = yend))
  
  if(edge_presentation == "width"){  # similarity score shown by edge width:
    
    if(any_edge){
      net_fig <- net_fig +
        geom_edges(aes(size = get(edge_var)), color = "black", alpha = 0.5)+
        scale_size_continuous(name = edge_title,
                              limits = edge_range,
                              range = c(0, 1.5))
    }
    
    if(interactive_plot){
      net_fig <- net_fig +
        geom_nodes(aes(color = get(color_var), text = text), size = 2)
    }else{
      net_fig <- net_fig +
        geom_nodes(aes(color = get(color_var)), size = 2)
    }
    
    net_fig <- net_fig  +
      scale_color_manual(name = color_title,
                         values = vertex_color_pal)
    
    warning("the size aesthetic mapping is used to adjust edge width, 
            so peptide enrichment frequence is not reflected by point size.")
    
  }else if(edge_presentation == "color"){  # similarity score shown by gray scale
    
    if(any_edge){
      net_fig <- net_fig +
        geom_edges(aes(color = get(edge_var)), size = 1, alpha = 0.6) +
        scale_color_gradient(name = edge_title, 
                             low = "gray99", high = "black",
                             limits = edge_range) 
    }
    
    if(interactive_plot){
      net_fig <- net_fig +
        geom_nodes(aes(fill = get(color_var), text = text, size = get(vertex_size_var)), 
                   shape = 21, stroke = 0.1)
    }else{
      net_fig <- net_fig +
        geom_nodes(aes(fill = get(color_var), size = get(vertex_size_var)), 
                   shape = 21, stroke = 0.1)
    }
    
    net_fig <- net_fig +
      scale_fill_manual(name = color_title, 
                        values = vertex_color_pal) + 
      scale_size_continuous(name = vertex_size_title,
                            limits = c(0, 1), 
                            range = c(2, 5))
    
  }else{
    stop("invalid choice of edge presentation. please select from 'width' and 'color'.")
  }
  
  net_fig <- net_fig +
    ggtitle(fig_title) +
    theme_blank() +
    theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 12),
          plot.margin  = unit(c(10, 10, 10, 10), units = "pt"),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 10)) +
    guides(size = FALSE, color = FALSE)
  
  # # extract the legend
  # net_legend <- cowplot::get_legend(net_fig)
  
  
  # plotly interactive figure
  if(interactive_plot){
    
    if(any_edge){
      net_fig_plotly <- {net_fig + 
          geom_point(aes(x = x1, y = y1, text = text), 
                     color = "gray75", alpha = 0.5, size = 1.5)} %>%
        ggplotly(tooltip = "text")
    }else{
      net_fig_plotly <- net_fig %>% 
        ggplotly(tooltip = "text")
    }
    
    return(net_fig_plotly)
    
  }else{
    return(net_fig)
  }
  
}





# 2. UI -------------------------------------------------------------------

ui <- fluidPage(
  
  # shiny theme
  theme = shinytheme("cerulean"),
  # shinythemes::themeSelector(),
  
  # title              
  titlePanel(h3("VirScan peptide network visualization", 
                h5("Developer: Siyang Xia, Daniel Monaco, H. Benjamin Larman")),
             windowTitle = "VirScan peptide network visualization tool"
  ),
  
  

  ## a) input panels ------------------------------------------------------------------------------

  fluidRow(
    
    ### upload and select peptides and sample Ab reactivity profiles -----
    column(3, 
           
      # allow the user to upload a file that contains peptide information
      fileInput(inputId = "peptide_upload",
                label = "Choose peptide information file",
                multiple = FALSE,
                accept = c(".csv", ".tsv")),
      
      # allow the user to select variables to perform further filtering of the peptides
      uiOutput("filter_variable"),
      
      uiOutput("filter1"), 
      
      br(),
            
      # allow the user to upload a file that contains PhIP-Seq antibody reactivity profile
      fileInput(inputId = "reactivity_upload",
                label = "Choose antibody reactivity file",
                multiple = FALSE,
                accept = c(".csv", ".tsv")),
      
      # select samples to look at
      selectizeInput(inputId = "select_sample", 
                     label = "Select samples",
                     choices = c("all", colnames(VRC_reactivity)[-1]),
                     multiple = TRUE, 
                     options = list(placeholder = "type sample namnes")),
      
      # threshold of z-scores or fold-changes that determines significant enrichment
      numericInput(inputId = "hit_thres", 
                   label = "Threshold of hit", 
                   min = 0, 
                   max = 1000, 
                   value = 1, 
                   step = 1),
      
      # filter of peptide frequency
      sliderInput(inputId = "frequency",
                  label = "Enrichment frequency",
                  min = 0,
                  max = 1,
                  value = c(0, 1), 
                  step = 0.01, 
                  ticks = FALSE, 
                  dragRange = TRUE),

      
    ),
    
    
    ### parameters for computing peptide pairwise alignments and correlations -----
    column(3,
      
      # input for the BLASTP alignment
      textInput(inputId = "local_dir",
                label = "Local directory for temporary files",
                value = "C:/Users/siyang_xia/R/", 
                placeholder = "local directory (do not allow space)"),
      
      textInput(inputId = "blastp_arg",
                label = "Argument for BLASTP",
                value = "-evalue 100 -max_hsps 1 -soft_masking false -word_size 7 -max_target_seqs 100000", 
                placeholder = "in the format of -arg1 value1 -arg2 value2"),
      
      br(),

      # input for the pairwiseAlignment with BiocGenerics (implemented with virlink)
      selectInput(inputId = "align_type", 
                  label = "Type of alignment",
                  choices = c("global", "local", "overlap", "global-local", "local-global"),
                  multiple = FALSE, 
                  selectize = FALSE,
                  selected = "local"),
      
      selectInput(inputId = "sub_matrix", 
                  label = "Amino acid substitution matrix",
                  choices = c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100",
                              "PAM30", "PAM40", "PAM70", "PAM120", "PAM250"),
                  multiple = FALSE, 
                  selectize = FALSE,
                  selected = "BLOSUM62"),
      
      numericInput(inputId = "gap_opening", 
                   label = "Panelty for gap opening", 
                   min = 0, 
                   max = 100, 
                   value = 10.0, 
                   step = 0.1),

      numericInput(inputId = "gap_extension", 
                   label = "Panelty for gap extension", 
                   min = 0, 
                   max = 100, 
                   value = 4.0, 
                   step = 0.1),
      
    ),
    
    
    ### filtering the pairwise calculation results -----
    column(3, 
      
      numericRangeInput(inputId   = "evalue", 
                        label     = "E-value ranges",
                        value     = c(0, 100),
                        separator = " - ",
                        min       = 0,
                        max       = 10000),
      
      br(),
           
      sliderInput(inputId = "nchar",
                  label = "Length of alignment",
                  min = 1,
                  max = 56,
                  value = c(7, 56), 
                  step = 1, 
                  ticks = FALSE, 
                  dragRange = TRUE),
      
      sliderInput(inputId = "matches",
                  label = "Number of matched amino acids",
                  min = 1,
                  max = 56,
                  value = c(5, 56),
                  step = 1, 
                  ticks = FALSE, 
                  dragRange = TRUE),
      
      sliderInput(inputId = "match_seq_length",
                  label = "Length of continuous matches",
                  min = 1,
                  max = 56,
                  value = c(3, 56), 
                  step = 1, 
                  ticks = FALSE, 
                  dragRange = TRUE),
      
      sliderInput(inputId = "sim_score",
                  label = "Sequence similarity",
                  min = 0,
                  max = 1,
                  value = c(0.1, 1), 
                  step = 0.01, 
                  ticks = FALSE, 
                  dragRange = TRUE),
      
      br(),
      
      sliderInput(inputId = "cor",
                  label = "Ab reactivity correlation",
                  min = -1,
                  max = 1,
                  value = c(-1, 1), 
                  step = 0.01, 
                  ticks = FALSE, 
                  dragRange = TRUE),
      
      sliderInput(inputId = "jaccard",
                  label = "Ab enrichment jaccard index",
                  min = 0,
                  max = 1,
                  value = c(0, 1), 
                  step = 0.01, 
                  ticks = FALSE, 
                  dragRange = TRUE),
      
      br(),
      
      # seed for determine the coordinate of vertices
      numericInput(inputId = "seed", 
                   label = "Seed", 
                   min = 1, 
                   max = 1000, 
                   value = 111, 
                   step = 1),
      
    ),
    
    
    ### plotting -----
    column(3, 
      
      # seed for determine the coordinate of vertices
      numericInput(inputId = "seed", 
                   label = "Seed", 
                   min = 1, 
                   max = 1000, 
                   value = 111, 
                   step = 1),
      
      selectInput(inputId = "color_var", 
                  label = "Color of nodes", 
                  choices = c("family", "genus", "species", "organism"), 
                  selected = "family",
                  multiple = FALSE)),
    
    ),
    
  
  fluidRow(
    
    column(3, offset = 1.5,
           actionButton(inputId = "load", label = "Submit")),
    
    column(3, offset = 1.5, 
           actionButton(inputId = "calculate", label = "Compute")),
    
    column(3, offset = 1.5,
           actionButton(inputId = "filter", label = "Filter")),
    
    column(3, offset = 1.5,
           actionButton(inputId = "plot", label = "Plot"))
  ),
  
  hr(),

  br(),
  

  
  # b) Output panels ----------------------------------------------------------------------------
  
  htmlOutput(outputId = "n_node"),
  
  br(),
  br(), 
  
  fluidRow(
    
    tabsetPanel(type = "tabs",
                tabPanel("peptide metadata",
                         dataTableOutput("peptide_info_table")),
                
                tabPanel("antibody reactivity",
                         dataTableOutput("ab_reactivity_table")),
                
                tabPanel("peptide pairwise sequence similarity",
                         dataTableOutput("pairwise_seq_sim")),
                
                tabPanel("peptide pairwise sequence BLASTP",
                         dataTableOutput("pairwise_blastp")),
                
                tabPanel("peptide pairwise antibody reactivity correlation",
                         dataTableOutput("pairwise_cor")),
                
                tabPanel("peptide pairwise jaccard index",
                         dataTableOutput("pairwise_jaccard")),
                
                tabPanel("peptide pairwise calculation",
                         dataTableOutput("pairwise_all")),
                
                tabPanel("sequence similarity", 
                         plotlyOutput("plotly_seq")),
                
                tabPanel("antibody reactivity correlation",
                         plotlyOutput("plotly_cor")),
                
                tabPanel("antibody hit jaccard index",
                         plotlyOutput("plotly_jaccard"))
    )
    
  )
  
)




# 3. Server ---------------------------------------------------------------


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # a) upload the peptide metadata and select variables to filter ---------

  # update peptide information data frame based on whether the user upload a file
  peptide_metadata <- reactive({
  
    if(is.null(input$peptide_upload)){  # no file upload, use the pre-loaded VRC data
      
      peptide_info <- VRC_peptide_info
      
    }else{                              # user updated a file
      
      req(input$peptide_upload)         
      
      peptide_upload_file <- input$peptide_upload
      ext <- tools::file_ext(peptide_upload_file$datapath)  # file extension
      
      # read the file according to the file extension
      if(ext == "csv"){
        peptide_info <- read_csv(file = peptide_upload_file$datapath,
                                 col_names = TRUE) 
      }else if(ext == "tsv"){
        peptide_info <- read_tsv(file = peptide_upload_file$datapath,
                                 col_names = TRUE)
      }else{
        stop("Please select from csv or tsv files.")
      }
      
    }
    return(peptide_info)
  })
  
  
  # update the UI depending on what columns exist in the peptide information data
  output$filter_variable <- renderUI({
    
    # a list of variables in peptide info
    filter_list <- names(peptide_metadata())
    filter_list <- filter_list[!(filter_list %in% c("u_pep_id", "pep_aa"))]
    
    default_filter <- grep(pattern = "taxon_", x = filter_list, value = TRUE)
    
    # create a multiple selection UI to allow the user to decide which variables to filter
    selectizeInput(inputId = "filter_var", 
                   label = "Vriables for further filtering",
                   choices = filter_list,
                   multiple = TRUE, 
                   options = list(placeholder = "type variable namnes"))
    
  })
  

  # create UI for filtering peptide
  output$filter1 <- renderUI({
    
    dynamic_ui <- lapply(
      
      X = input$filter_var, 
      
      FUN = function(x) {
        if(class(peptide_metadata()[[x]]) %in% c("integer", "numeric")){
          
          x_fullrange <- range(peptide_metadata()[[x]], na.rm = TRUE)
          
          # create a slider bar
          return(sliderInput(inputId = x,
                             label = x,
                             min = x_fullrange[1],
                             max = x_fullrange[2],
                             value = x_fullrange, 
                             ticks = FALSE, 
                             dragRange = TRUE))
          
        }else{
          
          x_fulllist <- sort(unique(as.character(peptide_metadata()[[x]])))
          
          # create a multiple selection UI
          return(selectizeInput(inputId = x, 
                                label = x,
                                choices = c("all", x_fulllist),
                                multiple = TRUE, 
                                options = list(placeholder = "type variable options",
                                               delimiter = " ", 
                                               create = T)))
          
        }
        
      })
    
    return(dynamic_ui)
    
  })
  
  
  
  # b) upload antibody reactivity profile and select samples ------------------------------------
  
  # update peptide information data frame based on whether the user upload a file
  ab_data <- reactive({
    
    if(is.null(input$reactivity_upload)){  # no file upload, use the pre-loaded VRC data
      
      ab <- VRC_reactivity
      
    }else{                              # user updated a file
      
      req(input$peptide_upload)         
      
      reactivity_upload_file <- input$reactivity_upload
      ext <- tools::file_ext(reactivity_upload_file$datapath)  # file extension
      
      # read the file according to the file extension
      if(ext == "csv"){
        ab <- read_csv(file = reactivity_upload_file$datapath,
                       col_names = TRUE) 
      }else if(ext == "tsv"){
        ab <- read_tsv(file = reactivity_upload_file$datapath,
                       col_names = TRUE)
      }else{
        stop("Please select from csv or tsv files.")
      }
      
    }
    
    return(ab)
  })
  
  
  observeEvent(input$reactivity_upload, 
               {
                 
                 temp <- names(ab_data())[-1]
                 updateSelectizeInput(session,
                                      inputId = "select_sample",
                                      choices = c("all", temp),
                                      selected = NULL,
                                      options = list(delimiter = " ", create = T))
               })
  
  

  # c) filter the peptides and samples ----------------------------------------------------------

  ### filter the peptide information based on the selected variables
  peptide_info_filtered <- reactiveValues()
  
  observeEvent(input$load, {
    
    # filter peptide metadata according to the selected variables
    var_list <- input$filter_var
    
    dt <- peptide_metadata()

    for(v in var_list){

      if(class(dt[[v]]) %in% c("integer", "numeric")){

        dt <- dt %>%
          dplyr::filter(between(get(v), input[[v]][1], input[[v]][2]))
        
      }else{
        
        if(!("all" %in% input[[v]])){
          dt <- dt %>%
            dplyr::filter(as.character(get(v)) %in% input[[v]])
        }
        
      }

    }
    
    # further filter peptides according to their prevalence in the selected samples
    sample_list    <- input$select_sample
    if("all" %in% sample_list | is.null(sample_list)) sample_list <- names(ab_data())[-1]
    
    hit_thres      <- input$hit_thres
    freq_threshold <- input$frequency
    
    hit_freq <- ab_data() %>% 
      select(u_pep_id, all_of(sample_list)) %>% 
      filter(u_pep_id %in% dt$u_pep_id) %>% 
      pivot_longer(cols = all_of(sample_list), names_to = "sample", values_to = "original") %>% 
      mutate(hit = as.integer(original > hit_thres)) %>% 
      group_by(u_pep_id) %>% 
      summarise(n_sample = n(),
                n_hit    = sum(hit),
                .groups  = "drop") %>% 
      mutate(frequency = n_hit / n_sample)
    
    dt <- dt %>% 
      left_join(hit_freq, by = "u_pep_id") %>% 
      filter(frequency >= freq_threshold[1], frequency <= freq_threshold[2])
    
    peptide_info_filtered$dt <- dt

  })
  
  
  ### filter the peptide information based on the selected variables
  ab_reactivity_filtered <- reactiveValues()
  
  observeEvent(input$load, {
    
    sample_list <- input$select_sample
    if("all" %in% sample_list | is.null(sample_list)) sample_list <- names(ab_data())[-1]
    
    dt <- ab_data() %>% 
      select(u_pep_id, all_of(sample_list))
    
    ab_reactivity_filtered$dt <- dt
    
  })
  
  
  

  # d) print the number of peptides and samples selected ----------------------------------------

  output$n_node <- renderUI({
    
    req(peptide_info_filtered$dt)
    req(ab_reactivity_filtered$dt)
    
    n_pep <- nrow(peptide_info_filtered$dt)
    n_pep2 <- nrow(ab_reactivity_filtered$dt)
    n_sam <- ncol(ab_reactivity_filtered$dt) - 1
    
    same_peptides <- setequal(peptide_info_filtered$dt$u_pep_id, ab_reactivity_filtered$dt$u_pep_id)
    
    if(n_pep == n_pep2 & same_peptides){
      text_to_print <- paste(paste("Number of peptides:", n_pep),
                             paste("Number of samples:", n_sam),
                             sep = '<br/>')
    }else{
      text_to_print <- paste(paste("Number of peptides with metadata:", n_pep),
                             paste("Number of peptides with Ab reactivity:", n_pep2), 
                             paste("Number of samples:", n_sam),
                             "Warning: Peptides in the metadata file and Ab reactivity file do not match!",
                             sep = '<br/>')
    }
    
    if(n_pep > 50 | n_pep2 > 50){
      text_to_print <- paste(text_to_print, 
                             "Warning: more than 50 peptides in the dataset.
                             Calculation may take a long time and the network may be too dense.",
                             sep = '<br/>')
    }
    
    
    return(HTML(text_to_print))
  })
  
  output$peptide_info_table <- renderDataTable(peptide_info_filtered$dt)
  
  output$ab_reactivity_table <- renderDataTable(ab_reactivity_filtered$dt)
  
  
  
  # e) calculate the sequence similarity and reactivity correlation -----------------------------
  peptide_pairwise <- reactiveValues()
  
  observeEvent(input$calculate, {
    
    withProgress(message = "calculating", value = 0, { 
    
    peptide_id_list <- peptide_info_filtered$dt$u_pep_id
    
    setProgress(value = 0.05, message = "start computing")
    
    ### sequence similarity analysis using the "virlink" package
    peptide_seq_sim <- peptide_pairwise_alignment(
      peptides        = peptide_info_filtered$dt,
      id_col          = "u_pep_id",
      seq_col         = "pep_aa",
      sub_matrix      = input$sub_matrix, 
      gap_opening     = input$gap_opening,
      gap_extension   = input$gap_extension,
      align_type      = input$align_type,
      self_comparison = TRUE,    # default
      full_align      = FALSE,   # default
      other_info      = TRUE,    # default
      parallel_ncore  = NULL,    # default
      output_str      = "tibble")
    
    self_align <- peptide_seq_sim %>% filter(id1 == id2)
    
    peptide_seq_sim <- peptide_seq_sim %>% 
      filter(id1 != id2) %>% 
      left_join({self_align %>% select(id1, id1_score = score)}, by = "id1") %>% 
      left_join({self_align %>% select(id2, id2_score = score)}, by = "id2") %>% 
      mutate(opt_score        = (id1_score + id2_score) / 2,
             sim_score        = score / opt_score,
             match_seq_length = sapply(X = string_compare,
                                       FUN = function(x){
                                         max(nchar(unlist(str_extract_all(string = x, pattern = "[A-Z]+"))))
                                       })) %>% 
      mutate(match_seq_length = ifelse(match_seq_length == -Inf, 0, match_seq_length)) %>% 
      select(-id1_score, -id2_score) %>% 
      arrange(id1, id2) %>% 
      as_tibble()
    
    peptide_pairwise$seq <- peptide_seq_sim
    
    setProgress(value = 0.25, message = "finished peptide alignment")
    
    ### sequence similarity analysis using BLAST
    peptide_blastp <- blastp_dm(pep_dt     = peptide_info_filtered$dt, 
                                fasta_dir  = input$local_dir, 
                                blastp_arg = input$blastp_arg,
                                other_info = TRUE)
    
    peptide_blastp <- peptide_blastp %>% filter(id1 != id2)
    
    peptide_pairwise$blastp <- peptide_blastp
    
    setProgress(value = 0.5, message = "finished BLASTP")
    
    
    ### antibody reactivity correlation
    ab_reactivity_formatted <- ab_reactivity_filtered$dt %>% 
      column_to_rownames(var = "u_pep_id") %>% 
      t() %>% 
      as.data.frame()
    
    ab_cor <- peptide_pairwise_correlation(
      d             = ab_reactivity_formatted, 
      analysis_type = "correlation",
      perform_test  = FALSE,
      cor_method    = "pearson",
      output_str    = "tibble")
    
    peptide_pairwise$cor <- ab_cor
    
    setProgress(value = 0.7, message = "finished Ab correlation calculation")
    
    
    ### antibody reactivity jaccard index
    ab_hit_formatted <- as.data.frame((ab_reactivity_formatted > input$hit_thres) * 1)
    
    ab_jaccard <- peptide_pairwise_correlation(
      d               = ab_hit_formatted, 
      analysis_type   = "cooccurrence",
      perform_test    = FALSE,
      occ_method      = "jaccard",
      hit_threshold   = 1,
      output_str      = "tibble")
    
    peptide_pairwise$jac <- ab_jaccard
    
    setProgress(value = 0.9, message = "finished Ab jaccard calculation")

    
    
    ### check if the peptides in the peptide metadata and the antibody reactivity are the same
    same_peptides <- setequal(peptide_info_filtered$dt$u_pep_id, ab_reactivity_filtered$dt$u_pep_id)
    peptide_pairwise$same <- same_peptides
    
    
    ### match the id1 and id2 between sequence alignment, ab reactivity correlation and jaccard index
    
    # if the peptides are different between the two inputs
    if(!same_peptides){
      peptide_id_list <- peptide_id_list[peptide_id_list %in% intersect(peptide_info_filtered$dt$u_pep_id, 
                                                                        ab_reactivity_filtered$dt$u_pep_id)]
      
      peptide_seq_sim <- peptide_seq_sim %>% 
        filter(id1 %in% peptide_id_list, id2 %in% peptide_id_list)
      
      ab_cor <- ab_cor %>% 
        filter(id1 %in% peptide_id_list, id2 %in% peptide_id_list)
      
      ab_jaccard <- ab_jaccard %>% 
        filter(id1 %in% peptide_id_list, id2 %in% peptide_id_list)
    }
    
    
    # peptide alignment:
    peptide_seq_sim <- peptide_seq_sim %>% 
      mutate(id1 = ordered(id1, levels = peptide_id_list),
             id2 = ordered(id2, levels = peptide_id_list))
    
    peptide_seq_sim[peptide_seq_sim$id1 > peptide_seq_sim$id2, c("id1", "id2")] <- 
      peptide_seq_sim[peptide_seq_sim$id1 > peptide_seq_sim$id2, c("id2", "id1")]
    
    peptide_seq_sim <- peptide_seq_sim %>% arrange(id1, id2)
    
    # ab reactivity correlation:
    ab_cor <- ab_cor %>% 
      mutate(id1 = ordered(id1, levels = peptide_id_list),
             id2 = ordered(id2, levels = peptide_id_list))
    
    ab_cor[ab_cor$id1 > ab_cor$id2, c("id1", "id2")] <- 
      ab_cor[ab_cor$id1 > ab_cor$id2, c("id2", "id1")]
    
    ab_cor <- ab_cor %>% arrange(id1, id2)
    
    # ab jaccard index: 
    ab_jaccard <- ab_jaccard %>% 
      mutate(id1 = ordered(id1, levels = peptide_id_list),
             id2 = ordered(id2, levels = peptide_id_list))
    
    ab_jaccard[ab_jaccard$id1 > ab_jaccard$id2, c("id1", "id2")] <- 
      ab_jaccard[ab_jaccard$id1 > ab_jaccard$id2, c("id2", "id1")]
    
    ab_jaccard <- ab_jaccard %>% arrange(id1, id2)
    
    
    ### combine the three measurement
    peptide_pairwise$all <- peptide_seq_sim %>% 
      dplyr::left_join(peptide_blastp, 
                       by = intersect(colnames(peptide_seq_sim), colnames(peptide_blastp))) %>% 
      dplyr::full_join(ab_cor, by = c("id1", "id2")) %>% 
      dplyr::full_join(ab_jaccard, by = c("id1", "id2"))
    
    })
    
  })
  


  # f) further filtering the pairwise calculations ----------------------------------------------
  
  pairwise_filter <- reactiveValues()
  
  observeEvent(input$filter, {
    
    req(peptide_pairwise)
    
    pairwise_filter$same <- peptide_pairwise$same
    
    # filter the BLASTP results
    pairwise_filter$blastp <- peptide_pairwise$blastp %>% 
      filter(between(evalue, input$evalue[1], input$evalue[2]))
    
    # filter the sequence alignment results
    pairwise_filter$seq <- peptide_pairwise$seq %>% 
      filter(between(nchar, input$nchar[1], input$nchar[2])) %>% 
      filter(between(matches, input$matches[1], input$matches[2])) %>% 
      filter(between(match_seq_length, input$match_seq_length[1], input$match_seq_length[2])) %>%
      filter(between(sim_score, input$sim_score[1], input$sim_score[2]))
    
    # filter the Ab correlation
    pairwise_filter$cor <- peptide_pairwise$cor %>% 
      filter(between(cor, input$cor[1], input$cor[2]))
      
    # filter the Ab jaccard index
    pairwise_filter$jac <- peptide_pairwise$jac %>% 
      filter(between(jaccard, input$jaccard[1], input$jaccard[2]))
    
    # combine the results
    pairwise_filter$all <- pairwise_filter$seq %>% 
      left_join(pairwise_filter$blastp, 
                by = intersect(colnames(pairwise_filter$seq),
                               colnames(pairwise_filter$blastp))) %>% 
      full_join(pairwise_filter$cor, by = c("id1", "id2")) %>% 
      full_join(pairwise_filter$jac, by = c("id1", "id2"))
    
  })
  
  
  # show the pairwise calculation results after filtering
  output$pairwise_seq_sim  <- renderDataTable({
    req(pairwise_filter)
    pairwise_filter$seq
  })
  
  output$pairwise_blastp  <- renderDataTable({
    req(pairwise_filter)
    pairwise_filter$blastp
  })
  
  output$pairwise_cor     <- renderDataTable({
    req(pairwise_filter)
    pairwise_filter$cor
  })
  
  output$pairwise_jaccard <- renderDataTable({
    req(pairwise_filter)
    pairwise_filter$jac
  })
  
  output$pairwise_all     <- renderDataTable({
    req(pairwise_filter)
    pairwise_filter$all
  })
  

  

  # g) network construction ---------------------------------------------------------------------

  network_dt <- reactive({

    # # weight to group species and genus together
    # weight_same_family <- 1
    # weight_same_genus <- 1
    # weight_same_species <- 1
    
    net_vertices <- peptide_info_filtered$dt
    
    ## sequence similarity network -----

    # construct an igraph object
    seq_sim_net <- igraph::graph_from_data_frame(d = peptide_pairwise$seq,
                                                 vertices = peptide_info_filtered$dt,
                                                 directed = FALSE)

    # # weight of edges
    # E(seq_sim_net)$weight <- E(seq_sim_net)$sim_score +
    #   E(seq_sim_net)$same_family  * weight_same_family +
    #   E(seq_sim_net)$same_genus   * weight_same_genus +
    #   E(seq_sim_net)$same_species * weight_same_species

    # fortify the igraph to a data frame suitable for ggplot2
    set.seed(input$seed)
    suppressWarnings({
      seq_sim_net_df <- ggnetwork::ggnetwork(x = seq_sim_net,
                                             layout = igraph::with_fr())
    })
    
    # filter the edges
    seq_sim_net_df <- seq_sim_net_df %>% 
      filter(between(nchar, input$nchar[1], input$nchar[2])) %>% 
      filter(between(matches, input$matches[1], input$matches[2])) %>% 
      filter(between(match_seq_length, input$match_seq_length[1], input$match_seq_length[2])) %>%
      filter(between(sim_score, input$sim_score[1], input$sim_score[2]))
    
    
    
    ## BLASTP network -----
    
    # all pairs of peptides
    blastp_edge <- peptide_pairwise$seq %>% 
      left_join(pairwise_filter$blastp, 
                by = intersect(colnames(pairwise_filter$seq),
                               colnames(pairwise_filter$blastp)))
    
    # construct an igraph object
    blastp_net <- igraph::graph_from_data_frame(d = blastp_edge,
                                                vertices = peptide_info_filtered$dt,
                                                directed = FALSE)
    
    # fortify the igraph to a data frame suitable for ggplot2
    set.seed(input$seed)
    suppressWarnings({
      blastp_net_df <- ggnetwork::ggnetwork(x = blastp_net,
                                            layout = igraph::with_fr())
    })
    
    # filter the edges
    blastp_net_df <- blastp_net_df %>% 
      filter(!is.na(evalue)) %>% 
      filter(between(evalue, input$evalue[1], input$evalue[2]))
    
    
    
    ## Ab reactivity correlation network -----
    
    # construct an igraph object
    seq_sim_net <- igraph::graph_from_data_frame(d = peptide_pairwise$cor,
                                                 vertices = peptide_info_filtered$dt,
                                                 directed = FALSE)
    
    # fortify the igraph to a data frame suitable for ggplot2
    set.seed(input$seed)
    suppressWarnings({
      seq_sim_net_df <- ggnetwork::ggnetwork(x = seq_sim_net,
                                             layout = igraph::with_fr())
    })
    
    # filter the edges
    seq_sim_net_df <- seq_sim_net_df %>% 
      filter(between(nchar, input$nchar[1], input$nchar[2])) %>% 
      filter(between(matches, input$matches[1], input$matches[2])) %>% 
      filter(between(match_seq_length, input$match_seq_length[1], input$match_seq_length[2])) %>%
      filter(between(sim_score, input$sim_score[1], input$sim_score[2]))

  })


  # observeEvent(input$filter_var,
  #              {
  #                print(input)
  #                # updateSelectizeInput(session,
  #                #                      inputId = input[[input$filter_var[1]]],
  #                #                      choices = c("all"),
  #                #                      selected = NULL,
  #                #                      server = TRUE)
  #              })
  
  # 
  # 
  # 
  # # select virus family
  # updateSelectizeInput(session, 
  #                      inputId = "family", 
  #                      choices = c("all", unique(taxa_protein$family)), 
  #                      selected = NULL,
  #                      server = TRUE,
  #                      options = list(delimiter = " ", create = T))
  # 
  # # select virus genus
  # observeEvent(input$family,
  #              {
  #                # update the genus choices according to virus family input
  #                temp <- taxa_protein %>% 
  #                {if(!("all" %in% input$family)){
  #                  dplyr::filter(., family %in% input$family)
  #                }else{
  #                  .
  #                }}
  #                
  #                updateSelectizeInput(session, 
  #                                     inputId = "genus", 
  #                                     choices = c("all", unique(temp$genus)), 
  #                                     selected = NULL,
  #                                     server = TRUE)
  #              })
  # 
  # # select virus species
  # observeEvent(input$genus,
  #              {
  #                # update the species choices according to virus genus input
  #                temp <- taxa_protein %>% 
  #                {if(!("all" %in% input$family)){
  #                  dplyr::filter(., family %in% input$family)
  #                }else{
  #                  .
  #                }} %>% 
  #                {if(!("all" %in% input$genus)){
  #                  dplyr::filter(., genus %in% input$genus)
  #                }else{
  #                  .
  #                }}
  #                
  #                updateSelectizeInput(session, 
  #                                     inputId = "species", 
  #                                     choices = c("all", unique(temp$species)), 
  #                                     selected = NULL,
  #                                     server = TRUE)
  #              })
  # 
  # # select virus organism
  # observeEvent(input$species,
  #              {
  #                # update the organism choices according to virus species input
  #                temp <- taxa_protein %>% 
  #                {if(!("all" %in% input$family)){
  #                  dplyr::filter(., family %in% input$family)
  #                }else{
  #                  .
  #                }} %>% 
  #                {if(!("all" %in% input$genus)){
  #                  dplyr::filter(., genus %in% input$genus)
  #                }else{
  #                  .
  #                }} %>% 
  #                {if(!("all" %in% input$species)){
  #                  dplyr::filter(., species %in% input$species)
  #                }else{
  #                  .
  #                }}
  #                
  #                updateSelectizeInput(session, 
  #                                     inputId = "organism", 
  #                                     choices = c("all", unique(temp$organism)), 
  #                                     selected = NULL,
  #                                     server = TRUE)
  #              })
  # 
  # # select protein
  # observeEvent(input$organism,
  #              {
  #                # update the protein choices according to virus organism input
  #                temp <- taxa_protein %>% 
  #                {if(!("all" %in% input$family)){
  #                  dplyr::filter(., family %in% input$family)
  #                }else{
  #                  .
  #                }} %>% 
  #                {if(!("all" %in% input$genus)){
  #                  dplyr::filter(., genus %in% input$genus)
  #                }else{
  #                  .
  #                }} %>% 
  #                {if(!("all" %in% input$species)){
  #                  dplyr::filter(., species %in% input$species)
  #                }else{
  #                  .
  #                }} %>% 
  #                {if(!("all" %in% input$organism)){
  #                  dplyr::filter(., organism %in% input$organism)
  #                }else{
  #                  .
  #                }}
  #                
  #                updateSelectizeInput(session, 
  #                                     inputId = "uniprot", 
  #                                     choices = c("all", unique(temp$UniProt_acc)), 
  #                                     selected = NULL,
  #                                     server = TRUE)
  #              })
  
  

  
  # # build the network
  # network_dt <- eventReactive(input$load, {
  #   
  #   var_list <- input$filter_var
  #   
  #   ### vextices
  #   if("all" %in% input$family){  # filter family
  #     vertex_d <- epitope_info
  #   }else{
  #     vertex_d <- epitope_info %>% 
  #       dplyr::filter(family %in% input$family)
  #   }
  #   
  #   if("all" %in% input$genus){  # filter genus
  #     vertex_d <- vertex_d
  #   }else{
  #     vertex_d <- vertex_d %>% 
  #       dplyr::filter(genus %in% input$genus)
  #   }
  #   
  #   if("all" %in% input$species){  # filter species
  #     vertex_d <- vertex_d
  #   }else{
  #     vertex_d <- vertex_d %>% 
  #       dplyr::filter(species %in% input$species)
  #   }
  #   
  #   if("all" %in% input$organism){  # filter virus organism
  #     vertex_d <- vertex_d
  #   }else{
  #     vertex_d <- vertex_d %>% 
  #       dplyr::filter(organism %in% input$organism)
  #   }
  #   
  #   if("all" %in% input$uniprot){  # filter protein UniProt accession number
  #     vertex_d <- vertex_d
  #   }else{
  #     vertex_d <- vertex_d %>% 
  #       dplyr::filter(UniProt_acc %in% input$uniprot)
  #   }
  #   
  #   if("freq" %in% names(epitope_info)){  # filter by peptide's enrichment frequency (column "freq")
  #     vertex_d <- vertex_d %>% 
  #       dplyr::filter(freq >= input$frequency[1], freq <= input$frequency[2])
  #   }else{
  #     vertex_d <- vertex_d %>% 
  #       dplyr::mutate(freq = 1)
  #   }
  #   
  #   
  #   ### edges
  #   edge_d <- epitope_pair %>% 
  #     dplyr::filter(subject_id %in% vertex_d$id,
  #                   pattern_id %in% vertex_d$id)
  #   
  #   # construct an igraph object
  #   net <- igraph::graph_from_data_frame(d = edge_d, 
  #                                        vertices = vertex_d, 
  #                                        directed = FALSE)
  #   
  #   # weight to group species and genus together
  #   weight_same_family <- 1
  #   weight_same_genus <- 1
  #   weight_same_species <- 1
  #   
  #   # weight of edges
  #   E(net)$weight <- E(net)$sim_score + 
  #     E(net)$same_family  * weight_same_family + 
  #     E(net)$same_genus   * weight_same_genus + 
  #     E(net)$same_species * weight_same_species
  #   
  #   # fortify the igraph to a data frame suitable for ggplot2
  #   set.seed(input$seed)
  #   suppressWarnings({
  #     net_fig_df <- ggnetwork::ggnetwork(x = net, 
  #                                        layout = igraph::with_fr())
  #   })
  #   
  #   return(net_fig_df)
  # })
  
  
  
  # filter the network
  network_flt <- eventReactive(input$filter, {
    
    if(nrow(network_dt()) > 0){
      net_fig_df <- network_dt() %>% 
        dplyr::filter((nchar >= input$nchar[1] & 
                       nchar <= input$nchar[2]) | 
                       is.na(nchar)) %>% 
        dplyr::filter((matches >= input$matches[1] & 
                       matches <= input$matches[2]) | 
                       is.na(matches)) %>% 
        dplyr::filter((match_seq_length >= input$match_seq_length[1] &
                       match_seq_length <= input$match_seq_length[2]) | 
                       is.na(match_seq_length)) %>% 
        dplyr::filter((sim_score >= input$sim_score[1] & 
                       sim_score <= input$sim_score[2]) | 
                       is.na(sim_score)) %>% 
        dplyr::filter((cor >= input$cor[1] & 
                       cor <= input$cor[2]) | 
                       is.na(cor)) %>% 
        dplyr::filter((jaccard >= input$jaccard[1] & 
                       jaccard <= input$jaccard[2]) | 
                       is.na(jaccard))
      
      # filter epitope pairs by their virus family
      if(input$same_family == "yes"){
        net_fig_df <- net_fig_df %>% dplyr::filter(same_family == TRUE | is.na(same_family))
      }else if(input$same_family == "no"){
        net_fig_df <- net_fig_df %>% dplyr::filter(same_family == FALSE | is.na(same_family))
      }
      
      # filter epitope pairs by their virus genus
      if(input$same_genus == "yes"){
        net_fig_df <- net_fig_df %>% dplyr::filter(same_genus == TRUE | is.na(same_genus))
      }else if(input$same_genus == "no"){
        net_fig_df <- net_fig_df %>% dplyr::filter(same_genus == FALSE | is.na(same_genus))
      }
      
      # filter epitope pairs by their virus species
      if(input$same_species == "yes"){
        net_fig_df <- net_fig_df %>% dplyr::filter(same_species == TRUE | is.na(same_species))
      }else if(input$same_species == "no"){
        net_fig_df <- net_fig_df %>% dplyr::filter(same_species == FALSE | is.na(same_species))
      }
      
      # filter epitope pairs by their virus organism
      if(input$same_organism == "yes"){
        net_fig_df <- net_fig_df %>% dplyr::filter(same_organism == TRUE | is.na(same_organism))
      }else if(input$same_organism == "no"){
        net_fig_df <- net_fig_df %>% dplyr::filter(same_organism == FALSE | is.na(same_organism))
      }
      
      # filter tiled epitopes
      if(input$tiling == "yes"){
        net_fig_df <- net_fig_df %>% dplyr::filter(tiling == TRUE | is.na(tiling))
      }else if(input$tiling == "no"){
        net_fig_df <- net_fig_df %>% dplyr::filter(tiling == FALSE | is.na(tiling))
      }

    }else{
      net_fig_df <- network_dt()
    }
    
    # # A message if the number of nodes and edges are too many
    # session$sendCustomMessage(type = 'testmessage',
    #                           message = 'Thank you for clicking')
    
    return(net_fig_df)
  })
  
  
  # visualize the network: sequence similarity
  output$plotly_seq <- renderPlotly({
    
    req(network_flt())
    
    if(nrow(network_flt()) > 0){
      
      epitope_network_visualization(net_df = network_flt(),
                                    color_var = input$color_var,
                                    color_title = "",
                                    color_pal = manualPalette,
                                    edge_var = "sim_score",
                                    edge_title = "",
                                    edge_presentation = "color",
                                    vertex_size_var = "freq",
                                    fig_title = "network of sequence similarity",
                                    interactive_plot = TRUE)
      
    }else{
      plotly_empty(type = "scatter", mode = "markers")
    }
    
  })
  
  
  # visualize the network: co-occurrence jaccard index
  output$plotly_jaccard <- renderPlotly({
    
    req(network_flt())
    
    if(nrow(network_flt()) > 0){
      
      epitope_network_visualization(net_df = network_flt(),
                                    color_var = input$color_var,
                                    color_title = "",
                                    color_pal = manualPalette,
                                    edge_var = "jaccard",
                                    edge_title = "",
                                    edge_presentation = "color",
                                    vertex_size_var = "freq",
                                    vertex_size_title = "",
                                    fig_title = "network of antibody hit co-occurrence",
                                    interactive_plot = TRUE)
    }else{
      plotly_empty(type = "scatter", mode = "markers")
    }
  })
  
  
  # visualize the network: correlation
  output$plotly_cor <- renderPlotly({
    
    req(network_flt())
    
    if(nrow(network_flt()) > 0){
      
      epitope_network_visualization(net_df = network_flt(),
                                    color_var = input$color_var,
                                    color_title = "",
                                    color_pal = manualPalette,
                                    edge_var = "cor",
                                    edge_title = "",
                                    edge_presentation = "color",
                                    vertex_size_var = "freq",
                                    vertex_size_title = "",
                                    fig_title = "network of antibody reactivity correlation",
                                    interactive_plot = TRUE)
    }else{
      plotly_empty(type = "scatter", mode = "markers")
    }
  })
  
}




# 4. Run the application --------------------------------------------------
shinyApp(ui = ui, server = server)
