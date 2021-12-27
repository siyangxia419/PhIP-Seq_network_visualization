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
load_lib <- c("shiny", "shinythemes", "shinyWidgets", "DT", 
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
manualPalette = c("#0072B2", "#009E73", "#F0E442", "#D55E00")

# virus taxa and protein UniProt accession numbers
taxa_protein <- VRC_peptide_info %>% 
  dplyr::select(starts_with("taxon_"), UniProt_acc) %>% 
  dplyr::distinct()

options(DT.options = list(
  paging = TRUE,
  searching = TRUE,
  fixedColumns = TRUE,
  autoWidth = TRUE,
  ordering = TRUE,
  dom = 'lf<"top"i>tpr<"bottom"B>RSPQ',
  buttons = c('copy', 'csv', 'excel')
))


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
  
  # remove duplicated calcuation
  predp <- predp %>% dplyr::distinct()
  
  # if evalue differ between the two calculation, use the smaller evalue
  predp2 <- predp %>% 
    dplyr::group_by(id1, id2) %>% 
    dplyr::summarise(evalue = min(evalue), .groups = "drop") %>% 
    dplyr::left_join(predp, by = c("id1", "id2", "evalue"))
  
  
  # add information of the peptides
  if(other_info){
    predp2 <- predp2 %>% 
      arrange(id1, id2) %>% 
      as_tibble() %>% 
      left_join({pep_dt %>% rename_with(~paste0("id1_", .))}, 
                by = c("id1" = "id1_u_pep_id")) %>% 
      left_join({pep_dt %>% rename_with(~paste0("id2_", .))}, 
                by = c("id2" = "id2_u_pep_id"))
  }
  
  return(predp2)
}



# d) function to plot the network -----------------------------------------

epitope_network_visualization <- function(net_df, 
                                          color_var = "taxon_genus",
                                          color_title = "",
                                          color_pal = c("#0072B2", "#009E73", "#F0E442", "#D55E00"),
                                          edge_var = "sim_score",
                                          vertex_size_var = "frequency",
                                          fig_title = "peptide network",
                                          interactive_plot = TRUE,
                                          var_to_show = ""){
  
  ### Colors of the vertices
  ncolor <- length(unique(net_df[, color_var]))
  vertex_color_pal <- colorRampPalette(color_pal)(ncolor)  # color palette
  
  
  ### Prepare the text shown when mouse hover (plotly)
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
                                no = (y + yend) / 2))
    
    # text to show when mouse hover
    net_df$edge_text <- ""
    for(v in var_to_show){
      net_df <- net_df %>% 
        mutate(edge_text = paste0(edge_text, "<br>", v, ": ", get(v)))
    }
    
    net_df <- net_df %>% 
      dplyr::mutate(text = ifelse(test = is.na(string_compare),
                                  yes  = paste0(paste0("u_pep_id: ", name), 
                                                edge_text), 
                                  no   = paste(paste0("peptide 1: ", name),
                                               paste0("peptide 2: ", nameend),
                                               paste0("E-value: ", evalue),
                                               paste0("alignment: ", string_compare), 
                                               paste0("sequence similarity: ", round(sim_score, 3)),
                                               paste0("Ab reactivity correlation: ", round(cor, 3)),
                                               paste0("Ab hit jaccard index: ", round(jaccard, 3)),
                                               sep = "<br>")))
  }
  
  
  ### Colors of the edges
  # check if any edge should be plotted
  any_edge <- any(!(is.na(net_df[, edge_var])))
  if(all(net_df[, edge_var] == 0 | is.na(net_df[, edge_var]))) any_edge <- FALSE
  
  # transform evalue to -log scale as it decrease exponsionally as the score of match increases
  net_df <- net_df %>% dplyr::mutate(evalue_log = -log(evalue))
  
  # hide edges with jaccard = 0
  net_df <- net_df %>% dplyr::mutate(jaccard = ifelse(test = jaccard == 0,
                                                      yes = NA, 
                                                      no = jaccard))
  
  # remove edges with missing values
  net_df <- net_df %>% 
    filter(if_any(.cols = all_of(edge_var), ~ ( !is.na(.) | is.na(string_compare))))
  
  # prepare the net_df for multiple edges
  if(length(edge_var) > 1){ 
    
    # ggplotly does not support geom_curve so plots with multiple edges cannot be converted to interactive plot
    interactive_plot <- FALSE
    
    # split correlation coefficients to positive and negative
    if("cor" %in% edge_var){
      net_df <- net_df %>% 
        dplyr::mutate(cor_pos = ifelse(test = cor >= 0,
                                       yes  = cor,
                                       no   = NA),
                      cor_neg = ifelse(test = cor < 0,
                                       yes = -cor, 
                                       no = NA))
      edge_var <- c(edge_var[edge_var != "cor"], "cor_pos", "cor_neg")
    }
    
    if("evalue" %in% edge_var){
      edge_var[edge_var == "evalue"] <- "evalue_log"
      net_df$evalue_log <- (net_df$evalue_log - (-log(100))) / ((-log(0.01)) - (-log(100)))
    }
    
    # sort the edge_var
    edge_var <- ordered(edge_var, levels = c("sim_score", "evalue_log", "cor_pos", "cor_neg", "jaccard"))
    edge_var <- as.character(sort(edge_var))
    
    
    # colors and curvatures for multiple edges
    edge_colors <- c(sim_score  = "#000000", 
                     evalue_log = "#332288", 
                     cor_pos    = "#E69F00", 
                     cor_neg    = "#56B4E9", 
                     jaccard    = "994F00")
    
    if((length(edge_var) %% 2) == 1){
      edge_curves <- c(0, -0.02, 0.02, -0.04, 0.04)
    }else{
      edge_curves <- c(-0.02, 0.02, -0.04, 0.04, -0.06, 0.06)
    }
  }
  
  
  
  ### Static ggplot
  net_fig <- ggplot(net_df, aes(x = x, y = y, xend = xend, yend = yend))
  
  if(any_edge){
    
    if(length(edge_var) == 1){
      
      if(edge_var == "sim_score"){
        net_fig <- net_fig +
          geom_edges(aes(color = sim_score), size = 1, alpha = 0.5) +
          scale_color_gradient(low = "white", high = "black",
                               limits = c(0, 1))
      }else if(edge_var == "evalue"){
        net_fig <- net_fig +
          geom_edges(aes(color = evalue_log), size = 1, alpha = 0.5) +
          scale_color_gradient(low = "gray95", high = "black",
                               limits = c(-log(100), -log(0.01)))
      }else if(edge_var == "cor"){
        net_fig <- net_fig +
          geom_edges(aes(color = cor), size = 1, alpha = 0.5, na.rm = TRUE) +
          scale_color_gradient2(low = "#56B4E9", mid = "gray95", high = "#E69F00",
                                limits = c(-1, 1))
      }else if(edge_var == "jaccard"){
        net_fig <- net_fig +
          geom_edges(aes(color = jaccard), size = 1, alpha = 0.5, na.rm = TRUE) +
          scale_color_gradient(low = "white", high = "black",
                               limits = c(0, 1))
      }
      
    }else{
      
      for(v in 1:length(edge_var)){
        temp_alpha <- net_df %>% filter(!is.na(string_compare)) %>% pull(get(edge_var[v]))
        temp_alpha[is.na(temp_alpha)] <- 0
        temp_color <- alpha(colour = edge_colors[edge_var[v]],
                            alpha = temp_alpha)
        
        net_fig <- net_fig +
          geom_edges(size = 0.5, color = temp_color, curvature = edge_curves[v], na.rm = TRUE)
      }
      
    }
    
  }
  
  # add nodes: color represents the variable of selection and size represents the frequency in the population
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
    scale_size_continuous(limits = c(0, 1), 
                          range = c(2, 5))
  
  # add ggplot theme elements
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
  
  
  ### Plotly interactive figure
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
                  dragRange = TRUE)
      
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
      
      uiOutput("cluster_variable"),
      
      uiOutput("color_variable"),
    
    )
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
                tabPanel("Table: peptide metadata",
                         dataTableOutput("peptide_info_table")),
                
                tabPanel("Table: antibody reactivity",
                         dataTableOutput("ab_reactivity_table")),
                
                tabPanel("Table: peptide pairwise sequence similarity",
                         dataTableOutput("pairwise_seq_sim")),
                
                tabPanel("Table: peptide pairwise sequence BLASTP",
                         dataTableOutput("pairwise_blastp")),
                
                tabPanel("Table: peptide pairwise antibody reactivity correlation",
                         dataTableOutput("pairwise_cor")),
                
                tabPanel("Table: peptide pairwise jaccard index",
                         dataTableOutput("pairwise_jaccard")),
                
                tabPanel("Table: peptide pairwise calculation",
                         dataTableOutput("pairwise_all")),
                
                tabPanel("Table: network data",
                         dataTableOutput("network_dt")),
                
                tabPanel("Network: sequence similarity", 
                         plotlyOutput("plotly_seq", height = "1000px")),
                
                tabPanel("Network: sequence BLASTP", 
                         plotlyOutput("plotly_blastp", height = "1000px")),
                
                tabPanel("Network: antibody reactivity correlation",
                         plotlyOutput("plotly_cor", height = "1000px")),
                
                tabPanel("Network: antibody hit jaccard index",
                         plotlyOutput("plotly_jaccard", height = "1000px")),
                
                tabPanel("Network: all measures combined",
                         plotOutput("plot_combined", height = "1000px"))
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
    
    # create a multiple selection UI to allow the user to decide which variables to filter
    selectizeInput(inputId = "filter_var", 
                   label = "Variables for further filtering",
                   choices = filter_list,
                   multiple = TRUE, 
                   options = list(placeholder = "type variable namnes"))
    
  })
  
  output$cluster_variable <- renderUI({
    
    # a list of variables in peptide info
    filter_list <- names(peptide_metadata())
    filter_list <- filter_list[!(filter_list %in% c("u_pep_id", "pep_aa"))]
    
    # create a multiple selection UI to allow the user to decide which variables to filter
    selectizeInput(inputId = "cluster_var", 
                   label = "Variables for clustering peptides",
                   choices = filter_list,
                   multiple = TRUE, 
                   options = list(placeholder = "type variable namnes"))
    
  })
  
  output$color_variable <- renderUI({
    
    # a list of variables in peptide info
    filter_list <- names(peptide_metadata())
    filter_list <- filter_list[!(filter_list %in% c("u_pep_id", "pep_aa"))]
    
    # create a multiple selection UI to allow the user to decide which variables to filter
    selectInput(inputId = "color_var", 
                label = "Variables for peptide colors",
                choices = filter_list,
                multiple = FALSE,
                selectize = FALSE)
    
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
  data_filtered <- reactiveValues()
  
  observeEvent(input$load, {
    
    ## filter peptide metadata according to the selected variables -----
    var_list <- input$filter_var
    
    dt <- peptide_metadata()

    for(v in var_list){

      if(class(dt[[v]]) %in% c("integer", "numeric")){

        dt <- dt %>%
          dplyr::filter(between(get(v), input[[v]][1], input[[v]][2]) | is.na(get(v)))
        
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
      filter(between(frequency, freq_threshold[1], freq_threshold[2]) | is.na(frequency)) %>% 
      mutate(u_pep_id = ordered(u_pep_id, levels = u_pep_id))
    
    data_filtered$pep <- dt
    
    
    
    ## filter Ab reactivity profile -----
    data_filtered$ab <- ab_data() %>% 
      select(u_pep_id, all_of(sample_list)) %>% 
      filter(u_pep_id %in% data_filtered$pep$u_pep_id) %>% 
      mutate(u_pep_id = ordered(u_pep_id, levels = levels(data_filtered$pep$u_pep_id))) %>% 
      arrange(u_pep_id)

  })
  
  

  # d) print the number of peptides and samples selected ----------------------------------------

  output$n_node <- renderUI({
    
    req(data_filtered$pep)
    req(data_filtered$ab)
    
    peptides_check <- all(data_filtered$pep$u_pep_id %in% data_filtered$ab$u_pep_id)
    
    n_pep <- nrow(data_filtered$pep)
    n_pep2 <- nrow(data_filtered$ab)
    n_sam <- ncol(data_filtered$ab) - 1
    
    if(n_pep == n_pep2 & peptides_check){
      text_to_print <- paste(paste("Number of peptides:", n_pep),
                             paste("Number of samples:", n_sam),
                             sep = '<br/>')
    }else{
      text_to_print <- paste(paste("Number of peptides with metadata:", n_pep),
                             paste("Number of peptides with Ab reactivity:", n_pep2), 
                             paste("Number of samples:", n_sam),
                             "Warning: Some peptides in the metadata file do not have Ab reactivity!",
                             sep = '<br/>')
    }
    
    if(n_pep > 50){
      text_to_print <- paste(text_to_print, 
                             "Warning: more than 50 peptides in the dataset.
                             Calculation may take a long time and the network may be too dense.",
                             sep = '<br/>')
    }
    
    return(HTML(text_to_print))
    
  })
  
  # output$peptide_info_table <- renderDataTable(data_filtered$pep)
  
  output$peptide_info_table <- DT::renderDataTable(
    DT::datatable(data = data_filtered$pep,
                  extensions = "Buttons",
                  filter = "top", 
                  class = "display"
    )
  )
  
  output$ab_reactivity_table <- DT::renderDataTable(
    DT::datatable(data = data_filtered$ab,
                  extensions = "Buttons",
                  filter = "top", 
                  class = "display"
    )
  )
  
  
  
  # e) calculate the sequence similarity and reactivity correlation -----------------------------
  peptide_pairwise <- reactiveValues()
  
  observeEvent(input$calculate, {
    
    withProgress(message = "calculating", value = 0, { 
    
    peptide_id_list <- levels(data_filtered$pep$u_pep_id)
    
    setProgress(value = 0.05, message = "start computing")
    
    ### sequence similarity analysis using the "virlink" package
    peptide_seq_sim <- peptide_pairwise_alignment(
      peptides        = data_filtered$pep,
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
    peptide_blastp <- blastp_dm(pep_dt     = data_filtered$pep, 
                                fasta_dir  = input$local_dir, 
                                blastp_arg = input$blastp_arg,
                                other_info = TRUE)
    
    peptide_blastp <- peptide_blastp %>% filter(id1 != id2)
    
    peptide_pairwise$blastp <- peptide_blastp
    
    setProgress(value = 0.5, message = "finished BLASTP")
    
    
    ### antibody reactivity correlation
    ab_reactivity_formatted <- data_filtered$ab %>% 
      mutate(u_pep_id = as.character(u_pep_id)) %>% 
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
    
    
    
    ### match the id1 and id2 between sequence alignment, ab reactivity correlation and jaccard index
    
    # peptide alignment:
    peptide_seq_sim <- peptide_seq_sim %>% 
      mutate(id1 = ordered(id1, levels = peptide_id_list),
             id2 = ordered(id2, levels = peptide_id_list))
    
    peptide_seq_sim[peptide_seq_sim$id1 > peptide_seq_sim$id2, c("id1", "id2")] <- 
      peptide_seq_sim[peptide_seq_sim$id1 > peptide_seq_sim$id2, c("id2", "id1")]
    
    peptide_seq_sim <- peptide_seq_sim %>% arrange(id1, id2)
    
    # peptide alignment:
    peptide_blastp <- peptide_blastp %>% 
      mutate(id1 = ordered(id1, levels = peptide_id_list),
             id2 = ordered(id2, levels = peptide_id_list))
    
    peptide_blastp[peptide_blastp$id1 > peptide_blastp$id2, c("id1", "id2")] <- 
      peptide_blastp[peptide_blastp$id1 > peptide_blastp$id2, c("id2", "id1")]
    
    peptide_blastp <- peptide_blastp %>% arrange(id1, id2)
    
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
      dplyr::left_join(ab_cor, by = c("id1", "id2")) %>% 
      dplyr::left_join(ab_jaccard, by = c("id1", "id2"))
    
    })
    
  })
  


  # f) further filtering the pairwise calculations ----------------------------------------------
  
  pairwise_filter <- reactiveValues()
  
  observeEvent(input$filter, {
    
    req(peptide_pairwise)
    
    # filter the BLASTP results
    pairwise_filter$blastp <- peptide_pairwise$blastp %>% 
      filter(between(evalue, input$evalue[1], input$evalue[2]) | is.na(evalue))
    
    # filter the sequence alignment results
    pairwise_filter$seq <- peptide_pairwise$seq %>% 
      filter(between(nchar,            input$nchar[1],            input$nchar[2])            | is.na(nchar)) %>% 
      filter(between(matches,          input$matches[1],          input$matches[2])          | is.na(matches)) %>% 
      filter(between(match_seq_length, input$match_seq_length[1], input$match_seq_length[2]) | is.na(match_seq_length)) %>%
      filter(between(sim_score,        input$sim_score[1],        input$sim_score[2])        | is.na(sim_score))
    
    # filter the Ab correlation
    pairwise_filter$cor <- peptide_pairwise$cor %>% 
      filter(between(cor, input$cor[1], input$cor[2]) | is.na(cor))
      
    # filter the Ab jaccard index
    pairwise_filter$jac <- peptide_pairwise$jac %>% 
      filter(between(jaccard, input$jaccard[1], input$jaccard[2]) | is.na(jaccard))
    
    # combine the results
    pairwise_filter$all <- peptide_pairwise$all %>% 
      filter(between(evalue,           input$evalue[1],           input$evalue[2])           | is.na(evalue)) %>% 
      filter(between(nchar,            input$nchar[1],            input$nchar[2])            | is.na(nchar)) %>% 
      filter(between(matches,          input$matches[1],          input$matches[2])          | is.na(matches)) %>% 
      filter(between(match_seq_length, input$match_seq_length[1], input$match_seq_length[2]) | is.na(match_seq_length)) %>%
      filter(between(sim_score,        input$sim_score[1],        input$sim_score[2])        | is.na(sim_score)) %>% 
      filter(between(cor,              input$cor[1],              input$cor[2])              | is.na(cor))%>% 
      filter(between(jaccard,          input$jaccard[1],          input$jaccard[2])          | is.na(jaccard))
      
    # pairwise_filter$all <- pairwise_filter$seq %>% 
    #   left_join(pairwise_filter$blastp, 
    #             by = intersect(colnames(pairwise_filter$seq),
    #                            colnames(pairwise_filter$blastp))) %>% 
    #   left_join(pairwise_filter$cor, by = c("id1", "id2")) %>% 
    #   left_join(pairwise_filter$jac, by = c("id1", "id2"))
    
  })
  
  
  # show the pairwise calculation results after filtering
  output$pairwise_seq_sim <- DT::renderDataTable(
    DT::datatable(data = pairwise_filter$seq,
                  extensions = "Buttons",
                  filter = "top", 
                  class = "display"
    )
  )
  
  output$pairwise_blastp <-  DT::renderDataTable(
    DT::datatable(data = pairwise_filter$blastp,
                  extensions = "Buttons",
                  filter = "top", 
                  class = "display"
    )
  )

  
  output$pairwise_cor <- DT::renderDataTable(
    DT::datatable(data = pairwise_filter$cor,
                  extensions = "Buttons",
                  filter = "top", 
                  class = "display"
    )
  )
  
  output$pairwise_jaccard <- DT::renderDataTable(
    DT::datatable(data = pairwise_filter$jac,
                  extensions = "Buttons",
                  filter = "top", 
                  class = "display"
    )
  )
  
  output$pairwise_all <- DT::renderDataTable(
    DT::datatable(data = pairwise_filter$all,
                  extensions = "Buttons",
                  filter = "top", 
                  class = "display"
    )
  )
  

  

  # g) network construction ---------------------------------------------------------------------

  network_dt <- reactiveValues()
  
  observeEvent(input$plot, {
    
    # edges
    peptide_edges <- peptide_pairwise$all
    
    # examine if the two peptides have the same value in the variables selected for clustering
    if(!is.null(input$cluster_var)){
      
      for(v in input$cluster_var){
        peptide_edges <- peptide_edges %>% 
          mutate(!!paste0("same_", v) := get(paste0("id1_", v)) == get(paste0("id2_", v))) %>% 
          mutate(!!paste0("same_", v) := ifelse(test = is.na(get(paste0("same_", v))), 
                                                yes = FALSE, 
                                                no = get(paste0("same_", v))))
      }
      
    }
    
    # total number of matching clustering variables
    peptide_edges <- peptide_edges %>% 
      rowwise() %>% 
      mutate(cluster_weight = sum(c_across(starts_with("same_")))) %>% 
      ungroup()
    
    
    # construct an igraph object
    peptide_net <- igraph::graph_from_data_frame(d        = peptide_edges,
                                                 vertices = data_filtered$pep,
                                                 directed = FALSE)
    
    # add weight to edges
    weight_per_var <- 1
    E(peptide_net)$weight <- E(peptide_net)$sim_score + 
      weight_per_var * E(peptide_net)$cluster_weight

    
    # fortify the igraph to a data frame suitable for ggplot2
    set.seed(input$seed)
    suppressWarnings({
      peptide_net_df <- ggnetwork::ggnetwork(x = peptide_net,
                                             layout = igraph::with_fr())
    })
    
    
    # filter the edges
    peptide_net_df <- peptide_net_df %>% 
      filter(between(evalue,           input$evalue[1],           input$evalue[2])           | is.na(evalue)) %>% 
      filter(between(nchar,            input$nchar[1],            input$nchar[2])            | is.na(nchar)) %>% 
      filter(between(matches,          input$matches[1],          input$matches[2])          | is.na(matches)) %>% 
      filter(between(match_seq_length, input$match_seq_length[1], input$match_seq_length[2]) | is.na(match_seq_length)) %>%
      filter(between(sim_score,        input$sim_score[1],        input$sim_score[2])        | is.na(sim_score)) %>% 
      filter(between(cor,              input$cor[1],              input$cor[2])              | is.na(cor))%>% 
      filter(between(jaccard,          input$jaccard[1],          input$jaccard[2])          | is.na(jaccard))
    
    network_dt$dt <- peptide_net_df

  })
  
  
  output$network_dt  <- DT::renderDataTable(
    DT::datatable(data = network_dt$dt,
                  extensions = "Buttons",
                  filter = "top", 
                  class = "display"
    )
  )
  


  # h) visualization ----------------------------------------------------------------------------

  # visualize the network: sequence similarity
  output$plotly_seq <- renderPlotly({
    
    req(network_dt$dt)
    
    var_list <- names(peptide_metadata())
    var_list <- var_list[!(var_list %in% c("u_pep_id"))]
    
    if(nrow(network_dt$dt) > 0){
      
      epitope_network_visualization(net_df = network_dt$dt,
                                    color_var = input$color_var,
                                    color_title = "",
                                    color_pal = manualPalette,
                                    edge_var = "sim_score",
                                    vertex_size_var = "frequency",
                                    fig_title = "network of sequence similarity",
                                    interactive_plot = TRUE, 
                                    var_to_show = var_list)
      
    }else{
      plotly_empty(type = "scatter", mode = "markers")
    }
    
  })
  
  
  # visualize the network: BLASTP
  output$plotly_blastp <- renderPlotly({
    
    req(network_dt$dt)
    
    var_list <- names(peptide_metadata())
    var_list <- var_list[!(var_list %in% c("u_pep_id"))]
    
    if(nrow(network_dt$dt) > 0){
      
      epitope_network_visualization(net_df = network_dt$dt,
                                    color_var = input$color_var,
                                    color_title = "",
                                    color_pal = manualPalette,
                                    edge_var = "evalue",
                                    vertex_size_var = "frequency",
                                    fig_title = "network of sequence BLASTP E-value",
                                    interactive_plot = TRUE, 
                                    var_to_show = var_list)
      
    }else{
      plotly_empty(type = "scatter", mode = "markers")
    }
    
  })
  
  
  # visualize the network: co-occurrence jaccard index
  output$plotly_jaccard <- renderPlotly({
    
    req(network_dt$dt)
    
    var_list <- names(peptide_metadata())
    var_list <- var_list[!(var_list %in% c("u_pep_id"))]
    
    if(nrow(network_dt$dt) > 0){
      
      epitope_network_visualization(net_df = network_dt$dt,
                                    color_var = input$color_var,
                                    color_title = "",
                                    color_pal = manualPalette,
                                    edge_var = "jaccard",
                                    vertex_size_var = "frequency",
                                    fig_title = "network of antibody hit jaccard index",
                                    interactive_plot = TRUE, 
                                    var_to_show = var_list)
    }else{
      plotly_empty(type = "scatter", mode = "markers")
    }
  })
  
  
  # visualize the network: correlation
  output$plotly_cor <- renderPlotly({
    
    req(network_dt$dt)
    
    var_list <- names(peptide_metadata())
    var_list <- var_list[!(var_list %in% c("u_pep_id"))]
    
    if(nrow(network_dt$dt) > 0){
      
      epitope_network_visualization(net_df = network_dt$dt,
                                    color_var = input$color_var,
                                    color_title = "",
                                    color_pal = manualPalette,
                                    edge_var = "cor",
                                    vertex_size_var = "frequency",
                                    fig_title = "network of antibody reactivity correlation",
                                    interactive_plot = TRUE, 
                                    var_to_show = var_list)
    }else{
      plotly_empty(type = "scatter", mode = "markers")
    }
  })
  
  
  # visualize the network: all measures combined (static plot)
  output$plot_combined <- renderPlot({
    
    req(network_dt$dt)
    
    var_list <- names(peptide_metadata())
    var_list <- var_list[!(var_list %in% c("u_pep_id"))]
    
    if(nrow(network_dt$dt) > 0){
      
      epitope_network_visualization(net_df = network_dt$dt,
                                    color_var = input$color_var,
                                    color_title = "",
                                    color_pal = manualPalette,
                                    edge_var = c("sim_score", "evalue", "cor", "jaccard"),
                                    vertex_size_var = "frequency",
                                    fig_title = "sequence similarity: gray, E-value: purple, correlation coefficient: orange-blue, jaccard index: brown",
                                    interactive_plot = FALSE, 
                                    var_to_show = var_list)
      
    }else{
      plot.new()
    }
    
  })
  
}




# 4. Run the application --------------------------------------------------
shinyApp(ui = ui, server = server)
