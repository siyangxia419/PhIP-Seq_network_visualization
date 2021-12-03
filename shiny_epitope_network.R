###
# Interactive network visualization of viral peptides in PhIP-Seq libraries
# show peptides' pairwise sequence similarity, correlation, and cooccurrence
#
# This is a R shiny app. To launch the app, click "Run App". 
#
# Contributor:
# Siyang Xia
# Daniel Monaco
# H. Benjamin Larman
#
# Version: 2021-12-01
###



# 1. Preparation ----------------------------------------------------------


# a) packages -------------------------------------------------------------

# install the packages below if not already installed
load_lib <- c("shiny", "shinythemes", 
              "tidyverse", "here", "devtools", "tools", 
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


# c) function to plot the network -----------------------------------------

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
  
  # sidebar with a slider input for number of bins 
  fluidRow(
    
    # input
    column(2, 
           
      # allow the user to upload a file that contains peptide information
      fileInput(inputId = "peptide_upload",
                label = "Choose peptide information file",
                multiple = FALSE,
                accept = c(".csv", ".tsv")),
      
      # allow the user to upload a file that contains PhIP-Seq antibody reactivity profile
      fileInput(inputId = "reactivity_upload",
                label = "Choose antibody reactivity file",
                multiple = FALSE,
                accept = c(".csv", ".tsv")),
      
      # allow the user to select variables to perform further filtering of the peptides
      uiOutput("filter_variable"),
      
      uiOutput("filter1"), 
      
      # select virus families
      selectizeInput(inputId = "family", 
                     label = "Virus family",
                     choices = NULL,
                     multiple = TRUE, 
                     options = list(placeholder = "type names of virus family")),
      
      # select virus genus
      selectizeInput(inputId = "genus", 
                     label = "Virus genus",
                     choices = NULL,
                     multiple = TRUE, 
                     options = list(placeholder = "type names of virus genus")),
      
      # select virus species
      selectizeInput(inputId = "species", 
                     label = "Virus species",
                     choices = NULL,
                     multiple = TRUE, 
                     options = list(placeholder = "type names of virus species")),
      
      # select virus organism
      selectizeInput(inputId = "organism", 
                     label = "Virus organism",
                     choices = NULL,
                     multiple = TRUE, 
                     options = list(placeholder = "type names of virus organism")),
      
      # select protein uniprot
      selectizeInput(inputId = "uniprot", 
                     label = "UniProt accession number",
                     choices = NULL,
                     multiple = TRUE, 
                     options = list(placeholder = "type UniProt accession numbers")),
      
      # filter of peptide frequency
      sliderInput(inputId = "frequency",
                  label = "Enrichment frequency",
                  min = 0,
                  max = 1,
                  value = c(0, 1), 
                  step = 0.01, 
                  ticks = FALSE, 
                  dragRange = TRUE),
      
      # seed for determine the coordinate of vertices
      numericInput(inputId = "seed", 
                   label = "Seed", 
                   min = 1, 
                   max = 1000, 
                   value = 111, 
                   step = 1),
      
      br(),
      
      # action button to generate the network figure
      actionButton(inputId = "go", label = "Submit")
      
    ),
    
    # filtering the network before visualization
    column(2, 
           
      sliderInput(inputId = "nchar",
                  label = "Length of alignment",
                  min = 1,
                  max = 56,
                  value = c(7, 56), 
                  step = 1, 
                  ticks = FALSE, 
                  dragRange = TRUE),
      
      sliderInput(inputId = "matches",
                  label = "# of matched amino acids",
                  min = 1,
                  max = 56,
                  value = c(5, 56),
                  step = 1, 
                  ticks = FALSE, 
                  dragRange = TRUE),
      
      sliderInput(inputId = "match_seq_length",
                  label = "Continuous matches",
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
      
      sliderInput(inputId = "cor",
                  label = "Ab reactivity correlation",
                  min = -1,
                  max = 1,
                  value = c(-1, 1), 
                  step = 0.01, 
                  ticks = FALSE, 
                  dragRange = TRUE),
      
      sliderInput(inputId = "jaccard",
                  label = "Ab hit jaccard index",
                  min = 0,
                  max = 1,
                  value = c(0, 1), 
                  step = 0.01, 
                  ticks = FALSE, 
                  dragRange = TRUE),
      
      selectInput(inputId = "same_family", 
                  label = "Same family", 
                  choices = c("yes", "no", "both"), 
                  selected = "both", 
                  multiple = FALSE),
      
      selectInput(inputId = "same_genus", 
                  label = "Same genus", 
                  choices = c("yes", "no", "both"), 
                  selected = "both", 
                  multiple = FALSE),
      
      selectInput(inputId = "same_species", 
                  label = "Same species", 
                  choices = c("yes", "no", "both"), 
                  selected = "both", 
                  multiple = FALSE),
      
      selectInput(inputId = "same_organism", 
                  label = "Same organism", 
                  choices = c("yes", "no", "both"), 
                  selected = "both", 
                  multiple = FALSE),
      
      selectInput(inputId = "tiling", 
                  label = "Tiling epitopes", 
                  choices = c("yes", "no", "both"), 
                  selected = "both", 
                  multiple = FALSE),
      
      br(),
      
      actionButton(inputId = "filter", label = "Filter")
    ),
    
    # Output: interactive network by plotly
    column(8,
           fluidRow(

             column(10, 
                    textOutput(outputId = "n_node")), 
             
             column(2, 
                    selectInput(inputId = "color_var", 
                                label = "Color of nodes", 
                                choices = c("family", "genus", "species", "organism"), 
                                selected = "family",
                                multiple = FALSE)),
             
             br(),
             
             tabsetPanel(type = "tabs",
                         tabPanel("sequence similarity", 
                                  plotlyOutput("plotly_seq")),
                         tabPanel("antibody reactivity correlation",
                                  plotlyOutput("plotly_cor")),
                         tabPanel("antibody hit jaccard index",
                                  plotlyOutput("plotly_jaccard"))
             )
             
           )
    )
  )
)



# 3. Server ---------------------------------------------------------------


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # update peptide information data frame based on whether the user upload a file
  peptide_info_dt <- reactive({
  
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
    filter_list <- names(peptide_info_dt())
    filter_list <- filter_list[!(filter_list %in% c("u_pep_id", "pep_aa"))]
    
    default_filter <- grep(pattern = "taxon_", x = filter_list, value = TRUE)
    
    # create a multiple selection UI to allow the user to decide which variables to filter
    selectInput(inputId = "filter_var", 
                label = "Vriables for further filtering",
                choices = filter_list,
                selectize = TRUE,
                multiple = TRUE, 
                selected = default_filter)
    
  })
  
  
  # UI for filtering peptide
  output$filter1 <- renderUI({
    
    dynamic_ui <- lapply(
      X = input$filter_var, 
      FUN = function(x) {
        if(class(peptide_info_dt()[[x]]) %in% c("integer", "numeric")){
          
          x_fullrange <- range(peptide_info_dt()[[x]], na.rm = TRUE)
          
          # create a slider bar
          return(sliderInput(inputId = x,
                             label = x,
                             min = x_fullrange[1],
                             max = x_fullrange[2],
                             value = x_fullrange, 
                             ticks = FALSE, 
                             dragRange = TRUE))
          
        }else{
          
          x_fulllist <- sort(unique(as.character(peptide_info_dt()[[x]])))
          
          # create a multiple selection UI
          return(selectizeInput(inputId = x, 
                                label = x,
                                choices = x_fulllist,
                                multiple = TRUE, 
                                options = list(placeholder = "type variable options")))
          
        }
        
      })
    
    return(dynamic_ui)
    
  })

  
  
  
  # select virus family
  updateSelectizeInput(session, 
                       inputId = "family", 
                       choices = c("all", unique(taxa_protein$family)), 
                       selected = NULL,
                       server = TRUE)
  
  # select virus genus
  observeEvent(input$family,
               {
                 # update the genus choices according to virus family input
                 temp <- taxa_protein %>% 
                 {if(!("all" %in% input$family)){
                   dplyr::filter(., family %in% input$family)
                 }else{
                   .
                 }}
                 
                 updateSelectizeInput(session, 
                                      inputId = "genus", 
                                      choices = c("all", unique(temp$genus)), 
                                      selected = NULL,
                                      server = TRUE)
               })
  
  # select virus species
  observeEvent(input$genus,
               {
                 # update the species choices according to virus genus input
                 temp <- taxa_protein %>% 
                 {if(!("all" %in% input$family)){
                   dplyr::filter(., family %in% input$family)
                 }else{
                   .
                 }} %>% 
                 {if(!("all" %in% input$genus)){
                   dplyr::filter(., genus %in% input$genus)
                 }else{
                   .
                 }}
                 
                 updateSelectizeInput(session, 
                                      inputId = "species", 
                                      choices = c("all", unique(temp$species)), 
                                      selected = NULL,
                                      server = TRUE)
               })
  
  # select virus organism
  observeEvent(input$species,
               {
                 # update the organism choices according to virus species input
                 temp <- taxa_protein %>% 
                 {if(!("all" %in% input$family)){
                   dplyr::filter(., family %in% input$family)
                 }else{
                   .
                 }} %>% 
                 {if(!("all" %in% input$genus)){
                   dplyr::filter(., genus %in% input$genus)
                 }else{
                   .
                 }} %>% 
                 {if(!("all" %in% input$species)){
                   dplyr::filter(., species %in% input$species)
                 }else{
                   .
                 }}
                 
                 updateSelectizeInput(session, 
                                      inputId = "organism", 
                                      choices = c("all", unique(temp$organism)), 
                                      selected = NULL,
                                      server = TRUE)
               })
  
  # select protein
  observeEvent(input$organism,
               {
                 # update the protein choices according to virus organism input
                 temp <- taxa_protein %>% 
                 {if(!("all" %in% input$family)){
                   dplyr::filter(., family %in% input$family)
                 }else{
                   .
                 }} %>% 
                 {if(!("all" %in% input$genus)){
                   dplyr::filter(., genus %in% input$genus)
                 }else{
                   .
                 }} %>% 
                 {if(!("all" %in% input$species)){
                   dplyr::filter(., species %in% input$species)
                 }else{
                   .
                 }} %>% 
                 {if(!("all" %in% input$organism)){
                   dplyr::filter(., organism %in% input$organism)
                 }else{
                   .
                 }}
                 
                 updateSelectizeInput(session, 
                                      inputId = "uniprot", 
                                      choices = c("all", unique(temp$UniProt_acc)), 
                                      selected = NULL,
                                      server = TRUE)
               })
  
  

  
  # build the network
  network_dt <- eventReactive(input$go, {
    
    ### vextices
    if("all" %in% input$family){  # filter family
      vertex_d <- epitope_info
    }else{
      vertex_d <- epitope_info %>% 
        dplyr::filter(family %in% input$family)
    }
    
    if("all" %in% input$genus){  # filter genus
      vertex_d <- vertex_d
    }else{
      vertex_d <- vertex_d %>% 
        dplyr::filter(genus %in% input$genus)
    }
    
    if("all" %in% input$species){  # filter species
      vertex_d <- vertex_d
    }else{
      vertex_d <- vertex_d %>% 
        dplyr::filter(species %in% input$species)
    }
    
    if("all" %in% input$organism){  # filter virus organism
      vertex_d <- vertex_d
    }else{
      vertex_d <- vertex_d %>% 
        dplyr::filter(organism %in% input$organism)
    }
    
    if("all" %in% input$uniprot){  # filter protein UniProt accession number
      vertex_d <- vertex_d
    }else{
      vertex_d <- vertex_d %>% 
        dplyr::filter(UniProt_acc %in% input$uniprot)
    }
    
    if("freq" %in% names(epitope_info)){  # filter by peptide's enrichment frequency (column "freq")
      vertex_d <- vertex_d %>% 
        dplyr::filter(freq >= input$frequency[1], freq <= input$frequency[2])
    }else{
      vertex_d <- vertex_d %>% 
        dplyr::mutate(freq = 1)
    }
    
    
    ### edges
    edge_d <- epitope_pair %>% 
      dplyr::filter(subject_id %in% vertex_d$id,
                    pattern_id %in% vertex_d$id)
    
    # construct an igraph object
    net <- igraph::graph_from_data_frame(d = edge_d, 
                                         vertices = vertex_d, 
                                         directed = FALSE)
    
    # weight to group species and genus together
    weight_same_family <- 1
    weight_same_genus <- 1
    weight_same_species <- 1
    
    # weight of edges
    E(net)$weight <- E(net)$sim_score + 
      E(net)$same_family  * weight_same_family + 
      E(net)$same_genus   * weight_same_genus + 
      E(net)$same_species * weight_same_species
    
    # fortify the igraph to a data frame suitable for ggplot2
    set.seed(input$seed)
    suppressWarnings({
      net_fig_df <- ggnetwork::ggnetwork(x = net, 
                                         layout = igraph::with_fr())
    })
    
    return(net_fig_df)
  })
  
  
  # number of epitopes selected
  output$n_node <- renderText({
    paste("Number of peptides:", sum(is.na(network_dt()$string_compare)))
  })
  
  
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
