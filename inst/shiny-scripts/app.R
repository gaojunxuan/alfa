library(shiny)

maxK <- 10

# Define UI for random distribution app ----
ui <- fluidPage(

  # App title ----
  titlePanel("alfa - Alignment-free analysis tool"),
  tags$p("alfa is an R package that allows for alignment-free comparative
         genomic analysis. It uses standard distance functions including
         Euclidean distance and standardized Euclidean distance to compute
         the distance between sequences instead of aligning them using
         dynamic programming. You can upload you sequence file in FASTA format
         and analyze it using alfa. A sample data can be found at
         https://raw.githubusercontent.com/gaojunxuan/alfa/main/inst/extdata/usflu.fasta

         You have to use standardized Euclidean distance for the example dataset
         because the pairwise distance would otherwise be too big to construct
         a phylogenetic tree using neighbor joining."),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
      tags$p("The file should be a FASTA file containing more than one
             possibly biologically related sequences."),
      fileInput("file", "Choose FASTA File",
                multiple = FALSE,
                accept = c(".txt",
                           ".fasta",
                           ".")),

      # Input: Select the random distribution type ----
      radioButtons("dist", "Distance function:",
                   c("Euclidean" = "euclidean",
                     "Standardized Euclidean" = "standardizedEuclidean")),

      radioButtons("phyloAlgo", "Phylogenetic tree construction method:",
                   c("UPGMA" = "upgma",
                     "Neighbor-joining" = "nj")),

      # br() element to introduce extra vertical spacing ----
      br(),

      # Input: Slider for the number of observations to generate ----
      numericInput("k",
                  "k-mer size:",
                  value = 5,
                  min = 1,
                  max = maxK),
      actionButton("submitBtn", "Submit job")

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", plotOutput("plot")),
                  tabPanel("Phylogenetic Tree", plotOutput("phylo")),
                  tabPanel("Pairwise Distance Matrix", tableOutput("mat"))
      )

    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {

  # Reactive expression to generate the requested distribution ----
  # This is called whenever the inputs change. The output functions
  # defined below then use the value computed from this expression

  submitJob <- observeEvent(input$submitBtn, {
    # submit job
    req(input$file)
    file <- input$file
    dnaSet <- alfa::readFASTA(file$datapath)
    seqLen <- max(stringr::str_length(dnaSet))
    maxK <- seqLen
    dist <- switch(input$dist,
                   euclidean = "euclidean",
                   standardizedEuclidean = "standardizedEuclidean",
                   "euclidean")

    distMat <- alfa::createDistanceMatrix(dnaSet, dist, k)
    output$plot <- renderPlot({
      if (exists("distMat") && !is.null(distMat)) {
        heatmap(distMat)
      }
    })

    output$mat <- renderTable({
      if (exists("distMat") && !is.null(distMat)) {
        distMat
      }
    })

    output$phylo <- renderPlot({
      if (exists("distMat") && !is.null(distMat)) {
        if (input$phyloAlgo == "nj") {
          alfa::neighborJoiningTree(distMat)
        } else {
          alfa::upgmaTree(distMat)
        }
      }
    })
  })

  # update max k value based on max sequence length of input
  updateMaxK <- observeEvent(input$file, {
    file <- input$file
    dnaSet <- alfa::readFASTA(file$datapath)
    seqLen <- max(stringr::str_length(dnaSet))
    maxK <- seqLen
    updateNumericInput(getDefaultReactiveDomain(), "k",
                       "k-mer size:",
                       value = 5,
                       min = 1,
                       max = maxK)
  })

}

# Create Shiny app ----
shinyApp(ui, server)
