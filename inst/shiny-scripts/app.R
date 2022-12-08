library(shiny)

k <- 1

# Define UI for random distribution app ----
ui <- fluidPage(

  # App title ----
  titlePanel("alfa - Alignment-free analysis tool"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      fileInput("file", "Choose FASTA File",
                multiple = FALSE,
                accept = c(".txt",
                           ".fasta",
                           ".")),

      # Input: Select the random distribution type ----
      radioButtons("dist", "Distance function:",
                   c("Euclidean" = "euclidean",
                     "Standardized Euclidean" = "stdEuclidean")),

      # br() element to introduce extra vertical spacing ----
      br(),

      # Input: Slider for the number of observations to generate ----
      numericInput("k",
                  "k-mer size:",
                  value = 5,
                  min = 1,
                  max = 100),
      actionButton("submitBtn", "Submit job")

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", plotOutput("plot")),
                  tabPanel("Phylogenetic Tree", verbatimTextOutput("phylo")),
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

  req(input$file)
  req(input$dist)
  req(input$k)

  observeEvent(input$submitBtn, {
    # submit job

  })



  d <- reactive({
    dist <- switch(input$dist,
                   norm = rnorm,
                   unif = runif,
                   lnorm = rlnorm,
                   exp = rexp,
                   rnorm)

    dist(input$n)
  })

  # Generate a plot of the data ----
  # Also uses the inputs to build the plot label. Note that the
  # dependencies on the inputs and the data reactive expression are
  # both tracked, and all expressions are called in the sequence
  # implied by the dependency graph.
  output$plot <- renderPlot({
    dist <- input$dist
    n <- input$n

    hist(d(),
         main = paste("r", dist, "(", n, ")", sep = ""),
         col = "#75AADB", border = "white")
  })

  # Generate a summary of the data ----
  output$summary <- renderPrint({
    summary(d())
  })

  # Generate an HTML table view of the data ----
  output$table <- renderTable({
    d()
  })

}

# Create Shiny app ----
shinyApp(ui, server)
