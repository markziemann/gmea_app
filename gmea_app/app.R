library(tidyverse)
library(ggplot2)
library(rmarkdown)
library(broom)
library(mitch)
options(shiny.maxRequestSize = 100 * 1024^2)

# Define UI
ui <- fluidPage(
  titlePanel("LAM - Infinium Enrichment Analysis Tool"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose limma CSV File",
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      fileInput("file2", "Choose GMT File",
                accept = c("text/tsv",
                           "text/tav-separated-values,text/plain",
                           ".gmt")),
      uiOutput("prioritisation"),
      uiOutput("minsetsize"),
      uiOutput("chip_selector"),
      downloadButton("download_report", "Download Report")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data", tableOutput("contents")),
        tabPanel("Summary", verbatimTextOutput("summary")),
        tabPanel("Structure", verbatimTextOutput("structure")),
        tabPanel("Structure GMT", verbatimTextOutput("structure2")),
        tabPanel("Enrichment Result", tableOutput("enrichment_data")),
        tabPanel("Enrichment Plot", plotOutput("enrichmentplot")),
        tabPanel("Regression Summary", verbatimTextOutput("regression_summary"))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  data <- reactive({
    req(input$file)
    df <- read.csv(input$file$datapath, header = TRUE, check.names = FALSE)
    return(df)
  })
  
  data2 <- reactive({
    req(input$file2)
    sets <- gmt_import(input$file2$datapath)
    return(sets)
  })
  
  output$contents <- renderTable({
    req(data())
    head(data(), rownames = TRUE)
  })
  
  output$summary <- renderPrint({
    summary(data())
  })
  
  output$structure <- renderPrint({
    str(data())
  })
  
  output$structure2 <- renderPrint({
    str(data2())
  })
  
  output$prioritisation <- renderUI({
    selectInput("prioritisation", "Prioritisation", choices = c("Effect","Significance"), selected = "Effect")
  })
  
  output$minsetsize <- renderUI({
    selectInput("minsetsize", "Minimum gene set size", choices = c("5","10"), selected = "5")
  })
  
  output$chip_selector <- renderUI({
    selectInput("chip", "Choose array chip", choices = c("450K","EPIC"), selected = "450K")
  })
  
  enrichment_data <- reactive({
    req(data(),data2(),input$chip)
    if (input$chip == "450K") {
      gt <- readRDS("450K.rds")
    }
    if (input$chip == "EPIC") {
      gt <- readRDS("EPIC.rds")
    }
    m <- mitch_import(x=data(),DEtype="limma",geneTable=gt,geneIDcol="Name")
    minsetsize <- as.numeric(input$minsetsize)
    mres <- mitch_calc(x=m,genesets=data2(),minsetsize=minsetsize, 
                       priority=tolower(input$prioritisation),cores=4)
    return(mres$enrichment_result)
  })
  
  output$enrichment_data <- renderTable({
    head(enrichment_data(),50)
  })

  plot_data <- reactive({
    edata <- head(enrichment_data(),20)
    edata <- edata[order(edata$s.dist),]
    vec <- edata$s.dist
    names(vec) <- edata$set
    names(vec) <- substr(names(vec), start = 1, stop = 60)
    names(vec) <- gsub("_"," ",names(vec))
    par(mar=c(5,40,3,3))
    barplot(abs(vec),col=sign(-vec)+3,horiz=TRUE,las=1,cex.names=1.1,xlab="Enrichment Score")
  })
  
  output$enrichmentplot <- renderPlot({
    plot_data()
  }, height = 800, width=800)
  
  regression_model <- reactive({
    req(data(), input$x_axis, input$y_axis)
    lm(formula = as.formula(paste(input$y_axis, "~", input$x_axis)), data = data())
  })
  
  regression_stats <- reactive({
    req(regression_model())
    model <- regression_model()
    summary_model <- summary(model)
    
    # Calculate correlation and its p-value
    cor_test <- cor.test(data()[[input$x_axis]], data()[[input$y_axis]])
    
    list(
      r_value = sqrt(summary_model$r.squared) * sign(summary_model$coefficients[2,1]),
      r_squared = summary_model$r.squared,
      adj_r_squared = summary_model$adj.r.squared,
      f_statistic = summary_model$fstatistic[1],
      p_value = summary_model$coefficients[2,4],
      correlation = cor_test$estimate,
      cor_p_value = cor_test$p.value
    )
  })
  
  output$regression_summary <- renderPrint({
    req(regression_model(), regression_stats())
    
    cat("Regression Summary:\n\n")
    cat("R value:", round(regression_stats()$r_value, 4), "\n")
    cat("R-squared:", round(regression_stats()$r_squared, 4), "\n")
    cat("Adjusted R-squared:", round(regression_stats()$adj_r_squared, 4), "\n")
    cat("F-statistic:", round(regression_stats()$f_statistic, 4), "\n")
    cat("P-value:", format.pval(regression_stats()$p_value, digits = 4), "\n\n")
    
    cat("Model Coefficients:\n")
    print(summary(regression_model())$coefficients)
    
    cat("\nAdditional Statistics:\n")
    cat("Correlation coefficient:", round(regression_stats()$correlation, 4), "\n")
    cat("Correlation p-value:", format.pval(regression_stats()$cor_p_value, digits = 4), "\n")
  })
  
  output$download_report <- downloadHandler(
    filename = function() {
      paste("report-", Sys.Date(), ".html", sep="")
    },
    content = function(file) {
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      params <- list(plot = plot_data(),
                     stats = regression_stats(),
                     x_var = input$x_axis,
                     y_var = input$y_axis,
                     coefficients = summary(regression_model())$coefficients)
      
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)