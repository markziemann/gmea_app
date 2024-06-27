library(markdown)
library(rmarkdown)
library(shiny)
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
        tabPanel("Instructions", includeMarkdown("intro.md")),
        tabPanel("limma CSV check", tableOutput("csvcheck1"),
                 textOutput("csvcheck2")),
        tabPanel("GMT check", "Peek at the gene sets",
                 verbatimTextOutput("gmtcheck1"),
                 textOutput("gmtcheck2"),
                 "Summary of gene set sizes",
                 verbatimTextOutput("gmtcheck3")),
        tabPanel("Enrichment Result", tableOutput("enrichment_data")),
        tabPanel("Enrichment Plot", 
                 "Pathways with FDR<0.05 with highest enrichment score",
                 plotOutput("enrichmentplot1", width = "600px", height = "500px"),
                 "Pathways with smallest p-values",
                 plotOutput("enrichmentplot2", width = "600px", height = "500px"))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  data <- reactive({
    req(input$file)
    df <- read.csv(input$file$datapath, header = TRUE, check.names = FALSE)
    colnames(df)[1] <- "Name"
    return(df)
  })
  
  data2 <- reactive({
    req(input$file2)
    sets <- gmt_import(input$file2$datapath)
    names(sets) <- gsub("_"," ",names(sets))
    return(sets)
  })
  
  output$csvcheck1 <- renderTable({
    req(data())
    head(data())
  })
  
  output$csvcheck2 <- renderText({
    req(data())
    paste("Limma CSV file has",nrow(data()),"rows of data")
  })
  
  output$gmtcheck1 <- renderPrint({
    req(data2())
    lapply(head(data2(),3),function(x) {x[1:5]} )
  })

  output$gmtcheck2 <- renderText({
    req(data2())
    paste("Number of gene sets:",length(data2()))
  })
  
  output$gmtcheck3 <- renderPrint({
    req(data2())
    summary(unlist(lapply(data2(),length)))
  })
  output$prioritisation <- renderUI({
    selectInput("prioritisation", "Prioritisation", 
                choices = c("Effect","Significance"), selected = "Effect")
  })
  
  output$minsetsize <- renderUI({
    selectInput("minsetsize", "Minimum gene set size", choices = c("5","10"), 
                selected = "5")
  })
  
  output$chip_selector <- renderUI({
    selectInput("chip", "Choose array chip", choices = c("450K","EPIC"), 
                selected = "EPIC")
  })
  
  edata <- reactive({
    req(data(),data2(),input$chip)
    if (input$chip == "450K") {
      gt <- readRDS("450K.rds")
    }
    if (input$chip == "EPIC") {
      gt <- readRDS("EPIC.rds")
    }
    if (input$chip == "EPIC") {
      gt <- readRDS("EPIC.rds")
    }
    m <- mitch_import(x=data(),DEtype="limma",geneTable=gt,geneIDcol="Name")
    minsetsize <- as.numeric(input$minsetsize)
    mres <- mitch_calc(x=m,genesets=data2(),minsetsize=minsetsize, 
                       priority=tolower(input$prioritisation),cores=4)
    return(mres)
  })
      
  enrichment_data <- reactive({
    mres <- edata() 
    return(mres$enrichment_result)
  })
  
  output$enrichment_data <- renderTable({
    head(enrichment_data(),50)
  })

  plot_data1 <- reactive({
    edata <- head(subset(enrichment_data(),p.adjustANOVA<0.05),20)
    if (nrow(edata)>1) {
      edata <- edata[order(edata$s.dist),]
      vec <- edata$s.dist
      names(vec) <- edata$set
      names(vec) <- substr(names(vec), start = 1, stop = 60)
      names(vec) <- gsub("_"," ",names(vec))
      par(mar=c(4,33,1,1))
      barplot(abs(vec),col=sign(-vec)+3,horiz=TRUE,las=1,cex.names=1,
              xlab="Enrichment Score")
      mtext("Red=higher,Blue=lower methylation",side=3,adj=1)
    } else {
      plot(1,1,cex=0,xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '')
      text(1,1,labels="None found")
    }
  })
  
  output$enrichmentplot1 <- renderPlot({
    plot_data1()
  }, height = 500, width=600)

  plot_data2 <- reactive({
    edata <- enrichment_data()
    edata <- edata[order(edata$p.adjustANOVA),]
    edata <- head(edata,20)
    edata <- edata[order(edata$s.dist),]
    vec <- edata$s.dist
    names(vec) <- edata$set
    names(vec) <- substr(names(vec), start = 1, stop = 60)
    names(vec) <- gsub("_"," ",names(vec))
    par(mar=c(4,33,1,1))
    barplot(abs(vec),col=sign(-vec)+3,horiz=TRUE,las=1,cex.names=1,
            xlab="Enrichment Score")
    mtext("Red=higher,Blue=lower methylation",side=3,adj=1)
  })
  
  output$enrichmentplot2 <- renderPlot({
    plot_data2()
  }, height = 500, width=600)
  
  output$download_report <- downloadHandler(
    filename = function() {
      paste("report-", Sys.Date(), ".html", sep="")
    },
    content = function(file) {
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      params <- list("res"=edata(),
                     "minsetsize"=input$minsetsize,
                     "prioritisation"=input$prioritisation)
      
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
