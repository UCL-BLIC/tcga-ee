#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(data.table)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)

data.dir <- "../data/"
cohorts <- NULL
genes <- fread(paste0("gunzip -c ", data.dir, "gencode.v22.annotation.gene.probeMap.gz"), header = T, data.table = T)
files <- list.files(data.dir, pattern = ".htseq_fpkm-uq")
cohort_names <- sub(".htseq_fpkm-uq.*", "", files, perl = T)

read_data <- function(session, data.dir) {
  these.cohorts <- list()
  progress <- Progress$new(session, min = 0, max = length(files))
  progress$set(value = 0, message = 'Loading TCGA data...')
  on.exit(progress$close())
  for (i in 1:length(files)) {
    this.file <- files[i]
    this.cohort_name <- sub(".htseq_fpkm-uq.*", "", this.file, perl = T)
    progress$set(value = i, message = paste0('Loading ', this.cohort_name, '...'))
    these.cohorts[[this.cohort_name]]$file <- this.file
    if (grepl(".gz$", this.file, perl = T)) {
      data <- fread(paste0("gunzip -c ", data.dir, "/", this.file), header = T, data.table = F)
    } else {
      data <- fread(paste0(data.dir, "/", this.file), header = T, data.table = F)
    }
    rownames(data) <- data[, 1]
    data <- data[, -1]

    # print(paste(this.cohort_name, paste(unique(sub("TCGA-..-....-", "", colnames(data))), collapse = " ")))
    
    ## Sample type codes are defined here: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
    
    ## Split into T and N samples (only solid tumours and only primaries)
    data.T <- data[, grepl("-01\\w$", colnames(data), perl = T), drop = F]
    data.N <- data[, grepl("-11\\w$", colnames(data), perl = T), drop = F]
    
    ## Subset to matching columns (i.e. Case) only
    cases <- intersect(sub("-\\d\\d\\w$", "", colnames(data.N), perl = T),
                          sub("-\\d\\d\\w$", "", colnames(data.T), perl = T))
    data.N <- data.N[, which(sub("-\\d\\d\\w$", "", colnames(data.N)) %in% cases), drop = F]
    data.T <- data.T[, which(sub("-\\d\\d\\w$", "", colnames(data.T)) %in% cases), drop = F]

    these.cohorts[[this.cohort_name]]$cases <- cases
    
    these.cohorts[[this.cohort_name]]$data.T <- data.T
    these.cohorts[[this.cohort_name]]$data.N <- data.N
    
    these.cohorts[[this.cohort_name]]$label <- paste0(this.cohort_name, " (", length(cases),")")
  }
  cohorts <<- these.cohorts
  return(cohorts)
}

gg_message <- function(message) {
  ggplot() +
    geom_text(label = message, aes(x = "x", y = "y")) +
    theme_void()
}


# Define UI for application that draws a histogram
ui <- navbarPage(
  title = "TCGA Expression Explorer",
  theme = shinytheme("sandstone"),

  # Tumour-Normal comparisons
  tabPanel("Tumour vs Normal",
  
   # Sidebar with a slider input for number of bins 
    sidebarLayout(

      sidebarPanel(
        selectInput("cohort",
                    "Cohort (num of matching T-N)",
                    choices = "Loading...",
                    # selected = ,
                    multiple = F,
                    selectize = F),
        selectInput("gene",
                    "Gene (select or type)",
                    choices = sort(genes$gene),
                    selected = "CDK9",
                    size = 10,
                    multiple = F,
                    selectize = F),
        selectInput("scale",
                    "Colourscale",
                    choices = c("Discrete", "Continuous"),
                    selected = "Continuous"),
        sliderInput("pointsize",
                    "Point-size",
                    min = 0.3, max = 5, step = 0.1, value = 2),
        hr(),
        p("Data downloaded from UCSC Xena",
          a(href = "http://xena.ucsc.edu", "http://xena.ucsc.edu"))
      ),
  
      # Show a plot of the generated distribution
      mainPanel(
        fluidRow(
          column(8, tabsetPanel(type = "tabs",
            tabPanel("Parallel Coordinate Plot",
                     plotOutput("TN.parallel", height = "400px", width = "400px")
            ),
            tabPanel("Scatter Plot",
                     plotOutput("TN.scatter", height = "400px", width = "400px")
            ),
            tabPanel("Density Plot",
                     plotOutput("TN.density", height = "400px", width = "400px")
            )
          )),
          column(4, wellPanel(htmlOutput("TN.ttest")))
        ),
        h3("Normal-Tumour FPKM-UQ data"),
        DT::dataTableOutput("TN.dataTable")
      )
    )
  ),
  # Tumour-Normal comparisons
  tabPanel(
    "About",
    h3("Data"),
    p("The data in this app come from the GDC portal and have been downlaoded from the UCSC Xena Browser",
      a(href = "http://xena.ucsc.edu", "http://xena.ucsc.edu")),
    p("The list of cohorts is:"),
    pre("
    GDC TCGA Acute Myeloid Leukemia (LAML) (11 datasets)
    GDC TCGA Adrenocortical Cancer (ACC) (11 datasets)
    GDC TCGA Bile Duct Cancer (CHOL) (11 datasets)
    GDC TCGA Bladder Cancer (BLCA) (11 datasets)
    GDC TCGA Breast Cancer (BRCA) (11 datasets)
    GDC TCGA Cervical Cancer (CESC) (11 datasets)
    GDC TCGA Colon Cancer (COAD) (11 datasets)
    GDC TCGA Endometrioid Cancer (UCEC) (11 datasets)
    GDC TCGA Esophageal Cancer (ESCA) (11 datasets)
    GDC TCGA Glioblastoma (GBM) (11 datasets)
    GDC TCGA Head and Neck Cancer (HNSC) (11 datasets)
    GDC TCGA Kidney Chromophobe (KICH) (11 datasets)
    GDC TCGA Kidney Clear Cell Carcinoma (KIRC) (11 datasets)
    GDC TCGA Kidney Papillary Cell Carcinoma (KIRP) (11 datasets)
    GDC TCGA Large B-cell Lymphoma (DLBC) (11 datasets)
    GDC TCGA Liver Cancer (LIHC) (11 datasets)
    GDC TCGA Lower Grade Glioma (LGG) (11 datasets)
    GDC TCGA Lung Adenocarcinoma (LUAD) (11 datasets)
    GDC TCGA Lung Squamous Cell Carcinoma (LUSC) (11 datasets)
    GDC TCGA Melanoma (SKCM) (11 datasets)
    GDC TCGA Mesothelioma (MESO) (11 datasets)
    GDC TCGA Ocular melanomas (UVM) (11 datasets)
    GDC TCGA Ovarian Cancer (OV) (11 datasets)
    GDC TCGA Pancreatic Cancer (PAAD) (11 datasets)
    GDC TCGA Pheochromocytoma & Paraganglioma (PCPG) (11 datasets)
    GDC TCGA Prostate Cancer (PRAD) (11 datasets)
    GDC TCGA Rectal Cancer (READ) (11 datasets)
    GDC TCGA Sarcoma (SARC) (11 datasets)
    GDC TCGA Stomach Cancer (STAD) (11 datasets)
    GDC TCGA Testicular Cancer (TGCT) (11 datasets)
    GDC TCGA Thymoma (THYM) (11 datasets)
    GDC TCGA Thyroid Cancer (THCA) (11 datasets)
    GDC TCGA Uterine Carcinosarcoma (UCS) (11 datasets)")
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  if (is.null(cohorts)) {
    cohorts <- read_data(session, data.dir)
    cohort_choices <- names(cohorts)
    names(cohort_choices) <- sapply(cohorts, function(x) {x$label})
    updateSelectInput(session, "cohort", "Cohort", 
                      choices = cohort_choices,
                      selected = cohort_choices[grep("(0)", names(cohort_choices), value = F, invert = T)[1]])

  }

  TN.data <- reactive({
    cohort_name <- input$cohort
    gene <- as.character(genes[genes$gene == input$gene, 1])

    if (!(cohort_name %in% names(cohorts))) {
      data <- data.frame(matrix(NA, ncol = 4, nrow = 0))
      colnames(data) <- c("Case", "Normal", "Tumour", "diff")
      return(data)
    }
    data.T <- cohorts[[cohort_name]]$data.T[gene, , drop = F]
    data.N <- cohorts[[cohort_name]]$data.N[gene, , drop = F]
    if (ncol(data.T) == 0) {
      data <- data.frame(matrix(NA, ncol = 4, nrow = 0))
      colnames(data) <- c("Case", "Normal", "Tumour", "diff")
      return(data)
    }
    data.N <- gather(data.N, key = "Case", value = "Normal")
    data.T <- gather(data.T, key = "Case", value = "Tumour")
    data.N$Case <- sub("-\\d\\d\\w$", "", data.N$Case)
    data.T$Case <- sub("-\\d\\d\\w$", "", data.T$Case)
    data <- inner_join(data.N, data.T, by = "Case")
    data$diff <- data$Tumour - data$Normal
    
    return(data)
  })

  output$TN.parallel <- renderPlot({
    cohort_name <- input$cohort
    gene <- as.character(genes[genes$gene == input$gene, 1])

    data <- TN.data()
    if (nrow(data) == 0) {
      return(gg_message(paste("No T-N data for gene\n", input$gene, "in", input$cohort)))
    }

    if (input$scale == "Discrete") {
      data$diff <- sign(data$diff)
    }
    limits <- max(abs(range(data$diff)))
    limits <- c(-limits, limits)

    ## id column is to sort out cases with more than one sample
    data$id <- 1:nrow(data)
    data <- gather(data, key = "Sample", value = "FPKM-UQ", -c("Case", "diff", "id"))

    ggplot(data = data, aes(x = `Sample`, y = `FPKM-UQ`, group = id, color = diff)) +
      geom_path(lineend = "round") + 
      geom_point(size = input$pointsize) +
      scale_color_gradientn(colors = c("blue", "darkblue", "black", "darkred", "red"),
                            values = c(0, 0.4, 0.5, 0.6, 1),
                            limits = limits,
                            guide = F) +
      theme_linedraw() +
      labs(title = paste(input$gene, "in", input$cohort))
  })

  output$TN.scatter <- renderPlot({
    cohort_name <- input$cohort
    gene <- as.character(genes[genes$gene == input$gene, 1])

    data <- TN.data()
    if (nrow(data) == 0) {
      return(gg_message(paste("No T-N data for gene\n", input$gene, "in", input$cohort)))
    }

    if (input$scale == "Discrete") {
      data$diff <- sign(data$diff)
    }
    limits <- max(abs(range(data$diff)))
    limits <- c(-limits, limits)
    
    ggplot(data = data, aes(x = `Normal`, y = `Tumour`, color = diff)) +
      geom_abline(linetype = 1, color = "black") +
      geom_point(size = input$pointsize) +
      coord_fixed(ratio = 1, xlim = range(data[, 2:3]), ylim = range(data[, 2:3])) +
      scale_color_gradientn(colors = c("blue", "darkblue", "black", "darkred", "red"),
                            values = c(0, 0.4, 0.5, 0.6, 1),
                            limits = limits,
                            guide = F) +
      theme_linedraw() +
      labs(title = paste(input$gene, "in", input$cohort),
           x = "Normal (FPKM-UQ)", y = "Tumour (FPKM-UQ)")
  })

  output$TN.density <- renderPlot({
    cohort_name <- input$cohort
    gene <- as.character(genes[genes$gene == input$gene, 1])

    data <- TN.data()
    if (nrow(data) == 0) {
      return(gg_message(paste("No T-N data for gene\n", input$gene, "in", input$cohort)))
    }
    
    ggplot(data = data) +
      geom_density(aes(x = `diff`), color = "black", fill = "lightblue") +
      geom_vline(linetype = 1, color = "black", xintercept = 0) +
      geom_vline(linetype = 2, color = "blue", xintercept = mean(data$diff)) +
      theme_linedraw() +
      labs(title = paste(input$gene, "in", input$cohort),
           x = "Tumour - Normal (FPKM-UQ)")

    # Bit of code to create a density plot with a gradient filling
    # y <- density(data$diff, n = 2000)
    # limits <- max(abs(range(y$x)))
    # limits <- c(-limits, limits)
    # 
    # ggplot(data.frame(x = y$x, y = y$y), aes(x, y)) +
    #   geom_segment(aes(xend = x, yend = 0, colour = x)) + 
    #   geom_line(col = "black") + 
    #   scale_color_gradientn(colors = c("blue", "darkblue", "black", "darkred", "red"),
    #                         values = c(0, 0.4, 0.5, 0.6, 1),
    #                         limits = limits,
    #                         guide = F) +
    #   labs(title = paste(input$gene, "in", input$cohort),
    #        x = "Tumour - Normal (FPKM-UQ)")
  })
  
  output$TN.ttest <- renderText({
    data <- TN.data()
    if (nrow(data) > 0) {
      t <- t.test(data$diff)
      out <- c("<p align='center'><b>One Sample t-test</b></p>\n",
               "<p>Is mean not equal to 0?</p>\n",
               "<p></p>\n",
               paste0("<p>mean: ", round(t$estimate, 3), "</p>\n"),
               "<p></p>\n",
               "<p>95% conf.int:</p>\n",
               paste0("&nbsp;&nbsp;&nbsp;&nbsp;", round(t$conf.int[1], 3), "&nbsp;&nbsp;", round(t$conf.int[2], 3), "</p>\n"),
               "<p></p>\n")
      if (!is.nan(t$p.value) & t$p.value < 0.05) {
          out <- c(out, paste("<p><b>p-value:", signif(t$p.value, 3), "</b></p>\n"))
      } else {
        out <- c(out, paste("<p>p-value:", signif(t$p.value, 3), "</p>\n"))
      }
    }
  })
  
  output$TN.ttest2 <- renderPrint({
    data <- TN.data()
    if (nrow(data) > 0) {
      t.test(data$diff)
    }
  })
  
  output$TN.dataTable <- renderDataTable({
    datatable(TN.data(),
              extensions = c('Buttons'),
              options = list(
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
              ))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

