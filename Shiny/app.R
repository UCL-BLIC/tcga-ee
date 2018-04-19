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

# Initialise variables
  data.dir <- "../data/"
  cohorts <- NULL
  cohorts.TumourOnly <- NULL

  # Read in the ensembl to HUGO mapping
  if (file.exists(paste0(data.dir, "gencode.v22.annotation.gene.probeMap"))) {
    genes <- fread(paste0(data.dir, "gencode.v22.annotation.gene.probeMap"), header = T, data.table = T)
  } else if (file.exists(paste0(data.dir, "gencode.v22.annotation.gene.probeMap.gz"))) {
    genes <- fread(paste0("gunzip -c ", data.dir, "gencode.v22.annotation.gene.probeMap.gz"), header = T, data.table = T)
  } else {
    stop("Cannot find gencode annotation file")
  }

  # Obtain the list of files containing the FPKM-UQ data
  files <- list.files(data.dir, pattern = ".htseq_fpkm-uq")
  cohort_names <- sub(".htseq_fpkm-uq.*", "", files, perl = T)

# Primary function for reading in data from disk
read_data <- function(session, data.dir)
{
  these.cohorts <- list()
  these.cohorts.TumourOnly <- list()

  # Progress bar on bottom-right
  progress <- Progress$new(session, min = 0, max = length(files))
  progress$set(value = 0, message = 'Loading TCGA data...')
  on.exit(progress$close())

  # Loop through each tsv.gz file and read in the values
  # Return a list object containing 2 elements: matched T-N; all T only
  for (i in 1:length(files))
  {
    this.file <- files[i]
    this.cohort_name <- sub(".htseq_fpkm-uq.*", "", this.file, perl = T)
    progress$set(value = i, message = paste0('Loading ', this.cohort_name, '...'))
    these.cohorts[[this.cohort_name]]$file <- this.file

    if (grepl(".gz$", this.file, perl = T))
    {
      data <- fread(paste0("gunzip -c ", data.dir, "/", this.file), header = T, data.table = F)
    }
    else
    {
      data <- fread(paste0(data.dir, "/", this.file), header = T, data.table = F)
    }
    rownames(data) <- data[, 1]
    data <- data[, -1]

    # Subset to matching columns (i.e. tumour-normal based on case ID) (unique number of cases reported in parentheses after cohort name)
    data.T <- data[, grepl("-01\\w$", colnames(data), perl = T), drop = F]
    data.N <- data[, grepl("-11\\w$", colnames(data), perl = T), drop = F]
    cases <- intersect(sub("-\\d\\d\\w$", "", colnames(data.N), perl = T),
                        sub("-\\d\\d\\w$", "", colnames(data.T), perl = T))
    data.N <- data.N[, which(sub("-\\d\\d\\w$", "", colnames(data.N)) %in% cases), drop = F]
    data.T <- data.T[, which(sub("-\\d\\d\\w$", "", colnames(data.T)) %in% cases), drop = F]
    these.cohorts[[this.cohort_name]]$cases <- cases
    these.cohorts[[this.cohort_name]]$data.T <- data.T
    these.cohorts[[this.cohort_name]]$data.N <- data.N
    these.cohorts[[this.cohort_name]]$label <- paste0(this.cohort_name, " (", length(cases[!duplicated(gsub("-01\\w$", "", cases))]),")")

    # Subset to select all tumour only (unique number of cases reported in parentheses after cohort name)
    data.TumourOnly <- data[, grepl("-01\\w$", colnames(data), perl = T), drop = F]
    cases.TumourOnly <- colnames(data.TumourOnly)
    these.cohorts.TumourOnly[[this.cohort_name]]$cases <- cases.TumourOnly
    these.cohorts.TumourOnly[[this.cohort_name]]$data.T <- data.TumourOnly
    these.cohorts.TumourOnly[[this.cohort_name]]$label <- paste0(this.cohort_name, " (", length(cases.TumourOnly[!duplicated(gsub("-01\\w$", "", cases.TumourOnly))]),")")
  }

  cohorts <<- these.cohorts
  cohorts.TumourOnly <<- these.cohorts.TumourOnly

  return(list(cohorts, cohorts.TumourOnly))
}

# Function for displaying error messgaes in a ggplot2 environment
gg_message <- function(message)
{
  ggplot() +
    geom_text(label = message, aes(x = "x", y = "y")) +
    theme_void()
}

# Define the UI (universal interface) for the application
ui <- navbarPage(
  title = "TCGA Expression Explorer (TCGA-EE)",
  theme = shinytheme("sandstone"),

  # New tab for tumour-normal comparisons
  tabPanel("Matched tumour vs normal",
    sidebarLayout(
                  # Side-bar with dataset selector, gene selector, discrete / continous colour shading bar, and a slider for plot point size
                  sidebarPanel(
                              selectInput("cohort",
                                          "Cohort (number of matching T-N)",
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
                              p("Data downloaded from UCSC Xena", a(href = "http://xena.ucsc.edu", "http://xena.ucsc.edu"))),

                  # Main panel with different plots (as sub-tabs) of the generated distributions (top-hald) and a data-table of the displayed values (bottom-half)
                  mainPanel(
                            fluidRow(column(8, tabsetPanel(type = "tabs",
                                                          tabPanel("Parallel Coordinate Plot",
                                                                  plotOutput("TN.parallel", height = "400px", width = "400px")),
                                                          tabPanel("Scatter Plot",
                                                                  plotOutput("TN.scatter", height = "400px", width = "400px")),
                                                          tabPanel("Density Plot",
                                                                  plotOutput("TN.density", height = "400px", width = "400px")))),
                                    column(4, wellPanel(htmlOutput("TN.ttest")))),
                            h3("Matched tumour-normal TCGA RNA-seq FPKM-UQ data"),
                            DT::dataTableOutput("TN.dataTable")))),

  # New tab: correlation analysis
  tabPanel("Tumour correlation analysis",
          sidebarLayout(
                      # Side-bar with dataset selector, 2 x gene selectors, and a slider for plot point size
                      sidebarPanel(
                                  selectInput("cohort.TumourOnly",
                                              "Cohort (num tumour)",
                                              choices = "Loading...",
                                              # selected = ,
                                              multiple = F,
                                              selectize = F),
                                  selectInput("gene1",
                                              "Gene 1 (select or type)",
                                              choices = sort(genes$gene),
                                              selected = "CDK9",
                                              size = 10,
                                              multiple = F,
                                              selectize = F),
                                  selectInput("gene2",
                                              "Gene 2 (select or type)",
                                              choices = sort(genes$gene),
                                              selected = "CDK9",
                                              size = 10,
                                              multiple = F,
                                              selectize = F),
                                  sliderInput("pointsize.TumourOnly",
                                              "Point-size",
                                              min = 0.3, max = 5, step = 0.1, value = 2),
                                  hr(),
                                  p("Data downloaded from UCSC Xena", a(href = "http://xena.ucsc.edu", "http://xena.ucsc.edu"))),

                      # Main panel with a plot of the correlation/regression (top half) and a data-table of the displayed values (bottom half)
                      mainPanel(
                                fluidRow(column(8, tabsetPanel(type = "tabs",
                                                              tabPanel("Scatter Correlation Plot",
                                                                      plotOutput("TT.scatter.cor", height = "400px", width = "400px")))),
                                column(4, wellPanel(htmlOutput("TT.cortest")))),
                      h3("Tumour TCGA RNA-seq FPKM-UQ data"),
                      DT::dataTableOutput("TT.dataTable")))),

  # New panel with information about the application
  tabPanel("About",
          h3("Data"),
            p("The data in this app come from the GDC portal and have been downlaoded from the UCSC Xena Browser",
              a(href = "http://xena.ucsc.edu", "http://xena.ucsc.edu")),
            p("The list of cohorts is:"),
              pre("
                GDC TCGA Acute Myeloid Leukemia (LAML)
                GDC TCGA Adrenocortical Cancer (ACC)
                GDC TCGA Bile Duct Cancer (CHOL)
                GDC TCGA Bladder Cancer (BLCA)
                GDC TCGA Breast Cancer (BRCA)
                GDC TCGA Cervical Cancer (CESC)
                GDC TCGA Colon Cancer (COAD)
                GDC TCGA Endometrioid Cancer (UCEC)
                GDC TCGA Esophageal Cancer (ESCA)
                GDC TCGA Glioblastoma (GBM)
                GDC TCGA Head and Neck Cancer (HNSC)
                GDC TCGA Kidney Chromophobe (KICH)
                GDC TCGA Kidney Clear Cell Carcinoma (KIRC)
                GDC TCGA Kidney Papillary Cell Carcinoma (KIRP)
                GDC TCGA Large B-cell Lymphoma (DLBC)
                GDC TCGA Liver Cancer (LIHC)
                GDC TCGA Lower Grade Glioma (LGG)
                GDC TCGA Lung Adenocarcinoma (LUAD)
                GDC TCGA Lung Squamous Cell Carcinoma (LUSC)
                GDC TCGA Melanoma (SKCM)
                GDC TCGA Mesothelioma (MESO)
                GDC TCGA Ocular melanomas (UVM)
                GDC TCGA Ovarian Cancer (OV)
                GDC TCGA Pancreatic Cancer (PAAD)
                GDC TCGA Pheochromocytoma & Paraganglioma (PCPG)
                GDC TCGA Prostate Cancer (PRAD)
                GDC TCGA Rectal Cancer (READ)
                GDC TCGA Sarcoma (SARC)
                GDC TCGA Stomach Cancer (STAD)
                GDC TCGA Testicular Cancer (TGCT)
                GDC TCGA Thymoma (THYM)
                GDC TCGA Thyroid Cancer (THCA)
                GDC TCGA Uterine Carcinosarcoma (UCS)"))
)

# Define server logic
server <- function(input, output, session)
{
  # If there are no datasets populated (this will be the first task to be executed when app is started)
  if (is.null(cohorts) || is.null(cohorts.TumourOnly))
  {
    # Call function to read in the data
    data <- read_data(session, data.dir)

    # read_data returns a list with 2 elements: 1, matched tumour-normal; 2, tumour only
    cohorts <- data[[1]]
    cohorts.TumourOnly <- data[[2]]

    # Update the dataset selectors on both panels ("Matched tumour vs normal" and "Tumour correlation analysis") in the app
    cohort_choices <- names(cohorts)
    names(cohort_choices) <- sapply(cohorts, function(x) {x$label})
    updateSelectInput(session, "cohort", "Cohort", 
                       choices = cohort_choices,
                       selected = cohort_choices[grep("(0)", names(cohort_choices), value = F, invert = T)[1]])

    cohort_choices.TumourOnly <- names(cohorts.TumourOnly)
    names(cohort_choices.TumourOnly) <- sapply(cohorts.TumourOnly, function(x) {x$label})
    updateSelectInput(session, "cohort.TumourOnly", "Cohort", 
                      choices = cohort_choices.TumourOnly,
                      selected = cohort_choices.TumourOnly[grep("(0)", names(cohort_choices.TumourOnly), value = F, invert = T)[1]])
  }

  # Function that will be called for reading in the matched T-N data when a cohort and gene is selected from "Matched tumour vs normal" panel
  TN.data <- reactive(
  {
    cohort_name <- input$cohort
    gene <- as.character(genes[genes$gene == input$gene, 1])

    # If the cohort name is not in the list, return an empty object
    if (!(cohort_name %in% names(cohorts)))
    {
      data <- data.frame(matrix(NA, ncol = 4, nrow = 0))
      colnames(data) <- c("Case", "Normal", "Tumour", "diff")
      return(data)
    }

    # Extract the expression values for the given gene in the selected cohort for tumour and then normal
    data.T <- cohorts[[cohort_name]]$data.T[gene, , drop = F]
    data.N <- cohorts[[cohort_name]]$data.N[gene, , drop = F]

    # If no tumour data for this gene
    if (ncol(data.T) == 0)
    {
      data <- data.frame(matrix(NA, ncol = 4, nrow = 0))
      colnames(data) <- c("Case", "Normal", "Tumour", "diff")
      return(data)
    }

    # Match the tumour with normal by TCGA barcode ('case')
    data.N <- gather(data.N, key = "Case", value = "Normal")
    data.T <- gather(data.T, key = "Case", value = "Tumour")
    data.N$Case <- sub("-\\d\\d\\w$", "", data.N$Case)
    data.T$Case <- sub("-\\d\\d\\w$", "", data.T$Case)
    data <- inner_join(data.N, data.T, by = "Case")
    data$diff <- data$Tumour - data$Normal

    # Remove duplicated samples
    data <- data[!duplicated(data$Case),]

    return(data)
  })

  # Function that will be called for reading in the tumour data when a cohort and 2 genes are selected from "Tumour correlation analysis" panel
  TT.data <- reactive(
  {
    cohort_name.TumourOnly <- input$cohort.TumourOnly
    gene1 <- as.character(genes[genes$gene == input$gene1, 1])
    gene2 <- as.character(genes[genes$gene == input$gene2, 1])

    # If the cohort name is not in the list, return an empty object
    if (!(cohort_name.TumourOnly %in% names(cohorts.TumourOnly)))
    {
      data <- data.frame(matrix(NA, ncol = 4, nrow = 0))
      colnames(data) <- c("Case", "Gene1", "Gene2", "diff")
      return(data)
    }

    # Extract the expression values for the given genes in the selected cohort for tumour
    data.T1 <- cohorts.TumourOnly[[cohort_name.TumourOnly]]$data.T[gene1, , drop = F]
    data.T2 <- cohorts.TumourOnly[[cohort_name.TumourOnly]]$data.T[gene2, , drop = F]

    # If no tumour data for either of the selected genes
    if (ncol(data.T1) == 0 || ncol(data.T2) == 0)
    {
      data <- data.frame(matrix(NA, ncol = 4, nrow = 0))
      colnames(data) <- c("Case", "Gene1", "Gene2", "diff")
      return(data)
    }

    # Match the values for gene1 with gene2 by TCGA barcode ('case')
    data.T1 <- gather(data.T1, key = "Case", value = "Gene1")
    data.T2 <- gather(data.T2, key = "Case", value = "Gene2")
    data.T1$Case <- sub("-\\d\\d\\w$", "", data.T1$Case)
    data.T2$Case <- sub("-\\d\\d\\w$", "", data.T2$Case)
    data <- inner_join(data.T1, data.T2, by = "Case")
    colnames(data) <- c("Case", "Gene1", "Gene2")
    data$diff <- data$Gene1 - data$Gene2

    # Remove duplicated samples
    data <- data[!duplicated(data$Case),]

    return(data)
  })

  # Function for generating the parallel / bi-partite plot for "Matched tumour vs normal" panel
  output$TN.parallel <- renderPlot(
  {
    cohort_name <- input$cohort
    gene <- as.character(genes[genes$gene == input$gene, 1])

    data <- TN.data()

    if (nrow(data) == 0)
    {
      return(gg_message(paste("No T-N data for gene\n", input$gene, "in", input$cohort)))
    }

    if (input$scale == "Discrete")
    {
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

  # Function for generating the scatter plot for "Matched tumour vs normal" panel
  output$TN.scatter <- renderPlot(
  {
    cohort_name <- input$cohort
    gene <- as.character(genes[genes$gene == input$gene, 1])

    data <- TN.data()

    if (nrow(data) == 0)
    {
      return(gg_message(paste("No T-N data for gene\n", input$gene, "in", input$cohort)))
    }

    if (input$scale == "Discrete")
    {
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

  # Function for generating the density plot for "Matched tumour vs normal" panel
  output$TN.density <- renderPlot(
  {
    cohort_name <- input$cohort
    gene <- as.character(genes[genes$gene == input$gene, 1])

    data <- TN.data()

    if (nrow(data) == 0)
    {
      return(gg_message(paste("No T-N data for gene\n", input$gene, "in", input$cohort)))
    }

    ggplot(data = data) +
      geom_density(aes(x = `diff`), color = "black", fill = "lightblue") +
      geom_vline(linetype = 1, color = "black", xintercept = 0) +
      geom_vline(linetype = 2, color = "blue", xintercept = mean(data$diff)) +
      coord_cartesian(xlim = c(-max(abs(data$diff)), max(abs(data$diff)))) + 
      theme_linedraw() +
      labs(title = paste(input$gene, "in", input$cohort),
           x = "Matched tumour - normal\n(FPKM-UQ)")
  })

  # Function for generating the scatter plot for "Matched tumour vs normal" panel
  output$TN.ttest <- renderText(
  {
    data <- TN.data()

    # Convert zeros to NAs - these will be later removed
    data$Tumour[data$Tumour==0] <- NA
    data$Normal[data$Normal==0] <- NA

    # Check for complete matched cases, i.e., where both T-N have values, and require at least 3 of these
    if (nrow(data[complete.cases(data$Tumour, data$Normal),])<3)
    {
      out <- c("<p>Error - not enough complete observations</p>")
    }
    else
    {
      # Paired student's t-test (parametric)
      t <- t.test(data$diff)

      # Wilcoxon Signed Rank test (non-parametric)
      wt <- wilcox.test(data$Tumour, data$Normal, paired=TRUE, exact=FALSE, conf.int=TRUE, conf.level=0.95)

      # Record the number of complete cases
      n.complete.obs <- nrow(data[complete.cases(data$Tumour, data$Normal),])

      options(scipen=999)

      out <- c("<p align='left'><b>One-sample student's t-test</b></p>\n",
              "<p><i>Is mean difference not equal to 0?</i></p>\n",
              paste0("<p>&nbsp;&nbsp;&nbsp;Estimated mean: ", round(t$estimate, 3), "</p>\n"),
              paste0("<p>&nbsp;&nbsp;&nbsp;95% CI: ", round(t$conf.int[1], 3), ", ", round(t$conf.int[2], 3), "</p>\n"),
              paste0("<p>&nbsp;&nbsp;&nbsp;p-value: ", signif(t$p.value, 3), "</p>\n"))

      out <- c(out, "<p><hr></p>\n")
      out <- c(out, "<p align='left'><b>Wilcoxon signed rank test</b></p>\n",
               paste0("<p>&nbsp;&nbsp;&nbsp;Estimated location parameter: ", round(wt$estimate, 3), "</p>\n"),
               paste0("<p>&nbsp;&nbsp;&nbsp;95% CI: ", round(wt$conf.int[1], 3), ", ", round(wt$conf.int[2], 3), "</p>\n"),
               paste0("<p>&nbsp;&nbsp;&nbsp;p-value: ", signif(wt$p.value, 3), "</p>\n"))
      out <- c(out, "<p><hr></p>\n")
      out <- c(out, paste0("<p>Number of complete non-zero observatons = ", n.complete.obs, "</p>\n"))
      out <- c(out, "<p>----------------------</p>")
      out <- c(out, "<p><b>Notes on statistical testing:</b></p>",
               "<ul><li>observations with a zero not included</li>",
               "<li>minimum 3 complete observations required</li></ul>")
    }
  })

  # Function for outputting the expression values as a table in the main panel for "Matched tumour vs normal" panel
  output$TN.dataTable <- renderDataTable(
  {
    datatable(TN.data(),
              extensions = c('Buttons'),
              options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })

  # Function for outputting the expression values as a table in the main panel for "Tumour correlation analysis" panel
  output$TT.dataTable <- renderDataTable(
  {
    datatable(TT.data(),
              extensions = c('Buttons'),
              colnames=c("Case", input$gene1, input$gene2, "Diff"),
              options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })

  # Function to generate the scatter plot for correlation analysis, with linear regression fit ("Tumour correlation analysis" panel)
  output$TT.scatter.cor <- renderPlot(
  {
    cohort_name.TumourOnly <- input$cohort.TumourOnly

    gene1 <- as.character(genes[genes$gene == input$gene1, 1])
    gene2 <- as.character(genes[genes$gene == input$gene2, 1])

    data <- TT.data()

    # Convert zeros to NAs - as NA values, ggplot2 will automatically remove them (and report as a warning message)
    data$Gene1[data$Gene1==0] <- NA
    data$Gene2[data$Gene2==0] <- NA

    # Check for complete cases, i.e., where both gene1 and gene2 have values
    if (nrow(data[complete.cases(data$Gene1, data$Gene2),])==0)
    {
      return(gg_message(paste("No complete cases for chosen genes in", input$cohort.TumourOnly)))
    }
    else
    {
      ggplot(data = data, aes(x = `Gene1`, y = `Gene2`)) +
        geom_point(size = input$pointsize.TumourOnly) +
        geom_smooth(method='lm',formula=y~x, colour="hotpink", size=1.5) +
        stat_smooth(method="lm", fullrange=TRUE, colour="hotpink")  +
        theme_linedraw() +
        labs(title = paste("Gene correlation in", input$cohort.TumourOnly),
              x = paste("Gene1 (", input$gene1, ")\n(FPKM-UQ)", sep=""),
              y = paste("Gene2 (", input$gene2, ")\n(FPKM-UQ)", sep=""))
    }
  })

  # Function to display correlation results in the "Tumour correlation analysis" panel
  output$TT.cortest <- renderText(
  {
    data <- TT.data()

    # Convert zeros to NAs - these will be later removed
    data$Gene1[data$Gene1==0] <- NA
    data$Gene2[data$Gene2==0] <- NA

    # Check for complete matched cases, i.e., where both T-N have values, and require at least 3 of these
    if (nrow(data[complete.cases(data$Gene1, data$Gene2),])<3)
    {
      out <- c("<p align='left'>Error - not enough complete observations</p>\n")
    }
    else
    {
      # Pearson r correlation
      r.p <- cor(data$Gene1, data$Gene2, method="pearson", use="complete.obs")
      p.p <- cor.test(data$Gene1, data$Gene2, method="pearson", use="complete.obs")

      # Spearman rho correlation
      r.sp <- cor(data$Gene1, data$Gene2, method="spearman", use="complete.obs")
      p.sp <- cor.test(data$Gene1, data$Gene2, method="spearman", use="complete.obs", exact=FALSE)

      # Kendall tau correlation
      r.k <- cor(data$Gene1, data$Gene2, method="kendall", use="complete.obs")
      p.k <- cor.test(data$Gene1, data$Gene2, method="kendall", use="complete.obs", exact=FALSE)

      # Record the number of complete cases
      n.complete.obs <- nrow(data[complete.cases(data$Gene1, data$Gene2),])

      options(scipen=999)
      
      out <- c("<p align='left'><b>Correlation stats</b></p>\n",
                "<p></p>\n",
                paste0("<p>Pearson r = ", round(r.p, 3), " (95% CI ", round(p.p$conf.int[1], 3), ", ", round(p.p$conf.int[2], 3), ")</p>\n"),
                paste0("<p>p-value = ", round(p.p$p.value, 9), "</p>\n"),
                "<p><hr></p>\n",
                paste0("<p>Spearman &rho; (rho) = ", round(r.sp, 3), "</p>\n"),
                paste0("<p>p-value = ", round(p.sp$p.value, 9), "</p>\n"),
                "<p><hr></p>\n",
                paste0("<p>Kendall &tau; (tau) = ", round(r.k, 3), "</p>\n"),
                paste0("<p>p-value = ", round(p.k$p.value, 9), "</p>\n"),
                "<p><hr></p>\n",
                paste0("<p>Number of complete non-zero observatons = ", n.complete.obs, "</p>\n"))
      out <- c(out, "<p>----------------------</p>")
      out <- c(out, "<p><b>Notes on statistical testing:</b></p>",
              "<ul><li>observations with a zero removed</li>",
              "<li>minimum 3 complete observations required</li></ul>")
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
