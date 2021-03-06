# TCGA Expression Explorer

Shiny app to explore the expression data from TCGA.

## Functionality

- [x] Explore paired Normal-Tumour expression data
- [x] Conduct correlation between 2 genes within the same tumour

## Description

This app loads the TCGA expression data in memory and allows the user to explore them. At the moment, the only functionality implemented is to look at the expression of the tumour and normal paired samples.

## Installation

1. Clone the repository on your local machine and open the code with RStudio.
2. Create a data/ folder
3. Obtain necessary data frm BLIC
4. Open the Shiny/app.R file
5. Install necessary packages
6. Run app by clicking on "Run App"

## Data

Originally, FPKM-UQ data from [UCSC Xena browser](http://xena.ucsc.edu) (see Appendix) was used; however, raw RSEM gene-summarised counts were then obtained from the [GDC Legacy Archive](https://portal.gdc.cancer.gov/legacy-archive) and re-processd at BLIC. New data files are available by request to Javier Herrero or Kevin Blighe.

<hr>

## Appendix

Original FPKM-UQ data used:

```bash
wget -N https://gdc.xenahubs.net/download/probeMaps/gencode.v22.annotation.gene.probeMap.gz
wget -N https://gdc.xenahubs.net/download/probeMaps/gencode.v22.annotation.gene.probeMap.gz
wget -N https://gdc.xenahubs.net/download/TCGA-LAML/Xena_Matrices/TCGA-LAML.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-ACC/Xena_Matrices/TCGA-ACC.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-CHOL/Xena_Matrices/TCGA-CHOL.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-BLCA/Xena_Matrices/TCGA-BLCA.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-BRCA/Xena_Matrices/TCGA-BRCA.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-CESC/Xena_Matrices/TCGA-CESC.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-COAD/Xena_Matrices/TCGA-COAD.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-UCEC/Xena_Matrices/TCGA-UCEC.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-ESCA/Xena_Matrices/TCGA-ESCA.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-GBM/Xena_Matrices/TCGA-GBM.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-HNSC/Xena_Matrices/TCGA-HNSC.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-KICH/Xena_Matrices/TCGA-KICH.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-KIRC/Xena_Matrices/TCGA-KIRC.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-KIRP/Xena_Matrices/TCGA-KIRP.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-DLBC/Xena_Matrices/TCGA-DLBC.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-LIHC/Xena_Matrices/TCGA-LIHC.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-LGG/Xena_Matrices/TCGA-LGG.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-LUAD/Xena_Matrices/TCGA-LUAD.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-LUSC/Xena_Matrices/TCGA-LUSC.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-SKCM/Xena_Matrices/TCGA-SKCM.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-MESO/Xena_Matrices/TCGA-MESO.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-UVM/Xena_Matrices/TCGA-UVM.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-OV/Xena_Matrices/TCGA-OV.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-PAAD/Xena_Matrices/TCGA-PAAD.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-PCPG/Xena_Matrices/TCGA-PCPG.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-PRAD/Xena_Matrices/TCGA-PRAD.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-READ/Xena_Matrices/TCGA-READ.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-SARC/Xena_Matrices/TCGA-SARC.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-STAD/Xena_Matrices/TCGA-STAD.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-TGCT/Xena_Matrices/TCGA-TGCT.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-THYM/Xena_Matrices/TCGA-THYM.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-THCA/Xena_Matrices/TCGA-THCA.htseq_fpkm-uq.tsv.gz
wget -N https://gdc.xenahubs.net/download/TCGA-UCS/Xena_Matrices/TCGA-UCS.htseq_fpkm-uq.tsv.gz
```
