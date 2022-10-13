
### Here, we will try to download some patient omics data for downstream machine learning applications.

### About GISTIC data:

# What is GISTIC score?
# GISTIC is a tool to identify genes targeted by somatic copy number variation (CNV). 
# The GISTIC algorithm defines CNV boundaries by a user-defined confidence level.

# For legacy data (data aligned to hg19) TCGAbiolinks is using GRCh37.p13 
# and for harmonized data (data aligned to hg38) now it is using Gencode version 36.

## I will download data for Lung adenocarcinoma from tcga database:


### I need to retrieve the sample TCGA barcodes. For this, I will use Firebrowser package in R:

getwd()

# [1] "C:/Users/ASUS/Documents"

library(FirebrowseR)

cohorts <- Metadata.Cohorts(format = "csv") ### Download all cohorts that are available.

cancer.Type = cohorts[grep("Lung adenocarcinoma", cohorts$description, ignore.case = T), 1] ### find the cohort with decription containing "pancreatic":

print(cancer.Type)

# [1] "LUAD"

# Now that we know that the breast cancer samples are identified be  "LUAD", we can retrieve a list
# of all patients associated with this identifier.

luad.Pats <- Samples.Clinical(cohort = cancer.Type, format="tsv")

dim(luad.Pats)

###[1] 150  96 ### there are 150 patients. However, there are more than this (631) according to (http://firebrowse.org/) or (https://gdac.broadinstitute.org/).

### The reason is that Firebrowse API returns 150 entries per page. So other entries are in next pages. Following script takes all the data from the following pages (for detail, see firebrowseR vignette):

all.Received = F                         ### this boolean variable is set initially FALSE and it enters the loop till it is False, and after we receive all data, it will be set to TRUE and it will exit from the while loop.
page.Counter = 1
page.size = 150                          ### number of maximum entry within each page. 
luad.Pats = list()                      ### Patients data from each page are stored in the list.

while(all.Received == F){
  luad.Pats[[page.Counter]] = Samples.Clinical(format = "csv",
                                                   cohort = cancer.Type,
                                                   page_size = page.size,
                                                   page = page.Counter)
  if(page.Counter > 1)
    colnames(luad.Pats[[page.Counter]]) = colnames(luad.Pats[[page.Counter-1]])  ### In the API, except first page,  pages do not have column name, so it adds the column name of from the data of first page stored as the first element of the list. 
  if(nrow(luad.Pats[[page.Counter]]) < page.size){    ### if the last page whose entries are received has less than 150 entries, it means that itis the last page. thus, all.Received boolean variable is set to TRUE and exits the loop.           
    all.Received = T
  } else{      ### if not true, it incerements the page.counter variable and continues with the retrieval of the patient data in the next page. 
    page.Counter = page.Counter + 1
  }
}

luad.Pats <- do.call(what = rbind , args = luad.Pats) ### do.call () function takes the function in "what" parameter and applies it to the list object assigned to the "args" parameter.

dim(luad.Pats)
# 522  81

### we will create a vector with the patients barcode:

sample_ids <- luad.Pats$tcga_participant_barcode

library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq",
  file.type = "normalized_results",
  experimental.strategy = "RNA-Seq",
  barcode = sample_ids,
  legacy = TRUE
)

library(dplyr)

query %>% head

### my current directory is "C:\Users\ASUS\Documents". Download worked only in this directory. I do not know the reason.

GDCdownload(
  query = query,
  method = "api",
  files.per.chunk = 10
)

data_gexp <- GDCprepare(query = query)

datatable(
  as.data.frame(colData(data)),
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
  rownames = FALSE)

### gene matrix:

datatable(
  assay(data)[1:20,],
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
  rownames = TRUE
)

## CNV data:


query_CNV <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Copy Number Variation",
  data.type = "Copy Number Segment",
  barcode = sample_ids
)

GDCdownload(query = query_CNV)
data_CNV <- GDCprepare(query = query_CNV)
table(data_CNV$Sample) %>% names() %>% unique()  ### more than 20 patients. I will check the command.

## SNV data:table(maf$Tumor_Sample_Barcode) %>% names() %>% length()

query_SNV <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open",barcode = sample_ids)

GDCdownload(query_SNV)

maf <- GDCprepare(query_SNV)

### check the number of unique samples:

table(maf$Tumor_Sample_Barcode) %>% names() %>% length()

session_info()

# Session info ──────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.2.1 (2022-06-23 ucrt)
# os       Windows 10 x64 (build 19043)
# system   x86_64, mingw32
# ui       RStudio
# language (EN)
# collate  Turkish_Turkey.utf8
# ctype    Turkish_Turkey.utf8
# tz       Europe/Istanbul
# date     2022-09-18
# rstudio  2022.07.1+554 Spotted Wakerobin (desktop)
# pandoc   NA
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# AnnotationDbi          1.58.0    2022-04-26 [1] Bioconductor
# assertthat             0.2.1     2019-03-21 [1] CRAN (R 4.1.2)
# Biobase              * 2.56.0    2022-04-26 [1] Bioconductor
# BiocFileCache          2.4.0     2022-04-26 [1] Bioconductor
# BiocGenerics         * 0.42.0    2022-04-26 [1] Bioconductor
# BiocManager          * 1.30.18   2022-05-18 [1] CRAN (R 4.2.1)
# biomaRt                2.52.0    2022-04-26 [1] Bioconductor
# Biostrings             2.64.1    2022-08-25 [1] Bioconductor
# bit                    4.0.4     2020-08-04 [1] CRAN (R 4.1.2)
# bit64                  4.0.5     2020-08-30 [1] CRAN (R 4.1.2)
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.1.1)
# blob                   1.2.3     2022-04-10 [1] CRAN (R 4.2.1)
# bslib                  0.4.0     2022-07-16 [1] CRAN (R 4.2.1)
# cachem                 1.0.6     2021-08-19 [1] CRAN (R 4.1.2)
# callr                  3.7.2     2022-08-22 [1] CRAN (R 4.2.1)
# cli                    3.4.0     2022-09-08 [1] CRAN (R 4.2.1)
# colorspace             2.0-3     2022-02-21 [1] CRAN (R 4.1.3)
# crayon                 1.5.1     2022-03-26 [1] CRAN (R 4.1.3)
# crosstalk              1.2.0     2021-11-04 [1] CRAN (R 4.1.2)
# curl                   4.3.2     2021-06-23 [1] CRAN (R 4.1.2)
# data.table             1.14.2    2021-09-27 [1] CRAN (R 4.1.3)
# DBI                    1.1.3     2022-06-18 [1] CRAN (R 4.2.1)
# dbplyr                 2.2.1     2022-06-27 [1] CRAN (R 4.2.1)
# DelayedArray           0.22.0    2022-04-26 [1] Bioconductor
# devtools             * 2.4.4     2022-07-20 [1] CRAN (R 4.2.1)
# digest                 0.6.29    2021-12-01 [1] CRAN (R 4.1.2)
# downloader             0.4       2015-07-09 [1] CRAN (R 4.2.1)
# dplyr                * 1.0.10    2022-09-01 [1] CRAN (R 4.2.1)
# DT                   * 0.24      2022-08-09 [1] CRAN (R 4.2.1)
# ellipsis               0.3.2     2021-04-29 [1] CRAN (R 4.1.2)
# fansi                  1.0.3     2022-03-24 [1] CRAN (R 4.1.3)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.2.1)
# fastmap                1.1.0     2021-01-25 [1] CRAN (R 4.1.2)
# filelock               1.0.2     2018-10-05 [1] CRAN (R 4.1.2)
# FirebrowseR          * 1.1.35    2022-09-17 [1] Github (mariodeng/FirebrowseR@c549f28)
# fs                     1.5.2     2021-12-08 [1] CRAN (R 4.1.2)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.32.4    2022-09-06 [1] Bioconductor
# GenomeInfoDbData       1.2.8     2022-09-12 [1] Bioconductor
# GenomicRanges        * 1.48.0    2022-04-26 [1] Bioconductor
# ggplot2              * 3.3.6     2022-05-03 [1] CRAN (R 4.2.1)
# glue                   1.6.2     2022-02-24 [1] CRAN (R 4.1.3)
# gtable                 0.3.1     2022-09-01 [1] CRAN (R 4.2.1)
# highr                  0.9       2021-04-16 [1] CRAN (R 4.1.3)
# hms                    1.1.2     2022-08-19 [1] CRAN (R 4.2.1)
# htmltools              0.5.3     2022-07-18 [1] CRAN (R 4.2.1)
# htmlwidgets            1.5.4     2021-09-08 [1] CRAN (R 4.1.2)
# httpuv                 1.6.6     2022-09-08 [1] CRAN (R 4.2.1)
# httr                   1.4.4     2022-08-17 [1] CRAN (R 4.2.1)
# IRanges              * 2.30.1    2022-08-25 [1] Bioconductor
# jquerylib              0.1.4     2021-04-26 [1] CRAN (R 4.1.2)
# jsonlite               1.8.0     2022-02-22 [1] CRAN (R 4.1.3)
# KEGGREST               1.36.3    2022-07-14 [1] Bioconductor
# knitr                  1.40      2022-08-24 [1] CRAN (R 4.2.1)
# labeling               0.4.2     2020-10-20 [1] CRAN (R 4.1.1)
# later                  1.3.0     2021-08-18 [1] CRAN (R 4.1.2)
# lattice                0.20-45   2021-09-22 [2] CRAN (R 4.2.1)
# lifecycle              1.0.2     2022-09-09 [1] CRAN (R 4.2.1)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.2.1)
# Matrix                 1.5-0     2022-09-10 [1] CRAN (R 4.2.1)
# MatrixGenerics       * 1.8.1     2022-06-30 [1] Bioconductor
# matrixStats          * 0.62.0    2022-04-19 [1] CRAN (R 4.2.1)
# memoise                2.0.1     2021-11-26 [1] CRAN (R 4.1.2)
# mime                   0.12      2021-09-28 [1] CRAN (R 4.1.1)
# miniUI                 0.1.1.1   2018-05-18 [1] CRAN (R 4.1.3)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.1.2)
# pillar                 1.8.1     2022-08-19 [1] CRAN (R 4.2.1)
# pkgbuild               1.3.1     2021-12-20 [1] CRAN (R 4.1.3)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.1.2)
# pkgload                1.3.0     2022-06-27 [1] CRAN (R 4.2.1)
# plyr                   1.8.7     2022-03-24 [1] CRAN (R 4.1.3)
# png                    0.1-7     2013-12-03 [1] CRAN (R 4.1.1)
# prettyunits            1.1.1     2020-01-24 [1] CRAN (R 4.1.2)
# processx               3.7.0     2022-07-07 [1] CRAN (R 4.2.1)
# profvis                0.3.7     2020-11-02 [1] CRAN (R 4.2.1)
# progress               1.2.2     2019-05-16 [1] CRAN (R 4.1.2)
# promises               1.2.0.1   2021-02-11 [1] CRAN (R 4.1.2)
# ps                     1.7.1     2022-06-18 [1] CRAN (R 4.2.1)
# purrr                  0.3.4     2020-04-17 [1] CRAN (R 4.1.2)
# R.methodsS3            1.8.2     2022-06-13 [1] CRAN (R 4.2.0)
# R.oo                   1.25.0    2022-06-12 [1] CRAN (R 4.2.0)
# R.utils                2.12.0    2022-06-28 [1] CRAN (R 4.2.1)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.1.2)
# rappdirs               0.3.3     2021-01-31 [1] CRAN (R 4.1.2)
# Rcpp                   1.0.9     2022-07-08 [1] CRAN (R 4.2.1)
# RCurl                  1.98-1.8  2022-07-30 [1] CRAN (R 4.2.1)
# readr                  2.1.2     2022-01-30 [1] CRAN (R 4.1.3)
# remotes                2.4.2     2021-11-30 [1] CRAN (R 4.1.3)
# rlang                  1.0.5     2022-08-31 [1] CRAN (R 4.2.1)
# rprojroot              2.0.3     2022-04-02 [1] CRAN (R 4.2.1)
# RSQLite                2.2.17    2022-09-10 [1] CRAN (R 4.2.1)
# rstudioapi             0.14      2022-08-22 [1] CRAN (R 4.2.1)
# rvest                  1.0.3     2022-08-19 [1] CRAN (R 4.1.3)
# S4Vectors            * 0.34.0    2022-04-26 [1] Bioconductor
# sass                   0.4.2     2022-07-16 [1] CRAN (R 4.2.1)
# scales                 1.2.1     2022-08-20 [1] CRAN (R 4.2.1)
# sessioninfo            1.2.2     2021-12-06 [1] CRAN (R 4.1.3)
# shiny                  1.7.2     2022-07-19 [1] CRAN (R 4.2.1)
# stringi                1.7.8     2022-07-11 [1] CRAN (R 4.2.1)
# stringr                1.4.1     2022-08-20 [1] CRAN (R 4.2.1)
# SummarizedExperiment * 1.26.1    2022-04-29 [1] Bioconductor
# TCGAbiolinks         * 2.24.3    2022-06-16 [1] Bioconductor
# TCGAbiolinksGUI.data   1.16.0    2022-04-28 [1] Bioconductor
# tibble                 3.1.8     2022-07-22 [1] CRAN (R 4.1.3)
# tidyr                  1.2.1     2022-09-08 [1] CRAN (R 4.2.1)
# tidyselect             1.1.2     2022-02-21 [1] CRAN (R 4.1.3)
# tzdb                   0.3.0     2022-03-28 [1] CRAN (R 4.1.3)
# urlchecker             1.0.1     2021-11-30 [1] CRAN (R 4.2.1)
# usethis              * 2.1.6     2022-05-25 [1] CRAN (R 4.2.1)
# utf8                   1.2.2     2021-07-24 [1] CRAN (R 4.1.2)
# vctrs                  0.4.1     2022-04-13 [1] CRAN (R 4.1.3)
# vroom                  1.5.7     2021-11-30 [1] CRAN (R 4.1.3)
# withr                  2.5.0     2022-03-03 [1] CRAN (R 4.1.3)
# xfun                   0.32      2022-08-10 [1] CRAN (R 4.2.1)
# XML                    3.99-0.10 2022-06-09 [1] CRAN (R 4.2.1)
# xml2                   1.3.3     2021-11-30 [1] CRAN (R 4.1.2)
# xtable                 1.8-4     2019-04-21 [1] CRAN (R 4.1.2)
# XVector                0.36.0    2022-04-26 [1] Bioconductor
# yaml                   2.3.5     2022-02-21 [1] CRAN (R 4.1.2)
# zlibbioc               1.42.0    2022-04-26 [1] Bioconductor
# 
# [1] C:/Users/ASUS/AppData/Local/R/win-library/4.2
# [2] C:/Program Files/R/R-4.2.1/library
