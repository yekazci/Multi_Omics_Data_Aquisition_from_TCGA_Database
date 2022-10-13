
### Here, we will try to download some patient omics data for downstream machine learning applications.

### About GISTIC data:

# What is GISTIC score?
# GISTIC is a tool to identify genes targeted by somatic copy number variation (CNV). 
# The GISTIC algorithm defines CNV boundaries by a user-defined confidence level.

# For legacy data (data aligned to hg19) TCGAbiolinks is using GRCh37.p13 
# and for harmonized data (data aligned to hg38) now it is using Gencode version 36.

## I will download data for Colorectal adenocarcinoma from tcga database.


### I need to retrieve the sample TCGA barcodes. For this, I will use Firebrowser package in R:


library(FirebrowseR)

cohorts <- Metadata.Cohorts(format = "csv") ### Download all cohorts that are available.

cancer.Type = cohorts[grep("colorectal", cohorts$description, ignore.case = T), 1] ### find the cohort with decription containing "pancreatic":

print(cancer.Type)

# [1] "COADREAD"

# Now that we know that the breast cancer samples are identified be  "COADREAD", we can retrieve a list
# of all patients associated with this identifier.

coadread.Pats <- Samples.Clinical(cohort = cancer.Type, format="tsv")

dim(coadread.Pats)

###[1] 150  96 ### there are 150 patients. However, there are more than this (631) according to (http://firebrowse.org/) or (https://gdac.broadinstitute.org/).

### The reason is that Firebrowse API returns 150 entries per page. So other entries are in next pages. Following script takes all the data from the following pages (for detail, see firebrowseR vignette):

all.Received = F                         ### this boolean variable is set initially FALSE and it enters the loop till it is False, and after we receive all data, it will be set to TRUE and it will exit from the while loop.
page.Counter = 1
page.size = 150                          ### number of maximum entry within each page. 
coadread.Pats = list()                      ### Patients data from each page are stored in the list.

while(all.Received == F){
  coadread.Pats[[page.Counter]] = Samples.Clinical(format = "csv",
                                               cohort = cancer.Type,
                                               page_size = page.size,
                                               page = page.Counter)
  if(page.Counter > 1)
    colnames(coadread.Pats[[page.Counter]]) = colnames(coadread.Pats[[page.Counter-1]])  ### In the API, except first page,  pages do not have column name, so it adds the column name of from the data of first page stored as the first element of the list. 
  if(nrow(coadread.Pats[[page.Counter]]) < page.size){    ### if the last page whose entries are received has less than 150 entries, it means that itis the last page. thus, all.Received boolean variable is set to TRUE and exits the loop.           
    all.Received = T
  } else{      ### if not true, it incerements the page.counter variable and continues with the retrieval of the patient data in the next page. 
    page.Counter = page.Counter + 1
  }
}

coadread.Pats <- do.call(what = rbind , args = coadread.Pats) ### do.call () function takes the function in "what" parameter and applies it to the list object assigned to the "args" parameter.

dim(coadread.Pats)

dim(coadread.Pats)
# [1] 629  96

#  table(coadread.Pats$primary_tumor_pathologic_spread)
# 
# t3 
# 1 
#  table(coadread.Pats$primary_lymph_node_presentation_assessment)
# 
# no yes 
# 19 593 
#  table(coadread.Pats$braf_gene_analysis_result)
# 
# abnormal   normal 
# 3       32 
#  table(coadread.Pats$anatomic_site_colorectal)
# 
# rectum 
# 1 
#  table(coadread.Pats$bcr_patient_canonical_status)      # can be used.
# 
# canonical canonical - plus                               
# 512              117 
#  table(coadread.Pats$perineural_invasion_present)  
# 
# no yes 
# 173  60 
# 
#  table(coadread.Pats$person_neoplasm_cancer_status)  can be used.
# 
# tumor free with tumor 
# 388        202 
#  table(coadread.Pats$colon_polyps_present)           ### 
# 
# no yes 
# 210  92 
#  table(coadread.Pats$followup_treatment_success)
# 
# complete remission/response  partial remission/response         progressive disease 
# 252                          16                          55 
# stable disease 
# 7 
#  table(coadread.Pats$histological_type)    ### can be used first.
# 
# colon adenocarcinoma  colon mucinous adenocarcinoma          rectal adenocarcinoma 
# 391                             62                            152 
# rectal mucinous adenocarcinoma 
# 13 
#  table(coadread.Pats$kras_mutation_found)
# 
# no yes 
# 32  31 
#  table(coadread.Pats$lymphatic_invasion)  ### can be used for prediction.
# 
# no yes 
# 335 232 
#  table(coadread.Pats$lymphnode_pathologic_spread)
# 
# n0 
# 1 
#  table(coadread.Pats$non_nodal_tumor_deposits)
# 
# no yes 
# 247  50 
#  table(coadread.Pats$number_of_lymphnodes_positive_by_he)
# 
# 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  17  18  20  21  23  24  27  30  31  32  50 
# 337  61  43  36  22  23  12  14   4   4   7   5   7   3   1   2   1   2   2   1   1   1   1   1   1   1   1 
#  table(coadread.Pats$number_of_lymphnodes_positive_by_ihc)
# 
# 0  1  2  9 10 12 
# 61  1  1  1  1  1 
#  table(coadread.Pats$pathologic_stage)
# 
# stage i   stage ia   stage ii  stage iia  stage iib  stage iic  stage iii stage iiia stage iiib stage iiic 
# 108          1         38        177         12          2         25         15         86         55 
# stage iv  stage iva  stage ivb 
# 64         24          2 
#  table(coadread.Pats$pathologic_t)
# 
# t1  t2  t3  t4 t4a t4b tis 
# 20 109 427  33  27  10   1 
#  table(coadread.Pats$perineural_invasion_present)
# 
# no yes 
# 173  60 
#  table(coadread.Pats$person_neoplasm_cancer_status)   ### can be used.
# 
# tumor free with tumor 
# 388        202 
#  table(coadread.Pats$postoperative_rx_tx)  ### can be used.
# 
# no yes 
# 317 212 
#  table(coadread.Pats$primary_lymph_node_presentation_assessment)
# 
# no yes 
# 19 593 
#  table(coadread.Pats$primary_therapy_outcome_success)
# 
# complete remission/response  partial remission/response         progressive disease 
# 248                          15                          33 
# stable disease 
# 5 
#  table(coadread.Pats$primary_tumor_pathologic_spread)
# 
# t3 
# 1 
#  table(coadread.Pats$residual_tumor)
# 
# r0  r1  r2  rx 
# 458   6  38  30 
#  table(coadread.Pats$synchronous_colon_cancer_present)
# 
# no yes 
# 535  28 
#  table(coadread.Pats$tissue_source_site)
# 
# 3l  4n  4t  5m  a6  aa  ad  af  ag  ah  am  au  ay  az  bm  ca  ci  ck  cl  cm  d5  dc  dm  dt  dy  ef  ei 
# 1   1   1   3  51 174  13  18  80   7   2   2  10  21   1  10   6  14   3  37  31  13  25   1   7   2  17 
# f4  f5  g4  g5  nh  qg  ql  ru  ss  t9  ws 
# 16  12  27   4   9   5   1   1   1   1   1 
#  table(coadread.Pats$tumor_stage)
# 
# stage iia 
# 1 
#  table(coadread.Pats$tumor_tissue_site)       ### can be used.
# 
# colon rectum 
# 457    168 
#  table(coadread.Pats$venous_invasion)           ### can be used.
# 
# no yes 
# 411 134 

### the following features will be predicted:

# lymphatic_invasion
# histological_type
# person_neoplasm_cancer_status
# bcr_patient_canonical_status

### we will create a vector with the patients barcode:

sample_ids <- coadread.Pats$tcga_participant_barcode

