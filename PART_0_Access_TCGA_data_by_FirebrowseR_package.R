
### I will use FirebrowseR package to download an analyzed data from Firebrowse API.

library(FirebrowseR)

cohorts <- Metadata.Cohorts(format = "csv") ### Download all cohorts that are available.


print(cohorts)

# cohort                                                      description
# 1       ACC                                         Adrenocortical carcinoma
# 2      BLCA                                     Bladder Urothelial Carcinoma
# 3      BRCA                                        Breast invasive carcinoma
# 4      CESC Cervical squamous cell carcinoma and endocervical adenocarcinoma
# 5      CHOL                                               Cholangiocarcinoma
# 6      COAD                                             Colon adenocarcinoma
# 7  COADREAD                                        Colorectal adenocarcinoma
# 8      DLBC                  Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
# 9      ESCA                                            Esophageal carcinoma 
# 10     FPPP                                              FFPE Pilot Phase II
# 11      GBM                                          Glioblastoma multiforme
# 12   GBMLGG                                                           Glioma
# 13     HNSC                            Head and Neck squamous cell carcinoma
# 14     KICH                                               Kidney Chromophobe
# 15    KIPAN                               Pan-kidney cohort (KICH+KIRC+KIRP)
# 16     KIRC                                Kidney renal clear cell carcinoma
# 17     KIRP                            Kidney renal papillary cell carcinoma
# 18     LAML                                           Acute Myeloid Leukemia
# 19      LGG                                         Brain Lower Grade Glioma
# 20     LIHC                                   Liver hepatocellular carcinoma
# 21     LUAD                                              Lung adenocarcinoma
# 22     LUSC                                     Lung squamous cell carcinoma
# 23     MESO                                                     Mesothelioma
# 24       OV                                Ovarian serous cystadenocarcinoma
# 25     PAAD                                        Pancreatic adenocarcinoma
# 26     PCPG                               Pheochromocytoma and Paraganglioma
# 27     PRAD                                          Prostate adenocarcinoma
# 28     READ                                            Rectum adenocarcinoma
# 29     SARC                                                          Sarcoma
# 30     SKCM                                          Skin Cutaneous Melanoma
# 31     STAD                                           Stomach adenocarcinoma
# 32     STES                                 Stomach and Esophageal carcinoma
# 33     TGCT                                      Testicular Germ Cell Tumors
# 34     THCA                                                Thyroid carcinoma
# 35     THYM                                                          Thymoma
# 36     UCEC                             Uterine Corpus Endometrial Carcinoma
# 37      UCS                                           Uterine Carcinosarcoma
# 38      UVM                                                   Uveal Melanoma


cancer.Type = cohorts[grep("pancreatic", cohorts$description, ignore.case = T), 1] ### find the cohort with decription containing "pancreatic":




print(cancer.Type)
# [1] "PAAD"

# Now that we know that the breast cancer samples are identified be  "PAAD", we can retrieve a list
# of all patients associated with this identifier.

paad.Pats <- Samples.Clinical(cohort = cancer.Type, format="tsv")

dim(paad.Pats)

###[1] 150  93 ### there are 150 patients. However, there are more than this according to (http://firebrowse.org/) or (https://gdac.broadinstitute.org/).

### The reason is that Firebrowse API returns 150 entries per page. So other entries are in next pages. Following script takes all the data from the following pages (for detail, see firebrowseR vignette):

all.Received = F                         ### this boolean variable is set initially FALSE and it enters the loop till it is False, and after we receive all data, it will be set to TRUE and it will exit from the while loop.
page.Counter = 1
page.size = 150                          ### number of maximum entry within each page. 
paad.Pats = list()                      ### Patients data from each page are stored in the list.

while(all.Received == F){
  paad.Pats[[page.Counter]] = Samples.Clinical(format = "csv",
                                               cohort = cancer.Type,
                                               page_size = page.size,
                                               page = page.Counter)
  if(page.Counter > 1)
    colnames(paad.Pats[[page.Counter]]) = colnames(paad.Pats[[page.Counter-1]])  ### In the API, except first page,  pages do not have column name, so it adds the column name of from the data of first page stored as the first element of the list. 
  if(nrow(paad.Pats[[page.Counter]]) < page.size){    ### if the last page whose entries are received has less than 150 entries, it means that itis the last page. thus, all.Received boolean variable is set to TRUE and exits the loop.           
    all.Received = T
  } else{      ### if not true, it incerements the page.counter variable and continues with the retrieval of the patient data in the next page. 
    page.Counter = page.Counter + 1
  }
}

paad.Pats <- do.call(what = rbind , args = paad.Pats) ### do.call () function takes the function in "what" parameter and applies it to the list object assigned to the "args" parameter.

dim(paad.Pats) ### we check if number of patients increased:
# [1] 185  93

# ### optional for high number of patients: "We now have collected all samples. Next we subset this data frame to patients being dead. 
# We only do this to keep the run time short, downloading mRNA expression data for a thousand patients would
# take a lot of time, later on.

### paad.Pats = paad.Pats[which(paad.Pats$vital_status == "dead"), ] ### filter the rows with dead status.

## We can do a gene-wise investigation of samples with specified gene names.

# ### "Here we define a vector containing genes known to be differential expressed in pancreatic cancer and
# download the mRNA expression data for these genes and our patients. Since there are a lots of PAAD
# samples available, we chunk this query into a gene-wise subset"

diff.Exp.Genes = c("TP53")

all.Found = F
page.Counter = 1
mRNA.Exp = list()
page.Size = 2000 # using a bigger page size is faster
while(all.Found == F){
  mRNA.Exp[[page.Counter]] = Samples.mRNASeq(format = "csv",
                                             gene = diff.Exp.Genes,
                                             cohort = "PAAD",
                                             tcga_participant_barcode =
                                               paad.Pats$tcga_participant_barcode,
                                             page_size = page.Size,
                                             page = page.Counter)
  if(nrow(mRNA.Exp[[page.Counter]]) < page.Size)
    all.Found = T
  else
    page.Counter = page.Counter + 1
}

Sample

mRNA.Exp = do.call(rbind, mRNA.Exp)

dim(mRNA.Exp)

###  [1] 183   8

# We only keep the samples having a primary tumor and corresponding normal tissue available. Normal
# tissue is encoded by `NT` and tumor tissue by `TP`. These identifiers can be decoded using the
# Metadata.SampleTypes("csv") function.
# 

Metadata.SampleTypes("csv")
# code short_letter_code                                        definition
# 1     1                TP                               Primary Solid Tumor
# 2     2                TR                             Recurrent Solid Tumor
# 3     3                TB   Primary Blood Derived Cancer - Peripheral Blood
# 4     4              TRBM      Recurrent Blood Derived Cancer - Bone Marrow
# 5     5               TAP                          Additional - New Primary
# 6     6                TM                                        Metastatic
# 7     7               TAM                             Additional Metastatic
# 8     8              THOC                        Human Tumor Original Cells
# 9     9               TBM        Primary Blood Derived Cancer - Bone Marrow
# 10   10                NB                              Blood Derived Normal
# 11   11                NT                               Solid Tissue Normal
# 12   12               NBC                                Buccal Cell Normal
# 13   13              NEBV                           EBV Immortalized Normal
# 14   14               NBM                                Bone Marrow Normal
# 15   20             CELLC                                   Control Analyte
# 16   40               TRB Recurrent Blood Derived Cancer - Peripheral Blood
# 17   50              CELL                                        Cell Lines
# 18   60                XP                          Primary Xenograft Tissue
# 19   61               XCL                Cell Line Derived Xenograft Tissue

# Patients with normal tissue
normal.Tissue.Pats = which(mRNA.Exp$sample_type == "NT")
# get the patients barcodes
patient.Barcodes = mRNA.Exp$tcga_participant_barcode[normal.Tissue.Pats]
# Subset the mRNA.Exp data frame, keeping only the pre-selected barcodes AND
# having a sample type of NT or TP
mRNA.Exp = mRNA.Exp[which(mRNA.Exp$tcga_participant_barcode %in% patient.Barcodes &
                            mRNA.Exp$sample_type %in% c("NT", "TP")), ]

dim(mRNA.Exp)

# [1] 8 8

###"Now we can use the famous ggplot2 package to plot the expression."

library(ggplot2)

p = ggplot(mRNA.Exp, aes(factor(gene), z.score))

p + geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")