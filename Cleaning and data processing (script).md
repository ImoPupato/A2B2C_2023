### Libraries to acces and data download
```R
library("TCGAbiolinks")
library("RTCGA")
library("RTCGA.clinical")
llibrary("SummarizedExperiment")
```
### To install the packages from [RTCGA](https://rtcga.github.io/RTCGA)
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RTCGA")
```
### Libraries to data cleaning and manipulation
```R
library(limma)
library(edgeR) # for cpm normalization
library(tidyverse)
library(org.Hs.eg.db) # for gene notation
```
### - Accesing to clinical data
```R
clinica<-as.data.frame(SKCM.clinical) 
clinica$bcr_patient_barcode<-toupper(clinica$patient.bcr_patient_barcode) # modifying the ID
id<-as.data.frame(clinica$bcr_patient_barcode)
x<-x[,c(1,9,11,13,16,28,30,31,438,439,448,449,450,770,775,778,779,780,785,917,935,939,940,941,944,945,846,949,956,964,1060,1867)] # selecting usefull information
```
### - Search and download the RNAseq (raw counts)
```R
query.exp <- GDCquery(
  project = "TCGA-SKCM", # I use the SKCM project, you may use BRCA or breast cancer or for MMRF bone marrow, etc.
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts") # For raw counts

GDCdownload(
  query = query.exp,
  files.per.chunk = 100)

skcm.exp <- GDCprepare(
  query = query.exp,
  save = TRUE,
  save.filename = "skcmExp.rda")

rnaseq <- assay(skcm.exp)
rnaseq <- rnaseq[which(!duplicated(rownames(rnaseq))),]   
write.table(rnaseq, "raw counts.txt") # saving the data set
```
### - Data cleaning
#### - Id changes
```R
id2<-as.data.frame(colnames(rnaseq))
colnames(id2)<-c("bcr_patient_barcode")
id2<-id%>%mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12))
colnames(rnaseq)<-id2$bcr_patient_barcode
```
#### - Selecting patients with complete information
```R
rnaseq<-select(rnaseq, by=c(x$id))
CPM<-rnaseq[,-1] # save it for later
id3<-as.data.frame(colnmaes(rnaseq))
```
#### - Add ENTREZ and GENE SYMBOL codification
```R
hs <- org.Hs.eg.db
my.symbols <- c(rnaseq$ENSEMBL)
IDS<-select(hs, 
            keys = my.symbols,
            columns = c("ENSEMBL", "SYMBOL","ENTREZID"),
            keytype = "ENSEMBL")
IDS<-na.omit(IDS)
rnaseq<-merge(rnaseq,IDS,by="ENSEMBL")
```
#### - Characterize the patients using the cutpoints extracted from the survival analisys (use the data created for the survival analisys)
```R
phenotype$PhenoV1<-ifelse(phenotype$VAV1<7.717641,"Low","High")
phenotype$PhenoV2<-ifelse(phenotype$VAV2<57.475259,"Low","High")
phenotype$PhenoV3<-ifelse(phenotype$VAV3<9.791567,"Low","High")
```
