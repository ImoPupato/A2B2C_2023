## Survival analisys  
### Libraries  
```R
library(survival)  
library(edgeR) # for the cpm normalization  
library(tidyverse) # for data cleaning and plotting  
library("TCGAbiolinks") # to acces to the project  
library("RTCGA")  
library("RTCGA.clinical") # to download de clinical relevant data, survival status and times  
llibrary("SummarizedExperiment")  
```
### Normalization and data selection
```R
VAVs_SKCM.surv<-survivalTCGA(SKCM.clinical)
colnames(VAVs_SKCM.surv)<-c("times","ID","status")
CPM<-cpm(CPM) 
CPM<-cbind(rnaseq$SYMBOL,CPM)
phenotype<-t(rbind(CPM[CPM$SYMBOL=="VAV1",],CPM[CPM$SYMBOL=="VAV2",],CPM[CPM$SYMBOL=="VAV3",]))
phenotype<-cbind(id3,phenotype)
colnames(phenotype)<-c("ID","VAV1","VAV2,"VAV3)
VAVs_SKCM.surv_rnaseq<- VAVs_SKCM.surv %>% left_join(CPM,by = "ID")
VAVs_SKCM.surv_rnaseq.cut<-surv_cutpoint(
  VAVs_SKCM.surv_rnaseq,
  time = "times",
  event = "status",
  variables = c("VAV1","VAV2", "VAV3")
)
summary(VAVs_SKCM.surv_rnaseq.cut) # to explore the cut points
```
### Fitting and ploting the curve
```R
fit <- survfit(Surv(times, status) ~ Vav*, #the * should be changed for each vav number
               data = VAVs_SKCM.surv_rnaseq.cat)
ggsurvplot(
   fit,                     
   risk.table = TRUE,       
   pval = TRUE,             
   conf.int = FALSE,                  
   xlim = c(0,5500),      
   break.time.by = 1000,    
   ggtheme = theme_bw(), 
   risk.table.y.text.col = T, 
  font.legend=8, 
  legend.title = "Expression",
  risk.table.y.text = FALSE,                          
)
```
