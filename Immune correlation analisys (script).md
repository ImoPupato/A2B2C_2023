### In this case we used the data downloaded from [TIMER2.0](http://timer.cistrome.org/)
```R
data<-immune.matrix[!duplicated(immune.matrix$ID),]
```
### Libraries
```R
library("corrplot")
library("ggstatsplot")
```
### Correlation analysis
```R
corrplot(cordata, tl.cex=0.1,tl.col = "black",number.cex = 0.1)
cor.test(data$VAV*,data$ESTIMATE_score,method="spearman") # the * should be change for the number of the Vav 
cor.test(data$VAV*,data$IMMUNE_score,method="spearman")
cor.test(data$VAV*,data$microenvironment.score_XCELL,method="spearman")
cor.test(data$VAV*,data$immune.score_XCELL,method="spearman")
```
