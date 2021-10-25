
### R: TMM and RPKM count for host/culture analysis

# R
library("limma")
library("edgeR")

#object with read counts from htseq output
counts = readDGE(RNAseq.data$HTSeqfile, header=F)$counts

cpms = cpm(counts) #count per million in file 'counts'
#read filtering
#matrix
counts.m <- data.matrix(cpms)
y <- DGEList(counts.m) #object used by edgeR to store count data


#normalization for composition bias
# Apply normalization to DGEList object
y <- calcNormFactors(y)

#gene length to DGEList
y$genes <- data.frame(gene_length)
#apply rpkm normalization
count.rpkm <- rpkm(y)


Linear model regression for protein conservation and expression variation.
#R
#dataset highly conserved genes (n=6070).
#libraries
library(dplyr)
library(GGally)
library(ggplot2)
library(ggfortify)
#log transformation of gene expression variation (CV).
dataT <- highly_conserved %>% mutate(CVp_log=log(CVplant), CVv_log=log(CVvitro).
# Relationship between expression and conservation
#null hypothesis
model_eVSc<-lm(CVp_log ~ amino_acid_identity, data=dataT)
intercept_null= model_eVSc[[1]][[1]]
slope_null= model_eVSc[[1]][[2]]
#subset by H3K27me3 marks
dataT_M<-dataT[dataT$H3K27m3==1,] #marked
dataT_NM<-dataT[dataT$H3K27m3==0,] #unmarked
#Plot H3K27me3 marked and unmarked regression models
ggplot(data = dataT, aes(y = CVv_log, x = amino_acid_identity, color=H3K27m3, shape=H3K27m3)) +
  geom_point(size = 1, alpha = 0.2) +
  geom_smooth(method = "lm", formula = y ~ x, fullrange = TRUE) +
  geom_abline(intercept = intercept_null, slope = slope_null, colour="black")
#Multiple regression relation between protein conservation and expression variation accounting for H3K27me3 marks.
lm_cVSe <- lm(CVp_log ~ amino_acid_identity * H3K27m3, data=dataT)
summary(lm_cVSe)
