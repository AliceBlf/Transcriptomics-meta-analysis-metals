library(GEOquery)
library(FactoMineR)
library(preprocessCore)
library(sva)
library(pheatmap)
library(limma)
library(fifer)

#Download GEO dataset

{
data_GSEXXXXX=getGEO("GSEXXXXX")   
  ##Replace 'XXXXX' by your GEO dataset number in the whole script
data_GSEXXXXX=data_GSEXXXXX$`GSEXXXXX_series_matrix.txt.gz`
expression_GSEXXXXX=exprs(data_GSEXXXXX)
annotation_GSEXXXXX=pData(data_GSEXXXXX)

chip_annotation_GSEXXXXX=data_GSEXXXXX@featureData@data
rownames(chip_annotation_GSEXXXXX)=chip_annotation_GSEXXXXX$ID

chip_annotation_GSEXXXXX=chip_annotation_GSEXXXXX[rownames(expression_GSEXXXXX),]
chip_annotation_GSEXXXXX=chip_annotation_GSEXXXXX[chip_annotation_GSEXXXXX$'Gene Symbol'!="",] 
  ##The typography of 'Gene Symbol' changes depending on the array

expression_GSEXXXXX=expression_GSEXXXXX[rownames(chip_annotation_GSEXXXXX),]
rownames(expression_GSEXXXXX)=chip_annotation_GSEXXXXX$'Gene Symbol'
}

#Treat probes that measure the expression of multiple gene (optional, depending on the array)

{
multiprobes = grepl("///", rownames(expression_GSEXXXXX))
  ##Replace '///' by the character that separate the gene symbol in multigene probes
multiprobes=unique(unlist(strsplit(rownames(expression_GSEXXXXX[multiprobes,]), split = " /// ", fixed = T)))
multiprobes=make.names(multiprobes)

single_probes=c()
for (k in 1:length(multiprobes)){
  b=expression_GSEXXXXX[grepl(multiprobes[k],rownames(expression_GSEXXXXX)),]
  if (length(b)>16) {b=colMeans(b)}
  single_probes=rbind(single_probes,b)}
rownames(single_probes)=multiprobes

expression_GSEXXXXX=rbind(expression_GSEXXXXX[!grepl("///", rownames(expression_GSEXXXXX)),],single_probes)
}

#Treat genes that are described by multiple probes (optional, depending on the array)

{
expression_GSEXXXXX = aggregate(expression_GSEXXXXX,by = list(rownames(expression_GSEXXXXX)),FUN=mean)
rownames(expression_GSEXXXXX)  = expression_GSEXXXXX$Group.1
expression_GSEXXXXX = expression_GSEXXXXX[,-1]
}

#Check gene expression distribution

{
boxplot(expression_GSEXXXXX,outline=F,col=rainbow(ncol(expression_GSEXXXXX)))

  ##QC: data should be in log fold change
  ##Optional, if data distribution is not gaussian-like
expression_GSEXXXXX=log2(expression_GSEXXXXX)

  ##QC: data should have similar distribution across conditions
  ##Optional, if data distribution vary from one sample to another
expression_norm_GSEXXXXX=as.data.frame(normalize.quantiles(as.matrix(expression_GSEXXXXX)))
colnames(expression_norm_GSEXXXXX)=colnames(expression_GSEXXXXX)
rownames(expression_norm_GSEXXXXX)=rownames(expression_GSEXXXXX)
expression_GSEXXXXX=expression_norm_GSEXXXXX
}

#Principal component analysis

{
  ##Selection of the genes with the higher variance to mean ratio
genes_by_VMR_GSEXXXXX=apply(expression_GSEXXXXX,MARGIN = 1,FUN = function(x) {var(x)/mean(x)})
hist(log10(genes_by_VMR_GSEXXXXX),n=100)
top_500_GSEXXXXX=genes_by_VMR_GSEXXXXX[order(genes_by_VMR_GSEXXXXX,decreasing = T)]
top_500_GSEXXXXX=names(top_500_GSEXXXXX[1:500])

PCA_GSEXXXXX=PCA(X = t(expression_GSEXXXXX[top_500_GSEXXXXX,]),scale.unit = T,graph = F)
plot.PCA(PCA_GSEXXXXX,choix = "ind",label = "none",col.ind = string.to.colors(annotation_GSEXXXXX$characteristics_ch1))
  ##the column used for color annotation has to be changed according to the considered dataset

  ##QC : identification of possible batch effect with PCA
  ##Optional, to correct batch effects
expression_GSEXXXXX = ComBat(as.matrix(expression_GSEXXXXX),batch = annotation_GSEXXXXX$characteristics_ch1)
}

#Creation of the condition matrix

{
Condition_GSEXXXXX=annotation_GSEXXXXX$title
Condition_GSEXXXXX=make.names(Condition_GSEXXXXX)
  ##Verify that the strings are indentical for replicates of a same condition and different between conditions
 
Condition_factor_GSEXXXXX=factor(Condition_GSEXXXXX, levels = unique(Condition_GSEXXXXX))
design_GSEXXXXX=model.matrix(~0+Condition_factor_GSEXXXXX)
colnames(design_GSEXXXXX)=levels(Condition_factor_GSEXXXXX)
}

#Comparison of the conditions of interest and quality control on data distribution

{
contrast_matrix_GSEXXXXX_A=makeContrasts("Condition_A-Control_Condition",levels=design_GSEXXXXX)
contrast_matrix_GSEXXXXX_B=makeContrasts("Condition_B-Control_Condition",levels=design_GSEXXXXX)
  ##Replace the condition names by the strings of the Condition_GSEXXXXX vector

fit_GSEXXXXX=lmFit(expression_GSEXXXXX,design_GSEXXXXX)

fit_GSEXXXXX_A=contrasts.fit(fit_GSEXXXXX, contrast_matrix_GSEXXXXX_A)
fit_GSEXXXXX_A=eBayes(fit_GSEXXXXX_A)
volcanoplot(fit_GSEXXXXX_A,coef=1)
  ##QC : Verify the p value and expression distribution

fit_GSEXXXXX_B=contrasts.fit(fit_GSEXXXXX, contrast_matrix_GSEXXXXX_B)
fit_GSEXXXXX_B=eBayes(fit_GSEXXXXX_B)
volcanoplot(fit_GSEXXXXX_B,coef=1)
  ##QC : Verify the p value and expression distribution
}

#Correlation between similar conditions

{
cor(fit_GSEXXXXX_A$coefficients,fit_GSEXXXXX_B$coefficients)

##To merge data if the condition A and B are sufficiently correlated (optional)
fit_GSEXXXXX_coeff=rowMeans(data.frame(fit_GSEXXXXX_24h$coefficients,fit_GSEXXXXX_30h$coefficients,fit_GSEXXXXX_36h$coefficients))
}

#Export results

write.table(fit_GSEXXXXX_coeff, file="GSEXXXXX.txt",quote = FALSE,sep = '\t')
