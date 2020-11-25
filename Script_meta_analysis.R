library(FactoMineR)
library(preprocessCore)
library(pheatmap)
library(wCorr)
library(MLmetrics)

#General table contruction

{

list_files = list.files("Path/to/directory/Gene_expression_tables")

unique_gene=c()
for (i in 1:length(list_files)){
  a=rownames(read.table(list_files[i], header=T, sep="\t"))
  unique_gene=unique(c(unique_gene,a))
}


all_data=matrix(data=NA, nrow=length(unique_gene), ncol=length(list_files))
for (i in 1:length(list_files)){
  b=read.table(list_files[i], header=T, sep="\t")
  for (k in 1:length(unique_gene)){
    all_data[k,i]=as.numeric(ifelse(unique_gene[k]%in%rownames(b), b[as.character(unique_gene[k]),], NA))
  }
}

rownames(all_data)=unique_gene
colnames(all_data)=gsub(".txt","",list_files)


Dataset_info=read.table("Path/to/directory/Dataset_information.txt", header=T, sep="\t")

}

#Gene and dataset filtering

{
NA_per_line = rowSums(is.na(all_data))
names(NA_per_line)=rownames(all_data)
genes_selected = names(NA_per_line[NA_per_line<length(list_files)/3])
genes_selected=genes_selected[genes_selected!=""]

all_data_select=all_data[genes_selected,]

rownames(all_data_select)=genes_selected
colnames(all_data_select)=as.character(Dataset_info[,1])

NA_per_column = colSums(is.na(all_data_select))
Dataset_info=Dataset_info[NA_per_column<(dim(all_data_select)[1]/3),]
Dataset_info=Dataset_info[Dataset_info$metal_type!="cobalt"&Dataset_info$metal_type!="nickel",]
all_data_select = all_data_select[,as.character(Dataset_info[,1])]
}

#Normalization

{

all_data_select_norm = normalize.quantiles(all_data_select)
rownames(all_data_select_norm)=rownames(all_data_select)
colnames(all_data_select_norm)=colnames(all_data_select)
}

#PCA

{
Var_logFC = apply(all_data_select_norm,MARGIN = 1,FUN = function(x) {var(x,na.rm=T)})
hist(sqrt(Var_logFC),100)
abline(v=0.38,col="red",lty=2,lwd=2)
Variable_genes = names(which(sqrt(Var_logFC)>0.38))

PCA_all_data = PCA(t(all_data_select_norm[Variable_genes,]),scale.unit = T, graph = F)
barplot(PCA_all_data$eig[,2])

plot.PCA(PCA_all_data,choix = "ind", col.ind = rainbow(ncol(all_data_select)), cex=0.6)
plot.PCA(PCA_all_data,choix = "ind", col.ind = rainbow(ncol(all_data_select)), axes=c(1,3), cex=0.6)

}

#Individual signatures

{
score_Pt = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="platinum",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Pt = rowMeans(score_Pt, na.rm = TRUE)
score_Pt_ordered = score_Pt[order(as.numeric(score_Pt))]
signature_platine = names(score_Pt_ordered[1:50])

score_Ti = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="titanium",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Ti = rowMeans(score_Ti, na.rm = TRUE)
score_Ti_ordered = score_Ti[order(as.numeric(score_Ti))]
signature_titane = names(score_Ti_ordered[1:50])

score_Fe_ion = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="iron"&Dataset_info$form=="ion",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Fe_ion = rowMeans(score_Fe_ion, na.rm = TRUE)
score_Fe_ion_ordered = score_Fe_ion[order(as.numeric(score_Fe_ion))]
signature_fer_ion = names(score_Fe_ion_ordered[1:50])

score_Fe_np = apply(as.data.frame(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="iron"&Dataset_info$form=="np",1])]),MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Fe_np = rowMeans(score_Fe_np, na.rm = TRUE)
score_Fe_np_ordered = score_Fe_np[order(as.numeric(score_Fe_np))]
signature_fer_np = names(score_Fe_np_ordered[1:50])

score_Au_ion = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="gold"&Dataset_info$form=="ion",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Au_ion = rowMeans(score_Au_ion, na.rm = TRUE)
score_Au_ion_ordered = score_Au_ion[order(as.numeric(score_Au_ion))]
signature_or_ion = names(score_Au_ion_ordered[1:50])

score_Au_np = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="gold"&Dataset_info$form=="np",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Au_np = rowMeans(score_Au_np, na.rm = TRUE)
score_Au_np_ordered = score_Au_np[order(as.numeric(score_Au_np))]
signature_or_np = names(score_Au_np_ordered[1:50])

score_Ag_ion = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="silver"&Dataset_info$form=="ion",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Ag_ion = rowMeans(score_Ag_ion, na.rm = TRUE)
score_Ag_ion_ordered = score_Ag_ion[order(as.numeric(score_Ag_ion))]
signature_argent_ion = names(score_Ag_ion_ordered[1:50])

score_Ag_np = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="silver"&Dataset_info$form=="np",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Ag_np = rowMeans(score_Ag_np, na.rm = TRUE)
score_Ag_np_ordered = score_Ag_np[order(as.numeric(score_Ag_np))]
signature_argent_np = names(score_Ag_np_ordered[1:50])

score_Cu_ion = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="copper"&Dataset_info$form=="ion",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Cu_ion = rowMeans(score_Cu_ion, na.rm = TRUE)
score_Cu_ion_ordered = score_Cu_ion[order(as.numeric(score_Cu_ion))]
signature_cuivre_ion = names(score_Cu_ion_ordered[1:50])

score_Cu_np = apply(as.data.frame(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="copper"&Dataset_info$form=="np",1])]),MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Cu_np = rowMeans(score_Cu_np, na.rm = TRUE)
score_Cu_np_ordered = score_Cu_np[order(as.numeric(score_Cu_np))]
signature_cuivre_np = names(score_Cu_np_ordered[1:50])


score_Zn_ion = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="zinc"&Dataset_info$form=="ion",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Zn_ion = rowMeans(score_Zn_ion, na.rm = TRUE)
score_Zn_ion_ordered = score_Zn_ion[order(as.numeric(score_Zn_ion))]
signature_zinc_ion = names(score_Zn_ion_ordered[1:50])

score_Zn_np = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="zinc"&Dataset_info$form=="np",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Zn_np = rowMeans(score_Zn_np, na.rm = TRUE)
score_Zn_np_ordered = score_Zn_np[order(as.numeric(score_Zn_np))]
signature_zinc_np = names(score_Zn_np_ordered[1:50])

score_Cd = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="cadmium",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Cd = rowMeans(score_Cd, na.rm = TRUE)
score_Cd_ordered = score_Cd[order(as.numeric(score_Cd))]
signature_cadmium = names(score_Cd_ordered[1:50])

score_C = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="carbon",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_C = rowMeans(score_C, na.rm = TRUE)
score_C_ordered = score_C[order(as.numeric(score_C))]
signature_carbone = names(score_C_ordered[1:50])

score_Si = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="silica",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Si = rowMeans(score_Si, na.rm = TRUE)
score_Si_ordered = score_Si[order(as.numeric(score_Si))]
signature_silice = names(score_Si_ordered[1:50])
}

#Normalization scores

{
  all_scores=cbind(score_Pt,score_Ti,score_Cd,score_Au_np,score_Au_ion,score_Ag_np,score_Ag_ion,
                   score_Cu_np,score_Cu_ion,score_Fe_np,score_Fe_ion,score_Zn_np,score_Zn_ion,
                   score_C, score_Si)
                   
  all_scores_norm=normalize.quantiles(all_scores)
  rownames(all_scores_norm)=rownames(all_scores)
  colnames(all_scores_norm)=colnames(all_scores)
}


#Weighted correlation between scores

{
  
correlation_matrix=matrix(NA, nrow=15, ncol=15)
pvalue_matrix=matrix(NA, nrow=15, ncol=15)
    
for (i in 1:15){
  for (j in 1:15){
    correlation_matrix[i,j]=weightedCorr(all_scores_norm[,i],all_scores_norm[,j], 
                                           weights = 1/((all_scores_norm[,i]^2)*(all_scores_norm[,j]^2)), method="Pearson")
  }
}

for (i in 1:15){
  for (j in 1:15){
    pvalue_matrix[i,j]=as.numeric(summary(lm(all_scores_norm[,i]~all_scores_norm[,j], 
                                             weight=1/(all_scores_norm[,i]^2*all_scores_norm[,j]^2)))[[5]][,4][2])
  }
}
    
correlation_matrix=as.data.frame(correlation_matrix)
pvalue_matrix=as.data.frame(pvalue_matrix)
    
colnames(correlation_matrix)=c("Platinum (ion)","Titanium (NP)", "Cadmium (ion)","Gold (NP)", "Gold (ion)",
                                 "Silver (NP)", "Silver (ion)","Copper (NP)", "Copper (ion)","Iron (NP)", 
                                 "Iron (ion)","Zinc (NP)", "Zinc (ion)", "Carbon (NP)", "Silicon (NP)")
rownames(correlation_matrix)=colnames(correlation_matrix)
colnames(pvalue_matrix)=colnames(correlation_matrix)
rownames(correlation_matrix)=colnames(correlation_matrix)
    
pheatmap(correlation_matrix, clustering_method = "ward.D", display_numbers=T)

  
}


#Response to Zn, Au, Ag, Cu and Cd

{
  
score_Zn = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="zinc",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Zn = rowMeans(score_Zn, na.rm = TRUE)
score_Au = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="gold",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Au = rowMeans(score_Au, na.rm = TRUE)
score_Ag = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="silver",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Ag = rowMeans(score_Ag, na.rm = TRUE)
score_Cu = apply(all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal_type=="copper",1])],MARGIN = 2,FUN = function(x) {rank(-x) } )
score_Cu = rowMeans(score_Cu, na.rm = TRUE)
  
Score_ZnAgAuCuCd= cbind(score_Zn,score_Ag,score_Au,score_Cu,score_Cd)
Score_ZnAgAuCuCd=normalize.quantiles(Score_ZnAgAuCuCd)
colnames(Score_ZnAgAuCuCd)=c("Zinc", "Silver", "Gold", "Copper", "Cadmium")
rownames(Score_ZnAgAuCuCd)=names(score_Cd)

Score_ZnAgAuCuCd=rowMeans(Score_ZnAgAuCuCd)
Score_ZnAgAuCuCd_ordered=Score_ZnAgAuCuCd[order(Score_ZnAgAuCuCd)]


  ##Gene selection

  {
  
  #Initialisation
  
  initial_genes=c(1:6)
  
  ideal_clust=ifelse(Dataset_info[,5]=="cadmium"|Dataset_info[,5]=="zinc"|Dataset_info[,5]=="copper"|Dataset_info[,5]=="gold"|Dataset_info[,5]=="silver",2,1)
  hierarc_clust_init=cutree(hclust(dist(t(all_data_select_norm[names(Score_ZnAgAuCuCd_ordered[initial_genes]),])), method="ward.D2"),k=2)
  F1_score_init=F1_Score(ideal_clust, hierarc_clust_init, positive="2")
  
  
  #Step forward
  
  F1_max=F1_score_init
  selected_genes_ZnAuAgCuCd=initial_genes
  
  for (i in 7:50)
    
  {
    hierarc_clust_boucle=cutree(hclust(dist(t(all_data_select_norm[names(Score_ZnAgAuCuCd_ordered[c(selected_genes_ZnAuAgCuCd,i)]),])), method="ward.D2"),k=2)
    F1_score_boucle=F1_Score(ideal_clust, hierarc_clust_boucle, positive="2")
    if (is.na(F1_score_boucle)) {F1_score_boucle=0}
    
    if (F1_score_boucle>F1_max) {selected_genes_ZnAuAgCuCd=c(selected_genes_ZnAuAgCuCd,i)
    F1_max=F1_score_boucle}
    
  }
  
  
  #Step backward
  
  final_genes_ZnAgAuCuCd=c()
  
  for (j in 1:length(selected_genes_ZnAuAgCuCd)){
    
    hierarc_clust_boucle_2=cutree(hclust(dist(t(all_data_select_norm[names(Score_ZnAgAuCuCd_ordered[selected_genes_ZnAuAgCuCd[-j]]),])), method="ward.D2"),k=2)
    F1_score_boucle_2=F1_Score(ideal_clust, hierarc_clust_boucle_2, positive="2")
    if (is.na(F1_score_boucle_2)) {F1_score_boucle_2=0}
    
    if (F1_score_boucle_2<=F1_max) {final_genes_ZnAgAuCuCd=c(final_genes_ZnAgAuCuCd, selected_genes_ZnAuAgCuCd[j])} }
  

  F1_max
  names(Score_ZnAgAuCuCd_ordered[c(final_genes_ZnAgAuCuCd)])
  
  
  
  }

data_heatmap=all_data_select_norm[c(names(Score_ZnAgAuCuCd_ordered[c(final_genes_ZnAgAuCuCd)])),as.character(Dataset_info[,1])]

annotation_heatmap=as.data.frame(Dataset_info[,c(5,6)])
rownames(annotation_heatmap)=colnames(data_heatmap)
colnames(annotation_heatmap)=c("Metal","Form")

pheatmap(data_heatmap, cluster_rows = F, cluster_cols = T, clustering_method = "ward.D2",show_colnames=F ,
         cutree_cols = 2, annotation_col = annotation_heatmap)
}


#Response to NPs

{
  score_NP=cbind(score_Ti, score_C, score_Au_np, score_Fe_np, score_Zn_np,score_Ag_np,score_Cu_np)
  
  score_NP=normalize.quantiles(score_NP)
  colnames(score_NP)=c("Titane", "Carbon", "Or", "Argent", "Cuivre","Fer","Zinc")
  rownames(score_NP)=names(score_Ti)
  
  
  score_NP_all=rowMeans(score_NP)
  score_NP_all_ordered=score_NP_all[order(score_NP_all)]
  score_NP_persist=rowMeans(score_NP[,c(1:3)])
  score_NP_persist_ordered=score_NP_persist[order(score_NP_persist)]
  
  
  ##Gene selection
  
  {
    
    #Initialisation
    
    initial_genes=c(1:6)
    
    ideal_clust=ifelse(Dataset_info[,6]=="np",2,1)
    hierarc_clust_init=cutree(hclust(dist(t(all_data_select_norm[names(score_NP_persist_ordered[initial_genes]),])), method="ward.D2"),k=2)
    F1_score_init=F1_Score(ideal_clust, hierarc_clust_init, positive="2")

    
    #Step forward
    
    F1_max=F1_score_init
    selected_genes_NP=initial_genes
    
    for (l in 7:50)
      
    {
      hierarc_clust_boucle=cutree(hclust(dist(t(all_data_select_norm[names(score_NP_persist_ordered[c(selected_genes_NP,l)]),])), method="ward.D2"),k=2)
      F1_score_boucle=F1_Score(ideal_clust, hierarc_clust_boucle, positive="2")
      if (is.na(F1_score_boucle)) {F1_score_boucle=0}
      
      if (F1_score_boucle>F1_max) {selected_genes_NP=c(selected_genes_NP,l)
      F1_max=F1_score_boucle}
      
    }
    
    
    #Step backward
    
    final_genes_NPs=c()
    
    for (j in 1:length(selected_genes_NP)){
      
      hierarc_clust_boucle_2=cutree(hclust(dist(t(all_data_select_norm[names(score_NP_persist_ordered[selected_genes_NP[-j]]),])), method="ward.D2"),k=2)
      F1_score_boucle_2=F1_Score(ideal_clust, hierarc_clust_boucle_2, positive="2")
      if (is.na(F1_score_boucle_2)) {F1_score_boucle_2=0}
      
      if (F1_score_boucle_2<=F1_max) {final_genes_NPs=c(final_genes_NPs, selected_genes_NP[j])} }
    

    F1_max
    names(score_NP_persist_ordered[c(final_genes_NPs)])
    
  }
  
  data_heatmap_NP=all_data_select_norm[names(score_NP_persist_ordered[c(final_genes_NPs)]),as.character(Dataset_info[,1])]
  
  pheatmap(data_heatmap_NP, cluster_rows = F, cluster_cols = T, clustering_method = "ward.D2", border_color="grey50",show_colnames=F ,
           cutree_cols = 2, annotation_col = annotation_heatmap)
  
}

#Response to essential metals

{
    score_essential=all_scores_norm[,c(9,11,13)]
    score_essential=rowMeans(score_essential)
    score_essential_ordered=score_essential[order(as.numeric(score_essential))]
    
    data_select_norm_metals=all_data_select_norm[,as.character(Dataset_info[Dataset_info$metal=="yes",1])]
    
    
    ##Gene selection
    
    
    {
      
      #Initialisation
      
      initial_genes=c(1:6)
      
      ideal_clust=ifelse(Dataset_info[Dataset_info$metal=="yes",3]=="essential",2,1)
      hierarc_clust_init=cutree(hclust(dist(t(data_select_norm_metals[names(score_essential_ordered[initial_genes]),])), method="ward.D2"),k=2)
      F1_score_init=F1_Score(ideal_clust, hierarc_clust_init, positive="2")
      
      
      #Step forward
      
      F1_max=F1_score_init
      selected_genes_essential=initial_genes
      
      for (l in 7:50)
        
      {
        hierarc_clust_boucle=cutree(hclust(dist(t(data_select_norm_metals[names(score_essential_ordered[c(selected_genes_essential,l)]),])), method="ward.D2"),k=2)
        F1_score_boucle=F1_Score(ideal_clust, hierarc_clust_boucle, positive="2")
        if (is.na(F1_score_boucle)) {F1_score_boucle=0}
        
        if (F1_score_boucle>F1_max) {selected_genes_essential=c(selected_genes_essential,l)
        F1_max=F1_score_boucle}
        
      }
      
      
      #Step backward
      
      final_genes_essential=c()
      
      for (j in 1:length(selected_genes_essential)){
        
        hierarc_clust_boucle_2=cutree(hclust(dist(t(data_select_norm_metals[names(score_essential_ordered[selected_genes_essential[-j]]),])), method="ward.D2"),k=2)
        F1_score_boucle_2=F1_Score(ideal_clust, hierarc_clust_boucle_2, positive="2")
        if (is.na(F1_score_boucle_2)) {F1_score_boucle_2=0}
        
        if (F1_score_boucle_2<=F1_max) {final_genes_essential=c(final_genes_essential, selected_genes_essential[j])} }
      
      
      F1_max
      names(score_essential_ordered[c(final_genes_essential)])
      
    }
    
    
    data_heatmap_essentiel=data_select_norm_metals[names(score_essential_ordered[c(final_genes_essential)]),]
    
    pheatmap(data_heatmap_essentiel, cluster_rows = F, cluster_cols = T, clustering_method = "ward.D2", border_color="grey50",show_colnames=F ,
             cutree_cols = 2, annotation_col = annotation_heatmap)
    
  }
  
#Response to highly toxic metals
  
{
  score_toxic=all_scores_norm[,c(1,3,7)]
    score_toxic=rowMeans(score_toxic)
    score_toxic_ordered=score_toxic[order(as.numeric(score_toxic))]
    
    
    ##Gene selection
    
    
    {
      
      #Initialisation
      
      initial_genes=c(1:6)
      
      ideal_clust=ifelse(Dataset_info[Dataset_info$metal=="yes",4]=="toxic",2,1)
      hierarc_clust_init=cutree(hclust(dist(t(data_select_norm_metals[names(score_toxic_ordered[initial_genes]),])), method="ward.D2"),k=2)
      F1_score_init=F1_Score(ideal_clust, hierarc_clust_init, positive="2")
      
      
      #Step forward
      
      F1_max=F1_score_init
      selected_genes_toxicity=initial_genes
      
      for (l in 7:50)
        
      {
        hierarc_clust_boucle=cutree(hclust(dist(t(data_select_norm_metals[names(score_toxic_ordered[c(selected_genes_toxicity,l)]),])), method="ward.D2"),k=2)
        F1_score_boucle=F1_Score(ideal_clust, hierarc_clust_boucle, positive="2")
        if (is.na(F1_score_boucle)) {F1_score_boucle=0}
        
        if (F1_score_boucle>F1_max) {selected_genes_toxicity=c(selected_genes_toxicity,l)
        F1_max=F1_score_boucle}
        
      }
      
      
      #Step backward
      
      final_genes_toxicity=c()
      
      for (j in 1:length(selected_genes_toxicity)){
        
        hierarc_clust_boucle_2=cutree(hclust(dist(t(data_select_norm_metals[names(score_toxic_ordered[selected_genes_toxicity[-j]]),])), method="ward.D2"),k=2)
        F1_score_boucle_2=F1_Score(ideal_clust, hierarc_clust_boucle_2, positive="2")
        if (is.na(F1_score_boucle_2)) {F1_score_boucle_2=0}
        
        if (F1_score_boucle_2<=F1_max) {final_genes_toxicity=c(final_genes_toxicity, selected_genes_toxicity[j])} }
      
      
      F1_max
      names(score_toxic_ordered[c(final_genes_toxicity)])
      
    }
    
    
    data_heatmap_toxicity=data_select_norm_metals[names(score_toxic_ordered[c(final_genes_toxicity)]),]
    
    pheatmap(data_heatmap_toxicity, cluster_rows = F, cluster_cols = T, clustering_method = "ward.D2", border_color="grey50",show_colnames=F ,
             cutree_cols = 2, annotation_col = annotation_heatmap)
    
    
  }
  