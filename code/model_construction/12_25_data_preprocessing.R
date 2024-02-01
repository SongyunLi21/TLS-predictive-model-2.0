library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(stringi)
library(ROCit)
library(png)
library(SeuratData)
library(patchwork)
####Data Preprocessing####
#RI model
#slide_list <- c("c_2",
#                "c_3",
#                "c_4",
#                "c_7",
#                "c_20",
#                "c_34",
#                "c_36",
#                "c_39",
#                "c_45",
#                "c_51")
#slide_list <- c("c_3",
#                "c_4",
#                "c_36")
#slide_list <- c("c_2",
#                "c_7",
#                "c_20",
#                "c_34",
#                "c_39",
#                "c_45",
#                "c_51")
#NRI model
#slide_list <- c("b_1",
#                "b_18",
#                "a_3",
#                "a_15")
slide_list <- c("b_1",
                "b_18",
                "a_3")


#slide_list <- c("a_15")
#slide_list <- c("c_8",
#                "c_5",
#                "c_23",
#                "c_57")
####inter####
library(readr)
inter <- read_csv("data/inter_NRI.csv")
inter <- inter$x
####preprocessing####
spatial_list <- sapply(slide_list,function(slide){
  print(slide)
  expr.url <- paste0("data/st/",slide,"/filtered_feature_bc_matrix.h5")
  spatial_object <- Seurat::Read10X_h5(filename =  expr.url )# Collect all genes coded on the mitochondrial genome
  spatial_object <- Seurat::CreateSeuratObject(counts = spatial_object, project = slide, assay = 'Spatial')
  mt.genes <- grep(pattern = "^MT-", x = rownames(spatial_object), value = TRUE)
  spatial_object$percent.mito <- (Matrix::colSums(spatial_object@assays$Spatial@counts[mt.genes, ])/Matrix::colSums(spatial_object@assays$Spatial@counts))*100
  #remove mt genes
  #genes_to_keep <- setdiff(names(which(Matrix::rowSums(spatial_object@assays$Spatial@counts )>5)),mt.genes)
  #spatial_object_subset <- subset(feature=genes_to_keep,spatial_object, subset = nFeature_Spatial > 300 & percent.mito < 30)
  spatial_object_subset <- subset(feature=inter,spatial_object, subset = nFeature_Spatial > 300 & percent.mito < 30)
  #cat("Spots removed: ", ncol(spatial_object) - ncol(spatial_object_subset), "\n")
  #cat("Genes kept: ", length(genes_to_keep),"from",nrow(spatial_object), "\n") 
  #spatial_object_subset <- SCTransform(spatial_object_subset, assay = "Spatial", verbose = T)
})
####intersection####
intersection <- c()
n = 0
for(i in slide_list){
  if(length(slide_list)!=3){
    break
  }
  gene <- rownames(spatial_list[[i]]@assays$Spatial@counts)
  if(n == 0){
    intersection <- gene
  }else{
    intersection <- intersect(x=intersection, y = gene)
  }
  n = n+1
}
write.csv(intersection,file="data/inter_RI.csv")

####annotation####
for(slide in slide_list){
  expr.url <- paste0("data/annotation/",slide,".csv")
  annot_table <- read.csv(file=expr.url,header=T)
  rownames(annot_table) <- annot_table[,1]
  annot_table$Barcode <- NULL
  spatial_list[[slide]] <- AddMetaData(object= spatial_list[[slide]],
                                       metadata = annot_table,
                                       col.name = paste0(colnames(annot_table),"_annot"))
}
#save
for(slide in slide_list){
  path <- paste0("data/result/",slide,"_matrix.csv")
  write.csv(spatial_list[[slide]]@assays$Spatial@counts,file=path)
  path <- paste0("data/result/",slide,"_annotation.csv")
  write.csv(spatial_list[[slide]]@meta.data["TLS_2_cat_annot"],file=path)
}

####DEG####
#merge data
library(harmony)
st <- merge(spatial_list[["a_3"]], y = spatial_list[c("b_1","b_18")],add.cell.ids = c("a_3","b_1","b_18"), merge.data = TRUE)
unique(sapply(X = strsplit(colnames(st), split = "_"), FUN = "[", 1))
st@meta.data[,"sample"] <- st@active.ident
#normalization
st <- NormalizeData(st, normalization.method = "LogNormalize", scale.factor = 1e4)
#dim reduction
st <- FindVariableFeatures(st,selection.method = 'vst', nfeatures = 2000)
st <- ScaleData(st)
st <- RunPCA(st, features = VariableFeatures(object = st))
#use harmony for integration
st <- RunHarmony(st, group.by.vars = "sample")
#DEG
gene<-FindMarkers(st,group.by="TLS_2_cat_annot",ident.1="TLS")
write.csv(gene,file="data/DEG_NRI.csv")
fcthreshold=2
pthreshold =0.05
a<-rownames(gene)
gene$change <- as.factor(ifelse(gene$p_val_adj<pthreshold & abs(gene$avg_log2FC)>log2(fcthreshold),
                                ifelse(gene$avg_log2FC>log2(fcthreshold),"Up","Down"),"Non"))
rownames(gene)<-a
gene$label <- ifelse(gene$p_val_adj<pthreshold & abs(gene$avg_log2FC)>log2(fcthreshold),as.character(rownames(gene)),"")

library(ggplot2)
library(ggrepel)

p.vol <- ggplot(data = gene,
                aes(x = avg_log2FC,y = -log10(p_val_adj),colour = change,fill = change))+
  scale_color_manual(values = c('green','grey','red'))+
  geom_point(alpha = 0.4,size = 3.5)+
  geom_text_repel(aes(x = avg_log2FC,y = -log10(p_val_adj),label = label),size = 5,
                  box.padding = unit(0.6,"lines"),point.padding = unit(0.7,"lines"),
                  segment.color = "black",show.legend = FALSE)+
  geom_vline(xintercept = c(-log2(fcthreshold),log2(fcthreshold)),lty = 4,col = "black",lwd = 0.8)+
  geom_hline(yintercept = -log10(pthreshold),lty = 4,col = "black",lwd = 0.8)+
  theme_bw()+
  labs(x = "log2(FoldChange)",y = "-log10(p-adjust)",title = "Volcano Plot of  Different Expression Proteins")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5,size = 20,face = "bold"), 
        legend.text = element_text(size = 20),legend.title = element_text(size = 20),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) 
ggsave(p.vol,filename = paste("Volcano_NRI.pdf"),width=10,height=10)
DEG <- gene[gene$p_val_adj<0.05,]
DEG <- DEG[abs(DEG$avg_log2FC)>1,]
write.csv(DEG,"DEG_NRI.csv")

