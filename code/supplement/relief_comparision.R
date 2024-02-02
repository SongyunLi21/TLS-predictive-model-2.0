####RI####
library(readr)
relief <- read_csv("RI_RELIEF_1.csv")
value <- as.numeric(relief$x)
table(cut(value,c(-20,-10,-1,0,1,10,20,max(value))))
#we use 20 as the thershold, consider the value relief>20 as important 
feature<-relief[relief$x>20,]
colnames(feature)<-c("gene","value")
print(paste0("The number of gene selected by relief is ",length(feature$gene)))
print("The gene selected by relief is")
print(as.character(feature$gene))
Chi_feature <- c('VIM','ENO1','TGFBI','SERPINE1','TMSB4X','IL32','MS4A1','IGHG3','UBC',
               'IGKV4-1','C7','BLK','CCL19','IGHG1','FCRL1','GLUL','MT2A','HSPB1','BANK1',
               'JCHAIN','IGKC','IGLV3-1','CTSD','FTL','MT1E','IGLC1','B2M','FN1','IGHG2','ACTB')
print("The gene selected by Chi-square test is")
print(Chi_feature)
inter <- intersect(x=Chi_feature, y = feature$gene)
print("The gene shared by CHi-square test and Relif is:")
print(inter)
print(paste0("Account for ",(length(inter)*100)/length(Chi_feature),"% of CHi-square feature"))
print(paste0("Account for ",(length(inter)*100)/length(feature$gene), "% of Relief feature "))

####NRI####
library(readr)
relief <- read_csv("NRI_RELIEF_1.csv")
value <- as.numeric(relief$x)
table(cut(value,c(-20,-10,-1,0,1,10,20,max(value))))
#we use 20 as the thershold, consider the value relief>20 as important 
feature<-relief[relief$x>20,]
colnames(feature)<-c("gene","value")
print(paste0("The number of gene selected by relief is ",length(feature$gene)))
print("The gene selected by relief is")
print(as.character(feature$gene))
Chi_feature <- c('VIM','RPL36','RPS27','SPP1','NDRG1','HLA-B','CD74','RPS21',
                 'GAPDH','CD24','B2M','RPL34','HLA-A','FTH1','RPS18','EEF1A1',
                 'RPL37A','RPS23','RPL41','RPL10','IGHG2','ITM2B','IGKC','MIF',
                 'RPL39','FTL','RPL13','TGFBI','IGLC3','TMSB10','RPS8','IGFBP7',
                 'TPT1','RPL37','IGLC2','RPLP1','IGHG1','RPS2')
print("The gene selected by Chi-square test is")
print(Chi_feature)
inter <- intersect(x=Chi_feature, y = feature$gene)
print("The gene shared by CHi-square test and Relif is:")
print(inter)
print(paste0("Account for ",(length(inter)*100)/length(Chi_feature),"% of CHi-square feature"))
print(paste0("Account for ",(length(inter)*100)/length(feature$gene), "% of Relief feature "))