###########################################################
###########################################################
##########EdgeR differential Gene expression analysis######
###########################################################
####AUthor: T. Hamdi Kitapci
##########################################################
##########################################################

setwd("C:/Users/thkit_000/Desktop/miRNA-Analysis/Merged_First_and_Second_run/edgeR_based_on_UMI_with_Lynda_Filtering")
library(FactoMineR)
#library(FactoInvestigate)
library(factoextra)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)

library("edgeR")

input_file="1466.all_samples.summary_only_UMI_piRNA_removed.csv"

table=read.csv(input_file,row.names = 1)
#table=read.csv("All_HTSeq_using_gencode_annotation.csv",row.names = 1)


# Do this trick to get samples names without the "X"
dummy_table=read.csv(input_file,row.names = 1,check.names=F)
sample_names=colnames(dummy_table)
rm(dummy_table)

#table=subset(table, select=-c(symbol)) #remove the "symbol field from the data frame"

#table=table[rownames(table) != "MIMAT0000063", ]

#################Filter very high counts#####################
#countCheck <- table < 2000
# It is recommended to keep only those genes with a cpm > 1 in at least a number of samples equal to the number of samples in the smallest group.
#keep <- which(rowSums(countCheck) >= 1)
#number_of_rows_kept=nrow(as.table(keep))
#table=table[keep,]

#############################################################



##NOR -> 24, 54, 57, 58    Based on Lynda's correction on 1/22/2018
##DOR -> 1A, 31, 46, 56
replicates=c("DOR","DOR","DOR","DOR","NOR","NOR","NOR","NOR")

#########PCA##############
table_for_PCA=cpm(table) #apply CPM filtering on the raw data


##########Different filtering methods####################################

#Filter such that in DOR each sample has at least Threshold UMI count or in NOR each sample has at least Threshold UMI count

Threshold=6

Check=table_for_PCA>Threshold
keep=which( (rowSums(Check[,1:4]) >= 4) | (rowSums(Check[,5:8]) >= 4) )
table_for_PCA=table_for_PCA[keep,]

##########################################################################


rev_table=t(table_for_PCA)

my.pca=PCA(rev_table,graph=FALSE)
#Investigate(my.pca,file="miRNA_min_threshold_1_in_at_least_4_samples_investigate",document=c("word_document"))

#Build in function to plot PCA
fviz_pca_ind(my.pca,title="",repel=TRUE) #Generates a very noisy plot for this data

####USe ggplot2 for PCA plot

my_PCA_coordinates=as.data.frame(my.pca$ind$coord)
PC1=my_PCA_coordinates$Dim.1
PC2=my_PCA_coordinates$Dim.2
PC3=my_PCA_coordinates$Dim.3
variation_explained=as.data.frame(my.pca$eig)
PC1_variation=round(variation_explained$`percentage of variance`[1],digits=2)
PC2_variation=round(variation_explained$`percentage of variance`[2],digits=2)
PC3_variation=round(variation_explained$`percentage of variance`[3],digits=2)

p1 <- ggplot() 
p1 <- p1+ theme_bw()+coord_fixed()
p1 <- p1 + labs(title="PCA gene expression from cpm",x=paste("Principal component 1 %",PC1_variation),y=paste("Principal component 2 %",PC2_variation))+theme(plot.title = element_text(hjust = 0.5))
p1 <- p1 +  geom_point(data= my_PCA_coordinates, aes(x=Dim.1, y=Dim.2,fill=replicates), shape=22, size=8)
p1 <- p1 + geom_text(data=my_PCA_coordinates,aes(x=Dim.1,y=Dim.2,label=sample_names))

p2 <- ggplot()
p2 <- p2+ theme_bw()+coord_fixed()
p2 <- p2 + labs(title="PCA gene expression from cpm",x=paste("Principal component 2 %",PC2_variation),y=paste("Principal component 3 %",PC3_variation))+theme(plot.title = element_text(hjust = 0.5))
p2 <- p2 +  geom_point(data= my_PCA_coordinates, aes(x=Dim.2, y=Dim.3,fill=replicates), shape=22, size=8)
p2 <- p2 + geom_text(data=my_PCA_coordinates,aes(x=Dim.2,y=Dim.3,label=sample_names))

p3 <- ggplot() 
p3 <- p3+ theme_bw()+coord_fixed()
p3 <- p3 + labs(title="PCA gene expression from cpm",x=paste("Principal component 1 %",PC1_variation),y=paste("Principal component 3 %",PC3_variation))+theme(plot.title = element_text(hjust = 0.5))
p3 <- p3 +  geom_point(data= my_PCA_coordinates, aes(x=Dim.1, y=Dim.3,fill=replicates), shape=22, size=8)
p3 <- p3 + geom_text(data=my_PCA_coordinates,aes(x=Dim.1,y=Dim.3,label=sample_names))

grid.arrange(p1,p2,p3,ncol=3)




###############DGE analysis###################################################

myEdgeRList <- DGEList(counts=table, group = replicates,samples = sample_names)

###Filter miRNA with low counts
# Filtering out genes with very low counts:

#Filter based on Lynda's suggestion
countCheck <- cpm(myEdgeRList) > Threshold
# It is recommended to keep only those genes with a cpm > 1 in at least a number of samples equal to the number of samples in the smallest group.
#keep <- which(rowSums(countCheck) >= 4)
keep <- which( (rowSums(countCheck[,1:4]) >= 4) | (rowSums(countCheck[,5:8]) >= 4) )

#number_of_rows_kept=nrow(as.table(keep))
myEdgeRList=myEdgeRList[keep, keep.lib.sizes=FALSE]

myEdgeRList=calcNormFactors(myEdgeRList)

########PLOT MDS##########################################################################

#Use ggplot2 to plot MDS from edgeR
my_MDS_coordinates=as.data.frame(plotMDS(myEdgeRList))
p4 <- ggplot() 
p4 <- p4+ theme_bw()+coord_fixed()
p4 <- p4 + labs(title="",x="Leading logFC dim 1",y="Leading logFC dim 2")+theme(plot.title = element_text(hjust = 0.5))
p4 <- p4 +  geom_point(data= my_MDS_coordinates, aes(x=x, y=y,fill=replicates), shape=22, size=9)
p4 <- p4 + geom_text(data=my_MDS_coordinates,aes(x=x,y=y,label=sample_names))
p4
######################################################################################

###Estimate dispersion and plot BCV###############

#designMat=model.matrix(~replicates)

myEdgeRList<-estimateDisp(myEdgeRList)
#myEdgeRList <- estimateGLMCommonDisp(myEdgeRList)
#myEdgeRList <- estimateGLMTrendedDisp(myEdgeRList)
#myEdgeRList <- estimateGLMTagwiseDisp(myEdgeRList)
plotBCV(myEdgeRList)  #plot biological coefficient of variation

#####################################################
##### GLM fit is better used for complex experimental designs (more than 2 groups)
#myGLMFit <- glmFit(myEdgeRList, designMat)
#myLrt <- glmLRT(myGLMFit)  #logFC is log2-fold-changes!!
######################################################


#Use exactTest() for comparing 2 groups 
myExactTest=exactTest(myEdgeRList,c("NOR","DOR")) #Positive FC is upregulated in NOR!!!



############Make volcano plot###############################

### USE glmLRT as input
#myData=as.data.frame(myLrt$table)  


####uSe exacttest() as input for volcano plot
myData=myExactTest$table   



#plot(myData$logFC,-log2(myData$PValue))
a=expression("log"[2]~"Fold Change DOR/NOR")
b=expression("-log"[10]~"(p-value)")

not_significant="black"
significant_but_small_effect_size="gray"
upregulated_in_small="brown"
upregulated_in_large="green"

ss1=subset(myData,myData$PValue<=0.05)
upregulated=subset(ss1,ss1$logFC>=1)
downregulated=subset(ss1,ss1$logFC<1)
number_of_genes_upregulated=dim(upregulated)[1]
number_of_genes_downregulated=dim(downregulated)[1]

p5 <- ggplot() + geom_point(aes(myData$logFC,-log10(myData$PValue)), color=ifelse(myData$PValue>0.05,not_significant,ifelse(myData$logFC<=-1,upregulated_in_large,ifelse(myData$logFC>=1,upregulated_in_small,significant_but_small_effect_size))))
p5 <- p5 + labs(title="Volcano plot for miRNA expression before multiple testing correction",x=a,y=b)+theme(plot.title = element_text(hjust = 0.5))
p5 <- p5 + annotate("text",x = 3, y = 3, label =paste("n=",number_of_genes_upregulated) ,color=upregulated_in_small)
p5 <- p5 + annotate("text", x=-5,y=3, label =paste("n=",number_of_genes_downregulated) ,color=upregulated_in_large)
p5

############################################################


#edgeR_result <- topTags(myLrt)

edgeR_result <-topTags(myExactTest,dim(myExactTest$table)[1],adjust.method="BH") #save everything 
write.csv(edgeR_result$table,"miRNA_DE_list_based_on_UMI_Lynda_Filtering.csv",quote=FALSE)

correctedData=edgeR_result$table
####Volcano plot with corrected p-values#########
a=expression("log"[2]~"Fold Change DOR/NOR")
b=expression("-log"[10]~"(p-value)")

not_significant="black"
significant_but_small_effect_size="gray"
upregulated_in_small="brown"
upregulated_in_large="green"

ss1=subset(correctedData,correctedData$FDR<=0.05)
upregulated=subset(ss1,ss1$logFC>=1)
downregulated=subset(ss1,ss1$logFC<1)
number_of_genes_upregulated=dim(upregulated)[1]
number_of_genes_downregulated=dim(downregulated)[1]

p6 <- ggplot() + geom_point(aes(correctedData$logFC,-log10(correctedData$FDR)), color=ifelse(correctedData$FDR>0.05,not_significant,ifelse(correctedData$logFC<=-1,upregulated_in_large,ifelse(correctedData$logFC>=1,upregulated_in_small,significant_but_small_effect_size))))
p6 <- p6 + labs(title="Volcano plot for miRNA expression after FDR correction",x=a,y=b)+theme(plot.title = element_text(hjust = 0.5))
#p6 <- p6 + annotate("text",x = 2, y = 2, label =paste("n=",number_of_genes_upregulated) ,color=upregulated_in_small)
#p6 <- p6 + annotate("text", x=-2.5,y=2, label =paste("n=",number_of_genes_downregulated) ,color=upregulated_in_large)
p6



#deGenes <- decideTestsDGE(myLrt, p=0.001)
#deGenes <- rownames(myLrt)[as.logical(deGenes)] 
#plotSmear(myLrt, de.tags=deGenes)
#abline(h=c(-1, 1), col=2)




