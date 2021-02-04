#Calcualte the RPKM value of all expressed genes, and calculate the mean RPKM of NSCs

library("EDASeq")

ensembl_list <-ExpDataCPM %>%
  mutate(ensemblID = gsub('[.][0-9]+$','',.$ensemblID)) %>% pull(1)

#ensembl_length<-getGeneLengthAndGCContent(ensembl_list, "hsa")
#write.csv(ensembl_length,'ensemble_gene_length.csv')

ensembl_length<-read.csv('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/metabolic_analysis/ensemble_gene_length.csv')

cpmMatrixFiltered_unlog <- data.frame(gene_id=as.character(countMatrixFiltered$gene_id), Count2CPM(countMatrixFiltered[,-1]))

cpm_matrix_Ens<-left_join(cpmMatrixFiltered_unlog, geneNames %>% dplyr::select(gene_id, gene_name)) %>% 
  mutate(gene_id = gsub('[.][0-9]+$','',.$gene_id))



cpm_matrix_Ense<-left_join(ensembl_length[,c(1,2)],cpm_matrix_Ens,by=c('X'='gene_id'))

## Group ID
NSC_ctrl_id <- Metadata %>% filter(cell_type == 'NSC', mutation == 'CTRL') %>%
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)

NSC_polg_id <- Metadata %>% filter(cell_type == 'NSC', mutation == 'POLG') %>%
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)

IPS_ctrl_id <- Metadata %>% filter(cell_type %in% c('ESC','IPSC'), mutation == 'CTRL') %>%
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)

IPS_polg_id <- Metadata %>% filter(cell_type %in% c('IPSC'), mutation == 'POLG') %>%
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)

## FPKM
#change count of fragments to RFKM: RFKM ( log transformed)= CPM*1000/gene_length  + 1
ExpData_FPKM<-data.frame(ensemblID=cpm_matrix_Ense$X,GeneSymbol=cpm_matrix_Ense$gene_name,
                         log2(cpm_matrix_Ense[,c(3:(length(colnames(cpm_matrix_Ense))-1))]*1000/cpm_matrix_Ense[,2]+1))

## FPKM not log transformed
ExpData_FPKM<-data.frame(ensemblID=cpm_matrix_Ense$X,GeneSymbol=cpm_matrix_Ense$gene_name,
                         cpm_matrix_Ense[,c(3:(length(colnames(cpm_matrix_Ense))-1))]*1000/cpm_matrix_Ense[,2])


##TPM
rpk<-Count_matrix_Ense[,c(3:(length(colnames(Count_matrix_Ense))-1))]*1000/Count_matrix_Ense[,2]
scale_factor<-colSums(rpk,na.rm = TRUE)/1000000
ExpData_TPM<-data.frame(ensemblID=Count_matrix_Ense$X,GeneSymbol=Count_matrix_Ense$gene_name,
                        rpk/scale_factor)

NSC_ctrl_mean_tpm<-ExpData_TPM %>% dplyr::select(ensemblID,GeneSymbol,NSC_ctrl_id) %>% mutate(mean_rpkm = rowMeans(.[,c(3:13)])) %>% dplyr::select(ensemblID,mean_rpkm)
IPS_ctrl_mean_tpm<-ExpData_TPM  %>% dplyr::select(ensemblID,GeneSymbol,IPS_ctrl_id) %>% mutate(mean_rpkm = rowMeans(.[,c(3:12)])) %>% dplyr::select(ensemblID,mean_rpkm)
NSC_polg_mean_tpm<-ExpData_TPM  %>% dplyr::select(ensemblID,GeneSymbol,NSC_polg_id) %>% mutate(mean_rpkm = rowMeans(.[,c(3:21)])) %>% dplyr::select(ensemblID,mean_rpkm)
IPS_polg_mean_tpm<-ExpData_TPM  %>% dplyr::select(ensemblID,GeneSymbol,IPS_polg_id) %>% mutate(mean_rpkm = rowMeans(.[,c(3:21)])) %>% dplyr::select(ensemblID,mean_rpkm)

table(NSC_ctrl_mean_tpm$mean_rpkm>1)

mean_tpm<-cbind(IPS_ctrl_mean_cpm, IPS_polg_mean_cpm[,2], NSC_ctrl_mean_cpm[,2], NSC_polg_mean_cpm[,2])
colnames(mean_tpm)<-c('ensemblID','IPS_ctrl','IPS_polg','NSC_ctrl','NSC_polg')

write.table(mean_tpm,'/home/yu/Postdoc_project/NSC_RNAseq/RNA_seq/metabolic_analysis/mean_tpm_metabolic_all_group.txt',row.names = F,quote = F, sep='\t')
write.csv(mean_tpm,'/home/yu/Postdoc_project/NSC_RNAseq/RNA_seq/metabolic_analysis/mean_tpm_metabolic_all_group.csv',row.names = F,quote = F)
## import the RPKM to matlab to generate GEM model and extract the metabolic-metabolic interaction tables
#......

# read the matlab HumanGEM analyzied data

library(dplyr)
metabolic_ID<-read.table('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/metabolic_analysis/meta_ID_nsc_ctrl.txt',header =T, stringsAsFactors = F, sep = '\t')
NSC_metabolic_ineraction<-read.table('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/metabolic_analysis/meta_interaction_nsc_ctrl.txt',header =F, sep = '\t')

metabolic_ID<-metabolic_ID %>% mutate(Met_name=gsub(' ','-',.$Met_name))
metabolic_ID$compartment<-substr(metabolic_ID[,2], nchar(metabolic_ID[,2]), nchar(metabolic_ID[,2]))

# library(regex)
# metabolic_ID<-metabolic_ID%>% mutate(met_Name=gsub('(\\([0-9](.*)[A-Z]\\)(-)*)','',.$met_Name))

metabolic_ID$Meta_name<-paste(metabolic_ID$Met_name,'[',metabolic_ID$compartment,']',sep ='')
  
# summary(duplicated(metabolic_ID$meta_name))


## remove duplicate metabolic pairs, for example 1--2, 2--1, keep only 1 of the above
## use a hash table 
library('hash')
h <-hash()
for (i in c(1:length(NSC_metabolic_ineraction[,1]))) {
  x <- paste(NSC_metabolic_ineraction[i,1],'+',NSC_metabolic_ineraction[i,2])
  y <- paste(NSC_metabolic_ineraction[i,2],'+',NSC_metabolic_ineraction[i,1])
  if (x == y) { next} else if (has.key(y, h)) {next}
  else (.set(h, keys=x, values=y))
}

?hash()

metabolic_interact<-ls(h)
metabolic_interact<-matrix(unlist(strsplit(metabolic_interact,' + ',fixed = T)),ncol=2,byrow=T)
metabolic_interact<-data.frame(metabolic_interact,stringsAsFactors = F)
metabolic_interact$X1<-strtoi(metabolic_interact$X1)
metabolic_interact$X2<-strtoi(metabolic_interact$X2)

## change ID to name
library(dplyr)
metabolic_interact_name<-metabolic_interact %>% left_join(metabolic_ID[,c(1,5)],by=c('X1'='Num') ) %>% 
  left_join(metabolic_ID[,c(1,5)],by=c('X2'='Num') ) %>% dplyr::select(Meta_name.x,Meta_name.y)



## gene and metabolic classification
met_model_gene_ID<-read.table('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/metabolic_analysis/Gene_ID_nsc_ctrl.txt',header =T, stringsAsFactors = F, sep = '\t')

Gene_metabolic_classification<-read.table('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/metabolic_analysis/meta_gene_classfication_nsc_ctrl.txt',header = F, sep = '\t')

#match gene name and ensemble

#library(org.Hs.eg.db)
#k<-keys(org.Hs.eg.db,keytype='ENSEMBL')
#en2ENSE<-AnnotationDbi::select(org.Hs.eg.db,keys=k,columns=c("ENTREZID",'SYMBOL'),keytype = 'ENSEMBL')


genes2transcripts <- read.csv("/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/Data/counts_all_samples.gene_name_mappings",stringsAsFactors = F,header = T) %>%
  dplyr::select(gene_id, transcript_id, gene_name, transcript_type, gene_type, seqnames) %>%
  filter(!is.na(transcript_id)) %>%
  unique %>%
  dplyr::select(gene_id, gene_name, gene_type, seqnames) %>%
  mutate(gene_id= gsub('[.][0-9]+$','',.$gene_id)) %>%
  unique

# Gene_metabolic Nun to name

Gene_metabolic<-Gene_metabolic_classification %>% 
  left_join(met_model_gene_ID,by=c('V1'='Num')) %>% 
  left_join(metabolic_ID[,c(1,5)],by=c('V2'='Num') ) %>% 
  left_join(genes2transcripts[,c(1,2)],by=c('Gene_ID'='gene_id'))  %>% 
  dplyr::select(Meta_name,Gene_ID,gene_name)

## gene p values import CP2A
gene_DE_cp2a<-read.csv('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/Yu_tables/DE_NSC-CP2A_vs_NSC-CTRL-all.tsv',
                       sep='\t',header=T,stringsAsFactors = F)

gene_DE_cp2a<-gene_DE_cp2a%>%
  mutate(EnsemblID = gsub('[.][0-9]+$','',.$EnsemblID))%>% 
  dplyr::select(EnsemblID,GeneSymbol,log2FoldChange,padj,pvalue) 

gene_p_value_cp2a<-met_model_gene_ID %>% left_join(gene_DE_cp2a,by=c('Gene_ID'='EnsemblID'))


##
#library(clusterProfiler)
#x1<-enricher(gene_p_value_DE_down$gene_ID, pvalueCutoff = 0.05,pAdjustMethod = "BH",
#             qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500,
#             TERM2GENE = Gene_metabolic[,c(1,2)], TERM2NAME = Gene_metabolic[,c(1,1)])


## Alpers 

gene_DE_alpers<-read.csv('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/Yu_tables/DE_NSC-ALPERSvs_NSC-CTRL-all.tsv',
                       sep='\t',header=T,stringsAsFactors = F) 

gene_DE_alpers<-gene_DE_alpers%>% 
  mutate(EnsemblID = gsub('[.][0-9]+$','',.$EnsemblID))%>% 
  dplyr::select(EnsemblID,GeneSymbol,log2FoldChange,padj,pvalue) 

gene_p_value_alpers<-met_model_gene_ID %>% left_join(gene_DE_alpers,by=c('gene_ID'='EnsemblID'))

## CP2A and Alpers 
gene_DE_cp2a_alpers<-read.csv('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/Yu_tables/DE_NSC-CP2A-Alpers_vs_NSC-CTRL-all.tsv',
                         sep='\t',header=T,stringsAsFactors = F) 

gene_DE_cp2a_alpers<-gene_DE_cp2a_alpers%>% 
  mutate(EnsemblID = gsub('[.][0-9]+$','',.$EnsemblID))%>%  
  dplyr::select(EnsemblID,GeneSymbol,log2FoldChange,padj,pvalue) 

gene_p_value_cp2a_alpers<-met_model_gene_ID %>% left_join(gene_DE_cp2a_alpers,by=c('Gene_ID'='EnsemblID'))


## WS5A, CP2A and Alpers 
gene_DE_polg_ctrl<-read.csv('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/Yu_tables/DE_NSC-All_polg_vs_NSC-CTRL-all.tsv',
                              sep='\t',header=T,stringsAsFactors = F) 

gene_DE_polg_ctrl<-gene_DE_polg_ctrl%>% 
  mutate(EnsemblID = gsub('[.][0-9]+$','',.$EnsemblID))%>%  
  dplyr::select(EnsemblID,GeneSymbol,log2FoldChange,padj,pvalue) 

gene_p_value_polg_ctrl<-met_model_gene_ID %>% left_join(gene_DE_polg_ctrl,by=c('Gene_ID'='EnsemblID'))



## use piano to conduct GSA analysis and draw metabolic network

library(piano)

#use piano, but first delete the genes without pvalues or fold changes
gene_DE_cp2a<-gene_DE_cp2a %>% na.omit(.[,c(3,4)])
gene_DE_alpers<-gene_DE_alpers %>% na.omit(.[,c(3,4)])
gene_DE_cp2a_alpers<-gene_DE_cp2a_alpers %>% na.omit(.[,c(3,4)])
gene_DE_polg_ctrl<-gene_DE_polg_ctrl %>% na.omit(.[,c(3,4)])

# pvalue data
gene_p_id_cp2a<-gene_DE_cp2a[,4]
names(gene_p_id_cp2a)<-gene_DE_cp2a$EnsemblID

gene_p_id_alpers<-gene_DE_alpers[,4]
names(gene_p_id_alpers)<-gene_DE_alpers$EnsemblID

gene_p_id_cp2a_alpers<-gene_DE_cp2a_alpers[,4]
names(gene_p_id_cp2a_alpers)<-gene_DE_cp2a_alpers$EnsemblID

gene_p_id_polg_ctrl<-gene_DE_polg_ctrl[,4]
names(gene_p_id_polg_ctrl)<-gene_DE_polg_ctrl$EnsemblID

# fold change data
gene_fc_id_cp2a<-gene_DE_cp2a[,3]
names(gene_fc_id_cp2a)<-gene_DE_cp2a$EnsemblID

gene_fc_id_alpers<-gene_DE_alpers[,3]
names(gene_fc_id_alpers)<-gene_DE_alpers$EnsemblID

gene_fc_id_cp2a_alpers<-gene_DE_cp2a_alpers[,3]
names(gene_fc_id_cp2a_alpers)<-gene_DE_cp2a_alpers$EnsemblID

gene_fc_id_polg_ctrl<-gene_DE_polg_ctrl[,3]
names(gene_fc_id_polg_ctrl)<-gene_DE_polg_ctrl$EnsemblID

# gene ID and metabolic data set
mygsc<-loadGSC(Gene_metabolic[,c(2,1)])

plot(c(1,2,3))

# run GSA analysis
gsaRes_cp2a <- runGSA(gene_p_id_cp2a,gene_fc_id_cp2a,gsc=mygsc, 
                 geneSetStat="reporter",
                 signifMethod="nullDist", 
                 nPerm=1000,
                 gsSizeLim=c(10,1000))

gsaRes_polg_ctrl <- runGSA(gene_p_id_polg_ctrl,gene_fc_id_polg_ctrl,gsc=mygsc, 
                      geneSetStat="reporter",
                      signifMethod="nullDist", 
                      nPerm=1000,
                      gsSizeLim=c(10,1000))
# GSAsummaryTable(gsaRes_cp2a, save=TRUE, file="gsaRes_cp2a_Tab.xls")
dev.on()

dev.new()
nw_cp2a <- networkPlot(gsaRes_cp2a, class="distinct", direction="both",
                       significance=0.0001,lay=layout_with_fr,edgeWidth = c(0, 15),overlap=3)

?networkPlot
networkPlot(gsaRes_cp2a, class="non")
?networkPlot
?networkPlot2

require(visNetwork)
cp2a<-visNetwork(nw_cp2a$x$nodes,nw_cp2a$x$edges)
View(nw_cp2a$x$nodes)
View(nw_cp2a$x$edges)
?visExport
visExport( cp2a,
  type = "pdf",
  name = "network"
)



dev.next()
gsaRes_alpers <- runGSA(gene_p_id_alpers,gene_fc_id_alpers,gsc=mygsc, 
                        geneSetStat="reporter",
                        signifMethod="nullDist", 
                        nPerm=1000,
                        gsSizeLim=c(5,1000),
                        adjMethod = "fdr")
nw_alpers <- networkPlot(gsaRes_alpers, class="distinct", direction="both",
                         significance=0.0001,lay=layout_with_kk,edgeWidth = c(0, 15),overlap=3)
?networkPlot
?runGSA
nw_alpers
addExport(nw_alpers, pdf = TRUE)
?addExport
getwd()



gsaRes_cp2a_alpers <- runGSA(gene_p_id_cp2a_alpers,gene_fc_id_cp2a_alpers,gsc=mygsc, 
                        geneSetStat="reporter",
                        signifMethod="nullDist", 
                        nPerm=1000,
                        gsSizeLim=c(10,1000),
                        adjMethod = "fdr")
nw_cp2a_alpers <- networkPlot(gsaRes_cp2a_alpers, class="distinct", direction="both",
                          significance=0.0001,lay=layout_with_fr,edgeWidth = c(0, 15),overlap=3)
nw_cp2a_alpers
nw_cp2a_alpers$x$nodes
View(nw_cp2a_alpers$x$nodes)

nw_cp2a_alpers$x$edges

View(nw_cp2a_alpers$x$edges)
nw_cp2a_alpers

gsaRes_polg_ctrl <- runGSA(gene_p_id_polg_ctrl,gene_fc_id_polg_ctrl,gsc=mygsc, 
                             geneSetStat="reporter",
                             signifMethod="nullDist", 
                             nPerm=1000,
                             gsSizeLim=c(10,1000),
                             adjMethod = "fdr")
nw_polg_ctrl <- networkPlot(gsaRes_polg_ctrl, class="distinct", direction="both",
                              significance=0.0001,lay = 2,edgeWidth = c(0, 15),overlap=3)

nw_polg_ctrl


library(igraph)

nw_cp2a_alpers_igra<-graph_from_data_frame(d=nw_cp2a_alpers$x$edges, vertices=nw_cp2a_alpers$x$nodes, directed=F)

dev.off()
plot(nw_cp2a_alpers_igra)

nw_cp2a_alpers_igra
V(nw_cp2a_alpers_igra)$size = nw_cp2a_alpers$x$nodes$size
V(nw_cp2a_alpers_igra)$color = nw_cp2a_alpers$x$nodes$color
V(nw_cp2a_alpers_igra)$label = nw_cp2a_alpers$x$nodes$label
V(nw_cp2a_alpers_igra)$label.color = 'black'

E(nw_cp2a_alpers_igra)$width = nw_cp2a_alpers$x$edges$width
E(nw_cp2a_alpers_igra)$color = nw_cp2a_alpers$x$edges$color

plot(nw_cp2a_alpers_igra,edge.curved=.1,vertex.label.cex=.7,layout=layout_with_gem)


write_graph(nw_cp2a_alpers_igra, 'nw_cp2a_alpers_igra.gml', format = c("gml"))

abc<-createNetworkFromIgraph(nw_cp2a_alpers_igra,"myIgraph")


?createNetworkFromIgraph
getwd()
library(RCy3)

nw_alpers
?
nw_alpers
?networkPlot2
nw_alpers$geneSets

GSAsummaryTable(gsaRes)
gene_p_value_cp2a
Gene_metabolic


## prepare the kiwi dataset
# gene set interaction

write.table(metabolic_interact_name,
            file= 'kiwi_meta_inter_name.txt',
            sep = '\t',quote = F, 
            row.names = F,col.names=F)


writeFilesForKiwi(gsaRes_cp2a,label='kiwi_cp2a_gs')
