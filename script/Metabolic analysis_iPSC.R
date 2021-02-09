#Calcualte the RPKM value of all expressed genes, and calculate the mean RPKM of NSCs


# read the matlab HumanGEM analyzied data

library(dplyr)
metabolic_ID<-read.table('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/metabolic_analysis/metabolic/meta_ID_ips_ctrl.txt',header =T, stringsAsFactors = F, sep = '\t')
IPS_metabolic_ineraction<-read.table('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/metabolic_analysis/metabolic/meta_interaction_ips_ctrl.txt',header =F, sep = '\t')

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
for (i in c(1:length(IPS_metabolic_ineraction[,1]))) {
  x <- paste(IPS_metabolic_ineraction[i,1],'+',IPS_metabolic_ineraction[i,2])
  y <- paste(IPS_metabolic_ineraction[i,2],'+',IPS_metabolic_ineraction[i,1])
  if (x == y) { next} else if (has.key(y, h)) {next}
  else (.set(h, keys=x, values=y))
}

metabolic_interact<-ls(h)
metabolic_interact<-matrix(unlist(strsplit(metabolic_interact,' + ',fixed = T)),ncol=2,byrow=T)
metabolic_interact<-data.frame(metabolic_interact,stringsAsFactors = F)
metabolic_interact$X1<-strtoi(metabolic_interact$X1)
metabolic_interact$X2<-strtoi(metabolic_interact$X2)

## change ID to name
metabolic_interact_name<-metabolic_interact %>% left_join(metabolic_ID[,c(1,5)],by=c('X1'='Num') ) %>% 
  left_join(metabolic_ID[,c(1,5)],by=c('X2'='Num') ) %>% dplyr::select(Meta_name.x,Meta_name.y)



## gene and metabolic classification
met_model_gene_ID<-read.table('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/metabolic_analysis/metabolic/Gene_ID_ips_ctrl.txt',header =T, stringsAsFactors = F, sep = '\t')

Gene_metabolic_classification<-read.table('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/metabolic_analysis/metabolic/meta_gene_classfication_ips_ctrl.txt',header = F, sep = '\t')

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

## gene p values import WS5A, CP2A compare control
gene_DE_polg_ctrl<-read.csv('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/Yu_tables/DE_IPS-CP2A-WS5A_vs_IPS-CTRL.tsv',
                            sep='\t',header=T,stringsAsFactors = F) 

gene_DE_polg_ctrl<-gene_DE_polg_ctrl%>% 
  mutate(EnsemblID = gsub('[.][0-9]+$','',.$EnsemblID))%>%  
  dplyr::select(EnsemblID,GeneSymbol,log2FoldChange,padj,pvalue) 

gene_p_value_polg_ctrl<-met_model_gene_ID %>% left_join(gene_DE_polg_ctrl,by=c('Gene_ID'='EnsemblID'))



gene_DE_alpers<-read.csv('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/Yu_tables/DE_IPSC-ALPERSvs_IPSC-CTRL-all.tsv',
                         sep='\t',header=T,stringsAsFactors = F) 

gene_DE_alpers<-gene_DE_alpers%>% 
  mutate(EnsemblID = gsub('[.][0-9]+$','',.$EnsemblID))%>% 
  dplyr::select(EnsemblID,GeneSymbol,log2FoldChange,padj,pvalue) 

gene_p_value_alpers<-met_model_gene_ID %>% left_join(gene_DE_alpers,by=c('Gene_ID'='EnsemblID'))


## use piano to conduct GSA analysis and draw metabolic network

library(piano)

#use piano, but first delete the genes without pvalues or fold changes

gene_DE_polg_ctrl<-gene_DE_polg_ctrl %>% na.omit(.[,c(3,4)])

gene_DE_alpers<-gene_DE_alpers %>% na.omit(.[,c(3,4)])

# pvalue data

gene_p_id_polg_ctrl<-gene_DE_polg_ctrl[,4]
names(gene_p_id_polg_ctrl)<-gene_DE_polg_ctrl$EnsemblID

gene_p_id_alpers<-gene_DE_alpers[,4]
names(gene_p_id_alpers)<-gene_DE_alpers$EnsemblID

# fold change data

gene_fc_id_polg_ctrl<-gene_DE_polg_ctrl[,3]
names(gene_fc_id_polg_ctrl)<-gene_DE_polg_ctrl$EnsemblID

gene_fc_id_alpers<-gene_DE_alpers[,3]
names(gene_fc_id_alpers)<-gene_DE_alpers$EnsemblID

# gene ID and metabolic data set
mygsc<-loadGSC(Gene_metabolic[,c(2,1)])

# run GSA analysis

gsaRes_polg_ctrl <- runGSA(gene_p_id_polg_ctrl,gene_fc_id_polg_ctrl,gsc=mygsc, 
                           geneSetStat="reporter",
                           signifMethod="nullDist", 
                           nPerm=1000,
                           gsSizeLim=c(10,1000),
                           adjMethod = "fdr")

nw_polg_ctrl <- networkPlot(gsaRes_polg_ctrl, class="distinct", direction="both",
                            significance=0.0001,lay = 2,edgeWidth = c(0, 15),overlap=3)
nw_polg_ctrl
nw_polg_ctrl <- networkPlot2(gsaRes_polg_ctrl, class="distinct", direction="both",
                            significance=0.0001,lay = 2,edgeWidth = c(0, 15),overlap=3)

nw_polg_ctrl


gsaRes_alpers <- runGSA(gene_p_id_alpers,gene_fc_id_alpers,gsc=mygsc, 
                        geneSetStat="reporter",
                        signifMethod="nullDist", 
                        nPerm=1000,
                        gsSizeLim=c(5,1000),
                        adjMethod = "fdr")
nw_alpers <- networkPlot(gsaRes_alpers, class="distinct", direction="both",
                         significance=0.0001,lay=5,edgeWidth = c(0, 15),overlap=3)


