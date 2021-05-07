library(magrittr)
library(tidyverse)
library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)



# extract ensemble and entrze gene ID mapping files
k<-keys(org.Hs.eg.db,keytype='ENSEMBL')
en2ENSE<-AnnotationDbi::select(org.Hs.eg.db,keys=k,
                               columns=c("ENTREZID",'SYMBOL'),
                               keytype = 'ENSEMBL')



# import the target gene list to be analyzed


#import DE results, import function
import_de<-function(a) {
  de_import<-read.table(a, sep= '\t',header =T,stringsAsFactors = F) %>%
    mutate(EnsemblID = gsub('[.][0-9]+$','',.$EnsemblID))%>% 
    dplyr::select(EnsemblID,log2FoldChange,padj,pvalue) %>% 
    left_join(en2ENSE[,1:2],by=(c('EnsemblID'='ENSEMBL')))
  return(de_import)
}

# import background genes

Gene_ID_ENTREZ<- ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  left_join(en2ENSE[,1:2],by=c('ensemblID'='ENSEMBL')) %>% dplyr::select(GeneSymbol,ENTREZID)

bg_genes<-Gene_ID_ENTREZ %>% dplyr::select(2) %>% pull() %>% unique()
bg_genes_ense<-ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  dplyr::select(ensemblID) %>% pull() %>% unique() 


# DE analysis

require("DESeq2")
require("colorspace")

GenericHumanAnno <- GetAnnoFiles("Generic_human")

#Group 1 WS5A and CP2A IPSC


contrast_id <- "IPS-CP2A_vs_IPS-WS5A"
contrast_str <- "subname_WS5A_vs_CP2A"
Mt <- Metadata %>% filter(cell_type %in% c('IPSC','ESC'), (subname_mut == 'WS5A' | subname_mut == 'CP2A')) %>%
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::rename(subname = subname_mut)
Mt %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

Ct <- Counts[,colnames(Counts) %in% Mt$sample_id]
rownames(Ct) <- Counts$gene_id
Md <- as.formula(paste0("~subname"))
dds <- DESeq2RUN(Ct, Mt, Md)

res_ws5a_cp2a_ips<-results(dds, alpha = 0.05, contrast=c('subname', 'CP2A','WS5A'),format = "DataFrame", independentFiltering = T)
res_ws5a_cp2a_ips$GeneSymbol <- geneNames$gene_name[match(rownames(res_ws5a_cp2a_ips), geneNames$gene_id)]
res_ws5a_cp2a_ips$EnsemblID <- rownames(res_ws5a_cp2a_ips)
res_ws5a_cp2a_ips %<>% data.frame %>% filter(GeneSymbol != "") %>% 
  dplyr::select(GeneSymbol, EnsemblID, log2FoldChange, pvalue, padj, everything()) %>%arrange(padj)


xlim <- c(1,1e5); ylim <- c(-15,15)
plotMA(res_ws5a_cp2a_ips, xlim=xlim, ylim=ylim, main="ipsc cp2a VS ws5a",alpha = 0.05)

ws5a_cp2a_ips_de_len<-length(res_ws5a_cp2a_ips[res_ws5a_cp2a_ips$padj<0.05,1]) 

de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_ws5a_cp2a_ips, de_file)

#Group 2 WS5A and CP2A NSC

contrast_id <- "NSC-CP2A_vs_NSC-WS5A"
contrast_str <- "subname_WS5A_vs_CP2A"
Mt <- Metadata %>% filter(cell_type %in% c('NSC'), (subname_mut == 'WS5A' | subname_mut == 'CP2A')) %>%
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::rename(subname = subname_mut)
Mt %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)
Ct <- Counts[,colnames(Counts) %in% Mt$sample_id]
rownames(Ct) <- Counts$gene_id
Md <- as.formula(paste0("~subname"))
dds <- DESeq2RUN(Ct, Mt, Md)

res_ws5a_cp2a_nsc<-results(dds, alpha = 0.05, contrast=c('subname', 'CP2A','WS5A'),format = "DataFrame", independentFiltering = T)

plotMA(res_ws5a_cp2a_nsc, xlim=xlim, ylim=ylim, main="nsc cp2a VS ws5a",alpha = 0.05)

res_ws5a_cp2a_nsc$GeneSymbol <- geneNames$gene_name[match(rownames(res_ws5a_cp2a_nsc), geneNames$gene_id)]
res_ws5a_cp2a_nsc$EnsemblID <- rownames(res_ws5a_cp2a_nsc)
res_ws5a_cp2a_nsc %<>% data.frame %>% filter(GeneSymbol != "") %>% 
  dplyr::select(GeneSymbol, EnsemblID, log2FoldChange, pvalue, padj, everything()) %>%arrange(padj)

de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_ws5a_cp2a_nsc, de_file)


##plot MA plot

par(mfrow=c(1,1), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-15,15)
plotMA(res_ws5a_cp2a_nsc_1, xlim=xlim, ylim=ylim, main="NSC",alpha = 0.05)
plotMA(res_ws5a_cp2a_ips_1, xlim=xlim, ylim=ylim, main="NSC",alpha = 0.05)
plotMA(res_cp2a_ips_1, xlim=xlim, ylim=ylim, main="CP2A iPSC",alpha = 0.05)
plotMA(res_CP2A_nsc1, xlim=xlim, ylim=ylim, main="CP2A NSC",alpha = 0.05)


## enrichment analysis

# import DE analysis results
ws5a_cp2a_nsc<-import_de('./Yu_tables/DE_NSC-CP2A_vs_NSC-WS5A.tsv')
ws5a_cp2a_ipsc<-import_de("./Yu_tables/DE_IPS-CP2A_vs_IPS-WS5A.tsv")

de_ws5a_cp2a_list<-list(nsc_up = ws5a_cp2a_nsc[,c(1,5)][ws5a_cp2a_nsc$padj<0.05 & ws5a_cp2a_nsc$log2FoldChange>0, ],
                        nsc_down = ws5a_cp2a_nsc[,c(1,5)][ws5a_cp2a_nsc$padj<0.05 & ws5a_cp2a_nsc$log2FoldChange<0, ],
                        ipsc_up = ws5a_cp2a_ipsc[,c(1,5)][ws5a_cp2a_ipsc$padj<0.05 & ws5a_cp2a_ipsc$log2FoldChange>0, ],
                        ipsc_down = ws5a_cp2a_ipsc[,c(1,5)][ws5a_cp2a_ipsc$padj<0.05 & ws5a_cp2a_ipsc$log2FoldChange<0, ])

de_ws5a_cp2a_list['nsc_up']
de_ws5a_cp2a_list['nsc_down']
de_ws5a_cp2a_list['ipsc_up']
de_ws5a_cp2a_list['ipsc_down']


## KEGG analysis

KEGG_ws5a_cp2a_enrich<-function(a) {
  enrich <- enrichKEGG(de_ws5a_cp2a_list[[a]]$ENTREZID, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05, pAdjustMethod = "BH",qvalueCutoff = 1,
                       universe = bg_genes)
  enrich <- setReadable(enrich, org.Hs.eg.db, 
                        keyType = "ENTREZID")
  return(enrich)
}

KEG_wc_n_up<-KEGG_ws5a_cp2a_enrich('nsc_up')
KEG_wc_n_down<-KEGG_ws5a_cp2a_enrich('nsc_down')
KEG_wc_i_up<-KEGG_ws5a_cp2a_enrich('ipsc_up')
KEG_wc_i_down<-KEGG_ws5a_cp2a_enrich('ipsc_down')

View(KEG_wc_i_up@result)

barplot(KEGG_ws5a_cp2a_enrich('nsc_up'),showCategory = 30)
heatplot(KEGG_ws5a_cp2a_enrich('nsc_down'),showCategory = 28)

heatplot(KEGG_ws5a_cp2a_enrich('ipsc_up'))
heatplot(KEGG_ws5a_cp2a_enrich('ipsc_down'))

## GO analysis
GO_ws5a_cp2a_enrich<-function(a) {
  ego <- enrichGO(gene          = de_ws5a_cp2a_list[[a]]$ENTREZID,
                  universe      = bg_genes,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  minGSSize = 20,
                  maxGSSize = 150,
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  return(ego)
}


GO_all_ws5a_cp2a_nsc_up<-GO_ws5a_cp2a_enrich('nsc_up')
GO_all_ws5a_cp2a_nsc_down<-GO_ws5a_cp2a_enrich('nsc_down')
GO_all_ws5a_cp2a_ipsc_up<-GO_ws5a_cp2a_enrich('ipsc_up')
GO_all_ws5a_cp2a_ipsc_down<-GO_ws5a_cp2a_enrich('ipsc_down')
View(GO_all_ws5a_cp2a_nsc_down@result)
View(GO_all_ws5a_cp2a_ipsc_down@result)


barplot(GO_ws5a_cp2a_enrich('nsc_up'))
heatplot(GO_ws5a_cp2a_enrich('nsc_down'))

heatplot(GO_ws5a_cp2a_enrich('ipsc_up'))
heatplot(GO_ws5a_cp2a_enrich('ipsc_down'))


## KEGG module analysis
MKEGG_ws5a_cp2a_enrich<-function(a) {
  enrich <- enrichMKEGG(de_ws5a_cp2a_list[[a]]$ENTREZID, organism = "hsa", minGSSize=1,
                        universe = bg_genes)
  enrich <- setReadable(enrich, org.Hs.eg.db, 
                        keyType = "ENTREZID")
  return(enrich)
}

heatplot(MKEGG_ws5a_cp2a_enrich('nsc_up'))
heatplot(MKEGG_ws5a_cp2a_enrich('nsc_down'))

heatplot(MKEGG_ws5a_cp2a_enrich('ipsc_up'))
heatplot(MKEGG_ws5a_cp2a_enrich('ipsc_down'))

##import wiki pathway annotation file 

wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%") %>% 
  filter(gene %in% Gene_ID_ENTREZ$ENTREZID)

wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME


## wiki enrichment analysis
wiki_ws5a_cp2a_enrich<-function(a) {
  enrich <- enricher(de_ws5a_cp2a_list[[a]]$ENTREZID, pvalueCutoff = 0.05,pAdjustMethod = "BH",
                     qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 200,
                     TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
  enrich <- setReadable(enrich, org.Hs.eg.db, 
                        keyType = "ENTREZID")
  return(enrich)
}

heatplot(wiki_ws5a_cp2a_enrich('nsc_up'),showCategory = 23)
barplot(wiki_ws5a_cp2a_enrich('nsc_down'),showCategory = 12)

wiki_ws5a_cp2a_enrich('ipsc_up')
wiki_ws5a_cp2a_enrich('ipsc_down')

## human mito carta 2.o database

## ID table
human_mito_carta<-read.table('./human_mito_carta.txt',sep='\t',header=TRUE)
colnames(human_mito_carta)<-c('ENTREZID','Ensembl')

## full_annotation table
human_mito_carta_full<-read.csv('./Human.MitoCarta2.0.csv',sep='>',header=TRUE)

View(human_mito_carta_full)


## instersect with DE up and down regulated gene list, and generate the mito_de gene information table, 
## note that de results is WS5A/CP2A, so in the table reverse the direction of log2foldchange

de_ws5a_cp2a_ipsc_up_mito<-as.character(intersect(de_ws5a_cp2a_list[['ipsc_up']]$EnsemblID,human_mito_carta$Ensembl))


ws5a_cp2a_ipsc_mito_down_table<-ws5a_cp2a_ipsc %>%  filter(EnsemblID %in% de_ws5a_cp2a_ipsc_up_mito) %>% 
  left_join(human_mito_carta_full[,c(4,6,7,10)],by = c("EnsemblID" = "EnsemblGeneID")) %>%
  mutate(log2FoldChange=-(.$log2FoldChange)) %>%
  dplyr::select(Symbol,EnsemblID,Description,MitoDomain_Score,log2FoldChange,padj)


de_ws5a_cp2a_ipsc_down_mito<-as.character(intersect(de_ws5a_cp2a_list[['ipsc_down']]$EnsemblID,human_mito_carta$Ensembl))

ws5a_cp2a_ipsc_mito_up_table<-ws5a_cp2a_ipsc %>%  filter(EnsemblID %in% de_ws5a_cp2a_ipsc_down_mito) %>% 
  left_join(human_mito_carta_full[,c(4,6,7,10)],by = c("EnsemblID" = "EnsemblGeneID")) %>%
  mutate(log2FoldChange=-(.$log2FoldChange)) %>%
  dplyr::select(Symbol,EnsemblID,Description,MitoDomain_Score,log2FoldChange,padj)

write.csv(ws5a_cp2a_ipsc_mito_up_table,'ws5a_cp2a_ipsc_mito_up_table.csv', col.names=T)


de_ws5a_cp2a_nsc_up_mito<-as.character(intersect(de_ws5a_cp2a_list[['nsc_up']]$EnsemblID,human_mito_carta$Ensembl))

ws5a_cp2a_nsc_mito_down_table<-ws5a_cp2a_nsc %>%  filter(EnsemblID %in% de_ws5a_cp2a_nsc_up_mito) %>% 
  left_join(human_mito_carta_full[,c(4,6,7,10)],by = c("EnsemblID" = "EnsemblGeneID")) %>%
  dplyr::select(Symbol,EnsemblID,Description,MitoDomain_Score,log2FoldChange,padj)

write.csv(ws5a_cp2a_nsc_mito_down_table,'ws5a_cp2a_nsc_mito_down_table.csv', col.names=T)


de_ws5a_cp2a_nsc_down_mito<-as.character(intersect(de_ws5a_cp2a_list[['nsc_down']]$EnsemblID,human_mito_carta$Ensembl))

ws5a_cp2a_nsc_mito_up_table<-ws5a_cp2a_nsc %>%  filter(EnsemblID %in% de_ws5a_cp2a_nsc_down_mito) %>% 
  left_join(human_mito_carta_full[,c(4,6,7,10)],by = c("EnsemblID" = "EnsemblGeneID")) %>%
  mutate(log2FoldChange=-(.$log2FoldChange)) %>%
  dplyr::select(Symbol,EnsemblID,Description,MitoDomain_Score,log2FoldChange,padj)

write.csv(ws5a_cp2a_nsc_mito_up_table,'ws5a_cp2a_nsc_mito_up_table.csv', col.names=T)

## export table with mito annotation data


## barplot of up and down regulated genes and mitochondrial related genes
length(de_ws5a_cp2a_ipsc_up_mito)
length(de_ws5a_cp2a_ipsc_down_mito)
length(de_ws5a_cp2a_nsc_up_mito)
length(de_ws5a_cp2a_nsc_down_mito)




length(na.omit(de_ws5a_cp2a_list[['ipsc_up']]$EnsemblID))
length(na.omit(de_ws5a_cp2a_list[['ipsc_down']]$EnsemblID))
length(na.omit(de_ws5a_cp2a_list[['nsc_up']]$EnsemblID))
length(na.omit(de_ws5a_cp2a_list[['nsc_down']]$EnsemblID))

 ## ipsc 
DE_number<-data.frame(Group=c('Up','Down','Up','Down'),
                      type =c('mito','mito','all','all'),
                      DE_gene_number=c(length(de_ws5a_cp2a_ipsc_up_mito),
                                       length(de_ws5a_cp2a_ipsc_down_mito),
                                       length(na.omit(de_ws5a_cp2a_list[['ipsc_up']]$EnsemblID)),
                                       length(na.omit(de_ws5a_cp2a_list[['ipsc_down']]$EnsemblID))))


p1 <- ggplot(DE_number, aes(x=Group, y=DE_gene_number,fill=type)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_classic() + 
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")) + 
  scale_fill_brewer(palette="Blues")

 ## nsc
DE_number<-data.frame(Group=c('Up','Down','Up','Down'),
                      type =c('mito','mito','all','all'),
                      DE_gene_number=c(length(de_ws5a_cp2a_nsc_up_mito),
                                       length(de_ws5a_cp2a_nsc_down_mito),
                                       length(na.omit(de_ws5a_cp2a_list[['nsc_up']]$EnsemblID)),
                                       length(na.omit(de_ws5a_cp2a_list[['nsc_down']]$EnsemblID))))

p2 <- ggplot(DE_number, aes(x=Group, y=DE_gene_number,fill=type)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_classic() + 
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")) +
  scale_fill_brewer(palette="Blues")


require(gridExtra)
grid.arrange(p1,p2, ncol=2)

## Venn plot
library(VennDiagram)

venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Set 1" , "Set 2 " , "Set 3"),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)


m1<-de_ws5a_cp2a_ipsc_up_mito
m1<-de_ws5a_cp2a_ipsc_down_mito
m1<-de_ws5a_cp2a_nsc_up_mito
m1<-de_ws5a_cp2a_nsc_down_mito

a1<-na.omit(de_ws5a_cp2a_list[['ipsc_up']]$EnsemblID)
a2<-na.omit(de_ws5a_cp2a_list[['ipsc_down']]$EnsemblID)
a3<-na.omit(de_ws5a_cp2a_list[['nsc_up']]$EnsemblID)
a4<-na.omit(de_ws5a_cp2a_list[['nsc_down']]$EnsemblID)


venn.diagram(
  x = list(m1, human_mito_carta$Ensembl),
  category.names = c("mito_up" , "all_up"),
  filename = 'ipsc',
  output=TRUE
)



## enrichment of up and down regulated mitochondrial related genes

#first change ID to entrze

de_ws5a_cp2a_ipsc_up_mito<-as.character(intersect(de_ws5a_cp2a_list[['ipsc_up']]$ENTREZID,human_mito_carta$ENTREZID))

de_ws5a_cp2a_ipsc_down_mito<-as.character(intersect(de_ws5a_cp2a_list[['ipsc_down']]$ENTREZID,human_mito_carta$ENTREZID))

de_ws5a_cp2a_nsc_up_mito<-as.character(intersect(de_ws5a_cp2a_list[['nsc_up']]$ENTREZID,human_mito_carta$ENTREZID))

de_ws5a_cp2a_nsc_down_mito<-as.character(intersect(de_ws5a_cp2a_list[['nsc_down']]$ENTREZID,human_mito_carta$ENTREZID))




enrich <- enrichMKEGG(de_ws5a_cp2a_nsc_up_mito, organism = "hsa", minGSSize=5,
                      universe = bg_genes,qvalueCutoff = 0.05, pvalueCutoff = 0.05,maxGSSize = 500)
enrich<-setReadable(enrich, org.Hs.eg.db, 
                    keyType = "ENTREZID")
View(enrich@result)
write.csv(enrich@result,'CP2A_vs_WS5A_down_regulat_MKEGG.csv')

?write.csv
heatplot(enrich,showCategory = 9)


enrich2 <- enrichMKEGG(de_ws5a_cp2a_nsc_down_mito, organism = "hsa", minGSSize=5,
                      universe = bg_genes,qvalueCutoff = 0.05, pvalueCutoff = 0.05,maxGSSize = 500)
enrich2<-setReadable(enrich2, org.Hs.eg.db, 
                    keyType = "ENTREZID")
heatplot(enrich2,showCategory = 50)

heatplot(enrich)
heatplot(enrich2)

heatplot(enrich,showCategory = 63)

enrich3 <- enrichMKEGG(de_ws5a_cp2a_ipsc_up_mito, organism = "hsa", minGSSize=5,
                      universe = bg_genes,qvalueCutoff = 0.05, pvalueCutoff = 0.05,maxGSSize = 500)
enrich3<-setReadable(enrich3, org.Hs.eg.db, 
                     keyType = "ENTREZID")
heatplot(enrich3,showCategory = 15)


enrich4 <- enrichMKEGG(de_ws5a_cp2a_ipsc_down_mito, organism = "hsa", minGSSize=5,
                       universe = bg_genes,qvalueCutoff = 0.05, pvalueCutoff = 0.05,maxGSSize = 500)
enrich4<-setReadable(enrich4, org.Hs.eg.db, 
                     keyType = "ENTREZID")
heatplot(enrich4)

a<-enricher(intersect(de_ws5a_cp2a_list[['nsc_up']]$ENTREZID,human_mito_carta$ENTREZID),
         pvalueCutoff = 0.05,pAdjustMethod = "BH",
                   qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 200,
                   TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

b<-enricher(intersect(de_ws5a_cp2a_list[['nsc_down']]$ENTREZID,human_mito_carta$ENTREZID), 
         pvalueCutoff = 0.05,pAdjustMethod = "BH",
         qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 200,
         TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
heatplot(b)
heatplot(a)



c<-enricher(intersect(de_ws5a_cp2a_list[['ipsc_up']]$ENTREZID,human_mito_carta$ENTREZID),
            pvalueCutoff = 1,pAdjustMethod = "BH",
            qvalueCutoff = 1, minGSSize = 1, maxGSSize = 200,
            TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

d<-enricher(intersect(de_ws5a_cp2a_list[['ipsc_down']]$ENTREZID,human_mito_carta$ENTREZID), 
            pvalueCutoff = 1,pAdjustMethod = "BH",
            qvalueCutoff = 1, minGSSize = 1, maxGSSize = 200,
            TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
heatplot(c,showCategory = 75)
heatplot(d)




### ipsc, NSC  markers

#ensembl_length<-getGeneLengthAndGCContent(ensembl_list, "hsa")
#write.csv(ensembl_length,'ensemble_gene_length.csv')

ensembl_length<-read.csv('/home/yu/Postdoc_project/NSC_RNAseq/RNA_seq/metabolic_analysis/ensemble_gene_length.csv')

cpmMatrixFiltered_unlog <- data.frame(gene_id=as.character(countMatrixFiltered$gene_id), Count2CPM(countMatrixFiltered[,-1]))

cpm_matrix_Ens<-left_join(cpmMatrixFiltered_unlog, geneNames %>% dplyr::select(gene_id, gene_name)) %>% 
  mutate(gene_id = gsub('[.][0-9]+$','',.$gene_id))
# 
cpm_matrix_Ense<-left_join(ensembl_length[,c(1,2)],cpm_matrix_Ens,by=c('X'='gene_id'))

#change count of fragments to RFKM: RFKM (log transformed)= log2(CPM*1000/gene_length  + 1)

ExpData_FPKM<-data.frame(ensemblID=cpm_matrix_Ense$X,GeneSymbol=cpm_matrix_Ense$gene_name,
                         log2(cpm_matrix_Ense[,c(3:(length(colnames(cpm_matrix_Ense))-1))]*1000/cpm_matrix_Ense[,2]+1))

ExpData_FPKM_nonlog<-data.frame(ensemblID=cpm_matrix_Ense$X,GeneSymbol=cpm_matrix_Ense$gene_name,
                         cpm_matrix_Ense[,c(3:(length(colnames(cpm_matrix_Ense))-1))]*1000/cpm_matrix_Ense[,2])


write.csv(ExpData_FPKM,'RFKM.csv',row.names=F)

## NSC or IPS sample ID

NSC_ctrl_WS5A_id <- Metadata %>% filter(cell_type == 'NSC', subname_mut %in% c('WS5A')) %>%
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)
NSC_ctrl_CP2A_id <- Metadata %>% filter(cell_type == 'NSC', subname_mut %in% c('CP2A')) %>%
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)
IPS_ctrl_WS5A_id <- Metadata %>% filter(cell_type == 'IPSC', subname_mut %in% c('WS5A')) %>%
  filter(!sample_outlier,  !possible_outlier) %>% pull(1) 
IPS_ctrl_CP2A_id <- Metadata %>% filter(cell_type == 'IPSC', subname_mut %in% c('CP2A')) %>%
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)

NSC_ctrl_WS5A_id <- Metadata %>% filter(cell_type == 'NSC', subname_mut %in% c('WS5A')) %>%
  filter(sample_outlier |  possible_outlier) %>% pull(1)
NSC_ctrl_CP2A_id <- Metadata %>% filter(cell_type == 'NSC', subname_mut %in% c('CP2A')) %>%
  filter(sample_outlier |  possible_outlier) %>% pull(1)
IPS_ctrl_WS5A_id <- Metadata %>% filter(cell_type == 'IPSC', subname_mut %in% c('WS5A')) %>%
  filter(sample_outlier |  possible_outlier) %>% pull(1) 
IPS_ctrl_CP2A_id <- Metadata %>% filter(cell_type == 'IPSC', subname_mut %in% c('CP2A')) %>%
  filter(sample_outlier | possible_outlier) %>% pull(1)



## astrocyte data

astrocyte_marker<-read.csv('/home/yu/Postdoc_project/NSC_RNAseq/WS5A_CP2A_paper/WS5A_CP2A_astrocyte_marker.txt',sep='\t')

astrocyte_gene_expression<-read.csv('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/Yu_Data/astrocyte_ws5a_cp2a.csv',sep=';', stringsAsFactors = F)

astrocyte_gene_expression[astrocyte_gene_expression$GeneSymbol == 'GFAP',]

astrocyte_rpkm<-data.frame(ensemblID = as.character(astrocyte_gene_expression[,4]),
                           GeneSymbol = astrocyte_gene_expression[,2],
                           log2(astrocyte_gene_expression[,c(9:23)]+1))

astrocyte_rpkm %>% filter (GeneSymbol %in% 'NANOG' ) %>% gather('sample','rpkm',3:length(.)) %>% 
  mutate(Group=substr(.$sample,1,7),CellType='Astrocyte') %>% 
  dplyr::select('ensemblID','GeneSymbol','Group','sample','rpkm','CellType') 

## main function to give cell type specific markers of specified gene
## all genes: levels=c('NANOG','POU5F1','SOX2','PAX6','MAP2','CDH2','SOX9','NFIX','ALDH1A1','GJA1','SLC1A3','GFAP')

Marker_Expr<-function(genename,celltype) {
  
  if (celltype == 'astrocyte') {
    astro_fpkm<-astrocyte_rpkm %>% filter (GeneSymbol %in% genename ) %>% gather('sample','rpkm',3:length(.)) %>% 
      mutate(Group=substr(.$sample,1,4),CellType='Astrocyte') %>% 
      dplyr::select('ensemblID','GeneSymbol','Group','sample','rpkm','CellType') 
    return(astro_fpkm)
  } 
  
  else if (celltype == 'NSC') {
    NSC_fpkm_WS5A<-ExpData_FPKM %>% mutate(Group= 'WS5A') %>% dplyr::select(ensemblID,GeneSymbol,Group, NSC_ctrl_WS5A_id) %>%
      filter(GeneSymbol %in% genename)  %>%
      gather('sample','rpkm',4:length(colnames(.)))
    
    NSC_fpkm_CP2A<-ExpData_FPKM %>% mutate(Group= 'CP2A') %>% dplyr::select(ensemblID,GeneSymbol,Group, NSC_ctrl_CP2A_id) %>%
      filter(GeneSymbol %in% genename)  %>%
      gather('sample','rpkm',4:length(colnames(.)))
    
    NSC_fpkm_merge<-rbind(NSC_fpkm_WS5A,NSC_fpkm_CP2A) %>% 
      mutate(GeneSymbol=factor(GeneSymbol,levels=genename),CellType='NSC') 
    
    return(NSC_fpkm_merge)
  }
  
  else if (celltype == 'IPS') {
    IPS_fpkm_WS5A<-ExpData_FPKM %>% mutate(Group= 'WS5A') %>% dplyr::select(ensemblID,GeneSymbol,Group, IPS_ctrl_WS5A_id) %>%
      filter(GeneSymbol %in% genename)  %>%
      gather('sample','rpkm',4:length(colnames(.)))
    
    IPS_fpkm_CP2A<-ExpData_FPKM %>% mutate(Group= 'CP2A') %>% dplyr::select(ensemblID,GeneSymbol,Group, IPS_ctrl_CP2A_id) %>%
      filter(GeneSymbol %in% genename)  %>%
      gather('sample','rpkm',4:length(colnames(.)))
    
    IPS_fpkm_merge<-rbind(IPS_fpkm_WS5A,IPS_fpkm_CP2A) %>% 
      mutate(GeneSymbol=factor(GeneSymbol,levels=genename),CellType='IPS') 
    
    return(IPS_fpkm_merge)
  }
  
  else {
    print('cell type shoule be one of the following: astrocyte,NSC,IPS.')
  }
  
}

Marker_Expr(c('NANOG'),'astrocyte')
Marker_Expr(c('NANOG'),'NSC')

plot_marker<-function(genename) {
  #combine gene expression of each cell type
  Marker_all_celltype<-rbind(Marker_Expr(genename,'astrocyte'),
                             Marker_Expr(genename,'NSC'),
                             Marker_Expr(genename,'IPS')) %>% 
    mutate(CellType= factor(CellType,levels=c('IPS','NSC','Astrocyte')))
  
  ## calculate the mean and sd of rpkm of each cell type in WS5A and CP2A respectively
  my_sum<-Marker_all_celltype %>%
    group_by(Group,CellType) %>%
    summarise( 
      n=n(),
      mean=mean(rpkm),
      sd=sd(rpkm)
    ) %>%
    mutate( se=sd/sqrt(n))  %>%
    mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
  
  my_sum_point <- Marker_all_celltype %>%
    group_by(Group,CellType)
  
  ## draw barplots
  
  ggplot(data=my_sum, aes(y=mean, x=CellType,fill=Group))+
    geom_bar(position="dodge", stat="identity",colour="black", alpha=1) + 
    #scale_fill_manual(values=c('#999999','#E69F00')) + 
    geom_errorbar(data=my_sum, aes(x=CellType, ymin=mean-se, ymax=mean+se), 
                  width=0.6, colour="black", alpha=1, size=0.6,
                  position=position_dodge(0.9))+
    theme_classic() + 
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5), legend.position = "none") + ggtitle(genename)
  
}


plot_marker_speci<-function(genename,cell_type) {
  ## calculate p values
  ## calculate the mean and sd of rpkm of each cell type in WS5A and CP2A respectively
  my_sum<-Marker_Expr(genename,cell_type) %>%
    group_by(Group) %>%
    summarise( 
      n=n(),
      mean=mean(rpkm),
      sd=sd(rpkm)
    ) %>%
    mutate( se=sd/sqrt(n))  %>%
    mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
  
  my_sum_point <- Marker_Expr(genename,cell_type) %>%
    group_by(Group)
  
  ## draw barplots
  a<-ggplot(data=my_sum, aes(y=mean, x=Group))+
    geom_bar(position="dodge", stat="identity",colour="black", fill=c('#588BAE','#E5FFCC'), alpha=1, width = 0.8) + 
    geom_errorbar(data=my_sum, aes(x=Group, ymin=mean-se, ymax=mean+se), 
                  width=0.6, colour="black", alpha=1, size=0.6,
                  position=position_dodge(0.9))+
    ylim(0, max(my_sum$mean+my_sum$sd)*1.03)+
    labs(y = "FPKM (log)") +
    theme_classic() + 
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5,size = 10,margin=margin(0,0,3,0)), legend.position = "none") + 
    ggtitle(genename) 
  
  
  ## add p values
  stat.test <- compare_means (rpkm~Group, data = Marker_Expr(genename,cell_type)) %>% 
    mutate(y.position = max(my_sum$mean+my_sum$sd))
  
  a + stat_pvalue_manual(stat.test,label = "p.signif") 
}


plot_marker_speci('NANOG','IPS')

## all genes: levels=c('NANOG','POU5F1','SOX2','PAX6','MAP2','CDH2','SOX9','NFIX','ALDH1A1','GJA1','SLC1A3','GFAP')

## IPS ('NANOG','POU5F1','SOX2','PAX6','MAP2','CDH2','SOX9','NFIX','ALDH1A1','SLC1A3','GFAP')

## NSC ('NANOG','POU5F1','SOX2','PAX6','MAP2','CDH2','SOX9','NFIX','ALDH1A1','SLC1A3','GFAP')

## astrocytes c('NANOG','PAX6','MAP2','CDH2','SOX9','NFIX','ALDH1A1','GJA1','SLC1A3','GFAP')

all_genes<-c('NANOG','POU5F1','SOX2','PAX6','MAP2','CDH2','SOX9','NFIX','ALDH1A1','GJA1','SLC1A3','GFAP')

plot_marker('NANOG')
plot_marker('SOX2')
plot_marker('POU5F1') # oct4
plot_marker('PODXL')
plot_marker('FUT4') # SSEA1

plot_marker('ISL1')
plot_marker('FABP7')

plot_marker('HES1')
plot_marker('SOX2')
plot_marker('PAX6')
plot_marker('NES')
plot_marker('SOX9')

plot_marker('MAP2')
plot_marker('CDH2')

plot_marker('ALDH1A1')
plot_marker('GJA1')
plot_marker('GFAP')
plot_marker('SLC1A3') # EAAT1
plot_marker('GLUL') # Glutamine synthetase
plot_marker('CD44')
plot_marker('NFIX')


library("ggpubr")
library(gtable)
library(gridExtra)
#IPSC
grid.arrange(plot_marker_speci('NANOG','IPS'),
             plot_marker_speci('POU5F1','IPS'),
             plot_marker_speci('SOX2','IPS'),
             plot_marker_speci('PODXL','IPS'),
             plot_marker_speci('FUT4','IPS'),
             ncol = 2)

wilcox.test(Marker_Expr('NANOG','IPS')$rpkm ~ Marker_Expr('NANOG','IPS')$Group)
wilcox.test(Marker_Expr('POU5F1','IPS')$rpkm ~ Marker_Expr('POU5F1','IPS')$Group)
wilcox.test(Marker_Expr('SOX2','IPS')$rpkm ~ Marker_Expr('SOX2','IPS')$Group)
wilcox.test(Marker_Expr('PODXL','IPS')$rpkm ~ Marker_Expr('PODXL','IPS')$Group)
wilcox.test(Marker_Expr('FUT4','IPS')$rpkm ~ Marker_Expr('FUT4','IPS')$Group)

# NSC
grid.arrange(plot_marker_speci('SOX9','NSC'),
             plot_marker_speci('SOX2','NSC'),
             plot_marker_speci('NES','NSC'),
             plot_marker_speci('PAX6','NSC'),
             plot_marker_speci('HES1','NSC'),
             ncol = 2)

wilcox.test(Marker_Expr('SOX9','NSC')$rpkm ~ Marker_Expr('SOX9','NSC')$Group)
wilcox.test(Marker_Expr('SOX2','NSC')$rpkm ~ Marker_Expr('SOX2','NSC')$Group)
wilcox.test(Marker_Expr('NES','NSC')$rpkm ~ Marker_Expr('NES','NSC')$Group)
wilcox.test(Marker_Expr('PAX6','NSC')$rpkm ~ Marker_Expr('PAX6','NSC')$Group)
wilcox.test(Marker_Expr('ISL1','NSC')$rpkm ~ Marker_Expr('ISL1','NSC')$Group)
wilcox.test(Marker_Expr('HES1','NSC')$rpkm ~ Marker_Expr('HES1','NSC')$Group)
wilcox.test(Marker_Expr('FABP7','NSC')$rpkm ~ Marker_Expr('FABP7','NSC')$Group)



NSC.list<-c('SOX2','SOX9','NES','PAX6','ISL1','HES1','FABP7')



res_ws5a_cp2a_nsc[res_ws5a_cp2a_nsc$GeneSymbol %in% NSC.list, ]


## FPKM, CPM and count of specific gene

# ExpData_FPKM_nonlog[ExpData_FPKM_nonlog$GeneSymbol == 'PAX6',]
# cpm_matrix_Ense[cpm_matrix_Ense$gene_name == 'PAX6',]
# countMatrixFiltered %>% mutate(gene_id = gsub('[.][0-9]+$','',.$gene_id)) %>% filter (gene_id == 'ENSG00000007372')

#astrocyte
grid.arrange(plot_marker_speci('CD44','astrocyte'),
             plot_marker_speci('NFIX','astrocyte'),
             plot_marker_speci('GLUL','astrocyte'),
             plot_marker_speci('GJA1','astrocyte'),
             plot_marker_speci('GFAP','astrocyte'),
             ncol = 2)

grid.arrange(plot_marker_speci('CD44','astrocyte'),
             plot_marker_speci('NFIX','astrocyte'),
             ncol = 2)

wilcox.test(Marker_Expr('CD44','astrocyte')$rpkm ~ Marker_Expr('CD44','astrocyte')$Group)
wilcox.test(Marker_Expr('NFIX','astrocyte')$rpkm ~ Marker_Expr('NFIX','astrocyte')$Group)
wilcox.test(Marker_Expr('GLUL','astrocyte')$rpkm ~ Marker_Expr('GLUL','astrocyte')$Group)
wilcox.test(Marker_Expr('GJA1','astrocyte')$rpkm ~ Marker_Expr('GJA1','astrocyte')$Group)
wilcox.test(Marker_Expr('GFAP','astrocyte')$rpkm ~ Marker_Expr('GFAP','astrocyte')$Group)


## plot function
plot_marker2<-function(x,y) 
{
  my_sum<-x %>%
  group_by(Group,GeneSymbol) %>%
  summarise( 
    n=n(),
    mean=mean(rpkm),
    sd=sd(rpkm)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

my_sum_point <- x %>%
  group_by(Group,GeneSymbol)


ggplot(data=my_sum, aes(y=mean, x=GeneSymbol,fill=Group))+
  geom_bar(position="dodge", stat="identity",colour="black", alpha=1) + 
  scale_fill_manual(values=c('#999999','#E69F00')) + 
  geom_errorbar(data=my_sum, aes(x=GeneSymbol, ymin=mean-se, ymax=mean+se), 
                width=0.6, colour="black", alpha=1, size=0.6,
                position=position_dodge(0.9))+
  theme_classic() + 
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5)) + ggtitle(y)
  
}
all_markers<-c('NANOG','POU5F1','SOX2','PAX6','MAP2','CDH2','SOX9','NFIX','ALDH1A1','GJA1','SLC1A3','GFAP','AQP4','S100B')

IPS<-c('NANOG','POU5F1','SOX2','PAX6','MAP2','CDH2','SOX9','NFIX','ALDH1A1','GFAP','AQP4','S100B')

NSC<-c('NANOG','POU5F1','SOX2','PAX6','MAP2','CDH2','SOX9','NFIX','ALDH1A1','GFAP','AQP4','S100B')

Astrocyte<-c('NANOG','PAX6','CDH2','SOX9','NFIX','ALDH1A1','GJA1','SLC1A3','GFAP','AQP4')


plot_marker2(Marker_Expr(IPS,'IPS'),'IPS')
plot_marker2(Marker_Expr(NSC,'NSC'),'NSC')
plot_marker2(Marker_Expr(Astrocyte,'astrocyte'),'astrocyte')
