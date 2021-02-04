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




## enrichment analysis

# import DE analysis results
polg_ctrl_nsc<-import_de('./Yu_tables/DE_NSC-All_polg_vs_NSC-CTRL-all.tsv')
polg_ctrl_ipsc<-import_de("./Yu_tables/DE_IPSC-All_polg_vs_IPSC-CTRL-all.tsv")

##pgc1
polg_ctrl_nsc[polg_ctrl_nsc$EnsemblID == 'ENSG00000109819',]
polg_ctrl_ipsc[polg_ctrl_ipsc$EnsemblID == 'ENSG00000109819',]

##PARKN
polg_ctrl_nsc[polg_ctrl_nsc$EnsemblID == 'ENSG00000158828',]
polg_ctrl_ipsc[polg_ctrl_ipsc$EnsemblID == 'ENSG00000158828',]

## mTOR
polg_ctrl_nsc[polg_ctrl_nsc$EnsemblID == 'ENSG00000198793',]

##PINK1
polg_ctrl_nsc[polg_ctrl_ipsc$EnsemblID == 'ENSG00000185345',]
polg_ctrl_ipsc[polg_ctrl_ipsc$EnsemblID == 'ENSG00000185345',]

ctrl_nsc_ipsc<-import_de("./Yu_tables/DE_NSC-All_CTRL_vs_IPS-CTRL-all.tsv")

de_polg_ctrl_list<-list(nsc_up = polg_ctrl_nsc[,c(1,5)][polg_ctrl_nsc$padj<0.05 & polg_ctrl_nsc$log2FoldChange>0, ],
                        nsc_down = polg_ctrl_ipsc[,c(1,5)][polg_ctrl_nsc$padj<0.05 & polg_ctrl_nsc$log2FoldChange<0, ],
                        ipsc_up = polg_ctrl_ipsc[,c(1,5)][polg_ctrl_nsc$padj<0.05 & polg_ctrl_ipsc$log2FoldChange>0, ],
                        ipsc_down = polg_ctrl_ipsc[,c(1,5)][polg_ctrl_ipsc$padj<0.05 & polg_ctrl_ipsc$log2FoldChange<0, ],
                        ctrl_nsc_down = ctrl_nsc_ipsc[,c(1,5)][ctrl_nsc_ipsc$padj<0.05 & ctrl_nsc_ipsc$log2FoldChange<0, ],
                        ctrl_nsc_up = ctrl_nsc_ipsc[,c(1,5)][ctrl_nsc_ipsc$padj<0.05 & ctrl_nsc_ipsc$log2FoldChange>0, ])

length(de_polg_ctrl_list[['nsc_up']][,1])
length(de_polg_ctrl_list[['nsc_down']][,1])
length(de_polg_ctrl_list[['ipsc_up']][,1])
length(de_polg_ctrl_list[['ipsc_down']][,1])

DE_number<-data.frame(Group=c('Up','Down','Up','Down'),
                      Cell_type=c('iPSC','iPSC','NSC','NSC'),
                      DE_gene_number=c(length(de_polg_ctrl_list[['ipsc_up']][,1]),length(de_polg_ctrl_list[['ipsc_down']][,1]),
                                       length(de_polg_ctrl_list[['nsc_up']][,1]),length(de_polg_ctrl_list[['nsc_down']][,1])))
p <- ggplot(DE_number, aes(x=Cell_type, y=DE_gene_number, fill=Group)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_classic() + 
  theme(axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black"))

View(de_polg_ctrl_list[['nsc_up']])


# Make our own gene index
# remove ensembl version number(dot) in ensemble ID and add ENTREZID gene ID from our own gene expression list
Gene_ID_ENTREZ<- ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  left_join(en2ENSE[,1:2],by=c('ensemblID'='ENSEMBL')) %>% dplyr::select(GeneSymbol,ENTREZID)



bg_genes<-Gene_ID_ENTREZ %>% dplyr::select(2) %>% pull() %>% unique()
bg_genes_ense<-ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  dplyr::select(ensemblID) %>% pull() %>% unique() 

## KEGG analysis

KEGG_polg_ctrl_enrich<-function(a) {
  enrich <- enrichKEGG(de_polg_ctrl_list[[a]]$ENTREZID, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05, pAdjustMethod = "none",qvalueCutoff = 1,
                       universe = bg_genes)
  enrich <- setReadable(enrich, org.Hs.eg.db, 
                        keyType = "ENTREZID")
  return(enrich)
}


View(KEGG_polg_ctrl_enrich('nsc_down')@result)
View(KEGG_polg_ctrl_enrich('nsc_up')@result)

write.csv(KEGG_polg_ctrl_enrich('nsc_up')@result, 'nsc_up.csv')
View(KEGG_polg_ctrl_enrich('ipsc_up')@result)
View(KEGG_polg_ctrl_enrich('ipsc_down')@result)


barplot(KEGG_polg_ctrl_enrich('nsc_up'))
barplot(KEGG_polg_ctrl_enrich('nsc_down'))
## GO analysis
GO_polg_ctrl_enrich<-function(a) {
  ego <- enrichGO(gene          = de_polg_ctrl_list[[a]]$ENTREZID,
                  universe      = bg_genes,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "none",
                  minGSSize = 5,
                  maxGSSize = 150,
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)
  return(ego)
}

de_polg_ctrl_list[['nsc_down']]$ENTREZID

## KEGG module analysis
MKEGG_polg_ctrl_enrich<-function(a) {
  enrich <- enrichMKEGG(de_polg_ctrl_list[[a]]$ENTREZID, organism = "hsa", minGSSize=1,
                        universe = bg_genes, pAdjustMethod= 'none',qvalueCutoff = 1)
  enrich <- setReadable(enrich, org.Hs.eg.db, 
                        keyType = "ENTREZID")
  return(enrich)
}

?enrichMKEGG
MKEGG_polg_ctrl_enrich('nsc_down')
MKEGG_polg_ctrl_enrich('nsc_up')

View(MKEGG_polg_ctrl_enrich('nsc_down'))

##import wiki pathway annotation file 

wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%") %>% 
  filter(gene %in% Gene_ID_ENTREZ$ENTREZID)

wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME


## wiki enrichment analysis
wiki_polg_ctrl_enrich<-function(a) {
  enrich <- enricher(de_polg_ctrl_list[[a]]$ENTREZID, pvalueCutoff = 0.05,pAdjustMethod = "BH",
                     qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 200,
                     TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
  enrich <- setReadable(enrich, org.Hs.eg.db, 
                        keyType = "ENTREZID")
  return(enrich)
}

## human mito carta 2.o database

## ID table
human_mito_carta<-read.table('./human_mito_carta.txt',sep='\t',header=TRUE)
colnames(human_mito_carta)<-c('ENTREZID','Ensembl')

## full_annotation table
human_mito_carta_full<-read.csv('./Human.MitoCarta2.0.csv',sep='>',header=TRUE)


## instersect with DE up and down regulated gene list, and generate the mito_de gene information table, 
## note that de results is WS5A/CP2A, so in the table reverse the direction of log2foldchange

de_polg_ctrl_ipsc_up_mito<-as.character(intersect(de_polg_ctrl_list[['ipsc_up']]$EnsemblID,human_mito_carta$Ensembl))


polg_ctrl_ipsc_mito_up_table<-polg_ctrl_ipsc %>%  filter(EnsemblID %in% de_polg_ctrl_ipsc_up_mito) %>% 
  left_join(human_mito_carta_full[,c(4,6,7,10)],by = c("EnsemblID" = "EnsemblGeneID")) %>%
  dplyr::select(Symbol,EnsemblID,Description,MitoDomain_Score,log2FoldChange,padj)


de_polg_ctrl_ipsc_down_mito<-as.character(intersect(de_polg_ctrl_list[['ipsc_down']]$EnsemblID,human_mito_carta$Ensembl))

polg_ctrl_ipsc_mito_down_table<-polg_ctrl_ipsc %>%  filter(EnsemblID %in% de_polg_ctrl_ipsc_down_mito) %>% 
  left_join(human_mito_carta_full[,c(4,6,7,10)],by = c("EnsemblID" = "EnsemblGeneID")) %>%
  dplyr::select(Symbol,EnsemblID,Description,MitoDomain_Score,log2FoldChange,padj)

#write.csv(ws5a_cp2a_ipsc_mito_up_table,'ws5a_cp2a_ipsc_mito_up_table.csv', col.names=T)


de_polg_ctrl_nsc_up_mito<-as.character(intersect(de_polg_ctrl_list[['nsc_up']]$EnsemblID,human_mito_carta$Ensembl))

polg_ctrl_nsc_mito_up_table<-polg_ctrl_nsc %>%  filter(EnsemblID %in% de_polg_ctrl_nsc_up_mito) %>% 
  left_join(human_mito_carta_full[,c(4,6,7,10)],by = c("EnsemblID" = "EnsemblGeneID")) %>%
  dplyr::select(Symbol,EnsemblID,Description,MitoDomain_Score,log2FoldChange,padj)

#write.csv(ws5a_cp2a_nsc_mito_down_table,'ws5a_cp2a_nsc_mito_down_table.csv', col.names=T)


de_polg_ctrl_nsc_down_mito<-as.character(intersect(de_polg_ctrl_list[['nsc_down']]$EnsemblID,human_mito_carta$Ensembl))

polg_ctrl_nsc_mito_down_table<-polg_ctrl_nsc %>%  filter(EnsemblID %in% de_polg_ctrl_nsc_down_mito) %>% 
  left_join(human_mito_carta_full[,c(4,6,7,10)],by = c("EnsemblID" = "EnsemblGeneID")) %>%
  dplyr::select(Symbol,EnsemblID,Description,MitoDomain_Score,log2FoldChange,padj)

#write.csv(ws5a_cp2a_nsc_mito_up_table,'ws5a_cp2a_nsc_mito_up_table.csv', col.names=T)

## export table with mito annotation data


## barplot of up and down regulated genes and mitochondrial related genes
length(de_polg_ctrl_ipsc_up_mito)
length(de_polg_ctrl_ipsc_down_mito)
length(de_polg_ctrl_nsc_up_mito)
length(de_polg_ctrl_nsc_down_mito)

length(na.omit(de_polg_ctrl_list[['ipsc_up']]$EnsemblID))
length(na.omit(de_polg_ctrl_list[['ipsc_down']]$EnsemblID))
length(na.omit(de_polg_ctrl_list[['nsc_up']]$EnsemblID))
length(na.omit(de_polg_ctrl_list[['nsc_down']]$EnsemblID))


## enrichment of up and down regulated mitochondrial related genes

#first change ID to entrze

de_polg_ctrl_ipsc_up_mito<-as.character(intersect(de_polg_ctrl_list[['ipsc_up']]$ENTREZID,human_mito_carta$ENTREZID))

de_polg_ctrl_ipsc_down_mito<-as.character(intersect(de_polg_ctrl_list[['ipsc_down']]$ENTREZID,human_mito_carta$ENTREZID))

de_polg_ctrl_nsc_up_mito<-as.character(intersect(de_polg_ctrl_list[['nsc_up']]$ENTREZID,human_mito_carta$ENTREZID))

de_polg_ctrl_nsc_down_mito<-as.character(intersect(de_polg_ctrl_list[['nsc_down']]$ENTREZID,human_mito_carta$ENTREZID))


bg_genes_polg_ctrl<-as.character(intersect(bg_genes,human_mito_carta$ENTREZID))


enrich <- enrichMKEGG(de_polg_ctrl_nsc_up_mito, organism = "hsa", minGSSize=5,
                      universe = bg_genes,qvalueCutoff = 1, pvalueCutoff = 1,maxGSSize = 500, pAdjustMethod='none')
enrich<-setReadable(enrich, org.Hs.eg.db, 
                    keyType = "ENTREZID")


enrich2 <- enrichMKEGG(de_polg_ctrl_nsc_down_mito, organism = "hsa", minGSSize=5,
                       universe = bg_genes,qvalueCutoff = 1, pvalueCutoff = 1,maxGSSize = 500, pAdjustMethod='none')
enrich2<-setReadable(enrich2, org.Hs.eg.db, 
                     keyType = "ENTREZID")


enrich3 <- enrichMKEGG(de_polg_ctrl_ipsc_up_mito, organism = "hsa", minGSSize=5,
                       universe = bg_genes,qvalueCutoff = 1, pvalueCutoff = 0.05,maxGSSize = 500, pAdjustMethod='none')
enrich3<-setReadable(enrich3, org.Hs.eg.db, 
                     keyType = "ENTREZID")


enrich4 <- enrichMKEGG(de_polg_ctrl_ipsc_down_mito, organism = "hsa", minGSSize=5,
                       universe = bg_genes,qvalueCutoff = 1, pvalueCutoff = 0.05,maxGSSize = 500, pAdjustMethod='none')
enrich4<-setReadable(enrich4, org.Hs.eg.db, 
                     keyType = "ENTREZID")


## mito de genes wiki analysis
a<-enricher(intersect(de_polg_ctrl_list[['nsc_up']]$ENTREZID,human_mito_carta$ENTREZID),
            pvalueCutoff = 0.05,pAdjustMethod = "BH",
            qvalueCutoff = 1, minGSSize = 10, maxGSSize = 200,
            TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
a<-setReadable(a, org.Hs.eg.db, 
               keyType = "ENTREZID")

b<-enricher(intersect(de_polg_ctrl_list[['nsc_down']]$ENTREZID,human_mito_carta$ENTREZID), 
            pvalueCutoff = 0.05,pAdjustMethod = "BH",
            qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 200,
            TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
b<-setReadable(b, org.Hs.eg.db, 
               keyType = "ENTREZID")

c<-enricher(intersect(de_polg_ctrl_list[['ipsc_up']]$ENTREZID,human_mito_carta$ENTREZID),
            pvalueCutoff = 0.05,pAdjustMethod = "BH",
            qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 200,
            TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
c<-setReadable(c, org.Hs.eg.db, 
               keyType = "ENTREZID")

d<-enricher(intersect(de_polg_ctrl_list[['ipsc_down']]$ENTREZID,human_mito_carta$ENTREZID), 
            pvalueCutoff = 0.05,pAdjustMethod = "BH",
            qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 200,
            TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
d<-setReadable(d, org.Hs.eg.db,
               keyType = "ENTREZID") 

##GSA analysis
MKEGG_gene<-read.table('./GSA_analysis/MKEGG_gene.txt',sep='\t',stringsAsFactors = F) %>% 
  mutate(ENTREZID=as.character(.$V1), MKEGG=substr(.[,2], 5, nchar(.[,2]))) %>% 
  dplyr::select(MKEGG, ENTREZID) %>% left_join(en2ENSE) 

MKEGG_ID<-read.table('./GSA_analysis/MKEGG_ID.txt',sep='\t',stringsAsFactors = F) %>% 
  dplyr::select(MKEGG=V1, MKEGG_name=V2)

MKEGG_gene_tbl<-MKEGG_gene %>% left_join(MKEGG_ID) %>% dplyr::select(MKEGG,MKEGG_name,ENSEMBL,ENTREZID,gene_name=SYMBOL)

polg_ctrl_nsc<-polg_ctrl_nsc %>% na.omit(.[,c(3,4)])
polg_ctrl_ipsc<-polg_ctrl_ipsc %>% na.omit(.[,c(3,4)])

nsc_p<-polg_ctrl_nsc$pvalue
names(nsc_p) <- polg_ctrl_nsc$EnsemblID
nsc_fc<-polg_ctrl_nsc$log2FoldChange
names(nsc_fc) <- polg_ctrl_nsc$EnsemblID

ips_p<-polg_ctrl_ipsc$pvalue
names(ips_p) <- polg_ctrl_ipsc$EnsemblID
ips_fc<-polg_ctrl_ipsc$log2FoldChange
names(ips_fc) <- polg_ctrl_ipsc$EnsemblID

library(piano)
polg_ctrl_gsc<-loadGSC(MKEGG_gene_tbl[,c(3,2)]) 
?loadGSC

# gene ID and metabolic data set 

library(igraph) 
library(circlize)
require(visNetwork)

col_fun_2 = colorRamp2(c(-100, 0, 100), c("#CC0000", "white","#005555"))

col<-c(col_fun_2(seq(-100, 100))[1:100],rev(col_fun_2(seq(-100, 100))[102:201]))

gsaRes_nsc_polg_ctrl <- runGSA(nsc_p,nsc_fc,gsc=polg_ctrl_gsc, 
                           geneSetStat="reporter",
                           signifMethod="nullDist", 
                           nPerm=1000,
                           adjMethod = "fdr",
                           gsSizeLim=c(1,1000))

?runGSA

#gsaRes_nsc_polg_ctrl <- runGSA(nsc_p,nsc_fc,gsc=polg_ctrl_gsc,gsSizeLim=c(1,1000))

networkPlot(gsaRes_nsc_polg_ctrl, class="distinct", direction="both",
            significance=0.01,lay=layout_with_fr,edgeWidth = c(0, 15),overlap=1,
            scoreColors=col)
gsaRes_nsc_polg_ctrl_tbl <- GSAsummaryTable(gsaRes_nsc_polg_ctrl)

gsaRes_ips_polg_ctrl <- runGSA(ips_p, ips_fc,gsc=polg_ctrl_gsc, 
                               geneSetStat="reporter",
                               signifMethod="nullDist", 
                               nPerm=1000,
                               adjMethod = "fdr",
                               gsSizeLim=c(1,1000))
#gsaRes_ips_polg_ctrl <- runGSA(ips_p,ips_fc,gsc=polg_ctrl_gsc,gsSizeLim=c(1,1000))
networkPlot(gsaRes_ips_polg_ctrl, class="distinct", direction="both",
            significance=0.01,lay=layout_with_fr,edgeWidth = c(0, 15),overlap=1,
            scoreColors=col)
gsaRes_ips_polg_ctrl_tbl <- GSAsummaryTable(gsaRes_ips_polg_ctrl)

