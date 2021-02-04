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



ExpDataCPM


# import the target gene list to be analyzed


#import DE results, import function
import_de<-function(a) {
  de_import<-read.table(a, sep= '\t',header =T,stringsAsFactors = F) %>%
    mutate(EnsemblID = gsub('[.][0-9]+$','',.$EnsemblID))%>% 
    dplyr::select(EnsemblID,log2FoldChange,padj,pvalue) %>% 
    left_join(en2ENSE[,1:2],by=(c('EnsemblID'='ENSEMBL')))
  return(de_import)
}


#CP2A NSC vs control NSC DE results
cp2a<-import_de('Yu_tables/DE_NSC-CP2A_vs_NSC-CTRL-all.tsv')
cp2a_IPS<-import_de('Yu_tables/DE_IPS-CP2A_vs_IPS-CTRL-all.tsv')

#WS5A NSC vs control NSC DE results
ws5a<-import_de('Yu_tables/DE_NSC-WS5A_vs_NSC-CTRL-all.tsv')
ws5a_IPS<-import_de('Yu_tables/DE_IPS-WS5A_vs_IPS-CTRL-all.tsv')
#alpers up and down regulating genes
alpers<-import_de('Yu_tables/DE_NSC-ALPERSvs_NSC-CTRL-all.tsv')
alpers_IPS<-import_de('Yu_tables/DE_IPS-ALPERSvs_IPS-CTRL-all.tsv')




# make a list for each mutation, included the DE genes, up and down regulated respectively
de_gene_list<-list(cp2a_up = cp2a[,c(1,5)][cp2a$padj<0.05 & cp2a$log2FoldChange>0, ],
                   cp2a_down = cp2a[,c(1,5)][cp2a$padj<0.05 &  cp2a$log2FoldChange<0, ],
                   ws5a_up = ws5a[,c(1,5)][ws5a$padj<0.05 &  ws5a$log2FoldChange>0, ],
                   ws5a_down = ws5a[,c(1,5)][ws5a$padj<0.05 &  ws5a$log2FoldChange<0, ],
                   alpers_up = alpers[,c(1,5)][alpers$padj<0.05 &  alpers$log2FoldChange>0, ],
                   alpers_down = alpers[,c(1,5)][alpers$padj<0.05 &  alpers$log2FoldChange<0, ],
                   cp2a_IPS_up = cp2a_IPS[,c(1,5)][cp2a_IPS$padj<0.05 & cp2a_IPS$log2FoldChange>0, ],
                   cp2a_IPS_down = cp2a_IPS[,c(1,5)][cp2a_IPS$padj<0.05 &  cp2a_IPS$log2FoldChange<0, ],
                   ws5a_IPS_up = ws5a_IPS[,c(1,5)][ws5a_IPS$padj<0.05 &  ws5a_IPS$log2FoldChange>0, ],
                   ws5a_IPS_down = ws5a_IPS[,c(1,5)][ws5a_IPS$padj<0.05 &  ws5a_IPS$log2FoldChange<0, ],
                   alpers_IPS_up = alpers_IPS[,c(1,5)][alpers_IPS$padj<0.05 &  alpers_IPS$log2FoldChange>0, ],
                   alpers_IPS_down = alpers_IPS[,c(1,5)][alpers_IPS$padj<0.05 &  alpers_IPS$log2FoldChange<0, ])



# Make our own gene index
# remove ensembl version number(dot) in ensemble ID and add ENTREZID gene ID from our own gene expression list
Gene_ID_ENTREZ<- ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  left_join(en2ENSE[,1:2],by=c('ensemblID'='ENSEMBL')) %>% dplyr::select(GeneSymbol,ENTREZID)



bg_genes<-Gene_ID_ENTREZ %>% dplyr::select(2) %>% pull() %>% unique()
bg_genes_ense<-ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  dplyr::select(ensemblID) %>% pull() %>% unique() 

#GSEA gene list ENTREZID
GSEA_list<-function(a) {
  gene_expression<- a %>% filter(!is.na(ENTREZID))%>% 
    mutate(ab_logFC = abs(log2FoldChange)) %>% 
    group_by(ENTREZID)%>% 
    filter(ab_logFC == max(ab_logFC)) %>%
    dplyr::select(ENTREZID, log2FoldChange) %>%
    arrange(desc(log2FoldChange)) %>% data.frame(.)
  gene_expression_2<-gene_expression[,2]
  names(gene_expression_2)<-gene_expression[,1]
  return(gene_expression_2)
}



Fc_combine<-data.frame(ID=names(GSEA_list(cp2a)),CP2A=GSEA_list(cp2a)) %>% 
  inner_join(data.frame(ID=names(GSEA_list(alpers)),alpers=GSEA_list(alpers)),by='ID') %>% 
  inner_join(data.frame(ID=names(GSEA_list(ws5a)),WS5A=GSEA_list(ws5a)),by='ID')
  
rownames(Fc_combine)<-Fc_combine[,1]  
Fc_combine<-Fc_combine[,2:4]


View(Fc_combine)
?cbind
##import wiki pathway annotation file 

wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%") %>% 
  filter(gene %in% Gene_ID_ENTREZ$ENTREZID)

wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME


## wiki enrichment analysis
wiki_enrich<-function(a) {
  enrich <- enricher(de_gene_list[[a]]$ENTREZID, pvalueCutoff = 0.05,pAdjustMethod = "BH",
                     qvalueCutoff = 0.05, minGSSize = 5, maxGSSize = 500,
                     TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
  enrich <- setReadable(enrich, org.Hs.eg.db, 
                                   keyType = "ENTREZID")
  return(enrich)
}

?enricher
?enricher
heatplot(wiki_enrich('cp2a_up'))
heatplot(wiki_enrich('cp2a_down'))
heatplot(wiki_enrich('ws5a_up'))
heatplot(wiki_enrich('ws5a_down'))
heatplot(wiki_enrich('alpers_down'))
heatplot(wiki_enrich('alpers_up'))


heatplot(wiki_enrich('cp2a_IPS_up'))
heatplot(wiki_enrich('cp2a_IPS_down'))
heatplot(wiki_enrich('ws5a_IPS_up'))
heatplot(wiki_enrich('ws5a_IPS_down'))
heatplot(wiki_enrich('alpers_IPS_down'))
heatplot(wiki_enrich('alpers_IPS_up'))



##wiki GSEA analysis
## 1 remove na entrezid ids, and
## 2 for duplicate in entrezid ids, 
##   keep the one with highest or lowest logFC values 'max(abs(logFC))'
## 3 ENTREZID id as columns
## 4 decrease the order of logFC

?GSEA
wiki_gsea<-function(a){
 wiki_GSEA <- GSEA(GSEA_list(a), TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)
 wiki_GSEA <- setReadable(wiki_GSEA, org.Hs.eg.db, keyType = "ENTREZID")
 return(wiki_GSEA)
}

wiki_gsea(ws5a)
wiki_gsea(cp2a)
wiki_gsea(alpers)

## Cell marker enrichment analysis
cell_markers <- read.table('Data/Human_cell_markers.txt', sep='\t',header=T,stringsAsFactors = F) %>%
  as_tibble(.) %>% 
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ','))

ws5a_nsc_ips_enricher<-import_de('Yu_tables/DE_NSC-WS5A_vs_IPS-WS5A.tsv') %>% 
  filter(padj<0.05 & log2FoldChange>0) %>% dplyr::select(4) %>% pull() %>%
  enricher(TERM2GENE=cell_markers, minGSSize=1, universe = bg_genes ) %>%
  setReadable(org.Hs.eg.db, keyType = "ENTREZID")

ws5a_nsc_ips_enricher

## KEGG pathway analysis
 #enrichment analysis
KEGG_enrich<-function(a) {
  enrich <- enrichKEGG(de_gene_list[[a]]$ENTREZID, organism = "hsa", keyType = "kegg",
                     pvalueCutoff = 0.05, pAdjustMethod = "BH",qvalueCutoff = 0.05,
                     universe = bg_genes)
  enrich <- setReadable(enrich, org.Hs.eg.db, 
                        keyType = "ENTREZID")
  return(enrich)
}

?enrichKEGG
cp2a_up_kegg<-KEGG_enrich('cp2a_up')
cp2a_down_kegg<-KEGG_enrich('cp2a_down')

alpers_up_kegg<-KEGG_enrich('alpers_up')
alpers_down_kegg<-KEGG_enrich('alpers_down')

# CP2A down regulated
KEGG_ID<-cp2a_down_kegg@result$Description[cp2a_down_kegg@result$p.adjust<0.05]
ggplot(cp2a_down_kegg@result[cp2a_down_kegg@result$p.adjust<0.05,],
       aes(x=Description,y=Count, fill = p.adjust,width=0.7)) +
  geom_bar(stat="identity",colour="black") + 
  #scale_fill_gradientn(colours = wes_palette("Rushmore1", 100, type = "continuous"))+
  #scale_fill_gradientn(colours = rev(wes_palette("Darjeeling1", 100, type = "continuous")))+
  scale_fill_continuous(low="#659EC7", high="white") +
  scale_x_discrete(limits = rev(KEGG_ID)) +
  theme(axis.text.x = element_text(color = "black", size = 11, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 11, angle = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 11, angle = 90, face = "plain")) +
  coord_flip()

# CP2A up regulated
KEGG_ID<-cp2a_up_kegg@result$Description[cp2a_up_kegg@result$p.adjust<0.05]
ggplot(cp2a_up_kegg@result[cp2a_down_kegg@result$p.adjust<0.05,],
       aes(x=Description,y=Count, fill = p.adjust,width=0.7)) +
  geom_bar(stat="identity",colour="black") + 
  #scale_fill_gradientn(colours = wes_palette("Rushmore1", 100, type = "continuous"))+
  #scale_fill_gradientn(colours = rev(wes_palette("Darjeeling1", 100, type = "continuous")))+
  scale_fill_continuous(low="#F75D59", high="white") +
  scale_x_discrete(limits = rev(KEGG_ID)) +
  theme(axis.text.x = element_text(color = "black", size = 11, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 11, angle = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 11, angle = 90, face = "plain")) +
  coord_flip()


# alpers down regulated
KEGG_ID<-alpers_down_kegg@result$Description[alpers_down_kegg@result$p.adjust<0.05]
ggplot(alpers_down_kegg@result[alpers_down_kegg@result$p.adjust<0.05,],
       aes(x=Description,y=Count, fill = p.adjust,width=0.7)) +
  geom_bar(stat="identity",colour="black") + 
  #scale_fill_gradientn(colours = wes_palette("Rushmore1", 100, type = "continuous"))+
  #scale_fill_gradientn(colours = rev(wes_palette("Darjeeling1", 100, type = "continuous")))+
  scale_fill_continuous(low="#659EC7", high="white") +
  scale_x_discrete(limits = rev(KEGG_ID)) +
  theme(axis.text.x = element_text(color = "black", size = 11, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 11, angle = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 11, angle = 90, face = "plain")) +
  coord_flip()

# alpers up regulated
KEGG_ID<-alpers_up_kegg@result$Description[alpers_up_kegg@result$p.adjust<0.05]
ggplot(alpers_up_kegg@result[alpers_up_kegg@result$p.adjust<0.05,],
       aes(x=Description,y=Count, fill = p.adjust,width=0.7)) +
  geom_bar(stat="identity",colour="black") + 
  #scale_fill_gradientn(colours = wes_palette("Rushmore1", 100, type = "continuous"))+
  #scale_fill_gradientn(colours = rev(wes_palette("Darjeeling1", 100, type = "continuous")))+
  scale_fill_continuous(low="#F75D59", high="white") +
  scale_x_discrete(limits = rev(KEGG_ID)) +
  theme(axis.text.x = element_text(color = "black", size = 11, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 11, angle = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 11, angle = 90, face = "plain")) +
  coord_flip()

barplot(cp2a_down_kegg,showCategory = 26)
barplot(cp2a_up_kegg,showCategory = 12)

cp2a_down_kegg
cp2a_up_kegg






cp2a_ipsc_up_kegg<-KEGG_enrich('cp2a_IPS_up')

alpers_ipsc_down_kegg<-KEGG_enrich('alpers_IPS_down')


KEGG_ID<-cp2a_ipsc_up_kegg@result$Description[cp2a_ipsc_up_kegg@result$p.adjust<0.05]
ggplot(cp2a_ipsc_up_kegg@result[cp2a_ipsc_up_kegg@result$p.adjust<0.05,],
       aes(x=Description,y=Count, fill = p.adjust,width=0.7)) +
  geom_bar(stat="identity",colour="black") + 
  #scale_fill_gradientn(colours = wes_palette("Rushmore1", 100, type = "continuous"))+
  #scale_fill_gradientn(colours = rev(wes_palette("Darjeeling1", 100, type = "continuous")))+
  scale_fill_continuous(low="#F75D59", high="white") +
  scale_x_discrete(limits = rev(KEGG_ID)) +
  theme(axis.text.x = element_text(color = "black", size = 11, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 11, angle = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 11, angle = 90, face = "plain")) +
  coord_flip()



KEGG_ID<-alpers_ipsc_down_kegg@result$Description[alpers_ipsc_down_kegg@result$p.adjust<0.05]
ggplot(alpers_ipsc_down_kegg@result[alpers_ipsc_down_kegg@result$p.adjust<0.05,],
       aes(x=Description,y=Count, fill = p.adjust,width=0.7)) +
  geom_bar(stat="identity",colour="black") + 
  #scale_fill_gradientn(colours = wes_palette("Rushmore1", 100, type = "continuous"))+
  #scale_fill_gradientn(colours = rev(wes_palette("Darjeeling1", 100, type = "continuous")))+
  scale_fill_continuous(low="#659EC7", high="white") +
  scale_x_discrete(limits = rev(KEGG_ID)) +
  theme(axis.text.x = element_text(color = "black", size = 11, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 11, angle = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 11, angle = 90, face = "plain")) +
  coord_flip()


## KEGG module analysis


MKEGG_enrich<-function(a) {
  enrich <- enrichMKEGG(de_gene_list[[a]]$ENTREZID, organism = "hsa", minGSSize=1,
                        universe = bg_genes,pvalueCutoff = 0.05, pAdjustMethod = "BH",qvalueCutoff = 0.05)
  enrich <- setReadable(enrich, org.Hs.eg.db, 
                        keyType = "ENTREZID")
  return(enrich)
}
?enrichMKEGG()

View(MKEGG_enrich('cp2a_down')[1:10])

MKEGG_enrich('alpers_IPS_up')
MKEGG_enrich('alpers_IPS_down')

MKEGG_enrich('cp2a_IPS_up')
MKEGG_enrich('cp2a_IPS_down')

MKEGG_enrich('ws5a_IPS_up')
MKEGG_enrich('ws5a_IPS_down')

heatplot(MKEGG_enrich('cp2a_down'))

heatplot(MKEGG_enrich('alpers_up'))
heatplot(MKEGG_enrich('alpers_down'))

MKEGG_enrich('cp2a_up')
MKEGG_enrich('cp2a_down')

MKEGG_enrich('alpers_down')
MKEGG_enrich('alpers_up')

## Pathway analysis visualization

library("pathview")
Path_view<-function(a,b){
  Path<-pathview(gene.data  = GSEA_list(a),
         pathway.id = b,
          species    = "hsa",
          limit      = list(gene=1, cpd=1))
}


Path<-pathview(gene.data  = Fc_combine,
               pathway.id = "hsa00010", #glycosis
               species    = "hsa",
               limit      = list(gene=1, cpd=1),
               same.layer = F, kegg.native = T,
               sign.pos = demo.paths$spos[1])


Path<-pathview(gene.data  = Fc_combine,
               pathway.id = "hsa00020", #TCA cycle
               species    = "hsa",
               limit      = list(gene=1, cpd=1),
               same.layer = F, kegg.native = T,
               sign.pos = demo.paths$spos[1])

Path<-pathview(gene.data  = Fc_combine,
               pathway.id = "hsa00190", # Oxidative phosphorylation
               species    = "hsa",
               limit      = list(gene=1, cpd=1),
               same.layer = F, kegg.native = T,
               sign.pos = demo.paths$spos[1])


#glycosis hsa00010
str(Path)
Path_view(Fc_combine,"hsa00020")
Path_view(alpers,"hsa00190")


##GO analysis
GO_enrich<-function(a) {
  ego <- enrichGO(gene          = de_gene_list[[a]]$ENTREZID,
                  universe      = bg_genes,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  return(ego)
}

cp2a_down_GO_enrich<-GO_enrich('cp2a_down')
cp2a_up_GO_enrich<-GO_enrich('cp2a_up')

##xx %>% ggplot(aes(x=factor(Description,level=level_order), y=GeneRatio)) +
##  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
## coord_flip() +
##  xlab("") +
##  theme_bw()
##plot_data(plot_data)


heatplot(cp2a_down_GO_enrich,showCategory = 30)
heatplot(cp2a_up_GO_enrich)

heatplot(cp2a_up_GO_enrich)
?heatplot
ws5a_down_GO_enrich<-GO_enrich('ws5a_down')
ws5a_up_GO_enrich<-GO_enrich('ws5a_up')
heatplot(ws5a_down_GO_enrich)

alpers_down_GO_enrich<-GO_enrich('alpers_down')

View(alpers_down_GO_enrich[1:322])

heatplot(alpers_down_GO_enrich)
alpers_up_GO_enrich<-GO_enrich('alpers_up')
View(alpers_down_GO_enrich[1:410])
barplot(alpers_down_GO_enrich,showCategory = 410)
heatplot(alpers_down_GO_enrich)
heatplot(alpers_up_GO_enrich)

## Since both alpers and cp2a have mitochondria enriched gene sets in down-regulated genes,
## take the interect of both gene groups
cp2a_alpers_down_ense<-intersect(de_gene_list[['cp2a_down']]$EnsemblID,
                                 de_gene_list[['alpers_down']]$EnsemblID)
cp2a_alpers_down<-intersect(de_gene_list[['cp2a_down']]$ENTREZID,
                                 de_gene_list[['alpers_down']]$ENTREZID)
ws5a_cp2a_down<-intersect(de_gene_list[['cp2a_down']]$ENTREZID,
                          de_gene_list[['ws5a_down']]$ENTREZID)

enrich <- enrichKEGG(cp2a_alpers_down, organism = "hsa", keyType = "kegg",
                     pvalueCutoff = 0.05, pAdjustMethod = "BH",
                     universe = bg_genes)

enrich[1:16]

enrich <- enrichMKEGG(cp2a_alpers_down, organism = "hsa", minGSSize=1,
                        universe = bg_genes)
  
heatplot(enrich)

enrich <- enricher(ws5a_cp2a_down, 
                   TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
heatplot(enrich)


ego <- enrichGO(gene          = cp2a_alpers_down_ense,
                universe      = bg_genes_ense,
                keyType = "ENSEMBL",
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ws5a_cp2a_up_ense<-intersect(de_gene_list[['cp2a_up']]$EnsemblID,
                          de_gene_list[['ws5a_up']]$EnsemblID)


ego_up <- enrichGO(gene          = ws5a_cp2a_up_ense,
                universe      = bg_genes_ense,
                keyType = "ENSEMBL",
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)


View(ego_up[1:47])

heatplot(ego)
View(ego[1:366])
?enrichGO
