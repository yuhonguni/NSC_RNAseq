## enrichment analysis
library(magrittr)
library(tidyverse)
library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)
library(piano)
library(enrichplot)


# 1.1 extract ensemble and entrze gene ID mapping files
k<-keys(org.Hs.eg.db,keytype='ENSEMBL')
en2ENSE<-AnnotationDbi::select(org.Hs.eg.db,keys=k,
                               columns=c("ENTREZID",'SYMBOL'),
                               keytype = 'ENSEMBL')

#1.2 Make our own gene index
# remove ensembl version number(dot) in ensemble ID and add ENTREZID gene ID from our own gene expression list
Gene_ID_ENTREZ<- ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  left_join(en2ENSE[,1:2],by=c('ensemblID'='ENSEMBL')) %>% dplyr::select(GeneSymbol,ENTREZID)

bg_genes<-Gene_ID_ENTREZ %>% dplyr::select(2) %>% pull() %>% unique()
bg_genes_ense<-ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  dplyr::select(ensemblID) %>% pull() %>% unique() 


# 1.3 import the target gene list and fold change list to be analyzed

# 1.3.1 import DE results, import function
import_de<-function(a) {
  de_import<-read.table(a, sep= '\t',header =T,stringsAsFactors = F) %>%
    mutate(EnsemblID = gsub('[.][0-9]+$','',.$EnsemblID))%>% 
    dplyr::select(EnsemblID,log2FoldChange,padj,pvalue) %>% 
    left_join(en2ENSE[,1:2],by=(c('EnsemblID'='ENSEMBL')))
  return(de_import)
}

# 1.3.2 import DE analysis results
alp_ctrl_nsc<-import_de('./Yu_tables/DE_NSC-ALPERSvs_NSC-CTRL-all.tsv')
alp_ctrl_ipsc<-import_de("./Yu_tables/DE_IPSC-ALPERSvs_IPSC-CTRL-all.tsv")




##mitphagy pathway


## 1.3.3 fold change
alp_ctrl_nsc<-alp_ctrl_nsc %>% na.omit(.[,c(3,4)])
alp_ctrl_ipsc<-alp_ctrl_ipsc %>% na.omit(.[,c(3,4)])

alp_nsc_fc<-alp_ctrl_nsc$log2FoldChange
names(alp_nsc_fc) <- alp_ctrl_nsc$ENTREZID
alp_nsc_fc = sort(alp_nsc_fc, decreasing = TRUE)

alp_ips_fc<-alp_ctrl_ipsc$log2FoldChange
names(alp_ips_fc) <- alp_ctrl_ipsc$ENTREZID
alp_ips_fc = sort(alp_ips_fc, decreasing = TRUE)

## 1.3.4 up and down regulated gene list

de_alp_ctrl_list<-list(nsc_up = alp_ctrl_nsc[,c(1,5)][alp_ctrl_nsc$padj<0.05 & alp_ctrl_nsc$log2FoldChange>0, ],
                        nsc_down = alp_ctrl_nsc[,c(1,5)][alp_ctrl_nsc$padj<0.05 & alp_ctrl_nsc$log2FoldChange<0, ],
                        ipsc_up = alp_ctrl_ipsc[,c(1,5)][alp_ctrl_ipsc$padj<0.05 & alp_ctrl_ipsc$log2FoldChange>0, ],
                        ipsc_down = alp_ctrl_ipsc[,c(1,5)][alp_ctrl_ipsc$padj<0.05 & alp_ctrl_ipsc$log2FoldChange<0, ])

length(de_alp_ctrl_list[['nsc_up']][,1])
length(de_alp_ctrl_list[['nsc_down']][,1])
length(de_alp_ctrl_list[['ipsc_up']][,1])
length(de_alp_ctrl_list[['ipsc_down']][,1])

## 1.3.5 plot number of up and down regulated genes
alp_DE_number<-data.frame(Group=c('Up','Down','Up','Down'),
                      Cell_type=c('iPSC','iPSC','NSC','NSC'),
                      DE_gene_number=c(length(de_alp_ctrl_list[['ipsc_up']][,1]),length(de_alp_ctrl_list[['ipsc_down']][,1]),
                                       length(de_alp_ctrl_list[['nsc_up']][,1]),length(de_alp_ctrl_list[['nsc_down']][,1])))
p <- ggplot(alp_DE_number, aes(x=Cell_type, y=DE_gene_number, fill=Group)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  scale_fill_manual(values=c( "blue","red")) +
  theme_minimal()+
  theme_classic() + 
  theme(axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black"))+
  labs(x="", y= "DE gene number" )

## 1.4 overexpression and enrichment and GSEA analysis

#### 1.4.1 GO analysis
GO_polg_ctrl_enrich<-function(a) {
  ego <- enrichGO(gene          = de_alp_ctrl_list[[a]]$ENTREZID,
                  universe      = bg_genes,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "fdr",
                  minGSSize = 20,
                  maxGSSize = 150,
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  return(ego)
}


alper_go_n_down <- GO_polg_ctrl_enrich('nsc_down')
alper_go_n_up <- GO_polg_ctrl_enrich('nsc_up')
alper_go_i_down <- GO_polg_ctrl_enrich('ipsc_down')
alper_go_i_up <- GO_polg_ctrl_enrich('ipsc_up')

barplot(alper_go_n_down, showCategory = 20)
barplot(alper_go_n_up, showCategory = 20)
barplot(alper_go_i_down, showCategory = 20)
barplot(alper_go_i_up, showCategory = 20)

View(alper_go_n_down@result)

## GO GSEA analysis
library(enrichplot)
nsc_gseGO <- gseGO(geneList = alp_nsc_fc,
              OrgDb        = org.Hs.eg.db,
              ont          = "ALL",
              nPerm        = 1000,
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 1,
              verbose      = FALSE,
              pAdjustMethod = "none")

dotplot(nsc_gseGO,split=".sign",showCategory = 12)+facet_grid(~.sign)

View(nsc_gseGO@result)

gseaplot2(nsc_gseGO,'GO:1903146',color="red",pvalue_table = T) ## autophage pathway
gseaplot2(nsc_gseGO,'GO:0048167',color="red",pvalue_table = T)

## mitochondrial disassembly, regulation of autophagy of mitochodrion, regulation of mitochondrial fusion
gseaplot2(nsc_gseGO,c('GO:1903146','GO:0010635','GO:0061726'),color="red",pvalue_table = F)

##mitochondrial metabolism 
gseaplot2(nsc_gseGO,c('GO:0090140','GO:1903146','GO:0010821','GO:0010635'),
          pvalue_table = T,
          subplots = 1)


View(nsc_gseGO@result)

ips_gseGO <- gseGO(geneList     = ips_fc,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "ALL",
                   nPerm        = 1000,
                   minGSSize    = 10,
                   maxGSSize    = 500,
                   pvalueCutoff = 1,
                   verbose      = FALSE,
                   pAdjustMethod = "none")

View(ips_gseGO@result)

dotplot(ips_gseGO,split=".sign",showCategory = 12)+facet_grid(~.sign)
gseaplot2(ips_gseGO,'GO:0098798',color="red",pvalue_table = F)
gseaplot2(ips_gseGO,c('GO:0098798','GO:0032981','GO:0006120','GO:1903146'),color="red",pvalue_table = T)

gseaplot2(ips_gseGO,c('GO:0090140','GO:1903146','GO:0010821','GO:0010635'),
          pvalue_table = T,
          subplots = 1# or 1:2,1:3
          )


gseaplot2(ips_gseGO,c('GO:1903146','GO:0010635','GO:0090140','GO:0010821'),color="red",pvalue_table = T)

## look into specific mitophagy and mitochondrial dynamic pathway genes

## mitochondrial dynamic pathway
alp_ctrl_nsc[alp_ctrl_nsc$ENTREZID == 55669,] #MFN1  significant up-regulate
alp_ctrl_nsc[alp_ctrl_nsc$ENTREZID == 9927,] #MFN2
alp_ctrl_nsc[alp_ctrl_nsc$ENTREZID == 10059,] #DRP1
alp_ctrl_nsc[alp_ctrl_nsc$ENTREZID == 56947,] #MFF near sifignificant Padj=0.06 down-regulate
alp_ctrl_nsc[alp_ctrl_nsc$ENTREZID == 4976,] #OPA1

alp_ctrl_ipsc[alp_ctrl_ipsc$ENTREZID == 55669,] #MFN1  near sifignificant P=0.08
alp_ctrl_ipsc[alp_ctrl_ipsc$ENTREZID == 9927,] #MFN2
alp_ctrl_ipsc[alp_ctrl_ipsc$ENTREZID == 10059,] #DRP1
alp_ctrl_ipsc[alp_ctrl_ipsc$ENTREZID == 56947,] #MFF 
alp_ctrl_ipsc[alp_ctrl_ipsc$ENTREZID == 4976,] #OPA1

## autophagy pathway 
alp_ctrl_nsc[alp_ctrl_nsc$ENTREZID == 84557,] ##LC3A MAP1LC3A significant up-regulate
alp_ctrl_nsc[alp_ctrl_nsc$ENTREZID == 81631,] ##LC3B MAP1LC3B significant up-regulate
alp_ctrl_nsc[alp_ctrl_nsc$ENTREZID == 8878,] ##SQSTM1/p62 P=0.09 up-regulate
alp_ctrl_nsc[alp_ctrl_nsc$ENTREZID == 65018 ,] ##PINK1
alp_ctrl_nsc[alp_ctrl_nsc$ENTREZID == 5071 ,] ##PAKN Padj = 0.07304728 up-regulate

alp_ctrl_nsc[alp_ctrl_nsc$ENTREZID == 9804 ,] ##TOMM20

alp_ctrl_ipsc[alp_ctrl_ipsc$ENTREZID == 84557,] ##LC3A MAP1LC3A significant up-regulate
alp_ctrl_ipsc[alp_ctrl_ipsc$ENTREZID == 81631,] ##LC3B MAP1LC3B
alp_ctrl_ipsc[alp_ctrl_ipsc$ENTREZID == 8878,] ##SQSTM1/p62 
alp_ctrl_ipsc[alp_ctrl_ipsc$ENTREZID == 65018 ,] ##PINK1
alp_ctrl_ipsc[alp_ctrl_ipsc$ENTREZID == 5071 ,] ##PAKN

alp_ctrl_ipsc[alp_ctrl_ipsc$ENTREZID == 9804 ,] ##TOMM20


## 1.4.2 KEGG analysis

KEGG_polg_ctrl_enrich<-function(a) {
  enrich <- enrichKEGG(de_polg_ctrl_list[[a]]$ENTREZID, organism = "hsa", keyType = "kegg",
                       pvalueCutoff = 0.05, pAdjustMethod = "fdr",qvalueCutoff = 0.05,
                       universe = bg_genes)
  enrich <- setReadable(enrich, org.Hs.eg.db, 
                        keyType = "ENTREZID")
  return(enrich)
}

KEG_a_n_down<-KEGG_polg_ctrl_enrich('nsc_down')
KEG_a_n_up<-KEGG_polg_ctrl_enrich('nsc_up')
KEG_a_i_down<-KEGG_polg_ctrl_enrich('ipsc_down')
KEG_a_i_up<-KEGG_polg_ctrl_enrich('ipsc_up')

## KEGG GSEA analysis

nsc_gseKEGG=gseKEGG(geneList     = nsc_fc,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 20,
                    pvalueCutoff = 0.5,
                    verbose      = FALSE,
                    pAdjustMethod = "none")

dotplot(nsc_gseKEGG,split=".sign",showCategory = 20)+facet_grid(~.sign)

View(nsc_gseKEGG@result)

gseaplot2(nsc_gseKEGG,'GO:1903146',color="red",pvalue_table = T)

ips_gseKEGG=gseKEGG(geneList     = ips_fc,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 120,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE,
                    pAdjustMethod = "none")

dotplot(ips_gseKEGG,split=".sign",showCategory = 20)+facet_grid(~.sign)


## 1.4.3 KEGG module analysis
MKEGG_polg_ctrl_enrich<-function(a) {
  enrich <- enrichMKEGG(de_polg_ctrl_list[[a]]$ENTREZID, organism = "hsa", minGSSize=1,
                        universe = bg_genes, pAdjustMethod= 'none')
  enrich <- setReadable(enrich, org.Hs.eg.db, 
                        keyType = "ENTREZID")
  return(enrich)
}


MKEGG_polg_ctrl_enrich('nsc_down')
MKEGG_polg_ctrl_enrich('nsc_up')
barplot(MKEGG_polg_ctrl_enrich('nsc_down'))
barplot(MKEGG_polg_ctrl_enrich('nsc_up'),showCategory = 15)

barplot(MKEGG_polg_ctrl_enrich('ipsc_down'))
barplot(MKEGG_polg_ctrl_enrich('ipsc_up'))

##MKEGG GSEA analysis

alp_nsc_gseMKEGG=gseMKEGG(geneList     = alp_nsc_fc,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 1,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE,
                    pAdjustMethod = "none")

View(alp_nsc_gseMKEGG@result)

dotplot(alp_nsc_gseMKEGG,split=".sign",showCategory = 10)+facet_grid(~.sign)

gseaplot2(nsc_gseMKEGG,c('M00142','M00146'),color="red",pvalue_table = T)

gseaplot2(nsc_gseMKEGG,c('M00001','M00049','M00146'),color="red",pvalue_table = T)

alp_ips_gseMKEGG=gseMKEGG(geneList     = alp_ips_fc,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 1,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE,
                    pAdjustMethod = "none")

View(alp_ips_gseMKEGG@result)

dotplot(alp_ips_gseMKEGG,split=".sign",showCategory = 4)+facet_grid(~.sign)
gseaplot2(ips_gseMKEGG,c('M00142','M00146'),color="red",pvalue_table = T)


## 1.5 GSA analysis, using piano package 
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

#gsaRes_nsc_polg_ctrl <- runGSA(nsc_p,nsc_fc,gsc=polg_ctrl_gsc,gsSizeLim=c(1,1000))

networkPlot(gsaRes_nsc_polg_ctrl, class="distinct", direction="both",
            significance=0.01,lay=layout_with_fr,edgeWidth = c(0, 15),overlap=1,
            scoreColors=col)
dev.off()
gsaRes_nsc_polg_ctrl_tbl <- GSAsummaryTable(gsaRes_nsc_polg_ctrl)

gsaRes_ips_polg_ctrl <- runGSA(ips_p, ips_fc,gsc=polg_ctrl_gsc, 
                               geneSetStat="reporter",
                               signifMethod="nullDist", 
                               nPerm=1000,
                               adjMethod = "fdr",
                               gsSizeLim=c(1,1000))
#gsaRes_ips_polg_ctrl <- runGSA(ips_p,ips_fc,gsc=polg_ctrl_gsc,gsSizeLim=c(1,1000))


networkPlot(gsaRes_ips_polg_ctrl, class="distinct", direction="both",
            significance=0.05,lay=layout_with_fr,edgeWidth = c(0, 15),overlap=1,
            scoreColors=col)


gsaRes_ips_polg_ctrl_tbl <- GSAsummaryTable(gsaRes_ips_polg_ctrl)
View(gsaRes_ips_polg_ctrl_tbl)
