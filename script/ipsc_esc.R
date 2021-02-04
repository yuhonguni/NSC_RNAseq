library(magrittr)
library(tidyverse)
library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)

#ipsc_nsc vs esc_nsc

Meta_ips_es<- Metadata %>% mutate(clonea=substr(Metadata$clone,1,nchar(Metadata$clone)-1))
  
contrast_id <- "NSC-ips_vs_NSC-es"
contrast_str <- "clonea_CTRL_vs_ESC"
Mt_ips_es_nsc <- Meta_ips_es %>% filter(cell_type %in% c('NSC'), mutation == 'CTRL') %>%
  filter(!sample_outlier,  !possible_outlier)
Mt_ips_es_nsc %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex,clonea) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

Ct_ips_es_nsc <- Counts[,colnames(Counts) %in% Mt_ips_es_nsc$sample_id]
rownames(Ct_ips_es_nsc) <- Counts$gene_id
Md_ips_es_nsc <- as.formula(paste0("~clonea"))
dds_ips_es_nsc <- DESeq2RUN(Ct_ips_es_nsc, Mt_ips_es_nsc, Md_ips_es_nsc)

res_ips_es_nsc<-results(dds_ips_es_nsc, alpha = 0.05, format = "DataFrame", independentFiltering = T)
res_ips_es_nsc$GeneSymbol <- geneNames$gene_name[match(rownames(res_ips_es_nsc), geneNames$gene_id)]
res_ips_es_nsc$EnsemblID <- rownames(res_ips_es_nsc)
res_ips_es_nsc %<>% data.frame %>% filter(GeneSymbol != "") %>% dplyr::select(GeneSymbol, EnsemblID, log2FoldChange, pvalue, padj, everything())

length(res_ips_es_nsc[res_ips_es_nsc$padj<0.05,1])

de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_ips_es_nsc, de_file)


#ips vs esc
View(Meta_ips_es)

contrast_id <- "ips_vs_es"
contrast_str <- "clonea_CTRL_vs_ESC"
Mt_ips_es <- Meta_ips_es %>% filter(cell_type %in% c('IPSC','ESC'), mutation == 'CTRL') %>%
  filter(!sample_outlier,  !possible_outlier)
Mt_ips_es %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex,clonea) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

Ct_ips_es <- Counts[,colnames(Counts) %in% Mt_ips_es$sample_id]
rownames(Ct_ips_es) <- Counts$gene_id
Md_ips_es <- as.formula(paste0("~clonea"))
dds_ips_es <- DESeq2RUN(Ct_ips_es, Mt_ips_es, Md_ips_es)

res_ips_es<-results(dds_ips_es, alpha = 0.05, format = "DataFrame", independentFiltering = T)
res_ips_es$GeneSymbol <- geneNames$gene_name[match(rownames(res_ips_es), geneNames$gene_id)]
res_ips_es$EnsemblID <- rownames(res_ips_es)
res_ips_es %<>% data.frame %>% filter(GeneSymbol != "") %>% dplyr::select(GeneSymbol, EnsemblID, log2FoldChange, pvalue, padj, everything())

length(res_ips_es[res_ips_es$padj<0.05,1])
de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_ips_es, de_file)



#import DE results, import function
import_de<-function(a) {
  de_import<-read.table(a, sep= '\t',header =T,stringsAsFactors = F) %>%
    mutate(EnsemblID = gsub('[.][0-9]+$','',.$EnsemblID))%>% 
    dplyr::select(EnsemblID,log2FoldChange,padj,pvalue) %>% 
    left_join(en2ENSE[,1:2],by=(c('EnsemblID'='ENSEMBL')))
  return(de_import)
}
###

ips_es_nsc<-import_de('./Yu_tables/DE_NSC-ips_vs_NSC-es.tsv')
ips_es<-import_de("./Yu_tables/DE_ips_vs_es.tsv")
ips_es_as<-read.table("./Yu_tables/upregulated_DE_ips_astrocyte.xls",sep='\t',header=T)

names(ips_es_as)[names(ips_es_as) == "Gene.ID"] <- "ENTREZID"
names(ips_es_as)[names(ips_es_as) == "Ensembl.Gene.ID"] <- "EnsemblID"

de_ips_es_list<-list(nsc_up = ips_es_nsc[,c(1:5)][ips_es_nsc$padj<0.05 & ips_es_nsc$log2FoldChange>0, ],
                        nsc_down = ips_es_nsc[,c(1:5)][ips_es_nsc$padj<0.05 & ips_es_nsc$log2FoldChange<0, ],
                        ipsc_up = ips_es[,c(1:5)][ips_es$padj<0.05 & ips_es$log2FoldChange>0, ],
                        ipsc_down = ips_es[,c(1:5)][ips_es$padj<0.05 & ips_es$log2FoldChange<0, ],
                     as_up=ips_es_as[,c(1,6)])



###
# Make our own gene index
# remove ensembl version number(dot) in ensemble ID and add ENTREZID gene ID from our own gene expression list
Gene_ID_ENTREZ<- ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  left_join(en2ENSE[,1:2],by=c('ensemblID'='ENSEMBL')) %>% dplyr::select(GeneSymbol,ENTREZID)



bg_genes<-Gene_ID_ENTREZ %>% dplyr::select(2) %>% pull() %>% unique()
bg_genes_ense<-ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  dplyr::select(ensemblID) %>% pull() %>% unique() 


MKEGG_ipsc_es_enrich<-function(a) {
  enrich <- enrichMKEGG(de_ips_es_list[[a]]$ENTREZID, organism = "hsa", minGSSize=1,pAdjustMethod = "fdr",pvalueCutoff = 1,qvalueCutoff = 1,
                        universe = bg_genes)
  enrich <- setReadable(enrich, org.Hs.eg.db, 
                        keyType = "ENTREZID")
  return(enrich)
}
?enrichMKEGG

heatplot(MKEGG_ipsc_es_enrich('nsc_down'))


heatplot(MKEGG_ipsc_es_enrich('ipsc_down'))

heatplot(MKEGG_ipsc_es_enrich('as_up'))

heatplot(MKEGG_ipsc_es_enrich('nsc_down'))
heatplot(MKEGG_ipsc_es_enrich('ipsc_down'))