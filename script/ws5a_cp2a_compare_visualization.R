#Group 1 WS5A and CP2A IPSC
contrast_id <- "IPS-CP2A_vs_IPS-WS5A"
contrast_str <- "subname_WS5A_vs_CP2A"
Mt_wc_ips <- Metadata %>% filter(cell_type %in% c('IPSC','ESC'), (subname_mut == 'WS5A' | subname_mut == 'CP2A')) %>%
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::rename(subname = subname_mut)
Mt_wc_ips %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)
?rename
Ct_wc_ips <- Counts[,colnames(Counts) %in% Mt_wc_ips$sample_id]
rownames(Ct_wc_ips) <- Counts$gene_id
Md_wc_ips <- as.formula(paste0("~subname"))
dds_wc_ips <- DESeq2RUN(Ct_wc_ips, Mt_wc_ips, Md_wc_ips)

res_ws5a_cp2a_ips_1<-results(dds_wc_ips, alpha = 0.05, format = "DataFrame", independentFiltering = T)

res_ws5a_cp2a_ips <- GetDESeq2Results(dds_wc_ips, coef=contrast_str, geneNames=geneNames) %>%arrange(padj)
ws5a_cp2a_ips_de_len<-length(na.omit(res_ws5a_cp2a_ips[res_ws5a_cp2a_ips$padj<0.05,2])) 




de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_ws5a_cp2a_ips, de_file)

#Group 2 WS5A and CP2A NSC

contrast_id <- "NSC-CP2A_vs_NSC-WS5A"
contrast_str <- "subname_WS5A_vs_CP2A"
Mt_wc_nsc <- Metadata %>% filter(cell_type %in% c('NSC'), (subname_mut == 'WS5A' | subname_mut == 'CP2A')) %>%
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::rename(subname = subname_mut)
Mt_wc_nsc  %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)
Ct_wc_nsc <- Counts[,colnames(Counts) %in% Mt_wc_nsc$sample_id]
rownames(Ct_wc_nsc) <- Counts$gene_id
Md_wc_nsc <- as.formula(paste0("~subname"))
dds_wc_nsc <- DESeq2RUN(Ct_wc_nsc, Mt_wc_nsc, Md_wc_nsc)

res_ws5a_cp2a_nsc_1<-results(dds_wc_nsc, alpha = 0.05, format = "DataFrame", independentFiltering = T)

res_ws5a_cp2a_nsc <- GetDESeq2Results(dds_wc_nsc, coef=contrast_str, geneNames=geneNames) %>%arrange(padj)
ws5a_cp2a_nsc_de_len<-length(res_ws5a_cp2a_nsc[res_ws5a_cp2a_nsc$padj<0.05,1]) 

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


ws5a_cp2a_nsc<-import_de('./Yu_tables/DE_NSC-CP2A_vs_NSC-WS5A.tsv')
ws5a_cp2a_ipsc<-import_de("./Yu_tables/DE_IPS-CP2A_vs_IPS-WS5A.tsv")


de_ws5a_cp2a_list<-list(nsc_up = ws5a_cp2a_nsc[,c(1:5)][ws5a_cp2a_nsc$padj<0.05 & ws5a_cp2a_nsc$log2FoldChange>0, ],
                        nsc_down = ws5a_cp2a_nsc[,c(1:5)][ws5a_cp2a_nsc$padj<0.05 & ws5a_cp2a_nsc$log2FoldChange<0, ],
                        ipsc_up = ws5a_cp2a_ipsc[,c(1:5)][ws5a_cp2a_ipsc$padj<0.05 & ws5a_cp2a_ipsc$log2FoldChange>0, ],
                        ipsc_down = ws5a_cp2a_ipsc[,c(1:5)][ws5a_cp2a_ipsc$padj<0.05 & ws5a_cp2a_ipsc$log2FoldChange<0, ])

##
human_mito_carta<-read.table('./human_mito_carta.txt',sep='\t',header=TRUE)
colnames(human_mito_carta)<-c('ENTREZID','Ensembl')



de_ws5a_cp2a_ipsc_up_mito<-as.character(intersect(de_ws5a_cp2a_list[['ipsc_up']]$ENTREZID,human_mito_carta$ENTREZID))

de_ws5a_cp2a_ipsc_up_mito2<-as.character(intersect(de_ws5a_cp2a_list[['ipsc_up']]$EnsemblID,human_mito_carta$Ensembl))

de_ws5a_cp2a_ipsc_down_mito<-as.character(intersect(de_ws5a_cp2a_list[['ipsc_down']]$ENTREZID,human_mito_carta$ENTREZID))

de_ws5a_cp2a_ipsc_down_mito2<-as.character(intersect(de_ws5a_cp2a_list[['ipsc_down']]$EnsemblID,human_mito_carta$Ensembl))



de_ws5a_cp2a_nsc_up_mito<-as.character(intersect(de_ws5a_cp2a_list[['nsc_up']]$ENTREZID,human_mito_carta$ENTREZID))

de_ws5a_cp2a_nsc_up_mito2<-as.character(intersect(de_ws5a_cp2a_list[['nsc_up']]$EnsemblID,human_mito_carta$Ensembl))

de_ws5a_cp2a_nsc_down_mito<-as.character(intersect(de_ws5a_cp2a_list[['nsc_down']]$ENTREZID,human_mito_carta$ENTREZID))

de_ws5a_cp2a_nsc_down_mito2<-as.character(intersect(de_ws5a_cp2a_list[['nsc_down']]$EnsemblID,human_mito_carta$Ensembl))

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
  x = list(a3,a4, human_mito_carta$Ensembl),
  category.names = c("nsc_up" ,"nsc_down" , "mito_cata"),
  filename = 'nsc',
  output = TRUE
)

venn.diagram(
  x = list(a1,a2, human_mito_carta$Ensembl),
  category.names = c("ips_up" ,"ips_down" , "mito_cata"),
  filename = 'ips',
  output=TRUE
)



## astrocyte WS5A vs CP2A

setwd('')

WS5A_CP2A_as_de<-read.table("/home/yu/Postdoc_project/NSC_RNAseq/WS5A_CP2A_paper/WS5A_CP2A_as.txt",sep='\t',header=TRUE,stringsAsFactors = F)



colnames(WS5A_CP2A_as_de)<-c('Gene_ID' ,'Gene_Symbol', 'log2FC', 'Qvalue', 'Pvalue', 'Type','Ensembl')
WS5A_CP2A_as_de$log2FC<-(WS5A_CP2A_as_de$log2FC)*(-1)
WS5A_CP2A_as_de$Gene_ID<-as.character(WS5A_CP2A_as_de$Gene_ID)


length(WS5A_CP2A_as_de$log2FC[WS5A_CP2A_as_de$log2FC<0])
length(WS5A_CP2A_as_de$log2FC[WS5A_CP2A_as_de$log2FC>0])


as_de_number<-data.frame(Direction=c('up','down'),
                         DE_number=c(16,25))
p2 <- ggplot(as_de_number, aes(x=Direction, y=DE_number)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_classic() + 
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")) +
  scale_fill_brewer(palette="Blues")

p2

intersect(WS5A_CP2A_as_de$Ensembl,human_mito_carta$Ensembl)

enrich_as <- enrichMKEGG(WS5A_CP2A_as_de$Gene_ID, organism = "hsa", minGSSize=1,
                       universe = bg_genes,qvalueCutoff = 1, pvalueCutoff = 1,maxGSSize = 500)

enrich <- enrichKEGG(WS5A_CP2A_as_de$Gene_ID, organism = "hsa", keyType = "kegg",
                     pvalueCutoff = 0.05, pAdjustMethod = "BH",qvalueCutoff = 1,
                     universe = bg_genes)


enrich_as<-enricher(WS5A_CP2A_as_de$Gene_ID, 
         pvalueCutoff = 1,pAdjustMethod = "BH",
         qvalueCutoff = 1, minGSSize = 1, maxGSSize = 200,
         TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

heatplot(enrich_as,showCategory = 33)


##heatmap for ipsc DE genes mito

## first, import the RPKM value of all genes

ensembl_length<-read.csv('/home/yu/Postdoc_project/NSC_RNAseq/RNA_seq/metabolic_analysis/ensemble_gene_length.csv')

Count_matrix_Ense<-left_join(cpmMatrixFiltered, geneNames %>% dplyr::select(gene_id, gene_name)) %>% 
  mutate(gene_id = gsub('[.][0-9]+$','',.$gene_id))

Count_matrix_Ense<-left_join(ensembl_length[,c(1,2)],Count_matrix_Ense,by=c('X'='gene_id'))

ExpData_RFKM<-data.frame(ensemblID=Count_matrix_Ense$X,GeneSymbol=Count_matrix_Ense$gene_name,
                         log2(Count2CPM(Count_matrix_Ense[,c(3:(length(colnames(Count_matrix_Ense))-1))])*1000/Count_matrix_Ense[,2]+1))

##ipsc a

de_ws5a_cp2a_ipsc_mito2<-c(de_ws5a_cp2a_ipsc_up_mito2,de_ws5a_cp2a_ipsc_down_mito2)

a<-na.omit(res_ws5a_cp2a_ips[res_ws5a_cp2a_ips$padj<0.05,c(2,3)]) %>% 
  mutate(EnsemblID = gsub('\\.[0-9]+$','',.$EnsemblID))

hp<-ExpData_RFKM[ExpData_RFKM$ensemblID %in% de_ws5a_cp2a_ipsc_mito2 ,c('ensemblID',Mt_wc_ips$sample_id)]  %>% 
  left_join(a,by=c('ensemblID'='EnsemblID')) %>% mutate(log2FoldChange=-(hp$log2FoldChange))

hp<-hp %>% mutate(log2FoldChange=-(.$log2FoldChange)) %>% mutate(FC=cut(.$log2FoldChange, 
                                                             breaks=c(-Inf, 0, Inf), 
                                                             labels=c("down_regulated","up_regulated")))

mat <- t(scale(t(data.matrix(hp[,c(2:(length(colnames(hp))-2))]))))
rownames(mat) <- paste0("row_", seq(nrow(FC)))

FC<-data.frame(hp$FC)
rownames(FC) <- paste0("row_", seq(nrow(FC)))

pheatmap(na.omit(mat), 
         clustering_method = "ward.D",
         color = colorRampPalette(c("green",'black', 'red'))(50),
         scale = 'row',
         cluster_cols=F,
         annotation_row =FC)



## heatmap for astrocyte DE genes

library(pheatmap)

WS5A_CP2A_as_de_count<-read.csv("/home/yu/Postdoc_project/NSC_RNAseq/WS5A_CP2A_paper/WS5A_vs_CP2A_astrocytes_DEGs.xls",header=TRUE)
rownames(WS5A_CP2A_as_de_count) <- paste0("row_", seq(nrow(WS5A_CP2A_as_de_count)))


View(WS5A_CP2A_as_de_count)

WS5A_CP2A_as_de_count[,c(7:21)]<-log2(WS5A_CP2A_as_de_count[,c(7:21)] + 1)
colnames(WS5A_CP2A_as_de_count)[7:21]<-gsub('.Expression','',colnames(WS5A_CP2A_as_de_count[,c(7:21)]))

logFC_as<-data.frame(c(rep('down_regulated',16),
            rep('up_regulated',25)))

rownames(logFC_as) <- paste0("row_", seq(nrow(logFC_as)))
colnames(logFC_as)<-'FC'

logFC_as2<-data.frame(WS5A_CP2A_as_de_count$log2..CP2A.WS5A.)
rownames(logFC_as2) <- paste0("row_", seq(nrow(logFC_as2)))
colnames(logFC_as2)<-'FC'

pheatmap(WS5A_CP2A_as_de_count[,c(7:21)],
         color = colorRampPalette(c("green",'black', 'red'))(50),
         scale = 'row',
         clustering_method = "ward.D",
         cluster_cols=F,
         annotation_row = logFC_as)

## heatmap for NSC mito genes


de_ws5a_cp2a_nsc_mito2<-c(de_ws5a_cp2a_nsc_up_mito2,de_ws5a_cp2a_nsc_down_mito2)


b<-na.omit(res_ws5a_cp2a_nsc[res_ws5a_cp2a_nsc$padj<0.05,2])
hp2<-ExpData_RFKM[ExpData_RFKM$ensemblID %in% de_ws5a_cp2a_nsc_mito2 ,c('ensemblID',Mt_wc_nsc$sample_id)]

res_ws5a_cp2a_nsc2<-res_ws5a_cp2a_nsc %>% mutate(EnsemblID = gsub('\\.[0-9]+$','',res_ws5a_cp2a_nsc$EnsemblID)) 

hp2<-hp2 %>% inner_join(res_ws5a_cp2a_nsc2[,c(2,3)],by = c('ensemblID'='EnsemblID')) %>% arrange(log2FoldChange)

mat <- t(scale(t(data.matrix(hp2[,c(2:(length(colnames(hp2))-1))]))))

pheatmap(na.omit(mat), clustering_method = "ward.D", scale = 'row',border_color = NA)

View(hp[,c(2:length(colnames(hp)))])