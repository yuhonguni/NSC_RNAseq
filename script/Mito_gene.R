library(magrittr)
library(tidyverse)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)


# extract ensemble and entrze gene ID mapping files
k<-keys(org.Hs.eg.db,keytype='ENSEMBL')
en2ENSE<-AnnotationDbi::select(org.Hs.eg.db,keys=k,
                               columns=c("ENTREZID",'SYMBOL'),
                               keytype = 'ENSEMBL')


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
alpers_IPS<-import_de('Yu_tables/DE_IPSC-ALPERSvs_IPSC-CTRL-all.tsv')
View(Metadata)

## sample ID from each group
## NSC
ctr_NSC_id <- Metadata %>% filter(cell_type %in% c('NSC'), (is.na(subname_mut))) %>%
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::select(sample_id) %>% pull()

ws5a_NSC_id <- Metadata %>% filter(cell_type %in% c('NSC'), subname_mut == 'WS5A') %>%
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::select(sample_id) %>% pull()

cp2a_NSC_id <- Metadata %>% filter(cell_type %in% c('NSC'), subname_mut == 'CP2A') %>%
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::select(sample_id) %>% pull()

alpers_NSC_id <- Metadata %>% filter(cell_type %in% c('NSC'), subname_mut == 'Alpers') %>%
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::select(sample_id) %>% pull()

polg_nsc_id<-c(ws5a_NSC_id,cp2a_NSC_id,alpers_NSC_id)

## IPS
ctr_IPS_id <- Metadata %>% filter(cell_type %in% c('IPSC','ESC'), (is.na(subname_mut))) %>%
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::select(sample_id) %>% pull()

ws5a_IPS_id <- Metadata %>% filter(cell_type %in% c('IPSC','ESC'), subname_mut == 'WS5A') %>%
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::select(sample_id) %>% pull()

cp2a_IPS_id <- Metadata %>% filter(cell_type %in% c('IPSC','ESC'), subname_mut == 'CP2A') %>%
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::select(sample_id) %>% pull()

alpers_IPS_id <- Metadata %>% filter(cell_type %in% c('IPSC','ESC'), subname_mut == 'Alpers') %>%
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::select(sample_id) %>% pull()

polg_ips_id<-c(ws5a_IPS_id,cp2a_IPS_id,alpers_IPS_id)


# mito gene symbol and entrezid gene ID

MitoGenes_id<-MitoGenes %>% mutate(gene_id = gsub('[.][0-9]+$','',.$gene_id)) %>% 
  filter(gene_type=='protein_coding') %>% dplyr::select(EnsemblID = gene_id) %>% pull()
mito_id<-en2ENSE %>% filter(ENSEMBL %in% MitoGenes_id) %>%  dplyr::select(ENTREZID,SYMBOL)
length(colnames(ExpDataCPM))




## mitochondria gene expression in all samples
ExpData_mito_tmp<-ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  left_join(en2ENSE[,1:2],by=c('ensemblID'='ENSEMBL')) %>%  filter(ENTREZID %in% mito_id$ENTREZID)

#mean in each sample group
 ## function
mito_mean_sum<-function(x) {
  a<-ExpData_mito_tmp %>% dplyr::select(x) %>% mutate (mean = rowMeans(.),sum=rowSums(.),ENTREZID=ExpData_mito_tmp$ENTREZID) %>% 
    dplyr::select(ENTREZID,mean,sum) 
  return(a)
}

 #mean in both case and control
mean_case_ctr<-function(a,b){
  mean_case_ctr_tmp<-(mito_mean_sum(a)[,3] + mito_mean_sum(b)[,3])/length(c(a,b))
  mean_case_ctr_tmp<-data.frame(ENTREZID=ExpData_mito_tmp$ENTREZID,mean=mean_case_ctr_tmp)
  return(mean_case_ctr_tmp)
}

#plot
require("ggrepel")
mito_gene_maplot<-function(a,b,c) {
  tmp<-a %>% 
    filter(EnsemblID %in% MitoGenes_id) %>% 
    inner_join(mean_case_ctr(b,c),by= 'ENTREZID') %>% 
    inner_join(mito_id,by= 'ENTREZID') %>% 
    mutate(padj_color=cut(.$padj,breaks = c(0, 0.05, 0.1, Inf), labels = c('p<0.05', "p<0.1", "p>0.1"),right  = FALSE))
  
  ggplot(tmp, aes(x=mean, y=log2FoldChange,colour=padj_color)) +
    geom_hline(yintercept=0, linetype="dashed", color = "black")+
    geom_point(size=3)+
    xlim(6, 12) + ylim(-2.2,1)+
    geom_text_repel(aes(label = SYMBOL),
                    size = 3.5,colour = 'black')+
    theme_classic() + 
    theme(
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) 
}


p1<-mito_gene_maplot(alpers_IPS,ctr_IPS_id,alpers_IPS_id)
p2<-mito_gene_maplot(cp2a_IPS,ctr_IPS_id,cp2a_IPS_id)
p3<-mito_gene_maplot(ws5a_IPS,ctr_IPS_id,ws5a_IPS_id)

p4<-mito_gene_maplot(alpers,ctr_NSC_id,alpers_NSC_id)
p5<-mito_gene_maplot(cp2a,ctr_NSC_id,cp2a_NSC_id)
p6<-mito_gene_maplot(ws5a,ctr_NSC_id,ws5a_NSC_id)


require(gridExtra)
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)

### mitochondria encoded genes 

 # control ips
ctr_IPS_mito_expr<-data.frame(cell_type = rep('IPS',12), muatation= rep('Control',12)) %>% 
  cbind(ExpData_mito_tmp %>% dplyr::select(ENTREZID,GeneSymbol,ctr_IPS_id)) %>%
  gather(sample,TPM,5:ncol(.)) %>% mutate(Group = 'Control')

 # polg ips
polg_IPS_mito_expr<-data.frame(cell_type = rep('IPS',12), muatation= rep('polg',12)) %>% 
  cbind(ExpData_mito_tmp %>% dplyr::select(ENTREZID,GeneSymbol,polg_ips_id)) %>%
  gather(sample,TPM,5:ncol(.)) %>% mutate(Group = NA)
polg_IPS_mito_expr[(polg_IPS_mito_expr$sample %in% ws5a_IPS_id),'Group'] = 'ws5a'
polg_IPS_mito_expr[(polg_IPS_mito_expr$sample %in% cp2a_IPS_id),'Group'] = 'cp2a'
polg_IPS_mito_expr[(polg_IPS_mito_expr$sample %in% alpers_IPS_id),'Group'] = 'alpers'

 # combine control ips and polg ips
IPS_mito_expr <- rbind(ctr_IPS_mito_expr,polg_IPS_mito_expr)


 #  control NSC
ctr_NSC_mito_expr<-data.frame(cell_type = rep('NSC',12), muatation= rep('Control',12)) %>% 
  cbind(ExpData_mito_tmp %>% dplyr::select(ENTREZID,GeneSymbol,ctr_NSC_id)) %>%
  gather(sample,TPM,5:ncol(.)) %>%mutate(Group = 'Control')

 #  polg NSC
polg_NSC_mito_expr<-data.frame(cell_type = rep('NSC',12), muatation= rep('polg',12)) %>% 
  cbind(ExpData_mito_tmp %>% dplyr::select(ENTREZID,GeneSymbol,polg_nsc_id)) %>%
  gather(sample,TPM,5:ncol(.)) %>% mutate(Group = NA)
polg_NSC_mito_expr[(polg_NSC_mito_expr$sample %in% ws5a_NSC_id),'Group'] = 'ws5a'
polg_NSC_mito_expr[(polg_NSC_mito_expr$sample %in% cp2a_NSC_id),'Group'] = 'cp2a'
polg_NSC_mito_expr[(polg_NSC_mito_expr$sample %in% alpers_NSC_id),'Group'] = 'alpers'

 # combine control ips and polg ips
NSC_mito_expr <- rbind(ctr_NSC_mito_expr,polg_NSC_mito_expr)


## wilcox test for compassion between all polg patient and control for mRNAs of mtdna encoded genes

  ## NSC
for (i in ExpData_mito_tmp$GeneSymbol) {
  a<-wilcox.test(ctr_NSC_mito_expr$TPM[ctr_NSC_mito_expr$GeneSymbol == i],
         polg_NSC_mito_expr$TPM[polg_NSC_mito_expr$GeneSymbol == i],alternative = 'greater')
  print(i)
  print(a)
}

for (i in ExpData_mito_tmp$GeneSymbol) {
  a<-wilcox.test(ctr_IPS_mito_expr$TPM[ctr_IPS_mito_expr$GeneSymbol == i],
            polg_IPS_mito_expr$TPM[polg_IPS_mito_expr$GeneSymbol == i],alternative = 'greater')
  print(i)
  print(a)
}

## plot 
NSC_mito_expr 
IPS_mito_expr

p1 <- ggboxplot(NSC_mito_expr, x = "muatation", y = "TPM",
               color = "muatation", palette = "jco",
               add = "jitter", facet.by = "GeneSymbol")

p1 + stat_compare_means(method = "t.test", 
                        method.args = list(alternative = "less"),
                        label =  "p.signif", label.x = 1.5, label.y = 11.9)

p2 <- ggboxplot(IPS_mito_expr, x = "muatation", y = "TPM",
               color = "muatation", palette = "jco",
               add = "jitter", facet.by = "GeneSymbol")
p2 + stat_compare_means(method = "t.test", method.args = list(alternative = "less"),
                        label =  "p.signif", 
                        label.x = 1.5, label.y = 12.5)

all_mito_expr<-rbind(cp2a_IPS_mito_expr,cp2a_NSC_mito_expr,alpers_IPS_mito_expr,alpers_NSC_mito_expr,ctr_IPS_mito_expr,ctr_NSC_mito_expr)

## mitochondria encoded genes expression in iPSC and NSC, respectively

ctr_mito_expr<-rbind(ctr_IPS_mito_expr,ctr_NSC_mito_expr)

p1<- ggplot(ctr_mito_expr, aes(x=geneSymbol, y=TPM, fill=cell_type)) +
  geom_boxplot(position=position_dodge(1)) + 
  scale_fill_brewer(palette="RdBu") #+ theme_classic()

    ## geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), size=0.5)

cp2a_mito_expr<-rbind(cp2a_IPS_mito_expr,cp2a_NSC_mito_expr)

p2<- ggplot(cp2a_mito_expr, aes(x=GeneSymbol, y=TPM, fill=cell_type)) +
  geom_boxplot(position=position_dodge(1)) + 
  scale_fill_brewer(palette="Blues") #+ theme_classic()

ws5a_mito_expr<-rbind(ws5a_IPS_mito_expr,ws5a_NSC_mito_expr)

p3<- ggplot(ws5a_mito_expr, aes(x=GeneSymbol, y=TPM, fill=cell_type)) +
  geom_boxplot(position=position_dodge(1)) + 
  scale_fill_brewer(palette="Blues") #+ theme_classic()

alpers_mito_expr<-rbind(alpers_IPS_mito_expr,alpers_NSC_mito_expr)

p4<- ggplot(alpers_mito_expr, aes(x=GeneSymbol, y=TPM, fill=cell_type)) +
  geom_boxplot(position=position_dodge(1)) + 
  scale_fill_brewer(palette="Blues") # + theme_classic()

require(gridExtra)
grid.arrange(p1,p2,p4, ncol=1)

##another way to illustrate

all_mito_expr$muatation<-factor(all_mito_expr$muatation,levels = c('Control', 'CP2A', 'alpers'))

mito_figure<-function(gene_name){
  p<- ggplot(all_mito_expr[all_mito_expr$GeneSymbol == gene_name ,], 
             aes(x=cell_type, y=TPM, fill=muatation)) + 
    geom_boxplot(position=position_dodge(1)) + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5,position=position_dodge(1))+ 
    scale_fill_brewer(palette="Blues") + 
    ggtitle(gene_name)+ 
    theme(plot.title = element_text(hjust = 0.5))
return(p)
}

mito_figure(ExpData_mito_tmp$GeneSymbol[1])

for (i in 1:length(ExpData_mito_tmp$GeneSymbol)){
  p<-mito_figure(ExpData_mito_tmp$GeneSymbol[i])
  assign(paste('x',i,sep=''),p)
}
grid.arrange(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12, ncol=4)
