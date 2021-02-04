library(tidyverse)
library(magrittr)
library(knitr)
library(kableExtra)
library(ggplot2)
library(ggfortify)
library(ggrepel)
source('functions.R')

# load("raw_data_input.RData")
# 1) Annotation file
genes2transcripts <- read_csv("/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/Data/counts_all_samples.gene_name_mappings") %>%
  dplyr::select(gene_id, transcript_id, gene_name, transcript_type, gene_type, seqnames) %>%
  filter(!is.na(transcript_id)) %>%
  unique

# 2) Sample metadata
Metadata <- read_csv("/home/yu/PostDocProject/NSC_RNAseq/RNA_seq//Data/metadata_xiao.csv") %>%
  mutate(biological_rep=as.factor(biological_rep), source_of_clone=as.factor(source_of_clone)) %>%
  dplyr::select(sample_id, dv200, cell_type, mutation, homo_het, subname_mut, clone, biological_rep, source_of_clone, age)

Data <- readRDS('/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/Data/countMatrix.Rds')

if (!exists("countMatrix")) {countMatrix <- readRDS("/home/yu/PostDocProject/NSC_RNAseq/RNA_seq//Data/countMatrix.Rds")}
geneNames <- genes2transcripts %>%
  dplyr::select(gene_id, gene_name, gene_type, seqnames) %>%
  unique

MitoGenes <- geneNames[geneNames$seqnames=="chrM",]
RiboGenes <- geneNames[geneNames$gene_type %in% c("rRNA", "Mt_rRNA"),]
ProtGenes <- geneNames[geneNames$gene_type=="protein_coding",]

# Now we remove the chr names from the geneNames because there are duplicate
# gene names in chrX and chrY 
geneNames %<>% dplyr::select(gene_id, gene_name, gene_type) %>%
  unique

## 2 Quality control
## 2.1 Transcript biotypes

# Counts for mitochondrial and rRNA genes
MitoCountSum <- colSums(countMatrix %>% filter(gene_id %in% MitoGenes$gene_id) %>% dplyr::select(-gene_id))
MitoCountSum <- data.frame(sample_id=names(MitoCountSum), mito_counts=round(MitoCountSum), stringsAsFactors=FALSE)
RiboCountSum <- colSums(countMatrix %>% filter(gene_id %in% RiboGenes$gene_id) %>% dplyr::select(-gene_id))
RiboCountSum <- data.frame(sample_id=names(RiboCountSum), ribo_counts=round(RiboCountSum), stringsAsFactors=FALSE)
ProtCountSum <- colSums(countMatrix %>% filter(gene_id %in% (ProtGenes %>% filter(seqnames!="chrM"))$gene_id) %>% dplyr::select(-gene_id))
ProtCountSum <- data.frame(sample_id=names(ProtCountSum), prot_counts=round(ProtCountSum), stringsAsFactors=FALSE)
TotCountSum  <- colSums(countMatrix %>% dplyr::select(-gene_id))
TotCountSum <- data.frame(sample_id=names(TotCountSum), counts=TotCountSum, stringsAsFactors=FALSE)

ReadFractions <- inner_join(MitoCountSum, RiboCountSum, by="sample_id") %>%
  inner_join(ProtCountSum, by="sample_id") %>%
  inner_join(TotCountSum, by="sample_id") %>%
  mutate(mito_fraction=mito_counts/counts, ribo_fraction=ribo_counts/counts, prot_fraction=prot_counts/counts,
         other_fraction=1-(mito_fraction+ribo_fraction+prot_fraction), lib_size=as.integer(round(counts))) %>%
  dplyr::select(sample_id, mito_fraction, ribo_fraction, prot_fraction, other_fraction, lib_size=counts)

Metadata %<>% left_join(ReadFractions)

ReadFractions %>% dplyr::select(-lib_size) %>% gather("Gene_biotype", "Proportion", -sample_id) %>%
  arrange(Proportion) %>% mutate(Gene_biotype=factor(Gene_biotype, levels=c("ribo_fraction","mito_fraction","other_fraction","prot_fraction"))) %>%
  ggplot(aes(x=sample_id, y=Proportion, fill=Gene_biotype)) +
  geom_bar(position="stack", stat="identity", colour="white") +
  labs(y="Proportion of reads", x="Sample") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks=element_blank())

##2.2 Highly expressed genes
## highly expressed 

# Counts excluding mitochondrial genes and mitochondrial reads fraction over
# total lib. size
minProp <- 0.01
MitoCountFiltered <- countMatrix %>% filter(!gene_id %in% MitoGenes$gene_id)
CountFreqs <- prop.table(as.matrix(MitoCountFiltered[,-1]), 2)
rownames(CountFreqs) <- MitoCountFiltered$gene_id

#CountFreqs[CountFreqs < minProp] <- NA

HighExpGenes <- CountFreqs[apply(CountFreqs, 1, function(x) any(x >= minProp)),] %>%
  as.data.frame %>% mutate(gene_id=rownames(.)) %>%
  left_join(geneNames) %>% 
  gather("sample_id", "Proportion", -gene_id, -gene_name, -gene_type)

s_order <- HighExpGenes %>% group_by(sample_id) %>% summarise(sum=sum(Proportion)) %>%
  arrange(sum) %>% pull(sample_id) 

g_order <- HighExpGenes %>% group_by(gene_name) %>% summarise(sum=sum(Proportion)) %>%
  arrange(sum) %>% pull(gene_name) 

HighExpGenes$sample_id <- factor(HighExpGenes$sample_id, levels=s_order)
HighExpGenes$gene_name <- factor(HighExpGenes$gene_name, levels=g_order)

ggplot(HighExpGenes, aes(sample_id, Proportion, fill=gene_name)) +
  geom_bar(stat = "identity", colour="white") +
  labs(y="Proportion of reads", x="Sample") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks=element_blank())


##Now that we have checked for highly abundant transcripts, we can transform the count data to log2(CPM) and do some more filtering on the non-mitochondrial genes. Genes are filtered out if:
  
##  no Gene Symbol (or duplicated)
##  expression is below the median expression for more than 20% of the samples

MitoGenes_exclude <- MitoGenes %>% filter(!(gene_type == 'protein_coding' & gene_name != 'MT-ND6'))

MitoCountFiltered<-countMatrix %>% filter(!gene_id %in% MitoGenes_exclude$gene_id)

countMatrixFiltered <- MitoCountFiltered

GenesKept <- rep(0, 3)
names(GenesKept) <- c("Total nuclear transcripts", "Trancripts above threshold", "With Gene Symbol")

# Create log2 CPM matrix after removal of mitochondria-encoded genes
cpmMatrixFiltered <- data.frame(gene_id=as.character(countMatrixFiltered$gene_id), log2(Count2CPM(countMatrixFiltered[,-1])+1))
GenesKept[1] <- nrow(cpmMatrixFiltered)

# Add gene symbols
ExpDataCPM <- left_join(cpmMatrixFiltered, geneNames %>% dplyr::select(gene_id, gene_name)) %>%
  dplyr::select(ensemblID=gene_id, GeneSymbol=gene_name, everything())

# Filter low-expressed genes similar to https://f1000research.com/articles/5-1438
# (although a bit stricter: keep genes that have at least 10 reads in 25% of
# the samples)
min_lib_size <- min(Metadata$lib_size)
min_cpm <- round(10/(min_lib_size/1e+6),1)
keep <- rowSums(cpmMatrixFiltered[,-1] > min_cpm) >= round(nrow(Metadata)*.25)
#keep %>% table

# Filter genes below noise level
ExpDataCPM <- ExpDataCPM[keep,]
GenesKept[2] <- nrow(ExpDataCPM)

# Filter out genes with no Gene Symbol or with duplicated Gene Symbol (highest
# expressed option of duplicates is kept)
ExpDataCPM %<>% filter(GeneSymbol != "", !is.na(GeneSymbol)) %>%
  mutate(Sum = do.call(pmax, select_if(., is.numeric))) %>%
  arrange(desc(Sum)) %>% 
  distinct(GeneSymbol, .keep_all=TRUE) %>%
  dplyr::select(-Sum) %>%
  arrange(ensemblID)

GenesKept[3] <- nrow(ExpDataCPM)
prot_coding_n <- c(geneNames %>% filter(gene_id %in% ExpDataCPM$ensemblID) %>%
                     pull(gene_type) %>% table %>% prop.table)["protein_coding"]


GenesKept %>% kable(caption="Number of genes passing filters") %>%
  kable_styling(bootstrap_options = "striped", full_width = F)

# Reorder samples in CPM matrix
if ( ! all(colnames(ExpDataCPM)[-c(1,2)] == Metadata$sample_id) ) {
  warning("Reordering CPM matrix...")
  ExpDataCPM <- ExpDataCPM[,c(1,2,(match(Metadata$sample_id, colnames(ExpDataCPM[,-c(1,2)]))+2))]
} else {
  message("No need to reorder CPM matrix")
}

# Filter low-expressed genes from counts matrix and reorder if necessary
Counts <- countMatrixFiltered %>%
  filter(gene_id %in% ExpDataCPM$ensemblID) %>%
  mutate_if(is.numeric, function(x) as.integer(round(x)))

if ( ! all(colnames(Counts)[-1] == Metadata$sample_id) ) {
  warning("Reordering count matrix...")
  Counts <- Counts[,c(1,(match(Metadata$sample_id, colnames(Counts[,-1]))+1))]
} else {
  message("No need to reorder count matrix")
}

##2.3 Assign genetic sex
##We will try to assign genetic sex to each sample based on the expression of X and Y markers.

sexMarkers <- read_csv(file="/home/yu/PostDocProject/NSC_RNAseq/RNA_seq/Data/sex_specific_genes.csv")

sex_marker_expression <- sexMarkers %>% dplyr::select(-gene_id) %>%
  left_join(ExpDataCPM) %>% select_if(is.numeric) %>%
  as.matrix
rownames(sex_marker_expression) <- sexMarkers$GeneSymbol

sex_pca <- prcomp(t(sex_marker_expression))
pc1 <- sex_pca$rotation[,1,drop=FALSE] %>% as.data.frame %>% mutate(GeneSymbol=rownames(.)) %>%
  left_join(sexMarkers) %>% mutate(pc1_dir=ifelse(sign(PC1)==1,1,0), sex_dir=ifelse(sex=="male", 0, 1))

if (!all(xor(pc1$pc1_dir,pc1$sex_dir)[1] == xor(pc1$pc1_dir, pc1$sex_dir)) ) {
  stop("Not all sex markers go in the expected direction")
} else {
  if (pc1$pc1_dir[1]!=pc1$sex_dir[1]) {
    # Positive is male
    sample_sex <- data.frame(sample_id=rownames(sex_pca$x), sex=ifelse(sign(sex_pca$x[,1])<0,"female","male"), stringsAsFactors=FALSE)
  } else {
    # Positive is female
    sample_sex <- data.frame(sample_id=rownames(sex_pca$x), sex=ifelse(sign(sex_pca$x[,1])>0,"female","male"), stringsAsFactors=FALSE)
  }
}

sex_pca$x
autoplot(sex_pca, loadings=TRUE, loadings.label=TRUE)

Metadata %<>% left_join(sample_sex) 

#add the sex marker expression data to the sample information table
#Metadata_with_sex_marker<-sex_marker_expression %>%
#  t(.) %>% data.frame() %>% mutate(sample_id=rownames(.)) %>%
#  inner_join(Metadata,by='sample_id')
#write.csv(Metadata_with_sex_marker,'Metadata_with_sex_marker.csv')


CTRL3_samples <- Metadata %>% filter(clone=="CTRL3") %>% pull(sample_id)

sex_marker_expression %>% as.data.frame %>%
  mutate(GeneSymbol=rownames(.)) %>%
  gather(sample_id, log_cpm, -GeneSymbol) %>%
  left_join(sample_sex) %>%
  ggplot(aes(sex, log_cpm)) +
  geom_boxplot(colour="grey40") +
  geom_point(aes(colour=sex)) +
  geom_text_repel(data=(. %>% filter(sample_id %in% CTRL3_samples)), aes(sex, log_cpm, label=sample_id)) +
  scale_colour_manual(values=c(male="darkturquoise", female="palevioletred")) +
  facet_grid(.~GeneSymbol)

##2.4 Sample outliers

## Now that we have filtered the genes, we can look at the samples and identify outliers. 
## To that end, we calculate the pairwise correlations in expression between samples. 
## We characterise each sample by its median correlation to the rest of the samples. 
## Based on this value, we mark as outliers the samples outside the interval Q1-1.5*IQR,R3+1.5*IQR 
## (Q1=first quartile, Q3=third quartile, IQR=interquartile range). 
## We restrict the correlation to protein-coding genes.

sc <- ExpDataCPM %>%
  filter(ensemblID %in% ProtGenes$gene_id) %>%
  dplyr::select(-ensemblID, -GeneSymbol) %>%
  as.matrix %>%
  cor
diag(sc) <- NA
medians <- matrixStats::rowMedians(sc, na.rm=TRUE)
names(medians) <- colnames(sc)

iqrange <- IQR(medians)
quartiles <- quantile(medians, c(0.25, 0.75))
lim_low <- quartiles[1]-(1.5*iqrange)
lim_high <- quartiles[2]+(1.5*iqrange)

sample.outliers <- (medians < lim_low) | (medians > lim_high)

# Update object
Metadata <- bind_cols(Metadata, sample_outlier=sample.outliers)
#sex mismatch
Metadata[Metadata$sample_id=="SL403669",]$sample_outlier <- TRUE
Metadata[Metadata$clone=="CTRL2",]$sample_outlier <- TRUE
Metadata[Metadata$clone=="ESC2",]$sample_outlier <- TRUE

## 3 Sample characterisation
## 3.1 Sample clustering
## We now plot all the samples and their correlation to further explore outliers.

library(ComplexHeatmap)

# Create temporal Metadata data.frame
dfMeta <- Metadata %>% dplyr::select(-sample_id, -biological_rep, -source_of_clone, -clone, -homo_het) %>%
  unite("mutation", mutation, subname_mut, remove=TRUE) %>%
  mutate(mutation=gsub("_NA", "", .$mutation)) %>%
  mutate_if(is.character, as.factor) %>%
  mutate_at(vars(matches("fraction$")), scales::rescale) %>%
  mutate(lib_size=scales::rescale(lib_size)) %>%
  as.data.frame
rownames(dfMeta) <- Metadata$sample_id

col_fun_type <- circlize::colorRamp2(c(0, 1), c("black", "orange3"))
col_fun_dv200 <- circlize::colorRamp2(c(95, 99), c("red3", "springgreen3"))
col_fun_libsize <- circlize::colorRamp2(c(0,1), c("red3", "springgreen3"))
library(RColorBrewer)
cols_age <- brewer.pal(length(levels(dfMeta$age)), "Set1")
names(cols_age) <- levels(dfMeta$age)

cols_mutation <- brewer.pal(length(levels(dfMeta$mutation)), "Set2")
names(cols_mutation) <- levels(dfMeta$mutation)

cols_celltype <- brewer.pal(length(levels(dfMeta$cell_type)), "Accent")
names(cols_celltype) <- levels(dfMeta$cell_type)
cols_sex <- c(male="darkturquoise", female="palevioletred")


column_ha = HeatmapAnnotation(GeneTypeFraction=cbind(mtDNA_Genes=dfMeta$mito_fraction,
                                                     Ribo_Genes=dfMeta$ribo_fraction,
                                                     Prot_Genes=dfMeta$prot_fraction,
                                                     Other_Genes=dfMeta$other_fraction),
                              DV200=dfMeta$dv200,
                              Library_Size=dfMeta$lib_size,
                              CellType=dfMeta$cell_type,
                              Mutation=dfMeta$mutation,
                              Homo_Het=dfMeta$homo_het,
                              Subname_Mutation=dfMeta$subname_mut,
                              Clone=dfMeta$clone,
                              Source_Clone=dfMeta$source_of_clone,
                              Age=dfMeta$age,
                              Sex=dfMeta$sex,
                              #Sample_Outlier=dfMeta$sample.outlier,
                              na_col="white",
                              col=list(DV200=col_fun_dv200,
                                       GeneTypeFraction=col_fun_type,
                                       Library_Size=col_fun_libsize,
                                       CellType=cols_celltype,
                                       Mutation=cols_mutation,
                                       Age=cols_age,
                                       Sex=cols_sex
                                       #Sample_Outlier=c("TRUE"="red", "FALSE"="grey")
                              ))
Heatmap(sc, column_title="Sample correlation in gene expression", col=viridis::viridis(10), top_annotation=column_ha,
        show_row_names=FALSE, show_column_names=TRUE, show_row_dend=FALSE, column_dend_height=unit(5,"cm"))

other_outliers <- c("SL403684", "SL403687")
outlier_cluster <- c("SL403678","SL403688","SL403648","SL403649")
Metadata[Metadata$sample_id %in% other_outliers,]$sample_outlier <- TRUE
Metadata$possible_outlier <- FALSE
Metadata[Metadata$sample_id %in% outlier_cluster,]$possible_outlier <- TRUE

## 3.2 Correlation between variables

## We can plot the pair-wise Pearson’s correlation between the variables to have 
## a general picture of which ones are more associated.

library(corrplot)

dfDummy <- fastDummies::dummy_cols(dfMeta) %>%
  select_if(is.numeric)
rownames(dfDummy) <- rownames(dfMeta)

M <- cor(dfDummy)
res <- cor.mtest(dfDummy, conf.level=0.95)
diag(M) <- NA

corrplot.mixed(M, tl.pos="lt", tl.col="black", tl.srt=45, order="original", na.label=" ")

##3.3 PCA over all samples and transcripts

## Now we plot the PCA of all samples.

sample_pca <- prcomp(t(ExpDataCPM[,-c(1,2)]))

pcameta<-dfMeta %>%
  unite(celltype_mutation,c('mutation','cell_type','age','sex'),remove = FALSE)
#

autoplot(sample_pca, data=pcameta, colour="celltype_mutation")


assocPCA(M=t(ExpDataCPM[,-c(1,2)]), variables=dfDummy, plot.return=TRUE, ncomp=5, plot.alpha=0.05, plot.colour="orangered") +
  ggtitle("Association of PC 1-5 with experimental variables")

# PCA over samples exclude outliers

ExpDataCPM
sample_id_no_outlier<-Metadata %>% filter(!sample_outlier, !possible_outlier) %>%
  dplyr::select(sample_id) %>% pull(.)

ExpDataCPM_no_outlier<-ExpDataCPM %>% dplyr::select(sample_id_no_outlier)

dfMeta_no_outlier<- dfMeta %>% mutate(possible_outlier = Metadata$possible_outlier,
                                      sample_outlier = Metadata$sample_outlier) %>%
  filter(!sample_outlier, !possible_outlier)
View(dfMeta_no_outlier)


sample_pca_no_outlier <- prcomp(t(ExpDataCPM_no_outlier))

pcameta_no_outlier<-dfMeta_no_outlier %>%
  unite(celltype_mutation,c('mutation','cell_type'),remove = FALSE)

#
autoplot(sample_pca_no_outlier, data = pcameta_no_outlier, colour="celltype_mutation")




## 4 Differential expression analysis

# Load libraries and get annotation
require("DESeq2")
require("ermineR")
require("colorspace")

GenericHumanAnno <- GetAnnoFiles("Generic_human")
?GetAnnoFiles

Meta_no_outlier<-Metadata %>% filter(!sample_outlier,  !possible_outlier)

write.csv(Meta_no_outlier, 'Metadata_no_outlier.csv',row.names = F )
length(rownames(Meta_no_outlier))

#Group_1
contrast_id <- "NSC-CP2A-Alpers_vs_NSC-CTRL-all"
contrast_str <- "mutation_POLG_vs_CTRL"
Mt <- Metadata %>% filter(cell_type %in% c('NSC'), (is.na(subname_mut)) | subname_mut %in% c('CP2A','Alpers')) %>%
  filter(!sample_outlier,  !possible_outlier)
Mt %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

Ct <- Counts[,colnames(Counts) %in% Mt$sample_id]
rownames(Ct) <- Counts$gene_id
Md <- as.formula(paste0("~mutation"))
dds <- DESeq2RUN(Ct, Mt, Md)

res_cp2a_alpers_nsc1<-results(dds, alpha = 0.05, format = "DataFrame", independentFiltering = T)

res_cp2a_alpers_nsc <- GetDESeq2Results(dds, coef=contrast_str, geneNames=geneNames) %>%arrange(padj)
cp2a_nsc_alpers_de_len<-length(res_cp2a_alpers_nsc[res_cp2a_alpers_nsc$padj<0.05,1])

de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_cp2a_alpers_nsc, de_file)





#Group 2
contrast_id <- "NSC-CP2A_vs_NSC-CTRL-all"
contrast_str <- "mutation_POLG_vs_CTRL"
Mt <- Metadata %>% filter(cell_type %in% c('NSC'), (is.na(subname_mut)) | subname_mut == 'CP2A') %>%
  filter(!sample_outlier,  !possible_outlier)
Mt %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

Ct <- Counts[,colnames(Counts) %in% Mt$sample_id]
rownames(Ct) <- Counts$gene_id
Md <- as.formula(paste0("~mutation"))
dds <- DESeq2RUN(Ct, Mt, Md)

res_cp2a_nsc1<-results(dds, alpha = 0.05, format = "DataFrame", independentFiltering = T)

res_cp2a_nsc <- GetDESeq2Results(dds, coef=contrast_str, geneNames=geneNames) %>%arrange(padj)
cp2a_nsc_de_len<-length(res_CP2A_nsc[res_CP2A_nsc$padj<0.05,1])

de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_cp2a_nsc, de_file)

#Group 3

contrast_id <- 'NSC-ALPERSvs_NSC-CTRL-all'
contrast_str <- "mutation_POLG_vs_CTRL"
Mt <- Metadata %>% filter(cell_type %in% c('NSC'), (is.na(subname_mut)) | subname_mut == 'Alpers') %>%
  filter(!sample_outlier,  !possible_outlier)
Mt %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

Ct <- Counts[,colnames(Counts) %in% Mt$sample_id]
rownames(Ct) <- Counts$gene_id
Md <- as.formula(paste0("~ mutation"))
dds <- DESeq2RUN(Ct, Mt, Md)

res_alpers_nsc1<-results(dds, alpha = 0.05, format = "DataFrame", independentFiltering = T)

res_alpers_nsc <- GetDESeq2Results(dds, coef=contrast_str, geneNames=geneNames) %>%arrange(padj)
alpers_nsc_de_len<-length(res_alpers_nsc[res_alpers_nsc$padj<0.05,1])

de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_alpers_nsc, de_file)



#Group 5
contrast_id <- "IPS-CP2A_vs_IPS-CTRL-all"
contrast_str <- "mutation_POLG_vs_CTRL"
Mt <- Metadata %>% filter(cell_type %in% c('IPSC','ESC'), (is.na(subname_mut)) | subname_mut == 'CP2A') %>%
  filter(!sample_outlier,  !possible_outlier)
Mt %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

Ct <- Counts[,colnames(Counts) %in% Mt$sample_id]
rownames(Ct) <- Counts$gene_id
Md <- as.formula(paste0("~mutation"))
dds <- DESeq2RUN(Ct, Mt, Md)
res_cp2a_ips_1<-results(dds, alpha = 0.05, format = "DataFrame", independentFiltering = T)

res_cp2a_ips <- GetDESeq2Results(dds, coef=contrast_str, geneNames=geneNames) %>%arrange(padj)
cp2a_ips_de_len<-length(res_cp2a_ips[res_cp2a_ips$padj<0.05,1]) #266 #10067

de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_cp2a_ips, de_file)

#Group 6

contrast_id <- 'IPSC-ALPERSvs_IPSC-CTRL-all'
contrast_str <- "mutation_POLG_vs_CTRL"
Mt <- Metadata %>% filter(cell_type %in% c('IPSC','ESC'), (is.na(subname_mut)) | subname_mut == 'Alpers') %>%
  filter(!sample_outlier,  !possible_outlier)
Mt %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

Ct <- Counts[,colnames(Counts) %in% Mt$sample_id]
rownames(Ct) <- Counts$gene_id
Md <- as.formula(paste0("~ mutation"))
dds <- DESeq2RUN(Ct, Mt, Md)

res_alpers_ips_1<-results(dds, alpha = 0.05, format = "DataFrame", independentFiltering = T)

res_alpers_ips <- GetDESeq2Results(dds, coef=contrast_str, geneNames=geneNames) %>%arrange(padj)
alpers_ips_de_len<-length(res_alpers_ips[res_alpers_ips$padj<0.05,1])

de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_alpers_ips, de_file)


#Group_7
contrast_id <- "NSC-CP2A-WS5A_vs_iPSC-CP2A-WS5A-all"
contrast_str <- "cell_type_NSC_vs_IPSC"
Mt <- Metadata %>% filter(subname_mut %in% c('CP2A','WS5A')) %>%
  filter(!sample_outlier,  !possible_outlier)
Mt %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

Ct <- Counts[,colnames(Counts) %in% Mt$sample_id]
rownames(Ct) <- Counts$gene_id
Md <- as.formula(paste0("~cell_type"))
dds <- DESeq2RUN(Ct, Mt, Md)

res_cp2a_ws5a_nsvsips1<-results(dds, alpha = 0.05, format = "DataFrame", independentFiltering = T)

res_cp2a_ws5a_nsvsips <- GetDESeq2Results(dds, coef=contrast_str, geneNames=geneNames) %>%arrange(padj)

down_regulation_list<-res_cp2a_ws5a_nsvsips  %>% filter(padj<0.05, log2FoldChange<0)  %>% select(GeneSymbol)
write_tsv(down_regulation_list, 'res_cp2a_ws5a_nsvsips_down_regulation_list.tsv')


de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_cp2a_alpers_nsc, de_file)


#Group_8
contrast_id <- "NSC-CP2A-WS5A_vs_NSC-CTRL"
contrast_str <- "mutation_POLG_vs_CTRL"

View(Metadata)
Mt <- Metadata %>% filter(cell_type %in% c('NSC'), (subname_mut %in% c('CP2A','WS5A') | mutation %in% c('CTRL'))) %>%
  filter(!sample_outlier,  !possible_outlier)
Mt %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

Ct <- Counts[,colnames(Counts) %in% Mt$sample_id]
rownames(Ct) <- Counts$gene_id
Md <- as.formula(paste0("~sex+mutation"))
dds <- DESeq2RUN(Ct, Mt, Md)

res_cp2a_ws5a_NSC_POLGvsCTRL1<-results(dds, alpha = 0.05, format = "DataFrame", independentFiltering = T)

res_cp2a_ws5a_NSC_POLGvsCTRL <- GetDESeq2Results(dds, coef=contrast_str, geneNames=geneNames) %>%arrange(padj)

res_cp2a_ws5a_NSC_POLGvsCTRL  %>% filter(padj<0.05)  %>% select(GeneSymbol)

de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_cp2a_ws5a_NSC_POLGvsCTRL, de_file)


#Group_9
contrast_id <- "IPS-CP2A-WS5A_vs_IPS-CTRL"
contrast_str <- "mutation_POLG_vs_CTRL"

Mt <- Metadata %>% filter(cell_type %in% c('IPSC','ESC'), (subname_mut %in% c('CP2A','WS5A') | mutation %in% c('CTRL'))) %>%
  filter(!sample_outlier,  !possible_outlier)
Mt %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

Ct <- Counts[,colnames(Counts) %in% Mt$sample_id]
rownames(Ct) <- Counts$gene_id
Md <- as.formula(paste0("~sex+mutation"))
dds <- DESeq2RUN(Ct, Mt, Md)

res_cp2a_ws5a_IPS_POLGvsCTRL1<-results(dds, alpha = 0.05, format = "DataFrame", independentFiltering = T)

res_cp2a_ws5a_IPS_POLGvsCTRL <- GetDESeq2Results(dds, coef=contrast_str, geneNames=geneNames) %>%arrange(padj)

res_cp2a_ws5a_IPS_POLGvsCTRL  %>% filter(padj<0.05)  %>% select(GeneSymbol)


de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_cp2a_ws5a_IPS_POLGvsCTRL, de_file)



getwd()


## DE number plot
DE_number<-data.frame(Group=c('Alpers','CP2A','Alpers','CP2A'),
                      Cell_type=c('iPSC','iPSC','NSC','NSC'),
                      DE_gene_number=c(alpers_ips_de_len,cp2a_ips_de_len,
                                       alpers_nsc_de_len,cp2a_nsc_de_len))
p <- ggplot(DE_number, aes(x=Group, y=DE_gene_number, fill=Cell_type)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_classic() + 
  theme(axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black")) +
  scale_fill_brewer(palette="Blues")



##plot MA plot

par(mfrow=c(2,2), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-15,15)
plotMA(res_alpers_ips_1, xlim=xlim, ylim=ylim, main="Alpers iPSC",alpha = 0.05)
plotMA(res_alpers_nsc1, xlim=xlim, ylim=ylim, main="Alpers NSC",alpha = 0.05)
plotMA(res_cp2a_ips_1, xlim=xlim, ylim=ylim, main="CP2A iPSC",alpha = 0.05)
plotMA(res_CP2A_nsc1, xlim=xlim, ylim=ylim, main="CP2A NSC",alpha = 0.05)

