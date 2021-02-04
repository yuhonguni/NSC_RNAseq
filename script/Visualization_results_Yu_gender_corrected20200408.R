# PCA over samples exclude outliers

ExpDataCPM
sample_id_no_outlier<-Metadata %>% filter(!sample_outlier, !possible_outlier) %>%
  select(sample_id) %>% pull(.)

ExpDataCPM_no_outlier<-ExpDataCPM %>% select(sample_id_no_outlier)

dfMeta_no_outlier<- dfMeta %>% mutate(possible_outlier = Metadata$possible_outlier,
                                      sample_outlier = Metadata$sample_outlier) %>%
  filter(!sample_outlier, !possible_outlier)


sample_pca_no_outlier <- prcomp(t(ExpDataCPM_no_outlier))

pcameta_no_outlier<-dfMeta_no_outlier %>%
  unite(celltype_mutation,c('mutation','cell_type'),remove = FALSE)

#
autoplot(sample_pca_no_outlier, data = pcameta_no_outlier, colour="celltype_mutation")


# correlation plot without outeliers
dfMeta <- Metadata %>% filter(!sample_outlier,  !possible_outlier)%>% 
  select(-sample_id, -biological_rep, -source_of_clone, -clone, -homo_het) %>%
  unite("mutation", mutation, subname_mut, remove=TRUE) %>%
  mutate(mutation=gsub("_NA", "", .$mutation)) %>%
  mutate_if(is.character, as.factor) %>%
  mutate_at(vars(matches("fraction$")), scales::rescale) %>%
  mutate(lib_size=scales::rescale(lib_size)) %>%
  as.data.frame
rownames(dfMeta) <- Metadata %>% filter(!sample_outlier,  !possible_outlier) %>%
  pull(sample_id)

View(dfMeta)

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


column_ha = HeatmapAnnotation(#GeneTypeFraction=cbind(mtDNA_Genes=dfMeta$mito_fraction,
                                                     #Ribo_Genes=dfMeta$ribo_fraction,
                                                     #Prot_Genes=dfMeta$prot_fraction,
                                                     #Other_Genes=dfMeta$other_fraction),
                              #DV200=dfMeta$dv200,
                              #Library_Size=dfMeta$lib_size,
                              CellType=dfMeta$cell_type,
                              Mutation=dfMeta$mutation,
                              Homo_Het=dfMeta$homo_het,
                              Subname_Mutation=dfMeta$subname_mut,
                              Clone=dfMeta$clone,
                              Source_Clone=dfMeta$source_of_clone,
                              #Age=dfMeta$age,
                              Sex=dfMeta$sex,
                              #Sample_Outlier=dfMeta$sample.outlier,
                              na_col="white",
                              col=list(#DV200=col_fun_dv200,
                                       #GeneTypeFraction=col_fun_type,
                                       #Library_Size=col_fun_libsize,
                                       CellType=cols_celltype,
                                       Mutation=cols_mutation,
                                       #Age=cols_age,
                                       Sex=cols_sex
                                       #Sample_Outlier=c("TRUE"="red", "FALSE"="grey")
                              ))

sc <- ExpDataCPM %>% select(c('ensemblID',rownames(dfMeta))) %>% 
  filter(ensemblID %in% ProtGenes$gene_id) %>%
  select(-ensemblID) %>%
  as.matrix %>%
  cor
diag(sc) <- NA

Heatmap(sc, column_title="Sample correlation in gene expression", col=viridis::viridis(10), top_annotation=column_ha,
        show_row_names=FALSE, show_column_names=TRUE, show_row_dend=FALSE, column_dend_height=unit(5,"cm"))


#number of differentially expressed genes
DE_number<-data.frame(cell_type = c(rep('IPS',3),rep('NSC',3)),
           mutation = c('WS5A','CP2A','Alpers','WS5A','CP2A','Alpers'),
           DE_number = c(2406,266,2442,3731,10062,5870))

p <- ggplot(DE_number, aes(x=mutation, y=DE_number, fill=cell_type)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  theme_minimal()
p + scale_fill_manual(values=c('#999999','#E69F00'))


#NSC and iPSC markers


ExpDataCPM_no_outlier_id<-ExpDataCPM %>% select( ensemblID, GeneSymbol,sample_id_no_outlier)

ExpDataCPM_no_outlier_id[ExpDataCPM_no_outlier_id$GeneSymbol =='GFAP',]
sample_info <- Metadata %>% filter(!sample_outlier, !possible_outlier) %>%
  select(sample_id,cell_type,mutation,clone,biological_rep,source_of_clone,age,sex) %>%
  data.frame(.)


# neural marker
Marker_expr<- ExpDataCPM_no_outlier_id %>% 
  filter(GeneSymbol %in% c('POU5F1','NANOG','SOX2','NES','PAX6','MAP2','CTNNB1','CDH2',
                           'POLG','POLG2')) %>% 
  select(-1)
Marker_expr <- setNames(data.frame(t(Marker_expr[,-1])),Marker_expr[,1])
Marker_expr<-cbind(sample_info,Marker_expr)

Marker_expr <- Marker_expr %>% mutate(individual = c(rep('Ctr_ES1',6),rep('Ctr_ES2',4),rep('Ctr_IPS',11),rep('WS5A',14),rep('CP2A',18),rep('ALPERS',6)),
                                      cell_stage = Marker_expr$cell_type,
                                      Group = c(rep('Ctr_ES',10),rep('Ctr_IPS',11),rep('WS5A',14),rep('CP2A',18),rep('ALPERS',6)))
Marker_expr[Marker_expr$cell_stage %in% c('ESC','IPSC'),]$cell_stage <-'PSC'


View(Marker_expr)



data <- Marker_expr %>% select(Group, cell_stage,TPM =POLG2) 
  # Calculates mean, sd, se and IC
my_sum <- data %>%
  group_by(Group, cell_stage) %>%
  summarise( 
    n=n(),
    mean=mean(TPM),
    sd=sd(TPM)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
my_sum_point <- data %>%
  group_by(Group, cell_stage)

level_order <- c('Ctr_ES', 'Ctr_IPS', 'WS5A','CP2A','ALPERS')

ggplot()+
  geom_bar(data=my_sum, aes(fill=cell_stage, y=mean, x=factor(Group,level=level_order)),
            position="dodge", stat="identity",colour='black',alpha=0.0)+
  geom_errorbar(data=my_sum, aes(x=factor(Group,level=level_order), ymin=mean-se, ymax=mean+se,fill=cell_stage), 
                width=0.6, colour="black", alpha=1, size=0.6,
                position=position_dodge(0.9))+
  geom_dotplot(data=my_sum_point, aes(fill=cell_stage, y=TPM, x=factor(Group,level=level_order)),
                colour="black",
                binaxis='y', stackdir='center', 
                position=position_dodge(0.9),binwidth = 0.2)+
  theme_classic() + 
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"))




#box plot  
ggplot(data=my_sum_point, aes(x=Group,y=PAX6, fill=cell_stage))+
  geom_boxplot(alpha=0.7)+
  geom_point(position=position_jitterdodge(),alpha=0.7)+
  labs(fill = "cell")


