
## heat map  of alpers together

heatmap_sample_id <- Metadata %>% filter((is.na(subname_mut)) | subname_mut == 'Alpers') %>%
  filter(!sample_outlier,  !possible_outlier) %>% .$sample_id

sc_alpers_ctrl <- ExpDataCPM %>%
  filter(ensemblID %in% ProtGenes$gene_id) %>%
  dplyr::select(heatmap_sample_id) %>%
  as.matrix %>%
  cor

library(ComplexHeatmap)

# Create temporal Metadata data.frame
dfMeta <- Metadata %>% dplyr::select(-sample_id, -biological_rep, -source_of_clone, -clone, -homo_het) %>%
  filter(is.na(subname_mut) | subname_mut == 'Alpers') %>% 
  filter(!sample_outlier,  !possible_outlier) %>%
  unite("mutation", mutation, subname_mut, remove=TRUE) %>%
  mutate(mutation=gsub("_NA", "", .$mutation)) %>%
  mutate_if(is.character, as.factor) %>%
  mutate_at(vars(matches("fraction$")), scales::rescale) %>%
  mutate(lib_size=scales::rescale(lib_size)) %>%
  as.data.frame
rownames(dfMeta) <- Metadata %>% filter(is.na(subname_mut) | subname_mut == 'Alpers')%>% 
  filter(!sample_outlier,  !possible_outlier)  %>%
  .$sample_id

library(RColorBrewer)

cols_celltype <- brewer.pal(length(levels(dfMeta$cell_type)), "Accent")
names(cols_celltype) <- levels(dfMeta$cell_type)

cols_mutation <- c("#7FC97F","#BEAED4")
names(cols_mutation) <- levels(dfMeta$mutation)


column_ha = HeatmapAnnotation(CellType=dfMeta$cell_type,
                              Mutation=dfMeta$mutation,
                              Subname_Mutation=dfMeta$subname_mut,
                              na_col="white",
                              col=list(CellType=cols_celltype,
                                       Mutation = cols_mutation
                              ))

Heatmap(sc_alpers_ctrl, 
        column_title="Sample correlation in gene expression", 
        col=viridis::viridis(10), 
        #col=colorRampPalette(brewer.pal(10, "RdYlBu"))(256),
        top_annotation=column_ha,
        show_row_names=FALSE, 
        show_column_names=FALSE,
        show_column_dend=FALSE, 
        #show_row_dend=FALSE, 
        #column_dend_height=unit(5,"cm"),
        row_dend_width=unit(5,"cm"))
dev.off()

# PCA over samples exclude outliers

ExpDataCPM_alpers_ctrl<-ExpDataCPM %>% dplyr::select(heatmap_sample_id )

dfMeta_alpers_ctrl<- dfMeta %>% filter(rownames(.) %in%  heatmap_sample_id)

  ## calculat PCA

sample_pca_alpers_ctrl <- prcomp(t(ExpDataCPM_alpers_ctrl))

pcameta_alpers_ctrl<-dfMeta_alpers_ctrl %>%
  unite(celltype_mutation,c('mutation','cell_type'),remove = FALSE)

df_pca<-data.frame(sample_pca_alpers_ctrl$x,Group = pcameta_alpers_ctrl$celltype_mutation)

  ## add PC percentage 
percentage<-round(sample_pca_alpers_ctrl$sdev^2/sum(sample_pca_alpers_ctrl$sdev^2)*100,2)
percentage<-paste(colnames(sample_pca_alpers_ctrl$x), '(',paste(as.character(percentage)),'%',')',sep='')

  # ggplot
ggplot(df_pca,aes(x=PC1,y=PC2,colour = Group))+
  geom_point()+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),axis.line=element_line(colour='black'))+
  labs(x=percentage[1],y=percentage[2])
  #stat_ellipse(level=0.9,show.legend=F)
  
  # another way of plot
autoplot(sample_pca_alpers_ctrl, 
         data = pcameta_alpers_ctrl, 
         colour="celltype_mutation")
