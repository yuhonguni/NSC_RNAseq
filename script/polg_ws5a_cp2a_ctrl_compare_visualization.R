
path_list<-gsaRes_nsc_polg_ctrl_tbl %>% filter(`p adj (dist.dir.dn)` < 0.05) %>% 
  dplyr::select(Name) %>% pull(.)

path_gene_fc<-MKEGG_gene_tbl %>% 
  filter(MKEGG_name %in% path_list) %>% 
  left_join(polg_ctrl_nsc[,c(1,2)],by=c('ENSEMBL'='EnsemblID')) %>%
  filter(is.na(.$log2FoldChange) == F)

gene_unqiue<-path_gene_fc %>% distinct(.$gene_name) %>% pull(.)

fc<-data.frame(matrix(nrow=length(path_list),ncol=length(gene_unqiue)))
colnames(fc)<-gene_unqiue
rownames(fc)<-path_list

for (i in c(1:length(path_gene_fc$log2FoldChange))) {
  path<-path_gene_fc$MKEGG_name[i]
  gene<-path_gene_fc$gene_name[i]
  logfc<-path_gene_fc$log2FoldChange[i]
  fc[path,gene]<-logfc
}

library(ComplexHeatmap)

library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#005555", "white", "#CC0000"))
col_fun = colorRamp2(c(-50, 0, 50), c("#005555", "white", "#CC0000"))
col_fun(seq(-50, 50))

Heatmap(data.matrix(fc),na_col='white',
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        border = TRUE,
        col = col_fun,
        rect_gp = gpar(col = "black", lwd = 0.3),
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 5)
        #cell_fun = function(j, i, x, y, width, height, fill) {
        #  if(is.na(data.matrix(fc)[i, j]) == FALSE)
        #    grid.text(sprintf("%.1f", data.matrix(fc)[i, j]), x, y, gp = gpar(fontsize = 4))}
        )
dev.off()


## heat map  of ws5a and cp2a together

heatmap_sample_id <- Metadata %>% filter(is.na(subname_mut) | subname_mut %in% c('CP2A','WS5A')) %>%
  filter(!sample_outlier,  !possible_outlier) %>% .$sample_id


sc_ws_cp_ctrl <- ExpDataCPM %>%
  filter(ensemblID %in% ProtGenes$gene_id) %>%
  dplyr::select(heatmap_sample_id) %>%
  as.matrix %>%
  cor

library(ComplexHeatmap)
str_replace

# Create temporal Metadata data.frame
dfMeta <- Metadata %>% dplyr::select(-sample_id, -biological_rep, -source_of_clone, -clone, -homo_het) %>%
  filter(is.na(subname_mut) | subname_mut %in% c('CP2A','WS5A')) %>% 
  filter(!sample_outlier,  !possible_outlier) %>%
  unite("mutation", mutation, subname_mut, remove=TRUE) %>%
  mutate(mutation=gsub("_NA", "", .$mutation)) %>%
  mutate(mutation2=str_replace(.$mutation, "POLG_WS5A","POLGCOMP")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate_at(vars(matches("fraction$")), scales::rescale) %>%
  mutate(lib_size=scales::rescale(lib_size)) %>%
  as.data.frame
rownames(dfMeta) <- Metadata %>% filter(is.na(subname_mut) | subname_mut %in% c('CP2A','WS5A'))%>% 
  filter(!sample_outlier,  !possible_outlier)  %>%
  .$sample_id

library(RColorBrewer)

cols_mutation <- brewer.pal(length(levels(dfMeta$mutation)), "Set2")
names(cols_mutation) <- levels(dfMeta$mutation)

cols_celltype <- brewer.pal(length(levels(dfMeta$cell_type)), "Accent")
names(cols_celltype) <- levels(dfMeta$cell_type)


column_ha = HeatmapAnnotation(CellType=dfMeta$cell_type,
                              Mutation=dfMeta$mutation,
                              Subname_Mutation=dfMeta$subname_mut,
                              na_col="white",
                              col=list(CellType=cols_celltype,
                                       Mutation=cols_mutation
                              ))

Heatmap(sc_ws_cp_ctrl, 
        column_title="Sample correlation in gene expression", 
        #col=viridis::viridis(10), 
        col=colorRampPalette(brewer.pal(10, "RdYlBu"))(256),
        top_annotation=column_ha,
        show_row_names=FALSE, 
        show_column_names=FALSE,
        show_column_dend=FALSE, 
        #show_row_dend=FALSE, 
        #column_dend_height=unit(5,"cm"),
        row_dend_width=unit(5,"cm"))


# PCA over samples exclude outliers

ExpDataCPM_ws_cp_ctrl<-ExpDataCPM %>% dplyr::select(heatmap_sample_id )

dfMeta_ws_cp_ctrl<- dfMeta %>% filter(rownames(.) %in%  heatmap_sample_id)

  ## calculat PCA

sample_pca_ws_cp_ctrl <- prcomp(t(ExpDataCPM_ws_cp_ctrl))

pcameta_ws_cp_ctrl<-dfMeta_ws_cp_ctrl %>%
  unite(celltype_mutation,c('mutation','cell_type'),remove = FALSE)

df_pca<-data.frame(sample_pca_ws_cp_ctrl$x,Group = pcameta_ws_cp_ctrl$celltype_mutation)

  ## add PC percentage 
percentage<-round(sample_pca_ws_cp_ctrl$sdev^2/sum(sample_pca_ws_cp_ctrl$sdev^2)*100,2)
percentage<-paste(colnames(sample_pca_ws_cp_ctrl$x), '(',paste(as.character(percentage)),'%',')',sep='')

  # ggplot
ggplot(df_pca,aes(x=PC1,y=PC2,colour = Group))+
  geom_point()+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),axis.line=element_line(colour='black'))+
  labs(x=percentage[1],y=percentage[2])
  #stat_ellipse(level=0.9,show.legend=F)
  
  # another way of plot
autoplot(sample_pca_ws_cp_ctrl, 
         data = pcameta_ws_cp_ctrl, 
         colour="celltype_mutation")
