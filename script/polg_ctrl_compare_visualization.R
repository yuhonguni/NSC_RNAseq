
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
