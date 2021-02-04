ensembl_length<-read.csv('/home/yu/Postdoc_project/NSC_RNAseq/RNA_seq/metabolic_analysis/ensemble_gene_length.csv')

cpmMatrixFiltered_unlog <- data.frame(gene_id=as.character(countMatrixFiltered$gene_id), Count2CPM(countMatrixFiltered[,-1]))

cpm_matrix_Ens<-left_join(cpmMatrixFiltered_unlog, geneNames %>% dplyr::select(gene_id, gene_name)) %>% 
  mutate(gene_id = gsub('[.][0-9]+$','',.$gene_id))



cpm_matrix_Ense<-left_join(ensembl_length[,c(1,2)],cpm_matrix_Ens,by=c('X'='gene_id'))

## Group ID
NSC_ctrl_id_anbin <- Metadata %>% filter(cell_type == 'NSC', mutation == 'CTRL') %>% filter(clone != 'ALPERS') %>% 
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)

NSC_polg_id_anbin <- Metadata %>% filter(cell_type == 'NSC', mutation == 'POLG') %>% filter(clone != 'ALPERS') %>% 
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)

IPS_ctrl_id_anbin <- Metadata %>% filter(cell_type %in% c('ESC','IPSC'), mutation == 'CTRL') %>% filter(clone != 'ALPERS') %>% 
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)

IPS_polg_id_anbin <- Metadata %>% filter(cell_type %in% c('IPSC'), mutation == 'POLG') %>% filter(clone != 'ALPERS') %>% 
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)

## FPKM
#change count of fragments to RFKM: FPKM ( log transformed)= CPM*1000/gene_length  + 1
ExpData_FPKM<-data.frame(ensemblID=cpm_matrix_Ense$X,GeneSymbol=cpm_matrix_Ense$gene_name,
                         log2(cpm_matrix_Ense[,c(3:(length(colnames(cpm_matrix_Ense))-1))]*1000/cpm_matrix_Ense[,2]+1))

AnbinList<-ExpData_FPKM[ExpData_FPKM$GeneSymbol %in% c('POLA1','POLB','POLG','POLD1','POLE',
                                            'REV3L','POLH','POLQ','POLI','POLK',
                                            'POLL','POLM','POLN','TENT4A',
                                            'REV1','DNTT','PRIMPOL'),]



Meta_no_outlier_anbin<-Metadata %>% filter(!sample_outlier,  !possible_outlier) %>% 
  filter(clone != 'ALPERS''') %>% 
  select(sample_id,cell_type,mutation,clone,biological_rep,source_of_clone)

write.csv(Meta_no_outlier_anbin,'sample_id_anbin.csv')


AnbinList <- AnbinList %>% select(ensemblID,GeneSymbol,Meta_no_outlier_anbin$sample_id)

write.csv(AnbinList,'gene_expression_anbin.csv')


Meta_no_outlier %>%  select(sample_id,)



