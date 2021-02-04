#ensembl_length<-getGeneLengthAndGCContent(ensembl_list, "hsa")
#write.csv(ensembl_length,'ensemble_gene_length.csv')

###CDR related genes Torbjorn

ensembl_length<-read.csv('/home/yu/Postdoc_project/NSC_RNAseq/RNA_seq/metabolic_analysis/ensemble_gene_length.csv')

cpmMatrixFiltered_unlog <- data.frame(gene_id=as.character(countMatrixFiltered$gene_id), Count2CPM(countMatrixFiltered[,-1]))

cpm_matrix_Ens<-left_join(cpmMatrixFiltered_unlog, geneNames %>% dplyr::select(gene_id, gene_name)) %>% 
  mutate(gene_id = gsub('[.][0-9]+$','',.$gene_id))



cpm_matrix_Ense<-left_join(ensembl_length[,c(1,2)],cpm_matrix_Ens,by=c('X'='gene_id'))

## Group ID
NSC_ctrl_id <- Metadata %>% filter(cell_type == 'NSC', mutation == 'CTRL') %>%
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)

NSC_polg_id <- Metadata %>% filter(cell_type == 'NSC', mutation == 'POLG') %>%
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)

IPS_ctrl_id <- Metadata %>% filter(cell_type %in% c('ESC','IPSC'), mutation == 'CTRL') %>%
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)

IPS_polg_id <- Metadata %>% filter(cell_type %in% c('IPSC'), mutation == 'POLG') %>%
  filter(!sample_outlier,  !possible_outlier) %>% pull(1)



## FPKM
#change count of fragments to RFKM: RFKM ( log transformed)= CPM*1000/gene_length  + 1
ExpData_FPKM<-data.frame(ensemblID=cpm_matrix_Ense$X,GeneSymbol=cpm_matrix_Ense$gene_name,
                         log2(cpm_matrix_Ense[,c(3:(length(colnames(cpm_matrix_Ense))-1))]*1000/cpm_matrix_Ense[,2]+1))

## FPKM not log transformed
ExpData_FPKM<-data.frame(ensemblID=cpm_matrix_Ense$X,GeneSymbol=cpm_matrix_Ense$gene_name,
                         cpm_matrix_Ense[,c(3:(length(colnames(cpm_matrix_Ense))-1))]*1000/cpm_matrix_Ense[,2])
#NSC_ctrl_mean_FPKM<-ExpData_FPKM %>% dplyr::select(ensemblID,GeneSymbol,NSC_ctrl_id) %>% mutate(mean_rpkm = rowMeans(.[,c(3:13)])) %>% dplyr::select(ensemblID,mean_rpkm)
#IPS_ctrl_mean_FPKM<-ExpData_FPKM %>% dplyr::select(ensemblID,GeneSymbol,IPS_ctrl_id) %>% mutate(mean_rpkm = rowMeans(.[,c(3:12)])) %>% dplyr::select(ensemblID,mean_rpkm)
#NSC_polg_mean_FPKM<-ExpData_FPKM %>% dplyr::select(ensemblID,GeneSymbol,NSC_polg_id) %>% mutate(mean_rpkm = rowMeans(.[,c(3:21)])) %>% dplyr::select(ensemblID,mean_rpkm)
#IPS_polg_mean_FPKM<-ExpData_FPKM %>% dplyr::select(ensemblID,GeneSymbol,IPS_polg_id) %>% mutate(mean_rpkm = rowMeans(.[,c(3:21)])) %>% dplyr::select(ensemblID,mean_rpkm)

#mean_FPKM<-cbind(IPS_ctrl_mean_FPKM, IPS_polg_mean_FPKM[,2], NSC_ctrl_mean_FPKM[,2], NSC_polg_mean_FPKM[,2])
#colnames(mean_FPKM)<-c('ensemblID','IPS_ctrl','IPS_polg','NSC_ctrl','NSC_polg')
#View(mean_FPKM)
#write.csv(mean_FPKM,'/home/yu/Postdoc_project/NSC_RNAseq/RNA_seq/metabolic_analysis/mean_rpkm_metabolic_all_group.txt',row.names = F,quote = F)

## FPKM for CDR1 CDR2 CDR2L
## CDR2 CDR2L
ExpData_FPKM[ExpData_FPKM$GeneSymbol %in% c('CDR2','CDR2L'),]
ExpData_FPKM[ExpData_FPKM$GeneSymbol %in% c('POLG','POLG2','PRIMPOL','TFAM','TWNK'),]
## CDR1, Gene length of CDR1 1299bp
CDR1_cpm<-cpmMatrixFiltered %>% mutate(gene_id = gsub('[.][0-9]+$','',.$gene_id)) %>% filter(gene_id == 'ENSG00000184258')
CDR1_FPKM<-data.frame(ensemblID='ENSG00000184258',GeneSymbol='CDR1',
                      CDR1_cpm[,c(2:(length(colnames(CDR1_cpm))))]*1000/1299)

cdr<-rbind(ExpData_FPKM[ExpData_FPKM$GeneSymbol %in% c('CDR2','CDR2L'),],CDR1_FPKM)
rownames(cdr)<-cdr$GeneSymbol

polg<-ExpData_FPKM[ExpData_FPKM$GeneSymbol %in% c('POLG','POLG2','PRIMPOL','TFAM','TWNK'),]
rownames(polg)<-polg$GeneSymbol

cdr<-t(cdr[,3:(length(colnames(cdr)))])
sample_id<-data.frame(rownames(cdr))
cdr<-cbind(cdr,sample_id)
colnames(cdr)[4]<-'sample_id'

polg<-t(polg[,3:(length(colnames(polg)))])
sample_id<-data.frame(rownames(polg))
polg<-cbind(polg,sample_id)
colnames(polg)[6]<-'sample_id'



polg<-Metadata %>% left_join(polg, by = c('sample_id','sample_id')) %>% 
  filter(!sample_outlier,  !possible_outlier) %>% 
  dplyr::select(sample_id,cell_type,mutation,clone,TWNK,TFAM,POLG,PRIMPOL,POLG2) %>% 
  as.data.frame(.) 

write.table(polg,'POLG.txt',sep=' ',row.names = F)

polg %>%  mutate(idvd=substr(.$clone,1,(nchar(.$clone)-1))) %>% 
  group_by(idvd) %>% summarise(TWNK_mean = mean(TWNK),
                               TFAM_mean = mean(TFAM),
                               POLG_mean = mean(POLG),
                               POLG2_mean = mean(POLG2),
                               PRIMPOL_mean = mean(PRIMPOL))

