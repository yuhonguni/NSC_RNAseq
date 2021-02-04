#Calcualte the RPKM value of all expressed genes, and calculate the mean RPKM of NSCs

library("EDASeq")

ensembl_list <-ExpDataCPM %>%
  mutate(ensemblID = gsub('[.][0-9]+$','',.$ensemblID)) %>% pull(1)

#ensembl_length<-getGeneLengthAndGCContent(ensembl_list, "hsa")
#write.csv(ensembl_length,'ensemble_gene_length.csv')

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


##TPM
rpk<-Count_matrix_Ense[,c(3:(length(colnames(Count_matrix_Ense))-1))]*1000/Count_matrix_Ense[,2]
scale_factor<-colSums(rpk,na.rm = TRUE)/1000000
ExpData_TPM<-data.frame(ensemblID=Count_matrix_Ense$X,GeneSymbol=Count_matrix_Ense$gene_name,
                        rpk/scale_factor)

NSC_ctrl_mean_tpm<-ExpData_TPM %>% dplyr::select(ensemblID,GeneSymbol,NSC_ctrl_id) %>% mutate(mean_rpkm = rowMeans(.[,c(3:13)])) %>% dplyr::select(ensemblID,mean_rpkm)
IPS_ctrl_mean_tpm<-ExpData_TPM  %>% dplyr::select(ensemblID,GeneSymbol,IPS_ctrl_id) %>% mutate(mean_rpkm = rowMeans(.[,c(3:12)])) %>% dplyr::select(ensemblID,mean_rpkm)
NSC_polg_mean_tpm<-ExpData_TPM  %>% dplyr::select(ensemblID,GeneSymbol,NSC_polg_id) %>% mutate(mean_rpkm = rowMeans(.[,c(3:21)])) %>% dplyr::select(ensemblID,mean_rpkm)
IPS_polg_mean_tpm<-ExpData_TPM  %>% dplyr::select(ensemblID,GeneSymbol,IPS_polg_id) %>% mutate(mean_rpkm = rowMeans(.[,c(3:21)])) %>% dplyr::select(ensemblID,mean_rpkm)

table(NSC_ctrl_mean_tpm$mean_rpkm>1)

mean_tpm<-cbind(IPS_ctrl_mean_cpm, IPS_polg_mean_cpm[,2], NSC_ctrl_mean_cpm[,2], NSC_polg_mean_cpm[,2])
colnames(mean_tpm)<-c('ensemblID','IPS_ctrl','IPS_polg','NSC_ctrl','NSC_polg')

write.table(mean_tpm,'/home/yu/Postdoc_project/NSC_RNAseq/RNA_seq/metabolic_analysis/mean_tpm_metabolic_all_group.txt',row.names = F,quote = F, sep='\t')
write.csv(mean_tpm,'/home/yu/Postdoc_project/NSC_RNAseq/RNA_seq/metabolic_analysis/mean_tpm_metabolic_all_group.csv',row.names = F,quote = F)
## import the RPKM to matlab to generate GEM model and extract the metabolic-metabolic interaction tables
#......
