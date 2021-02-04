try1<-Metadata %>% filter(cell_type %in% c('IPSC','ESC'), (is.na(subname_mut))) %>% 
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::select(sample_id,cell_type, sex)

Mito_express_ips_ips<-ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  left_join(en2ENSE[,1:2],by=c('ensemblID'='ENSEMBL')) %>% dplyr::select(ENTREZID,try1$sample_id) %>% 
  filter(ENTREZID %in% mito_id$ENTREZID) %>% mutate(mean = rowSums(.[,2:11])/10)


try2<-Metadata %>% filter(cell_type %in% c('IPSC','ESC'), subname_mut=='Alpers') %>% 
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::select(sample_id,cell_type, sex)

Mito_express_alpers_ips<-ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  left_join(en2ENSE[,1:2],by=c('ensemblID'='ENSEMBL')) %>% dplyr::select(ENTREZID,try2$sample_id) %>% 
  filter(ENTREZID %in% mito_id$ENTREZID) %>% mutate(mean = rowSums(.[,2:4])/3,sdd = sd(.[,2:4]))


try3<-Metadata %>% filter(cell_type %in% c('NSC'), is.na(subname_mut)) %>% 
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::select(sample_id,cell_type, sex)

Mito_express_ips_nsc<-ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  left_join(en2ENSE[,1:2],by=c('ensemblID'='ENSEMBL')) %>% dplyr::select(ENTREZID,try3$sample_id) %>% 
  filter(ENTREZID %in% mito_id$ENTREZID) %>% mutate(mean = rowSums(.[,2:12])/11,sdd = sd(.[,2:12]))


try4<-Metadata %>% filter(cell_type %in% c('NSC'), subname_mut=='Alpers') %>% 
  filter(!sample_outlier,  !possible_outlier) %>% dplyr::select(sample_id,cell_type, sex)

Mito_express_alpers_nsc<-ExpDataCPM %>% mutate(ensemblID = gsub('\\.[0-9]+$','',ExpDataCPM$ensemblID)) %>%
  left_join(en2ENSE[,1:2],by=c('ensemblID'='ENSEMBL')) %>% dplyr::select(ENTREZID,try4$sample_id) %>% 
  filter(ENTREZID %in% mito_id$ENTREZID) %>% mutate(mean = rowSums(.[,2:4])/3)