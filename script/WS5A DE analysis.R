## 4 Differential expression analysis

# Load libraries and get annotation
require("DESeq2")
require("ermineR")
require("colorspace")

GenericHumanAnno <- GetAnnoFiles("Generic_human")
View(GenericHumanAnno)

#Group 1
contrast_id <- "NSC-WS5A_vs_NSC-CTRL-all"
contrast_str <- "mutation_POLG_vs_CTRL"
Mt <- Metadata %>% filter(cell_type %in% c('NSC'), (is.na(subname_mut)) | subname_mut == 'WS5A') %>%
  filter(!sample_outlier,  !possible_outlier)
Mt %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

Ct <- Counts[,colnames(Counts) %in% Mt$sample_id]
rownames(Ct) <- Counts$gene_id
Md <- as.formula(paste0("~mutation"))
dds <- DESeq2RUN(Ct, Mt, Md)

res_WS5A_1<-results(dds, alpha = 0.05, format = "DataFrame", independentFiltering = T)
res_WS5A <- GetDESeq2Results(dds, coef=contrast_str, geneNames=geneNames) %>%arrange(padj)
length(res_WS5A[res_WS5A$padj<0.05,1]) #2408 #3726

de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_WS5A, de_file)


#Group 4
contrast_id <- "IPS-WS5A_vs_IPS-CTRL-all"
contrast_str <- "mutation_POLG_vs_CTRL"
Mt <- Metadata %>% filter(cell_type %in% c('IPSC','ESC'), (is.na(subname_mut)) | subname_mut == 'WS5A') %>%
  filter(!sample_outlier,  !possible_outlier)
Mt %>% dplyr::select(sample_id:age, -biological_rep,source_of_clone,sex) %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

Ct <- Counts[,colnames(Counts) %in% Mt$sample_id]
rownames(Ct) <- Counts$gene_id
Md <- as.formula(paste0("~ mutation"))
dds <- DESeq2RUN(Ct, Mt, Md)

res_WS5A_1<-results(dds, alpha = 0.05, format = "DataFrame", independentFiltering = T)
res_WS5A <- GetDESeq2Results(dds, coef=contrast_str, geneNames=geneNames) %>%arrange(padj)
length(res_WS5A[res_WS5A$padj<0.05,1]) #2408 #3726

de_file <- paste0("./Yu_tables/DE_", contrast_id, ".tsv")
write_tsv(res_WS5A, de_file)