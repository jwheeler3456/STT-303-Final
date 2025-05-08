if (!requireNamespace("TCGAbiolinks")) install.packages("BiocManager"); BiocManager::install("TCGAbiolinks")
if (!requireNamespace("tidyverse")) install.packages("tidyverse")
if (!requireNamespace("survival")) install.packages("survival")
if (!requireNamespace("survminer")) install.packages("survminer")
install.packages("pheatmap")

library(TCGAbiolinks)
library(tidyverse)
library(survival)
library(survminer)
library(pheatmap)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
data <- GDCprepare(query)

clinical <- GDCquery_clinic("TCGA-BRCA", type = "clinical")

counts <- as.data.frame(SummarizedExperiment::assay(data))

gene_info <- as.data.frame(SummarizedExperiment::rowData(data))

gene_info$gene_id_clean <- gsub("\\..*", "", gene_info$gene_id)
counts$gene_id_clean <- gsub("\\..*", "", rownames(counts))

counts_with_names <- counts %>%
  left_join(gene_info[, c("gene_id_clean", "gene_name")], by = "gene_id_clean")

genes_of_interest <- c("BRCA1", "TP53", "ESR1")
filtered_counts <- counts_with_names %>%
  filter(gene_name %in% genes_of_interest)

expression_cols <- colnames(counts)

expression_data <- as.data.frame(t(filtered_counts[, expression_cols]))

colnames(expression_data) <- filtered_counts$gene_name

expression_data$sample <- substr(rownames(expression_data), 1, 12)

clinical_clean <- clinical %>%
  select(submitter_id, vital_status, days_to_death, days_to_last_follow_up, ajcc_pathologic_stage) %>%
  mutate(
    survival_time = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death),
    status = ifelse(vital_status == "Alive", 0, 1)
  ) %>%
  rename(tumor_stage = ajcc_pathologic_stage)

merged_data <- expression_data %>%
  left_join(clinical_clean, by = c("sample" = "submitter_id"))

merged_data <- merged_data %>%
  mutate(across(all_of(genes_of_interest), as.numeric))

summary(merged_data)

for (gene in genes_of_interest) {
  cat("\nT-test for", gene, ":\n")
  print(t.test(merged_data[[gene]] ~ merged_data$status))
}

for (gene in genes_of_interest) {
  cat("\nANOVA for", gene, "by tumor stage:\n")
  model <- aov(merged_data[[gene]] ~ merged_data$tumor_stage)
  print(summary(model))
}

clustering_data <- merged_data %>%
  select(sample, all_of(genes_of_interest)) %>%
  na.omit()

set.seed(123)
kmeans_result <- kmeans(clustering_data %>% select(-sample), centers = 2)

clustering_data$kmeans_cluster <- kmeans_result$cluster

merged_data <- merged_data %>%
  distinct(sample, .keep_all = TRUE)

clustering_data <- clustering_data %>%
  distinct(sample, .keep_all = TRUE)

merged_data <- merged_data %>%
  left_join(clustering_data %>% select(sample, kmeans_cluster), by = "sample")

clustering_matrix_scaled <- scale(clustering_data %>% select(-sample, -kmeans_cluster))

pheatmap(clustering_matrix_scaled,
         clustering_distance_rows = dist(clustering_matrix_scaled),
         clustering_method = "complete",
         main = "Standardized Heatmap of Gene Expression")

logit_model <- glm(status ~ BRCA1 + TP53 + ESR1 + tumor_stage,
                   data = merged_data, family = "binomial")
summary(logit_model)

surv_obj <- Surv(time = merged_data$survival_time, event = merged_data$status)
km_fit <- survfit(surv_obj ~ kmeans_cluster, data = merged_data)

ggsurvplot(km_fit, data = merged_data, pval = TRUE,
           title = "Survival by Clusters")

write.csv(merged_data, "TCGA_BRCA_Processed.csv", row.names = FALSE)
