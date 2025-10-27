library(tidyverse)
library(matrixStats)  # 用于高效矩阵计算

# ====== 1. 定义核心处理函数 ======
process_proteomics_data <- function(file_path, group_name) {
  # 1. 读取数据并处理ID
  df <- read_tsv(file_path) %>%
    mutate(idx = str_replace(idx, "\\.[0-9]+$", ""))
  
  # 2. 添加基因符号
  df <- df %>%
    left_join(gene_annotation, by = c("idx" = "ID")) %>%
    mutate(
      gene_symbol = coalesce(symbol, idx),
      .keep = "unused"
    ) %>%
    distinct(gene_symbol, .keep_all = TRUE)
  # 3. 设置行名为基因符号
  
  
  gene_symbols <- df$gene_symbol
  expr_data <- df %>% dplyr::select(, -gene_symbol)
  expr_data<-as.data.frame(expr_data)
  rownames(expr_data) <- gene_symbols
  
  # 6. 过滤低覆盖度蛋白质(少于70%的人有数据)
  # 7. NA转为0
  numeric_cols <- expr_data  %>% dplyr::select(where(is.numeric)) %>% names()
  filtered_data <- expr_data %>%
    dplyr::filter(rowSums(!is.na(.)) >= ncol(expr_data)*0.7) %>%
    mutate(across(all_of(numeric_cols), ~ replace_na(.x, 0)))
  # 返回结果并验证行名
  
  return(as.data.frame(filtered_data))
}
# ====== 2. 预加载基因注释 ======
gene_annotation <- read_tsv("gencode.v36.annotation.gtf.gene.probemap") %>%
  transmute(
    ID = str_replace(id, "\\.[0-9]+$", ""),
    symbol = gene  # 根据实际列名调整（gene/gene_name/symbol）
  ) %>%
  distinct(ID, .keep_all = TRUE)  # 确保ID唯一性[6](@ref)

# ====== 3. 并行处理肿瘤和对照组 ======
tumor_data <- process_proteomics_data(
  "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
  "Tumor"
)

control_data <- process_proteomics_data(
  "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt",
  "Control"
)
mad_control <- apply(control_data, 1, mad, na.rm=TRUE)
mad_tumor <- apply(tumor_data, 1, mad, na.rm=TRUE)
tumor_data <- tumor_data[which(mad_tumor > quantile(mad_tumor, probs=seq(0, 1, 0.25))[3]),]
control_data <- control_data[which(mad_control > quantile(mad_control, probs=seq(0, 1, 0.25))[3]),]
# ====== 4. 数据整合与保存 ======
# 合并数据集（取基因交集）
common_genes <- intersect(rownames(tumor_data), rownames(control_data))
tumor_data<-tumor_data[common_genes,]
control_data<-control_data[common_genes,]
combined_data <- cbind(
  tumor_data[common_genes, ],
  control_data[common_genes, ]
)

# 保存结果（保留基因符号为行名）
write.csv(combined_data, "processed data/LSCC_combined_proteomics_processed.csv")
write.csv(tumor_data, "processed data/LSCC_tumor_proteomics_processed.csv")
write.csv(control_data, "processed data/LSCC_control_proteomics_processed.csv")

# ====== 5. 质量验证 ======
cat("Final dataset dimensions:\n")
cat("  Tumor:", dim(tumor_data), "\n")
cat("  Control:", dim(control_data), "\n")
cat("  Combined:", dim(combined_data), "\n")

# 验证归一化效果
check_normalization <- function(data) {
  medians <- colMedians(as.matrix(data), na.rm = TRUE)
  cat("Median values range:", range(medians), "\n")
}
check_normalization(combined_data)
