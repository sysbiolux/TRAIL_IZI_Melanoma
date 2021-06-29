# Aurélien Ginolhac DLSM, University of Luxembourg
# 2020 MIT License.

#' @param  .data is a SummarizedExperiment
#' @param   meta is the meta data from the PCA pre-processing
#' @param   fc feactureCounts object, list of counts and annotations
#' @param   .cell a string (relevant is A375 or WM1346)
#' @return a subset of SummarizedExperiment
subset_dds <- function(.data, meta, fc, .cell = NULL) {
  stopifnot(!is.null(.cell))
  meta %>%
    filter(cell == .cell) %>%
    pull(name) -> to_keep

  meta %>%
    filter(cell == .cell) %>%
    remove_rownames() %>%
    column_to_rownames(var = "name") %>%
    select(condition) %>%
    mutate(condition = factor(str_replace(condition, "-", "_"))) -> col_data

  sub_dds <- DESeqDataSetFromMatrix(countData = counts(.data)[, to_keep],
                                    colData   = col_data,
                                    design    = ~ condition)
  # annotate as in https://support.bioconductor.org/p/62374/ Mike Love
  select(fc$annotation, GeneID, gene_name) %>%
    semi_join(tibble(GeneID = rownames(sub_dds))) -> features
  stopifnot(all(rownames(sub_dds) == features$GeneID))
  mcols(sub_dds) <- cbind(mcols(sub_dds), gene_name = features$gene_name)
  sub_dds <- estimateSizeFactors(sub_dds)
  # remove uninformative genes
  sub_dds <- sub_dds[rowSums(counts(sub_dds, normalized = FALSE)) > 10, ]
  # normalization and preprocessing
  sub_dds <- DESeq(sub_dds, parallel = TRUE)
  sub_dds
}

sub_plotPCA <- function(.data, .cell) {
  vst_counts <- vst(.data, blind = FALSE)
  df <- plotPCA(vst_counts, intgroup = c("condition"), returnData = TRUE)
  percentVar <- round(100 * attr(df, "percentVar"))
  df %>%
    mutate(IZI = if_else(str_detect(condition, "IZI"), TRUE, FALSE),
           TRAIL = case_when(
             str_detect(condition, "^p") ~ "parental",
             str_detect(condition, "^c") ~ "conditioned",
             TRUE ~ "NA"),
           cell = case_when(
             str_detect(condition,"A375") ~ "A375",
             str_detect(condition, "Malme") ~ "Malme3M",
             str_detect(condition, "WM1346") ~ "WM1346")) %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(aes(color = TRAIL, shape = IZI), size = 3, alpha = .8) +
    ggrepel::geom_text_repel(aes(label = name),
                             min.segment.length = .1, show.legend = FALSE) +
    labs(title = .cell,
         x = paste0("PC1: ", percentVar[1], "% variance"),
         y = paste0("PC2: ", percentVar[2], "% variance")) +
    theme_ipsum_rc(12)
}

sub_results <- function(.data, .coef) {
  stopifnot(!is.null(.coef))
  # if we don't have the right contrast
  if (!.coef %in% resultsNames(.data)) {
    # relevel by the desired reference
    my_ref <- stringr::str_split(.coef, "_vs_", simplify = TRUE)[1, 2]
    .data$condition <- relevel(.data$condition, my_ref)
    # redo the test as https://support.bioconductor.org/p/123247/
    .data <- nbinomWaldTest(.data, quiet = TRUE)
    # now it should work, otherwise fail
    stopifnot(.coef %in% resultsNames(.data))
  }
  message(paste("contrasting", stringr::str_split(.coef, "_vs_", simplify = TRUE)[1, 1],
                "versus",
                stringr::str_split(.coef, "_vs_", simplify = TRUE)[1, 2]))
  res <- results(.data, name = .coef, parallel = parallel)
  # shrink fold changes for low expressed genes
  res <- lfcShrink(.data, coef = .coef, res = res, type = "apeglm")
  # annotate with gene names
  stopifnot(all(row.names(res) == row.names(rowData(.data))))
  res$gene_name <- rowData(.data)$gene_name
  # sort by p-value
  res <- res[order(res$padj), ]
  res
}


#' @param  .data is a tibble of DESeq2 results
two_volcanoes <- function(.data, .title = NULL, fc_thr = 3) {
  stopifnot(!is.null(.title))
  .data %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(contrast = fct_relevel(factor(contrast), "pIZI50_vs_p")) -> df
  df %>%
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_hex(bins = 90) +
    geom_point(data = subset(df, padj < 0.1 & abs(log2FoldChange) >= 1),
               colour = "tomato1") +
    # ggrepel::geom_text_repel(data = subset(df,
    #                                        !is.na(gene_name) &
    #                                          abs(log2FoldChange) >= fc_thr & -log10(padj) > 10),
    #                          colour = "tomato3", aes(label = gene_name)) +
    scale_fill_gradient(low = "grey80", high = "grey30") +
    theme_minimal(14) +
    labs(title = .title,
         y = bquote(-log[10] * " ("*italic(padj)*")"),
         caption = paste0("red dots: adj-pval < 0.1 and abs(lFC) >= 1
       gene names: abs(log2FoldChange) >= ", fc_thr, " and -log10(adj-pval) > 10")) +
    facet_wrap(~ contrast)
}

