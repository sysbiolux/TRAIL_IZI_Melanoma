subset_dds <- function(.data, .cell = NULL) {
  stopifnot(!is.null(.cell))
  pca_df_meta %>%
    filter(cell == .cell) %>%
    pull(name) -> to_keep

  pca_df_meta %>%
    filter(cell == .cell) %>%
    column_to_rownames(var = "name") %>%
    select(condition) %>%
    mutate(condition = factor(str_replace(condition, "-", "_"))) -> col_data

  sub_dds <- DESeqDataSetFromMatrix(countData = count_mat[, to_keep],
                                    colData   = col_data,
                                    design    = ~ condition)
  # annotate as in https://support.bioconductor.org/p/62374/ Mike Love
  stopifnot(all(rownames(sub_dds) == fc$annotation$GeneID))
  mcols(sub_dds) <- cbind(mcols(sub_dds), gene_name = fc$annotation$gene_name)
  sub_dds <- estimateSizeFactors(sub_dds)
  # remove uninformative genes
  sub_dds <- sub_dds[ rowSums(counts(sub_dds, normalized = FALSE)) > 10, ]
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

two_volcanoes <- function(.data, .title = NULL, fc_thr = 3) {
  stopifnot(!is.null(.title))
  .data %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_hex(bins = 90) +
    geom_point(data = subset(.data, padj < 0.1 & abs(log2FoldChange) >= 1),
               colour = "tomato1") +
    ggrepel::geom_text_repel(data = subset(.data,
                                           !is.na(gene_name) &
                                             abs(log2FoldChange) >= fc_thr & -log10(padj) > 10),
                             colour = "tomato3", aes(label = gene_name)) +
    scale_fill_gradient(low = "grey80", high = "grey30") +
    theme_minimal(14) +
    labs(title = .title,
         y = bquote(-log[10] * " ("*italic(padj)*")"),
         caption = paste0("red dots: adj-pval < 0.1 and abs(lFC) >= 1
       gene names: abs(log2FoldChange) >= ", fc_thr, " and -log10(adj-pval) > 10")) +
    facet_wrap(~ contrast) +
    theme_ipsum_rc()
}

two_lfc <- function(.data, .title = NULL) {
  stopifnot(!is.null(.title))
  .data %>%
    pivot_wider(id_cols = c(gene_id, baseMean, gene_name),
                names_from = contrast,
                values_from = c(log2FoldChange:padj)) %>%
    ggplot(aes(x = log2FoldChange_conditioned_vs_parental,
               y = log2FoldChange_IZIc_vs_IZIp)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggpointdensity::geom_pointdensity() +
    ggrepel::geom_text_repel(data = function(x) filter(x, abs(log2FoldChange_conditioned_vs_parental - log2FoldChange_IZIc_vs_IZIp) > 5),
                             aes(label = gene_name),
                             max.overlaps = 8,
                             min.segment.length = 0.1,
                             segment.size = 0.5,
                             segment.curvature = -0.1,
                             segment.ncp = 3, segment.angle = 20,
                             segment.linetype = 1, size = 3,
                             color = "white",     # text color
                             bg.color = "grey30", # shadow color
                             bg.r = 0.1,         # shadow radius
                             arrow = arrow(length = unit(0.015, "npc"))) +
    theme_ipsum_rc() +
    scale_colour_viridis_b() +
    labs(title = .title,
         caption = "gene_names indicated when absolute difference
         between lFC is greater than 5")
}

select_gene_unite <- function(.data, column) {

  select(.data, gene = {{column}}) %>%
    filter(!is.na(gene)) %>%
    left_join(select(IZI, gene_id, gene = gene_name)) %>%
    distinct() %>%
    filter(!is.na(gene_id)) %>%
    rowwise() %>%
    mutate(counts = list(plotCounts(dds, gene = gene_id, returnData = TRUE) %>%
                           filter(str_detect(condition, "Malme", negate = TRUE)))) %>%
    unnest(counts) %>%
    mutate(IZI = if_else(str_detect(condition, "IZI"), TRUE, FALSE),
           TRAIL = case_when(
             str_detect(condition, "^p") ~ "parental",
             str_detect(condition, "^c") ~ "conditioned",
             TRUE ~ "NA"),
           cell = case_when(
             str_detect(condition,"A375") ~ "A375",
             str_detect(condition, "Malme") ~ "Malme3M",
             str_detect(condition, "WM1346") ~ "WM1346",
             TRUE ~ "NA")) %>%
    unite("id", c(TRAIL, IZI), remove = FALSE)
}


gg_games_howell <- function(.data) {
  label_gene <- unique(.data$gen)
  label_cell <- unique(.data$cellule)
  # if all counts are 0, add a jittering
  if (unique(.data$log10_count) == log10(0.5)) {
    .data$log10_count <-  .data$log10_count + runif(nrow(.data), max = 0.1)
  }
  ggstatsplot::ggbetweenstats(.data,
                              x = condition,
                              y = log10_count,
                              plot.type = "box",
                              type = "p",
                              title = paste(label_cell, label_gene),
                              results.subtitle = NULL,
                              grouping.var = id,
                              sample.size.label = FALSE,
                              p.adjust.method = "fdr") +
    labs(caption = NULL, x = NULL) +
    theme_bw(10) +
    theme(legend.position = "none") +
    scale_y_continuous(expand = expand_scale(mult = c(0.05, 0.15)))
}

plot_multipage <- function(gglist, nrow = 4L, ncol = 2L, output_file = NULL) {
  npages <- ceiling(length(gglist) / (nrow * ncol))
  npanels <- nrow * ncol
  res <- vector(mode = "list", length = npages)
  start <- 1L
  for (i in seq(1, length(gglist), npanels)) {
    end <- ifelse(i + (npanels - 1) > length(gglist), length(gglist), i + (npanels - 1))
    print(paste("page", sprintf("%2d", start), ":", sprintf("%2d", i), sprintf("%2d", end)))
    res[[start]] <- cowplot::plot_grid(plotlist = gglist[i:end], ncol = ncol)
    start <- start + 1
  }
  # print page as pdf in a temp folder
  plot_dir <- fs::path_temp()
  fs::dir_create(plot_dir)
  iwalk(res, ~ ggsave(fs::file_temp(pattern = paste0(sprintf("%02i", .y), "cowplot"),
                                    tmp_dir = plot_dir,
                                    ext = "pdf"),
                      plot = .x,
                      width = 210, height = 297, units = "mm")
  )
  qpdf::pdf_combine(input = fs::dir_ls(plot_dir, glob = "*.pdf"),
                    output = output_file)
  fs::file_delete(plot_dir)
  #unlink(plot_dir, recursive = TRUE)
}

plot_heatmap <- function(mat, dk_mat, gene_size = 6, ...) {
  stl <- RColorBrewer::brewer.pal(4, "Set1")
  Heatmap(mat,
          name = "log2FC",
          show_column_dend = FALSE,
          show_row_dend = FALSE,
          col = circlize::colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
          row_names_gp = gpar(fontsize = gene_size),
          column_names_gp = gpar(fontsize = 6),
          column_names_centered = FALSE,
          top_annotation = HeatmapAnnotation(cells = factor(str_extract(colnames(dk_mat), "[:alnum:]+(?=_)")),
                                             contrast = factor(str_extract(colnames(dk_mat), "[:alnum:]+_vs_[:alnum:]+")),
                                             col = list(cells = c("A375" = "blueviolet", "WM1346" = "gold2"),
                                                        contrast = c(c_vs_p = stl[1], cIZI50_vs_c = stl[2],
                                                                     IZIc_vs_IZIp = stl[3], pIZI50_vs_p = stl[4]))),
          column_title_gp = gpar(fontsize = 10, fontface = "bold"), ...)
}