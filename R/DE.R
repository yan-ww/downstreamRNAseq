DE <- function(a, condition, FC = 1.5, pval = 0.05, metadata) {
  require(DESeq2)

  # ============ 1. metadata ============
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }


  if (is.null(rownames(metadata)) || any(rownames(metadata) == "")) {
    warning("metadata row names are missing or empty. Attempting to use first column as row names...")
    if (ncol(metadata) >= 1) {
      rownames(metadata) <- metadata[,1]
      metadata <- metadata[,-1, drop = FALSE]
      message("Used first column as row names")
    } else {
      stop("metadata has no row names and no columns to use as row names")
    }
  }


  if (ncol(metadata) != 1) {
    warning("metadata has ", ncol(metadata), " columns. Only the first column will be used as 'group'")
    metadata <- metadata[,1, drop = FALSE]
  }


  colnames(metadata) <- "group"


  if (!all(colnames(a) %in% rownames(metadata))) {
    stop("Sample names in count matrix do not match row names in metadata")
  }


  a <- a[, rownames(metadata), drop = FALSE]

  # ============ 2. count matrix============

  if (!is.matrix(a) && !is.data.frame(a)) {
    stop("a must be a matrix or data frame")
  }


  if (!is.null(rownames(a)) && all(rownames(a) == a[,1])) {
    warning("Row names appear to be duplicated in first column. Removing first column...")
    a <- a[,-1, drop = FALSE]
  }


  if (is.null(rownames(a)) || any(rownames(a) == "")) {
    warning("Row names are missing. Attempting to use first column as row names...")
    if (ncol(a) >= 1) {
      rownames(a) <- a[,1]
      a <- a[,-1, drop = FALSE]
      message("Used first column as row names")
    } else {
      stop("Count matrix has no row names")
    }
  }


  if (any(duplicated(rownames(a)))) {
    warning("Duplicate gene names found. Adding unique suffixes...")
    rownames(a) <- make.unique(rownames(a))
  }


  if (!all(abs(a - round(a)) < 1e-6)) {
    warning("Non-integer values detected in count matrix. Rounding to nearest integer...")
    a <- round(a)
    message("Values have been rounded to integers")
  }


  if (any(a < 0)) {
    stop("Negative values found in count matrix. Raw counts cannot be negative.")
  }


  zero_ratio <- sum(a == 0) / length(a)
  if (zero_ratio > 0.5) {
    warning("More than 50% of counts are zero. Consider using DESeq2's filters or alternative methods.")
  }


  if (is.data.frame(a)) {
    a <- as.matrix(a)
  }


  lib_sizes <- colSums(a)
  if (any(lib_sizes < 1e5)) {
    warning("Some samples have very small library sizes (< 100,000). This may affect DESeq2 results.")
  }

  # ============ 3. DESeq2 ============

  metadata$group <- as.factor(metadata$group)


  if (length(levels(metadata$group)) != 2) {
    stop("Group variable must have exactly 2 levels. Currently has: ",
         paste(levels(metadata$group), collapse = ", "))
  }


  dds <- DESeqDataSetFromMatrix(countData = a,
                                colData = metadata,
                                design = ~ group)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group",
                                   levels(metadata$group)[2],
                                   levels(metadata$group)[1]))
  res <- res[order(res$pvalue), ]
  res$change <- as.factor(ifelse(!is.na(res$padj) & res$padj < pval &
                                   res$log2FoldChange > log2(FC), "up",
                                 ifelse(!is.na(res$padj) & res$padj < pval &
                                          res$log2FoldChange < (-log2(FC)), "down", "not")))

  # ============ 4. output ============
  attr(res, "design_info") <- list(
    contrast = paste(levels(metadata$group)[2], "vs", levels(metadata$group)[1]),
    fc_threshold = FC,
    pval_threshold = pval,
    n_up = sum(res$change == "up", na.rm = TRUE),
    n_down = sum(res$change == "down", na.rm = TRUE),
    n_not = sum(res$change == "not", na.rm = TRUE)
  )


  message("\n========== DESeq2 Analysis Summary ==========")
  message("Contrast: ", attr(res, "design_info")$contrast)
  message("Up-regulated genes: ", attr(res, "design_info")$n_up)
  message("Down-regulated genes: ", attr(res, "design_info")$n_down)
  message("Not significant: ", attr(res, "design_info")$n_not)
  message("============================================\n")

  return(res)
}
