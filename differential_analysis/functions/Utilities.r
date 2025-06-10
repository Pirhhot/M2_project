#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR3/bin/R
library(parallelly)
library(ArchR)
library(ggplot2)
library(dplyr)

CODE_PATH <- '/data/users/tberthom/Project/code'
source(file.path(CODE_PATH, "low_level", "custom_libra_functions.r"))

getLibraDiffAnalysis <- function(proj_archr, diff_annotation, group_forward, group_ref,
                   cell_type = "cell.type",
                   paired_factor = "Patient.number",
                   sample_col = "Sample.number",
                   useMatrix = "PeakMatrix",
                   de_family = "pseudobulk",
                   de_method = "DESeq2",
                   de_type = "Wald",
                   normalization = "log_tp10k",
                   binarization = FALSE,
                   n_threads = ceiling(parallelly::availableCores() / 2)
                   ) {
    
    idx_condition <- BiocGenerics::which(proj_archr@cellColData[,diff_annotation] %in% c(group_ref, group_forward))
    cellsSample <- proj_archr$cellNames[idx_condition]
    proj_subset <- proj_archr[cellsSample, ]

    meta <- proj_subset@cellColData[, c(paired_factor, cell_type, diff_annotation)]
    meta[, paired_factor] <- as.vector(meta[, paired_factor])
    meta[, sample_col] <- as.vector(meta[, sample_col])
    meta[, cell_type] <- as.vector(meta[, cell_type])
    meta[, diff_annotation] <- relevel(factor(meta[, diff_annotation]), ref = group_ref)

    if (useMatrix == "PeakMatrix") {
        type_analysis <- 'scATAC'
    } else if (useMatrix == "GeneScoreMatrix") {
        type_analysis <- 'scRNA' 
    } else {
        message(paste("Matrix type", useMatrix, "unknown. Defaulting type of analysis to scATAC.", sep = " "))
        type_analysis <- 'scATAC' 
    }
    
    peakmatrix <- getMatrixFromProject(
        ArchRProj = proj_subset,
        useMatrix = useMatrix,
        binarize = FALSE,
        verbose = TRUE
    )
    
    if (useMatrix == "PeakMatrix") {
      peaks <- rowRanges(peakmatrix)
      peak_labels <- paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks))
      rownames(peakmatrix) <- peak_labels
    }

    DE <- run_de(assay(peakmatrix),
        meta = meta,
        replicate_col = paired_factor,
        cell_type_col = cell_type,
        label_col = diff_annotation,
        de_family = de_family,
        de_method = de_method,
        de_type = de_type,
        input_type = type_analysis,
        normalization = normalization,
        binarization = binarization,
        n_threads = n_threads
    )
    return(DE)
}

gr_from_gene <- function(gene_names) {
    # Try to parse expected format (chr1:1000-2000)
    parsed <- tryCatch({
      strcapture("chr([0-9XYM]+):([0-9]+)-([0-9]+)",
                 gene_names,
                 proto = list(chr = character(), start = integer(), end = integer()))
    }, error = function(e) {
      stop("Gene names not in expected 'chr1:1000-2000' format")
    })
    
    # Check if parsing was successful
    if(all(is.na(parsed$chr))) {
      stop("Failed to parse genomic coordinates from gene names")
    }
    
    GRanges(seqnames = paste0("chr", parsed$chr),
            ranges = IRanges(start = parsed$start, end = parsed$end))
}
getLibraPseudoBulkAnalysis <- function(proj_archr, diff_annotation, 
                   group_ref = NULL, 
                   group_control = NULL,
                   cell_type = "cell.type",
                   replicate_col = "Patient.number",
                   control_col = NULL,
                   sample_col = "Sample.number",
                   paired_by_replicate = FALSE,
                   useMatrix = "PeakMatrix",
                   de_family = "pseudobulk",
                   de_method = "DESeq2",
                   de_type = "Wald",
                   minCell = 150,
                   impute_matrix = FALSE, 
                   n_threads = ceiling(parallelly::availableCores() / 2)
                   ) { 

  if (!is.numeric(proj_archr@cellColData[,diff_annotation])) {
    if (is.null(group_ref) || is.null(group_control)) {
      stop("Control group or condition group is not stated")
    }
    idx_condition <- BiocGenerics::which(proj_archr@cellColData[,diff_annotation] %in% c(group_ref, group_control))
    cellsSample <- proj_archr$cellNames[idx_condition]
    proj_subset <- proj_archr[cellsSample, ]
  } else {
    idx_condition <- BiocGenerics::which(!is.na(proj_archr@cellColData[,diff_annotation]))
    cellsSample <- proj_archr$cellNames[idx_condition]
    proj_subset <- proj_archr[cellsSample, ]
  }
  if (!is.null(control_col)){
    meta <- proj_subset@cellColData[, c(replicate_col, cell_type, diff_annotation, control_col, sample_col)]
    meta[, control_col] <- as.vector(meta[, control_col])
  } else {
    meta <- proj_subset@cellColData[, c(replicate_col, cell_type, diff_annotation, sample_col)]
  }
  meta[, replicate_col] <- as.vector(meta[, replicate_col]) 
  meta[, sample_col] <- as.vector(meta[, sample_col]) 
  meta[, cell_type] <- as.vector(meta[, cell_type])
  if (!is.numeric(proj_archr@cellColData[,diff_annotation])) {
    meta[, diff_annotation] <- relevel(factor(meta[, diff_annotation]), ref = group_ref)
  } else {
    meta[, diff_annotation] <- as.vector(meta[, diff_annotation])
  }

  if (useMatrix == "PeakMatrix") {
    type_analysis <- 'scATAC'
  } else if (useMatrix == "GeneScoreMatrix") {
    type_analysis <- 'scRNA' 
  } else {
    message(paste("Matrix type", useMatrix, "unknown. Defaulting type of analysis to scRNA.", sep = " "))
    type_analysis <- 'scRNA'
  }
  
  peakmatrix <- getMatrixFromProject(
    ArchRProj = proj_subset,
    useMatrix = useMatrix,
    binarize = FALSE,
    verbose = TRUE
  )

  if (useMatrix == "PeakMatrix") {
    peaks <- rowRanges(peakmatrix)
    peak_labels <- paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks))
    rownames(peakmatrix) <- peak_labels
  } else if (useMatrix == "GeneScoreMatrix") {
    rownames(peakmatrix) <- rowData(peakmatrix)$name
  }

  if (impute_matrix) {
    iW <- getImputeWeights(proj_archr)
    # Impute
    assay(peakmatrix) <- imputeMatrix(assay(peakmatrix), iW) # TO CHECK !!!!!!!!!
  }

  DE <- custom_run_de(assay(peakmatrix),
    meta = meta,
    replicate_col = replicate_col,
    cell_type_col = cell_type,
    label_col = diff_annotation,
    control_col = control_col,
    sample_col = sample_col,
    paired_by_replicate = paired_by_replicate,
    min_cells = minCell,
    de_family = "pseudobulk",
    de_method = de_method,
    de_type = de_type,
    input_type = type_analysis,
    n_threads = n_threads
  )
  return(DE)
}

#peak gr is obtained with getPeakSet
convertLibraOutToSE <- function(df, peak_gr) {
  # Check for required columns first
  required_cols <- c("cell_type", "gene", "avg_logFC", "p_val", "p_val_adj")
  stopifnot(all(required_cols %in% colnames(df)))
  
  # Extract cell type names before transformation
  col_names <- unique(df$cell_type)
  
  # Add index if not already present
  df <- df %>% mutate(idx = row_number())
  
  # Create assay matrices: wide format with cell_type as columns
  df_wide <- df %>%
    select(cell_type, gene, Log2FC = avg_logFC, FDR = p_val_adj, Pval = p_val) %>%
    pivot_wider(names_from = cell_type, values_from = c(Log2FC, FDR, Pval))
  
  # Reconstruct GRanges from gene column
  
  gene_gr <- gr_from_gene(df_wide$gene)
  
  # Find overlaps with peak_gr
  overlaps <- findOverlaps(gene_gr, peak_gr, select = "first")
  
  # Keep only valid matches
  valid_idx <- which(!is.na(overlaps))
  if (length(valid_idx) < length(overlaps)) {
    warning(sprintf("%d ranges did not overlap and were removed.", 
                  length(overlaps) - length(valid_idx)))
  }
  
  # Filter df_wide and gene_gr accordingly
  df_wide <- df_wide[valid_idx, ]
  row_gr <- peak_gr[overlaps[valid_idx]]
  
  # Consistent row names across all objects
  row_names <- as.character(seq_along(row_gr))
  names(row_gr) <- row_names
  
  # Create assay matrices
  make_matrix <- function(prefix) {
    mat <- df_wide %>%
      select(starts_with(prefix)) %>%
      as.matrix()
    colnames(mat) <- sub(paste0("^", prefix, "_"), "", colnames(mat)) # Strip prefix consistently
    rownames(mat) <- row_names  # Use consistent row names
    mat
  }
  
  assays_list <- list(
    Log2FC = make_matrix("Log2FC"),
    FDR    = make_matrix("FDR"),
    Pval   = make_matrix("Pval")
  )
  
  # Construct rowData with consistent row names
  row_df <- as.data.frame(row_gr)[, c("seqnames", "start", "end")]
  rownames(row_df) <- row_names
  
  # Final check: all parts must have the same number of rows
  stopifnot(
    length(row_gr) == nrow(df_wide),
    nrow(assays_list$Log2FC) == nrow(df_wide),
    nrow(row_df) == nrow(df_wide),
    all(dim(assays_list$Log2FC)[2] == sapply(assays_list, ncol)) # Check column counts too
  )
  
  # Build SummarizedExperiment
  se <- SummarizedExperiment(
    assays = assays_list,
    rowData = row_df,
    colData = DataFrame(row.names = col_names)
  )
  
  # Set metadata for ArchR
  metadata(se)$Params$useMatrix <- "PeakMatrix"
  
  return(se)
}


plotVolcanoLibra <- function(df, logfc_col = "avg_logFC", pval_col = "p_val_adj",
                         label_col = "gene", title = "Volcano Plot",
                         logfc_thresh = 0.25, pval_thresh = 0.05) {
  
  # Add significance category
  df <- df %>%
    mutate(
      significance = case_when(
        !!sym(pval_col) < pval_thresh & !!sym(logfc_col) > logfc_thresh ~ "Upregulated",
        !!sym(pval_col) < pval_thresh & !!sym(logfc_col) < -logfc_thresh ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  # Base plot
  p <- ggplot(df, aes(x = !!sym(logfc_col), y = -log10(!!sym(pval_col)))) +
    geom_point(aes(color = significance), alpha = 0.8, size = 2) +
    scale_color_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "Not significant" = "grey70")) +
    theme_minimal(base_size = 14) +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-Value",
      color = "Significance"
    ) +
    geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed", color = "grey50")
  
  # Optional: Label top significant points
  top_genes <- df %>%
    filter(!!sym(pval_col) < pval_thresh & abs(!!sym(logfc_col)) > logfc_thresh) %>%
    arrange(!!sym(pval_col)) %>%
    slice_head(n = 5)
  
  p <- p + ggrepel::geom_text_repel(
    data = top_genes,
    aes(label = !!sym(label_col)),
    size = 3.5,
    max.overlaps = 20
  )
  
  return(p)
}
addCombinedAnnotation <- function(project, annotation1, annotation2){
    combined_anno = paste(annotation1, "by", annotation2, sep = "_")
    project@cellColData[[combined_anno]] <- mapply(paste, 
                                                     as.vector(project@cellColData[[annotation1]]), 
                                                     as.vector(project@cellColData[[annotation2]]), 
                                                     sep = "_")
    unname(project@cellColData[[combined_anno]])
    return(project)
}