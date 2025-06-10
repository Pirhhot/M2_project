#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR4/bin/R
library(Seurat)
library(dplyr)
library(purrr)
library(ArchR)
library(DESeq2)
library(edgeR)
# the purpose of this file is to recreate Libra package functions for my goals
# I have modified the functions to add the following features 
# 1) control_col parameter : (optional) indicates a column in meta 
# corresponding to a factor we want to control for our data
# in the GLM
# 2) paired option : to precise if we want the analysis to be paired on the replicate factor or not

#adapted from https://github.com/neurorestore/Libra
check_inputs = function(input,
                        meta,
                        replicate_col = 'replicate',
                        cell_type_col = 'cell_type', 
                        label_col = 'label',  
                        control_col = NULL, 
                        sample_col = "Sample.number") {

  # extract cell types and label from metadata
  if ("Seurat" %in% class(input)) {
    # confirm Seurat is installed
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("install \"Seurat\" R package for Libra compatibility with ",
           "input Seurat object", call. = FALSE)
    }
    meta = input@meta.data %>%
      droplevels()
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    if (!is.factor(meta[[label_col]])) {
      labels = meta[[label_col]]
    } else {
      labels = as.character(meta[[label_col]])
    }
    cell_types = as.character(meta[[cell_type_col]])
    expr = Seurat::GetAssayData(input, slot = 'counts')
  } else if ("cell_data_set" %in% class(input)) {
    # confirm monocle3 is installed
    if (!requireNamespace("monocle3", quietly = TRUE)) {
      stop("install \"monocle3\" R package for Libra compatibility with ",
           "input monocle3 object", call. = FALSE)
    }
    meta = monocle3::pData(input) %>%
      droplevels() %>%
      as.data.frame()
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    if (!is.factor(meta[[label_col]])) {
      labels = meta[[label_col]]
    } else {
      labels = as.character(meta[[label_col]])
    }
    cell_types = as.character(meta[[cell_type_col]])
    expr = monocle3::exprs(input)
  } else if ("SingleCellExperiment" %in% class(input)){
    # confirm SingleCellExperiment is installed
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("install \"SingleCellExperiment\" R package for Libra ",
           "compatibility with input SingleCellExperiment object",
           call. = FALSE)
    }
    meta = SummarizedExperiment::colData(input) %>%
      droplevels() %>%
      as.data.frame()
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    if (!is.factor(meta[[label_col]])) {
      labels = meta[[label_col]]
    } else {
      labels = as.character(meta[[label_col]])
    }
    cell_types = as.character(meta[[cell_type_col]])
    expr = SummarizedExperiment::assay(input)
  } else if ("Signac" %in% class(input)){
    # confirm Signac is installed
    if (!requireNamespace("Signac", quietly = TRUE)) {
      stop("install \"Signac\" R package for Libra compatibility with ",
           "input Signac object", call. = FALSE)
    }
    meta = input@meta.data %>%
      droplevels()
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    if (!is.factor(meta[[label_col]])) {
      labels = meta[[label_col]]
    } else {
      labels = as.character(meta[[label_col]])
    }
    cell_types = as.character(meta[[cell_type_col]])
    expr = Signac::GetAssayData(input, slot = 'peaks')
  } else if ("ArchR" %in% class(input)){
    # confirm ArchR is installed
    if (!requireNamespace("ArchR", quietly = TRUE)) {
      stop("install \"ArchR\" R package for Libra compatibility with ",
           "input ArchR object", call. = FALSE)
    }
    meta = data.frame(getCellColData(input)) %>%
      droplevels()
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    if (!is.factor(meta[[label_col]])) {
      labels = meta[[label_col]]
    } else {
      labels = as.character(meta[[label_col]])
    }
    cell_types = as.character(meta[[cell_type_col]])
    expr = getMatrixFromProject(input, useMatrix='PeakMatrix')
  } else if ("snap" %in% class(input)){
    # confirm ArchR is installed
    if (!requireNamespace("SnapATAC", quietly = TRUE)) {
      stop("install \"SnapATAC\" R package for Libra compatibility with ",
           "input SnapATAC object", call. = FALSE)
    }
    meta = data.frame(input@meta.data) %>%
      droplevels()
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    if (!is.factor(meta[[label_col]])) {
      labels = meta[[label_col]]
    } else {
      labels = as.character(meta[[label_col]])
    }
    cell_types = as.character(meta[[cell_type_col]])
    expr = input@pmat
  } else {
    # check if input is sparse matrix or numberic matrix/df
    valid_input = is(input, 'sparseMatrix') ||
      tester::is_numeric_matrix(input) ||
      tester::is_numeric_dataframe(input)
    if (!valid_input)
      stop("input must be Seurat, monocle, sparse matrix, numeric matrix, or ",
           "numeric data frame")
    if (is.null(meta))
      stop("input matrix must be accompanied by a metadata table")
    expr = input
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    if (!is.factor(meta[[label_col]])) {
      labels = meta[[label_col]]
    } else {
      labels = as.character(meta[[label_col]])
    }
    cell_types = as.character(meta[[cell_type_col]])
  }
  
  # check dimensions are non-zero
  if (length(dim(expr)) != 2 || !all(dim(expr) > 0)) {
    stop("expression matrix has at least one dimension of size zero")
  }

  # check dimensions match
  n_cells1 = nrow(meta)
  n_cells2 = ncol(expr)
  if (n_cells1 != n_cells2) {
    stop("number of cells in metadata (", n_cells1, ") does not match number ",
         "of cells in expression (", n_cells2, ")")
  }

  # check at least two labels
  if (n_distinct(labels) == 1) {
    stop("only one label provided: ", unique(labels))
  }
  if (n_distinct(labels) > 2) {
    warning("more than 2 labels provided: ", unique(labels))
  }

  # check for missing labels or cell types
  if (any(is.na(labels))) {
    stop("labels contain ", sum(is.na(labels)), "missing values")
  }
  if (any(is.na(cell_types))) {
    stop("cell types contain ", sum(is.na(cell_types)), "missing values")
  }
  if (!is.null(replicate_col) && any(is.na(replicates))) {
    stop("replicates contain ", sum(is.na(replicates)), "missing values")
  }
  if (!is.null(control_col) && any(is.na(meta[[control_col]]))) {
    stop("Control factor column contain ", sum(is.na(meta[[control_col]])), " missing values")
  }
  

  # check for missing replicates
  if (!is.null(replicate_col) && is.null(replicates)) {
    stop("metadata does not contain replicate information")
  }

  # remove missing values
  missing = is.na(expr)
  if (any(missing)) {
    stop("matrix contains", sum(missing), "missing values")
  }
  
  ## clean up the meta data
  if (!is.null(replicate_col) && !is.null(control_col)) {
    meta %<>% as.data.frame() %>%
      mutate(cell_barcode = rownames(meta),
             replicate = meta[[replicate_col]],
             cell_type = meta[[cell_type_col]], 
             control = meta[[control_col]],
             label = meta[[label_col]],
             sample = meta[[sample_col]])
  } else if (!is.null(replicate_col)) {
    meta %<>% as.data.frame() %>%
      mutate(cell_barcode = rownames(meta),
             cell_type = meta[[cell_type_col]],
             label = meta[[label_col]], 
             replicate = meta[[replicate_col]],
             sample = meta[[sample_col]]) 
  } else if (!is.null(control_col)) {
    meta %<>% as.data.frame() %>%
      mutate(cell_barcode = rownames(meta),
             cell_type = meta[[cell_type_col]],
             label = meta[[label_col]], 
             control = meta[[control_col]],
             sample = meta[[sample_col]]) 
  } else {
    meta %<>% as.data.frame() %>%
      mutate(cell_barcode = rownames(meta),
             cell_type = meta[[cell_type_col]], 
             label = meta[[label_col]],
             sample = meta[[sample_col]])
  }

  # keep factors if present in meta
  if (!is.factor(meta$label) && !is.numeric(meta$label)) {
    meta$label <- as.factor(meta$label)
  }

  if (!is.factor(meta$replicate)) {
    meta$replicate <- as.factor(meta$replicate)
  }
  if (!is.factor(meta$cell_type)) {
    meta$cell_type <- as.factor(meta$cell_type)
  }
  
  # make sure meta contains row names and is a data frame
  rownames(meta) <- colnames(expr)
  meta <- as.data.frame(meta)
  to_return <- list(
    expr = expr,
    meta = meta
  )
  return(to_return)
}

custom_to_pseudobulk <- function(input,
                         meta = NULL,
                         replicate_col = 'replicate',
                         cell_type_col = 'cell_type',
                         label_col = 'label',
                         sample_col = "Sample.number",
                         min_cells = 150,
                         min_reps = 2,
                         min_features = 0,
                         external = T) {
  if (external) {
    # first, make sure inputs are correct
    inputs = check_inputs(
      input,
      meta = meta,
      replicate_col = replicate_col,
      cell_type_col = cell_type_col,
      label_col = label_col, 
      control_col = control_col,
      sample_col = sample_col)
    expr = inputs$expr
    meta = inputs$meta
  } else {
    expr = input
  }

  # convert to characters
  meta %<>% mutate(replicate = as.character(replicate),
                   cell_type = as.character(cell_type),
                   label = as.character(label),
                   sample = as.character(sample))
  
  # keep only cell types with enough cells
  keep = meta %>%
    dplyr::count(cell_type, label) %>%
    group_by(cell_type) %>%
    filter(all(n >= min_cells)) %>%
    pull(cell_type) %>%
    unique()
  
  #filter meta to remove samples/cell typee without enough cells
  keep_sample <- meta %>%
    dplyr::count(cell_type, sample) %>%
    filter(n >= min_cells) %>%
    select(cell_type, sample)

  meta <- meta %>%
    semi_join(keep_sample, by = c("cell_type", "sample"))
  #create colnames
  meta$group <- paste(meta$replicate, meta$label, meta$sample, sep = ":")

  global_group_info <- meta %>%
    select(group, replicate, label, sample) %>%
    distinct() %>%
    mutate(replicate_label = paste(replicate, label, sep = ":"))

  rep_label_counts <- global_group_info %>%
    count(replicate_label, name = "n")
      
  global_group_info <- global_group_info %>%
    left_join(rep_label_counts, by = "replicate_label") %>%
    mutate(final_name = ifelse(n == 1, replicate_label, paste(group, sep = ""))) %>%
    select(group, final_name)

  # process data into gene x replicate x cell_type matrices
  pseudobulks <- keep %>%
    map(~ {
      print(.)
      cell_type = .
      meta0 = meta %>% filter(cell_type == !!cell_type)
      expr0 = expr %>% magrittr::extract(, meta0$cell_barcode)
      # catch cell types without replicates or conditions
      if (n_distinct(meta0$label) < 2)
        return(NA)
      replicate_counts = distinct(meta0, label, replicate) %>%
        group_by(label) %>%
        summarise(replicates = n_distinct(replicate)) %>%
        pull(replicates)
      if (any(replicate_counts < min_reps))
        return(NA)
      
      # process data into gene X replicate X cell_type matrice (X sample if there are 2 samples in one replicate X cell_type label combination)

      mm = model.matrix(~ 0 + group, data = meta0)
      mat_mm = expr0 %*% mm
      keep_genes = rowSums(mat_mm > 0) >= min_features 
      mat_mm = mat_mm[keep_genes, ] %>% as.data.frame()
      mat_mm %<>% as.data.frame()

      group_info <- global_group_info %>% filter(group %in% meta0$group)
      
      colnames(mat_mm) <- group_info$final_name[match(colnames(mat_mm), paste0("group", group_info$group))]

      # drop empty columns
      keep_samples = colSums(mat_mm) > 0
      mat_mm %<>% magrittr::extract(, keep_samples)
      return(mat_mm)
    }) %>%
    setNames(keep)
  
  # drop NAs
  pseudobulks %<>% magrittr::extract(!is.na(.))
  
  # also filter out cell types with no retained genes
  min_dim = map(pseudobulks, as.data.frame) %>% map(nrow)
  pseudobulks %<>% magrittr::extract(min_dim > 1)
  
  # also filter out types without replicates
  min_repl = map_int(pseudobulks, ~ {
    # make sure we have a data frame a not a vector
    tmp = as.data.frame(.)
    targets = data.frame(group_sample = colnames(tmp)) %>%
      mutate(group_parts = strsplit(group_sample, ":")) %>%
      mutate(group = sapply(group_parts, function(x) x[2])) %>%
      select(-group_parts) 
    if (n_distinct(targets$group) == 1)
      return(as.integer(0))
    min(table(targets$group))
  })
  pseudobulks %<>% magrittr::extract(min_repl >= min_reps)
  return(pseudobulks)
}


custom_pseudobulk_de <- function(input,
                                 meta = NULL,
                                 replicate_col = "replicate",
                                 cell_type_col = "cell_type",
                                 paired_by_replicate = FALSE,
                                 control_col = NULL,
                                 label_col = "label",
                                 sample_col = "Sample.number",
                                 min_cells = 150,
                                 min_reps = 2,
                                 min_features = 0,
                                 de_family = "pseudobulk",
                                 de_method = "edgeR",
                                 de_type = "LRT") {
  # check args
  if (de_method == "limma") {
    if (de_type != "voom") {
      # change default type to use
      de_type <- "trend"
    }
  }

  # get the pseudobulks list
  pseudobulks <- custom_to_pseudobulk(
    input = input,
    meta = meta,
    replicate_col = replicate_col,
    cell_type_col = cell_type_col,
    label_col = label_col,
    sample_col = sample_col,
    min_cells = min_cells,
    min_reps = min_reps,
    min_features = min_features,
    external = F
  )

  results <- purrr::imap(pseudobulks, function(x, current_cell_type) {
    # create targets matrix
    targets <- data.frame(group_sample = colnames(x)) %>%
      mutate(group_parts = strsplit(group_sample, ":")) %>%
      mutate(group = sapply(group_parts, function(x) x[2])) %>%
      mutate(replicate = sapply(group_parts, function(x) x[1])) %>%
      select(-group_parts)
    ## optionally, carry over factor levels from entire dataset
    if (is.factor(meta[[label_col]])) {
      targets$group %<>% factor(levels = levels(meta[[label_col]]))
    }
    if (is.numeric(meta[[label_col]])) {
      targets$group <- as.numeric(targets$group)
    }
    if (n_distinct(targets$group) > 2 && !is.numeric(meta[[label_col]])) {
      warning("more than two groups in categorical label column")
      return(NULL)
    }
    # create design
    if (!is.null(control_col)) { 
      grouped_data_for_join <- meta %>%
        mutate(
          replicate = as.character(replicate),
          cell_type = as.character(cell_type),
          label     = as.character(label),
          sample    = as.character(sample),
          replicate_label = paste(replicate, label, sep = ":")
        ) %>%
        group_by(replicate_label) %>%
        mutate(
            group_sample = if (n_distinct(sample) == 1) replicate_label else paste(replicate_label, sample, sep = ":")
          ) %>%
        ungroup() %>%
        select(group_sample, cell_type, control) %>%
        distinct(group_sample, cell_type, control, .keep_all = TRUE)

      tmp_test <- grouped_data_for_join %>%
        select(group_sample, cell_type) %>%
        distinct(group_sample, cell_type, .keep_all = TRUE)
        
      if (nrow(grouped_data_for_join) == nrow(tmp_test)) {
        grouped_data_for_join %<>% 
          filter(cell_type == current_cell_type) %>%
          select(group_sample, control) %>%
          distinct(group_sample, control)
        targets %<>% left_join(grouped_data_for_join, by = join_by(group_sample))
      } else {
        stop("Control column should have a unique value for each sample/diff_annotation/cell_type combination.")
      }

      if (!is.factor(targets$control) && !is.numeric(targets$control)) {
        targets$control <- factor(targets$control)
      }
      if (paired_by_replicate) {
        design <- model.matrix(~ replicate + control + group, data = targets)
      } else {
        design <- model.matrix(~ control + group, data = targets) 
      }
      
    } else {
      if (paired_by_replicate) {
        design <- model.matrix(~ replicate + group, data = targets)
      } else {
        design <- model.matrix(~ group, data = targets)
      }
    }
    
    message("Design used :")
    message(paste(capture.output(print(design)), collapse = "\n"))
    DE <- switch(de_method,
      edgeR = {
        tryCatch(
          {
            y <- DGEList(counts = x, group = targets$group) %>%
              calcNormFactors(method = "TMM") %>%
              estimateDisp(design)
            test <- switch(de_type,
              QLF = {
                fit <- glmQLFit(y, design)
                test <- glmQLFTest(fit, coef = -1)
              },
              LRT = {
                fit <- glmFit(y, design = design)
                test <- glmLRT(fit)
              }
            )
            mat_out <- t(t(y$counts) / (y$samples$norm.factors * y$samples$lib.size)) * mean(y$samples$lib.size)
            res <- topTags(test, n = Inf) %>%
              as.data.frame() %>%
              tibble::rownames_to_column("gene") %>%
              # flag metrics in results
              mutate(
                de_family = "pseudobulk",
                de_method = de_method,
                de_type = de_type
              )
            return(list(res = res, mat_out = mat_out))
          },
          error = function(e) {
            message(e)
            data.frame()
          }
        )
      },
      DESeq2 = {
        tryCatch(
          {
            if (!is.null(control_col)){
              if (paired_by_replicate) {
                dds <- DESeqDataSetFromMatrix(
                  countData = x,
                  colData = targets,
                  design = ~ replicate + control + group
                )
              } else {
                dds <- DESeqDataSetFromMatrix(
                  countData = x,
                  colData = targets,
                  design = ~ control + group
                )
              }
              
            } else {
              if (paired_by_replicate) {
                dds <- DESeqDataSetFromMatrix(
                  countData = x,
                  colData = targets,
                  design = ~ replicate + group
                )
              } else {
                dds <- DESeqDataSetFromMatrix(
                  countData = x,
                  colData = targets,
                  design = ~ group
                )
              }
            }
            
            dds <- switch(de_type,
              Wald = {
                dds <- try(DESeq(dds,
                  test = "Wald",
                  fitType = "parametric",
                  sfType = "poscounts",
                  betaPrior = F
                ))
              },
              LRT = {
                dds <- try(DESeq(dds,
                  test = "LRT",
                  reduced = ~1,
                  fitType = "parametric",
                  sfType = "poscounts",
                  betaPrior = F
                ))
              }
            )
            res <- results(dds)
            mat_out <- counts(dds, normalized = TRUE)
            # write
            res <- as.data.frame(res) %>%
              mutate(gene = rownames(x)) %>%
              # flag metrics in results
              mutate(
                de_family = "pseudobulk",
                de_method = de_method,
                de_type = de_type
              )
            return(list(res = res, mat_out = mat_out))
          },
          error = function(e) {
            message(e)
            data.frame()
          }
        )
      },
      limma = {
        tryCatch(
          {
            x <- switch(de_type,
              trend = {
                trend_bool <- T
                dge <- DGEList(as.matrix(x), group = targets$group)
                dge <- calcNormFactors(dge)
                x <- new("EList")
                x$E <- cpm(dge, log = TRUE, prior.count = 3)
                x
              },
              voom = {
                counts <- all(as.matrix(x) %% 1 == 0)
                if (counts) {
                  trend_bool <- F
                  x <- voom(as.matrix(x), design)
                  x
                }
              }
            )
            # get fit
            fit <- lmFit(x, design) %>%
              eBayes(trend = trend_bool, robust = trend_bool)
            # format the results
            res <- fit %>%
              # extract all coefs except intercept
              topTable(number = Inf, coef = -1) %>%
              rownames_to_column("gene") %>%
              # flag metrics in results
              mutate(
                de_family = "pseudobulk",
                de_method = de_method,
                de_type = de_type
              )
          },
          error = function(e) {
            message(e)
            data.frame()
          }
        )
      }
    )
  })

  if (de_method %in% c("DESeq2", "edgeR")) {
    res_list <- lapply(results, `[[`, "res")
    mat_list <- lapply(results, `[[`, "mat_out")
    result_df <- res_list %<>% bind_rows(.id = "cell_type")
    return(list(results = result_df, matrices = mat_list))
  } else {
    result_df <- results %<>% bind_rows(.id = "cell_type")
    return(list(results = result_df, matrices = pseudobulks))
  }
}

custom_run_de <- function(input,
                  meta = NULL,
                  replicate_col = 'replicate',
                  cell_type_col = 'cell_type',
                  label_col = 'label',
                  control_col = NULL,
                  sample_col = "Sample.number",
                  paired_by_replicate = FALSE,
                  min_cells = 150,
                  min_reps = 2,
                  min_features = 0,
                  de_family = 'pseudobulk',
                  de_method = 'edgeR',
                  de_type = 'LRT',
                  input_type = 'scRNA',
                  normalization = 'log_tp10k',
                  binarization = FALSE,
                  latent_vars = NULL,
                  n_threads = 2) {
  # first, make sure inputs are correct
  inputs <- check_inputs(
    input = input,
    meta = meta,
    replicate_col = replicate_col,
    cell_type_col = cell_type_col,
    label_col = label_col,
    control_col = control_col,
    sample_col = sample_col
  )

  input <- inputs$expr
  meta <- inputs$meta
  if (is.factor(meta$label)) {
    label_levels <- levels(meta$label)
    cell_type_levels <- levels(meta$cell_type)
    sc <- CreateSeuratObject(input, meta.data = meta) %>% NormalizeData(verbose = F)
    out_stats <- data.frame()
    for (ct in cell_type_levels) {
      label1_barcodes <- meta %>%
        filter(cell_type == ct, label == label_levels[1]) %>%
        rownames(.)
      label2_barcodes <- meta %>%
        filter(cell_type == ct, label == label_levels[2]) %>%
        rownames(.)
      label1_mean_expr <- rowMeans(input[, label1_barcodes])
      label2_mean_expr <- rowMeans(input[, label2_barcodes])
      tmp_stats <- Seurat::FoldChange(sc, label1_barcodes, label2_barcodes, base = exp(1)) %>%
        mutate(gene = rownames(.)) %>%
        magrittr::set_rownames(NULL) %>%
        dplyr::select(gene, avg_logFC, pct.1, pct.2)
      mean_expr <- data.frame(
        gene = names(label1_mean_expr),
        exp1 = label1_mean_expr,
        exp2 = label2_mean_expr
      )
      out_stats %<>% rbind(tmp_stats %>%
        dplyr::left_join(mean_expr, by = "gene") %>%
        mutate(cell_type = ct) %>%
        dplyr::relocate(cell_type, .before = gene))
    }
  }

  # run differential expression
  DE <- switch(de_family,
    pseudobulk = custom_pseudobulk_de(
      input = input,
      meta = meta,
      replicate_col = replicate_col,
      cell_type_col = cell_type_col,
      label_col = label_col,
      control_col = control_col,
      sample_col = sample_col,
      paired_by_replicate = paired_by_replicate,
      min_cells = min_cells,
      min_reps = min_reps,
      min_features = min_features,
      de_family = "pseudobulk",
      de_method = de_method,
      de_type = de_type
    ),
    mixedmodel = mixedmodel_de(
      input = input,
      meta = meta,
      replicate_col = replicate_col,
      cell_type_col = cell_type_col,
      label_col = label_col,
      min_features = min_features,
      de_family = "mixedmodel",
      de_method = de_method,
      de_type = de_type,
      n_threads = n_threads
    ),
    singlecell = singlecell_de(
      input = input,
      meta = meta,
      cell_type_col = cell_type_col,
      label_col = label_col,
      min_features = min_features,
      de_method = de_method,
      normalization = normalization,
      binarization = binarization,
      latent_vars = latent_vars,
      input_type = input_type
    ),
    snapatac_findDAR = singlecell_de(
      input = input,
      meta = meta,
      cell_type_col = cell_type_col,
      label_col = label_col,
      min_features = min_features,
      de_method = de_method
    )
  )

  # clean up the output
  suppressWarnings(
    colnames(DE$results) %<>%
      stringr::str_replace("^logFC\\.group.*", "avg_logFC") %>% # necessary for controlled edgeR output
      forcats::fct_recode(
        "p_val" = "p.value", ## DESeq2
        "p_val" = "pvalue", ## DESeq2
        "p_val" = "p.value", ## t/wilcox
        "p_val" = "P.Value", ## limma
        "p_val" = "PValue", ## edgeR
        "p_val_adj" = "padj", ## DESeq2/t/wilcox
        "p_val_adj" = "adj.P.Val", ## limma
        "p_val_adj" = "FDR", ## edgeR
        "avg_logFC" = "log2FoldChange", ## DESEeq2
        "avg_logFC" = "logFC", ## limma/edgeR
        "avg_logFC" = "avg_log2FC" # Seurat V4
      )
  ) %>%
    as.character()

  DE$results %<>% 
    # calculate adjusted p values
    group_by(cell_type) %>%
    mutate(p_val_adj = p.adjust(p_val, method = "BH")) %>%
    # make sure gene is a character not a factor
    mutate(gene = as.character(gene)) %>%
    dplyr::select(cell_type, gene, avg_logFC, p_val, p_val_adj, de_family, de_method, de_type)

  if (is.factor(meta$label)) {
    DE$results %<>% 
      select(!c(avg_logFC)) %>%
      dplyr::left_join(out_stats, by = c("gene", "cell_type")) %>%
      dplyr::select(
        cell_type,
        gene,
        avg_logFC,
        pct.1,
        pct.2,
        exp1,
        exp2,
        p_val,
        p_val_adj,
        de_family,
        de_method,
        de_type
      ) %>%
      ungroup() %>%
      arrange(cell_type, gene) %>%
      dplyr::rename(
        !!paste0(label_levels[1], ".exp") := exp1,
        !!paste0(label_levels[2], ".exp") := exp2,
        !!paste0(label_levels[1], ".pct") := pct.1,
        !!paste0(label_levels[2], ".pct") := pct.2
      )
  }
  if (input_type == "scATAC") {
    DE$results %<>%
      dplyr::rename(
        da_family = de_family,
        da_method = de_method,
        da_type = de_type
      )
  }
  DE
}
