#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR4/bin/R
library(SingleCellExperiment)
library(mgcv)
library(foreach)
library(doSNOW) # For progress bar with foreach + parallel

fit_gene <- function(gene_expr, pseudotime, sample, gene_name = NULL, verbose = FALSE) {
    sample_factor <- as.factor(sample)
    df <- data.frame(
        Expression = as.numeric(gene_expr),
        Pseudotime = as.numeric(pseudotime),
        Sample.number = sample_factor,
        stringsAsFactors = FALSE
    )

    model_fit <- tryCatch({
        mgcv::bam(Expression ~ s(Pseudotime) +
                              s(Sample.number, bs = "re") +
                              s(Pseudotime, Sample.number, bs = "fs"),
                  data = df,
                  family = gaussian(link = "identity"),
                  discrete = TRUE)
    }, error = function(e) {
        warning(paste0("Full interaction model failed for gene ", gene_name, ": ", e$message, ". Falling back to no-interaction model."), call. = FALSE)
        mgcv::bam(Expression ~ s(Pseudotime) + s(Sample.number, bs = "re"),
                  data = df,
                  family = gaussian(link = "identity"),
                  discrete = TRUE)
    })

    # --- AGGRESSIVE STRIPPING OF HEAVY COMPONENTS ---
    if (!is.null(model_fit) && inherits(model_fit, "gam")) {
        model_fit$model <- NULL      # Remove the model frame (redundant if fit=FALSE, but good practice)
        model_fit$data <- NULL       # Remove data slot (if any remains)
        model_fit$X <- NULL          # Remove the overall design matrix
        model_fit$y <- NULL          # Remove response vector
        model_fit$residuals <- NULL  # Remove residuals
        model_fit$fitted.values <- NULL # Fitted values can be re-predicted
        model_fit$weights <- NULL    # Remove observation weights
        model_fit$prior.weights <- NULL # Remove prior weights
        model_fit$call <- NULL       # <--- CRUCIAL: Breaks links to potentially large calling environment

        # Also strip design matrices from individual smooth terms
        if (!is.null(model_fit$smooth)) {
            for (j in seq_along(model_fit$smooth)) {
                if (!is.null(model_fit$smooth[[j]]$X)) {
                    model_fit$smooth[[j]]$X <- NULL # Remove design matrix bits for smooths
                }
                # Remove 'by.y' if present, sometimes it holds data references
                if (!is.null(model_fit$smooth[[j]]$by.y)) {
                    model_fit$smooth[[j]]$by.y <- NULL
                }
            }
        }
    }

    return(model_fit)
}


fit_pseudotime_gam <- function(sce,
                               expr_assay = "counts",
                               pseudotime_col = "pseudotime",
                               sample_col = "Sample.number",
                               n_cores = 4,
                               verbose = TRUE) {

    # Check input
    if (!pseudotime_col %in% colnames(colData(sce))) {
        stop(paste0("'", pseudotime_col, "' not found in colData(sce)"))
    }
    if (!sample_col %in% colnames(colData(sce))) {
        stop(paste0("'", sample_col, "' not found in colData(sce)"))
    }
    if (!expr_assay %in% assayNames(sce)) {
        stop(paste0("'", expr_assay, "' assay not found in sce"))
    }

    # Extract raw data
    pseudotime_raw <- colData(sce)[[pseudotime_col]]
    sample_raw <- colData(sce)[[sample_col]]
    expr_matrix_raw <- assay(sce, expr_assay)
    
    na_pseudotime <- is.na(pseudotime_raw)
    na_sample <- is.na(sample_raw) # Check for NAs in sample_raw here as well.
    if (sum(na_sample) > 0) {
        stop("Missing values for sample.")
    }
    # Determine which cells (columns) to keep
    cells_to_keep <- !(na_pseudotime) # Filter if either pseudotime OR sample is NA

    if (verbose) {
        n_removed_cells <- sum(!cells_to_keep)
        if (n_removed_cells > 0) {
            cat("Removing", n_removed_cells, "cells due to missing pseudotime.\n")
        } else {
            cat("No cells removed due to missing pseudotime or sample information.\n")
        }
    }

    # Filter the data
    pseudotime_filtered <- pseudotime_raw[cells_to_keep]
    sample_filtered <- as.factor(sample_raw[cells_to_keep]) 
    expr_matrix_filtered <- expr_matrix_raw[, cells_to_keep, drop = FALSE]
    gene_names_filtered <- rownames(expr_matrix_filtered)

    if (verbose) {
        n_genes_removed <- nrow(expr_matrix_raw) - nrow(expr_matrix_filtered)
        if (n_genes_removed > 0) {
            cat("Removed", n_genes_removed, "genes that had no non-NA/non-zero expression after cell filtering.\n")
        }
        cat("Proceeding with", ncol(expr_matrix_filtered), "cells and", nrow(expr_matrix_filtered), "genes for modeling.\n")
    }

    # Check if there are any genes left to process
    if (nrow(expr_matrix_filtered) == 0) {
        warning("No genes remaining after filtering. Returning empty list.")
        return(list())
    }

    # --- Set up parallel processing with doSNOW and foreach ---
    if (verbose) cat("Fitting GAMs for", nrow(expr_matrix_filtered), "genes using foreach + doSNOW...\n")

    cl <- makeCluster(n_cores)
    registerDoSNOW(cl) # Register the SNOW backend for foreach

    # Setup progress bar
    iterations <- nrow(expr_matrix_filtered)
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    # Use foreach with %dopar% for parallel processing and progress bar
    results <- foreach(i = seq_len(iterations),
                       .options.snow = opts,
                       .export = c("fit_gene", "pseudotime_filtered", "sample_filtered"), 
                       .packages = c("mgcv")) %dopar% {
        fit_gene(expr_matrix_filtered[i, ], pseudotime_filtered, sample_filtered, gene_names_filtered[i], verbose = FALSE)
    }

    close(pb) # Close the progress bar
    stopCluster(cl) # Stop the cluster
    print("Fitting complete !")
    names(results) <- gene_names_filtered
    
    df <- data.frame(t(expr_matrix_filtered))
    df$Pseudotime <- pseudotime_filtered
    df$Sample <- sample_filtered
    return(list(results = results, df = df, pseudotime_filtered = pseudotime_filtered))
}
