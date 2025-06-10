#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR4/bin/R
library(SingleCellExperiment)
library(mgcv)
library(foreach)
library(doSNOW) # For progress bar with foreach + parallel

extractFormatGAMoutput <- function(fit_gam_output_list, p.adjust.method = "bonferroni") {
    gam_results_list <- fit_gam_output_list$results
    pseudotime_filtered <- fit_gam_output_list$pseudotime_filtered
    sample_filtered <- fit_gam_output_list$df$Sample
    gene_names <- names(gam_results_list)
    
    p_values_pseudotime <- numeric(length(gene_names))
    max_fitted_original <- numeric(length(gene_names))
    min_fitted_original <- numeric(length(gene_names))
    range_of_change <- numeric(length(gene_names))
    variance_along_pseudotime <- numeric(length(gene_names))
    
    # List to temporarily store ordered fitted vectors before forming matrix
    all_fitted_vectors_ordered <- vector("list", length(gene_names))
    order_idx <- order(pseudotime_filtered)
    ordered_pseudotime_points <- pseudotime_filtered[order_idx]
    num_pseudotime_points <- length(ordered_pseudotime_points)

    original_data_for_predict <- data.frame(
        Pseudotime = pseudotime_filtered,
        Sample.number = sample_filtered,
        stringsAsFactors = FALSE
    )

    for (i in seq_along(gam_results_list)) {
        gene_name <- gene_names[i]
        model_object <- gam_results_list[[i]]
        
        # Default values for this gene if model is invalid or extraction fails
        p_values_pseudotime[i] <- NA
        max_fitted_original[i] <- NA
        min_fitted_original[i] <- NA
        range_of_change[i] <- NA
        variance_along_pseudotime[i] <- NA
        all_fitted_vectors_ordered[[i]] <- rep(NA, num_pseudotime_points) 

        if (is.null(model_object) || !inherits(model_object, "gam")) {
            warning(paste("Model for gene", gene_name, "is not a valid 'gam' object or is NULL. Skipping full extraction for this gene."))
            next
        }
        
        model_summary <- summary(model_object)
        
        if ("s(Pseudotime)" %in% rownames(model_summary$s.table)) {
            p_values_pseudotime[i] <- model_summary$s.table["s(Pseudotime)", "p-value"]
        } else {
            warning(paste("s(Pseudotime) term not found in summary for gene", gene_name, ". Assigning NA p-value."))
        }
        
        fitted_values_current_gene <- tryCatch({
            predict(model_object, newdata = original_data_for_predict, type = "response")
        }, error = function(e) {
            warning(paste("Could not predict fitted values for gene", gene_name, ":", e$message))
            return(NA)
        })
        
        if (!any(is.na(fitted_values_current_gene)) && length(fitted_values_current_gene) == num_pseudotime_points) {
            # Order the fitted values according to the sorted pseudotime indices
            fitted_values_current_gene_ordered <- fitted_values_current_gene[order_idx]
            all_fitted_vectors_ordered[[i]] <- fitted_values_current_gene_ordered

            current_max_expr <- max(fitted_values_current_gene_ordered)
            current_min_expr <- min(fitted_values_current_gene_ordered)
            
            max_fitted_original[i] <- current_max_expr
            min_fitted_original[i] <- current_min_expr
            range_of_change[i] <- current_max_expr - current_min_expr
            variance_along_pseudotime[i] <- var(fitted_values_current_gene_ordered)
        } else {
            warning(paste("Fitted values for gene", gene_name, "are missing or length mismatch. Assigning NAs for strength metrics."))
        }
    }
    
    adjusted_p_values <- p_values_pseudotime
    valid_p_values_idx <- !is.na(p_values_pseudotime)
    
    if (sum(valid_p_values_idx) > 0) {
        adjusted_p_values[valid_p_values_idx] <- p.adjust(p_values_pseudotime[valid_p_values_idx], method = p.adjust.method)
    }
    
    gene_stats_df <- data.frame(
        Gene = gene_names,
        P_Value_Pseudotime = p_values_pseudotime,
        Adjusted_P_Value_Pseudotime = adjusted_p_values,
        Max_Fitted = max_fitted_original,
        Min_Fitted = min_fitted_original,
        Range_of_Change = range_of_change,
        Variance_Predicted_Values = variance_along_pseudotime,
        stringsAsFactors = FALSE
    )
    
    # Create the fitted values matrix
    fitted_values_matrix <- do.call(rbind, all_fitted_vectors_ordered)
    rownames(fitted_values_matrix) <- gene_names
    colnames(fitted_values_matrix) <- paste0("Cell_PT_", round(ordered_pseudotime_points, 4))
    
    return(list(
        gene_stats = gene_stats_df,
        fitted_values_matrix = fitted_values_matrix,
        pseudotime_points = ordered_pseudotime_points
    ))
}

generate_overall_pseudotime_predictions <- function(fit_gam_output_list, num_pseudotime_points = 100, ci_level = 0.95) {
    gam_results_list <- fit_gam_output_list$results
    pseudotime_filtered <- fit_gam_output_list$pseudotime_filtered
    sample_filtered <- as.factor(fit_gam_output_list$df$Sample)
    min_pt <- min(pseudotime_filtered, na.rm = TRUE)
    max_pt <- max(pseudotime_filtered, na.rm = TRUE)
    all_predictions_list = list()
    # Create a prediction grid that only contains Pseudotime, as sample terms are excluded
    prediction_grid <- data.frame(
        Pseudotime = seq(min_pt, max_pt, length.out = num_pseudotime_points),
        Sample.number = factor(levels(sample_filtered)[1], levels = levels(sample_filtered)) # Add placeholder
    )
    gene_names <- names(gam_results_list)

    cat("Generating overall (average) pseudotime predictions with CIs...\n")
    pb_pred <- txtProgressBar(max = length(gam_results_list), style = 3)

    z_score <- qnorm(1 - (1 - ci_level) / 2)

    for (i in seq_along(gam_results_list)) {
        gene_name <- gene_names[i]
        model_object <- gam_results_list[[i]]

        if (is.null(model_object) || !inherits(model_object, "gam")) {
            warning(paste("Model for gene", gene_name, "is not valid. Skipping predictions."))
            all_predictions_list[[gene_name]] <- NULL
            next
        }

        # Predict, EXCLUDING the sample-specific random effect and interaction terms
        # Ensure the term names match those in summary(model_object)$smooth.labels
        predictions_se <- tryCatch({
            predict(model_object,
                    newdata = prediction_grid,
                    type = "response", # Predict on the response scale
                    se.fit = TRUE,
                    # Exclude the specific smooth terms related to Sample.number
                    # The names should match summary(model_object)$smooth.labels
                    exclude = c("s(Sample.number)", "s(Pseudotime,Sample.number)"))
        }, error = function(e) {
            warning(paste("Overall prediction failed for gene", gene_name, ":", e$message))
            return(NULL)
        })

        if (!is.null(predictions_se)) {
            predicted_expression <- predictions_se$fit
            se_predicted_expression <- predictions_se$se.fit
            
            # Calculate confidence interval bounds
            lower_ci <- predicted_expression - z_score * se_predicted_expression
            upper_ci <- predicted_expression + z_score * se_predicted_expression

            gene_predictions_df <- data.frame(
                Gene = gene_name,
                Pseudotime = prediction_grid$Pseudotime,
                Predicted_Expression = predicted_expression,
                SE_Predicted_Expression = se_predicted_expression,
                Lower_CI = lower_ci,
                Upper_CI = upper_ci,
                stringsAsFactors = FALSE
            )
            all_predictions_list[[gene_name]] <- gene_predictions_df
        }
        setTxtProgressBar(pb_pred, i)
    }
    close(pb_pred)

    return(all_predictions_list)
}


generate_per_sample_predictions <- function(fit_gam_output_list, num_pseudotime_points = 100, ci_level = 0.95) {
    gam_results_list <- fit_gam_output_list$results
    pseudotime_filtered <- fit_gam_output_list$pseudotime_filtered
    sample_filtered <- fit_gam_output_list$df$Sample

    unique_samples <- levels(sample_filtered)
    min_pt <- min(pseudotime_filtered, na.rm = TRUE)
    max_pt <- max(pseudotime_filtered, na.rm = TRUE)
    
    prediction_grid_df_list <- lapply(unique_samples, function(s_num) {
        data.frame(
            Pseudotime = seq(min_pt, max_pt, length.out = num_pseudotime_points),
            Sample.number = factor(s_num, levels = unique_samples)
        )
    })
    prediction_grid <- do.call(rbind, prediction_grid_df_list)

    all_predictions_list <- list()
    gene_names <- names(gam_results_list)

    cat("Generating per-sample predictions with CIs...\n")
    pb_pred <- txtProgressBar(max = length(gam_results_list), style = 3)

    z_score <- qnorm(1 - (1 - ci_level) / 2)

    for (i in seq_along(gam_results_list)) {
        gene_name <- gene_names[i]
        model_object <- gam_results_list[[i]]

        if (is.null(model_object) || !inherits(model_object, "gam")) {
            warning(paste("Model for gene", gene_name, "is not valid. Skipping predictions."))
            all_predictions_list[[gene_name]] <- NULL
            next
        }

        predictions_se <- tryCatch({
            predict(model_object, newdata = prediction_grid, se.fit = TRUE)
        }, error = function(e) {
            warning(paste("Prediction failed for gene", gene_name, ":", e$message))
            return(NULL)
        })

        if (!is.null(predictions_se)) {
            predicted_expression <- predictions_se$fit
            se_predicted_expression <- predictions_se$se.fit
            
            lower_ci <- predicted_expression - z_score * se_predicted_expression
            upper_ci <- predicted_expression + z_score * se_predicted_expression

            gene_predictions_df <- data.frame(
                Gene = gene_name,
                Pseudotime = prediction_grid$Pseudotime,
                Sample = prediction_grid$Sample.number,
                Predicted_Expression = predicted_expression,
                SE_Predicted_Expression = se_predicted_expression,
                Lower_CI = lower_ci,
                Upper_CI = upper_ci,
                stringsAsFactors = FALSE
            )
            all_predictions_list[[gene_name]] <- gene_predictions_df
        }
        setTxtProgressBar(pb_pred, i)
    }
    close(pb_pred)

    return(all_predictions_list)
}