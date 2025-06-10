#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR4/bin/R
plotHeatmapTopMotifs <- function(overall_prediction_output, format_output, number_genes = 100, ranking_metric = "Range_of_Change") {
    require(dplyr)
    require(pheatmap)

    named_list_predict = lapply(overall_prediction_output, function(i) i$Predicted_Expression)
    mat_predicted_overall <- do.call(rbind, named_list_predict)
    colnames(mat_predicted_overall) <- overall_prediction_output[[1]]$Pseudotime

    top_n_var_genes <- format_output$gene_stats %>%
        filter(Adjusted_P_Value_Pseudotime < 0.05) %>%
        arrange(dplyr::desc(.data[[ranking_metric]])) %>%
        head(number_genes) %>%
        pull(Gene)

    numeric_vals <- as.numeric(colnames(mat_predicted_overall))
    annotation_col <- data.frame(Pseudotime = numeric_vals)
    rownames(annotation_col) <- colnames(mat_predicted_overall)

    # Use a color function for continuous annotation
    ann_colors <- list(
        Pseudotime = colorRampPalette(c("white", "darkblue"))(100)
    )

    pheatmap(
        mat_predicted_overall[top_n_var_genes,],
        cluster_col = FALSE,
        annotation_col = annotation_col,
        annotation_colors = ann_colors,
        show_colnames = FALSE
    )
}

plotChromVARmotifEmbedding <- function(sce, motif_name) {
    require(dplyr)
    require(ggplot2)
    plot_df <- data.frame(colData(sce), t(assay(sce, "counts")))
    plot_df <- plot_df[sample(nrow(plot_df)),]
    plot_df %>%
      ggplot(aes(x = DR1, y = DR2, color = .data[[motif_name]])) +
      geom_point(size = 1) +
      scale_color_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0
      ) +
      theme_minimal()
}

plotSampleTrajectoryMotif <- function(motif, overall_prediction_output, sample_prediction_output) {
  # Load necessary libraries
  require(dplyr)
  require(ggplot2)
  require(RColorBrewer) 

  overall_plot_df <- overall_prediction_output[[motif]]
  overall_plot_df$Sample <- "Global Prediction"
  target_col_order <- colnames(sample_prediction_output[[motif]])
  overall_plot_df <- overall_plot_df[, target_col_order]

  plot_df <- do.call(rbind, list(overall_plot_df, sample_prediction_output[[motif]]))

  # --- Color Palette Setup ---
  unique_samples <- unique(plot_df$Sample)
  sample_specific_names <- setdiff(unique_samples, "Global Prediction")
  n_sample_specific <- length(sample_specific_names)

  # Generate a distinct color palette for sample-specific lines
  if (n_sample_specific > 8) {
    sample_specific_colors <- scales::hue_pal()(n_sample_specific)
  } else {
    sample_specific_colors <- brewer.pal(max(3, n_sample_specific), "Set1")[1:n_sample_specific]
  }

  # Create a named vector for all colors, explicitly assigning "black" to "Global Prediction"
  all_colors <- c("Global Prediction" = "black")
  names(sample_specific_colors) <- sample_specific_names
  all_colors <- c(all_colors, sample_specific_colors)

  # Ensure the levels of the 'Sample' factor are ordered correctly for the legend,
  # placing "Global Prediction" first for clarity.
  plot_df$Sample <- factor(plot_df$Sample, levels = names(all_colors))

  # --- Plot Generation ---
  p <- ggplot(plot_df, aes(x = Pseudotime, y = Predicted_Expression, color = Sample)) +
    # Add the confidence ribbon for the global prediction, slightly transparent grey
    geom_ribbon(data = subset(plot_df, Sample == "Global Prediction"),
                aes(x = Pseudotime, ymin = Lower_CI, ymax = Upper_CI), # <<< FIX: Added x = Pseudotime
                fill = "grey", alpha = 0.45,
                inherit.aes = FALSE) +

    # Add lines for individual sample trajectories
    geom_line(data = subset(plot_df, Sample != "Global Prediction"),
              linewidth = 1) + # Standard line width for samples

    # Add a prominent line for the global prediction (always black and thicker)
    geom_line(data = subset(plot_df, Sample == "Global Prediction"),
              linewidth = 1.8) + # Thicker line for the global trajectory

    # Apply the custom color scale
    scale_color_manual(values = all_colors) +

    # --- Labels and Theme Customization ---
    labs(title = paste("Per-Sample Trajectory for", motif, "with 95% CI"),
         x = "Pseudotime",
         y = "Predicted Expression",
         color = "Sample") + # Legend title
    # Apply a refined minimal theme for a clean, modern look
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold", margin = margin(b = 15)), # Centered, larger, bold title
      axis.title = element_text(size = 14, margin = margin(t = 10, r = 10)), # Larger axis titles
      axis.text = element_text(size = 12), # Larger axis text
      legend.title = element_text(size = 14, face = "bold"), # Larger, bold legend title
      legend.text = element_text(size = 12), # Larger legend text
      legend.position = "right", # Place legend on the right for better use of space
      panel.grid.major = element_line(color = "grey90", linewidth = 0.6), # Finer major grid lines
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.3) # Finer minor grid lines
    )

  return(p)
}
plotMultipleGenesPseudotime <- function(overall_prediction_output, motifs_to_plot, label = NULL) {
  # Load necessary library
  require(dplyr)
  require(ggplot2)
  # Combine data frames
  full_df <- do.call(rbind, overall_prediction_output)
  if (is.null(label)) {
    label = paste(motifs_to_plot, collapse = "/")
  }
  # Generate the plot
  p <- ggplot(subset(full_df, Gene %in% motifs_to_plot),
              aes(x = Pseudotime, y = Predicted_Expression)) +
    geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI, group = Gene), fill = "lightgrey", alpha = 0.2) + # The confidence band
    labs(title = paste("Average Trajectory of",label, "with 95% CI"),
         x = "Pseudotime",
         y = "Predicted Expression",
         color = "Motif") +
    geom_line(linewidth = 1, aes(color = Gene)) + # The predicted trajectory line

    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}