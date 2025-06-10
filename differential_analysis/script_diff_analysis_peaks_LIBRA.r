#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR4/bin/R

print("###################################################################")
print('############### Starting Peak Libra diff analysis #################')
print("###################################################################")
# libraries
suppressMessages(suppressWarnings(library(ArchR)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(argparse)))
suppressMessages(suppressWarnings(library(BSgenome.Hsapiens.UCSC.hg38)))

# parser

parser <- ArgumentParser(description = "Do differential analysis")

parser$add_argument("--archr_path", type = "character", help = "archr project directory path")
parser$add_argument("--out_path", type = "character", default = None, help = "path to output data directory")
parser$add_argument("--code_path", type = "character", default = None, help = "path to code directory")
parser$add_argument("--sample_included", nargs = "+", default = c("All"), help = "list of samples to include (use 'All' to include all)")
parser$add_argument("--diff_annotation", type = "character", default = "State", help = "annotation column to use for differential analysis")
parser$add_argument("--group_forward", type = "character", default = "rec1", help = "forward comparison group name")
parser$add_argument("--group_background", type = "character", default = "ini", help = "background comparison group name")
parser$add_argument("--cell_type_column", type = "character", default = "cell.type", help = "column name for cell type")
parser$add_argument("--replicate_column", type = "character", default = "Patient.number", help = "column name for identifying replicates")
parser$add_argument("--control_column", type = "character", default = "No_control", help = "column name for controlling batch")
parser$add_argument("--paired", type = "logical", default = FALSE, help = "should i pair on replicate")
parser$add_argument("--de_method", type = "character", default = "DESeq2", help = "differential expression method")
parser$add_argument("--de_type", type = "character", default = "Wald", help = "differential expression type")
parser$add_argument("--binarization", type = "logical", default = FALSE, help = "whether to binarize input matrix")
parser$add_argument("--min_cell_threshold", type = "integer", default = 150, help = "minimum number of cells per group/sample")
parser$add_argument("--max_cells", type = "integer", default = 5000, help = "maximum number of cells per group")
parser$add_argument("--reproducibility", type = "character", default = "max(ceiling(n/3),2)", help = "peak reproducibility strategy")
parser$add_argument("--seed", type = "integer", default = 123, help = "random seed")

parser$add_argument("--motifdb", type = "character", default = "JASPAR2020", help = "Motif DB used in motif enrichment")
parser$add_argument("--threshold_FDR", type = "integer", default = 0.05, help = "Threshold of adjusted pvalue for peaks")
parser$add_argument("--threshold_Log2FC", type = "integer", default = 0.5, help = "Threshold of Log2FC for peaks")
parser$add_argument("--peak_call", type = "logical", default = FALSE, help = "whether I should do the peak calling")
parser$add_argument("--use_peak", type = "character", default = NULL, help = "path to peakset i should use if i don't do peak calling")

# Parse arguments
args <- parser$parse_args()

proj_global_dir <- args$archr

OUT_LIBRA <- args$out_path

sample_included <- args$sample_included
diff_annotation <- args$diff_annotation
group_forward <- args$group_forward
group_background <- args$group_background
cell_type_column <- args$cell_type_column
replicate_column <- args$replicate_column
control_column <-  args$control_column
paired <- args$paired
de_method <- args$de_method
de_type <- args$de_type
binarization <- args$binarization
minCellThreshold <- args$min_cell_threshold
maxCells <- args$max_cells
reproducibility_strategy <- args$reproducibility
SEED <- args$seed
motifdb <- args$motifdb
threshold_FDR <- args$threshold_FDR
threshold_LogFC <- args$threshold_Log2FC
peak_call <- args$peak_call
use_peak <- args$use_peak

if (motifdb %in% c("JASPAR2020", "JASPAR2018")) {
  library(motifdb, character.only=TRUE)
}

# paths
OUT_ARROWS <- file.path(DATA_PATH, "Arrows", project_name)
OUT_LOGS <- file.path(DATA_PATH, "outputs", project_name, "ArchRLogs")
OUT_ARCHR <- file.path(DATA_PATH, "outputs", project_name, "ArchROutput")
OUT_LOGS_LARGE <- file.path(DATA_PATH, "outputs", project_name, "ArchRLogs_Large")
OUT_ARCHR_LARGE <- file.path(DATA_PATH, "outputs", project_name, "ArchRProject_Large")
OUT_MSG <- file.path(DATA_PATH, "outputs", project_name, "Prints")
OUT_PEAKS <- file.path(DATA_PATH, "outputs", project_name, "Peaks")
OUT_LIBRA <- file.path(DATA_PATH, "outputs", project_name, "Libra")

set.seed(SEED)
datetime <- strftime(Sys.time(), format = "%Y:%m:%d_%H:%M")

# functions scripts
source_all_r_files <- function(path_list) {
  origin_wd <- getwd()
  lapply(
    path_list, function(x) {
      setwd(x)
      file.sources <- list.files(pattern = "*.r$")
      sapply(file.sources, source, .GlobalEnv)
    }
  )
  setwd(origin_wd)
}
source_all_r_files("PATH TO CHANGE") ################################### TO CHANGE ################### change this path to the folder containing the used functions
check_a_path(OUT_LIBRA)
# set up
addArchRThreads(ceiling(parallelly::availableCores() / 2))
addArchRGenome("hg38")
addArchRLocking(locking = T)

# load project
proj_global <- loadArchRProject(proj_global_dir)

# run
if (peak_call) {
  proj_global <- addCombinedAnnotation(proj_global, diff_annotation, cell_type_column)
  combined_name <- paste(diff_annotation, "by", cell_type_column, sep = "_")

  count_cells <- proj_global@cellColData[, c("Sample.number", combined_name)] %>%
    as.data.frame() %>%
    group_by(.data[[combined_name]], .data[["Sample.number"]]) %>%
    summarize(count = n(), .groups = "drop")
  
  count_per_condition <- count_cells %>% filter(count >= minCellThreshold) %>%
    group_by(.data[[combined_name]]) %>%
    summarize(num_samples = n(), .groups = "drop")

  message("Number passing conditions :")
  message(count_per_condition)
  

  minReplicate <- count_per_condition %>%
    filter(
      !grepl("undetermined", .data[[combined_name]]) &
      (grepl(group_forward, .data[[combined_name]]) | grepl(group_background, .data[[combined_name]]))
    ) %>%
    summarize(min(num_samples)) %>% as.numeric()

maxReplicate <- count_per_condition %>%
  filter(
    !grepl("undetermined", .data[[combined_name]]) &
    (grepl(group_forward, .data[[combined_name]]) | grepl(group_background, .data[[combined_name]]))
  ) %>%
  summarize(max(num_samples)) %>% as.numeric() + 2

  proj_global <- addGroupCoverages(
    ArchRProj = proj_global,
    groupBy = combined_name, 
    useLabels = TRUE,
    sampleLabels = "Sample",
    minCells = minCellThreshold,
    maxCells = maxCells,
    minReplicates = minReplicate,
    maxReplicates = maxReplicate,
    force = TRUE
  )

  proj_global <- addReproduciblePeakSet(
    ArchRProj = proj_global,
    groupBy = combined_name,
    reproducibility = reproducibility_strategy,
    force = TRUE
  )

  proj_global <- addPeakMatrix(proj_global)
  proj_global <- saveArchRProject(proj_global)
} else if (!is.null(use_peak)) {
  proj_global@peakSet <- readRDS(use_peak)
  proj_global <- addPeakMatrix(proj_global)
  proj_global <- saveArchRProject(proj_global)
}

if (!identical(sample_included, c("All"))) {
  idx_condition <- BiocGenerics::which(proj_global@cellColData[,"Sample.number"] %in% sample_included)
  cellsSample <- proj_global$cellNames[idx_condition]
  proj_global <- proj_global[cellsSample, ]
}


if(control_column == "No_control") {
  control_col <- NULL
} else {
  control_col <- control_column
}

if (is.numeric(proj_global@cellColData[, diff_annotation])) {
  group_forward <- NULL
  group_background <- NULL
}

out_libra_pseudobulk <- getLibraPseudoBulkAnalysis(proj_global, diff_annotation, group_forward, group_background,
                   cell_type = cell_type_column,
                   replicate_col = replicate_column,
                   control_col = control_col,
                   paired_by_replicate = paired,
                   useMatrix = "PeakMatrix",
                   de_family = "pseudobulk",
                   de_method = de_method,
                   de_type = de_type,
                   minCell= minCellThreshold,
                   n_threads = ceiling(parallelly::availableCores() / 2)
                   )


peak_global <- getPeakSet(proj_global)
results_se <- convertLibraOutToSE(out_libra_pseudobulk$results, peak_global)

# motif analysis
proj_global <- addMotifAnnotations(proj_global, motifSet = motifdb, annoName = "Motif", force = TRUE)


threshold_up = paste('FDR <=', threshold_FDR, '& Log2FC >=', threshold_LogFC, sep = " ")

motifsUp <- peakAnnoEnrichment(
  seMarker = results_se,
  ArchRProj = proj_global,
  peakAnnotation = "Motif",
  cutOff = threshold_up
  )

threshold_down = paste('FDR <= ', threshold_FDR, ' & Log2FC <= -', threshold_LogFC, sep = "")

motifsDown <- peakAnnoEnrichment(
  seMarker = results_se,
  ArchRProj = proj_global,
  peakAnnotation = "Motif",
  cutOff = threshold_down
  )


#combining results into one rds
results_fin <- list(Peaks =  results_se, MotifsUp = motifsUp, MotifsDown = motifsDown, matrices = out_libra_pseudobulk$matrices, Libra = out_libra_pseudobulk$results)

saveRDS(results_fin, file = file.path(OUT_LIBRA, paste(datetime, basename(proj_global_dir), "libra_diff", diff_annotation, control_column, de_method, de_type, sep = "_")))

print("###################################################################")
print('############# Successful analysis for peaks analysis ##############')
print("###################################################################")

