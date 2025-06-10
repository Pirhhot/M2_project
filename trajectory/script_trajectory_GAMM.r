#! !!!!!!!!!!!!!!!change to your r bin !!!!!!!!!!!

#libraries
library(ArchR)
library(dplyr)
library(tidyr)
library(argparse)

SEED <- 123
set.seed(SEED)

#functions scripts
source_all_r_files <- function(path_list){
    origin_wd = getwd()
    lapply(
        path_list, function(x) {
            setwd(x)
            file.sources = list.files(pattern =  "*.r$")
            sapply(file.sources,source,.GlobalEnv)
        }
    )
    setwd(origin_wd)
}
source_all_r_files("PATH TO CHANGE") ################################### TO CHANGE ################### change this path to the folder containing the used functions
parser <- ArgumentParser(description = "Do trajectory analysis")

parser$add_argument("--sce_path", type = "character", help = "sce input path")
parser$add_argument("--out_path", type = "character", default = None, help = "path to output directory")
parser$add_argument("--name_assay", type = "character", default = None, help = "name of stored assay to use")
parser$add_argument("--sample_col", type = "character", default = None, help = "column used as sample")
parser$add_argument("--pseudotime_col", type = "character", default = None, help = "column used as pseudotime")

#load project
sce_path <- args$sce_path

OUT_DIR <- args$out_path

name_assay <- args$name_assay
sample_col <- args$sample_col
pseudotime_col <- args$pseudotime_col

my_sce <- readRDS(sce_path)

print("##### Starting Fitting model #######")

fit_output <- fit_pseudotime_gam(my_sce,
                               expr_assay = name_assay,
                               pseudotime_col = pseudotime_col,
                               sample_col = sample_col,
                               n_cores = 4,
                               verbose = TRUE)

gc()
print("##### Fitting succeded ! #######")

print("Starting exctracting")

format_output <- extractFormatGAMoutput(fit_output) 
print("Finished exctracting")

saveRDS(format_output, file.path(OUT_DIR, "formatted_output.rds"))
print("Saved!")
rm(format_output)
gc()

print("Starting overall")

overall_predictions <- generate_overall_pseudotime_predictions(fit_output)
print("Finished overall")

saveRDS(overall_predictions, file.path(OUT_DIR, "overall_predictions_output.rds"))
print("Saved!")

rm(overall_predictions)
gc()

print("Starting sample")

sample_predictions <- generate_per_sample_predictions(fit_output)
print("Finished sample")

saveRDS(sample_predictions, file.path(OUT_DIR, "sample_prediction_output.rds"))
print("Saved!")