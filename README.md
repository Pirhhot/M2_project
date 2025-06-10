# Heterogeneity Analysis of IDH-Mutant Gliomas

This repository contains the code for my master's project investigating the heterogeneity of IDH-mutant gliomas using single-cell ATAC-seq data.
It can both be used in practice to perform the analysis or, looking at the provided scripts and functions, better understand how they are performed in my report.

## Repository Overview

The analysis is structured around two complementary approaches presented in the thesis:

### Paired Differential Analysis
This analysis examines differential chromatin accessibility patterns across conditions. The folder contains:
- A main processing script that handles the differential accessibility analysis workflow
- A Jupyter notebook with plotting functions and visualization code
- An HTML rendering of the notebook demonstrating figure generation as presented in the report
- A `functions/` subdirectory with utility functions specific to this analysis

### Trajectory Analysis
This analysis employs chromVAR values with Generalized Additive Mixed Models (GAMMs) to model temporal trajectories. The folder includes:
- A main processing script for trajectory modeling and analysis
- A Jupyter notebook containing plotting functions and visualization methods
- An HTML rendering showing the figure generation process
- A `functions/` subdirectory with utility functions specific to this analysis
### Additional Components
- **`environment.yml`**: Complete conda environment specifications

## Data Requirements

- **Differential Analysis**: Requires a preprocessed ArchR project object
- **Trajectory Analysis**: Requires a SingleCellExperiment object with chromVAR values as an assay

## Installation

The repository includes comprehensive conda environment specifications in `environment.yml`. This file is optimized for CentOS7 Linux systems and contains extensive package dependencies. It should serve as a reference for environment setup rather than direct installation.

### Core Dependencies
The primary packages can be installed using:
```bash
conda install -c conda-forge bioconda r-archr
conda install -c r r-mgcv  
conda install -c bioconda pheatmap
```

Additional dependencies are documented in the YML specification file.

## Usage
1. Modify the scripts to your file storage as explained in the scripts' comments
2. Execute the main processing scripts in each analysis directory to generate intermediary results
3. Use the Jupyter notebooks to reproduce the figures and visualizations shown in the HTML outputs
4. Utility functions in each `functions/` subdirectory provide analysis-specific support for the corresponding scripts and notebooks
   
## Output

The analyses generate the figures and results presented in the master's thesis, demonstrating heterogeneity patterns in IDH-mutant gliomas through differential accessibility and trajectory modeling approaches.
