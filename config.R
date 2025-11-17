# Configuration File for scRNA-seq Analysis
# This file contains all configurable paths and parameters
# Edit this file to match your local setup

# ============================================================================
# DATA PATHS - UPDATE THESE BEFORE RUNNING
# ============================================================================

# Directory containing 10x Genomics files (barcodes.tsv, genes.tsv, matrix.mtx)
DATA_DIR <- "PATH-TO-FILES/"
# Example Windows: "C:/Users/YourName/Documents/scRNA-seq/data/"
# Example Mac/Linux: "/Users/YourName/scRNA-seq/data/"

# Path to reference dataset for cell type annotation
REFERENCE_PATH <- "PATH-TO-REFERENCE/reference.rds"
# Example: "data/reference.rds"

# ============================================================================
# OUTPUT PATHS
# ============================================================================

# Directory for saving results
RESULTS_DIR <- "results/"

# Directory for figures
FIGURE_DIR <- "results/figures/"

# Directory for tables
TABLE_DIR <- "results/tables/"

# ============================================================================
# ANALYSIS PARAMETERS
# ============================================================================

# Project name
PROJECT_NAME <- "pbmc3k"

# Quality Control Thresholds
QC_PARAMS <- list(
  min_counts = 1000,        # Minimum UMI counts per cell
  max_counts = 4000,        # Maximum UMI counts per cell
  min_features = 500,       # Minimum genes per cell
  max_features = 1200,      # Maximum genes per cell
  max_mito = 5,            # Maximum mitochondrial percentage
  min_ribo = 20            # Minimum ribosomal percentage
)

# Normalization
NORM_PARAMS <- list(
  method = "LogNormalize",
  scale_factor = 10000
)

# Variable Features
HVG_PARAMS <- list(
  method = "vst",
  n_features = 2000
)

# PCA
PCA_PARAMS <- list(
  n_pcs = 50,              # Total PCs to compute
  pcs_use = 10             # PCs to use for downstream analysis
)

# Clustering
CLUSTER_PARAMS <- list(
  resolution_test = seq(0.1, 2.0, 0.1),  # Range of resolutions to test
  resolution_final = 0.5                  # Final resolution to use
)

# UMAP/tSNE
REDUCTION_PARAMS <- list(
  seed = 1001,
  pcs_use = 10
)

# Differential Expression
DEG_PARAMS <- list(
  only_pos = FALSE,
  min_pct = 0.25,
  logfc_threshold = 0.5,
  p_adj_cutoff = 0.05
)

# ============================================================================
# OPTIONAL ANALYSIS PARAMETERS
# ============================================================================

# SCTransform parameters
SCT_PARAMS <- list(
  n_features = 3000,
  pcs_use = 25,
  resolution = 0.7
)

# Gene Set Enrichment
ENRICHMENT_PARAMS <- list(
  ontology = c("BP", "CC", "MF"),  # GO categories
  show_category = 10,
  p_adj_cutoff = 0.05
)

# ============================================================================
# MARKER GENES FOR MANUAL ANNOTATION
# ============================================================================

MARKER_GENES <- list(
  B_cells = c("CD19", "MS4A1", "CD79A", "CD79B"),
  T_cells = c("CD3D", "CD3E", "CD3G"),
  CD4_T_cells = c("CD4", "IL7R"),
  CD8_T_cells = c("CD8A", "CD8B"),
  NK_cells = c("GNLY", "NKG7", "NCAM1"),
  Monocytes = c("CD14", "LYZ", "FCGR3A"),
  Dendritic_cells = c("FCER1A", "CST3")
)

# ============================================================================
# VISUALIZATION PARAMETERS
# ============================================================================

VIZ_PARAMS <- list(
  point_size = 1.5,
  label_size = 5,
  colors_divergent = c("lightgrey", "red"),
  figure_width = 10,
  figure_height = 6
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Create output directories
create_output_dirs <- function() {
  dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(FIGURE_DIR, "qc_plots"), showWarnings = FALSE)
  dir.create(file.path(FIGURE_DIR, "clustering"), showWarnings = FALSE)
  dir.create(file.path(FIGURE_DIR, "cell_types"), showWarnings = FALSE)
  dir.create(file.path(FIGURE_DIR, "degs"), showWarnings = FALSE)
  message("Output directories created successfully!")
}

# Validate data paths
validate_paths <- function() {
  errors <- c()

  if (!dir.exists(DATA_DIR)) {
    errors <- c(errors, paste("Data directory not found:", DATA_DIR))
  }

  required_files <- c("barcodes.tsv", "genes.tsv", "matrix.mtx")
  for (file in required_files) {
    full_path <- file.path(DATA_DIR, file)
    if (!file.exists(full_path)) {
      errors <- c(errors, paste("Required file not found:", full_path))
    }
  }

  if (length(errors) > 0) {
    stop("Path validation failed:\n", paste(errors, collapse = "\n"))
  } else {
    message("All required data files found!")
  }
}

# Print configuration summary
print_config <- function() {
  cat("\n")
  cat("========================================\n")
  cat("scRNA-seq Analysis Configuration\n")
  cat("========================================\n")
  cat("Project:", PROJECT_NAME, "\n")
  cat("Data directory:", DATA_DIR, "\n")
  cat("Results directory:", RESULTS_DIR, "\n")
  cat("\n")
  cat("Quality Control:\n")
  cat("  - Count range:", QC_PARAMS$min_counts, "-", QC_PARAMS$max_counts, "\n")
  cat("  - Feature range:", QC_PARAMS$min_features, "-", QC_PARAMS$max_features, "\n")
  cat("  - Max mito %:", QC_PARAMS$max_mito, "\n")
  cat("\n")
  cat("Analysis Parameters:\n")
  cat("  - HVGs:", HVG_PARAMS$n_features, "\n")
  cat("  - PCs to use:", PCA_PARAMS$pcs_use, "\n")
  cat("  - Clustering resolution:", CLUSTER_PARAMS$resolution_final, "\n")
  cat("========================================\n")
  cat("\n")
}
