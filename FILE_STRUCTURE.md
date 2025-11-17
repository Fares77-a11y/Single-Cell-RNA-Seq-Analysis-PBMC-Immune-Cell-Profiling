# Repository File Structure and Summary

## Generated Files Overview

This document provides a complete overview of all files generated for the scRNA-seq analysis repository.

---

## 📁 Core Analysis Files

### 1. `scRNA_seq_analysis.Rmd` 
**Type:** RMarkdown Notebook  
**Purpose:** Interactive, reproducible analysis with embedded documentation  
**Key Features:**
- Complete workflow from data loading to differential expression
- Embedded plots and results
- Generates HTML report when knitted
- Best for: Learning, exploration, sharing results

**How to Use:**
```r
# In RStudio
rmarkdown::render("scRNA_seq_analysis.Rmd")
```

---

### 2. `scRNA_seq_analysis_main.R`
**Type:** R Script  
**Purpose:** Core analysis pipeline as standalone script  
**Key Features:**
- All essential analysis steps
- Well-commented and annotated
- Modular structure with clear sections
- Best for: Batch processing, automation

**Sections:**
1. Setup and Installation
2. Load Libraries
3. Load Data
4. Quality Control
5. Normalization
6. Variable Feature Selection
7. Scaling
8. PCA
9. Clustering
10. t-SNE/UMAP
11. Cell Type Prediction (Manual)
12. Cell Type Prediction (Reference-based)
13. Differential Expression
14. Hierarchical Clustering

**How to Use:**
```r
source("scRNA_seq_analysis_main.R")
```

---

### 3. `scRNA_seq_analysis_optional.R`
**Type:** R Script  
**Purpose:** Advanced and optional analyses  
**Key Features:**
- Gene set overrepresentation (GO, KEGG)
- SCTransform normalization workflow
- Additional visualization options

**Sections:**
1. Gene Set Overrepresentation Analysis
   - GO: Biological Process, Cellular Component, Molecular Function
   - KEGG Pathways
2. SCTransform Workflow
   - Alternative normalization approach
   - Enhanced clustering with more PCs
   - Marker validation

**How to Use:**
```r
# After running main analysis
source("scRNA_seq_analysis_optional.R")
```

---

## 📋 Documentation Files

### 4. `README.md`
**Type:** Markdown  
**Purpose:** Main repository documentation  
**Contents:**
- Project overview
- Dataset description
- Installation instructions (brief)
- Usage guidelines
- Analysis workflow summary
- Expected results
- Troubleshooting tips
- References and resources

**Best Practices:**
- First file users should read
- Links to all other documentation
- Badges for dependencies and status

---

### 5. `INSTALLATION.md`
**Type:** Markdown  
**Purpose:** Comprehensive installation guide  
**Contents:**
- System requirements
- Platform-specific instructions (Windows, macOS, Linux)
- Step-by-step R and RStudio installation
- Package installation (interactive and batch)
- Verification scripts
- Troubleshooting common errors

**Key Sections:**
- Prerequisites
- R installation for each OS
- RStudio setup
- System dependencies
- Package installation methods
- Verification procedures

---

### 6. `DATA_MANIFEST.md`
**Type:** Markdown  
**Purpose:** Data file documentation and acquisition  
**Contents:**
- Description of required data files
- File formats and structure
- Download instructions
- Alternative data sources
- Data validation procedures
- Troubleshooting data issues

**Data Sources:**
- 10x Genomics PBMC dataset
- Reference dataset for annotation
- Alternative datasets (optional)

---

## ⚙️ Configuration Files

### 7. `config.R`
**Type:** R Configuration Script  
**Purpose:** Centralized parameter and path management  
**Key Features:**
- All paths in one place
- Adjustable analysis parameters
- QC thresholds
- Helper functions for setup
- Easy customization without editing main scripts

**Configuration Categories:**
- Data paths (input/output)
- QC parameters (thresholds)
- Normalization settings
- Clustering parameters
- Visualization options
- Marker gene lists

**How to Use:**
```r
# At the start of your analysis
source("config.R")

# Modify parameters
QC_PARAMS$max_mito <- 10  # Adjust as needed

# Use in analysis
pbmc <- subset(pbmc, subset = percent.mito < QC_PARAMS$max_mito)
```

---

## 📊 Expected Output Structure

After running the complete analysis, your repository will have:

```
scRNA-seq-analysis/
│
├── README.md
├── INSTALLATION.md
├── DATA_MANIFEST.md
├── scRNA_seq_analysis.Rmd
├── scRNA_seq_analysis_main.R
├── scRNA_seq_analysis_optional.R
├── config.R
│
├── data/                           # User creates and populates
│   ├── barcodes.tsv
│   ├── genes.tsv
│   ├── matrix.mtx
│   └── reference.rds
│
├── results/                        # Auto-generated during analysis
│   ├── pbmc_filt.rds              # Processed Seurat object
│   ├── pbmc_sct.rds               # SCT-processed object (optional)
│   │
│   ├── figures/
│   │   ├── qc_plots/
│   │   │   ├── hemoglobin_vln.png
│   │   │   ├── mito_ribo_vln.png
│   │   │   ├── mito_ribo_scatter.png
│   │   │   └── count_feature_vln.png
│   │   │
│   │   ├── clustering/
│   │   │   ├── pca_plot.png
│   │   │   ├── elbow_plot.png
│   │   │   ├── cluster_tree.png
│   │   │   ├── tsne_clusters.png
│   │   │   └── umap_clusters.png
│   │   │
│   │   ├── cell_types/
│   │   │   ├── marker_features.png
│   │   │   ├── marker_violins.png
│   │   │   ├── predicted_celltypes.png
│   │   │   └── celltype_barplot.png
│   │   │
│   │   └── degs/
│   │       ├── deg_heatmap.png
│   │       ├── hierarchical_tree.png
│   │       ├── go_bp_dotplot.png
│   │       ├── go_cc_dotplot.png
│   │       ├── go_mf_dotplot.png
│   │       └── kegg_dotplot.png
│   │
│   └── tables/
│       ├── qc_metrics.csv
│       ├── cluster_markers.csv
│       ├── cell_type_predictions.csv
│       ├── upregulated_degs.csv
│       ├── downregulated_degs.csv
│       └── enrichment_results.csv
│
└── scRNA_seq_analysis.html        # Knitted HTML report
```

---

## 🔄 Workflow Relationships

```
config.R
   ↓
scRNA_seq_analysis_main.R
   ↓
scRNA_seq_analysis_optional.R

OR

config.R
   ↓
scRNA_seq_analysis.Rmd → scRNA_seq_analysis.html
```

---

## 📝 File Usage Guide

### For Beginners:
1. Read `README.md`
2. Follow `INSTALLATION.md`
3. Get data from `DATA_MANIFEST.md`
4. Use `scRNA_seq_analysis.Rmd` in RStudio

### For Experienced Users:
1. Review `README.md`
2. Edit `config.R` with your paths
3. Run `scRNA_seq_analysis_main.R`
4. Optionally run `scRNA_seq_analysis_optional.R`

### For Customization:
1. Modify parameters in `config.R`
2. Edit scripts as needed
3. Create new analysis sections
4. Update documentation

---

## 🔍 Key Path Variables

All scripts expect these paths to be defined:

### Input Data:
- **10x files:** `DATA_DIR` in config.R
  - Default: `"data/"`
  - Contains: barcodes.tsv, genes.tsv, matrix.mtx

- **Reference:** `REFERENCE_PATH` in config.R
  - Default: `"data/reference.rds"`
  - Optional but recommended

### Output:
- **Results:** `RESULTS_DIR` in config.R
  - Default: `"results/"`

- **Figures:** `FIGURE_DIR` in config.R
  - Default: `"results/figures/"`

- **Tables:** `TABLE_DIR` in config.R
  - Default: `"results/tables/"`

---

## 📦 Required Data Files Summary

| File | Location | Required? | Description |
|------|----------|-----------|-------------|
| barcodes.tsv | data/ | ✓ | Cell barcodes from 10x |
| genes.tsv | data/ | ✓ | Gene annotations |
| matrix.mtx | data/ | ✓ | Expression matrix |
| reference.rds | data/ | ○ | Annotated reference (optional) |

---

## 🎯 Quick Start Checklist

- [ ] Read README.md
- [ ] Follow INSTALLATION.md to set up R environment
- [ ] Download data per DATA_MANIFEST.md
- [ ] Create `data/` directory
- [ ] Place 10x files in `data/`
- [ ] (Optional) Download reference.rds
- [ ] Edit `config.R` with correct paths
- [ ] Run `scRNA_seq_analysis.Rmd` OR `scRNA_seq_analysis_main.R`
- [ ] Review results in `results/`

---

## 📚 Additional Resources

### Tutorial Origin:
- **Author:** Kristoffer Nilsson Grimstad
- **Institution:** University of Skövde
- **Date:** October 2024
- **Original File:** scRNA-seq_computer_lab_2024.pdf

### External Links:
- Seurat: https://satijalab.org/seurat/
- 10x Genomics: https://www.10xgenomics.com/
- Bioconductor: https://bioconductor.org/

---

## 🔧 Troubleshooting Quick Reference

| Issue | File to Check | Solution |
|-------|--------------|----------|
| Can't find data | config.R | Update DATA_DIR path |
| Package errors | INSTALLATION.md | Follow installation steps |
| Wrong output | config.R | Check parameters |
| Missing reference | DATA_MANIFEST.md | Download or skip reference step |
| Memory error | INSTALLATION.md | See troubleshooting section |

---

## 📞 Support

For issues or questions:
1. Check relevant .md documentation file
2. Review troubleshooting sections
3. Consult Seurat documentation
4. Open GitHub issue

---

**Repository Version:** 1.0  
**Last Updated:** November 2024  
**Compatibility:** R ≥ 4.1.0, Seurat < 5.0

---

## License

MIT License - See repository for details
