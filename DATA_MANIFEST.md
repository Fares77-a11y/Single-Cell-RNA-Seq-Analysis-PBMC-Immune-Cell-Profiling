# Data Manifest for scRNA-seq Analysis

This document describes all data files required for the analysis.

---

## Required Data Files

### 1. 10x Genomics PBMC Dataset

**Source:** 10x Genomics Public Datasets  
**Dataset:** 3k PBMCs from a Healthy Donor  
**URL:** https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k

#### Files:
```
data/
├── barcodes.tsv       # 2,700 cell barcodes
├── genes.tsv          # 32,738 gene annotations (filtered to ~13,700)
└── matrix.mtx         # Sparse gene expression matrix
```

#### File Formats:

**barcodes.tsv:**
- Tab-separated file
- One barcode per line
- No header
- Example: `AAACATACAACCAC-1`

**genes.tsv:**
- Tab-separated file
- Two columns: Gene ID (Ensembl) and Gene Symbol
- No header
- Example: `ENSG00000243485    MIR1302-10`

**matrix.mtx:**
- MatrixMarket format (sparse matrix)
- Header line with dimensions
- Entries: gene_index cell_index count
- Example:
  ```
  %%MatrixMarket matrix coordinate integer general
  32738 2700 2286884
  32709 1 4
  32707 1 1
  ```

---

### 2. Reference Dataset for Cell Type Annotation

**Filename:** `reference.rds`  
**Format:** R Data Serialization (RDS)  
**Contents:** Pre-processed Seurat object with annotated cell types

#### Required Metadata:
- `cell_type` column with cell type annotations
- Normalized expression data
- Variable features identified

#### Expected Cell Types in Reference:
- CD4+ T cells
- CD8+ T cells
- B cells
- NK cells
- Classical monocytes
- Non-classical monocytes
- Conventional dendritic cells
- Plasma cells

#### Creating Your Own Reference:

If you don't have the reference file, you can:

1. **Use a public annotated dataset:**
   - Download from cellxgene: https://cellxgene.cziscience.com/
   - Use Human Cell Atlas data: https://www.humancellatlas.org/

2. **Process it into reference format:**
```r
# Example code to prepare reference
library(Seurat)

# Load your annotated dataset
ref_data <- Read10X(data.dir = "path/to/reference/")
reference <- CreateSeuratObject(counts = ref_data)

# Add cell type annotations (from your metadata)
reference$cell_type <- your_cell_type_annotations

# Process
reference <- reference %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)

# Save
saveRDS(reference, "data/reference.rds")
```

---

## Alternative: Run Without Reference

You can skip reference-based annotation and use only manual annotation based on known markers:

**Marker Genes:**
- **B cells:** CD19, MS4A1 (CD20), CD79A, CD79B
- **T cells:** CD3D, CD3E, CD3G
  - **CD4+ T cells:** CD4, IL7R
  - **CD8+ T cells:** CD8A, CD8B
- **NK cells:** GNLY, NKG7, NCAM1 (CD56)
- **Monocytes:** CD14, LYZ
  - **Classical:** CD14+, FCGR3A (CD16)-
  - **Non-classical:** CD14+, FCGR3A (CD16)+
- **Dendritic cells:** FCER1A, CST3

---

## Data Placement

Ensure all files are placed in the `data/` directory:

```bash
scRNA-seq-analysis/
└── data/
    ├── barcodes.tsv
    ├── genes.tsv
    ├── matrix.mtx
    └── reference.rds  # Optional
```

---

## Download Instructions

### Option 1: 10x Genomics Website

1. Visit: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k
2. Download "Feature / cell matrix (filtered)"
3. Extract files to `data/` directory

### Option 2: Command Line (Linux/Mac)

```bash
# Create data directory
mkdir -p data

# Download dataset
cd data
wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

# Extract
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz

# Move files to correct location
mv filtered_gene_bc_matrices/hg19/* .
rm -rf filtered_gene_bc_matrices pbmc3k_filtered_gene_bc_matrices.tar.gz
```

### Option 3: R Download

```r
# Download directly in R
dir.create("data", showWarnings = FALSE)

url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
download.file(url, destfile = "data/pbmc3k.tar.gz")

# Extract (requires external tools)
untar("data/pbmc3k.tar.gz", exdir = "data")

# Clean up
file.remove("data/pbmc3k.tar.gz")
```

---

## Data Validation

After downloading, verify your data files:

```r
library(Seurat)

# Test loading data
pbmc.data <- Read10X(data.dir = "data/")

# Check dimensions
dim(pbmc.data)
# Expected: ~32,000 genes x 2,700 cells

# Check for reference
if (file.exists("data/reference.rds")) {
  ref <- readRDS("data/reference.rds")
  print("Reference loaded successfully")
  print(table(ref$cell_type))
} else {
  print("Reference file not found - will use manual annotation only")
}
```

Expected output:
```
[1] 32738  2700
```

---

## Troubleshooting Data Issues

### Issue 1: Files not found
**Solution:** Check file paths and ensure files are in `data/` directory

### Issue 2: Wrong file format
**Solution:** Ensure files are uncompressed (.tsv and .mtx, not .gz)

### Issue 3: Permission errors
**Solution:** Check file permissions: `chmod 644 data/*`

### Issue 4: Large file size
**Solution:** 
- Matrix file is sparse format (~20 MB compressed)
- If much larger, you may have the wrong dataset
- Expected total size: ~50-100 MB uncompressed

---

## Additional Datasets (Optional)

For extended analysis or comparison, you can download additional datasets:

1. **PBMC 6k:** https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc6k
2. **PBMC 8k:** https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k
3. **Frozen PBMCs:** https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3

---

Last Updated: November 2024
