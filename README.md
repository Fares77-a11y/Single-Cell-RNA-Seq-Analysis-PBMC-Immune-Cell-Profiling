# Single-Cell RNA-Seq Analysis: PBMC Immune Cell Profiling

[![R](https://img.shields.io/badge/R-4.3+-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-4.3.0-orange.svg)](https://satijalab.org/seurat/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **Comprehensive computational pipeline for single-cell transcriptomic profiling of peripheral blood mononuclear cells (PBMCs) with automated cell type annotation and functional characterization**

---

## 🔬 Project Overview

This repository presents a production-ready bioinformatics pipeline for analyzing single-cell RNA sequencing (scRNA-seq) data from human peripheral blood. The workflow implements industry-standard computational methods to identify and characterize immune cell populations through unsupervised clustering, reference-based annotation, and pathway enrichment analysis.

### Scientific Context

Peripheral blood mononuclear cells (PBMCs) comprise diverse immune cell populations critical for immune surveillance, inflammation, and disease response. Single-cell transcriptomics enables high-resolution characterization of cellular heterogeneity, revealing distinct functional states and rare cell populations invisible to bulk RNA-seq approaches.

### Biological Significance

This analysis identifies **7 transcriptionally distinct immune cell populations** from 2,700 cells, including:
- **T lymphocytes** (CD4+ helper T cells, CD8+ cytotoxic T cells)
- **B lymphocytes** (antibody-producing cells)
- **Monocytes** (classical CD14+ and non-classical CD16+ subsets)
- **Natural Killer (NK) cells** (cytotoxic lymphocytes)

Each cluster is validated through:
- Canonical cell surface marker expression (CD3E, CD19, CD14, etc.)
- Reference-based label transfer from annotated datasets
- Functional pathway enrichment (Gene Ontology, KEGG)

---

## 📊 Key Results

### Cell Population Distribution

| Cluster | Cell Type | Count | Key Markers | Biological Function |
|---------|-----------|-------|-------------|---------------------|
| 0, 1, 5 | CD4+ T cells | ~1,200 | CD3E, IL7R, CD4 | Adaptive immunity, T helper response |
| 4 | CD8+ T cells | ~300 | CD3E, CD8A | Cytotoxic T cell response |
| 3 | B cells | ~350 | CD19, MS4A1 (CD20) | Antibody production |
| 2, 6 | Monocytes | ~700 | CD14, LYZ | Innate immunity, phagocytosis |
| 7 | NK cells | ~150 | GNLY, NKG7 | Cytotoxic innate immunity |

### Analysis Highlights

- **Quality Control**: Implemented multi-parameter filtering removing low-quality cells (<5% mitochondrial RNA, >20% ribosomal RNA)
- **Dimensionality Reduction**: PCA identified 10 principal components capturing biological variation while filtering technical noise
- **Cluster Optimization**: Tested 20 resolution parameters (0.1-2.0) with cluster tree stability analysis
- **Validation**: Reference-based annotation achieved >95% concordance with manual marker-based classification
- **Functional Insights**: Pathway enrichment revealed cluster-specific immune activation signatures

---

## 🛠️ Methodology

### Computational Pipeline

```
Raw 10x Data (13,714 genes × 2,700 cells)
    ↓
Quality Control & Filtering
    ↓
Normalization (LogNormalize) + Scaling
    ↓
Feature Selection (2,000 HVGs)
    ↓
Dimensionality Reduction (PCA → UMAP/t-SNE)
    ↓
Graph-Based Clustering (Louvain algorithm)
    ↓
Cell Type Annotation (Manual + Reference-based)
    ↓
Differential Expression Analysis
    ↓
Pathway Enrichment (GO, KEGG)
```

### Technical Implementation

**Core Analysis:**
- **Platform**: R 4.3+ with Seurat 4.3.0
- **Normalization**: Log-normalization with scaling factor 10,000
- **Clustering**: Shared Nearest Neighbor (SNN) graph with Louvain modularity optimization
- **Visualization**: UMAP and t-SNE for 2D projection
- **Statistics**: Wilcoxon rank-sum test for differential expression (Bonferroni-adjusted p < 0.05)

**Advanced Features:**
- SCTransform normalization for variance stabilization
- ClusterProfiler for Gene Ontology enrichment
- Reference-guided annotation with label transfer
- Hierarchical clustering for cluster relationship inference

---

## 📁 Repository Structure

```
.
├── scRNA_seq_analysis.Rmd          # Interactive RMarkdown notebook
├── scRNA_seq_analysis_main.R       # Core analysis script
├── scRNA_seq_analysis_optional.R   # Enrichment & SCTransform workflow
├── config.R                        # Configuration parameters
├── README.md                       # This file
├── INSTALLATION.md                 # Setup instructions
├── DATA_MANIFEST.md                # Data documentation
│
├── data/                           # Input data (user-provided)
│   ├── barcodes.tsv
│   ├── genes.tsv
│   ├── matrix.mtx
│   └── reference.rds
│
└── results/                        # Analysis outputs (auto-generated)
    ├── figures/
    │   ├── qc_plots/               # Quality control visualizations
    │   ├── clustering/             # PCA, UMAP, t-SNE, cluster trees
    │   ├── cell_types/             # Marker expression, annotations
    │   └── enrichment/             # GO/KEGG pathway plots
    │
    └── tables/
        ├── cluster_markers.csv     # Differentially expressed genes
        ├── cell_annotations.csv    # Cell type predictions
        └── enrichment_results.csv  # Pathway analysis results
```

---

## 🚀 Quick Start

### Prerequisites

- R ≥ 4.1.0
- RStudio (recommended)
- 8 GB RAM minimum (16 GB recommended)

### Installation

```r
# Install required packages
install.packages("devtools")
devtools::install_version("Seurat", version = "4.3.0")
install.packages(c("dplyr", "ggplot2", "clustree"))

# Optional: Enrichment analysis
BiocManager::install(c("clusterProfiler", "enrichplot", "org.Hs.eg.db"))
```

### Usage

**Interactive Analysis (Recommended):**
```r
# Open in RStudio
rmarkdown::render("scRNA_seq_analysis.Rmd")
```

**Scripted Execution:**
```r
source("scRNA_seq_analysis_main.R")
source("scRNA_seq_analysis_optional.R")  # Optional
```

**Data Requirements:**
- Download PBMC 3k dataset from [10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k)
- Place files in `data/` directory
- Update paths in `config.R`

See [INSTALLATION.md](INSTALLATION.md) for detailed setup instructions.

---

## 📈 Visualizations

### Dimensionality Reduction
![UMAP Clustering](results/figures/clustering/umap_clusters.png)
*UMAP projection showing 7 transcriptionally distinct clusters*

### Cell Type Annotation
![Cell Type Distribution](results/figures/cell_types/celltype_barplot.png)
*Distribution of predicted cell types across clusters*

### Differential Expression
![Marker Heatmap](results/figures/enrichment/marker_heatmap.png)
*Heatmap of cluster-specific marker genes*

### Pathway Enrichment
![GO Enrichment](results/figures/enrichment/go_biological_process.png)
*Gene Ontology enrichment for cluster-specific genes*

---

## 🔍 Key Findings

### 1. **Immune Cell Diversity**
Identified canonical PBMC populations with expected proportions:
- T cells (60%): Dominant population with CD4/CD8 subset separation
- Monocytes (25%): Classical and intermediate phenotypes
- B cells (13%): Mature B lymphocytes
- NK cells (2%): Small but distinct cytotoxic population

### 2. **Functional Characterization**
Pathway analysis revealed:
- **T cells**: Enriched for T cell activation, cytokine signaling
- **B cells**: Immunoglobulin production, B cell receptor signaling
- **Monocytes**: Inflammatory response, antigen presentation
- **NK cells**: Cytotoxicity, interferon gamma response

### 3. **Technical Validation**
- >95% concordance between manual and reference-based annotation
- Stable clustering across resolution parameters (0.4-0.7)
- SCTransform validation confirmed cluster robustness

---

## 💻 Skills Demonstrated

This project showcases expertise in:

- **Computational Biology**: Single-cell data analysis, dimensionality reduction, clustering
- **Bioinformatics**: RNA-seq processing, differential expression, pathway enrichment
- **Statistical Analysis**: Multiple testing correction, parameter optimization, validation
- **Programming**: R/Bioconductor, reproducible research, version control
- **Data Visualization**: Multi-panel figures, publication-quality graphics
- **Immunology**: PBMC biology, immune cell markers, functional annotation
- **Scientific Communication**: Technical documentation, result interpretation

---

## 📚 References

### Key Methods

1. **Seurat Framework**  
   Stuart et al. (2019). *Comprehensive Integration of Single-Cell Data.* Cell.  
   https://doi.org/10.1016/j.cell.2019.05.031

2. **UMAP Visualization**  
   McInnes et al. (2018). *UMAP: Uniform Manifold Approximation and Projection.*  
   https://arxiv.org/abs/1802.03426

3. **Best Practices**  
   Luecken & Theis (2019). *Current best practices in single-cell RNA-seq analysis.* Molecular Systems Biology.  
   https://doi.org/10.15252/msb.20188746

### Resources

- **Seurat Documentation**: https://satijalab.org/seurat/
- **10x Genomics Datasets**: https://support.10xgenomics.com/
- **Bioconductor**: https://bioconductor.org/books/release/OSCA/

---

## 🤝 Contributing

Contributions, issues, and feature requests are welcome! Feel free to:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Commit changes (`git commit -am 'Add new feature'`)
4. Push to branch (`git push origin feature/improvement`)
5. Open a Pull Request

---

## 📝 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

---

## 👤 Author

**Fares Ibrahim**

- GitHub: [@Fares77-a11y](https://github.com/Fares77-a11y)
- Email: Ibrahimfares825@gmail.com

---

## 🙏 Acknowledgments

- Seurat development team (Satija Lab, New York Genome Center)
- 10x Genomics for public dataset access
- Bioconductor community for computational tools

---

**Last Updated**: November 2025  
**Status**: Production-ready  
**Version**: 1.0.0
