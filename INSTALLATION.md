# Installation Guide for scRNA-seq Analysis

Complete step-by-step installation instructions for all platforms.

---

## Table of Contents

1. [System Requirements](#system-requirements)
2. [Install R and RStudio](#install-r-and-rstudio)
3. [Install R Packages](#install-r-packages)
4. [Verify Installation](#verify-installation)
5. [Troubleshooting](#troubleshooting)
6. [Platform-Specific Notes](#platform-specific-notes)

---

## System Requirements

### Minimum Requirements:
- **RAM:** 8 GB (16 GB recommended)
- **Storage:** 10 GB free space
- **Processor:** 64-bit processor
- **OS:** Windows 10+, macOS 10.13+, or Linux (Ubuntu 18.04+)

### Software Requirements:
- **R:** Version 4.1.0 or higher
- **RStudio:** Latest version (optional but recommended)
- **Internet connection** for package installation

---

## Install R and RStudio

### Windows

#### Step 1: Install R
1. Visit: https://cran.r-project.org/bin/windows/base/
2. Download the latest R installer (e.g., `R-4.3.x-win.exe`)
3. Run the installer with default settings
4. Verify installation:
   ```cmd
   R --version
   ```

#### Step 2: Install RStudio
1. Visit: https://www.rstudio.com/products/rstudio/download/
2. Download RStudio Desktop (free version)
3. Run the installer
4. Launch RStudio

#### Step 3: Install Rtools (Required for package compilation)
1. Visit: https://cran.r-project.org/bin/windows/Rtools/
2. Download Rtools matching your R version
3. Install with default settings
4. Restart RStudio

---

### macOS

#### Step 1: Install R
1. Visit: https://cran.r-project.org/bin/macosx/
2. Download the appropriate .pkg file for your macOS version
3. Open the .pkg file and follow installation instructions
4. Verify installation in Terminal:
   ```bash
   R --version
   ```

#### Step 2: Install RStudio
1. Visit: https://www.rstudio.com/products/rstudio/download/
2. Download RStudio Desktop for macOS
3. Open the .dmg file and drag RStudio to Applications
4. Launch RStudio

#### Step 3: Install XCode Command Line Tools (Required)
```bash
xcode-select --install
```

---

### Linux (Ubuntu/Debian)

#### Step 1: Install R
```bash
# Update package list
sudo apt update

# Install dependencies
sudo apt install -y software-properties-common dirmngr

# Add CRAN repository
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

# Install R
sudo apt install -y r-base r-base-dev

# Verify installation
R --version
```

#### Step 2: Install RStudio
```bash
# Download RStudio (check website for latest version)
wget https://download1.rstudio.org/electron/jammy/amd64/rstudio-2023.12.0-369-amd64.deb

# Install
sudo dpkg -i rstudio-2023.12.0-369-amd64.deb

# Fix dependencies if needed
sudo apt-get install -f
```

#### Step 3: Install System Dependencies
```bash
sudo apt install -y \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libfontconfig1-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev \
  libhdf5-dev
```

---

## Install R Packages

### Method 1: Interactive Installation (Recommended)

Open RStudio and run the following code:

```r
# ============================================================================
# STEP 1: Install BiocManager
# ============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# ============================================================================
# STEP 2: Install devtools
# ============================================================================
install.packages("devtools")

# ============================================================================
# STEP 3: Install Seurat (specific version for compatibility)
# ============================================================================
# This may take 10-20 minutes
devtools::install_version("Seurat", version = "4.3.0", 
                          repos = "http://cran.us.r-project.org")

# ============================================================================
# STEP 4: Install Core Packages
# ============================================================================
install.packages(c(
  "dplyr",
  "ggplot2",
  "clustree",
  "rmarkdown",
  "knitr"
))

# ============================================================================
# STEP 5: Install Optional Packages for Enrichment Analysis
# ============================================================================
BiocManager::install(c(
  "clusterProfiler",
  "enrichplot",
  "org.Hs.eg.db"
))

install.packages("ggtree")

# ============================================================================
# STEP 6: Verify Installation
# ============================================================================
# Load all packages to verify
library(Seurat)
library(dplyr)
library(ggplot2)
library(clustree)

# Print versions
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("dplyr version:", as.character(packageVersion("dplyr")), "\n")
cat("ggplot2 version:", as.character(packageVersion("ggplot2")), "\n")
```

### Method 2: Batch Installation Script

Save the following as `install_packages.R` and run it:

```r
#!/usr/bin/env Rscript

# Batch installation script
cat("Starting package installation...\n")

# Function to install if not already installed
install_if_missing <- function(pkg, source = "CRAN") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "from", source, "...\n")
    if (source == "CRAN") {
      install.packages(pkg, dependencies = TRUE)
    } else if (source == "Bioconductor") {
      BiocManager::install(pkg, update = FALSE)
    }
  } else {
    cat(pkg, "already installed.\n")
  }
}

# Install BiocManager first
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Core packages
core_packages <- c("devtools", "dplyr", "ggplot2", "clustree", 
                   "rmarkdown", "knitr")
for (pkg in core_packages) {
  install_if_missing(pkg, "CRAN")
}

# Install specific Seurat version
if (packageVersion("Seurat") != "4.3.0") {
  devtools::install_version("Seurat", version = "4.3.0")
}

# Bioconductor packages
bioc_packages <- c("clusterProfiler", "enrichplot", "org.Hs.eg.db")
for (pkg in bioc_packages) {
  install_if_missing(pkg, "Bioconductor")
}

# Additional packages
install_if_missing("ggtree", "CRAN")

cat("Installation complete!\n")
```

Run in terminal:
```bash
Rscript install_packages.R
```

---

## Verify Installation

### Complete Verification Script

```r
# ============================================================================
# Installation Verification Script
# ============================================================================

cat("\n=== scRNA-seq Analysis Environment Check ===\n\n")

# Check R version
r_version <- R.Version()$version.string
cat("R Version:", r_version, "\n")
if (as.numeric(R.Version()$major) >= 4 && as.numeric(R.Version()$minor) >= 1) {
  cat("✓ R version OK\n\n")
} else {
  cat("✗ R version too old. Please upgrade to R >= 4.1.0\n\n")
}

# Required packages
required_packages <- c(
  "Seurat", "dplyr", "ggplot2", "clustree", 
  "rmarkdown", "knitr"
)

# Optional packages
optional_packages <- c(
  "clusterProfiler", "enrichplot", "org.Hs.eg.db", "ggtree"
)

# Check required packages
cat("Required Packages:\n")
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- as.character(packageVersion(pkg))
    cat("  ✓", pkg, "(v", version, ")\n")
  } else {
    cat("  ✗", pkg, "NOT FOUND\n")
  }
}

# Check optional packages
cat("\nOptional Packages:\n")
for (pkg in optional_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- as.character(packageVersion(pkg))
    cat("  ✓", pkg, "(v", version, ")\n")
  } else {
    cat("  ○", pkg, "not installed (optional)\n")
  }
}

# Test Seurat functionality
cat("\nTesting Seurat functionality...\n")
tryCatch({
  library(Seurat)
  test_data <- matrix(rnbinom(1000, mu = 5, size = 2), ncol = 10)
  rownames(test_data) <- paste0("Gene", 1:100)
  colnames(test_data) <- paste0("Cell", 1:10)
  test_obj <- CreateSeuratObject(counts = test_data)
  cat("✓ Seurat is working correctly\n")
}, error = function(e) {
  cat("✗ Seurat test failed:", e$message, "\n")
})

# Memory check
memory_limit_mb <- memory.limit()  # Windows only
cat("\nSystem Memory:\n")
if (!is.na(memory_limit_mb)) {
  cat("  Available:", memory_limit_mb, "MB\n")
} else {
  cat("  Check manually with system tools\n")
}

cat("\n=== Verification Complete ===\n")
```

Expected output:
```
=== scRNA-seq Analysis Environment Check ===

R Version: R version 4.3.x (YYYY-MM-DD)
✓ R version OK

Required Packages:
  ✓ Seurat (v 4.3.0)
  ✓ dplyr (v 1.x.x)
  ✓ ggplot2 (v 3.x.x)
  ✓ clustree (v 0.x.x)
  ✓ rmarkdown (v 2.x.x)
  ✓ knitr (v 1.x.x)

Optional Packages:
  ✓ clusterProfiler (v 4.x.x)
  ✓ enrichplot (v 1.x.x)
  ✓ org.Hs.eg.db (v 3.x.x)
  ✓ ggtree (v 3.x.x)

Testing Seurat functionality...
✓ Seurat is working correctly

=== Verification Complete ===
```

---

## Troubleshooting

### Common Installation Errors

#### Error 1: "Package compilation failed"

**Windows:**
```r
# Install pre-compiled binaries
install.packages("Seurat", type = "binary")
```

**macOS/Linux:**
```bash
# Install system dependencies
# See platform-specific sections above
```

#### Error 2: "Non-zero exit status"

**Solution:**
```r
# Update all packages first
update.packages(ask = FALSE, checkBuilt = TRUE)

# Then retry installation
```

#### Error 3: "Unable to install Seurat version 4.3.0"

**Solution:**
```r
# Remove existing Seurat
remove.packages("Seurat")

# Clear cache
.libPaths()  # Check library paths
# Manually delete Seurat folder from library path

# Reinstall
devtools::install_version("Seurat", version = "4.3.0")
```

#### Error 4: "openssl/ssl.h not found" (Linux)

**Solution:**
```bash
sudo apt install libssl-dev libcurl4-openssl-dev
```

#### Error 5: Memory allocation error

**Solution:**
```r
# Increase memory limit (Windows)
memory.limit(size = 16000)

# On Linux/Mac, close other applications
```

---

## Platform-Specific Notes

### Windows

- **Administrator rights** may be required for installation
- Install **Rtools** before installing packages
- If using corporate network, may need to configure proxy:
  ```r
  Sys.setenv(http_proxy = "http://proxy.company.com:port")
  Sys.setenv(https_proxy = "https://proxy.company.com:port")
  ```

### macOS

- **XCode Command Line Tools** required
- On Apple Silicon (M1/M2):
  ```r
  # Some packages may need Rosetta
  # Install R for ARM64 architecture
  ```
- May need to allow apps from "unidentified developers" in Security settings

### Linux

- System dependencies vary by distribution
- Ubuntu/Debian: See commands in Linux section
- CentOS/RHEL:
  ```bash
  sudo yum install openssl-devel libcurl-devel libxml2-devel
  ```
- Arch Linux:
  ```bash
  sudo pacman -S r gcc-fortran
  ```

---

## Next Steps

After successful installation:

1. **Download data:** See `DATA_MANIFEST.md`
2. **Configure paths:** Edit `config.R`
3. **Run analysis:** Open `scRNA_seq_analysis.Rmd` in RStudio

---

## Getting Help

If you encounter issues:

1. **Check R version:** Must be >= 4.1.0
2. **Check package versions:** Especially Seurat (must be < 5.0)
3. **Google the error message:** Often finds solutions
4. **Seurat GitHub issues:** https://github.com/satijalab/seurat/issues
5. **Bioconductor support:** https://support.bioconductor.org/

---

Last Updated: November 2024
