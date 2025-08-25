#named for: Yggdrasil, the World Tree. Sacred in Norse cosmology. Holds all
#Nine Worlds explored in Norse belief, so seemed fitting to call the script
#that holds all helper functions that let the code run after it.
#install_if_missing: download packages if not already present
#usage example: 
#install_if_missing(cran_pkgs = "dplyr")
#install_if_missing(bioc_pkgs = "DESeq2")
#install_if_missing(
# cran_pkgs = c("ggplot2", "data.table"),
# bioc_pkgs = c("edgeR", "tximport")
# )

install_if_missing <- function(cran_pkgs = NULL, bioc_pkgs = NULL) {
  #install BiocManager first if any Bioconductor packages are requested
  if (!is.null(bioc_pkgs)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      message("Installing BiocManager...")
      install.packages("BiocManager")
    }
  }
  
  #handle CRAN packages
  if (!is.null(cran_pkgs)) {
    for (pkg in as.character(cran_pkgs)) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing CRAN package: ", pkg)
        install.packages(pkg)
      } else {
        message("CRAN package already installed: ", pkg)
      }
    }
  }
  
  #handle bioconductor packages
  if (!is.null(bioc_pkgs)) {
    for (pkg in as.character(bioc_pkgs)) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing Bioconductor package: ", pkg)
        BiocManager::install(pkg)
      } else {
        message("Bioconductor package already installed: ", pkg)
      }
    }
  }
}


################################################################################
#load_packages: loads multiple packages, both CRAN and BioC
#usage example:
#load_packages(cran_pkgs = "dplyr)
#load_packages(bioc_pkgs = "DESeq2")
#load_packages(
# cran_pkgs = c("ggplot2", "data.table"),
# bioc_pkgs = c("edgeR", "tximport")
# )
load_packages <- function(cran_pkgs = NULL, bioc_pkgs = NULL) {
  #make sure we have BiocManager if Bioconductor packages are requested
  if (!is.null(bioc_pkgs)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  }
  
  #handle CRAN packages
  if (!is.null(cran_pkgs)) {
    for (pkg in cran_pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing CRAN package: ", pkg)
        install.packages(pkg)
      }
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
  
  #handle Bioconductor packages
  if (!is.null(bioc_pkgs)) {
    for (pkg in bioc_pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing Bioconductor package: ", pkg)
        BiocManager::install(pkg)
      }
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
}