#role: harmonizing metadata and selecting which runs to use for downstream analysis
#named after: Forseti, one of the gods of justice in Norse belief. His primary role
#was that of a judge/mediator, capable of resolving disputes. Since here is where
#samples are selected and metadata is harmonized, I thought this was a decent name. 
bioc <- c("SRAdb","dplyr")
cran <- c("rentrez", "xml2")
source("~/Yggdrasil_helper_functions.R")

install_if_missing(bioc_pkgs = bioc, cran_pkgs = cran)
load_packages(bioc_pkgs = bioc, cran_pkgs = cran)
