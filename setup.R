# setup.R for First Time
# Project-wide package setup

# Ensure renv is active (already happens via .Rprofile, but explicit is fine)
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::activate()

# Define required packages
packages <- readLines("requirements_R.txt")
packages <- trimws(packages)
packages <- packages[packages != "" & substr(packages,1,1) != "#"]

cran_pkgs <- c()
bioc_pkgs <- c()

for (pkg in packages) {
  if (startsWith(pkg, "CRAN:")) {
    cran_pkgs <- c(cran_pkgs, sub("CRAN:", "", pkg))
  } else if (startsWith(pkg, "BIOC:")) {
    bioc_pkgs <- c(bioc_pkgs, sub("BIOC:", "", pkg))
  }
}

# Install missing packages
install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing)) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    if (bioc && !requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    for (pkg in missing) {
      if (bioc) {
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }
  }
}

# Install Any Missing Packages
install_if_missing(cran_pkgs)
install_if_missing(bioc_pkgs, bioc = TRUE)


# Load quietly
suppressPackageStartupMessages(
  lapply(packages, library, character.only = TRUE)
)

# Save to renv.lock
renv::snapshot()
# Notify Completion
message("Project environment is ready.")})


# Setup for All Other Times:
# renv::restore()

# Reference genome setup

download_and_verify <- function(url, dest, checksum = NULL, algo = "sha256") {
  if (!dir.exists(dirname(dest))) {
    dir.create(dirname(dest), recursive = TRUE)
  }

  if (!file.exists(dest)) {
    message("Downloading: ", basename(dest))
    download.file(url, dest, mode = "wb")
  } else {
    message("File already exists: ", basename(dest))
  }

  # Verify checksum if provided
  if (!is.null(checksum)) {
    # Compute local hash
    computed <- tools::md5sum(dest)
    if (computed != checksum) {
      stop("Checksum mismatch for ", dest, 
           "\nExpected: ", checksum, 
           "\nGot: ", computed)
    } else {
      message("Checksum OK for ", basename(dest))
    }
  }
}

library(jsonlite)
ref_gen_sources <- fromJSON("config_reference_genomes.json")

for (sp in names(ref_gen_sources)) {
  for (type in names(ref_gen_sources[[sp]])) {
    url <- ref_gen_sources[[sp]][[type]]$url
    checksum <- ref_gen_sources[[sp]][[type]]$checksum
    dest <- file.path("reference_genome", sp, basename(url))
    download_and_verify(url, dest, checksum)
  }
}



