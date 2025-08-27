#role: lazy load metadata for sample selection
#named for: Ratatoskr, the squirrel in the World Tree (Yggdrasil). Responsible
#for running messages from the Eagles in its branches to the serpent (Nidhoggr) 
#at its roots. Up to the end user to decide who they are in this case. 
#check for duplicated column names
#check on failed samples
#will want to do this for not just cell lines but from tissue as well and from other dbs (TCGA, GTEx?, CCLE? let's stick with one at a time for now and remember that some of these dbs will be controlled)


# -----------------------------
# 0. Setup
# -----------------------------
bioc <- c("SRAdb","dplyr")
cran <- c("rentrez", "xml2","yaml")
source("~/Children_of_Loki/RNA-seq/helper_code/Yggdrasil_helper_functions.R")

install_if_missing(bioc_pkgs = bioc, cran_pkgs = cran)
load_packages(bioc_pkgs = bioc, cran_pkgs = cran)

#Read the YAML
config <- yaml::read_yaml("~/Children_of_Loki/RNA-seq/helper_code/Mimir_config.yaml")

options(timeout = 120)  # seconds

# -----------------------------
# 1. Define variables
# -----------------------------
#make sure you check the number of IDs by passing the query value to the SRA
#website located here: https://www.ncbi.nlm.nih.gov/sra
query <- paste0(config$query, "[", names(config$query), "]", collapse = " AND ")
study_name <- gsub(" ", "_", config$indication)
search_batch <- config$search_batch
retstart <- config$retstart
batch_size <- config$batch_size
output_file_run_meta  <- file.path(config$paths$output_dir, 
                                   paste0(study_name, "_", config$paths$run_metadata))
skipped_file_sra <- file.path(config$paths$output_dir, 
                              paste0(study_name, "_", config$paths$skipped_sra_ids))
output_file_bio_meta  <- file.path(config$paths$output_dir, 
                                   paste0(study_name, "_", config$paths$biosample_metadata))
skipped_file_bio <- file.path(config$paths$output_dir, 
                              paste0(study_name, "_", config$paths$skipped_biosamples))
output_file_expr_meta <- file.path(config$paths$output_dir, 
                                   paste0(study_name, "_", config$paths$experiment_metadata))
skipped_file_expr     <- file.path(config$paths$output_dir, 
                                   paste0(study_name, "_", config$paths$skipped_experiments))


# -----------------------------
# 2. Search all matching IDs
# -----------------------------
all_ids <- c()

start_time <- Sys.time()

repeat {
  res <- try(entrez_search(db="sra",
                           term=query,
                           retmax=search_batch,
                           retstart=retstart,
                           timeout=60),
             silent = TRUE)
  
  if (inherits(res, "try-error")) { Sys.sleep(5); next }
  
  all_ids <- c(all_ids, res$ids)
  cat("Retrieved", length(all_ids), "of", res$count, "IDs\n")
  
  if (length(all_ids) >= res$count) break
  retstart <- retstart + search_batch
}

cat("Total IDs retrieved:", length(all_ids), "\n")

end_time <- Sys.time()
cat("Time to retrieve all relevant IDs in minutes: ",
    as.numeric(
      round(
        difftime(end_time,
                 start_time,
                 units = "mins"),
        2)
      ), 
    "\n")

# -----------------------------
# 3. Fetch run metadata
# -----------------------------
# Create output directory if missing
if (!dir.exists(config$paths$output_dir)) {
  dir.create(config$paths$output_dir, recursive = TRUE)
}

# Initialize skipped list
if (file.exists(skipped_file_sra)) {
  skipped_sra_ids <- read.csv(skipped_file_sra, stringsAsFactors = FALSE)$skipped_sra
} else {
  skipped_sra_ids <- character()
}

# Determine IDs to process
if (file.exists(output_file_run_meta)) {
  processed <- read.csv(output_file_run_meta, stringsAsFactors = FALSE)
  processed_ids <- unique(processed$run_acc)
  remaining_ids <- setdiff(all_ids, processed_ids)
  message("Resuming run metadata collection. Already processed: ", length(processed_ids),
          ", Remaining: ", length(remaining_ids))
} else {
  remaining_ids <- all_ids
}

start_time <- Sys.time()

#lazy batch processing
for (start in seq(1, length(remaining_ids), by = batch_size)) {
  batch_ids <- remaining_ids[start:min(start + batch_size - 1, length(remaining_ids))]
  
  # Retry logic
  retries <- 3
  wait_sec <- 5
  res <- NULL
  for (i in seq_len(retries)) {
    res <- try(entrez_summary(db = "sra", id = batch_ids, timeout = 60), silent = TRUE)
    if (!inherits(res, "try-error")) break
    message("Retry ", i, " for IDs: ", paste(batch_ids, collapse = ", "))
    Sys.sleep(wait_sec)
  }
  
  # If still failing, record skipped IDs and continue
  if (inherits(res, "try-error")) {
    message("Failed to fetch IDs: ", paste(batch_ids, collapse = ", "))
    skipped_sra_ids <- c(skipped_sra_ids, batch_ids)
    write.csv(data.frame(skipped_sra = skipped_sra_ids), skipped_file_sra, row.names = FALSE)
    next
  }
  
  # Parse current batch
  batch_df <- bind_rows(lapply(res, function(x) {
    if (is.null(x$expxml)) return(NULL)
    xml_doc <- read_xml(paste0("<root>", x$expxml, "</root>"))
    
    data.frame(
      run_acc        = xml_attr(xml_find_first(xml_doc, ".//Run"), "acc"),
      exp_acc        = xml_attr(xml_find_first(xml_doc, ".//Experiment"), "acc"),
      study_acc      = xml_attr(xml_find_first(xml_doc, ".//Study"), "acc"),
      title          = xml_text(xml_find_first(xml_doc, ".//Title")),
      library_name   = xml_text(xml_find_first(xml_doc, ".//LIBRARY_NAME")),
      library_strategy = xml_text(xml_find_first(xml_doc, ".//LIBRARY_STRATEGY")),
      library_source   = xml_text(xml_find_first(xml_doc, ".//LIBRARY_SOURCE")),
      library_layout   = ifelse(length(xml_find_all(xml_doc, ".//PAIRED")) > 0, "PAIRED", "SINGLE"),
      instrument       = xml_attr(xml_find_first(xml_doc, ".//Platform"), "instrument_model"),
      total_spots      = xml_attr(xml_find_first(xml_doc, ".//Statistics"), "total_spots"),
      total_bases      = xml_attr(xml_find_first(xml_doc, ".//Statistics"), "total_bases"),
      createdate       = x$createdate,
      updatedate       = x$updatedate,
      submitter        = xml_attr(xml_find_first(xml_doc, ".//Submitter"), "acc"),
      bioproject       = xml_text(xml_find_first(xml_doc, ".//Bioproject")),
      biosample        = xml_text(xml_find_first(xml_doc, ".//Biosample")),
      lib_protocol     = xml_text(xml_find_first(xml_doc, ".//LIBRARY_CONSTRUCTION_PROTOCOL")),
      stringsAsFactors = FALSE
    )
  }))
  
  # Append batch to CSV incrementally
  if (nrow(batch_df) > 0) {
    write.table(
      batch_df,
      output_file_run_meta,
      sep = ",",
      row.names = FALSE,
      col.names = !file.exists(output_file_run_meta),
      append = TRUE
    )
    message("Processed batch ", start, " to ", min(start + batch_size - 1, length(remaining_ids)),
            " — rows: ", nrow(batch_df))
  }
  
  Sys.sleep(0.34) # Respect NCBI rate limits
}

#get summary
end_time <- Sys.time()
cat("Time to fetch run metadata in minutes:",
    round(as.numeric(difftime(end_time, start_time, units = "mins")), 2), "\n")

if (length(skipped_sra_ids) > 0) {
  cat("Skipped SRA IDs saved to:", skipped_file_sra, "\n")
}


# -----------------------------
# 4. Fetch BioSample metadata (batch, resume-safe, vector-validated)
# -----------------------------
# Ensure output directory exists
if (!dir.exists(config$paths$output_dir)) {
  dir.create(config$paths$output_dir, recursive = TRUE)
}

# Load run metadata
run_metadata_df <- read.csv(output_file_run_meta, stringsAsFactors = FALSE)
biosample_ids <- toupper(trimws(run_metadata_df$biosample))
biosample_ids <- biosample_ids[!is.na(biosample_ids) & biosample_ids != ""]

# Load previously fetched metadata
if (file.exists(output_file_bio_meta)) {
  biosample_dedup <- read.csv(output_file_bio_meta, stringsAsFactors = FALSE) %>%
    distinct(BioSample_accession, .keep_all = TRUE)
  processed_ids <- toupper(trimws(biosample_dedup$BioSample_accession))
} else {
  processed_ids <- character()
}

# Load skipped IDs if any
if (file.exists(skipped_file_bio)) {
  skipped_biosamples <- read.csv(skipped_file_bio, stringsAsFactors = FALSE)$skipped_biosample
} else {
  skipped_biosamples <- character()
}

# Determine remaining IDs
remaining_ids <- setdiff(biosample_ids, processed_ids)
total_samples <- length(remaining_ids)
message("BioSample metadata: Already processed ", length(processed_ids),
        ", Remaining ", total_samples)

start_time <- Sys.time()
batch_size <- 50  # adjust for API limits

# -----------------------------
# Function: fetch and dynamically parse BioSample XML
# -----------------------------
fetch_biosample_dynamic <- function(bs_id) {
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
                "db=biosample&id=", bs_id, "&retmode=xml")
  
  xml_raw <- tryCatch({
    httr::GET(url) %>% httr::content(as="text", encoding="UTF-8") %>% read_xml()
  }, error = function(e) NULL)
  
  if (is.null(xml_raw)) return(NULL)
  
  # Extract core fields
  core_fields <- data.frame(
    BioSample_accession = xml_text(xml_find_first(xml_raw, ".//Ids/Id[@db='BioSample']")),
    Organism           = xml_text(xml_find_first(xml_raw, ".//Organism/OrganismName")),
    Title              = xml_text(xml_find_first(xml_raw, ".//Description")),
    Package            = xml_text(xml_find_first(xml_raw, ".//Package")),
    stringsAsFactors   = FALSE
  )
  
  if (is.na(core_fields$BioSample_accession) | core_fields$BioSample_accession == "") {
    core_fields$BioSample_accession <- bs_id
  }
  core_fields$BioSample_accession <- toupper(trimws(core_fields$BioSample_accession))
  
  # Extract all <Attribute> nodes dynamically
  attr_nodes <- xml_find_all(xml_raw, ".//Attributes/Attribute")
  if (length(attr_nodes) > 0) {
    attr_df <- data.frame(
      name  = xml_attr(attr_nodes, "attribute_name"),
      value = xml_text(attr_nodes),
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(name) & name != "") %>%
      pivot_wider(names_from = name, values_from = value, values_fn = function(x) paste(x, collapse="; "))
    
    # Combine core + dynamic attributes
    df <- bind_cols(core_fields, attr_df)
  } else {
    df <- core_fields
  }
  
  return(df)
}

# -----------------------------
# Batch processing loop
# -----------------------------
for (start in seq(1, total_samples, by = batch_size)) {
  batch_ids <- remaining_ids[start:min(start + batch_size - 1, total_samples)]
  message("Processing batch ", start, " to ", min(start + batch_size - 1, total_samples))
  
  batch_results <- lapply(seq_along(batch_ids), function(i) {
    bs_id <- batch_ids[i]
    message("  Sample ", start + i - 1, " of ", total_samples, ": ", bs_id)
    
    df <- tryCatch(fetch_biosample_dynamic(bs_id), error = function(e) NULL)
    
    if (is.null(df)) {
      skipped_biosamples <<- unique(c(skipped_biosamples, bs_id))
      return(NULL)
    } else {
      return(df)
    }
  })
  
  batch_df <- dplyr::bind_rows(batch_results)
  
  if (nrow(batch_df) > 0) {
    write.table(
      batch_df,
      output_file_bio_meta,
      sep = ",",
      row.names = FALSE,
      col.names = !file.exists(output_file_bio_meta),
      append = TRUE
    )
  }
  
  # Save skipped IDs incrementally
  write.csv(data.frame(skipped_biosample = skipped_biosamples),
            skipped_file_bio, row.names = FALSE)
  
  Sys.sleep(0.34)  # NCBI rate limit
}

end_time <- Sys.time()
cat("Time to fetch BioSample metadata in hours: ",
    round(as.numeric(difftime(end_time, start_time, units = "hours")), 2), "\n")
cat("Skipped BioSample IDs saved to: ", skipped_file_bio, "\n")


# -----------------------------
# 5. Fetch Experiment metadata
# -----------------------------
# Ensure output directory exists
if (!dir.exists(config$paths$output_dir)) {
  dir.create(config$paths$output_dir, recursive = TRUE)
}

# Load list of unique experiment IDs
experiment_ids <- unique(run_metadata_df$exp_acc)
experiment_ids <- experiment_ids[!is.na(experiment_ids) & experiment_ids != ""]

# Load skipped experiments if present
if (file.exists(skipped_file_expr)) {
  skipped_experiments <- read.csv(skipped_file_expr, stringsAsFactors = FALSE)$skipped_exp
} else {
  skipped_experiments <- character()
}

# Determine IDs to resume if metadata already exists
if (file.exists(output_file_expr_meta)) {
  processed <- read.csv(output_file_expr_meta, stringsAsFactors = FALSE)
  processed_ids <- unique(processed$exp_acc)
  remaining_ids <- setdiff(experiment_ids, processed_ids)
  message("Resuming Experiment metadata. Already processed: ", length(processed_ids),
          ", Remaining: ", length(remaining_ids))
} else {
  remaining_ids <- experiment_ids
}

start_time <- Sys.time()

# Lazy batch processing
for (start in seq(1, length(remaining_ids), by = batch_size)) {
  batch_ids <- remaining_ids[start:min(start + batch_size - 1, length(remaining_ids))]
  
  exp_batch <- lapply(batch_ids, function(exp_id) {
    tryCatch({
      xml_raw <- entrez_fetch(db = "sra", id = exp_id, rettype = "xml", parsed = FALSE)
      xml_doc <- xml2::read_xml(xml_raw)
      parse_experiment(xml_doc, exp_id)
    }, error = function(e) {
      message("Skipping Experiment ID ", exp_id, ": ", e$message)
      skipped_experiments <<- c(skipped_experiments, exp_id)
      write.csv(data.frame(skipped_exp = skipped_experiments),
                skipped_file_expr, row.names = FALSE)
      return(NULL)
    })
  })
  
  exp_batch_df <- dplyr::bind_rows(exp_batch)
  if (nrow(exp_batch_df) > 0) {
    write.table(
      exp_batch_df,
      output_file_expr_meta,
      sep = ",",
      row.names = FALSE,
      col.names = !file.exists(output_file_expr_meta),
      append = TRUE
    )
    message("Processed Experiment batch ", start, " to ",
            min(start + batch_size - 1, length(remaining_ids)),
            " — collected rows: ", nrow(exp_batch_df))
  }
  
  Sys.sleep(0.34)  # Respect NCBI API rate limit
}

# Summary
end_time <- Sys.time()
cat("Time to fetch Experiment metadata in hours:",
    round(as.numeric(difftime(end_time, start_time, units = "hours")), 2), "\n")

if (length(skipped_experiments) > 0) {
  cat("Skipped Experiment IDs saved to:", skipped_file_expr, "\n")
}
      