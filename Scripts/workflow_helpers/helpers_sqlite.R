# SQLite import and target-record filtering helpers.
# The raw database preserves selected occurrence-level ALA fields so that
# downstream taxonomic decisions remain traceable to source records.

# Select the ALA fields required for taxonomy, filtering, and provenance.
selected_ala_columns <- function(reference_header) {
  wanted <- c(
    "raw_scientificName", "scientificName", "scientificNameAuthorship",
    "scientificNameID", "taxonID", "taxonRank", "taxonomicStatus",
    "acceptedNameUsage", "acceptedNameUsageID", "parentNameUsage",
    "parentNameUsageID", "kingdom", "phylum", "class", "order",
    "superfamily", "family", "subfamily", "genus", "species",
    "specificEpithet", "infraspecificEpithet", "verbatimTaxonRank",
    "originalNameUsage", "country", "stateProvince", "decimalLatitude",
    "decimalLongitude", "coordinateUncertaintyInMeters", "locality",
    "minimumDepthInMeters", "maximumDepthInMeters", "occurrenceID",
    "dataResourceUid", "dataResourceName", "datasetName", "basisOfRecord",
    "eventDate", "identifiedBy", "identificationQualifier",
    "identificationVerificationStatus"
  )
  unique(wanted[wanted %in% reference_header])
}

# Import both NSW and QLD occurrence archives into a new SQLite database.
import_ala_sqlite <- function(config, output_path, selected_columns) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(output_path)) unlink(output_path)
  archive <- configured_archive(config)
  zipfile <- here::here(archive$zip)
  members <- resolve_zip_members(zipfile, archive$expected_members)
  con <- DBI::dbConnect(RSQLite::SQLite(), output_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbExecute(con, "PRAGMA journal_mode = WAL")
  DBI::dbExecute(con, "PRAGMA synchronous = NORMAL")
  for (i in seq_along(members)) {
    message("Importing ", i, "/", length(members), ": ", basename(members[[i]]))
    records <- read_zip_csv(zipfile, members[[i]], select = selected_columns)
    DBI::dbWriteTable(con, "ala_occurrences_raw", records, append = i > 1L, row.names = FALSE)
    rm(records)
    gc()
  }
  DBI::dbExecute(con, "CREATE INDEX idx_raw_phylum ON ala_occurrences_raw (phylum)")
  DBI::dbExecute(con, "CREATE INDEX idx_raw_occurrence ON ala_occurrences_raw (occurrenceID)")
  DBI::dbExecute(con, "ANALYZE")
  output_path
}

# Read a complete SQLite table into R when a smaller table is appropriate.
read_sqlite_table <- function(database, table) {
  con <- DBI::dbConnect(RSQLite::SQLite(), database)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbReadTable(con, table)
}

# Retain rank-resolved records of the configured kingdom (taxon.ala_kingdom)
# and exclude the configured phyla, classes, orders and genera.
# ALA phylum/class values are not trusted as sole inclusion criteria.
filter_target_records <- function(database, accepted_ranks, exclude_phyla, exclude_classes, output_path,
                                  exclude_orders = character(), exclude_genera = character(), config = NULL) {
  con <- DBI::dbConnect(RSQLite::SQLite(), database)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  rank_sql <- paste(DBI::dbQuoteString(con, as.character(unlist(accepted_ranks))), collapse = ", ")
  phylum_sql <- paste(DBI::dbQuoteString(con, as.character(unlist(exclude_phyla))), collapse = ", ")
  class_sql <- paste(DBI::dbQuoteString(con, as.character(unlist(exclude_classes))), collapse = ", ")
  sql <- paste0(
    "SELECT * FROM ala_occurrences_raw WHERE taxonRank IN (", rank_sql,
    ")", kingdom_sql(con, config, prefix = ""),  # taxon.ala_kingdom (default Plantae)
    " AND (phylum IS NULL OR TRIM(phylum) = '' OR phylum NOT IN (", phylum_sql, "))",
    " AND (class IS NULL OR TRIM(class) = '' OR class NOT IN (", class_sql, "))",
    seed_plant_sql(con, exclude_orders, exclude_genera, prefix = "")
  )
  records <- DBI::dbGetQuery(con, sql)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(records, output_path)
  output_path
}

# SQL clauses that drop ferns and lycophytes. ALA places all land plants in
# Charophyta / Equisetopsida, so these groups can only be removed by order,
# or by genus when the order field is blank. `order` is quoted because it is
# an SQL keyword.
seed_plant_sql <- function(con, exclude_orders, exclude_genera, prefix = "") {
  out <- ""
  if (length(exclude_orders)) {
    out <- paste0(out, " AND (", prefix, "\"order\" IS NULL OR TRIM(", prefix, "\"order\") = '' OR ",
                  prefix, "\"order\" NOT IN (", paste(DBI::dbQuoteString(con, as.character(unlist(exclude_orders))), collapse = ", "), "))")
  }
  if (length(exclude_genera)) {
    out <- paste0(out, " AND (", prefix, "genus IS NULL OR TRIM(", prefix, "genus) = '' OR ",
                  prefix, "genus NOT IN (", paste(DBI::dbQuoteString(con, as.character(unlist(exclude_genera))), collapse = ", "), "))")
  }
  out
}

# Read the fern/lycophyte genus list named in the config (one genus per line).
read_excluded_genera <- function(config) {
  path <- config$target$exclude_genera_file
  if (is.null(path) || !nzchar(path)) return(character())
  genera <- trimws(readLines(here::here(path), warn = FALSE))
  unique(genera[nzchar(genera) & !startsWith(genera, "#")])
}

# Produce basic raw-database and target-file quality-control metrics.
qc_summary <- function(database, target_file) {
  con <- DBI::dbConnect(RSQLite::SQLite(), database)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  raw_n <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM ala_occurrences_raw")$n
  target <- data.table::fread(target_file, nrows = 0)
  data.table::data.table(
    metric = c("raw_records", "target_columns"),
    value = c(raw_n, ncol(target))
  )
}
