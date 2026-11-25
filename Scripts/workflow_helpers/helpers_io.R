# Input and configuration helpers for the ALA-to-CRABS workflow.
# The production input is one combined NSW/QLD ALA archive (config$input).
# CSV members are streamed out of the ZIP with `unzip -p`, so the multi-GB
# archive never has to be extracted to disk.

# Read the YAML configuration. Stops immediately if the barcode marker is not
# configured, so a missing setting fails in seconds instead of after the import.
read_config <- function(path) {
  stopifnot(file.exists(path))
  config <- yaml::read_yaml(path)
  marker_settings(config)  # helpers_entrez.R; stops if marker.name / marker.entrez_query are missing
  config
}

# Which organisms the list is for (config `taxon:` block). Defaults reproduce
# the plant setup, so older configs without the block behave as before.
#   label:                text for log messages and the report
#   ala_kingdom:          ALA kingdom filter for records and ALA name matching ("" = no filter)
#   ncbi_organism_filter: optional extra NCBI nucleotide term ANDed to the marker query ("" = none)
#   ncbi_root_taxid:      NCBI subtree for name matching (old key: taxonomy.ncbi_root_taxid)
#   seed_plant_taxid:     subtree for the ncbi_seed_plant audit column; ~ to skip
taxon_settings <- function(config) {
  tx <- if (is.list(config$taxon)) config$taxon else list()
  str1 <- function(x, default) if (is.null(x) || !length(x) || is.na(x[[1]])) default else trimws(as.character(x[[1]]))
  seed <- if ("seed_plant_taxid" %in% names(tx)) tx$seed_plant_taxid else 58024
  list(
    label = str1(tx$label, "Plantae"),
    ala_kingdom = str1(tx$ala_kingdom, "Plantae"),
    ncbi_organism_filter = str1(tx$ncbi_organism_filter, ""),
    ncbi_root_taxid = str1(tx$ncbi_root_taxid %||% config$taxonomy$ncbi_root_taxid, "35493"),
    seed_plant_taxid = if (is.null(seed) || !length(seed) || is.na(seed[[1]]) || !nzchar(trimws(as.character(seed[[1]])))) NA_character_ else trimws(as.character(seed[[1]]))
  )
}

# SQL clause for the ALA kingdom filter; records with a blank kingdom are kept
# (ALA leaves it empty for some records). Empty ala_kingdom = no filter.
kingdom_sql <- function(con, config, prefix = "") {
  k <- taxon_settings(config)$ala_kingdom
  if (!nzchar(k)) return("")
  col <- paste0(prefix, "kingdom")
  paste0(" AND (", col, " = ", DBI::dbQuoteString(con, k), " OR ", col, " IS NULL OR TRIM(", col, ") = '')")
}

# Absolute path of a project-relative path (project root = folder with .here).
project_path <- function(path) here::here(path)

# File names inside a ZIP archive.
zip_members <- function(zipfile) {
  if (!is.character(zipfile) || length(zipfile) != 1L || is.na(zipfile) || !nzchar(zipfile)) {
    stop("Invalid combined ZIP path received: ", paste(dput(zipfile), collapse = ""))
  }
  zipfile <- normalizePath(zipfile, winslash = "/", mustWork = TRUE)
  utils::unzip(zipfile, list = TRUE)$Name
}

# Map each expected CSV name (config$input$expected_members) to its path inside
# the ZIP; stops unless each one occurs exactly once.
resolve_zip_members <- function(zipfile, expected) {
  actual <- zip_members(zipfile)
  expected <- as.character(expected)
  vapply(expected, function(name) {
    matches <- actual[basename(actual) == name]
    if (length(matches) != 1L) stop("Expected one ZIP member named ", name, " in ", zipfile)
    as.character(matches[[1L]])
  }, character(1))
}

# Shell command that streams one ZIP member to stdout (used by fread(cmd = ...)).
zip_command <- function(zipfile, member) {
  unzip_exe <- Sys.which("unzip")
  if (!nzchar(unzip_exe)) stop("The unzip executable is required in PATH.")
  quote_type <- if (.Platform$OS.type == "windows") "cmd" else "sh"
  paste(
    shQuote(unzip_exe, type = quote_type), "-p",
    shQuote(normalizePath(zipfile, winslash = "/", mustWork = TRUE), type = quote_type),
    shQuote(member, type = quote_type)
  )
}

# Read one CSV member straight from the ZIP (optionally only some columns or rows).
read_zip_csv <- function(zipfile, member, select = NULL, nrows = Inf, tmpdir = NULL) {
  if (is.null(tmpdir) || !nzchar(tmpdir)) tmpdir <- file.path(tempdir(), "ala_crabs_fread")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  data.table::fread(
    cmd = zip_command(zipfile, member), select = select, nrows = nrows,
    check.names = FALSE, tmpdir = normalizePath(tmpdir, winslash = "/", mustWork = TRUE)
  )
}

# Archive path and expected members from the config.
configured_archive <- function(config) {
  list(zip = config$input$combined_zip, expected_members = config$input$expected_members)
}

# One-row manifest of the archive: path, existence, size and MD5 checksum.
build_input_manifest <- function(config) {
  archive <- configured_archive(config)
  path <- here::here(archive$zip)
  data.table::data.table(
    archive = "NSW_Qld_combined", path = path, exists = file.exists(path),
    size_bytes = if (file.exists(path)) file.info(path)$size else NA_real_,
    md5 = if (file.exists(path)) as.character(tools::md5sum(path)) else NA_character_
  )
}

# Stop if the archive is missing or lacks an expected member; returns the manifest.
validate_inputs <- function(config) {
  manifest <- build_input_manifest(config)
  if (!isTRUE(manifest$exists[[1L]])) stop("Configured combined ALA archive is missing.")
  archive <- configured_archive(config)
  resolve_zip_members(here::here(archive$zip), archive$expected_members)
  manifest
}

# Read only the header row of every member and stop unless all are identical,
# so the parts can be stacked into one SQLite table.
validate_headers <- function(config) {
  archive <- configured_archive(config)
  zipfile <- here::here(archive$zip)
  members <- resolve_zip_members(zipfile, archive$expected_members)
  headers <- lapply(members, function(member) names(read_zip_csv(zipfile, member, nrows = 0)))
  reference <- headers[[1L]]
  if (!all(vapply(headers, identical, logical(1), y = reference))) {
    stop("Occurrence CSV headers are not identical in the combined archive.")
  }
  list(reference_header = reference, n_files = length(headers))
}

# Render an analysis page to HTML. rmarkdown needs pandoc, which ships with
# RStudio and Quarto but is not found when R runs from Rscript or Positron.
# Look in the usual install folders first; without pandoc, knit to Markdown
# so the run still finishes, and say how to get the HTML.
render_report <- function(input, output_dir) {
  if (!rmarkdown::pandoc_available("2.8")) {
    dirs <- c(
      Sys.getenv("RSTUDIO_PANDOC"),
      "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
      "C:/Program Files/RStudio/bin/quarto/bin/tools",
      "C:/Program Files/Quarto/bin/tools",
      file.path(Sys.getenv("LOCALAPPDATA"), "Programs/Positron/resources/app/quarto/bin/tools"),
      "C:/Program Files/Positron/resources/app/quarto/bin/tools",
      file.path(Sys.getenv("LOCALAPPDATA"), "Pandoc")
    )
    dirs <- dirs[nzchar(dirs) & dir.exists(dirs)]
    if (length(dirs)) rmarkdown::find_pandoc(dir = dirs)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (rmarkdown::pandoc_available("2.8")) {
    return(rmarkdown::render(input, output_dir = output_dir, envir = new.env(parent = globalenv()), quiet = TRUE))
  }
  warning("pandoc not found: wrote a Markdown report instead of HTML. Install pandoc or Quarto for the HTML version.")
  out <- normalizePath(file.path(output_dir, sub("\\.Rmd$", ".md", basename(input))), mustWork = FALSE)
  old <- setwd(dirname(input)); on.exit(setwd(old), add = TRUE)
  knitr::knit(basename(input), output = out, quiet = TRUE, envir = new.env(parent = globalenv()))
}
