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

# QC summary (target qc_file): one row per metric with a section and a short
# description, from the files the pipeline has already written. Counts are
# names or species, not records, unless the metric says records.
count_file_rows <- function(path, chunk = 64 * 1024^2) {
  if (is.null(path) || !length(path) || !file.exists(path)) return(NA_real_)
  con <- file(path, "rb"); on.exit(close(con)); n <- 0
  repeat { b <- readBin(con, "raw", chunk); if (!length(b)) break; n <- n + sum(b == as.raw(10L)) }
  max(n - 1, 0)  # minus header
}

build_qc_summary <- function(raw_qc, target_records, candidates, normalised_names, taxonomy_matches,
                             synonym_file, genus_expansion_file, entrez_checks, config, output_path) {
  rd <- function(p) if (!is.null(p) && length(p) && file.exists(p)) data.table::fread(p, na.strings = c("", "NA"), encoding = "UTF-8") else data.table::data.table()
  out_path <- function(key) here::here(config$output[[key]])
  col <- function(dt, name) if (name %in% names(dt)) dt[[name]] else rep(NA, nrow(dt))
  rows <- list()
  add <- function(section, metric, value, description) {
    rows[[length(rows) + 1L]] <<- data.table::data.table(section = section, metric = metric,
                                                         value = as.numeric(value), description = description)
  }
  raw_n <- raw_qc$value[raw_qc$metric == "raw_records"]
  add("records", "raw_records", raw_n, "ALA records imported from the ZIP")
  add("records", "target_records", count_file_rows(target_records), "records left after the taxon, rank and genus filters")
  add("records", "target_columns", raw_qc$value[raw_qc$metric == "target_columns"], "columns kept in the target records")

  norm <- rd(normalised_names)
  add("names", "candidate_names", count_file_rows(candidates), "distinct name strings taken from the ALA name fields")
  add("names", "normalised_names", data.table::uniqueN(col(norm, "normalised_name")), "distinct names after cleaning")
  m <- rd(taxonomy_matches)
  dec <- col(m, "decision"); inc_values <- include_decision_values(config)
  included <- !is.na(dec) & dec %chin% inc_values
  add("names", "canonical_names", nrow(m), "canonical names (authors removed, rank markers standardised)")
  add("names", "canonical_species", data.table::uniqueN(col(m, "canonical_binomial")[included & !is.na(col(m, "canonical_binomial"))]),
      "distinct binomials among the included names")

  src <- col(m, "decision_source")
  add("decisions", "included_automatic", sum(included & src %chin% "automatic"), "clean names included without review")
  add("decisions", "decided_manually", sum(src %chin% "manual"), "names with a decision in taxonomy_decisions.csv")
  add("decisions", "included_total", sum(included), "names with an include decision")
  add("decisions", "excluded", sum(!is.na(dec) & !included), "names with exclude or insufficient_evidence (exclusions.csv)")
  add("decisions", "undecided", count_file_rows(out_path("crabs_review_queue")), "names still without a decision (crabs_review_queue.csv)")
  add("decisions", "genus_level_included", sum(included & col(m, "canonical_rank") %chin% "genus_or_uncertain"),
      "included names at genus level (not written as genera; see genus expansion)")

  st <- col(m, "authority_status")
  add("authority", "authority_accepted", sum(grepl("^accepted", st)), paste0("names accepted in ", authority_name(config)))
  add("authority", "authority_synonym", sum(grepl("synonym|misapplied", st)), "names that are synonyms or misapplications, replaced by the accepted name")
  add("authority", "not_in_authority", sum(st %chin% "not_in_authority"), "names not found in the authority (ALA name used)")
  add("authority", "ncbi_matched", sum(!is.na(col(m, "ncbi_taxid"))), "names matched to an NCBI Taxonomy ID")
  add("authority", "ncbi_unmatched", count_file_rows(out_path("ncbi_unmatched")), "names without an NCBI match (ncbi_unmatched.csv)")

  ge <- rd(genus_expansion_file)
  gsrc <- col(ge, "source")
  add("genus_expansion", "authority_species_added", sum(gsrc %chin% "authority_genus_expansion"), "authority species of included genus-level names")
  add("genus_expansion", "ncbi_phrase_names_added", sum(gsrc %chin% "ncbi_phrase_name"), "phrase names that NCBI lists as taxa")

  e <- rd(entrez_checks)
  hits <- suppressWarnings(as.numeric(col(e, "ncbi_marker_hits")))
  marker <- if (!is.null(config$marker$name)) config$marker$name else "marker"
  add("ncbi_marker", "species_checked", nrow(e), "names sent to the NCBI nucleotide check")
  add("ncbi_marker", "with_marker_sequences", sum(!is.na(hits) & hits > 0), paste0("names with at least one ", marker, " sequence"))
  add("ncbi_marker", "without_marker_sequences", sum(!is.na(hits) & hits == 0), paste0("names with no ", marker, " sequence"))
  add("ncbi_marker", "failed_lookups", sum(!col(e, "entrez_status") %chin% "queried"), "lookups that failed (retried on the next run)")

  crabs_file <- out_path("crabs_species_list")
  lines <- if (file.exists(crabs_file)) trimws(readLines(crabs_file, warn = FALSE, encoding = "UTF-8")) else character()
  lines <- lines[nzchar(lines)]
  hybrid <- grepl("^[^ ]+ [x×] |^[x×] ", lines)
  phrase <- !hybrid & grepl(" sp\\. ", lines)
  add("crabs_list", "species_lines", length(lines), "lines in crabs_species_list.txt")
  add("crabs_list", "binomials", sum(!hybrid & !phrase), "plain binomials")
  add("crabs_list", "named_hybrids", sum(hybrid), "hybrids NCBI lists as species (crabs.keep_ncbi_hybrids)")
  add("crabs_list", "phrase_names", sum(phrase), "phrase names NCBI lists as taxa (crabs.keep_ncbi_phrase_names)")
  add("crabs_list", "with_marker_sequences", sum(lines %chin% col(e, "species_name")[!is.na(hits) & hits > 0]),
      paste0("list lines with at least one ", marker, " sequence"))

  rf <- rd(out_path("region_flagged_species"))
  regions <- sub("^n_inside_", "", grep("^n_inside_", names(rf), value = TRUE))
  for (r in regions) {
    flagged <- if (r %in% names(rf)) as.logical(rf[[r]]) else rep(FALSE, nrow(rf))
    flagged[is.na(flagged)] <- FALSE
    in_list <- as.logical(col(rf, "in_crabs_list")); in_list[is.na(in_list)] <- FALSE
    add("regions", paste0(r, "_flagged"), sum(flagged), paste0("species with all georeferenced records in region ", r))
    add("regions", paste0(r, "_in_list"), sum(flagged & in_list), paste0("of these, species in the CRABS list"))
  }

  qc <- data.table::rbindlist(rows)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(qc, output_path)
  output_path
}
