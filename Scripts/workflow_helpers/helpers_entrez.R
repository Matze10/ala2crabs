# Entrez helpers for NCBI taxonomy and marker sequence counts (marker: block;
# optional taxon.ncbi_organism_filter).
# Credentials are read from the ignored Scripts/config/secrets.env file.
# Results are cached by species name so reruns do not repeat queries.

# Read simple KEY=value entries without printing credentials.
read_secrets_env <- function(path) {
  path <- scalar_secret(path)
  if (!isTRUE(nzchar(path)) || !isTRUE(file.exists(path))) return(list())
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- trimws(lines[!grepl("^\\s*#|^\\s*$", lines)])
  entries <- strsplit(lines, "=", fixed = TRUE)
  values <- lapply(entries, function(x) {
    value <- paste(x[-1], collapse = "=")
    # Strip surrounding quotes. (The old pattern had no second group, so a
    # quoted value became an empty string and the API key was silently lost.)
    sub('^(["\'])(.*)\\1$', "\\2", trimws(value), perl = TRUE)
  })
  stats::setNames(values, vapply(entries, function(x) trimws(x[[1]]), character(1)))
}

# Convert missing or zero-length values to one safe scalar string.
scalar_secret <- function(value) {
  if (isTRUE(length(value) == 0L) || isTRUE(is.na(value[[1]]))) {
    ""
  } else {
    as.character(value[[1]])
  }
}

# First element of a config value as a non-negative number, else `default`.
scalar_number <- function(value, default) {
  if (isTRUE(length(value) == 0L) || isTRUE(is.na(value[[1]]))) {
    return(default)
  }
  value <- suppressWarnings(as.numeric(value[[1]]))
  if (isTRUE(length(value) == 0L) || isTRUE(is.na(value)) || value < 0) default else value
}

# Retrieve credentials from the configured secrets file or environment.
entrez_credentials <- function(config) {
  taxonomy <- if (is.list(config) && is.list(config$taxonomy)) config$taxonomy else list()
  secret_file <- scalar_secret(taxonomy$secrets_file)
  if (isTRUE(nzchar(secret_file))) secret_file <- here::here(secret_file)
  secrets <- read_secrets_env(secret_file)
  email <- scalar_secret(Sys.getenv("NCBI_EMAIL", unset = ""))
  api_key <- scalar_secret(Sys.getenv("NCBI_API_KEY", unset = ""))
  if (!isTRUE(nzchar(email))) email <- scalar_secret(secrets$NCBI_EMAIL)
  if (!isTRUE(nzchar(api_key))) api_key <- scalar_secret(secrets$NCBI_API_KEY)
  list(email = email, api_key = api_key)
}

# Execute one Entrez search with optional NCBI credentials.
entrez_get <- function(term, db, email = NULL, api_key = NULL) {
  if (!requireNamespace("rentrez", quietly = TRUE)) {
    stop("Install the rentrez package before running Entrez checks.")
  }
  args <- list(db = db, term = term, retmode = "json")
  email <- scalar_secret(email)
  api_key <- scalar_secret(api_key)
  if (nzchar(email)) args$email <- email
  if (nzchar(api_key)) args$api_key <- api_key
  do.call(rentrez::entrez_search, args)
}

# Retrieve scientific names for NCBI Taxonomy identifiers.
entrez_taxonomy_name <- function(ids, email = NULL, api_key = NULL) {
  if (length(ids) == 0) return(NA_character_)
  args <- list(db = "taxonomy", id = ids)
  email <- scalar_secret(email)
  api_key <- scalar_secret(api_key)
  if (nzchar(email)) args$email <- email
  if (nzchar(api_key)) args$api_key <- api_key
  summaries <- tryCatch(do.call(rentrez::entrez_summary, args), error = function(e) NULL)
  if (is.null(summaries)) return(NA_character_)
  names_found <- vapply(summaries, function(x) {
    if (is.list(x) && !is.null(x$scientificname)) as.character(x$scientificname) else NA_character_
  }, character(1))
  paste(unique(names_found[!is.na(names_found)]), collapse = "; ")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# Barcode marker used for the NCBI sequence count (config `marker:` block).
# `entrez_query` is the nuccore search term ANDed with "<species>"[Organism].
# There is no default: the pipeline stops if the marker is not configured
# (checked first, in read_config()).
marker_settings <- function(config) {
  m <- if (is.list(config) && is.list(config$marker)) config$marker else list()
  name <- scalar_secret(m$name)
  query <- scalar_secret(m$entrez_query)
  if (!nzchar(trimws(name)) || !nzchar(trimws(query))) {
    stop("No barcode marker configured. Set both marker.name and marker.entrez_query in ",
         "Scripts/config/ala_crabs.yml (see the comment there and the README chapter ",
         "'Barcode marker and the NCBI query').", call. = FALSE)
  }
  # taxon.ncbi_organism_filter (default "") is an optional extra term, e.g.
  # Insecta[Organism]; empty means the query is exactly as before.
  list(name = trimws(name), entrez_query = trimws(query),
       organism_filter = taxon_settings(config)$ncbi_organism_filter)
}

# Full nuccore term for one species: "<species>"[Organism] AND (<marker query>)
# [AND (<organism filter>)].
entrez_marker_term <- function(species_name, marker) {
  term <- paste0('"', species_name, '"[Organism] AND (', marker$entrez_query, ')')
  filt <- if (is.null(marker$organism_filter)) "" else marker$organism_filter
  if (nzchar(filt)) term <- paste0(term, " AND (", filt, ")")
  term
}

# Cache file per species AND marker query. Results cached before the marker
# block existed were made with this ITS query and a name-only key; reusing
# that key for the identical query keeps them valid.
legacy_its_query <- "ITS[Title] OR internal transcribed spacer[All Fields]"

# Cache file name (MD5) for one species under the configured marker. Without an
# organism filter the key is unchanged from before, so existing cache files are
# reused; with one, the filter is part of the key.
entrez_cache_key <- function(species_name, marker) {
  filt <- if (is.null(marker$organism_filter)) "" else marker$organism_filter
  if (nzchar(filt)) return(digest::digest(paste(species_name, marker$entrez_query, filt, sep = "\r")))
  if (identical(marker$entrez_query, legacy_its_query)) digest::digest(species_name)
  else digest::digest(paste(species_name, marker$entrez_query, sep = "\r"))
}

# Query NCBI Taxonomy and nucleotide records for one reviewed species.
entrez_taxonomy_marker <- function(species_name, marker, credentials) {
  tax <- tryCatch(
    entrez_get(paste0('"', species_name, '"[Scientific Name]'), "taxonomy", credentials$email, credentials$api_key),
    error = function(e) NULL
  )
  hits <- tryCatch(
    entrez_get(entrez_marker_term(species_name, marker), "nuccore",
               credentials$email, credentials$api_key),
    error = function(e) NULL
  )
  data.table::data.table(
    species_name = species_name,
    ncbi_taxonomy_ids = if (is.null(tax)) NA_character_ else paste(tax$ids, collapse = ";"),
    ncbi_taxonomy_hits = if (is.null(tax)) NA_integer_ else tax$count,
    ncbi_scientific_names = if (is.null(tax)) NA_character_ else entrez_taxonomy_name(tax$ids, credentials$email, credentials$api_key),
    marker = marker$name,
    ncbi_marker_hits = if (is.null(hits)) NA_integer_ else hits$count,
    entrez_status = if (is.null(tax) || is.null(hits)) "query_error" else "queried",
    query_date = as.character(Sys.Date())
  )
}

# Cached results from before the marker block used the column ncbi_its_hits.
read_entrez_cache <- function(cache_file, marker) {
  res <- data.table::as.data.table(readRDS(cache_file))
  if ("ncbi_its_hits" %in% names(res) && !"ncbi_marker_hits" %in% names(res)) {
    data.table::setnames(res, "ncbi_its_hits", "ncbi_marker_hits")
  }
  if (!"marker" %in% names(res)) res[, marker := marker$name]
  res
}

# Run or load cached Entrez checks for every reviewed species.
run_entrez_checks <- function(species_file, config, output_path) {
  if (!isTRUE(file.exists(species_file))) {
    stop("Reviewed species list does not exist: ", species_file)
  }
  credentials <- entrez_credentials(config)
  species_lines <- readLines(species_file, warn = FALSE, encoding = "UTF-8")
  species <- unique(trimws(as.character(species_lines)))
  species <- species[!is.na(species) & nzchar(species)]
  if (length(species) == 0L) {
    results <- data.table::data.table(
      species_name = character(), ncbi_taxonomy_ids = character(), ncbi_taxonomy_hits = integer(),
      ncbi_scientific_names = character(), marker = character(), ncbi_marker_hits = integer(),
      entrez_status = character(), query_date = character()
    )
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(results, output_path)
    return(output_path)
  }
  taxonomy <- if (is.list(config) && is.list(config$taxonomy)) config$taxonomy else list()
  cache_relative <- scalar_secret(taxonomy$entrez_cache)
  if (!isTRUE(nzchar(cache_relative))) cache_relative <- "Scripts/data/derived/entrez_cache"
  cache_dir <- here::here(cache_relative)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  interval <- scalar_number(config$taxonomy$request_interval_seconds, 0.34)
  marker <- marker_settings(config)
  query_one <- function(name) {
    cache_file <- file.path(cache_dir, paste0(entrez_cache_key(name, marker), ".rds"))
    if (isTRUE(file.exists(cache_file))) return(read_entrez_cache(cache_file, marker))
    result <- entrez_taxonomy_marker(name, marker, credentials)
    # Failed queries are not cached, so the next run retries them.
    if (identical(result$entrez_status, "queried")) saveRDS(result, cache_file)
    Sys.sleep(interval)
    result
  }
  cached <- file.exists(file.path(cache_dir, paste0(vapply(species, entrez_cache_key, character(1), marker = marker), ".rds")))
  todo <- species[!cached]
  # Each name needs about three NCBI requests and one query takes ~2 s, so a
  # single worker manages ~0.5 names/s. Several workers stay under NCBI's
  # 10 requests/s limit with an API key (3/s without one: one worker).
  # Optional taxonomy.entrez_workers in the config; default 4 with an API key.
  workers <- as.integer(scalar_number(taxonomy$entrez_workers, 4))
  if (!nzchar(scalar_secret(credentials$api_key))) workers <- 1L
  message("Entrez (", marker$name, ", ", taxon_settings(config)$label, "): ", sum(cached), " names cached, ", length(todo), " to query with ", workers, " worker(s).")
  if (length(todo) > 0L && workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    worker_env <- environment()
    parallel::clusterExport(cl, c(
      "query_one", "cache_dir", "config", "credentials", "interval", "marker",
      "entrez_taxonomy_marker", "entrez_get", "entrez_taxonomy_name", "scalar_secret",
      "entrez_cache_key", "read_entrez_cache", "legacy_its_query", "entrez_marker_term"
    ), envir = worker_env)
    parallel::clusterEvalQ(cl, { library(data.table); NULL })
    invisible(parallel::parLapplyLB(cl, todo, query_one))
  }
  results <- data.table::rbindlist(lapply(species, query_one), fill = TRUE)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(results, output_path)
  output_path
}
