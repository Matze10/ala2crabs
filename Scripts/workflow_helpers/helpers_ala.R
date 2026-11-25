# ALA name-matching check for genera, used to drop ferns and lycophytes.
# ALA files every land plant under Charophyta / Equisetopsida, so ferns can only
# be recognised from order, family or ALA's "Ferns and allies" species group.
# Pteridaceae (Adiantum, Cheilanthes, Pteris ...) has a blank order in ALA, which
# is why the species group is checked as well. Lookups are cached per genus.

fern_lycophyte_orders <- c(
  "Lycopodiales", "Isoetales", "Selaginellales", "Equisetales", "Psilotales",
  "Ophioglossales", "Marattiales", "Osmundales", "Hymenophyllales",
  "Gleicheniales", "Schizaeales", "Salviniales", "Cyatheales", "Polypodiales"
)

# Distinct genus values plus the first word of each raw name in the database.
raw_genera <- function(database) {
  con <- DBI::dbConnect(RSQLite::SQLite(), database)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  sql <- paste(
    "SELECT DISTINCT TRIM(genus) AS g FROM ala_occurrences_raw WHERE genus IS NOT NULL AND TRIM(genus) <> ''",
    "UNION SELECT DISTINCT CASE WHEN INSTR(TRIM(raw_scientificName) || ' ', ' ') > 1",
    "THEN SUBSTR(TRIM(raw_scientificName), 1, INSTR(TRIM(raw_scientificName) || ' ', ' ') - 1) END",
    "FROM ala_occurrences_raw WHERE raw_scientificName IS NOT NULL"
  )
  g <- DBI::dbGetQuery(con, sql)$g
  g <- sub("^!+", "", trimws(g))
  sort(unique(g[!is.na(g) & grepl("^[A-Z][a-z-]+$", g)]))
}

# Look up one genus with ALA name matching, restricted to the configured
# kingdom (taxon.ala_kingdom, default Plantae; "" = no kingdom restriction).
ala_genus_lookup <- function(genus, base_url, kingdom = "Plantae") {
  empty <- data.table::data.table(
    genus = genus, ala_success = NA, ala_rank = NA_character_, ala_order = NA_character_,
    ala_family = NA_character_, ala_species_group = NA_character_, ala_issues = NA_character_,
    ala_status = "query_error", query_date = as.character(Sys.Date())
  )
  resp <- tryCatch(
    httr::GET(base_url, query = c(list(scientificName = genus), if (nzchar(kingdom)) list(kingdom = kingdom)),
              httr::timeout(30)),
    error = function(e) NULL
  )
  if (is.null(resp) || httr::status_code(resp) != 200L) return(empty)
  j <- tryCatch(jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8")), error = function(e) NULL)
  if (is.null(j)) return(empty)
  val <- function(x) if (is.null(x) || length(x) == 0L) NA_character_ else paste(as.character(x), collapse = "|")
  data.table::data.table(
    genus = genus,
    ala_success = isTRUE(j$success),
    ala_rank = val(j$rank),
    ala_order = val(j$order),
    ala_family = val(j$family),
    ala_species_group = val(j$speciesGroup),
    ala_issues = val(j$issues),
    ala_status = "queried",
    query_date = as.character(Sys.Date())
  )
}

# Check every genus against ALA (cached) and flag ferns and lycophytes.
run_ala_genus_checks <- function(genera, config, output_path) {
  taxonomy <- if (is.list(config$taxonomy)) config$taxonomy else list()
  base_url <- taxonomy$ala_namematching_url %||%
    "https://api.ala.org.au/namematching/api/searchByClassification"
  cache_dir <- here::here(taxonomy$ala_cache %||% "Scripts/data/derived/ala_genus_cache")
  interval <- scalar_number(taxonomy$ala_request_interval_seconds, 0.05)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  fern_orders <- as.character(unlist(config$target$exclude_orders %||% fern_lycophyte_orders))
  # ALA species groups that exclude a genus (default: "Ferns and allies").
  groups <- as.character(unlist(config$target$exclude_species_groups %||% "Ferns and allies"))
  kingdom <- taxon_settings(config)$ala_kingdom
  # Cache key: the genus alone for Plantae (keeps the existing cache valid),
  # genus + kingdom otherwise, so results for another kingdom are not mixed in.
  genus_key <- function(g) if (identical(kingdom, "Plantae")) digest::digest(g) else digest::digest(paste(g, kingdom, sep = "\r"))

  results <- data.table::rbindlist(lapply(genera, function(g) {
    cache_file <- file.path(cache_dir, paste0(genus_key(g), ".rds"))
    if (file.exists(cache_file)) return(readRDS(cache_file))
    res <- ala_genus_lookup(g, base_url, kingdom)
    if (identical(res$ala_status, "queried")) saveRDS(res, cache_file)  # retry errors next run
    Sys.sleep(interval)
    res
  }), fill = TRUE)

  if (nrow(results) == 0L) {
    results <- data.table::data.table(genus = character(), ala_order = character(),
                                      ala_species_group = character(), ala_status = character())
  }
  results[, ala_fern_lycophyte := (!is.na(ala_order) & ala_order %in% fern_orders) |
            (!is.na(ala_species_group) & Reduce(`|`, lapply(groups, function(x) grepl(x, ala_species_group, fixed = TRUE)), FALSE))]
  static <- read_excluded_genera(config)
  results[, in_static_list := genus %in% static]
  results[, exclude := ala_fern_lycophyte | in_static_list]
  results[, exclude_reason := data.table::fcase(
    ala_fern_lycophyte & in_static_list, "ALA + static list",
    ala_fern_lycophyte, "ALA",
    in_static_list, "static list",
    default = ""
  )]
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(results, output_path)
  n_err <- sum(results$ala_status != "queried", na.rm = TRUE)
  if (n_err > 0L) warning(n_err, " ALA genus lookups failed; the static genus list still applies to them. Rerun to retry.")
  output_path
}

# Genera to exclude: anything ALA flags as fern/lycophyte, plus the static list
# (which covers ALA homonyms such as Alsophila, Polypodium and Pteridium).
excluded_genera_from_checks <- function(check_file, config) {
  checks <- data.table::fread(check_file)
  sort(unique(c(checks[exclude == TRUE, genus], read_excluded_genera(config))))
}
