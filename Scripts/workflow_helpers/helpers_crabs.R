# CRABS export helpers.
# The final species file is plain text with one explicitly accepted name per line.

# Export names with explicit include decisions for Entrez checking.
export_reviewed_species <- function(matches_file, output_path, config = NULL, expansion_file = NULL) {
  matches <- data.table::fread(matches_file)
  inc <- include_decision_values(config)  # crabs.include_decisions (default include)
  species <- matches[
    decision %chin% inc &
      !is.na(accepted_name) & accepted_name != "" &
      !is.na(accepted_rank) & accepted_rank %chin% c("species", "subspecies", "variety", "form", "infraspecific"),
    unique(accepted_name)
  ]
  # Strip authorities and standardise rank markers so Entrez sees one name per taxon.
  parsed <- canonicalise_names(as.character(species))
  species <- parsed[canonical_rank != "genus_or_uncertain", unique(canonical_name)]
  species <- as.character(species)
  # Included names whose NCBI name is a phrase name (written with
  # crabs.keep_ncbi_phrase_names) need marker counts under that NCBI name.
  if (isTRUE(config$crabs$keep_ncbi_phrase_names) && "ncbi_name" %in% names(matches)) {
    species <- c(species, matches[decision %chin% inc & !is.na(ncbi_name) & grepl(" sp[.] ", ncbi_name), unique(ncbi_name)])
  }
  # Genus expansion and NCBI phrase names (helpers_expansion.R) get marker
  # counts too; their names are used as they are (NCBI or authority spelling).
  if (!is.null(expansion_file) && file.exists(expansion_file[1])) {
    exp <- data.table::fread(expansion_file[1], na.strings = c("", "NA"))
    if (nrow(exp)) species <- c(species, as.character(exp$crabs_name))
  }
  species <- unique(species[!is.na(species) & nzchar(species)])
  species <- sort(species)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(as.character(species), output_path, useBytes = TRUE)
  output_path
}

# Rows for build_crabs_outputs() from the genus expansion table, with the same
# columns and flags as the included decisions. NCBI phrase names are written
# exactly as NCBI spells them (they are the only non-binomial lines besides
# hybrids), so they bypass the binomial parsing.
expansion_include_rows <- function(expansion_file) {
  if (is.null(expansion_file) || !file.exists(expansion_file[1])) return(NULL)
  exp <- data.table::fread(expansion_file[1], na.strings = c("", "NA"))
  if (!nrow(exp)) return(NULL)
  parsed <- canonicalise_names(exp$crabs_name)
  phrase <- exp$source == "ncbi_phrase_name"
  data.table::data.table(
    accepted_name = exp$crabs_name,
    accepted_name_input = exp$crabs_name,
    accepted_binomial = ifelse(phrase, exp$crabs_name, parsed$canonical_binomial),
    accepted_rank = ifelse(phrase, "phrase_name", parsed$canonical_rank),
    has_binomial = phrase | !is.na(parsed$canonical_binomial),
    has_infraspecific_marker = !phrase & parsed$canonical_rank %chin% c("subspecies", "variety", "form", "nothosubspecies", "nothovariety", "infraspecific_unmarked"),
    has_hybrid_marker = !phrase & parsed$canonical_rank == "hybrid",
    has_uncertainty_marker = !phrase & parsed$canonical_rank == "genus_or_uncertain",
    decision = "include",
    decision_source = exp$source,
    ncbi_name = exp$crabs_name
  )
}

# Build the final CRABS list and companion exclusion/review audit files.
build_crabs_outputs <- function(matches_file, entrez_file, config, species_path, exclusions_path, review_path,
                                latitude_flags_file = NULL, region_path = NULL, expansion_file = NULL) {
  flags <- read_latitude_flags(latitude_flags_file)
  matches <- data.table::fread(matches_file)
  inc <- include_decision_values(config)  # crabs.include_decisions (default include)
  include <- matches[
    decision %chin% inc &
      !is.na(accepted_name) & accepted_name != ""
  ]
  # Canonicalise the reviewer's accepted name: authorities removed, rank taken
  # from the marker in the name (fixes var. names labelled as subspecies).
  parsed <- canonicalise_names(include$accepted_name)
  include[, `:=`(
    accepted_name_input = accepted_name,
    accepted_name = parsed$canonical_name,
    accepted_binomial = parsed$canonical_binomial,
    accepted_rank = parsed$canonical_rank,
    has_binomial = !is.na(parsed$canonical_binomial),
    has_infraspecific_marker = parsed$canonical_rank %chin% c("subspecies", "variety", "form", "nothosubspecies", "nothovariety", "infraspecific_unmarked"),
    has_hybrid_marker = parsed$canonical_rank == "hybrid",
    has_uncertainty_marker = parsed$canonical_rank == "genus_or_uncertain"
  )]
  # Species from the genus expansion and NCBI phrase names (helpers_expansion.R)
  # join here and pass the same filters below (hybrids, marker hit, region
  # exclusion; they have no ALA records, so no region flags).
  exp_rows <- expansion_include_rows(expansion_file)
  if (!is.null(exp_rows)) include <- data.table::rbindlist(list(include, exp_rows), use.names = TRUE, fill = TRUE)
  # Named hybrids that NCBI lists as species (e.g. Citrus x aurantiifolia,
  # Tilia x europaea) are kept under their NCBI spelling when
  # crabs.keep_ncbi_hybrids is true; hybrid formulas and unmatched hybrids stay out.
  hybrid_pattern <- "^[A-Z][a-z-]+ [x×] [a-z-]+$|^[x×] [A-Z][a-z-]+ [a-z-]+$"
  keep_hybrids <- isTRUE(config$crabs$keep_ncbi_hybrids) && "ncbi_name" %in% names(include)
  if (keep_hybrids) {
    include[, ncbi_hybrid := !is.na(ncbi_name) & grepl(hybrid_pattern, ncbi_name)]
  } else {
    include[, ncbi_hybrid := FALSE]
  }
  include[ncbi_hybrid == TRUE, `:=`(accepted_name = ncbi_name, accepted_binomial = ncbi_name,
                                     has_binomial = TRUE, has_uncertainty_marker = FALSE)]
  # Included names whose NCBI name is a phrase name (e.g. an ALA name that the
  # synonym step matched to "Olearia sp. Rhizomatica (I.R.Telford 11549)") are
  # written under that NCBI name when crabs.keep_ncbi_phrase_names is true.
  if (isTRUE(config$crabs$keep_ncbi_phrase_names) && "ncbi_name" %in% names(include)) {
    include[, ncbi_phrase := !is.na(ncbi_name) & grepl(" sp[.] ", ncbi_name) & !ncbi_hybrid]
    include[ncbi_phrase == TRUE, `:=`(accepted_name = ncbi_name, accepted_binomial = ncbi_name, has_binomial = TRUE,
                                      has_uncertainty_marker = FALSE, has_infraspecific_marker = FALSE)]
  }
  entrez <- data.table::fread(entrez_file)
  # Older cache files / outputs used the column name ncbi_its_hits.
  if ("ncbi_its_hits" %in% names(entrez) && !"ncbi_marker_hits" %in% names(entrez)) {
    data.table::setnames(entrez, "ncbi_its_hits", "ncbi_marker_hits")
  }
  # crabs.require_marker_hit (old name: require_its_hit): keep only species with
  # at least one NCBI sequence for the configured marker.
  require_hit <- isTRUE(config$crabs$require_marker_hit) || isTRUE(config$crabs$require_its_hit)
  if (require_hit && nrow(entrez) > 0) {
    with_hits <- entrez[!is.na(ncbi_marker_hits) & ncbi_marker_hits > 0, unique(species_name)]
    include <- include[accepted_name %chin% with_hits]
  }
  if (isTRUE(config$crabs$exclude_uncertain)) include <- include[!include$has_uncertainty_marker]
  if (isTRUE(config$crabs$exclude_hybrids)) include <- include[!include$has_hybrid_marker | include$ncbi_hybrid]
  if (isTRUE(config$crabs$exclude_genus_only)) include <- include[include$has_binomial]
  if (isTRUE(config$crabs$exclude_infraspecific)) include <- include[!include$has_infraspecific_marker]
  # Region-flagged species (geography.region_flags) stay in the list unless
  # their region has exclude: true.
  regions <- region_flag_settings(config)
  if (!is.null(flags)) {
    for (reg in regions) {
      if (reg$exclude && reg$name %in% names(flags)) {
        dropped <- flags[get(reg$name) %in% TRUE, species_binomial]
        include <- include[!(accepted_binomial %chin% dropped)]
      }
    }
  }
  collapse <- !isFALSE(config$crabs$collapse_to_binomial)
  species <- as.character(unique(if (collapse) include$accepted_binomial else include$accepted_name))
  species <- species[!is.na(species) & nzchar(species)]
  species <- sort(species)
  if (nrow(entrez) > 0) {
    species_audit <- data.table::data.table(species_name = species)[entrez, on = "species_name"]
  } else {
    species_audit <- data.table::data.table(species_name = species)
  }
  exclusions <- matches[decision %chin% c("exclude", "insufficient_evidence")]
  review <- matches[is.na(decision) | decision == "needs_review"]
  if (!is.null(flags)) review <- add_region_flags(review, flags, "normalised_name")
  dir.create(dirname(species_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(exclusions_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(review_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(species, species_path, useBytes = TRUE)
  data.table::fwrite(exclusions, exclusions_path)
  data.table::fwrite(data.table::merge.data.table(review, species_audit, by.x = "accepted_name", by.y = "species_name", all.x = TRUE), review_path)
  out <- c(species_path, exclusions_path, review_path)
  # Audit list of every species flagged by any region, and whether it is in the CRABS list.
  if (!is.null(flags) && !is.null(region_path)) {
    flagged <- flags[nzchar(region_flags) & !is.na(region_flags)]
    region_cols <- region_columns(flags)
    count_cols <- c("n_records", "n_georeferenced", paste0("n_inside_", region_cols), paste0("n_outside_", region_cols),
                    "southernmost_latitude", "northernmost_latitude", "westernmost_longitude", "easternmost_longitude")
    flagged <- flagged[, intersect(c("species_binomial", "region_flags", region_cols, count_cols), names(flagged)), with = FALSE]
    flagged[, in_crabs_list := species_binomial %chin% (if (collapse) species else canonicalise_names(species)$canonical_binomial)]
    dir.create(dirname(region_path), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(flagged, region_path)
    out <- c(out, region_path)
  }
  out
}
