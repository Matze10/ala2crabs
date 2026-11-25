# Candidate-name extraction and normalisation helpers.
# Candidate names are aggregated directly in SQLite to avoid loading all
# NSW/QLD occurrence rows into R memory.

# ALA name fields used for candidate names (target.candidate_fields; default
# raw_scientificName, scientificName, acceptedNameUsage, originalNameUsage).
# Stops if a listed field is not a column of the imported table.
candidate_fields <- function(con, config) {
  fields <- as.character(unlist(config$target$candidate_fields))
  if (!length(fields)) fields <- c("raw_scientificName", "scientificName", "acceptedNameUsage", "originalNameUsage")
  cols <- DBI::dbListFields(con, "ala_occurrences_raw")
  bad <- setdiff(fields, cols)
  if (length(bad)) stop("target.candidate_fields: not a column of the ALA records: ", paste(bad, collapse = ", "))
  fields
}

# One row per distinct name string among the target records. Names are taken
# from the ALA name fields in target.candidate_fields plus genus + epithet(s)
# rebuilt from the rank fields (always on), with the number of records that
# carry each name. Uses the same record filter as filter_target_records(), and
# drops names whose first word is an excluded genus.
extract_name_candidates <- function(database, config, output_path, excluded_genera = read_excluded_genera(config)) {
  con <- DBI::dbConnect(RSQLite::SQLite(), database)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  rank_sql <- paste(DBI::dbQuoteString(con, as.character(unlist(config$target$accepted_ranks))), collapse = ", ")
  phylum_sql <- paste(DBI::dbQuoteString(con, as.character(unlist(config$target$exclude_phyla))), collapse = ", ")
  class_sql <- paste(DBI::dbQuoteString(con, as.character(unlist(config$target$exclude_classes))), collapse = ", ")

  # Name fields from target.candidate_fields (checked against the table).
  fields <- candidate_fields(con, config)
  field_sql <- vapply(seq_along(fields), function(i) {
    col <- DBI::dbQuoteIdentifier(con, fields[i]); lab <- DBI::dbQuoteString(con, fields[i])
    if (i == 1L) paste0("SELECT rowid AS source_record_rowid, ", lab, " AS candidate_source, ", col, " AS original_name FROM ala_occurrences_raw ")
    else paste0("UNION ALL SELECT rowid, ", lab, ", ", col, " FROM ala_occurrences_raw ")
  }, character(1))

  sql <- paste0(
    "WITH name_observations AS (",
    paste(field_sql, collapse = ""),
    "UNION ALL SELECT rowid, 'reconstructed_binomial', TRIM(genus || ' ' || specificEpithet) FROM ala_occurrences_raw ",
    "WHERE genus IS NOT NULL AND TRIM(genus) <> '' AND specificEpithet IS NOT NULL AND TRIM(specificEpithet) <> '' ",
    "UNION ALL SELECT rowid, 'reconstructed_infraspecific', TRIM(genus || ' ' || specificEpithet || ' ' || ",
    "CASE LOWER(taxonRank) WHEN 'subspecies' THEN 'subsp. ' WHEN 'variety' THEN 'var. ' WHEN 'form' THEN 'f. ' ",
    "WHEN 'nothosubspecies' THEN 'nothosubsp. ' WHEN 'nothovariety' THEN 'nothovar. ' ELSE '' END || ",
    "infraspecificEpithet) FROM ala_occurrences_raw ",
    "WHERE genus IS NOT NULL AND TRIM(genus) <> '' AND specificEpithet IS NOT NULL AND TRIM(specificEpithet) <> '' AND infraspecificEpithet IS NOT NULL AND TRIM(infraspecificEpithet) <> ''",
    "), filtered_observations AS (",
    "SELECT o.* FROM name_observations o WHERE o.original_name IS NOT NULL AND TRIM(o.original_name) <> ''",
    "AND EXISTS (SELECT 1 FROM ala_occurrences_raw r WHERE r.rowid = o.source_record_rowid",
    " AND r.taxonRank IN (", rank_sql, ")",
    kingdom_sql(con, config, prefix = "r."),  # taxon.ala_kingdom (default Plantae)
    " AND (r.phylum IS NULL OR TRIM(r.phylum) = '' OR r.phylum NOT IN (", phylum_sql, "))",
    " AND (r.class IS NULL OR TRIM(r.class) = '' OR r.class NOT IN (", class_sql, "))",
    seed_plant_sql(con, config$target$exclude_orders, excluded_genera, prefix = "r."),
    ")) SELECT TRIM(original_name) AS normalised_name,",
    "GROUP_CONCAT(DISTINCT candidate_source) AS candidate_sources,",
    "COUNT(DISTINCT source_record_rowid) AS source_observations",
    " FROM filtered_observations GROUP BY TRIM(original_name)"
  )

  candidates <- data.table::setDT(DBI::dbGetQuery(con, sql))
  candidates <- candidates[!is.na(normalised_name) & normalised_name != ""]
  # Second pass on the name itself, for records whose genus field is blank or
  # disagrees with the name (e.g. raw 'Lindsaya' vs matched 'Lindsaea').
  first_word <- sub("^[!\\s]*(\\S+).*$", "\\1", candidates$normalised_name, perl = TRUE)
  candidates <- candidates[!(first_word %in% excluded_genera)]
  candidates[, original_taxon_rank := NA_character_]
  candidates[, original_genus := NA_character_]
  candidates[, original_name := normalised_name]

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(candidates, output_path)
  output_path
}

# Add the canonical name, binomial and rank (helpers_canonical_names.R) and the
# review flags used by build_review_queue(): uncertain, hybrid, infraspecific,
# authority text. Drops strings that do not start like "Genus epithet".
normalise_name_candidates <- function(candidate_file, output_path) {
  candidates <- data.table::fread(candidate_file)
  # has_leading_exclamation must be read before the "!" is stripped.
  candidates[, has_leading_exclamation := grepl("^\\s*!+", normalised_name)]
  candidates[, normalised_name := trimws(sub("^\\s*!+\\s*", "", normalised_name))]
  # Flags come from the canonical-name parser (helpers_canonical_names.R);
  # the earlier \\b-after-dot regexes missed most rank markers and authorities.
  parsed <- canonicalise_names(candidates$normalised_name)
  candidates[, `:=`(
    canonical_name = parsed$canonical_name,
    canonical_binomial = parsed$canonical_binomial,
    canonical_rank = parsed$canonical_rank,
    # Kept permissive so genus-only and phrase names stay reviewable (flagged uncertain).
    looks_like_binomial = grepl("^[A-Z][[:alpha:]-]+\\s+[a-z][[:alpha:]-]+", normalised_name),
    contains_infraspecific_marker = parsed$canonical_rank %chin% c("subspecies", "variety", "form", "nothosubspecies", "nothovariety", "infraspecific_unmarked"),
    contains_uncertainty_marker = parsed$canonical_rank == "genus_or_uncertain" |
      grepl("(^|\\s)(sp|spp|cf|aff|indet)\\.?(\\s|$)|\\?", normalised_name, ignore.case = TRUE),
    contains_hybrid_marker = parsed$canonical_rank == "hybrid",
    contains_authority_like_text = parsed$authority_removed
  )]
  candidates <- candidates[
    !is.na(normalised_name) & normalised_name != "" &
      looks_like_binomial
  ]
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(candidates, output_path)
  output_path
}
