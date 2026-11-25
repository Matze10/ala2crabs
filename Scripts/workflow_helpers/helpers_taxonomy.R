# Human-in-the-loop taxonomy review helpers.
# Automated flags generate proposals; only explicit manual decisions can
# promote a name into the reviewed species list.

# ALA phrase names such as Grevillea sp. "Baldy Range" contain quote marks.
# Every CSV round trip (fwrite/fread, Excel) doubles them, so the same name
# can appear with 16 or 64 quote characters. Match decisions to queue names on
# a key with repeated quotes collapsed, then use the queue's spelling.
name_match_key <- function(x) gsub('"+', '"', trimws(as.character(x)))

# Replace each decision's normalised_name by the queue's spelling when they
# match on name_match_key(); unmatched names are left unchanged.
align_decision_names <- function(decisions, queue_names) {
  if (!nrow(decisions)) return(decisions)
  key_queue <- name_match_key(queue_names)
  idx <- match(name_match_key(decisions$normalised_name), key_queue)
  hit <- !is.na(idx)
  decisions[hit, normalised_name := queue_names[idx[hit]]]
  decisions
}

# Read the version-controlled manual decision table.
read_manual_decisions <- function(config) {
  path <- here::here(config$taxonomy$manual_decisions)
  required <- c(
    "normalised_name", "accepted_name", "accepted_rank", "decision",
    "reason", "authority", "authority_version"
  )
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.table::as.data.table(setNames(
      lapply(required, character), required
    )))
  }
  decisions <- data.table::fread(path, na.strings = c("", "NA"))
  missing <- setdiff(required, names(decisions))
  if (length(missing) > 0) {
    stop("Manual decision file is missing: ", paste(missing, collapse = ", "))
  }
  decisions[, (required) := lapply(.SD, as.character), .SDcols = required]
  # Template rows left blank in the spreadsheet are not decisions.
  decisions <- decisions[!is.na(normalised_name) & nzchar(trimws(normalised_name))]
  decisions[decision == "ok", decision := "include"]
  decisions[, ..required]
}

# Generate one review row per normalized candidate name.
# Occurrence-level detail remains in name_candidates_normalised.csv; this queue
# contains only the evidence needed for a taxonomic decision.
build_review_queue <- function(normalised_file, output_path) {
  candidates <- data.table::fread(normalised_file)
  queue <- candidates[
    , .(
      original_names = paste(sort(unique(normalised_name)), collapse = " | "),
      candidate_sources = paste(sort(unique(unlist(strsplit(candidate_sources, ",", fixed = TRUE)))), collapse = ", "),
      original_taxon_ranks = "from-filtered-SQL",
      genera = "",
      occurrence_records = sum(source_observations, na.rm = TRUE),
      source_observations = sum(source_observations, na.rm = TRUE),
      has_binomial = any(looks_like_binomial, na.rm = TRUE),
      has_infraspecific_marker = any(contains_infraspecific_marker, na.rm = TRUE),
      has_uncertainty_marker = any(contains_uncertainty_marker, na.rm = TRUE),
      has_hybrid_marker = any(contains_hybrid_marker, na.rm = TRUE),
      has_authority_like_text = any(contains_authority_like_text, na.rm = TRUE)
    ),
    by = .(normalised_name)
  ]
  # Canonical form of each queue name (authority stripped, rank markers
  # standardised; see helpers_canonical_names.R). A name that parses to a clean
  # species or infraspecific canonical name is included automatically, so
  # author variants and marked infraspecific names no longer need a manual
  # decision. Only uncertain (sp./cf./aff./?/genus-only), hybrid and unparsed
  # names go to the manual review queue.
  parsed <- canonicalise_names(queue$normalised_name)
  clean_ranks <- c("species", "subspecies", "variety", "form", "nothosubspecies", "nothovariety", "infraspecific_unmarked")
  queue[, `:=`(
    canonical_name = parsed$canonical_name,
    canonical_binomial = parsed$canonical_binomial,
    canonical_rank = parsed$canonical_rank
  )]
  queue[, auto_ok := canonical_rank %chin% clean_ranks & !is.na(canonical_binomial) &
          !has_uncertainty_marker & !has_hybrid_marker]
  queue[, `:=`(
    proposed_accepted_name = data.table::fifelse(is.na(canonical_name), normalised_name, canonical_name),
    proposed_rank = canonical_rank,
    proposed_status = data.table::fifelse(auto_ok, "candidate", "needs_review"),
    decision = data.table::fifelse(auto_ok, "include", NA_character_),
    accepted_name = data.table::fifelse(auto_ok, canonical_name, NA_character_),
    accepted_rank = data.table::fifelse(
      auto_ok,
      data.table::fifelse(canonical_rank == "infraspecific_unmarked", "infraspecific", canonical_rank),
      NA_character_
    ),
    reason = NA_character_,
    authority = NA_character_,
    authority_version = NA_character_
  )]
  data.table::setcolorder(queue, c(
    "normalised_name", "proposed_accepted_name", "proposed_rank", "proposed_status",
    "original_names", "candidate_sources", "original_taxon_ranks", "genera",
    "occurrence_records", "source_observations", "has_binomial",
    "has_infraspecific_marker", "has_uncertainty_marker", "has_hybrid_marker",
    "has_authority_like_text", "decision", "accepted_name", "accepted_rank",
    "reason", "authority", "authority_version",
    "canonical_name", "canonical_binomial", "canonical_rank"
  ))
  queue[, auto_ok := NULL]
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(queue, output_path)

  review_only_path <- sub("taxonomy_review_queue\\.csv$", "manual_review_queue.csv", output_path)
  data.table::fwrite(queue[decision != "include" | is.na(decision)], review_only_path)
  output_path
}

# Join validated manual decisions back to the generated review queue.
apply_manual_decisions <- function(review_file, config, output_path, synonym_file = NULL) {
  review <- data.table::fread(review_file)
  decisions <- read_manual_decisions(config)
  required <- c("normalised_name", "accepted_name", "accepted_rank", "decision", "reason", "authority", "authority_version")
  missing <- setdiff(required, names(decisions))
  if (length(missing) > 0) stop("Manual decision file is missing: ", paste(missing, collapse = ", "))
  decisions <- decisions[, ..required]
  decisions <- align_decision_names(decisions, review$normalised_name)
  if (anyDuplicated(decisions$normalised_name)) stop("Manual decisions must have one row per normalised name.")
  matches <- merge(review, decisions, by = "normalised_name", all.x = TRUE, suffixes = c("_proposal", ""))
  # Manual decisions are only made for doubtful names. Every other name keeps
  # the automatic proposal from build_review_queue() (clean binomials are
  # proposed as include); a manual row always overrides the proposal.
  matches[, decision_source := data.table::fifelse(!is.na(decision) & nzchar(decision), "manual", "automatic")]
  fill_cols <- c("decision", "accepted_name", "accepted_rank", "reason", "authority", "authority_version")
  for (cc in fill_cols) {
    prop <- paste0(cc, "_proposal")
    if (prop %in% names(matches)) {
      matches[decision_source == "automatic", (cc) := as.character(get(prop))]
    }
  }
  matches[decision_source == "automatic" & decision == "include" & is.na(reason), reason := "automatic: clean binomial"]
  matches[is.na(decision) | decision == "", decision := "needs_review"]
  # APC/AFD accepted name and NCBI name (helpers_synonyms.R): automatic includes
  # take the NCBI name; a blank or placeholder authority is filled from the sources.
  if (!is.null(synonym_file)) matches <- add_synonym_matches(matches, synonym_file)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(matches, output_path)
  output_path
}
