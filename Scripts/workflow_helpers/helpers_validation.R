# Validate reviewer decisions before they affect downstream outputs.
# The table must contain one row per queue name and use an allowed decision.
validate_manual_decisions <- function(review_file, decision_file, output_file, latitude_flags_file = NULL) {
  review <- data.table::fread(review_file)
  required <- c("normalised_name", "accepted_name", "accepted_rank", "decision", "reason", "authority", "authority_version")
  if (!file.exists(decision_file) || file.info(decision_file)$size == 0) {
    decisions <- data.table::as.data.table(setNames(lapply(required, character), required))
  } else {
    decisions <- data.table::fread(decision_file, na.strings = c("", "NA"))
    missing <- setdiff(required, names(decisions))
    if (length(missing) > 0) stop("Manual decision file is missing: ", paste(missing, collapse = ", "))
    decisions[, (required) := lapply(.SD, as.character), .SDcols = required]
    # Template rows left blank in the spreadsheet are not decisions.
    decisions <- decisions[!is.na(normalised_name) & nzchar(trimws(normalised_name))]
  }
  missing <- setdiff(required, names(decisions))
  if (length(missing) > 0) stop("Manual decision file is missing: ", paste(missing, collapse = ", "))
  # Treat the reviewer shorthand `ok` as an explicit inclusion.
  decisions[decision == "ok", decision := "include"]
  allowed <- c("include", "exclude", "needs_review", "insufficient_evidence")
  invalid <- setdiff(unique(decisions$decision), allowed)
  invalid <- invalid[!is.na(invalid)]
  if (length(invalid) > 0) stop("Invalid manual decisions: ", paste(invalid, collapse = ", "))
  decisions <- align_decision_names(decisions, review$normalised_name)  # helpers_taxonomy.R
  if (anyDuplicated(decisions$normalised_name)) stop("Manual decisions must have one row per normalised name.")
  unknown <- setdiff(decisions$normalised_name, review$normalised_name)
  if (length(unknown) > 0) stop("Manual decisions contain names absent from the review queue: ",
                                paste(head(unknown, 10), collapse = "; "))
  # Warn when the stated rank disagrees with the name, or an include has no species epithet.
  inc <- decisions[decision == "include" & !is.na(accepted_name)]
  if (nrow(inc) > 0) {
    parsed <- canonicalise_names(inc$accepted_name)
    rank_from_name <- ifelse(parsed$canonical_rank == "infraspecific_unmarked", NA_character_, parsed$canonical_rank)
    bad_rank <- !is.na(rank_from_name) & !is.na(inc$accepted_rank) & rank_from_name != "genus_or_uncertain" &
      inc$accepted_rank != rank_from_name
    if (any(bad_rank)) warning("accepted_rank disagrees with the name's rank marker for: ",
      paste0(inc$accepted_name[bad_rank], " (", inc$accepted_rank[bad_rank], " -> ", rank_from_name[bad_rank], ")", collapse = "; "))
    no_species <- parsed$canonical_rank == "genus_or_uncertain"
    if (any(no_species)) message(sum(no_species), " included names are genus-level (",
      length(unique(inc$accepted_name[no_species])), " genera). They are not written to the CRABS list as genera;",
      " the genus_expansion step adds their species instead (see output.genus_expansion_audit).")
  }
  # Review sheet: EVERY review candidate (names build_review_queue() did not
  # include automatically) plus any other name with a manual decision. Existing
  # decisions from taxonomy_decisions.csv are pre-filled so the reviewer can
  # edit them here and copy the decision columns back into
  # taxonomy_decisions.csv. review_status says which rows are still open.
  # (Earlier versions listed only names WITHOUT a decision, so the sheet went
  # empty once every candidate had one.)
  candidates <- review[(is.na(decision) | decision != "include") | normalised_name %in% decisions$normalised_name]
  decision_cols <- c("decision", "accepted_name", "accepted_rank", "reason", "authority", "authority_version",
                     "reviewer", "review_date")
  for (cc in setdiff(decision_cols, names(decisions))) decisions[, (cc) := NA_character_]
  manual <- decisions[, c("normalised_name", decision_cols), with = FALSE]
  manual[, (decision_cols) := lapply(.SD, as.character), .SDcols = decision_cols]
  proposal_cols <- intersect(decision_cols, names(candidates))
  data.table::setnames(candidates, proposal_cols, paste0("auto_", proposal_cols))
  # Decisions live only in taxonomy_decisions.csv. The sheet leaves the
  # decision columns blank and only says whether a row already exists there.
  unresolved <- data.table::copy(candidates)
  unresolved[, in_taxonomy_decisions := normalised_name %chin% manual$normalised_name]
  unresolved[, (decision_cols) := NA_character_]
  # Why the name needs review, in words instead of five TRUE/FALSE columns.
  why <- function(flag, label) data.table::fifelse(flag %in% TRUE, label, NA_character_)
  unresolved[, why_review := apply(cbind(
    why(has_uncertainty_marker, "uncertain (sp./cf./?/genus only)"),
    why(has_hybrid_marker, "hybrid"),
    why(!has_binomial %in% TRUE, "no binomial"),
    why(has_infraspecific_marker, "infraspecific"),
    why(has_authority_like_text, "authority or extra text")
  ), 1, function(x) paste(x[!is.na(x)], collapse = "; "))]
  unresolved[!nzchar(why_review), why_review := "clean name, listed because it has a manual decision"]
  # Repeated quote marks from CSV round trips are collapsed for display;
  # decisions still match (align_decision_names()).
  unresolved[, normalised_name := name_match_key(normalised_name)]
  unresolved[, suggested_name := canonical_name]
  keep <- c(
    # first nine columns = taxonomy_decisions.csv, copy across as a block
    "normalised_name", "accepted_name", "accepted_rank", "decision", "reason",
    "authority", "authority_version", "reviewer", "review_date",
    # evidence for the decision
    "in_taxonomy_decisions", "why_review", "suggested_name", "canonical_rank", "occurrence_records"
  )
  unresolved <- unresolved[, intersect(keep, names(unresolved)), with = FALSE]
  data.table::setorder(unresolved, in_taxonomy_decisions, normalised_name)  # undecided rows first
  # Region flags (helpers_geography.R), appended after the evidence columns:
  # informational, nothing is removed here.
  flags <- read_latitude_flags(latitude_flags_file)
  if (!is.null(flags)) unresolved <- add_region_flags(unresolved, flags, "normalised_name")
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  # This file is regenerated on every run. If someone has typed decisions into
  # it, keep a timestamped copy before overwriting so no edits are lost.
  if (file.exists(output_file)) {
    old <- tryCatch(data.table::fread(output_file, na.strings = c("", "NA")), error = function(e) NULL)
    if (!is.null(old) && "decision" %in% names(old) && any(!is.na(old$decision) & nzchar(as.character(old$decision)))) {
      backup_dir <- here::here("Scripts/data/manual/review_queue_backups")
      dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
      backup <- file.path(backup_dir, paste0("review_queue_", format(file.mtime(output_file), "%Y-%m-%d_%H%M%S"), ".csv"))
      file.copy(output_file, backup, overwrite = FALSE, copy.date = TRUE)
      message("review_queue.csv contained decisions; copy kept at ", backup)
    }
  }
  # BOM so Excel opens the file as UTF-8 (diacritics, x signs in hybrids).
  data.table::fwrite(unresolved, output_file, bom = TRUE)
  list(decisions = decisions, unresolved = unresolved, path = output_file)
}
