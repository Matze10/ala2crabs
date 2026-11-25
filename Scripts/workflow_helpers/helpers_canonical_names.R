# Canonical scientific-name parsing for the ALA -> CRABS workflow.
# Strips authorities, standardises infraspecific rank markers, and returns
# both the full canonical name and the binomial used for the CRABS list.
# Base R + data.table only.

.rank_markers <- c(
  "subsp." = "subsp.", "subsp" = "subsp.", "ssp." = "subsp.", "ssp" = "subsp.",
  "var." = "var.", "var" = "var.", "f." = "f.", "fo." = "f.", "forma" = "f.",
  "nothosubsp." = "nothosubsp.", "nothovar." = "nothovar.", "nvar." = "nothovar."
)
.marker_rank <- c(
  "subsp." = "subspecies", "var." = "variety", "f." = "form",
  "nothosubsp." = "nothosubspecies", "nothovar." = "nothovariety"
)
# ALA taxonRank -> marker, used to repair unmarked reconstructed trinomials.
.ala_rank_marker <- c(
  subspecies = "subsp.", variety = "var.", form = "f.",
  nothosubspecies = "nothosubsp.", nothovariety = "nothovar."
)
.diacritics <- c("\u00eb" = "e", "\u00ef" = "i", "\u00f6" = "o", "\u00fc" = "u",
                 "\u00e4" = "a", "\u00e9" = "e", "\u00e8" = "e")
.name_particles <- c("de", "van", "von", "der", "den", "du", "la", "le", "da", "di", "ex", "et", "in", "&")
.uncertain_epithets <- c("sp", "sp.", "spp", "spp.", "cf", "cf.", "aff", "aff.", "indet", "indet.", "complex", "hybrid")

# TRUE for words that look like a species or infraspecific epithet: lowercase
# letters/hyphens, not an author particle (de, van, ex ...) or rank marker.
.is_epithet <- function(w) grepl("^[a-z][a-z-]+$", w) & !(w %in% .name_particles) & !(tolower(w) %in% names(.rank_markers))

# Parse one name. rank_hint is the ALA taxonRank for the record, if known.
.parse_one <- function(s, rank_hint = NA_character_) {
  out <- list(canonical_name = NA_character_, canonical_binomial = NA_character_,
              canonical_rank = "unparsed", infraspecific_epithet = NA_character_,
              authority_removed = FALSE)
  s <- sub("^\\s*!+\\s*", "", s)
  # Non-breaking and other Unicode spaces (common in ALA raw names, e.g.
  # "Backhousia myrtifolia Hook. & Harv.") would otherwise glue the
  # epithet to the authority and make the name unparseable.
  s <- gsub("[  -​  　]", " ", s, perl = TRUE)
  # Common name first, scientific name after a dash
  # ("Cape weed. - Arctotheca calendula", "Bidgee- widgee. -Acaena anserinifolia").
  # The common name ends in a full stop before the dash; "A b - C d"
  # intergrade names (no full stop) are left alone.
  s <- sub("^[^-]*\\.\\s+-\\s*(?=[A-Z][a-z-]+\\s+[a-z])", "", s, perl = TRUE)
  # Trailing introduced-species asterisk, full stop or stray bracket on the
  # epithet ("Lantana camara*", "Amyema cambagei.", "Neptunia monosperma)").
  s <- sub("^(\\S+\\s+[a-z][a-z-]+)[*.)]+(\\s|$)", "\\1\\2", s, perl = TRUE)
  s <- gsub("×", " × ", s)
  for (k in seq_along(.diacritics)) s <- gsub(names(.diacritics)[k], .diacritics[[k]], s, fixed = TRUE)  # e.g. Hierochloe
  w <- strsplit(trimws(gsub("\\s+", " ", s)), " ", fixed = TRUE)[[1]]
  if (!length(w) || !grepl("^[A-Z][a-z-]+$", w[1])) return(out)
  genus <- w[1]; i <- 2L; hybrid <- FALSE
  if (i <= length(w) && w[i] %in% c("×", "x", "X")) { hybrid <- TRUE; i <- i + 1L }
  if (i > length(w) || !.is_epithet(w[i]) || w[i] %in% .uncertain_epithets) {
    out$canonical_rank <- "genus_or_uncertain"; return(out)
  }
  sp <- w[i]; rest <- w[-seq_len(i)]
  binomial <- if (hybrid) paste(genus, "×", sp) else paste(genus, sp)
  out$canonical_binomial <- binomial
  out$canonical_name <- binomial
  out$canonical_rank <- if (hybrid) "hybrid" else "species"
  if (!length(rest)) return(out)

  # 1. Explicit rank marker anywhere after the epithet (authority may sit between).
  mk <- which(tolower(rest) %in% names(.rank_markers))
  mk <- mk[mk < length(rest) & .is_epithet(rest[pmin(mk + 1L, length(rest))])]
  if (length(mk)) {
    j <- mk[1]
    marker <- unname(.rank_markers[tolower(rest[j])]); ie <- rest[j + 1L]
    out$authority_removed <- length(rest) > 2L
    if (ie == sp) return(out)                       # autonym -> treat as the species
    out$canonical_name <- paste(binomial, marker, ie)
    out$canonical_rank <- unname(.marker_rank[marker])
    out$infraspecific_epithet <- ie
    return(out)
  }
  # 2. Unmarked trinomial (ALA reconstructed genus + species + infraspecificEpithet).
  if (.is_epithet(rest[1])) {
    ie <- rest[1]
    if (ie == sp) return(out)                       # autonym
    if (ie == "hybrid") { out$canonical_rank <- "hybrid"; return(out) }
    out$infraspecific_epithet <- ie
    marker <- if (!is.na(rank_hint)) unname(.ala_rank_marker[tolower(rank_hint)]) else NA_character_
    if (!is.na(marker)) {
      out$canonical_name <- paste(binomial, marker, ie)
      out$canonical_rank <- unname(.marker_rank[marker])
    } else {
      out$canonical_name <- paste(binomial, ie)
      out$canonical_rank <- "infraspecific_unmarked"
    }
    return(out)
  }
  # 3. Anything else after the epithet is authorship.
  out$authority_removed <- TRUE
  out
}

# Vectorised wrapper. Unmarked trinomials without a rank hint are resolved
# from marked names elsewhere in the same input when exactly one rank exists.
canonicalise_names <- function(names, rank_hint = NA_character_) {
  names <- as.character(names)
  if (!length(names)) {
    return(data.table::data.table(
      input_name = character(), canonical_name = character(), canonical_binomial = character(),
      canonical_rank = character(), infraspecific_epithet = character(), authority_removed = logical()
    ))
  }
  rank_hint <- rep_len(rank_hint, length(names))
  res <- data.table::rbindlist(Map(.parse_one, names, rank_hint))
  res[, input_name := names]
  marked <- unique(res[canonical_rank %chin% unname(.marker_rank),
    .(canonical_binomial, infraspecific_epithet, canonical_name, canonical_rank)])
  marked <- marked[, if (.N == 1L) .SD, by = .(canonical_binomial, infraspecific_epithet)]
  un <- res$canonical_rank == "infraspecific_unmarked"
  if (any(un) && nrow(marked) > 0L) {
    hit <- marked[res[un], on = .(canonical_binomial, infraspecific_epithet),
                  .(canonical_name = x.canonical_name, canonical_rank = x.canonical_rank)]
    ok <- !is.na(hit$canonical_name)
    idx <- which(un)[ok]
    data.table::set(res, idx, "canonical_name", hit$canonical_name[ok])
    data.table::set(res, idx, "canonical_rank", hit$canonical_rank[ok])
  }
  data.table::setcolorder(res, "input_name")
  res[]
}
