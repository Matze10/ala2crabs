# APC and NCBI name resolution for the ALA -> CRABS workflow.
#
# Every canonical name in the review queue (helpers_canonical_names.R) is
#   1. looked up in the national name authority (taxonomy.authority): the
#      Australian Plant Census (APC) for plants, the Australian Faunal Directory
#      (AFD) for animals, or none (NCBI only). For APC: Accepted
#      names are kept; synonyms and orthographic variants are replaced by the
#      APC accepted name. Misapplied, pro parte, doubtful and excluded APC
#      entries are not followed (the ALA name is used and the status flagged).
#   2. looked up in the NCBI Taxonomy dump below taxon.ncbi_root_taxid
#      (35493 land plants, 33208 animals), trying in order: the accepted name, the ALA canonical name, the other APC
#      synonyms of the accepted taxon, and Latin gender-ending variants
#      (-us/-a/-um, -is/-e, -er/-ra/-rum) of the name. The first hit wins and
#      the NCBI scientific name of that taxid becomes the CRABS name, because
#      CRABS queries NCBI with "<name>"[Organism].
#
# Inputs (config$taxonomy): authority, authority_file (manual download, see
# ALA_Records/README.md), taxon.ncbi_root_taxid (helpers_io.R),
# ncbi_taxdump_dir (taxdump.tar.gz downloaded once into this folder).
# Base R + data.table only.

`%||%` <- function(x, y) if (is.null(x)) y else x

ncbi_taxdump_url <- "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# APC statuses whose accepted name is followed.
apc_synonym_statuses <- c(
  "taxonomic synonym", "nomenclatural synonym", "orthographic variant",
  "basionym", "replaced synonym", "alternative name"
)
apc_ranks <- c("Species", "Subspecies", "Variety", "Form")

# ---- NCBI taxdump ----------------------------------------------------------

# Download and unpack names.dmp and nodes.dmp once. Delete the folder (or set
# taxonomy.ncbi_taxdump_refresh to a new date) to fetch a newer dump.
fetch_ncbi_taxdump <- function(config) {
  dir <- here::here(config$taxonomy$ncbi_taxdump_dir)
  files <- file.path(dir, c("names.dmp", "nodes.dmp"))
  if (all(file.exists(files))) return(files)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  archive <- file.path(dir, "taxdump.tar.gz")
  old <- options(timeout = max(1800, getOption("timeout")))
  on.exit(options(old), add = TRUE)
  utils::download.file(ncbi_taxdump_url, archive, mode = "wb", quiet = TRUE)
  utils::untar(archive, files = c("names.dmp", "nodes.dmp"), exdir = dir)
  unlink(archive)
  writeLines(format(Sys.Date()), file.path(dir, "download_date.txt"))
  if (!all(file.exists(files))) stop("NCBI taxdump download did not produce names.dmp and nodes.dmp in ", dir)
  files
}

# Read a .dmp file: fields are separated by "\t|\t", so splitting on tabs puts
# the values in the odd columns.
read_dmp <- function(path, cols, col_names) {
  x <- data.table::fread(path, sep = "\t", header = FALSE, quote = "", select = cols,
                         colClasses = "character", encoding = "UTF-8", showProgress = FALSE)
  data.table::setnames(x, col_names)
  x
}

# All taxids under `root` (35493 = Streptophyta, 58024 = Spermatophyta).
ncbi_descendants <- function(nodes, root) {
  children <- split(nodes$taxid, nodes$parent)
  found <- root; frontier <- root
  while (length(frontier)) {
    frontier <- unlist(children[frontier], use.names = FALSE)
    frontier <- frontier[!frontier %chin% found]
    found <- c(found, frontier)
  }
  found
}

# Name index for the configured NCBI subtree (taxon.ncbi_root_taxid): one row per (name, taxid) for scientific names and
# synonyms, plus the scientific name and rank of each taxid.
read_ncbi_names <- function(taxdump_files, root = "35493", seed_root = "58024") {
  names_file <- taxdump_files[basename(taxdump_files) == "names.dmp"]
  nodes_file <- taxdump_files[basename(taxdump_files) == "nodes.dmp"]
  nodes <- read_dmp(nodes_file, c(1, 3, 5), c("taxid", "parent", "ncbi_rank"))
  nodes <- nodes[taxid != parent]
  plants <- ncbi_descendants(nodes, as.character(root))
  # taxon.seed_plant_taxid (default 58024 Spermatophyta) only fills the audit
  # column ncbi_seed_plant; NA skips it (column NA).
  seed <- if (is.na(seed_root)) NULL else ncbi_descendants(nodes, as.character(seed_root))
  nm <- read_dmp(names_file, c(1, 3, 7), c("taxid", "name", "name_class"))
  nm <- nm[taxid %chin% plants & name_class %chin% c(
    "scientific name", "synonym", "equivalent name", "authority", "misspelling", "genbank synonym"
  )]
  sci <- nm[name_class == "scientific name", .(taxid, ncbi_name = name)]
  sci <- merge(sci, nodes[, .(taxid, ncbi_rank)], by = "taxid", all.x = TRUE)
  sci[, ncbi_seed_plant := if (is.null(seed)) NA else taxid %chin% seed]
  list(names = unique(nm[, .(name, taxid)]), taxa = sci)
}

# ---- Name authority ---------------------------------------------------------

authority_name <- function(config) {
  a <- toupper(trimws(as.character(config$taxonomy$authority %||% "none")))
  if (!a %in% c("APC", "AFD", "NONE")) stop("taxonomy.authority must be APC, AFD or none, not: ", a)
  a
}
# NCBI subtree for name matching: taxon.ncbi_root_taxid, else the old
# taxonomy.ncbi_root_taxid, else 35493 (see taxon_settings() in helpers_io.R).
ncbi_root <- function(config) taxon_settings(config)$ncbi_root_taxid

# File tracked by the authority_file target. With authority: none there is no
# file, so the config itself is tracked instead.
authority_input <- function(config) {
  if (authority_name(config) == "NONE") return(here::here("Scripts/config/ala_crabs.yml"))
  path <- here::here(config$taxonomy$authority_file)
  if (!file.exists(path)) stop("taxonomy.authority_file not found: ", path,
                               ". Download it by hand, see 'Name authority file' in README.md.")
  path
}

# Every reader returns list(names, synonyms, species, version, label):
#   names:    canonical_name, authority_status, authority_accepted_name, authority_taxon_id
#   synonyms: authority_accepted_name, synonym
#   species:  accepted species for the genus expansion (helpers_expansion.R):
#             authority_accepted_name (canonical, species rank), genus,
#             authority_taxon_id, distribution (raw text) and
#             distribution_parsed ("code|qualifier" entries separated by ";").
#             NULL if the authority has no distribution data (authority: none).
read_authority <- function(path, config) {
  a <- authority_name(config)
  res <- switch(a,
    APC = read_apc(path),
    AFD = read_afd(path),
    NONE = list(
      names = data.table::data.table(canonical_name = character(), authority_status = character(),
                                     authority_accepted_name = character(), authority_taxon_id = character()),
      synonyms = data.table::data.table(authority_accepted_name = character(), synonym = character()),
      species = NULL,
      version = NA_character_
    )
  )
  res$label <- if (a == "NONE") NA_character_ else a
  res
}

# Australian Faunal Directory checklist export (the CSV download on an AFD taxon
# page, e.g. https://biodiversity.org.au/afd/taxa/VESPIDAE). Columns used:
# GENUS, SPECIES, SUBSPECIES (canonical name, without subgenus and authors),
# NAME_TYPE ("Valid Name" = accepted; other types such as "Synonym" or
# "Generic Combination" are synonyms), RANK (Species / Subspecies), TAXON_GUID
# (shared by a valid name and its synonyms; used as authority_taxon_id) and
# NAME_LAST_UPDATE / TAXON_LAST_UPDATE (version). The file is Latin-1 encoded.
# Same contract as read_apc(). The AFD checklist has no distribution column, so
# the species table has no state codes: unless a DISTRIBUTION or STATES column
# is present (comma-separated, codes mapped as in afd_state_code()), the genus
# expansion has nothing to filter by and is skipped for AFD.
read_afd <- function(path) {
  afd <- data.table::fread(path, colClasses = "character", na.strings = "", encoding = "Latin-1")
  need <- c("GENUS", "SPECIES", "NAME_TYPE", "RANK", "TAXON_GUID")
  miss <- setdiff(need, names(afd))
  if (length(miss)) stop("AFD file is missing columns: ", paste(miss, collapse = ", "), " (", path, ")")
  if (!"SUBSPECIES" %in% names(afd)) afd[, SUBSPECIES := NA_character_]
  afd <- afd[!is.na(GENUS) & !is.na(SPECIES) & RANK %chin% c("Species", "Subspecies")]
  afd[, canonical_name := trimws(paste(GENUS, SPECIES, data.table::fifelse(is.na(SUBSPECIES), "", SUBSPECIES)))]
  afd[, accepted := NAME_TYPE == "Valid Name"]
  valid <- unique(afd[accepted == TRUE, .(TAXON_GUID, authority_accepted_name = canonical_name)], by = "TAXON_GUID")
  afd <- merge(afd, valid, by = "TAXON_GUID", all.x = TRUE, sort = FALSE)
  # Same resolution as for APC: an accepted entry wins, then a single synonym target.
  res <- afd[, {
    acc <- unique(authority_accepted_name[accepted])
    syn <- unique(authority_accepted_name[!accepted & !is.na(authority_accepted_name)])
    if (length(acc) == 1L) {
      list(authority_status = "accepted", authority_accepted_name = acc, authority_taxon_id = TAXON_GUID[accepted][1])
    } else if (length(syn) == 1L) {
      list(authority_status = "synonym", authority_accepted_name = syn, authority_taxon_id = TAXON_GUID[!accepted][1])
    } else if (length(acc) + length(syn) > 1L) {
      list(authority_status = "ambiguous", authority_accepted_name = NA_character_, authority_taxon_id = NA_character_)
    } else {
      list(authority_status = tolower(NAME_TYPE[1]), authority_accepted_name = NA_character_, authority_taxon_id = NA_character_)
    }
  }, by = canonical_name]
  syn <- unique(afd[!accepted & !is.na(authority_accepted_name) & canonical_name != authority_accepted_name,
                    .(authority_accepted_name, synonym = canonical_name)])
  dist_col <- intersect(c("DISTRIBUTION", "STATES"), names(afd))
  species <- afd[accepted == TRUE & RANK == "Species",
                 .(authority_accepted_name = canonical_name, genus = GENUS, authority_taxon_id = TAXON_GUID,
                   distribution = if (length(dist_col)) get(dist_col[1]) else NA_character_)]
  species <- unique(species, by = "authority_accepted_name")
  species[, distribution_parsed := vapply(distribution, function(d) {
    if (is.na(d) || !nzchar(d)) return("")
    p <- strsplit(parse_distribution(d), ";", fixed = TRUE)[[1]]
    paste(paste0(afd_state_code(sub("[|].*$", "", p)), "|", sub("^[^|]*[|]?", "", p)), collapse = ";")
  }, character(1), USE.NAMES = FALSE)]
  # Without a distribution column the species table is kept with empty
  # distributions; the genus expansion then uses it only with
  # crabs.expand_genus_without_distribution: true (all Australian species).
  if (!any(nzchar(species$distribution_parsed))) {
    message("AFD file has no distribution column: no state data for the genus expansion.")
  }
  dates <- c(afd$NAME_LAST_UPDATE, afd$TAXON_LAST_UPDATE)
  dates <- dates[!is.na(dates)]
  version <- if (length(dates)) {
    d <- max(substr(dates, 1, 8)); paste(substr(d, 1, 4), substr(d, 5, 6), substr(d, 7, 8), sep = "-")
  } else sub("^AFD-([0-9]{4})([0-9]{2})([0-9]{2}).*$", "\\1-\\2-\\3", basename(path))
  list(names = res, synonyms = syn, species = species, version = version)
}

# AFD state names -> the APC codes used in taxonomy.authority_states.
afd_state_code <- function(x) {
  map <- c("new south wales" = "NSW", "queensland" = "Qld", "victoria" = "Vic", "tasmania" = "Tas",
           "south australia" = "SA", "western australia" = "WA", "northern territory" = "NT",
           "australian capital territory" = "ACT", "lord howe island" = "LHI", "norfolk island" = "NI",
           "christmas island" = "ChI", "nsw" = "NSW", "qld" = "Qld", "vic" = "Vic", "tas" = "Tas",
           "sa" = "SA", "wa" = "WA", "nt" = "NT", "act" = "ACT")
  k <- tolower(trimws(x))
  ifelse(k %in% names(map), unname(map[k]), trimws(x))
}

# ---- APC -------------------------------------------------------------------

# APC export -> one row per canonical name with its resolution.
read_apc <- function(apc_file) {
  apc <- data.table::fread(apc_file, colClasses = "character", na.strings = "", encoding = "UTF-8",
                           select = c("taxonID", "acceptedNameUsage", "taxonomicStatus", "scientificName",
                                      "canonicalName", "taxonRank", "modified", "nameType", "taxonDistribution"))
  # Accepted species with their distribution, for the genus expansion.
  species <- apc[taxonomicStatus == "accepted" & taxonRank == "Species" & nameType == "scientific" & !is.na(canonicalName),
                 .(authority_accepted_name = canonicalName, genus = sub(" .*", "", canonicalName),
                   authority_taxon_id = taxonID, distribution = taxonDistribution)]
  species <- unique(species, by = "authority_accepted_name")
  species[, distribution_parsed := parse_distribution(distribution)]
  apc <- apc[taxonRank %chin% apc_ranks & !is.na(canonicalName)]
  accepted <- apc[taxonomicStatus == "accepted", .(acceptedNameUsage = scientificName,
                                                   authority_accepted_name = canonicalName, authority_taxon_id = taxonID)]
  accepted <- unique(accepted, by = "acceptedNameUsage")
  apc <- merge(apc, accepted, by = "acceptedNameUsage", all.x = TRUE, sort = FALSE)
  apc[taxonomicStatus == "accepted", `:=`(authority_accepted_name = canonicalName, authority_taxon_id = taxonID)]
  apc[, kind := data.table::fcase(
    taxonomicStatus == "accepted", "accepted",
    taxonomicStatus %chin% apc_synonym_statuses & !is.na(authority_accepted_name), "synonym",
    default = "other"
  )]
  # Resolve each name: an accepted entry wins, then a single synonym target.
  res <- apc[, {
    acc <- unique(authority_accepted_name[kind == "accepted"])
    syn <- unique(authority_accepted_name[kind == "synonym"])
    if (length(acc) == 1L) {
      list(authority_status = "accepted", authority_accepted_name = acc, authority_taxon_id = authority_taxon_id[kind == "accepted"][1])
    } else if (length(syn) == 1L) {
      list(authority_status = "synonym", authority_accepted_name = syn, authority_taxon_id = authority_taxon_id[kind == "synonym"][1])
    } else if (length(acc) + length(syn) > 1L) {
      list(authority_status = "ambiguous", authority_accepted_name = NA_character_, authority_taxon_id = NA_character_)
    } else {
      list(authority_status = taxonomicStatus[1], authority_accepted_name = NA_character_, authority_taxon_id = NA_character_)
    }
  }, by = .(canonical_name = canonicalName)]
  # Synonyms of every accepted taxon, used as extra NCBI search names.
  syn <- unique(apc[kind == "synonym", .(authority_accepted_name, synonym = canonicalName)])
  version <- max(substr(apc$modified, 1, 10), na.rm = TRUE)
  list(names = res, synonyms = syn, species = species, version = version)
}

# APC distribution text, e.g. "WA, NT, Qld, NSW (naturalised), LHI", as
# "code|qualifier" entries separated by ";" (qualifier empty if none):
# "WA|;NT|;Qld|;NSW|naturalised;LHI|".
parse_distribution <- function(x) {
  vapply(x, function(d) {
    if (is.na(d) || !nzchar(trimws(d))) return("")
    e <- trimws(strsplit(d, ",", fixed = TRUE)[[1]])
    e <- e[nzchar(e)]
    code <- trimws(sub("[(].*$", "", e))
    qual <- ifelse(grepl("[(]", e), trimws(sub("[)].*$", "", sub("^[^(]*[(]", "", e))), "")
    paste(paste0(code, "|", qual), collapse = ";")
  }, character(1), USE.NAMES = FALSE)
}

# ---- Matching --------------------------------------------------------------

# Latin gender-ending variants of the last epithet (ALA reconstructed names
# often disagree with the genus, e.g. Styphelia lanceolatus).
ending_variants <- function(x) {
  swaps <- list(us = c("a", "um"), a = c("us", "um"), um = c("us", "a"),
                is = "e", e = "is", er = c("ra", "rum"), ra = c("er", "rum"))
  out <- lapply(x, function(n) {
    if (is.na(n)) return(character())
    for (end in names(swaps)) {
      if (grepl(paste0(end, "$"), n)) return(paste0(sub(paste0(end, "$"), "", n), swaps[[end]]))
    }
    character()
  })
  out
}

# Build the search order for every queue name and keep the first NCBI hit.
match_names <- function(queue_file, authority_file, taxdump_files, config, output_path, unmatched_path) {
  queue <- data.table::fread(queue_file, colClasses = "character",
                             select = c("normalised_name", "canonical_name", "canonical_binomial", "canonical_rank"))
  names_dt <- unique(queue[!is.na(canonical_name) & nzchar(canonical_name),
                           .(canonical_name, canonical_binomial, canonical_rank)], by = "canonical_name")
  apc <- read_authority(authority_file, config)
  names_dt <- merge(names_dt, apc$names, by = "canonical_name", all.x = TRUE, sort = FALSE)
  names_dt[is.na(authority_status), authority_status := "not_in_authority"]
  # Names missing from the authority: try their gender-ending variants in APC as well.
  miss <- names_dt[authority_status == "not_in_authority", canonical_name]
  if (length(miss)) {
    ev <- ending_variants(miss)
    alt <- data.table::data.table(canonical_name = rep(miss, lengths(ev)), variant = unlist(ev))
    alt <- merge(alt, apc$names[authority_status %chin% c("accepted", "synonym")],
                 by.x = "variant", by.y = "canonical_name")
    alt <- unique(alt, by = "canonical_name")
    names_dt[alt, on = "canonical_name", `:=`(
      authority_status = paste(i.authority_status, "(ending variant)"),
      authority_accepted_name = i.authority_accepted_name, authority_taxon_id = i.authority_taxon_id
    )]
  }

  tries <- rbind(
    names_dt[!is.na(authority_accepted_name), .(canonical_name, step = 1L, via = "authority_accepted", search = authority_accepted_name)],
    names_dt[, .(canonical_name, step = 2L, via = "ala_name", search = canonical_name)],
    merge(names_dt[!is.na(authority_accepted_name), .(canonical_name, authority_accepted_name)], apc$synonyms,
          by = "authority_accepted_name", allow.cartesian = TRUE)[, .(canonical_name, step = 3L, via = "authority_synonym", search = synonym)]
  )
  base <- names_dt[, .(canonical_name, base = data.table::fifelse(is.na(authority_accepted_name), canonical_name, authority_accepted_name))]
  ev <- ending_variants(base$base)
  tries <- rbind(tries, data.table::data.table(
    canonical_name = rep(base$canonical_name, lengths(ev)), step = 4L, via = "ending_variant", search = unlist(ev)
  ))

  ncbi <- read_ncbi_names(taxdump_files, ncbi_root(config), taxon_settings(config)$seed_plant_taxid)
  hits <- merge(tries, ncbi$names, by.x = "search", by.y = "name", allow.cartesian = TRUE)
  hits <- merge(hits, ncbi$taxa, by = "taxid", all.x = TRUE)
  data.table::setorder(hits, canonical_name, step, search)
  hits <- hits[, .SD[step == step[1] & search == search[1]], by = canonical_name]
  best <- hits[, .(
    ncbi_matched_on = search[1], ncbi_match_via = via[1],
    ncbi_taxid = paste(unique(taxid), collapse = ";"),
    ncbi_name = ncbi_name[1], ncbi_rank = ncbi_rank[1],
    ncbi_seed_plant = all(ncbi_seed_plant), ncbi_multiple_taxa = data.table::uniqueN(taxid) > 1L
  ), by = canonical_name]

  out <- merge(names_dt, best, by = "canonical_name", all.x = TRUE, sort = TRUE)
  out[is.na(ncbi_match_via), ncbi_match_via := "no_match"]
  # Name written to the CRABS list: NCBI name, else authority accepted name, else ALA canonical.
  out[, crabs_name := data.table::fcoalesce(ncbi_name, authority_accepted_name, canonical_name)]
  dump_date <- file.path(dirname(taxdump_files[1]), "download_date.txt")
  out[, `:=`(
    authority_label = apc$label,
    authority_version_date = apc$version,
    ncbi_taxdump_date = if (file.exists(dump_date)) readLines(dump_date, n = 1L) else NA_character_
  )]
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(unmatched_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(out, output_path)
  data.table::fwrite(out[ncbi_match_via == "no_match"], unmatched_path)
  c(output_path, unmatched_path)
}

# Add APC/NCBI columns to the decision table. Automatic includes take the
# CRABS name as accepted_name; manual decisions keep the reviewer's
# accepted_name, and only a blank or placeholder ("Author") authority is filled.
add_synonym_matches <- function(matches, synonym_file) {
  syn <- data.table::fread(synonym_file, colClasses = "character", na.strings = "")
  keep <- c("canonical_name", "authority_label", "authority_status", "authority_accepted_name", "authority_taxon_id", "authority_version_date",
            "ncbi_taxid", "ncbi_name", "ncbi_match_via", "ncbi_matched_on", "ncbi_seed_plant",
            "ncbi_taxdump_date", "crabs_name")
  matches <- merge(matches, syn[, ..keep], by = "canonical_name", all.x = TRUE, sort = FALSE)
  in_apc <- !is.na(matches$authority_accepted_name) & nzchar(matches$authority_accepted_name)
  in_ncbi <- !is.na(matches$ncbi_name) & nzchar(matches$ncbi_name)
  matches[, syn_authority := data.table::fcase(
    in_apc & in_ncbi, paste0(authority_label, "; NCBI Taxonomy"), in_apc, authority_label, in_ncbi, "NCBI Taxonomy", default = "ALA")]
  matches[, syn_version := data.table::fcase(
    in_apc & in_ncbi, paste0(authority_label, " ", authority_version_date, "; NCBI taxdump ", ncbi_taxdump_date),
    in_apc, paste0(authority_label, " ", authority_version_date), in_ncbi, paste0("NCBI taxdump ", ncbi_taxdump_date),
    default = "ALA download")]
  matches[decision_source == "automatic" & decision == "include" & !is.na(crabs_name), accepted_name := crabs_name]
  matches[(is.na(authority) | !nzchar(authority) | authority %chin% c("Author", "author")) & !is.na(crabs_name),
          `:=`(authority = syn_authority, authority_version = syn_version)]
  matches[, c("syn_authority", "syn_version") := NULL]
  matches
}
