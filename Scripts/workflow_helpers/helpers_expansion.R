# Genus expansion and NCBI phrase names for included genus-level decisions.
#
# Genus-level decisions (e.g. "Acacia sp. New England (J.B.Williams 97011)" ->
# Acacia, rank genus) do not go into the CRABS species list, because a genus
# line would download every sequence of the genus worldwide. Two generic steps
# recover what is usable from them:
#   1. NCBI phrase names (crabs.keep_ncbi_phrase_names): the phrase name itself
#      is looked up in the local NCBI Taxonomy dump (below taxon.ncbi_root_taxid).
#      If NCBI has it (e.g. "Samadera sp. Mt Nardi"), the exact NCBI scientific
#      name goes into the list.
#   2. Genus expansion (crabs.expand_genus_from_authority): the name authority's
#      accepted species of the genus whose distribution includes one of
#      taxonomy.authority_states are added, if they are not in the list yet.
#      They are mapped to NCBI like match_names(): accepted name, authority
#      synonyms, Latin ending variants. Needs an authority with a species table
#      (APC now; AFD once its reader exists); with authority: none it is skipped.
# Expanded species have no ALA records of their own, so they get no region flags.

# Decisions counted as included (crabs.include_decisions, default include).
include_decision_values <- function(config) {
  x <- as.character(unlist(config$crabs$include_decisions))
  if (!length(x)) "include" else x
}

# Comparison key for phrase names: lower case, letters/digits only, single spaces.
phrase_key <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9 ]", " ", x)
  trimws(gsub(" +", " ", x))
}

# Phrase key without the voucher in brackets and a trailing "<state> Herbarium".
phrase_key_short <- function(x) {
  x <- sub("[(].*$", "", x)
  x <- sub(" *(WA|NSW|Qld|QLD|NT|Vic|VIC|SA|Tas|TAS|ACT)? *Herbarium *'?$", "", x)
  phrase_key(x)
}

# Look up included genus-level phrase names in the NCBI name index.
ncbi_phrase_matches <- function(genus_rows, ncbi) {
  if (!nrow(genus_rows)) return(data.table::data.table())
  nm <- ncbi$names[grepl(" sp[.] ", name)]
  nm[, `:=`(genus = sub(" .*", "", name), k_full = phrase_key(name), k_short = phrase_key_short(name))]
  q <- unique(rbind(
    genus_rows[, .(normalised_name, genus, query = normalised_name)],
    genus_rows[!is.na(canonical_name), .(normalised_name, genus, query = canonical_name)]
  ))
  q <- q[grepl(" sp[.] ", query)]
  q[, `:=`(k_full = phrase_key(query), k_short = phrase_key_short(query))]
  # The short key needs at least one word after "sp", so "Acacia sp." alone never matches.
  q <- q[grepl(" sp [a-z0-9]", k_short) | grepl(" sp [a-z0-9]", k_full)]
  full <- merge(q, nm[, .(genus, k_full, taxid)], by = c("genus", "k_full"))[, via := "ncbi_phrase_exact"]
  short <- merge(q[!normalised_name %chin% full$normalised_name], nm[, .(genus, k_short, taxid)],
                 by = c("genus", "k_short"))[, via := "ncbi_phrase_without_voucher"]
  hits <- rbind(full[, .(normalised_name, genus, taxid, via)], short[, .(normalised_name, genus, taxid, via)])
  if (!nrow(hits)) return(data.table::data.table())
  hits <- merge(hits, ncbi$taxa[, .(taxid, ncbi_name)], by = "taxid")
  hits[, .(
    source = "ncbi_phrase_name", genus = genus[1], source_names = paste(sort(unique(normalised_name)), collapse = "; "),
    authority_accepted_name = NA_character_, authority_taxon_id = NA_character_, distribution = NA_character_,
    ncbi_taxid = paste(unique(taxid), collapse = ";"), ncbi_match_via = via[1]
  ), by = .(crabs_name = ncbi_name)]
}

# Authority species of the genera, in the configured states, mapped to NCBI.
authority_genus_species <- function(genus_rows, auth, ncbi, config, covered) {
  if (is.null(auth$species) || !nrow(auth$species)) {
    message("genus expansion skipped: no name authority species table (taxonomy.authority: ",
            authority_name(config), ")")
    return(data.table::data.table())
  }
  states <- as.character(unlist(config$taxonomy$authority_states))
  bad <- tolower(as.character(unlist(config$taxonomy$authority_state_qualifiers_excluded)))
  src <- genus_rows[, .(source_names = paste(sort(unique(normalised_name)), collapse = "; ")), by = genus]
  sp <- data.table::copy(auth$species[genus %chin% src$genus])
  if (!nrow(sp)) return(data.table::data.table())
  has_dist <- any(nzchar(auth$species$distribution_parsed))
  if (has_dist) {
    if (!length(states)) stop("taxonomy.authority_states is empty; set the state codes used in the authority's distribution field.")
    # distribution_parsed: "code|qualifier" entries separated by ";".
    sp[, in_states := vapply(strsplit(distribution_parsed, ";", fixed = TRUE), function(e) {
      if (!length(e)) return(FALSE)
      code <- sub("[|].*$", "", e); qual <- tolower(sub("^[^|]*[|]?", "", e))
      any(code %chin% states & !(qual %chin% bad))
    }, logical(1))]
  } else if (isTRUE(config$crabs$expand_genus_without_distribution)) {
    # Authority without distributions (e.g. the AFD checklist): every accepted
    # species of the genus in the authority counts (Australia-wide).
    message("genus expansion without distribution data: all ", authority_name(config), " species of the genera are used.")
    sp[, in_states := TRUE]
  } else {
    message("genus expansion skipped: ", authority_name(config), " file has no distribution data ",
            "(set crabs.expand_genus_without_distribution: true to add all its species of the genus).")
    return(data.table::data.table())
  }
  sp <- sp[in_states == TRUE & !(authority_accepted_name %chin% covered)]
  if (!nrow(sp)) return(data.table::data.table())
  tries <- rbind(
    sp[, .(authority_accepted_name, step = 1L, via = "authority_accepted", search = authority_accepted_name)],
    merge(sp[, .(authority_accepted_name)], auth$synonyms, by = "authority_accepted_name",
          allow.cartesian = TRUE)[, .(authority_accepted_name, step = 3L, via = "authority_synonym", search = synonym)]
  )
  ev <- ending_variants(sp$authority_accepted_name)
  tries <- rbind(tries, data.table::data.table(authority_accepted_name = rep(sp$authority_accepted_name, lengths(ev)),
                                               step = 4L, via = "ending_variant", search = unlist(ev)))
  hits <- merge(tries, ncbi$names, by.x = "search", by.y = "name", allow.cartesian = TRUE)
  hits <- merge(hits, ncbi$taxa[, .(taxid, ncbi_name)], by = "taxid", all.x = TRUE)
  data.table::setorder(hits, authority_accepted_name, step, search)
  best <- hits[, .SD[step == step[1] & search == search[1]], by = authority_accepted_name][
    , .(ncbi_name = ncbi_name[1], ncbi_taxid = paste(unique(taxid), collapse = ";"), ncbi_match_via = via[1]),
    by = authority_accepted_name]
  out <- merge(sp, best, by = "authority_accepted_name", all.x = TRUE)
  out <- merge(out, src, by = "genus")
  out[is.na(ncbi_match_via), ncbi_match_via := "no_match"]
  out[, .(source = "authority_genus_expansion", genus, source_names, authority_accepted_name, authority_taxon_id,
          distribution, crabs_name = data.table::fcoalesce(ncbi_name, authority_accepted_name), ncbi_taxid, ncbi_match_via)]
}

# Target genus_expansion (build_genus_expansion): writes the derived table and its audit copy.
build_genus_expansion <- function(matches_file, authority_file, taxdump_files, config, output_path, audit_path) {
  inc <- include_decision_values(config)
  m <- data.table::fread(matches_file, na.strings = c("", "NA"), colClasses = "character")
  included <- m[decision %chin% inc & !is.na(accepted_name) & nzchar(accepted_name)]
  genus_rows <- included[(tolower(accepted_rank) == "genus" | canonical_rank == "genus_or_uncertain")]
  genus_rows[, genus := sub(" .*", "", sub("^x ", "", trimws(accepted_name)))]
  genus_rows <- genus_rows[grepl("^[A-Z][a-z-]+$", genus)]
  # Names already covered by the list (accepted, CRABS, authority and binomial forms).
  covered <- unique(na.omit(c(included$accepted_name, included$crabs_name, included$authority_accepted_name,
                              included$canonical_binomial)))
  # Phrase names among them (an automatic include whose NCBI name is a phrase
  # name) are written by build_crabs_outputs() only via keep_ncbi_phrase_names,
  # so they do not count as covered here.
  covered <- covered[!grepl(" sp[.] ", covered)]
  ncbi <- read_ncbi_names(taxdump_files, ncbi_root(config), NA_character_)
  parts <- list()
  if (isTRUE(config$crabs$keep_ncbi_phrase_names)) parts$phrase <- ncbi_phrase_matches(genus_rows, ncbi)
  if (isTRUE(config$crabs$expand_genus_from_authority)) {
    if (authority_name(config) == "NONE") {
      message("genus expansion skipped: no name authority (taxonomy.authority: none)")
    } else {
      auth <- read_authority(authority_file, config)
      parts$authority <- authority_genus_species(genus_rows, auth, ncbi, config, covered)
    }
  }
  cols <- c("source", "genus", "source_names", "authority_accepted_name", "authority_taxon_id", "distribution",
            "crabs_name", "ncbi_taxid", "ncbi_match_via")
  out <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  if (!nrow(out)) out <- data.table::as.data.table(stats::setNames(replicate(length(cols), character(), simplify = FALSE), cols))
  out <- out[, ..cols]
  out <- out[!is.na(crabs_name) & !(crabs_name %chin% covered)]
  out <- unique(out, by = "crabs_name")
  out[, already_in_list := FALSE]
  data.table::setorder(out, source, genus, crabs_name)
  message("Genus expansion: ", sum(out$source == "authority_genus_expansion"), " authority species from ",
          data.table::uniqueN(genus_rows$genus), " genera; ", sum(out$source == "ncbi_phrase_name"), " NCBI phrase names.")
  for (p in c(output_path, audit_path)) {
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(out, p)
  }
  c(output_path, audit_path)
}
