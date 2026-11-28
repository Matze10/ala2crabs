# AFD name authority for the COI run (Example 2: arthropods without crustaceans),
# built from the ALA species index in the AFD checklist CSV format read_afd()
# expects (helpers_synonyms.R). Counterpart of prepare_apc.R for plants.
#
# Why: the AFD website exports at most 1,000 taxa per download (Arthropoda has
# ~120,000), and the full AFD copy on ChecklistBank (2022) lacks whole families.
# ALA builds its animal taxonomy from the current AFD: every taxon carries its
# AFD GUID (https://biodiversity.org.au/afd/taxa/<GUID>), the AFD classification
# and its AFD synonyms (field additionalNames).
#
# Source: ALA species index (BIE) bulk download (always the CURRENT index;
# older versions cannot be downloaded again),
#   https://bie-ws.ala.org.au/ws/download?q=*&fq=idxtype:TAXON
#     &fq=rk_phylum:ARTHROPODA&fq=-rk_subphylum:CRUSTACEA
#     &fq=rank:(species OR subspecies)
# Only rows with datasetName "AFD" are kept; ALA's own additions
# ("Taxonomic Additions - ... 2025 - Internal Lists", ~1,400 names) are dropped.
#
# Files (in ALA_Records/):
#   AFD-ALA-arthropoda-<YYYYMMDD>/bie_download.csv   raw download (not tracked)
#   AFD-ALA-arthropoda-<YYYYMMDD>.csv               authority file (not tracked)
#   AFD-ALA-arthropoda-<YYYYMMDD>_download_info.txt URL, time, MD5 checksums, row
#                                                    counts and structure of both CSVs (tracked)
# Columns: GENUS, SPECIES, SUBSPECIES, NAME_TYPE ("Valid Name" / "Synonym"),
#   RANK ("Species" / "Subspecies"), TAXON_GUID (AFD GUID of the accepted taxon,
#   shared by its synonyms), SCIENTIFIC_NAME, FAMILY, ORDER, CLASS, SUBPHYLUM,
#   NAME_LAST_UPDATE / TAXON_LAST_UPDATE (download date).
#
# Version pinning (helpers_authority_files.R), as in prepare_apc.R:
#   Rscript Scripts/prepare_afd_from_ala.R
#       Never changes the version. Checks the pinned file (taxonomy.authority_file
#       of the config). If it is missing, rebuilds it from the raw download, or
#       downloads the current index and keeps it only if its md5 equals the pinned
#       raw_md5 (otherwise stops: copy the pinned files from another machine).
#   Rscript Scripts/prepare_afd_from_ala.R --update
#       Downloads the current index and compares its structure (columns and the
#       values of taxonRank, taxonomicStatus) with the pinned download. Same
#       structure: writes AFD-ALA-arthropoda-<today> with its tracking file and
#       sets taxonomy.authority_file in the config (and in Scripts/config/ala_crabs.yml
#       if that pointed to the same file). Then commit the config(s) and the
#       tracking file and push. Different structure: stops; nothing is changed.
#   --config <file>   config with the pinned file (default Scripts/config/ala_crabs_coi_arthropods.yml)
#
# Usage (from the project root): Rscript Scripts/prepare_afd_from_ala.R [--update] [--config <file>]

suppressMessages(library(data.table))
source("Scripts/workflow_helpers/helpers_authority_files.R")

url <- paste0(
  "https://bie-ws.ala.org.au/ws/download?q=*",
  "&fq=idxtype:TAXON",
  "&fq=rk_phylum:ARTHROPODA",
  "&fq=-rk_subphylum:CRUSTACEA",
  "&fq=rank:(species%20OR%20subspecies)"
)
raw_value_cols <- c("taxonRank", "taxonomicStatus")  # values the conversion depends on
active_config <- "Scripts/config/ala_crabs.yml"

args <- commandArgs(trailingOnly = TRUE)
update <- "--update" %in% args
config_file <- if ("--config" %in% args) args[which(args == "--config") + 1L] else "Scripts/config/ala_crabs_coi_arthropods.yml"

paths_for <- function(stamp) {
  raw_dir <- file.path("ALA_Records", paste0("AFD-ALA-arthropoda-", stamp))
  list(stamp = stamp, raw_dir = raw_dir, raw = file.path(raw_dir, "bie_download.csv"),
       out = file.path("ALA_Records", paste0("AFD-ALA-arthropoda-", stamp, ".csv")))
}

download_raw <- function() {
  tmp <- tempfile(fileext = ".csv")
  message("downloading ", url)
  if (download.file(url, tmp, mode = "wb", quiet = TRUE) != 0) stop("download failed")
  list(tmp = tmp, md5 = file_md5(tmp), time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
}

# ---- Conversion: raw ALA download -> AFD checklist format ---------------------
convert <- function(raw_file, out_file, stamp) {
  bie <- fread(raw_file, colClasses = "character", na.strings = "", encoding = "UTF-8")
  need <- c("taxonID", "taxonRank", "scientificName", "taxonomicStatus", "datasetName",
            "additionalNames", "genus", "family", "order", "class", "subphylum", "phylum")
  miss <- setdiff(need, names(bie))
  if (length(miss)) stop("ALA download is missing columns: ", paste(miss, collapse = ", "))
  message("ALA rows: ", nrow(bie))

  # Checks on the filter: arthropods only, no crustaceans, accepted AFD taxa.
  stopifnot(all(bie$phylum == "ARTHROPODA"), !any(bie$subphylum %chin% "CRUSTACEA", na.rm = TRUE))
  n_other <- bie[datasetName != "AFD" | is.na(datasetName), .N]
  bie <- bie[datasetName == "AFD" & taxonomicStatus == "accepted" & grepl("biodiversity.org.au/afd/taxa/", taxonID)]
  message("dropped ", n_other, " non-AFD rows (ALA taxonomic additions); AFD rows: ", nrow(bie))

  # Canonical name: without subgenus in parentheses and without authors
  # (scientificName has no authors in the BIE download).
  canonical <- function(x) trimws(gsub("[[:space:]]+", " ", gsub("[(][^)]*[)]", " ", x)))
  split_name <- function(cn) {
    p <- strsplit(cn, " ", fixed = TRUE)
    data.table(GENUS = vapply(p, `[`, "", 1), SPECIES = vapply(p, `[`, "", 2),
               SUBSPECIES = vapply(p, function(z) if (length(z) >= 3) z[3] else NA_character_, ""),
               n_words = lengths(p))
  }

  bie[, TAXON_GUID := sub("^.*/afd/taxa/", "", taxonID)]
  bie[, cn := canonical(scientificName)]
  higher <- bie[, .(TAXON_GUID, FAMILY = family, ORDER = order, CLASS = class, SUBPHYLUM = subphylum)]

  valid <- cbind(bie[, .(TAXON_GUID, SCIENTIFIC_NAME = scientificName, cn)], split_name(bie$cn))
  valid[, RANK := fifelse(n_words >= 3, "Subspecies", "Species")]
  valid[, NAME_TYPE := "Valid Name"]
  bad <- valid[n_words < 2 | n_words > 3]
  if (nrow(bad)) message("skipped ", nrow(bad), " accepted names that are not binomials/trinomials, e.g. ",
                         paste(head(bad$SCIENTIFIC_NAME, 3), collapse = "; "))
  valid <- valid[n_words %in% 2:3]

  # Synonyms: additionalNames is a comma-separated list of the AFD synonyms of
  # the taxon (original combinations, junior synonyms). Only species-group names
  # are kept (two or three words after removing the subgenus).
  syn <- bie[!is.na(additionalNames) & nzchar(additionalNames),
             .(SCIENTIFIC_NAME = trimws(unlist(strsplit(additionalNames, ",", fixed = TRUE)))), by = TAXON_GUID]
  syn <- syn[nzchar(SCIENTIFIC_NAME)]
  syn[, cn := canonical(SCIENTIFIC_NAME)]
  syn <- cbind(syn, split_name(syn$cn))
  syn <- syn[n_words %in% 2:3 & grepl("^[A-Z][a-z-]+$", GENUS) & grepl("^[a-z-]+$", SPECIES) &
               (is.na(SUBSPECIES) | grepl("^[a-z-]+$", SUBSPECIES))]
  syn[, RANK := fifelse(n_words == 3, "Subspecies", "Species")]
  syn[, NAME_TYPE := "Synonym"]
  syn <- syn[!valid, on = .(TAXON_GUID, cn)]  # same canonical name as the valid name

  cols <- c("GENUS", "SPECIES", "SUBSPECIES", "NAME_TYPE", "RANK", "TAXON_GUID", "SCIENTIFIC_NAME")
  out <- unique(rbind(valid[, ..cols], syn[, ..cols]))
  out <- merge(out, higher, by = "TAXON_GUID", all.x = TRUE, sort = FALSE)
  out[, `:=`(NAME_LAST_UPDATE = stamp, TAXON_LAST_UPDATE = stamp)]
  setcolorder(out, c("GENUS", "SPECIES", "SUBSPECIES", "NAME_TYPE", "RANK", "TAXON_GUID", "SCIENTIFIC_NAME",
                     "FAMILY", "ORDER", "CLASS", "SUBPHYLUM"))
  setorder(out, GENUS, SPECIES, SUBSPECIES, -NAME_TYPE)
  fwrite(out, out_file)
  message("written: ", out_file)
  message("  valid names: ", out[NAME_TYPE == "Valid Name", .N], " (",
          out[NAME_TYPE == "Valid Name" & RANK == "Species", .N], " species, ",
          out[NAME_TYPE == "Valid Name" & RANK == "Subspecies", .N], " subspecies)")
  message("  synonyms:    ", out[NAME_TYPE == "Synonym", .N])
  print(out[NAME_TYPE == "Valid Name", .N, by = CLASS][order(-N)])
  invisible(out)
}

# ---- Pinned version ----------------------------------------------------------
pinned <- yaml::read_yaml(config_file)$taxonomy$authority_file
if (is.null(pinned) || !nzchar(pinned)) stop("no taxonomy.authority_file in ", config_file)
info <- read_download_info(pinned)
if (is.null(info) || any(is.na(info[c("md5", "raw_md5")]))) stop("no tracking file with md5 and raw_md5 lines for ", pinned)
p <- paths_for(sub("^AFD-ALA-arthropoda-([0-9]+)[.]csv$", "\\1", basename(pinned)))

# Run as a function so that source() returns instead of ending the R session.
main <- function() {
if (!update) {
  if (file.exists(pinned) && file_md5(pinned) == info["md5"]) {
    message("ok: ", pinned, " is the pinned version (md5 ", info["md5"], ")"); return(invisible(NULL))
  }
  if (file.exists(pinned)) stop(pinned, " does not match its pinned md5 ", info["md5"], ". Restore the pinned file.")
  if (!file.exists(p$raw)) {
    d <- download_raw()
    if (d$md5 != info["raw_md5"]) {
      stop("the current ALA index (md5 ", d$md5, ") is not the pinned download (raw_md5 ", info["raw_md5"],
           "); ALA no longer provides it. Copy ", p$raw, " or ", pinned, " from another machine, or ",
           "switch to the current index with --update.")
    }
    dir.create(p$raw_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(d$tmp, p$raw)
  }
  if (file_md5(p$raw) != info["raw_md5"]) stop(p$raw, " does not match its pinned raw_md5 ", info["raw_md5"])
  convert(p$raw, pinned, p$stamp)
  if (file_md5(pinned) != info["md5"]) stop("rebuilt ", pinned, " does not match its pinned md5 ", info["md5"])
  message("rebuilt pinned version: ", pinned)
  return(invisible(NULL))
}

# ---- --update ----------------------------------------------------------------
pinned_sig <- info[grepl("^raw_(columns|values_)", names(info))]
if (!length(pinned_sig)) {
  if (!file.exists(p$raw)) stop("no raw_ structure lines in the tracking file of ", pinned, " and ", p$raw, " is missing")
  pinned_sig <- structure_signature(p$raw, raw_value_cols, prefix = "raw_")
  pinned_sig <- pinned_sig[names(pinned_sig) != "raw_rows"]
}
d <- download_raw()
if (d$md5 == info["raw_md5"]) { message("the current ALA index is the pinned download; nothing to update"); return(invisible(NULL)) }
new_sig <- structure_signature(d$tmp, raw_value_cols, prefix = "raw_")
diffs <- compare_signature(pinned_sig, new_sig)
if (length(diffs)) {
  stop("the current ALA index has a different structure than the pinned download; nothing changed. ",
       "Differences:\n  ", paste(diffs, collapse = "\n  "), "\nCheck the conversion before using it.")
}
n <- paths_for(format(Sys.Date(), "%Y%m%d"))
if (file.exists(n$raw_dir) || file.exists(n$out)) stop(n$raw_dir, " or ", n$out, " exists already; not overwritten")
dir.create(n$raw_dir, recursive = TRUE)
file.copy(d$tmp, n$raw)
out <- convert(n$raw, n$out, n$stamp)
write_download_info(n$out, c(url = url, downloaded = d$time,
                             raw_file = paste0(basename(n$raw_dir), "/", basename(n$raw)), raw_md5 = d$md5, new_sig,
                             file = basename(n$out), md5 = file_md5(n$out), rows = as.character(nrow(out))))
changed <- config_file
set_authority_file(config_file, n$out)
if (file.exists(active_config) && normalizePath(active_config) != normalizePath(config_file) &&
    identical(yaml::read_yaml(active_config)$taxonomy$authority_file, pinned)) {
  set_authority_file(active_config, n$out); changed <- c(changed, active_config)
}
message("same structure as the pinned download; new version: ", n$out)
message("taxonomy.authority_file set to ", n$out, " in: ", paste(changed, collapse = ", "))
message("commit and push: git add ", paste(c(changed, download_info_path(n$out)), collapse = " "))
}

main()
