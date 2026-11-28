# APC name authority for plant runs (Example 1): the "APC Taxon" CSV export of
# the Australian Plant Census, in the format read_apc() expects
# (helpers_synonyms.R). Counterpart of prepare_afd_from_ala.R for fauna.
#
# Source: National Species List export service, the same file as "APC Taxon" on
#   https://biodiversity.org.au/nsl/services/export/index
#   download URL: https://biodiversity.org.au/nsl/services/export/taxonCsv
# The service always returns the CURRENT export; older versions cannot be
# downloaded again.
#
# Version pinning (helpers_authority_files.R). The version in use is the file in
# taxonomy.authority_file of the config, pinned by the MD5 checksum and the
# structure in its tracking file ALA_Records/<file>_download_info.txt (in Git).
#
#   Rscript Scripts/prepare_apc.R
#       Never changes the version. Checks the pinned file; if it is missing,
#       downloads the current export and keeps it only if its md5 equals the
#       pinned one (otherwise stops: copy the pinned file from another machine).
#   Rscript Scripts/prepare_apc.R --update
#       Downloads the current export and compares its structure (columns and
#       the values of taxonomicStatus, taxonRank, nameType) with the pinned
#       version. Same structure: writes it under the server's file name (e.g.
#       APC-taxon-2026-11-28-1134.csv) with its tracking file, and sets
#       taxonomy.authority_file in the config (and in Scripts/config/ala_crabs.yml
#       if that pointed to the same file). Then commit the config(s) and the
#       tracking file and push; git pull on another machine and run this script
#       there (without --update) to fetch the same file. Different structure: stops
#       and lists the differences; nothing is changed.
#   --config <file>   config with the pinned file (default Scripts/config/ala_crabs_its_plants.yml)
#
# Usage (from the project root): Rscript Scripts/prepare_apc.R [--update] [--config <file>]

source("Scripts/workflow_helpers/helpers_authority_files.R")

url <- "https://biodiversity.org.au/nsl/services/export/taxonCsv"
out_dir <- "ALA_Records"
value_cols <- c("taxonomicStatus", "taxonRank", "nameType")  # values read_apc() depends on
active_config <- "Scripts/config/ala_crabs.yml"

args <- commandArgs(trailingOnly = TRUE)
update <- "--update" %in% args
config_file <- if ("--config" %in% args) args[which(args == "--config") + 1L] else "Scripts/config/ala_crabs_its_plants.yml"

pinned <- yaml::read_yaml(config_file)$taxonomy$authority_file
if (is.null(pinned) || !nzchar(pinned)) stop("no taxonomy.authority_file in ", config_file)
info <- read_download_info(pinned)
if (is.null(info) || is.na(info["md5"])) stop("no tracking file with an md5 line for ", pinned, ": version not pinned")

# Download the current export into a temporary file.
download_current <- function() {
  tmp <- tempfile(fileext = ".csv")
  message("downloading ", url)
  resp <- curl::curl_fetch_disk(url, tmp, handle = curl::new_handle(followlocation = TRUE))
  if (resp$status_code != 200) stop("download failed: HTTP ", resp$status_code)
  headers <- curl::parse_headers_list(resp$headers)
  type <- headers[["content-type"]] %||% ""
  if (!grepl("text/csv", type)) {
    stop("the server did not return a CSV (content type '", type, "'), probably the human-verification ",
         "page. Download 'APC Taxon' by hand from https://biodiversity.org.au/nsl/services/export/index ",
         "(see 'Name authority file' in README.md).")
  }
  fname <- sub('.*filename="?([^";]+)"?.*', "\\1", headers[["content-disposition"]] %||% "")
  if (!grepl("^APC-taxon-.*[.]csv$", fname)) fname <- paste0("APC-taxon-", format(Sys.Date(), "%Y-%m-%d"), ".csv")
  list(tmp = tmp, fname = fname, md5 = file_md5(tmp), time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
}

# Run as a function so that source() returns instead of ending the R session.
main <- function() {
if (!update) {
  if (file.exists(pinned)) {
    if (file_md5(pinned) == info["md5"]) { message("ok: ", pinned, " is the pinned version (md5 ", info["md5"], ")"); return(invisible(NULL)) }
    stop(pinned, " does not match its pinned md5 ", info["md5"], ". Restore the pinned file.")
  }
  d <- download_current()
  if (d$md5 != info["md5"]) {
    stop("the current APC export (", d$fname, ", md5 ", d$md5, ") is not the pinned version ", basename(pinned),
         " (md5 ", info["md5"], "); the server no longer provides it. Copy the pinned file from another ",
         "machine, or switch to the current export with --update.")
  }
  if (!file.copy(d$tmp, pinned)) stop("could not write ", pinned)
  message("restored pinned version: ", pinned)
  return(invisible(NULL))
}

# --update ---------------------------------------------------------------------
pinned_sig <- info[grepl("^(columns|values_)", names(info))]
if (!length(pinned_sig)) {
  if (!file.exists(pinned)) stop("no structure lines in the tracking file of ", pinned, " and the file is missing")
  pinned_sig <- structure_signature(pinned, value_cols)
  pinned_sig <- pinned_sig[names(pinned_sig) != "rows"]
}
d <- download_current()
if (d$md5 == info["md5"]) { message("the current export is the pinned version; nothing to update"); return(invisible(NULL)) }
new_sig <- structure_signature(d$tmp, value_cols)
diffs <- compare_signature(pinned_sig, new_sig)
if (length(diffs)) {
  stop("the current APC export (", d$fname, ") has a different structure than the pinned version ",
       basename(pinned), "; nothing changed. Differences:\n  ", paste(diffs, collapse = "\n  "),
       "\nCheck read_apc() against the new export before using it.")
}
if (as.numeric(new_sig["rows"]) < 100000) stop("APC export has only ", new_sig["rows"], " rows")

dest <- file.path(out_dir, d$fname)
if (file.exists(dest)) stop(dest, " exists already; not overwritten")
if (!file.copy(d$tmp, dest)) stop("could not write ", dest)
write_download_info(dest, c(url = url, downloaded = d$time, file = d$fname, md5 = d$md5, new_sig))
changed <- config_file
set_authority_file(config_file, dest)
if (file.exists(active_config) && normalizePath(active_config) != normalizePath(config_file) &&
    identical(yaml::read_yaml(active_config)$taxonomy$authority_file, pinned)) {
  set_authority_file(active_config, dest); changed <- c(changed, active_config)
}
message("same structure as ", basename(pinned), "; new version: ", dest, " (", new_sig["rows"], " rows)")
message("taxonomy.authority_file set to ", dest, " in: ", paste(changed, collapse = ", "))
message("commit and push: git add ", paste(c(changed, download_info_path(dest)), collapse = " "))
}

main()
