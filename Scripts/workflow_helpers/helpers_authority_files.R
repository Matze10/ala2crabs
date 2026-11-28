# Version pinning of the name authority CSVs (APC, AFD).
#
# Each authority CSV in ALA_Records/ has a tracking file
# <file>_download_info.txt, tracked in Git, with lines "key: value":
#   url, downloaded, file, md5, rows       source, time, checksum, row count
#   columns                                 column names, in file order
#   values_<column>                         sorted distinct values (joined by " | ")
#                                           of the columns the readers depend on
# (AFD from ALA: raw_file, raw_md5, raw_rows and raw_* structure lines for the
# raw ALA download as well.)
#
# The md5 pins the version: tar_make() stops if the file does not match it
# (check_pinned_md5(), used by authority_input() in helpers_synonyms.R).
# The structure lines let Scripts/prepare_apc.R and prepare_afd_from_ala.R
# accept a newer download (--update) only if it has the same structure.

download_info_path <- function(csv) sub("[.]csv$", "_download_info.txt", csv)

read_download_info <- function(csv) {
  f <- download_info_path(csv)
  if (!file.exists(f)) return(NULL)
  l <- readLines(f, warn = FALSE)
  l <- l[grepl(":", l, fixed = TRUE)]
  setNames(trimws(sub("^[^:]*:", "", l)), sub(":.*$", "", l))
}

write_download_info <- function(csv, info) {
  writeLines(paste0(names(info), ": ", unname(info)), download_info_path(csv))
}

file_md5 <- function(f) unname(tools::md5sum(f))

# Structure of a CSV: column names and the distinct values of value_cols.
structure_signature <- function(csv, value_cols, prefix = "") {
  dt <- data.table::fread(csv, colClasses = "character", na.strings = "", encoding = "UTF-8",
                          showProgress = FALSE)
  sig <- c(columns = paste(names(dt), collapse = ", "))
  for (v in value_cols) {
    vals <- if (v %in% names(dt)) sort(unique(stats::na.omit(dt[[v]]))) else "<column missing>"
    sig[paste0("values_", v)] <- paste(vals, collapse = " | ")
  }
  sig <- c(sig, rows = as.character(nrow(dt)))
  if (nzchar(prefix)) names(sig) <- paste0(prefix, names(sig))
  sig
}

# Differences between two signatures (rows are ignored: they may change).
compare_signature <- function(pinned, new) {
  keys <- setdiff(union(names(pinned), names(new)), grep("rows$", names(new), value = TRUE))
  bad <- keys[is.na(pinned[keys]) | is.na(new[keys]) | pinned[keys] != new[keys]]
  vapply(bad, function(k) sprintf("%s\n    pinned: %s\n    new:    %s", k, pinned[k], new[k]), "")
}

# tar_make() check: the authority file must match the md5 of its tracking file.
check_pinned_md5 <- function(path) {
  info <- read_download_info(path)
  if (is.null(info) || is.na(info["md5"])) {
    warning("No tracking file with an md5 line for ", basename(path), ": its version is not pinned.")
    return(invisible(NA))
  }
  md5 <- file_md5(path)
  if (!identical(md5, unname(info["md5"]))) {
    stop(basename(path), " does not match its tracking file ", basename(download_info_path(path)),
         " (md5 ", md5, ", pinned ", info["md5"], "). Restore the pinned file, or make a new version with ",
         "--update (see 'Name authority file' in README.md).")
  }
  invisible(TRUE)
}

# Set taxonomy.authority_file in a config, keeping its comments and layout.
set_authority_file <- function(config_file, new_path) {
  l <- readLines(config_file, warn = FALSE)
  i <- grep("^  authority_file:", l)
  if (length(i) != 1L) stop("expected one 'authority_file:' line in ", config_file)
  old <- trimws(sub("^  authority_file:", "", l[i]))
  l[i] <- paste0("  authority_file: ", new_path)
  writeLines(l, config_file)
  old
}
