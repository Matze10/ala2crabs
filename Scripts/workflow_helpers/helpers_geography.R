# Region flags per species.
# The download covers all of NSW and QLD, but the study area may be smaller.
# Each region in geography.region_flags is a latitude/longitude box: a species
# gets the region's flag when it has at least one georeferenced record and ALL
# its georeferenced records lie inside the box (e.g. fnq_only: every record
# north of Cairns). Flagged species are KEPT in the CRABS export unless the
# region has exclude: true, because close relatives in the reference database
# protect against confident misassignment; the flag marks them so that a hit
# gets checked before it is trusted. Every target record counts, including
# cultivated/garden records.

# Normalised, validated regions from the config: a list of
# list(name, description, boxes, clear_states, exclude), where boxes is a list
# of list(lat_min, lat_max, lon_min, lon_max). An area gives its bounds either
# directly (one box) or as a boxes: list (one flag for several areas).
# A missing bound (left out or ~) means no limit on that side. latitude_min and
# longitude_min are strict (">"), the max bounds inclusive ("<="), so a record
# exactly at the FNQ boundary (-16.92) counts as outside, as before.
# Backward compatibility: without geography.region_flags, the old
# geography.fnq_boundary_latitude (default -16.92) gives the single region
# fnq_only; the deprecated crabs.exclude_fnq_only: true sets its exclude.
region_flag_settings <- function(config) {
  geo <- config$geography
  regions <- geo$region_flags
  if (is.null(regions)) {
    b <- if (is.null(geo$fnq_boundary_latitude)) -16.92 else as.numeric(geo$fnq_boundary_latitude)
    regions <- list(list(
      name = "fnq_only", description = "Far North Queensland only (all records north of Cairns)",
      latitude_min = b, clear_states = list("New South Wales"), exclude = FALSE
    ))
  }
  num <- function(x) if (is.null(x) || (length(x) == 1L && is.na(x))) NA_real_ else as.numeric(x)
  out <- lapply(regions, function(r) {
    name <- if (is.null(r$name)) "" else as.character(r$name)
    if (!nzchar(name)) stop("geography.region_flags: every region needs a name.")
    # The name becomes a column in the flag and audit tables.
    if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", name)) {
      stop("geography.region_flags: name '", name, "' must start with a letter and contain only ",
           "letters, digits and underscores (it becomes a column name), e.g. fnq_only. ",
           "Use description for display text.")
    }
    # Bounds are given either directly on the area (one box) or as a boxes:
    # list (several boxes, one flag; a record is inside if it is in any box).
    bound_keys <- c("latitude_min", "latitude_max", "longitude_min", "longitude_max")
    as_box <- function(b, where) {
      box <- list(lat_min = num(b$latitude_min), lat_max = num(b$latitude_max),
                  lon_min = num(b$longitude_min), lon_max = num(b$longitude_max))
      if (all(is.na(unlist(box)))) {
        stop("geography.region_flags: ", where, " has no latitude or longitude bound.")
      }
      box
    }
    has_direct <- any(bound_keys %in% names(r))
    has_boxes <- !is.null(r$boxes)
    if (has_direct && has_boxes) {
      stop("geography.region_flags: region '", name, "' has both direct bounds and boxes; use one or the other.")
    }
    if (!has_direct && !has_boxes) {
      stop("geography.region_flags: region '", name, "' has no latitude or longitude bound and no boxes.")
    }
    boxes <- if (has_boxes) {
      if (!length(r$boxes)) stop("geography.region_flags: region '", name, "' has an empty boxes list.")
      lapply(seq_along(r$boxes), function(i) as_box(r$boxes[[i]], paste0("region '", name, "', box ", i)))
    } else {
      list(as_box(r, paste0("region '", name, "'")))
    }
    list(
      name = name,
      description = if (is.null(r$description)) name else as.character(r$description),
      boxes = boxes,
      clear_states = as.character(unlist(r$clear_states)),
      exclude = isTRUE(r$exclude)
    )
  })
  region_names <- vapply(out, `[[`, "", "name")
  if (anyDuplicated(region_names)) {
    stop("geography.region_flags: duplicate region name(s): ",
         paste(unique(region_names[duplicated(region_names)]), collapse = ", "))
  }
  # Deprecated key, kept working: crabs.exclude_fnq_only: true -> exclude fnq_only.
  if (isTRUE(config$crabs$exclude_fnq_only)) {
    for (i in seq_along(out)) if (out[[i]]$name == "fnq_only") out[[i]]$exclude <- TRUE
  }
  out
}

# SQL: the record is georeferenced for this region: it has a latitude, plus a
# longitude if any of the region's boxes has a longitude bound.
region_georef_sql <- function(reg) {
  uses_lon <- any(vapply(reg$boxes, function(b) !is.na(b$lon_min) || !is.na(b$lon_max), logical(1)))
  if (uses_lon) "(g.lat IS NOT NULL AND g.lon IS NOT NULL)" else "g.lat IS NOT NULL"
}

# SQL: the record is georeferenced for this region and lies inside any of its boxes.
region_inside_sql <- function(reg) {
  box_sql <- vapply(reg$boxes, function(b) {
    parts <- c(
      if (!is.na(b$lat_min)) paste0("g.lat > ", b$lat_min),
      if (!is.na(b$lat_max)) paste0("g.lat <= ", b$lat_max),
      if (!is.na(b$lon_min)) paste0("g.lon > ", b$lon_min),
      if (!is.na(b$lon_max)) paste0("g.lon <= ", b$lon_max)
    )
    paste0("(", paste(parts, collapse = " AND "), ")")
  }, "")
  paste0("(", region_georef_sql(reg), " AND (", paste(box_sql, collapse = " OR "), "))")
}

# SQL: the record counts against the region: georeferenced but outside all its
# boxes, or in one of the clear_states (also without coordinates).
region_outside_sql <- function(reg, con) {
  outside <- paste0("(", region_georef_sql(reg), " AND NOT ", region_inside_sql(reg), ")")
  if (length(reg$clear_states)) {
    outside <- paste0("(", outside, " OR g.stateProvince IN (",
                      paste(DBI::dbQuoteString(con, reg$clear_states), collapse = ", "), "))")
  }
  outside
}

# One row per canonical binomial with record counts, coordinate extremes and,
# per region, n_inside_<name>, n_outside_<name> and the TRUE/FALSE column
# <name>. region_flags lists the regions a species is flagged for ("" if none).
# Records are matched to a species through any of their name fields, like the
# candidate-name extraction.
species_latitude_flags <- function(database, config, output_path, excluded_genera = read_excluded_genera(config)) {
  regions <- region_flag_settings(config)
  geo <- config$geography
  val <- function(x, default) if (is.null(x)) default else as.numeric(x)
  lat_min <- val(geo$valid_latitude_min, -45)
  lat_max <- val(geo$valid_latitude_max, -9)
  lon_min <- val(geo$valid_longitude_min, 112)
  lon_max <- val(geo$valid_longitude_max, 160)  # includes Lord Howe Island (~159 E)

  con <- DBI::dbConnect(RSQLite::SQLite(), database)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  rank_sql <- paste(DBI::dbQuoteString(con, as.character(unlist(config$target$accepted_ranks))), collapse = ", ")
  phylum_sql <- paste(DBI::dbQuoteString(con, as.character(unlist(config$target$exclude_phyla))), collapse = ", ")
  class_sql <- paste(DBI::dbQuoteString(con, as.character(unlist(config$target$exclude_classes))), collapse = ", ")

  # Same target filter as extract_name_candidates(). Coordinates outside the
  # plausible band (e.g. 0,0 placeholders) are treated as missing.
  coord_sql <- function(col, lo, hi) paste0(
    "CASE WHEN CAST(r.", col, " AS REAL) BETWEEN ", lo, " AND ", hi,
    " AND TRIM(CAST(r.", col, " AS TEXT)) <> '' THEN CAST(r.", col, " AS REAL) END"
  )
  # Same name fields as the candidate names (target.candidate_fields, helpers_names.R).
  fields <- candidate_fields(con, config)
  qf <- as.character(DBI::dbQuoteIdentifier(con, fields))
  DBI::dbExecute(con, "DROP TABLE IF EXISTS temp.geo_records")
  DBI::dbExecute(con, paste0(
    "CREATE TEMP TABLE geo_records AS SELECT r.rowid AS rid, ",
    paste0("r.", qf, collapse = ", "), ", r.genus, r.specificEpithet, ",
    coord_sql("decimalLatitude", lat_min, lat_max), " AS lat, ",
    coord_sql("decimalLongitude", lon_min, lon_max), " AS lon, ",
    "r.stateProvince FROM ala_occurrences_raw r WHERE r.taxonRank IN (", rank_sql, ")",
    kingdom_sql(con, config, prefix = "r."),  # taxon.ala_kingdom (default Plantae)
    " AND (r.phylum IS NULL OR TRIM(r.phylum) = '' OR r.phylum NOT IN (", phylum_sql, "))",
    " AND (r.class IS NULL OR TRIM(r.class) = '' OR r.class NOT IN (", class_sql, "))",
    seed_plant_sql(con, config$target$exclude_orders, excluded_genera, prefix = "r.")
  ))

  obs_sql <- paste0(
    paste0(ifelse(seq_along(qf) == 1L, "SELECT rid, TRIM(", "UNION SELECT rid, TRIM("), qf,
           ifelse(seq_along(qf) == 1L, ") AS name FROM geo_records ", ") FROM geo_records "), collapse = ""),
    "UNION SELECT rid, TRIM(genus || ' ' || specificEpithet) FROM geo_records ",
    "WHERE genus IS NOT NULL AND TRIM(genus) <> '' AND specificEpithet IS NOT NULL AND TRIM(specificEpithet) <> ''"
  )

  # Canonicalise each distinct name once in R, then join back in SQLite so that
  # every record is counted once per species.
  names <- DBI::dbGetQuery(con, paste0("SELECT DISTINCT name FROM (", obs_sql, ") WHERE name IS NOT NULL AND name <> ''"))$name
  parsed <- canonicalise_names(names)
  name_map <- data.table::data.table(name = names, species_binomial = parsed$canonical_binomial)
  name_map <- name_map[!is.na(species_binomial) & nzchar(species_binomial)]
  DBI::dbWriteTable(con, "name_binomial_map", name_map, temporary = TRUE, overwrite = TRUE)
  DBI::dbExecute(con, "CREATE INDEX temp.idx_name_map ON name_binomial_map (name)")

  region_sql <- vapply(regions, function(reg) paste0(
    "SUM(CASE WHEN ", region_inside_sql(reg), " THEN 1 ELSE 0 END) AS n_inside_", reg$name, ", ",
    "SUM(CASE WHEN ", region_outside_sql(reg, con), " THEN 1 ELSE 0 END) AS n_outside_", reg$name
  ), "")

  flags <- data.table::setDT(DBI::dbGetQuery(con, paste0(
    "WITH obs AS (", obs_sql, "), ",
    "rb AS (SELECT DISTINCT o.rid, m.species_binomial FROM obs o JOIN name_binomial_map m ON m.name = o.name) ",
    "SELECT rb.species_binomial, COUNT(*) AS n_records, ",
    "SUM(g.lat IS NOT NULL) AS n_georeferenced, ",
    "MIN(g.lat) AS southernmost_latitude, MAX(g.lat) AS northernmost_latitude, ",
    "MIN(g.lon) AS westernmost_longitude, MAX(g.lon) AS easternmost_longitude, ",
    paste(region_sql, collapse = ", "),
    " FROM rb JOIN geo_records g ON g.rid = rb.rid GROUP BY rb.species_binomial"
  )))

  for (reg in regions) {
    n_in <- flags[[paste0("n_inside_", reg$name)]]
    n_out <- flags[[paste0("n_outside_", reg$name)]]
    flags[, (reg$name) := flags$n_georeferenced > 0 & n_out == 0 & n_in > 0]
  }
  region_names <- vapply(regions, `[[`, "", "name")
  hit <- as.matrix(flags[, region_names, with = FALSE])
  flags[, region_flags := apply(hit, 1, function(x) paste(region_names[x %in% TRUE], collapse = ";"))]
  data.table::setorder(flags, species_binomial)

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(flags, output_path)
  output_path
}

# Region flag columns present in a flags table (from its n_inside_<name> columns).
region_columns <- function(flags) sub("^n_inside_", "", grep("^n_inside_", names(flags), value = TRUE))

# Add region_flags and each region's flag column to a table, matched on the
# binomial of `name_col`.
add_region_flags <- function(dt, flags, name_col) {
  dt <- data.table::copy(data.table::as.data.table(dt))
  cols <- c("region_flags", region_columns(flags))
  if (!nrow(dt)) {
    for (cc in cols) dt[, (cc) := NA]
    return(dt)
  }
  key <- canonicalise_names(as.character(dt[[name_col]]))$canonical_binomial
  dt[, .region_key := key]
  f <- flags[, c("species_binomial", cols), with = FALSE]
  out <- f[dt, on = .(species_binomial = .region_key)]
  data.table::setnames(out, "species_binomial", ".region_key")
  out[, .region_key := NULL]
  data.table::setcolorder(out, c(setdiff(names(out), cols), cols))
  out[]
}

# Read species_latitude_flags.csv; NULL if no path is given or the file is missing.
read_latitude_flags <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  data.table::fread(path)
}
