# ALA-to-CRABS targets pipeline, version 1.0.3
#
# Builds a reviewed species list from an ALA occurrence download for a CRABS
# reference database (taxon and barcode marker set in the config):
#
#   config -> archive checks -> SQLite import -> ALA genus check ->
#   target records -> region flags -> candidate names -> canonical names ->
#   review queue (clean names auto-included) -> APC/AFD + NCBI names ->
#   manual decisions -> genus expansion / NCBI phrase names ->
#   NCBI taxonomy + marker check -> CRABS species list and audit files
#
# Run from the project root (the ALA folder):
#   targets::tar_make(script = "Scripts/_targets.R")
# Stop at the manual-review pause (writes the review sheet output.review_queue_export):
#   targets::tar_make(names = validated_decisions, script = "Scripts/_targets.R")
#
# All settings are in Scripts/config/ala_crabs.yml; output paths are under
# config$output. See README.md for the full workflow.

library(targets)

source("Scripts/workflow_helpers/helpers_io.R")              # config, ZIP archive access, input checks
source("Scripts/workflow_helpers/helpers_sqlite.R")          # SQLite import, seed-plant record filter, QC
source("Scripts/workflow_helpers/helpers_canonical_names.R") # name parser: authorities, rank markers, binomials
source("Scripts/workflow_helpers/helpers_names.R")           # candidate-name extraction and review flags
source("Scripts/workflow_helpers/helpers_taxonomy.R")        # review queue, automatic proposals, applying decisions
source("Scripts/workflow_helpers/helpers_synonyms.R")        # APC/AFD accepted names and NCBI names per canonical name
source("Scripts/workflow_helpers/helpers_authority_files.R") # authority CSV version pinning (tracking files)
source("Scripts/workflow_helpers/helpers_validation.R")      # decision-file checks, review sheet
source("Scripts/workflow_helpers/helpers_entrez.R")          # NCBI taxonomy and marker sequence counts (cached)
source("Scripts/workflow_helpers/helpers_ala.R")             # ALA name-matching genus check for ferns/lycophytes
source("Scripts/workflow_helpers/helpers_geography.R")       # region flags by latitude/longitude (e.g. FNQ-only)
source("Scripts/workflow_helpers/helpers_expansion.R")       # genus expansion from the authority, NCBI phrase names
source("Scripts/workflow_helpers/helpers_crabs.R")           # CRABS species list and audit exports

# Packages used by target expressions and helper functions.
tar_option_set(
  packages = c("data.table", "DBI", "RSQLite", "here", "yaml", "digest", "rentrez", "httr", "jsonlite"),
  format = "rds"
)

list(
  # ---- Configuration and inputs -------------------------------------------
  # The YAML is tracked as a file, so any edit reruns everything downstream.
  # read_config() also stops at once if the barcode marker is not configured.
  tar_target(config_file, "Scripts/config/ala_crabs.yml", format = "file"),
  tar_target(config, read_config(config_file)),
  # Static fern/lycophyte genus list (manual input, tracked in Git).
  tar_target(excluded_genera_file, "Scripts/data/manual/fern_lycophyte_genera.txt", format = "file"),
  # The combined ALA ZIP exists and contains the expected CSV members.
  tar_target(input_validation, validate_inputs(config)),
  # Archive path, size and MD5 written to the audit folder.
  tar_target(
    input_manifest_file,
    {
      path <- here::here(config$output$input_manifest)
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      data.table::fwrite(input_validation, path)
      path
    },
    format = "file"
  ),
  # All CSV parts must have identical headers before they are combined.
  tar_target(header_validation, validate_headers(config)),
  tar_target(selected_columns, selected_ala_columns(header_validation$reference_header)),

  # ---- Raw import (slow, ~10 min) -----------------------------------------
  # Every CSV part is streamed into one SQLite table (ala_occurrences_raw).
  tar_target(
    raw_database,
    {
      header_validation
      import_ala_sqlite(config, here::here(config$output$database), selected_columns)
    },
    format = "file"
  ),

  # ---- Fern and lycophyte exclusion ---------------------------------------
  # ALA places all land plants under Charophyta / Equisetopsida, so ferns are
  # found by genus: every genus in the raw data is checked against ALA name
  # matching (order, "Ferns and allies" group; cached per genus) and combined
  # with the static genus list.
  tar_target(raw_genus_list, raw_genera(raw_database)),
  tar_target(
    ala_genus_checks,
    run_ala_genus_checks(raw_genus_list, config, here::here(config$output$ala_genus_checks)),
    format = "file"
  ),
  tar_target(
    excluded_genera,
    {
      excluded_genera_file
      excluded_genera_from_checks(ala_genus_checks, config)
    }
  ),

  # ---- Seed-plant records --------------------------------------------------
  # Species-level and infraspecific records, without excluded phyla, classes,
  # fern/lycophyte orders and genera. Garden/cultivated records are kept.
  tar_target(
    target_records,
    {
      filter_target_records(
        raw_database,
        config$target$accepted_ranks,
        config$target$exclude_phyla,
        config$target$exclude_classes,
        here::here(config$output$target_records),
        exclude_orders = config$target$exclude_orders,
        exclude_genera = excluded_genera,
        config = config
      )
    },
    format = "file"
  ),
  tar_target(raw_qc, qc_summary(raw_database, target_records)),
  # Region flags (geography.region_flags): a species is flagged for a region
  # when all its georeferenced records lie inside the region's lat/lon box
  # (e.g. fnq_only). Flagged species are kept unless the region has exclude: true.
  tar_target(
    latitude_flags,
    species_latitude_flags(raw_database, config, here::here(config$output$latitude_flags), excluded_genera),
    format = "file"
  ),

  # ---- Names ---------------------------------------------------------------
  # One row per distinct name string from all ALA name fields, aggregated in SQLite.
  tar_target(
    candidates,
    extract_name_candidates(raw_database, config, here::here(config$output$candidates), excluded_genera),
    format = "file"
  ),
  # Canonical name, binomial and rank per candidate, plus review flags
  # (uncertain, hybrid, infraspecific, authority text).
  tar_target(
    normalised_names,
    normalise_name_candidates(candidates, here::here(config$output$normalised_names)),
    format = "file"
  ),
  # Proposal for every name: clean names -> include automatically; doubtful
  # names -> no decision (also written to manual_review_queue.csv).
  tar_target(
    review_queue,
    build_review_queue(normalised_names, here::here(config$output$review_queue)),
    format = "file"
  ),

  # ---- Name authority and NCBI names ----------------------------------------
  # The authority CSV (APC for plants, AFD for animals) is made beforehand with
  # Scripts/prepare_apc.R or Scripts/prepare_afd_from_ala.R,
  # path in taxonomy.authority_file (see README.md). The NCBI
  # taxonomy dump is downloaded once. Every canonical name gets its accepted
  # name and the NCBI name used in the CRABS list; no NCBI match -> audit file.
  tar_target(authority_file, authority_input(config), format = "file"),
  tar_target(ncbi_taxdump, fetch_ncbi_taxdump(config), format = "file"),
  tar_target(
    synonym_matches,
    match_names(review_queue, authority_file, ncbi_taxdump, config,
                here::here(config$output$synonym_matches),
                here::here(config$output$ncbi_unmatched)),
    format = "file"
  ),

  # ---- Manual review (the pause) ------------------------------------------
  # Checks the decision file (taxonomy.manual_decisions) against the queue and writes the review sheet
  # output.review_queue_export (all review candidates, blank decision
  # columns; a backup is kept first if it contains typed decisions).
  tar_target(
    validated_decisions,
    validate_manual_decisions(
      review_queue,
      here::here(config$taxonomy$manual_decisions),
      here::here(config$output$review_queue_export),
      latitude_flags
    )
  ),
  # Final decision per name: manual decision if present, automatic proposal otherwise.
  tar_target(
    taxonomy_matches,
    {
      validated_decisions
      apply_manual_decisions(review_queue, config, here::here(config$output$taxonomy_matches), synonym_matches[1])
    },
    format = "file"
  ),
  # Included genus-level names: authority species of the genus in the configured
  # states (crabs.expand_genus_from_authority) and phrase names NCBI knows
  # (crabs.keep_ncbi_phrase_names), mapped to NCBI with the local taxdump.
  tar_target(
    genus_expansion,
    build_genus_expansion(taxonomy_matches, authority_file, ncbi_taxdump, config,
                          here::here(config$output$genus_expansion),
                          here::here(config$output$genus_expansion_audit)),
    format = "file"
  ),
  # Included species and infraspecific names (canonical form) plus the genus
  # expansion, for NCBI.
  tar_target(
    reviewed_species,
    export_reviewed_species(taxonomy_matches, here::here(config$output$reviewed_species), config, genus_expansion[1]),
    format = "file"
  ),

  # ---- NCBI (slow on first run, cached per species and marker) ------------
  # NCBI Taxonomy ID and the number of nucleotide sequences for the marker:
  # "<species>"[Organism] AND (<marker.entrez_query>).
  tar_target(
    entrez_checks,
    run_entrez_checks(reviewed_species, config, here::here(config$output$entrez_checks)),
    format = "file"
  ),

  # ---- Outputs -------------------------------------------------------------
  # CRABS species list (one binomial per line), exclusions, CRABS review queue
  # and the region-flagged species audit list.
  tar_target(
    crabs_outputs,
    build_crabs_outputs(
      taxonomy_matches,
      entrez_checks,
      config,
      here::here(config$output$crabs_species_list),
      here::here(config$output$exclusions),
      here::here(config$output$crabs_review_queue),
      latitude_flags,
      here::here(config$output$region_flagged_species),
      genus_expansion[1]
    ),
    format = "file"
  ),
  # HTML summary of record and name counts through the filters. Depends on the
  # CRABS outputs so it is re-rendered whenever they change. Needs pandoc
  # (bundled with RStudio/Quarto); without it a Markdown report is written.
  tar_target(
    filtering_report,
    {
      crabs_outputs
      render_report("Scripts/analysis/05_filtering_summary.Rmd",
                    here::here(config$output$report_dir %||% "Scripts/output/reports"))
    },
    format = "file"
  ),
  # QC summary (output.qc): record, name, decision,
  # authority, NCBI, genus-expansion, CRABS-list and region counts, one row per
  # metric with a short description (build_qc_summary() in helpers_io.R).
  tar_target(
    qc_file,
    {
      crabs_outputs
      build_qc_summary(raw_qc, target_records, candidates, normalised_names, taxonomy_matches,
                       synonym_matches[1], genus_expansion[1], entrez_checks, config,
                       here::here(config$output$qc))
    },
    format = "file"
  )
)
