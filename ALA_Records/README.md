# ALA_Records

Input data for the pipeline. Everything else (setup, configuration, review, outputs) is in [`../README.md`](../README.md).

This folder holds:

- the ALA occurrence download of each example as one ZIP (set `input.combined_zip` and `input.expected_members` in `Scripts/config/ala_crabs.yml`):
  - Example 1 (plants): `plant_ala_records-2026-11-21/records-2026-11-21.zip`;
  - Example 2 (arthropods): `insect_ala_records-2026-11-28/records-2026-11-28.zip` (ALA species group Arthropods, crustaceans included; they are removed by the config);
- the ALA metadata of each download (`citation.csv`, `doi.txt`, `headings.csv`, `README.html`) and, for Example 1, `search_settings_ala.txt`;
- the name authority CSV (set `taxonomy.authority_file`):
  - Example 1: `APC-taxon-2026-11-22-4420.csv`, downloaded by hand on 2026-11-22; `Scripts/prepare_apc.R` checks it against its tracking file, and `--update` gets a newer export;
  - Example 2: `AFD-ALA-arthropoda-20261128.csv`, built by `Scripts/prepare_afd_from_ala.R` from the ALA species index (raw download `AFD-ALA-arthropoda-20261128/bie_download.csv`);
  - `AFD-20261125T151017.csv`: AFD website export of Vespidae (used to check the reader; not used by the examples).

Each authority CSV has a tracking file `<file>_download_info.txt` (source URL, download time, MD5 checksum, row count, structure) that pins its version; see [Version pinning](../README.md#version-pinning-tracking-files).

The ZIPs, the authority CSVs and the raw species-index download are large and not tracked in Git; the metadata and the tracking files are. The example data can be downloaded from their DOIs, see [Example 1](../README.md#example-1-pollen-its-list-for-nsw-and-qld-seed-plants) and [Example 2](../README.md#example-2-coi-list-for-nsw-and-qld-arthropods-without-crustaceans).

The sections that pipeline messages and config comments refer to are now in the main README: [Name authority file](../README.md#name-authority-file-apc-or-afd) and [Barcode marker and the NCBI query](../README.md#marker-barcode-marker-and-the-ncbi-query).
