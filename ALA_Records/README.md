# ALA_Records

Input data for the pipeline. Everything else (setup, configuration, review, outputs) is in [`../README.md`](../README.md).

This folder holds:

- the ALA occurrence download as one ZIP, e.g. `plant_ala_records-2026-11-21/records-2026-11-21.zip` (set `input.combined_zip` and `input.expected_members` in `Scripts/config/ala_crabs.yml`);
- the ALA metadata of each download (`citation.csv`, `doi.txt`, `headings.csv`, `README.html`) and `search_settings_ala.txt`;
- the manually downloaded name authority CSV, e.g. `APC-taxon-2026-11-22-4420.csv` (set `taxonomy.authority_file`).

The ZIP and the authority CSV are large and not tracked in Git; the metadata is. The example data can be downloaded from their DOIs, see [Example 1](../README.md#example-1-pollen-its-list-for-nsw-and-qld-seed-plants).

The sections that pipeline messages and config comments refer to are now in the main README: [Name authority file](../README.md#name-authority-file-apc-or-afd-manual-download-required) and [Barcode marker and the NCBI query](../README.md#marker-barcode-marker-and-the-ncbi-query).
