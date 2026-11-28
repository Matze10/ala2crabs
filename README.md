---
version: 1.0.4
modified: 2026-11-25
---
# ALA to CRABS: species lists for reference databases from Atlas of Living Australia records

Version 1.0.4 (2026-11-25). See [Version history](#version-history).

This pipeline turns an Atlas of Living Australia (ALA) occurrence download into a reviewed species list for a CRABS (Creating Reference databases for Amplicon-Based Sequencing) reference database. It filters the records to the taxon you configure, cleans and canonicalises the names, resolves synonyms against a national name authority and NCBI Taxonomy, sends only doubtful names to manual review, counts the NCBI sequences of your barcode marker for every species, and writes one species per line in the spelling CRABS needs. Taxon, marker, filters and region flags are all set in one YAML file.

It works for **Australian flora and fauna**: plants are resolved against the Australian Plant Census (APC, ALA kingdom Plantae, e.g. ITS, rbcL, matK), animals against the Australian Faunal Directory (AFD, ALA kingdom Animalia, e.g. COI). The same code runs for both; only the `taxon:`, `taxonomy:` and `marker:` settings change ([Name authority file](#name-authority-file-apc-or-afd)). The repository has two worked examples, each with its own config file: a pollen ITS list of NSW and Queensland seed plants ([Example 1](#example-1-pollen-its-list-for-nsw-and-qld-seed-plants), `Scripts/config/ala_crabs_its_plants.yml`) and a COI list of NSW and Queensland arthropods without crustaceans ([Example 2](#example-2-coi-list-for-nsw-and-qld-arthropods-without-crustaceans), `Scripts/config/ala_crabs_coi_arthropods.yml`). The pipeline always reads `Scripts/config/ala_crabs.yml`, so the example to run is copied to that name ([Config file](#config-file-the-active-config-must-be-named-ala_crabsyml)).

All commands are run from the project root, the folder that contains this README, `Scripts/` and `renv.lock`.

## Quick start

1. Restore the R packages:

   ```r
   renv::restore()
   ```

2. Choose the config. The pipeline reads only `Scripts/config/ala_crabs.yml`. Copy an example to that name, e.g. for the arthropod COI example:

   ```r
   file.copy("Scripts/config/ala_crabs_coi_arthropods.yml", "Scripts/config/ala_crabs.yml", overwrite = TRUE)
   ```

   See [Config file](#config-file-the-active-config-must-be-named-ala_crabsyml).
3. Put your ALA download into `ALA_Records/`: one ZIP with the occurrence CSVs, e.g. `ALA_Records/plant_ala_records-2026-11-21/records-2026-11-21.zip`. Set `input.combined_zip` to its path and list the CSV names in `input.expected_members` in `Scripts/config/ala_crabs.yml`. To reproduce an example, download its datasets from their DOIs ([Example 1](#example-1-pollen-its-list-for-nsw-and-qld-seed-plants), [Example 2](#example-2-coi-list-for-nsw-and-qld-arthropods-without-crustaceans)). The data are not in this repository.
4. Get the name authority file and set `taxonomy.authority_file`: `Scripts/prepare_apc.R` downloads the APC for plants, `Scripts/prepare_afd_from_ala.R` builds the AFD for arthropods. See [Name authority file](#name-authority-file-apc-or-afd). For a group without an Australian authority, set `taxonomy.authority: none`.
5. Check the `marker:` and `taxon:` blocks in the config ([Configuration](#configuration)).
6. Copy `Scripts/config/secrets.env.example` to `Scripts/config/secrets.env` and fill in `NCBI_EMAIL` and `NCBI_API_KEY`.
7. Run the pipeline. In R:

   ```r
   targets::tar_make(script = "Scripts/_targets.R")
   ```

   From PowerShell on Windows (see [Running from PowerShell](#running-from-powershell)):

   ```powershell
   & "C:/Program Files/R/R-4.6.2/bin/Rscript.exe" -e "targets::tar_make(script = 'Scripts/_targets.R')"
   ```

For the example data, the first run takes about 6 hours: roughly 2.5 hours of local steps (SQLite import, filtering, region flags, names) and about 3 hours for the NCBI checks of ~26,000 names. A full rerun with the NCBI results cached takes about 3.5 hours. Later runs rerun only the steps whose code, configuration or inputs changed. Any edit to `ala_crabs.yml` reruns almost everything, because every step reads it.

## How it works

The steps are defined in [`Scripts/_targets.R`](Scripts/_targets.R) with the `targets` package. Functions are in [`Scripts/workflow_helpers/`](Scripts/workflow_helpers/), all settings in [`Scripts/config/ala_crabs.yml`](Scripts/config/ala_crabs.yml), and the report page in [`Scripts/analysis/05_filtering_summary.Rmd`](Scripts/analysis/05_filtering_summary.Rmd). Target names are given in brackets.

```text
config
 -> archive and header checks
 -> SQLite import of the ALA records
 -> genus check against ALA name matching (e.g. ferns/lycophytes)
 -> target records (kingdom, ranks, excluded groups)
 -> region flags (e.g. FNQ-only)
 -> candidate names
 -> canonical names
 -> review queue (clean names included automatically)
 -> name authority (APC/AFD) + NCBI Taxonomy synonyms
 -> manual decisions
 -> genus expansion + NCBI phrase names
 -> NCBI marker counts
 -> CRABS species list + audit files + filtering report
```

1. **Configuration** (`config_file`, `config`). The YAML is read and the `marker:` block checked; the pipeline stops at once if the marker is missing.
2. **Archive checks** (`input_validation`, `input_manifest_file`, `header_validation`). The ZIP must exist and contain every file in `input.expected_members`; its size and MD5 go to `Scripts/output/audit/input_manifest.csv`, and all CSV headers must match.
3. **SQLite import** (`raw_database`). The CSVs are streamed out of the ZIP into `Scripts/data/derived/ala_raw.sqlite` (about 10 minutes for the example).
4. **Genus check** (`ala_genus_checks`, `excluded_genera`). Every genus is looked up in ALA name matching (`taxon.ala_kingdom`). Genera in an `exclude_orders` order or an `exclude_species_groups` group are excluded, plus the genera in `exclude_genera_file`. Audit: `Scripts/output/audit/ala_genus_checks.csv`.
5. **Target records** (`target_records`). Species-level and infraspecific records (`target.accepted_ranks`) of the configured kingdom, without the excluded phyla, classes, orders and genera, are written to `Scripts/data/derived/target_occurrences.csv`. Garden and cultivated records count.
6. **Region flags** (`latitude_flags`). Species whose records all lie inside a configured latitude/longitude area are flagged ([Region flags](#geography-region-flags)).
7. **Candidate names** (`candidates`). Every distinct name string from the ALA name fields, counted in SQLite.
8. **Canonical names** (`normalised_names`). Authorities are stripped and rank markers standardised (`helpers_canonical_names.R`), including names with invisible spaces, a trailing `*`, or a common name in brackets or before a dash.
9. **Review queue** (`review_queue`). Clean species and infraspecific names are proposed `include`; doubtful names (sp./phrase names, cf./aff./?, hybrids, a/b names, genus only) get no decision.
10. **Synonyms** (`authority_file`, `ncbi_taxdump`, `synonym_matches`). Every name is resolved against the name authority and then matched to a local copy of NCBI Taxonomy ([What the synonym step does](#what-the-synonym-step-does)).
11. **Manual decisions** (`validated_decisions`, `taxonomy_matches`, `reviewed_species`). `taxonomy_decisions.csv` is checked, the review sheet written, and the final decision per name applied: manual first, automatic otherwise ([Manual review](#manual-review)).
12. **Genus expansion** (`genus_expansion`). Included genus-level names add the phrase names NCBI lists and the authority's species of the genus in the configured states ([Genus expansion](#genus-expansion-and-ncbi-phrase-names)).
13. **NCBI marker counts** (`entrez_checks`). For every included species, the NCBI taxonomy ID and the number of nucleotide sequences for the marker ([`marker`](#marker-barcode-marker-and-the-ncbi-query)).
14. **Outputs** (`crabs_outputs`, `qc_file`, `filtering_report`). The CRABS list, audit files, QC summary and the HTML report ([Outputs](#outputs)).

## Configuration

All settings are in [`Scripts/config/ala_crabs.yml`](Scripts/config/ala_crabs.yml), in the blocks below (in file order). Steps that only make sense for plants are marked; for another group they are replaced or emptied ([Plant-specific steps](#plant-specific-steps)).

### Config file: the active config must be named `ala_crabs.yml`

**The pipeline reads only `Scripts/config/ala_crabs.yml`.** The file name is fixed in `Scripts/_targets.R` (target `config_file`), in the report `Scripts/analysis/05_filtering_summary.Rmd` and in `authority_input()`; a config under any other name is ignored. The examples are kept under their own names and copied to `ala_crabs.yml` to run them:

| File | What it is |
|---|---|
| `Scripts/config/ala_crabs.yml` | **The active config**: the one `tar_make()` runs. Currently a copy of Example 2. |
| [`Scripts/config/ala_crabs_its_plants.yml`](Scripts/config/ala_crabs_its_plants.yml) | [Example 1](#example-1-pollen-its-list-for-nsw-and-qld-seed-plants): ITS, NSW and QLD seed plants. |
| [`Scripts/config/ala_crabs_coi_arthropods.yml`](Scripts/config/ala_crabs_coi_arthropods.yml) | [Example 2](#example-2-coi-list-for-nsw-and-qld-arthropods-without-crustaceans): COI, NSW and QLD arthropods without crustaceans. |

To switch, copy the example over `ala_crabs.yml`:

```r
file.copy("Scripts/config/ala_crabs_its_plants.yml", "Scripts/config/ala_crabs.yml", overwrite = TRUE)      # Example 1
file.copy("Scripts/config/ala_crabs_coi_arthropods.yml", "Scripts/config/ala_crabs.yml", overwrite = TRUE)  # Example 2
```

If you change settings for a run, make the change in the example file as well (or save your own copy under a new name), so the active config can always be restored from a named file. The two examples write to different `output:` paths (Example 2 under `Scripts/data/derived/coi/` and `Scripts/output/coi/`), so switching does not overwrite the other example's results. There is one `targets` store (`_targets/`) for both: after a switch, `tar_make()` reruns every step for the new config.

### `project`

`name`, `version` and `target_rank` are project information only and are not read by the code. In particular, `project.target_rank` does not filter anything; the taxon filter is `taxon.ala_kingdom`.

### `input`

| Key | What it does |
|---|---|
| `combined_zip` | Path of the ZIP with the ALA occurrence CSVs. |
| `expected_members` | File names of the CSVs inside the ZIP. Members are found by file name, so a folder inside the ZIP does not matter; any names work. |

Reading the ZIP needs the `unzip` program on the PATH (on Windows e.g. from Rtools). Do not modify the ZIP; the pipeline writes all derived data elsewhere. For a new download, save it in a new `ALA_Records/records-<date>/` folder and update both keys.

### `marker`: Barcode marker and the NCBI query

The NCBI step is marker independent. The `marker:` block is required and has no default: if `name` or `entrez_query` is missing or empty, `tar_make()` stops at the first step.

```yaml
marker:
  name: ITS
  entrez_query: ITS[Title] OR internal transcribed spacer[All Fields]
```

For every species in `Scripts/data/derived/reviewed_species.txt` the pipeline runs two NCBI searches:

1. NCBI Taxonomy: `"<species>"[Scientific Name]` gives the taxonomy ID and NCBI name (independent of the marker and the taxon).
2. NCBI Nucleotide: `"<species>"[Organism] AND (<entrez_query>)`, plus `AND (<taxon.ncbi_organism_filter>)` if that is set, gives the number of sequences for the marker. For *Acacia aneura* with the ITS setting:

   ```text
   "Acacia aneura"[Organism] AND (ITS[Title] OR internal transcribed spacer[All Fields])
   ```

Results go to `Scripts/data/derived/entrez_marker_checks.csv` (columns `ncbi_taxonomy_ids`, `ncbi_scientific_names`, `marker`, `ncbi_marker_hits`, `entrez_status`). To use another marker, change both lines. `entrez_query` is plain NCBI search syntax without the `[Organism]` part, which the pipeline adds; test it on the NCBI Nucleotide website first.

| Marker | `entrez_query` |
|---|---|
| ITS (current) | `ITS[Title] OR internal transcribed spacer[All Fields]` |
| rbcL | `rbcL[Gene] OR rbcL[Title]` |
| matK | `matK[Gene] OR matK[Title]` |
| trnL | `trnL[Gene] OR trnL[Title]` |
| COI (untested) | `COI[Gene] OR COX1[Gene] OR cytochrome oxidase subunit 1[Title]` |

A new marker means one fresh NCBI run (about 3 hours for ~26,000 species); switching back to a marker used before reuses its cache ([NCBI caches](#ncbi-caches)).

### `taxon`

Which organisms the list is for. The defaults reproduce the plant setup.

```yaml
taxon:
  label: Plantae
  ala_kingdom: Plantae
  ncbi_organism_filter: ""
  ncbi_root_taxid: 35493
  seed_plant_taxid: 58024
```

| Key | Default | What it does |
|---|---|---|
| `label` | `Plantae` | Name of the group in log messages, e.g. `Entrez (ITS, Plantae): ...`. |
| `ala_kingdom` | `Plantae` | Keeps ALA records of this kingdom (records with a blank kingdom are kept too) in the record filter, the name extraction and the region flags, and restricts ALA name matching of genera to it. `""` turns the filter off. |
| `ncbi_organism_filter` | `""` | Extra term added to every NCBI nucleotide search, e.g. `Insecta[Organism]`. Empty means the query is exactly as above and the cache is reused; a non-empty filter is part of the cache key, so it gets its own cache files (one fresh NCBI run). |
| `ncbi_root_taxid` | `35493` | NCBI Taxonomy subtree searched by the synonym step: 35493 land plants, 33208 animals, 6656 arthropods. Replaces `taxonomy.ncbi_root_taxid`, which is still read if `taxon:` has no value. |
| `seed_plant_taxid` | `58024` | Subtree used only for the audit column `ncbi_seed_plant` in `synonym_matches.csv` (58024 Spermatophyta); `~` leaves the column empty. |

### `target`

| Key | What it does |
|---|---|
| `accepted_ranks` | ALA `taxonRank` values kept (species, subspecies, variety, form, nothovariety). |
| `exclude_phyla`, `exclude_classes` | Records in these ALA phyla/classes are dropped. *Plant-specific*: mosses, liverworts, hornworts, algae, lycophytes, ferns. |
| `exclude_orders` | Records in these orders are dropped, and genera in them are excluded by the genus check. *Plant-specific*: fern and lycophyte orders. ALA files every land plant under Charophyta / Equisetopsida, so the class list alone never removes ferns. |
| `exclude_genera_file` | Static genus list, one per line (`Scripts/data/manual/fern_lycophyte_genera.txt`); covers records with a blank order and ALA homonyms such as *Alsophila*, *Polypodium* and *Pteridium*. *Plant-specific.* |
| `exclude_species_groups` | Genera whose ALA species group contains one of these texts are excluded (default `Ferns and allies`). *Plant-specific.* |
| `candidate_fields` | ALA name fields that candidate names are taken from (default `raw_scientificName`, `scientificName`, `acceptedNameUsage`, `originalNameUsage`); also used for the region flags. Each must be a column of the download, otherwise the pipeline stops. Genus + epithet(s) rebuilt from the rank fields are always added. |

#### Plant-specific steps

For another group, set `exclude_phyla`, `exclude_classes`, `exclude_orders` and `exclude_species_groups` to `[]` (or to that group's exclusions) and `exclude_genera_file` to `""`. Set them to `[]` rather than deleting them: if `exclude_orders` or `exclude_species_groups` is left out entirely, the built-in fern defaults are used. The name authority (APC) and `taxon.seed_plant_taxid` are the other plant choices. ALA genus checks are cached per genus, and per kingdom when it is not Plantae, so plant results are not reused for another kingdom.

### `taxonomy`

| Key | What it does |
|---|---|
| `authority`, `authority_file` | Name authority: `APC`, `AFD` or `none`, and the path of its manually downloaded CSV ([Name authority file](#name-authority-file-apc-or-afd)). |
| `manual_decisions` | `Scripts/data/manual/taxonomy_decisions.csv` ([Manual review](#manual-review)). |
| `secrets_file` | `Scripts/config/secrets.env` with `NCBI_EMAIL` and `NCBI_API_KEY`. |
| `request_interval_seconds` | Pause between NCBI requests per worker. |
| `entrez_workers` | Optional; parallel NCBI workers (default 4 with an API key, 1 without). |
| `entrez_cache` | Folder of the NCBI cache. |
| `ala_namematching_url`, `ala_cache`, `ala_request_interval_seconds` | ALA name-matching service, its cache folder and pause. |
| `ncbi_taxdump_dir` | Folder of the NCBI Taxonomy dump (downloaded once). |
| `authority_states` | State codes, as written in the authority's distribution field, for the [genus expansion](#genus-expansion-and-ncbi-phrase-names) (default `[NSW, Qld, LHI]`; APC writes Lord Howe Island as `LHI`, separate from NSW). |
| `authority_state_qualifiers_excluded` | Distribution entries with these qualifiers do not count (default `[presumed extinct]`). Naturalised, doubtfully naturalised and other entries count. |

With an API key the NCBI step uses 4 parallel workers (about 2 names per second, within NCBI's limit of 10 requests per second); without a key it uses one worker (about 0.5 names per second). Quotes around the values in `secrets.env` are allowed. Failed lookups are not cached and are retried on the next run. `Scripts/config/secrets.env` is in `.gitignore`; never commit it.

### `geography`: region flags

Region flags mark species whose records all lie in one part of the range, so that a hit on such a species can be checked before it is trusted, or the species dropped.

```yaml
geography:
  valid_latitude_min: -45
  valid_latitude_max: -9
  valid_longitude_min: 112
  valid_longitude_max: 160   # includes Lord Howe Island (~159 E)
  region_flags:
    - name: fnq_only
      description: Far North Queensland only (all records north of Cairns)
      latitude_min: -16.92
      clear_states: [New South Wales]
      exclude: false
```

- `valid_*`: coordinates outside these bands (e.g. 0,0 placeholders) count as missing.
- Each entry in `region_flags` is one area with up to four bounds: `latitude_min`, `latitude_max`, `longitude_min`, `longitude_max`. Leave out (or set `~`) any bound you don't need, but each area needs at least one. `latitude_min` and `longitude_min` are strict (`>`), the max bounds inclusive (`<=`), so a record exactly at -16.92 counts as outside `fnq_only`.
- A species gets the flag when it has at least one georeferenced record and all its georeferenced records lie inside the area. An area with a longitude bound needs both coordinates.
- `clear_states`: records with one of these `stateProvince` values always count as outside, also without coordinates (for `fnq_only`, any NSW record clears the flag).
- `exclude: false` keeps flagged species in the CRABS list; `exclude: true` drops them. A species flagged by several areas is dropped if any of them has `exclude: true`.
- `name` becomes a column name: letters, digits and underscores only (`fnq_only`, not `FNQ-only`); put display text in `description`. The pipeline stops on an invalid or duplicate name, or an area without any bound.

More areas, each with only the bounds it needs (commented out in the YAML; remove the `#` to use them):

```yaml
    - name: cape_york_only
      latitude_min: -14
      exclude: true
    - name: west_of_divide_only      # longitude cut only
      longitude_max: 148
      exclude: false
    - name: some_box_only            # full bounding box
      latitude_min: -28
      latitude_max: -26
      longitude_min: 152
      longitude_max: 153.5
      exclude: true
```

A species can be flagged by several areas; the column `region_flags` then lists them, e.g. `fnq_only;cape_york_only`.

**Several boxes as one flag.** Instead of direct bounds, an area can have a `boxes:` list. Each box takes any subset of the four bounds (at least one), and a record is inside the area if it is inside any of its boxes. An area has either direct bounds or `boxes:`, not both.

```yaml
    - name: remote_islands_only       # several boxes = one flag
      description: Lord Howe or Torres Strait islands only
      boxes:
        - {latitude_min: -31.8, latitude_max: -31.4, longitude_min: 158.9, longitude_max: 159.3}
        - {latitude_min: -10.8, latitude_max: -9.0, longitude_min: 141.5, longitude_max: 144.0}
      exclude: true
```

Use `boxes:` for one flag and one column for species found only in some combination of places (here: only on Lord Howe and/or in the Torres Strait). Use separate areas for a separate flag per place; then a species recorded in both places gets neither flag, because its records are not all inside one area.

Outputs: `Scripts/data/derived/species_latitude_flags.csv` (every species with `n_records`, `n_georeferenced`, the extreme coordinates, and per area `n_inside_<name>`, `n_outside_<name>` and the TRUE/FALSE column `<name>`), `Scripts/output/audit/region_flagged_species.csv` (every flagged species, with `in_crabs_list`), and the columns `region_flags` plus one per area in `review_queue.csv` and `crabs_review_queue.csv`. Changing an area reruns the flag step (about 20 minutes for the example) and the steps after it. Without `region_flags`, the old `geography.fnq_boundary_latitude` (default -16.92) gives the single area `fnq_only`.

### `crabs`: CRABS list options

Changing these reruns only the CRABS export and the report.

| Key | Default | What it does |
|---|---|---|
| `require_marker_hit` | `false` | `true` writes only species with at least one NCBI sequence for the marker. |
| `exclude_hybrids` | `true` | Drops names whose rank parses as hybrid. |
| `keep_ncbi_hybrids` | `true` | Keeps named hybrids that NCBI lists as species (`Genus x epithet` or `x Genus epithet`, e.g. *Citrus* x *aurantiifolia*) under the NCBI spelling, even with `exclude_hybrids: true`. Hybrid formulas and hybrids without an NCBI match stay excluded. |
| `exclude_uncertain` | `true` | Drops names still uncertain or genus-only after canonicalisation. |
| `exclude_genus_only` | `true` | Genus-level decisions are not written; CRABS would need the plain genus name or `--query '"Acacia"[Organism] AND ...'`. |
| `exclude_infraspecific` | `false` | `true` drops subspecies and varieties. |
| `collapse_to_binomial` | `true` | One line per species: infraspecific names and author variants collapse to the binomial. |
| `include_decisions` | `[include]` | Decisions that count as included (for the NCBI check, the list and the genus expansion). |
| `expand_genus_from_authority` | `true` | Adds the authority's species of included genus-level names ([Genus expansion](#genus-expansion-and-ncbi-phrase-names)). |
| `expand_genus_without_distribution` | `false` | For an authority without distributions (the AFD checklist export): `true` adds every accepted species of the genus in the authority (Australia-wide), `false` skips the expansion. Ignored for APC. |
| `keep_ncbi_phrase_names` | `true` | Writes phrase names that NCBI Taxonomy lists (e.g. *Samadera* sp. Mt Nardi) under the exact NCBI name. |

The old `exclude_fnq_only` is deprecated; if set to `true` it still excludes the area `fnq_only`.

### Genus expansion and NCBI phrase names

Genus-level decisions (e.g. "Acacia sp. New England (J.B.Williams 97011)" decided as *Acacia*, rank genus) are not written as genera, because a genus line would download every sequence of the genus worldwide. The step `genus_expansion` recovers what is usable from them, using only the local files (no web queries):

1. **NCBI phrase names** (`keep_ncbi_phrase_names`). The phrase name is looked up in the NCBI Taxonomy dump below `taxon.ncbi_root_taxid`, first as written and then without the voucher in brackets and a trailing "<state> Herbarium". If NCBI lists it (e.g. *Samadera* sp. Mt Nardi, *Olearia* sp. Rhizomatica (I.R.Telford 11549)), the exact NCBI scientific name goes into the list. Included names that the synonym step matched to an NCBI phrase name are written the same way. These are the only lines of the list that are not "Genus epithet", apart from named hybrids.
2. **Authority species** (`expand_genus_from_authority`). The name authority's accepted species of the genus whose distribution lists one of `taxonomy.authority_states` (without an excluded qualifier) are added if they are not in the list yet. They are mapped to NCBI like the synonym step: accepted name, authority synonyms, Latin ending variants. This works the same for plants (APC) and fauna (AFD). The AFD checklist export has no distribution column, so for AFD the state filter cannot be applied: with `crabs.expand_genus_without_distribution: true` every accepted AFD species of the genus is added (the AFD covers Australian fauna, so this is Australia-wide, not NSW/Qld only); with `false` (default) the expansion is skipped with a message. With `authority: none` it is always skipped. If an AFD file has a `DISTRIBUTION` or `STATES` column, the state filter is used as for APC.

Expanded species have no ALA records of their own, so they get no region flags; they still pass the hybrid and `require_marker_hit` filters and get NCBI marker counts. Everything added is listed in `Scripts/output/audit/genus_expansion_species.csv` (`source` = `ncbi_phrase_name` or `authority_genus_expansion`, the genus, the genus-level names it came from, authority name, ID and distribution, CRABS name, NCBI taxid and `ncbi_match_via`). In the example, 19 authority species and 8 NCBI phrase names are added.

**`require_marker_hit`: `false` is recommended.** CRABS downloads sequences per name anyway: a species without a sequence contributes nothing, but keeping it makes the list a complete inventory of the target flora or fauna, and the count can still be used to filter later. Use `true` only if the list itself should show what the reference database can contain. In the example about half of the checked names (8,975 of 19,243) had an ITS sequence.

### `output`

Paths of all derived files and outputs, including the report folder (`report_dir`, default `Scripts/output/reports`). The code and the report take every path from this block. Example 1 keeps intermediates in `Scripts/data/derived/` (not tracked) and the final files in `Scripts/output/`; Example 2 uses `Scripts/data/derived/coi/` and `Scripts/output/coi/`. The paths in [Outputs](#outputs) are those of Example 1.

## Name authority file (APC or AFD)

The pipeline checks every name against a national name authority and replaces synonyms with the accepted name before matching to NCBI. The authority file is not part of the ALA download and is not fetched by `tar_make()`. It is made with a script, and its version is **pinned**: a new version changes the results, so it is only used after a deliberate update that is committed to Git.

| Authority | Script | Pinned file (example) |
|---|---|---|
| APC (plants, Example 1) | [`Scripts/prepare_apc.R`](Scripts/prepare_apc.R) | `ALA_Records/APC-taxon-2026-11-22-4420.csv` |
| AFD (arthropods without crustaceans, Example 2) | [`Scripts/prepare_afd_from_ala.R`](Scripts/prepare_afd_from_ala.R) | `ALA_Records/AFD-ALA-arthropoda-20261128.csv` |

| Taxon group | `taxonomy.authority` | Source | `taxon.ncbi_root_taxid` |
|---|---|---|---|
| Plants (ITS, rbcL, matK, trnL) | `APC` | Australian Plant Census | `35493` (land plants) |
| Animals (e.g. COI) | `AFD` | Australian Faunal Directory | `33208` (animals) or `6656` (arthropods) |
| Anything, without an Australian authority | `none` | ALA names go straight to NCBI | as above |

### Version pinning: tracking files

The version in use is the file in `taxonomy.authority_file` of the config. Next to it is its **tracking file** `<file>_download_info.txt`, which is tracked in Git (the large CSVs are not). It holds the source URL, the download time, the **MD5 checksum** and row count, and the **structure** of the file: its column names and the values of the columns the readers depend on (`values_taxonomicStatus`, `values_taxonRank`, ...; for the AFD file also of the raw ALA download, `raw_*`). The functions are in [`helpers_authority_files.R`](Scripts/workflow_helpers/helpers_authority_files.R).

- **`tar_make()`** checks the authority file against the `md5` of its tracking file and stops if they differ (`check_pinned_md5()`), so a run can only use the pinned version.
- **The scripts without options never change the version.** If the pinned file is there and matches, they report `ok`. If it is missing, they download the current file and keep it only if its checksum equals the pinned one (the AFD file is rebuilt from its raw download if that is still there). Both servers only provide the current version, so once it has changed upstream the pinned file cannot be downloaded again: copy it from the machine that made it.
- **`--update`** downloads the current file and compares its structure with the pinned version. If the structure is the same, it writes the new file under a new name with its tracking file and sets `taxonomy.authority_file` in the example config (and in `ala_crabs.yml` if that pointed to the same file). If the structure differs, it stops, lists the differences, and changes nothing; then check the reader (`read_apc()` / `read_afd()`) against the new file first.

To update the authority and share it:

```r
# on the machine that updates (example: APC)
#   Rscript Scripts/prepare_apc.R --update
# then commit the changed config(s) and the new tracking file, e.g.
#   git add Scripts/config/ala_crabs_its_plants.yml ALA_Records/APC-taxon-<date>-<n>_download_info.txt
#   git commit -m "APC version <date>"; git push
# on another machine
#   git pull
#   Rscript Scripts/prepare_apc.R      # fetches the same file if the server still has it, else copy it
```

An update reruns the name matching and everything after it.

### Plants: Australian Plant Census (APC)

`Rscript Scripts/prepare_apc.R` (or `source("Scripts/prepare_apc.R")`) checks or restores the pinned **APC Taxon** export (CSV, about 100 MB, 114,000 rows); `--update` gets the current one. It downloads <https://biodiversity.org.au/nsl/services/export/taxonCsv>, the same file as the "APC Taxon" link on <https://biodiversity.org.au/nsl/services/export/index>, and saves it under the server's file name, e.g. `APC-taxon-2026-11-28-1134.csv` (download date and time). If the site answers with its human-verification page instead of the CSV, the script stops; then download **APC Taxon** by hand from the export page, save it in `ALA_Records/` without renaming, and write its tracking file with `write_download_info()` (as was done for the Example 1 file, downloaded by hand on 2026-11-22). The APC version is taken from the file's `modified` column and written to `authority_version` (Example 1: APC 2026-11-17).

```yaml
taxonomy:
  authority: APC
  authority_file: ALA_Records/APC-taxon-2026-11-22-4420.csv
```

Structure compared on `--update`: column names and the values of `taxonomicStatus`, `taxonRank` and `nameType`.

### Animals: Australian Faunal Directory (AFD)

The AFD website exports at most **1,000 taxa per download** ("Bulk data limit exceeded" above that), so it only works for small groups. Arthropoda has about 120,000 taxa.

**Small group (up to 1,000 taxa), from the AFD website:** open the AFD page of the group, e.g. <https://biodiversity.org.au/afd/taxa/VESPIDAE>, download the group's checklist as CSV (a file like `AFD-20261125T151017.csv`), save it in `ALA_Records/` without renaming, and set `taxonomy.authority: AFD`, `taxonomy.authority_file` and `taxon.ncbi_root_taxid` (33208 animals, 6656 arthropods). Such a file has no tracking file unless you write one; `tar_make()` then warns that its version is not pinned.

**Arthropods without crustaceans (Example 2), from the ALA species index:** ALA builds its animal taxonomy from the current AFD; every taxon carries its AFD GUID, the AFD classification and its AFD synonyms. `Rscript Scripts/prepare_afd_from_ala.R` (or `source(...)`) checks or rebuilds the pinned file; `--update` gets the current index. It downloads
`https://bie-ws.ala.org.au/ws/download?q=*&fq=idxtype:TAXON&fq=rk_phylum:ARTHROPODA&fq=-rk_subphylum:CRUSTACEA&fq=rank:(species OR subspecies)`
into `ALA_Records/AFD-ALA-arthropoda-<YYYYMMDD>/bie_download.csv` and converts it to the AFD checklist format in `AFD-ALA-arthropoda-<YYYYMMDD>.csv`; the conversion is deterministic, so the same raw download always gives the same file (same md5). Only rows with `datasetName` `AFD` are kept: ALA's own additions (lists "Taxonomic Additions - ... 2025", 1,396 names, e.g. exotic biosecurity species) are dropped. Synonyms come from the ALA field `additionalNames`. The pinned download of 2026-11-28 gave 82,915 accepted AFD names (81,069 species, 1,846 subspecies; Insecta 72,486, Arachnida 8,598, the rest Diplopoda, Pycnogonida, Collembola, Chilopoda, Pauropoda, Entognatha, Symphyla) and 59,938 synonyms. All 239 names of the AFD website export for Vespidae (`AFD-20261125T151017.csv`) are in it with the same accepted names. Structure compared on `--update`: the raw download's column names and the values of `taxonRank` and `taxonomicStatus`.

Not suitable: the AFD copy on ChecklistBank (dataset 127379, version 2022-10-12) lacks whole groups, e.g. Vespidae.

The reader (`read_afd()` in `helpers_synonyms.R`) uses the columns `GENUS`, `SPECIES`, `SUBSPECIES` (canonical name without subgenus and authors), `NAME_TYPE` (`Valid Name` = accepted; other types, e.g. synonyms, are linked to the valid name), `RANK` (Species, Subspecies), `TAXON_GUID` (shared by a valid name and its synonyms; written as `authority_taxon_id`) and `NAME_LAST_UPDATE` / `TAXON_LAST_UPDATE` (version). The file is Latin-1 encoded. The checklist export has no distribution column, so the genus expansion cannot filter AFD species by state (see `crabs.expand_genus_without_distribution`). A checklist of a whole family can contain only valid names (the Vespidae example has 239 valid names and no synonyms); then the synonym step uses the AFD accepted names and the NCBI search order as for plants. AFD state names (e.g. "New South Wales") would be mapped to the APC codes (NSW, Qld, Vic, Tas, SA, WA, NT, ACT, LHI, NI, ChI) if a distribution column is present.

### What the synonym step does

For every canonical name, the authority's accepted name is used for synonyms and orthographic variants (misapplied, pro parte, doubtful and excluded entries are not followed). NCBI Taxonomy, below `taxon.ncbi_root_taxid`, is then searched for, in order: the accepted name, the ALA name, the other synonyms of the accepted taxon, and Latin gender-ending variants (`-us`/`-a`/`-um` ...). The NCBI scientific name of the first hit becomes the name in the CRABS list, because CRABS queries NCBI with `"<name>"[Organism]`. The route is recorded in `ncbi_match_via` (`authority_accepted`, `ala_name`, `authority_synonym`, `ending_variant` or `no_match`).

- `Scripts/data/derived/synonym_matches.csv`: one row per canonical name with `authority_status`, `authority_accepted_name`, `authority_taxon_id`, `ncbi_taxid`, `ncbi_name`, `ncbi_match_via` and `crabs_name`.
- `Scripts/output/audit/ncbi_unmatched.csv`: names with no NCBI match. They stay in the list under their authority or ALA name; CRABS will find no sequences for them.
- In `taxonomy_matches.csv`, `authority` / `authority_version` read e.g. `APC; NCBI Taxonomy` / `APC 2026-11-17; NCBI taxdump 2026-11-24`. Manual decisions keep their accepted name; only a blank or `Author` authority is filled.

Replacing the authority file, or changing `taxonomy.authority` or `taxon.ncbi_root_taxid`, reruns the matching and the steps after it. For a newer APC export, save it in `ALA_Records/` and change `taxonomy.authority_file`.

## Manual review

### What needs review

Names that parse to a clean species or infraspecific name are included automatically, so most names (~57,000 in the example) do not need a decision. Only doubtful names need one: sp./phrase names and sp. nov., cf./aff./?, hybrids, a/b names, complexes and intergrades, and names without a species epithet. A row in `taxonomy_decisions.csv` always overrides the automatic proposal, so a clean name can also be excluded or corrected there.

### The two files

- **`Scripts/output/audit/review_queue.csv`: the review sheet.** You edit it; the pipeline rewrites it on every run. It lists every review candidate, plus any other name that has a row in `taxonomy_decisions.csv`. The decision columns are always blank. The first nine columns are in the same order as `taxonomy_decisions.csv`; after them come the evidence columns `in_taxonomy_decisions`, `why_review`, `suggested_name`, `canonical_rank`, `occurrence_records` and the region-flag columns. Rows without a decision are sorted to the top. The file is UTF-8 with a BOM so Excel shows special characters correctly.
- **`Scripts/data/manual/taxonomy_decisions.csv`: the decisions.** The pipeline only reads it. One row per name.

`Scripts/data/derived/manual_review_queue.csv` is a plain list of the same doubtful names without decision columns. `Scripts/data/derived/taxonomy_review_queue.csv` is the full proposal table for all names; it is not meant for editing.

### Steps

1. Run to the review pause. This writes `review_queue.csv` and stops before NCBI:

   ```r
   targets::tar_make(names = validated_decisions, script = "Scripts/_targets.R")
   ```

2. Fill in the first nine columns of `review_queue.csv` for the rows you decide.
3. Copy those nine columns of the decided rows into `taxonomy_decisions.csv`.
4. Close both files in Excel (Windows blocks writing to an open file) and run `tar_make()`.

If `review_queue.csv` contains typed decisions when the pipeline rewrites it, a timestamped copy is first saved in `Scripts/data/manual/review_queue_backups/`.

### Two passes

1. **First run.** `tar_make()` runs completely with the automatically included names and any decisions already made. Doubtful names without a decision are not sent to NCBI.
2. **After review.** Add the decisions and run `tar_make()` again. Only the downstream steps rerun, and NCBI is queried only for newly included names.

### Decision file format

Allowed decisions: `include`, `exclude`, `needs_review`, `insufficient_evidence`. Columns of each row:

```text
normalised_name
accepted_name
accepted_rank
decision
reason
authority
authority_version
reviewer
review_date
```

Blank template rows are ignored, and `ok` is read as `include`. Decision matching tolerates the doubled quote marks that CSV round trips add to phrase names such as `Grevillea sp. "Baldy Range"`. Changing this file reruns the downstream steps but not the SQLite import. Keep old decisions separate if a new download changes the candidate names substantially.

## Outputs

A run writes many files, but you only need a few. Everything else is either an audit file you can open when you want to know why a name was kept or dropped, or an intermediate file the pipeline keeps so that a rerun can skip finished steps.

### 1. The result

| File | What it is |
|---|---|
| `Scripts/output/crabs/crabs_species_list.txt` | **The CRABS input**: one species per line (NCBI spelling). This is the file you give to CRABS. |
| `Scripts/output/reports/05_filtering_summary.html` | Report of the counts through the filters, the decisions and the region flags (`.md` if pandoc is not found). |

### 2. The file you edit

| File | What it is |
|---|---|
| `Scripts/output/audit/review_queue.csv` | Review sheet for doubtful names (rewritten every run, not tracked). Fill in the decision columns and copy them into `Scripts/data/manual/taxonomy_decisions.csv` ([Manual review](#manual-review)), which is tracked. |

### 3. Audit files (optional)

Open these to check or explain the list; nothing else needs them. The ones marked "yes" are tracked in Git, so the published repository shows what the example run kept and dropped.

| File | What it is | In Git |
|---|---|---|
| `Scripts/output/audit/qc_summary.csv` | The numbers of the run in one table (`section`, `metric`, `value`, `description`): records, names, decisions, authority and NCBI matches, genus expansion, NCBI marker check, CRABS list (binomials, hybrids, phrase names, lines with marker sequences) and region flags. | yes |
| `Scripts/output/audit/exclusions.csv` | Names with decision `exclude` or `insufficient_evidence`. | yes |
| `Scripts/output/audit/ncbi_unmatched.csv` | Names with no NCBI Taxonomy match (they cannot reach CRABS). | yes |
| `Scripts/output/audit/region_flagged_species.csv` | Species flagged by any area, with counts and whether each is in the CRABS list. | yes |
| `Scripts/output/audit/genus_expansion_species.csv` | Species added by the genus expansion and NCBI phrase names, with their source genus-level names. | yes |
| `Scripts/output/audit/crabs_review_queue.csv` | Names still without a decision, with their NCBI results. | no |
| `Scripts/output/audit/ala_genus_checks.csv` | ALA name-matching result per genus and why a genus was excluded. | no |
| `Scripts/output/audit/input_manifest.csv` | The input files with size and MD5 checksum. | no |

### 4. Intermediate files (not tracked, rebuilt by the pipeline)

`Scripts/data/derived/` holds the step-by-step files and caches. You do not need to open them, but `targets` uses them to skip steps that have not changed. The two large ones are `ala_raw.sqlite` (~7 GB) and `target_occurrences.csv` (~5 GB); see [Cleaning up disk space](#cleaning-up-disk-space) before deleting anything.

| File or folder | Step |
|---|---|
| `ala_raw.sqlite` | All ALA records imported from the ZIP. |
| `target_occurrences.csv` | Records left after the taxon, rank and region filters. |
| `species_latitude_flags.csv` | Per-species region flags. |
| `name_candidates.csv`, `name_candidates_normalised.csv` | Names taken from the records, then cleaned and canonicalised. |
| `taxonomy_review_queue.csv`, `manual_review_queue.csv` | All names with their automatic proposal; the second holds only the names that need review. |
| `synonym_matches.csv` | Per canonical name: authority status, accepted name, NCBI taxid and name, match route (`ncbi_match_via`), CRABS name. |
| `taxonomy_matches.csv` | Names after the manual decisions are applied. |
| `genus_expansion.csv` | Species added by the genus expansion (copied to the audit file above). |
| `reviewed_species.txt` | Species sent to the NCBI marker check. |
| `entrez_marker_checks.csv` | Per species: NCBI taxonomy IDs and number of marker sequences (`ncbi_marker_hits`). |
| `entrez_cache/`, `ncbi_taxdump/`, `ala_genus_cache/` | Caches of the NCBI and ALA web queries ([NCBI caches](#ncbi-caches)). |

### Confirm outputs

```r
readLines("Scripts/output/crabs/crabs_species_list.txt")            # one species per line
data.table::fread("Scripts/data/derived/synonym_matches.csv")      # authority and NCBI names
data.table::fread("Scripts/output/audit/ncbi_unmatched.csv")
data.table::fread("Scripts/data/derived/entrez_marker_checks.csv") # marker sequence counts
data.table::fread("Scripts/output/audit/review_queue.csv")         # review sheet
data.table::fread("Scripts/output/audit/crabs_review_queue.csv")
data.table::fread("Scripts/output/audit/exclusions.csv")
data.table::fread("Scripts/output/audit/region_flagged_species.csv")
data.table::fread("Scripts/output/audit/genus_expansion_species.csv")
data.table::fread("Scripts/output/audit/ala_genus_checks.csv")
data.table::fread("Scripts/output/audit/qc_summary.csv")
targets::tar_meta(fields = c(name, bytes, time, seconds, warnings, error))
```

The CRABS file has one canonical binomial per line. CRABS turns each line into an NCBI `"name"[Organism]` query without cleaning it, so the file must not contain authorities or `sp.` placeholders. The only exceptions are named hybrids (`Genus x epithet`) and phrase names that NCBI lists under exactly that name (e.g. `Samadera sp. Mt Nardi`), which CRABS can query as they are.

## Running the pipeline

### Rerun behaviour

`targets::tar_make(script = "Scripts/_targets.R")` runs only the steps invalidated by changed inputs, configuration or code. To run up to one step, e.g. only the report: `targets::tar_make(names = filtering_report, script = "Scripts/_targets.R")`. The dependency graph: `targets::tar_visnetwork(script = "Scripts/_targets.R")`.

To rerun everything, e.g. after replacing the ALA ZIP or changing the database rules:

```r
targets::tar_destroy(destroy = "all", script = "Scripts/_targets.R")
targets::tar_make(script = "Scripts/_targets.R")
```

The NCBI and ALA genus caches survive `tar_destroy()`, so names already checked are not queried again.

### Running from PowerShell

`Rscript` is often not on the PATH on Windows; on the development PC it is `C:/Program Files/R/R-4.6.2/bin/Rscript.exe`. PowerShell strips double quotes inside an argument, so use single quotes inside the R code:

```powershell
& "C:/Program Files/R/R-4.6.2/bin/Rscript.exe" -e "targets::tar_make(script = 'Scripts/_targets.R')"
```

Writing `-e 'targets::tar_make(script = "Scripts/_targets.R")'` fails in PowerShell with `unexpected symbol`.

### Report needs pandoc

The HTML report needs pandoc. pandoc ships with RStudio and Quarto but is not found when R runs from Rscript or Positron, so `render_report()` (in `helpers_io.R`) also looks in the RStudio, Quarto and Positron install folders. If none has pandoc, it writes `05_filtering_summary.md` instead, with a warning.

### NCBI caches

The caches live in `Scripts/data/derived/`, outside the `targets` store. When `targets` reruns `entrez_checks` (for example after any edit to the YAML), the step still reuses the cache, and only names without a cached result are sent to NCBI.

- **`entrez_cache/`** (about 7 MB for ~26,000 names): one `.rds` file per species and marker query. The file name is a hash of the species name and `marker.entrez_query`, plus `taxon.ncbi_organism_filter` if set (`entrez_cache_key()` in `helpers_entrez.R`). Only successful queries are saved, so failed ones are retried on the next run. Nothing expires automatically: to get fresh counts, e.g. to include sequences added to GenBank since, rename the folder before a run (e.g. to `entrez_cache_old`); the next run then queries all names again (about 3 hours).
- **`ncbi_taxdump/`** (about 530 MB): the NCBI Taxonomy dump (`taxdump.tar.gz`, about 80 MB download) used by the synonym step, downloaded on the first run. Its date is recorded in `ncbi_taxdump_date`. Delete the folder to fetch a newer dump.
- **`ala_genus_cache/`** (about 2 MB): ALA name matching per genus (and kingdom, if not Plantae).

### Cleaning up disk space

Sizes measured for the example after the 1.0.3 and 1.0.4 runs (November 2026).

| File or folder | Size | What it is for | Safe to delete? | Cost on the next run |
|---|---|---|---|---|
| `Scripts/data/derived/target_occurrences.csv` | 5.0 GB | Filtered records (result of `target_records`) | Yes | Filtering and all later steps rerun (about 2 h) |
| `Scripts/data/derived/ala_raw.sqlite` (+ `-wal`, `-shm`) | 4.5 GB + 0.9 GB (up to ~7 GB during an import) | SQLite copy of the ALA ZIP | Yes | Import (about 10 min) plus everything after it |
| `ALA_Records/plant_ala_records-2026-11-21/records-2026-11-21.zip` | 3.0 GB | The ALA download (input) | No, unless you can download it again | Needed for any rerun |
| `Scripts/data/derived/ncbi_taxdump/` | 530 MB | Local NCBI Taxonomy for the synonym step | Yes | Dump downloaded again (minutes) |
| `ALA_Records/APC-taxon-*.csv` | 96 MB | APC name authority (input) | No | Download again (`Scripts/prepare_apc.R`) |
| `Scripts/data/derived/entrez_cache/` | 7 MB | NCBI marker counts per species | Keep | All names queried again (about 3 h) |
| `Scripts/data/derived/ala_genus_cache/` | 2 MB | ALA name matching per genus | Keep | All genera checked again |
| `_targets/` | < 1 MB | `targets` metadata; results are stored as the files above | Only to force a full rerun | Everything reruns |

The large files can all be rebuilt from the ZIP and the authority file; deleting them after a final run is safe, but the next run then takes the full time. Not used by the current code, and safe to delete: `ALA_Records/ALA_terrestrial_inventory.sqlite` (+ `-wal`, `-shm`; 8.1 GB + 1.1 GB, an earlier database).

## Example 1: pollen ITS list for NSW and QLD seed plants

Config: [`Scripts/config/ala_crabs_its_plants.yml`](Scripts/config/ala_crabs_its_plants.yml) (copy it to `ala_crabs.yml` to run it). Purpose: a reference list of pollen-producing seed plants for ITS pollen metabarcoding in central QLD and SEQ. All NSW and QLD records are used, and garden and cultivated records count.

**Data.** Two ALA occurrence downloads made on 2026-11-21 (lifeform Plants, ALA quality filters as in [`search_settings_ala.txt`](ALA_Records/plant_ala_records-2026-11-21/search_settings_ala.txt)). They are not in the repository because they are about 19 GB uncompressed (3 GB zipped):

| State | DOI | Download date | CSV parts in the ZIP | Metadata folder |
|---|---|---|---|---|
| NSW | [10.26197/ala.9f02afb0-8210-4967-b28b-1fb8ce206601](https://doi.org/10.26197/ala.9f02afb0-8210-4967-b28b-1fb8ce206601) | 2026-11-21 | `records-2026-11-21_nsw.csv`, `records-2026-11-21_part2_nsw.csv` to `_part5_nsw.csv` | [`records-2026-11-21_nsw_metadata/`](ALA_Records/plant_ala_records-2026-11-21/records-2026-11-21_nsw_metadata/) |
| QLD | [10.26197/ala.06897d05-4b11-41b2-a84a-fdd2ff8d9f7c](https://doi.org/10.26197/ala.06897d05-4b11-41b2-a84a-fdd2ff8d9f7c) | 2026-11-21 | `records-2026-11-21_qld.csv`, `records-2026-11-21_part2_qld.csv` | [`records-2026-11-21_qld_metadata/`](ALA_Records/plant_ala_records-2026-11-21/records-2026-11-21_qld_metadata/) |

The CSVs of both downloads were renamed by hand to the part names above and zipped together as `ALA_Records/plant_ala_records-2026-11-21/records-2026-11-21.zip`. The DOIs give the exact 2026-11-21 snapshots; a new ALA search made later returns different records.

**Settings.** Marker ITS; `taxon` defaults (Plantae, NCBI root 35493); authority APC (`APC-taxon-2026-11-22-4420.csv`, APC version 2026-11-17); fern and lycophyte exclusion by phylum/class, order, genus list and ALA species group; area `fnq_only` (species recorded only north of Cairns, latitude > -16.92; any NSW record clears it), flagged but kept, because close relatives in the reference database reduce misassignment.

**Manual review conventions** (reviewer MR):

- Uncertain subspecies (`?ssp. x`, `subsp. ? x`): `Genus species x`, rank `subspecies`, include.
- Uncertain variety (`?var. x`), `Species ?`, `(?)`, `indet.`, `spp. agg.`: the binomial, rank `species`, include.
- `Species hybrid`: the parent species, rank `species`, include.
- sp. phrase names, sp. nov., `Genus subsp. x`: the plain genus name (`Acacia`, not `Acacia sp`), rank `genus`, include. Genus-level names are not written to the CRABS list as genera; the [genus expansion](#genus-expansion-and-ncbi-phrase-names) adds the NCBI phrase names and the APC species of those genera instead.
- a/b names, complexes, intergrades, `A ? B`, unusable labels: exclude, reason `rejected`.

**Results.** Final 1.0.4 run of 2026-11-25 (clean start, 2 h 24 min with the NCBI results of the stopped first attempt cached; a fully fresh run takes about 6 h): 12,119,773 raw records, 11,495,771 target records, 57,399 normalised names; **17,423 species in the CRABS list**, including 36 named hybrids and 11 phrase names that NCBI lists as taxa; 27 lines added by the genus expansion (19 APC species and 8 NCBI phrase names, see `genus_expansion_species.csv`); 3,949 names without an NCBI match; 1,199 species flagged `fnq_only` (kept); 8,987 of the 19,273 checked names with at least one ITS sequence at NCBI. No target ended with a warning or error. All numbers are in `Scripts/output/audit/qc_summary.csv`.

## Example 2: COI list for NSW and QLD arthropods without crustaceans

Config: [`Scripts/config/ala_crabs_coi_arthropods.yml`](Scripts/config/ala_crabs_coi_arthropods.yml) (copy it to `ala_crabs.yml` to run it). Purpose: a COI reference list of the NSW and QLD arthropods without crustaceans (insects, arachnids, myriapods, springtails, sea spiders and the other non-crustacean classes). Not yet run to a CRABS list.

**Data.** Two ALA occurrence downloads made on 2026-11-28 with ALA species group **Arthropods**, one per state, with the same filters: no fossil or eDNA records, no duplicates across datasets, no cultivated or captive records, no taxon-match problems, no absence records, no records needing identification, no records before 1700, no records with open user assertions or coordinate problems (0,0, swapped or negated, state or country centroid), spatially valid only, fewer than 3 environmental outlier layers, coordinate uncertainty up to 10 km. The full query is in each metadata folder's `README.html`.

| State | DOI | Download date | CSV in the ZIP | Metadata folder |
|---|---|---|---|---|
| NSW | [10.26197/ala.f8908ddd-bd16-4e91-9132-d6dfc692d837](https://doi.org/10.26197/ala.f8908ddd-bd16-4e91-9132-d6dfc692d837) | 2026-11-28 | `records-2026-11-28_nsw.csv` | [`records-2026-11-28_nsw_metadata/`](ALA_Records/insect_ala_records-2026-11-28/records-2026-11-28_nsw_metadata/) |
| QLD | [10.26197/ala.fd899d05-4447-4b23-bcc5-6cf8db5d0e3d](https://doi.org/10.26197/ala.fd899d05-4447-4b23-bcc5-6cf8db5d0e3d) | 2026-11-28 | `records-2026-11-28_qld.csv` | [`records-2026-11-28_qld_metadata/`](ALA_Records/insect_ala_records-2026-11-28/records-2026-11-28_qld_metadata/) |

Both CSVs are zipped together as `ALA_Records/insect_ala_records-2026-11-28/records-2026-11-28.zip` (4.5 GB uncompressed, 0.8 GB zipped). The download contains all arthropods; the crustaceans are removed by the config, not by the download. The download's classes (NSW + QLD records): Insecta 2,687,301, Arachnida 305,949, Malacostraca 140,961, Copepoda 22,688, Diplopoda 8,312, Chilopoda 6,815, Thecostraca 4,692, Collembola 4,490, Branchiopoda 4,420, Ostracoda 3,600, Pycnogonida 1,607, Symphyla 97, Ichthyostraca 71, Entognatha 22, Pauropoda 6, blank 981.

**Name authority.** AFD names of the non-crustacean arthropods built from the ALA species index on 2026-11-28 (`ALA_Records/AFD-ALA-arthropoda-20261128.csv`, [AFD](#animals-australian-faunal-directory-afd)).

**Settings.** Differences from Example 1:

```yaml
marker:
  name: COI
  entrez_query: COI[Gene] OR COX1[Gene] OR cytochrome oxidase subunit 1[Title]
taxon:
  label: Arthropoda
  ala_kingdom: Animalia
  ncbi_organism_filter: ""
  ncbi_root_taxid: 6656               # arthropods
  seed_plant_taxid: ~
target:
  accepted_ranks: [species, subspecies]
  exclude_phyla: []
  exclude_classes: [Malacostraca, Copepoda, Thecostraca, Branchiopoda, Ostracoda, Ichthyostraca,
                    Hexanauplia, Maxillopoda, Cephalocarida, Remipedia, Mystacocarida, Tantulocarida]
  exclude_orders: []
  exclude_genera_file: ""
  exclude_species_groups: []
taxonomy:
  manual_decisions: Scripts/data/manual/taxonomy_decisions_coi.csv
  authority: AFD
  authority_file: ALA_Records/AFD-ALA-arthropoda-20261128.csv
output:                               # all paths under Scripts/data/derived/coi/ and Scripts/output/coi/
```

- **Crustaceans** (subphylum Crustacea) are excluded by their classes: the first six occur in the download, the others are listed in case ALA uses them. The 981 records with a blank class are not removed by this filter; crustaceans among them show up in the review queue.
- **Manual decisions** go to their own file, `Scripts/data/manual/taxonomy_decisions_coi.csv`, so the plant decisions are not mixed in. It does not exist before the first review; the first run treats every name as undecided.
- **Region flag** `fnq_only` as in Example 1, flagged but kept.
- **Genus expansion:** the AFD file has no distributions, so `crabs.expand_genus_without_distribution: false` skips it (with `true`, genus-level names would add all AFD species of the genus, Australia-wide). NCBI phrase names work either way.
- The COI query has not been tested on the NCBI Nucleotide website yet. A new marker and taxon mean a fresh NCBI run and fresh ALA genus checks; the plant caches are not reused.

## Repository

Tracked in Git: the code in `Scripts/`, the configuration, the manual inputs (`Scripts/data/manual/taxonomy_decisions.csv`, `fern_lycophyte_genera.txt`), the READMEs, `LICENSE`, `.Rprofile`, `.here`, `renv.lock` with `renv/activate.R` and `renv/settings.json`, the ALA download metadata (citation, DOI, headings, search settings), the `*_download_info.txt` tracking files of the authority CSVs, and the small final outputs (`crabs_species_list.txt`, `exclusions.csv`, `ncbi_unmatched.csv`, `qc_summary.csv`, `region_flagged_species.csv`, `genus_expansion_species.csv`) of both examples (`Scripts/output/` and `Scripts/output/coi/`).

Not tracked (see [`.gitignore`](.gitignore)): the ALA data (cited by DOI), the APC/AFD authority CSVs and the raw ALA species-index download, the SQLite databases, the `_targets/` store, `Scripts/data/derived/` (intermediates and caches), the review backups, `Scripts/config/secrets.env` with the NCBI API key, and Excel lock files.

```text
ALA/
├── README.md                          (this file)
├── renv.lock
├── ALA_Records/                       (ALA downloads and the authority CSVs; see its README)
├── Scripts/
│   ├── _targets.R                     (pipeline definition)
│   ├── config/ala_crabs.yml           (the active config: the only one the pipeline reads)
│   ├── config/ala_crabs_its_plants.yml      (Example 1)
│   ├── config/ala_crabs_coi_arthropods.yml  (Example 2)
│   ├── config/secrets.env.example     (template for config/secrets.env)
│   ├── prepare_afd_from_ala.R         (builds the arthropod AFD authority file)
│   ├── workflow_helpers/              (helpers_*.R functions)
│   ├── analysis/05_filtering_summary.Rmd
│   ├── data/manual/                   (taxonomy_decisions*.csv, fern_lycophyte_genera.txt)
│   ├── data/derived/                  (intermediates and caches, not tracked; coi/ for Example 2)
│   └── output/                        (crabs/, audit/, reports/; coi/ for Example 2)
└── _targets/                          (targets store, not tracked)
```

## Licence and citation

This repository is licensed under CC BY 4.0 (see [`LICENSE`](LICENSE)). When using a list made with it, cite your ALA download by its DOI and the data providers in its `citation.csv`. For the example list, cite:

- Atlas of Living Australia occurrence download (NSW), 2026-11-21, <https://doi.org/10.26197/ala.9f02afb0-8210-4967-b28b-1fb8ce206601>, and (QLD), 2026-11-21, <https://doi.org/10.26197/ala.06897d05-4b11-41b2-a84a-fdd2ff8d9f7c>. The data are licensed by their providers, mostly CC-BY;
- the Australian Plant Census (APC), Council of Heads of Australasian Herbaria, with the APC version given in `authority_version` of `taxonomy_matches.csv`;
- NCBI Taxonomy and GenBank (NCBI), with the taxonomy dump date given in `ncbi_taxdump_date`.

## Version history

- **Unreleased**: Example 2 (COI, NSW and QLD arthropods without crustaceans).
  - Two example configs, `ala_crabs_its_plants.yml` (Example 1) and `ala_crabs_coi_arthropods.yml` (Example 2); the active config must be `ala_crabs.yml`.
  - `Scripts/prepare_afd_from_ala.R`: arthropod AFD names and synonyms from the ALA species index (the AFD website exports at most 1,000 taxa). `Scripts/prepare_apc.R`: the APC Taxon export by script (no manual download needed).
  - Authority versions pinned by tracking files `<file>_download_info.txt` in Git (md5 and structure); `tar_make()` stops if the authority file does not match; `--update` takes a newer version only if its structure is unchanged (`helpers_authority_files.R`).
  - Report folder configurable (`output.report_dir`); the report reads all paths from `output:` instead of fixed paths.
  - `.gitignore` tracks the final outputs of Example 2 under `Scripts/output/coi/`.
- **1.0.4 (2026-11-25)**: genus expansion, Australian fauna (AFD) support, config cleanup; first published version.
  - Genus expansion for plants and fauna (`crabs.expand_genus_from_authority`, `crabs.expand_genus_without_distribution`, `taxonomy.authority_states`, `taxonomy.authority_state_qualifiers_excluded`) and NCBI phrase names (`crabs.keep_ncbi_phrase_names`) for included genus-level names; audit `genus_expansion_species.csv`.
  - `target.candidate_fields` and `crabs.include_decisions` now used by the code; unused keys removed (`taxonomy.primary_source`, `taxonomy.primary_source_version`, `target.retain_ambiguous_names`, `target.retain_infraspecific`).
  - Seven autonyms in `taxonomy_decisions.csv` (e.g. *Acacia victoriae victoriae*) given rank `species`, matching their name; *Acaena* sp. Thredbo River Gorge decided as *Acaena* (was *Acacia*); `target_pha_count` removed from `qc_summary.csv`.
  - Removed: the unused analysis pages `01_`–`04_*.Rmd`, `_workflowr.yml`, the old template `Scripts/secrets.env.example` (the template is now `Scripts/config/secrets.env.example`) and the pre-1.0.2 file `entrez_its_checks.csv`.
  - AFD reader (`read_afd()`) for Australian fauna, tested on an AFD Vespidae checklist; README states the flora/fauna setups.
  - README: outputs grouped into result, review file, audit files and intermediate files.
  - `qc_summary.csv` extended from 5 to about 30 metrics with sections and descriptions (`build_qc_summary()`), including species counts for the CRABS list.
- **1.0.3 (2026-11-25)**: synonym step, configurable taxon and region flags.
  - Names resolved against APC (plants) or AFD (animals; reader tested on a Vespidae checklist) from a manually downloaded CSV (`taxonomy.authority`, `taxonomy.authority_file`), then matched to the NCBI Taxonomy dump below `taxon.ncbi_root_taxid`; the CRABS list uses the NCBI name. New files `synonym_matches.csv` and `ncbi_unmatched.csv`.
  - `authority` / `authority_version` filled from the sources instead of the placeholder.
  - Compared with 1.0.2, almost all names that left the list were renamed to their accepted name, which is in the list; the rest were hybrids and phrase names.
  - Named NCBI hybrids kept (`crabs.keep_ncbi_hybrids`).
  - User-set `taxon:` block (`label`, `ala_kingdom`, `ncbi_organism_filter`, `ncbi_root_taxid`, `seed_plant_taxid`) replaces the hard-coded kingdom Plantae and NCBI roots; ALA species groups for the genus check configurable (`target.exclude_species_groups`). Defaults reproduce the plant setup and reuse the NCBI cache.
  - Configurable region flags by latitude/longitude (`geography.region_flags`) replace the fixed FNQ flag, including multi-box areas (`boxes:`); `fnq_only` is the configured area, and `fnq_only_species.csv` became `region_flagged_species.csv`. Deprecated but still read: `geography.fnq_boundary_latitude` (used only if `region_flags` is absent) and `crabs.exclude_fnq_only` (`true` excludes the area `fnq_only`).
  - Filtering report: looks for pandoc in the RStudio, Quarto and Positron folders and writes Markdown if none is found; counts the occurrence file in chunks instead of loading it.
  - One general README; `ALA_Records/README.md` is now a short note on that folder.
- **1.0.2 (2026-11-22)**: first working end-to-end version.
  - Fern and lycophyte filter by order and genus, with an ALA name-matching check of every genus.
  - Canonical name parser (authorities stripped, rank markers standardised; handles invisible spaces, trailing `*`, common names in brackets or before a dash). CRABS list collapsed to binomials.
  - FNQ-only flag (species kept, flagged).
  - Clean names included automatically; manual review only for doubtful names, via `review_queue.csv` (blank decision columns, same order as `taxonomy_decisions.csv`) with backups before each rewrite. Decision matching tolerates doubled quote marks.
  - NCBI check: marker-independent (`marker:` block, required with no default; `crabs.require_marker_hit`), API key read correctly from quoted values, failed lookups retried instead of cached, 4 parallel workers.
