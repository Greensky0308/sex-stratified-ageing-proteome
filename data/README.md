# Input data

None of the input files are redistributed with this code. All are public or
controlled-access summary statistics and should be obtained from the sources below, then
placed at the paths listed. Files that require registration or an access agreement are
marked as controlled access.

## `data/instruments_ukbppp_grch37.tsv`

Top cis-pQTL per protein from the UK Biobank Pharma Proteomics Project (Olink), one row per
protein. Required columns: `protein`, `rsid`, `chromosome`, `position`, `effect_allele`,
`beta`. The GRCh38 version is produced by `liftover_instruments.py`.

## `data/sumstats/<ACCESSION>.h.tsv.gz`

Harmonised GWAS summary statistics, one file per accession (a matching `.tbi` enables the
ranged reads used by some scripts).

| trait group | accession range / IDs |
|---|---|
| brain imaging-derived phenotypes, sex-stratified | GWAS Catalog, GCST90429852 – GCST90432054 |
| appendicular lean mass | GCST90000027 (women), GCST90000026 (men) |
| grip-based muscle weakness | GCST90007527 (women), GCST90007528 (men) |

## `data/neale_panel/`

Body-composition and related GWAS summary statistics from UK Biobank round 2 (Neale lab),
file names `<field>_<sex>.tsv.bgz`, position-sorted BGZF without rsIDs. Fields used
include 23101 (whole-body fat-free mass), 23102 (whole-body water mass), 23105 (basal
metabolic rate), 23113 (leg fat-free mass), 23125 (arm fat-free mass), together with the
specificity controls listed in the analysis scripts.

## `data/decode/`

SomaScan summary data for individual proteins, one file per protein (downloads from the
deCODE summary-data portal). Used for the cross-platform checks.

## `data/kp4cd/`

External lean-mass summary statistics for the replication cohort, downloaded from the
Knowledge Portal Network (nodes for women and men).

## `data/igf1/`

Summary statistics for circulating IGF-1, used for the colocalization and mediation steps.

## `data/*.json` and `data/BrainIDP_SupplementaryData22.xlsx`

Supporting tables: brain phenotype list and sex-matched phenotype pairs
(`brain_idp_list.json`, `brain_pairs.json`), and gene coordinates for the instrument set
(`gene_coords_ukbppp.json`, `gene_coords_grch38.json`). The supplementary table lists the
brain imaging phenotypes and their accession numbers.

## Controlled access

UK Biobank individual-level data are not used here; the summary statistics listed above are
public, with the exception of any source that requires registration at its own portal.
Follow each source's terms of use.
