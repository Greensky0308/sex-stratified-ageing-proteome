# Sex-stratified cis-pQTL Mendelian randomization of the ageing proteome

Analysis code for a sex-stratified cis-pQTL Mendelian randomization screen of ageing
phenotypes. For each circulating protein the causal effect on an ageing phenotype is
estimated separately in women and men, the two effects are tested for a difference, and
candidate signals are carried through colocalization, cross-platform and external-cohort
checks.

## Repository layout

```
code/      analysis scripts (this directory)
data/      input summary statistics (not distributed here; see data/README.md)
out/       intermediate and derived tables written by the scripts
tables/    the phenotype summary table
figures/   the figures
```

Each script resolves the repository root as the parent of `code/` and reads and writes
relative to it. Output directories are created automatically on first run.

## Requirements

Python 3.10 or later:

```
pip install -r requirements.txt
```

`curl` must be available on the PATH; it is used for ranged reads of remote GWAS files.
If your network requires a proxy, export `HTTPS_PROXY` before running the scripts that read
remote files.

## Data

Input files are public or controlled-access summary statistics and are not redistributed
with this code. `data/README.md` lists each input, its accession or download route, and
where it is expected to sit under `data/`.

## Running

Run the scripts from the repository root, in the order below. Each step writes a table to
`out/` that later steps read.

| order | script | purpose |
|---|---|---|
| 1 | `liftover_instruments.py` | lift cis-pQTL instruments from GRCh37 to GRCh38 |
| 2 | `sex_diff_mr.py` | sex-differential MR against brain imaging phenotypes |
| 3 | `alm_sexdiff.py` | sex-differential MR against appendicular lean mass |
| 4 | `alm_locus_dump.py` | regional z-score dumps around the lead locus |
| 5 | `coloc_abf.py` | colocalization (ABF) utilities |
| 6 | `sult2a1_coloc_windows.py` | colocalization across window widths |
| 7 | `sult2a1_crossplatform.py` | cross-platform check of the lead cis-pQTL |
| 8 | `igf1_coloc.py` | colocalization of the locus with IGF-1 |
| 9 | `mediation_mr.py` | two-step mediation estimates |
| 10 | `reverse_mr.py` | reverse-direction MR |
| 11 | `neale_ffm_sexdiff.py` | sex-differential MR against body-composition traits |
| 12 | `neale_panel_phewas.py` | variant-level trait panel |
| 13 | `fetch_neale_snp.py` | ranged read of a single variant from a remote GWAS file |
| 14 | `su2025_screen.py` | screen in the external lean-mass cohort |
| 15 | `cst3_crossplatform.py` | analysis of the creatinine / cystatin C index |
| 16 | `cst3_locus_inspect.py` | shared readers for the cystatin C analyses |
| 17 | `cst3_locus_dump.py` | regional dumps at the cystatin C locus |
| 18 | `make_figures.py` | figures |
| 19 | `make_phenotype_table.py` | phenotype-wise summary table |
| 20 | `igf1_coloc_n.py` | number of SNPs entering each IGF-1 colocalization window |
| 21 | `pag1_coloc_reproduce.py` | colocalization of the PAG1 pQTL with the strongest brain phenotypes |
| 22 | `sult2a1_locus_dump.py` | regional dumps at the lead locus |
| 23 | `mw_locus_dump.py` | regional dumps for the grip-based weakness trait |
| 24 | `make_figS1.py` | supplementary figure for the creatinine / cystatin C index |

Example:

```
python code/liftover_instruments.py
python code/alm_sexdiff.py
```

## Notes

- No credentials are stored in this code. Network access, where used, reads the proxy from
  the `HTTPS_PROXY` environment variable and is otherwise direct.
- All paths are resolved relative to the repository root; no absolute paths are hard-coded.
- The scripts are deterministic given the input files; intermediate outputs in `out/` are
  overwritten on each run.

## License

MIT — see `LICENSE`.
