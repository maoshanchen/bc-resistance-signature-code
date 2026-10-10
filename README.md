# Reproduction bundle

Code for the manuscript:

*An overall-survival-derived gene-expression score reads proliferation and
prognosis, not acquired drug resistance: multi-cohort audit and repositioning
in HR+/HER2− disease.*

Every input is a public accession (see the manuscript's Data Availability
Statement), so all tables and figures regenerate end to end once the raw
public data are placed in the expected locations.

## Repository layout

```
├── README.md                 this file
├── LICENSE                   MIT
├── requirements.txt          exact package versions (Python 3.14.7)
├── MANIFEST.txt              SHA256 + size for every file in this bundle
├── scripts/                  all analysis code, numbered in pipeline order
│   └── stats_toolkit.py      shared statistics library (imported by the audit
│                             scripts; intentionally left unnumbered)
├── example/                  minimal runnable example (no large data)
│   ├── cohort_example.csv    tiny synthetic hotspot-panel cohort (7 patients)
│   ├── resistance_genes.csv  the resistance-panel gene list (29 genes)
│   └── expected_output/      output produced by scripts/04_score_cohort.py
├── config/                   (not used — each script carries its configuration
│                             constants at the top of the file)
└── results/                  output directory (see results/README.md);
    └── tables/               the derived submission tables (Table 1–3, S1–S5)
```

## Environment

- Python 3.14.7
- Install exact versions: `pip install -r requirements.txt`
- Key packages: numpy 2.5.3, pandas 2.3.3, scipy 1.18.1, scikit-learn 1.9.1,
  statsmodels 0.15.0, lifelines 0.30.3, scanpy 1.12.4, anndata 0.13.3.post0,
  harmonypy 0.2.0, matplotlib 3.11.2, openpyxl 3.1.5, python-docx 1.2.0

## Quick start (the example)

```bash
pip install -r requirements.txt
python scripts/04_score_cohort.py
```

This scores the tiny example cohort in `example/` and writes
`example/expected_output/cohort_ras_scored.csv`. To run it on your own
hotspot-panel cohort, edit the three path constants at the top of
`scripts/04_score_cohort.py` and prepare a cohort table with the column
layout documented there.

## Data (not included; all public)

Fetch the accessions below with `scripts/01_download_gse.py` /
`scripts/02_download_gpl.py` or the cBioPortal API, and arrange them as:

```
data/
├── 01_core/                 C1 GSE26459 (in vitro); C3 TCGA-BRCA; C4 METABRIC (cBioPortal brca_metabric)
├── 02_acquired_resistance/  C2 GSE45255 (treatment-dose series)
├── 03_single_cell/          C5 GSE264205 (scRNA-seq); GSE249361; GSE161529 (derived synthesis S1)
├── 04_api_reference/        cBioPortal API reference downloads (brca_aurora_2023, brca_mbcproject_2022)
├── 05_gpl/                  GEO platform annotation files
└── 06_metastatic/           C6 AURORA US (GSE209998); C7 MET500; C8 MBCproject/Count Me In
```

Cohort codes C1–C8 and the derived synthesis S1 follow the manuscript's
Table 1 (`results/tables/Table1_cohort_inventory.csv`).

## Pipeline order (main analyses)

1. **Data acquisition** — `scripts/01_download_gse.py`,
   `scripts/02_download_gpl.py`, `scripts/03_download_fix.py`;
   cBioPortal queries (`scripts/05_survey_cbioportal_cohorts.py`),
   cohort probing (`06`) and panel coverage (`07`).
2. **In-vitro / early-cohort preparation** — `08`–`11` (4-OHT labels,
   DE-gene mapping, GSE resistance comparison, gene frequencies).
3. **Score construction** — `scripts/12_rep_build.py` (Cox-weighted REP30
   assembly; METABRIC = fitting cohort, TCGA-BRCA held out),
   `scripts/14_rep30_resistance_validation.py` and
   `scripts/15_rep30_resistance_within_er.py` (the 461-comparison audit);
   `scripts/13`, `16`–`19` (HR+/HER2− refit, GEO validation, panel trimming).
4. **Statistical toolkit** — `scripts/stats_toolkit.py` is imported by all
   audit scripts (Firth penalised logistic regression, permutation tests,
   Cliff's delta, Hodges–Lehmann, Newcombe CI). Run nothing without it.
5. **Prognostic layer** — `scripts/20`–`25` (backtest, per-cohort
   performance, subtype analysis, HR+/HER2− metrics).
6. **Null benchmarks** — `scripts/26_audit_venet_null.py` …
   `30_audit_venet_null5.py` (random-gene-set nulls A–D),
   `scripts/31`–`34` (C-index cost, transfer stability, proliferation
   control), `scripts/35`–`36` (sensitivity demo + plots).
7. **Single-cell layer** — `scripts/37_sc_pipeline.py` …
   `45_sc_integrated_report.py` (Harmony integration, per-study DE,
   resistance-evolution trajectories).
8. **Metastatic / clinical layer** — `scripts/46`–`57` (MET500, AURORA,
   MBCproject, CDK4/6 genomic resistance, IGR2015 landscape, RAS scoring).
9. **Tables and figures** — `scripts/59_build_tables.py` (submission
   tables → `results/tables/`), `scripts/60_build_manuscript_figures.py` +
   `scripts/61_export_panels.py` (figures → `results/figures/`);
   `scripts/58` (evidence table) and `scripts/62` (submission docx,
   requires the manuscript source, not included) are internal builders.
   The ready-made derived tables are committed under `results/tables/`.

## Important note on paths

All path constants in the scripts are repository-relative: each script
resolves the repository root via `Path(__file__).resolve().parents[n]`, so
the repository can be cloned anywhere and run as-is, provided the input
layout under `data/` is recreated (see the download scripts).
`scripts/04_score_cohort.py` is the generic entry point for applying the
resistance-alteration score to any prepared cohort table.
