# Phosphoprotein essentiality & sub-cellular localisation

Relates phosphosite quality tiers (Gold/Silver/Bronze, from a
false-localisation-rate-thresholded meta-analysis of *T. gondii*
phosphoproteomics datasets) to (a) CRISPR essentiality phenotype scores and
(b) hyperLOPIT sub-cellular localisation.

All three scripts use only relative filenames — no hardcoded paths — and run
from inside this folder.

## Required input files (all included in this folder)

| File | Source |
|---|---|
| `G2S1B_0.05_protein_pos_all_prot_mapping_without_contam.csv` | This project's core phosphosite table: one row per (protein, position, residue), FLR-thresholded at 0.05, with a `PTM_FLR_category` (Gold/Silver/Bronze) call per site. Same file used in `Toxo_syntenic_ortho_mapped/`. |
| `GT1_phenotype_data.tsv` | ToxoDB export of the GT1 genome-wide CRISPR screen phenotype scores (Sidik et al. 2016, *Cell*; PMID 27594428). Columns: `Gene ID`, `source_id`, `Organism`, `Genomic Location (Gene)`, `Product Description`, plus `T.gondii GT1 CRISPR Phenotype - Mean Phenotype` and four tissue composite scores (Liver/Lung/Peritoneum/Spleen). |
| `MCMC_LOPIT_data.tsv` | hyperLOPIT sub-cellular proteome data for ME49, from Barylyuk et al. 2020, *Cell Host & Microbe* (PMID 33053376). Columns used: `Accession` (ME49 gene ID), `tagm.mcmc.allocation` (TAGM-MCMC compartment call). |
| `ME49_gene_table.tsv` | ToxoDB gene table export for ME49. Columns used: `Gene ID`, `Product Description`, `Gene Name or Symbol`. |

## Pipeline / run order

Run from inside this folder (`python <script>.py`, no arguments).

### 1. `phosphoprotein_phenotype_boxplots.py`
Reads `G2S1B_0.05_...csv` (filtered to `TGGT1_` proteins) and
`GT1_phenotype_data.tsv`. For each GT1 gene, takes its best (highest-confidence)
phosphosite tier — Gold > Silver > Bronze — or "Unmodified" if none, and merges
this onto the phenotype table. Runs Mann-Whitney U (phosphoprotein vs not) and
Kruskal-Wallis (across phosphosite-count bins) tests.
Output: `GT1_phosphoprotein_phenotype_merged.csv`, `boxplot_Mean_Phenotype.png`.

### 2. `phosphoprotein_localisation.py`
Reads `G2S1B_0.05_...csv` (filtered to `TGME49_` proteins) and
`MCMC_LOPIT_data.tsv`. Same best-tier logic as above, merged onto hyperLOPIT
compartment calls; compartments with fewer than 30 assigned proteins are
dropped from the plots (not from the output CSV).
Output: `ME49_phosphoprotein_localisation_merged.csv` plus six summary plots
(`stacked_tier_by_location.png`, `stacked_binary_by_location.png`,
`heatmap_tier_by_location.png`, `dotplot_phosphoprotein_rate_by_location.png`,
`bar_top5_compartments_GSBN.png`, `bar_all_compartments_GSBN.png`).

### 3. `add_gene_annotations.py`
Reads `ME49_phosphoprotein_localisation_merged.csv` (output of step 2) and
`ME49_gene_table.tsv`; left-joins `Product Description` and
`Gene Name or Symbol` onto it by Gene ID, overwriting the file in place.

## Dependencies

Python 3, `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`.

```
pip install pandas numpy matplotlib seaborn scipy
```
