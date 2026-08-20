# Cross-strain (syntenic ortholog) phosphosite mapping

Maps *T. gondii* phosphosites detected in any of four strains (ME49, GT1,
RH88, VEG) onto shared syntenic-orthogroup alignment positions, so the same
biological site can be compared/counted across strains despite differing
gene models, protein lengths, or residue numbering.

All 5 scripts use `os.path.dirname(os.path.abspath(__file__))` to locate
their own folder — **no hardcoded paths**, run from anywhere.

## What's included vs. what was skipped

The original pipeline ran MUSCLE alignment of ~3,860 orthogroups as SLURM
array jobs on an HPC cluster — that step is **not** reproduced here. Instead,
the pre-computed alignment output (`orthogroup_alignments/`,
`orthogroup_alignments_same_length/` — ~3,860 small `.afa` files, ~28MB
total) is bundled directly as an input. Every other step is real code run
against real (non-derived) inputs.

## Required input files (all included in this folder)

| File | Source |
|---|---|
| `AllME49_genes_from_other_strains.tsv`, `AllGT1_...`, `AllRH88_...`, `AllVEG_...` | ToxoDB "genes by ortholog" bulk export, one per strain. Columns used: `Gene ID`, `Ortholog Group`, `Product Description`, `Input Ortholog(s)`. |
| `protein_fasta_files/ToxoDB-68_Tgondii{ME49,GT1,RH88,VEG}_AnnotatedProteins.fasta` | ToxoDB release 68 annotated protein FASTAs, one per strain. Headers must carry `gene=<gene_id>`. |
| `G2S1B_0.05_protein_pos_all_prot_mapping_without_contam.csv` | This project's core phosphosite table (same file used in `Essentiality_sub_cellular/`). |
| `all_datasets_merged_Site_Peptidoform_centric_UniProt_tgondii.tsv` | This project's peptidoform-centric merged search-engine output across all 10 datasets in the meta-analysis (~105MB). Columns used: `Proteins`, `Protein Modification Positions`, `Site Passes Threshold [0.05]`, `Decoy Modification Site`, `Peptidoform Site ID`, `PSM Count Passing Threshold [0.05]`, `Source Dataset Identifier`. |
| `orthogroup_alignments/*.afa`, `orthogroup_alignments_same_length/*.afa` | Pre-computed MUSCLE v5.3 alignments (see above — MUSCLE/SLURM step not reproduced here). |

## Pipeline / run order

Run from inside this folder (`python <script>.py`, no arguments).

### 1. `build_orthogroup_summary.py`
Inputs: the four `AllXXX_genes_from_other_strains.tsv` files.
Groups genes into syntenic orthogroups via connected components on the
pairwise `Input Ortholog(s)` links (deliberately *not* grouped by the
broader OrthoMCL/OG family ID, so paralogs sharing an OG but with distinct
syntenic partners — e.g. MIC17A/B/C — are correctly kept separate).
Output: `syntenic_orthogroup_summary_v2.tsv`

### 2. `map_phosphosites_to_orthogroups.py`
Inputs: the four protein FASTAs, `G2S1B_0.05_...csv`, `syntenic_orthogroup_summary_v2.tsv`
For each orthogroup, picks one representative protein per gene (the isoform
actually carrying phosphosites, if any). Orthogroups where every member is
the same length get `Aligned_pos = Protein_pos` directly (no indels
assumed); groups with differing lengths are flagged `needs_MSA` — the
per-orthogroup FASTAs this step would write for MUSCLE are **not** needed
here since the alignments already exist in `orthogroup_alignments/`.
Output: `phosphosites_orthogroup_annotated.tsv`

### 3. `build_aligned_site_table.py`
Inputs: `phosphosites_orthogroup_annotated.tsv`, `orthogroup_alignments/*.afa`
Resolves every phosphosite to a final `(Group_ID, Aligned_pos)` — directly
for same-length groups, via alignment-column lookup (parsing the `.afa`
files) for MSA groups.
Output: `phosphosites_crossgenome_sites_v2.tsv` — one row per
(orthogroup, aligned position), with per-strain hit counts, site details,
and best FLR category.

### 4. `check_orthogroup_splits.py`
Inputs: `phosphosites_crossgenome_sites_v2.tsv`, `all_datasets_merged_Site_Peptidoform_centric_UniProt_tgondii.tsv`
QC cross-check: finds peptidoforms observed identically across ≥2 strains
whose matched protein positions land in *different* orthogroup/aligned-pos
entries (i.e. cases where the syntenic mapping likely mis-split a gene model
relative to direct peptidoform evidence), and flags them.
Outputs: `ambiguous_orthogroup_splits.tsv`, `ambiguous_site_pairs_summary.tsv`,
**`phosphosites_crossgenome_sites_v3.tsv`** (v2 + `Ambiguous_split_with`
column) — **this is `Data_S4`**.

### 5. `recover_phosphosites_from_alignment.py`
Inputs: `phosphosites_crossgenome_sites_v3.tsv`, `syntenic_orthogroup_summary_v2.tsv`, both `orthogroup_alignments*/` directories
For aligned positions where ≥1 strain has a detected site but another
strain doesn't, checks whether that other strain's residue at the same
aligned column is S/T and its "+1" motif residue matches a strain that does
have the site — if so, flags it as a recoverable site (a real but
MS-undetected site, e.g. due to a differing tryptic peptide from sequence
divergence).
Output: **`phosphosites_crossgenome_sites_v4_recovered.tsv`** — **this is `Data_S3`**.



## Dependencies

Python 3, `pandas` (only used indirectly via `csv`/stdlib in these 5
scripts — no third-party packages required for this subset of the
pipeline).
