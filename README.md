<h1 align="center"><b>A comprehensive meta-analysis of phosphosites in the protozoan parasite <i>Toxoplasma gondii</i></b></h1>

## Project Description
This analysis is based on a reanalysis of <i>Toxoplasma gondii</i> (<i>T. gondii</i>) phosphorylation datasets (PRIDE accession: [PXD078724](https://www.ebi.ac.uk/pride/archive/projects/PXD078724)) which used [alanine as a decoy](https://github.com/PGB-LIV/mzidFLR) to control for the false localisation rate.

## Contents
| File | Description |
|------|-------------|
| `00_preprocess_and_summary/scripts/01_ME49_GSB.py` | Code to generate single protein mapping file (ME49 only). |
| `00_preprocess_and_summary/scripts/02_summary_figure.py` | Code to generate `Figure 1` in manuscript (GSB summary). |
| `00_preprocess_and_summary/scripts/03_upset_plot.R` | Code to determine overlap of peptidoform sites in *T. gondii* strains. |
| `01_Toxo_syntenic_ortho_mapped/` | Maps phosphosites onto syntenic orthogroup alignments across four *T. gondii* strains (ME49/GT1/RH88/VEG) so sites can be compared across strains. Reproduces `Data_S3`/`Data_S4` (`SupplementaryFile_sites_recovered_sites_across_strains.tsv`, `SuppFile_phosphosites_crossgenome_sites.tsv`). See its own README for input files and run order. |
| `02_motif_analysis/scripts/01_15mer_background_v3.py` | Code to generate background for motif analysis - all ST sites in ME49 proteins. |
| `02_motif_analysis/scripts/02_overall_background.py` | Code to aggregate all ST sites identified across datasets. |
| `02_motif_analysis/scripts/03_motif_foreground.py` | Code to generate foreground of ME49 GSB ST sites for motif analysis. |
| `02_motif_analysis/scripts/04_motif.R` | Code for rmotifx analysis (ME49). |
| `03_metapredict/scripts/01_metapredict_processing_ME49.py` | Code for processing the output file from metapredict. |
| `03_metapredict/scripts/02_disorder_plot.R` | Code to generate boxplots from the disorder analysis scores for ST sites versus phosphorylated ST sites (ME49 single mapping). |
| `04_Essentiality_sub_cellular/` | Relates phosphosite quality tiers (Gold/Silver/Bronze) to CRISPR essentiality phenotypes (GT1) and hyperLOPIT sub-cellular localisation (ME49). Reproduces `Data_S7`/`Data_S8` (`Supplement_essentiality_scores.csv`, `Supplement_Localisation_data.csv`). See its own README for input files and run order. |
