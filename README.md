<h1 align="center"><b>A comprehensive meta-analysis of phosphosites in the protozoan parasite <i>Toxoplasma gondii</i></b></h1>

# Project Description
This analysis is based on a reanalysis of <i>Toxoplasma gondii</i> (T gondii) phosphorylation datasets (PRIDE accession: [PXD078724](https://www.ebi.ac.uk/pride/archive/projects/PXD078724)) which used [alanine as a decoy](https://github.com/PGB-LIV/mzidFLR) to control for the false localisation rate.


| File | Description |
|------|-------------|
| `00_preprocess_and_summary/scripts/01_ME49_GSB.py` | Code to generate single protein mapping file (ME49 only) |
| `00_preprocess_and_summary/scripts/02_summary_figure.py` | Code to generate Figure 1 in manuscript (GSB summary) |
| `00_preprocess_and_summary/scripts/03_upset_plot.R` | Code to determine overlap of peptidoform sites in T gondii strains |
| `01_motif_analysis/scripts/01_15mer_background_v3.py` | Code to generate background for motif analysis - all ST sites in ME49 proteins |
| `01_motif_analysis/scripts/02_overall_background.py` | Code to aggregate all ST sites identified across datasets |
| `01_motif_analysis/scripts/03_motif_foreground.py` | Code to generate foreground of ME49 GSB ST sites for motif analysis |
| `01_motif_analysis/scripts/04_motif.R` | Code for rmotifx analysis (ME49) |
| `02_metapredict/scripts/01_metapredict_processing_ME49.py` | Code for processing the output file from metapredict |
| `02_metapredict/scripts/02_disorder_plot.R` | Code to generate boxplots from the disorder analysis scores for ST sites versus phosphorylated ST sites (ME49 single mapping) |
