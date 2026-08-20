"""
Add product description and gene name/symbol to ME49_phosphoprotein_localisation_merged.csv
by joining on Gene ID against ME49_gene_table.tsv.
"""

import pandas as pd

MERGED_FILE = "ME49_phosphoprotein_localisation_merged.csv"
GENE_TABLE_FILE = "ME49_gene_table.tsv"

merged = pd.read_csv(MERGED_FILE)
genes = pd.read_csv(
    GENE_TABLE_FILE,
    sep="\t",
    usecols=["Gene ID", "Product Description", "Gene Name or Symbol"],
).drop_duplicates(subset="Gene ID")

out = merged.merge(genes, on="Gene ID", how="left")
out.to_csv(MERGED_FILE, index=False)
print(f"Wrote {len(out)} rows to {MERGED_FILE}")
