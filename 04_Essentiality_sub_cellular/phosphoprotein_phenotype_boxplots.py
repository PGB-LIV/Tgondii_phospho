"""
Build phosphoprotein quality tiers (Gold/Silver/Bronze) for T. gondii GT1 strain,
merge with CRISPR phenotype scores, and plot boxplots:
  Left:  Not phosphoprotein vs Phosphoprotein (any tier)
  Right: Unmodified / Bronze / Silver / Gold
for the Mean Phenotype score.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, kruskal

SITE_FILE = "G2S1B_0.05_protein_pos_all_prot_mapping_without_contam.csv"
PHENO_FILE = "GT1_phenotype_data.tsv"

# Colours for Gold / Silver / Bronze / Unmodified (and Phospho/Not-phospho)
TIER_COLORS = {
    "Gold": "#D4AF37",
    "Silver": "#A8A9AD",
    "Bronze": "#AD6F3B",
    "Unmodified": "#6E91C7",
}
BINARY_COLORS = {
    "Phosphoprotein": "#7B3F9E",
    "Not phosphoprotein": "#6E91C7",
}

TIER_ORDER = ["Unmodified", "Bronze", "Silver", "Gold"]
BINARY_ORDER = ["Not phosphoprotein", "Phosphoprotein"]

SITE_BINS = [-0.5, 0.5, 1.5, 4.5, 10.5, np.inf]
SITE_BIN_LABELS = ["0", "1", "2-4", "5-10", ">10"]
SITE_BIN_COLORS = {
    "0": "#6E91C7",
    "1": "#9DBF9E",
    "2-4": "#D9C25C",
    "5-10": "#D98E4B",
    ">10": "#B23A48",
}

PHENOTYPE_COLS = [
    "T.gondii GT1 CRISPR Phenotype - Mean Phenotype",
]


def load_protein_tiers():
    df = pd.read_csv(SITE_FILE, usecols=["Protein", "PTM_FLR_category"])
    df = df[df["Protein"].str.startswith("TGGT1_")].copy()
    df["Gene ID"] = df["Protein"].str.split("-", n=1).str[0]

    tier_rank = {"Gold": 3, "Silver": 2, "Bronze": 1}
    df["rank"] = df["PTM_FLR_category"].map(tier_rank)

    best = df.groupby("Gene ID")["rank"].max().reset_index()
    rank_to_tier = {3: "Gold", 2: "Silver", 1: "Bronze"}
    best["Tier"] = best["rank"].map(rank_to_tier)

    site_counts = df.groupby("Gene ID").size().reset_index(name="Site_count")
    best = best.merge(site_counts, on="Gene ID", how="left")
    return best[["Gene ID", "Tier", "Site_count"]]


def label_with_n(df, col, order):
    counts = df[col].value_counts()
    return [f"{cat}\n(n={counts.get(cat, 0)})" for cat in order]


def main():
    tiers = load_protein_tiers()
    pheno = pd.read_csv(PHENO_FILE, sep="\t")

    merged = pheno.merge(tiers, on="Gene ID", how="left")
    merged["Tier"] = merged["Tier"].fillna("Unmodified")
    merged["Site_count"] = merged["Site_count"].fillna(0).astype(int)
    merged["Phospho_status"] = merged["Tier"].apply(
        lambda t: "Not phosphoprotein" if t == "Unmodified" else "Phosphoprotein"
    )

    print("Protein tier counts:")
    print(merged["Tier"].value_counts().reindex(TIER_ORDER))
    print()
    print("Phospho status counts:")
    print(merged["Phospho_status"].value_counts().reindex(BINARY_ORDER))

    merged.to_csv("GT1_phosphoprotein_phenotype_merged.csv", index=False)

    for pheno_col in PHENOTYPE_COLS:
        plot_df = merged.dropna(subset=[pheno_col]).copy()
        plot_df["Site_bin"] = pd.cut(
            plot_df["Site_count"], bins=SITE_BINS, labels=SITE_BIN_LABELS
        )
        short_name = pheno_col.split(" - ")[-1]

        fig, axes = plt.subplots(1, 3, figsize=(13, 5), sharey=True)

        sns.boxplot(
            data=plot_df,
            x="Phospho_status",
            y=pheno_col,
            order=BINARY_ORDER,
            hue="Phospho_status",
            palette=BINARY_COLORS,
            legend=False,
            ax=axes[0],
            showfliers=False,
        )
        sns.stripplot(
            data=plot_df,
            x="Phospho_status",
            y=pheno_col,
            order=BINARY_ORDER,
            color="black",
            alpha=0.06,
            size=2,
            ax=axes[0],
        )
        axes[0].set_xlabel("")
        axes[0].set_ylabel(short_name)
        axes[0].tick_params(axis="x", rotation=20)
        axes[0].set_xticks(range(len(BINARY_ORDER)))
        axes[0].set_xticklabels(
            label_with_n(plot_df, "Phospho_status", BINARY_ORDER)
        )
        axes[0].text(
            -0.05, 1.05, "A", transform=axes[0].transAxes,
            fontsize=14, fontweight="bold", va="bottom", ha="right",
        )

        group_a = plot_df.loc[
            plot_df["Phospho_status"] == BINARY_ORDER[0], pheno_col
        ]
        group_b = plot_df.loc[
            plot_df["Phospho_status"] == BINARY_ORDER[1], pheno_col
        ]
        stat, pval = mannwhitneyu(group_a, group_b, alternative="two-sided")
        p_text = "p < 0.001" if pval < 0.001 else f"p = {pval:.3g}"
        axes[0].text(
            0.5,
            0.98,
            f"Mann-Whitney U {p_text}",
            transform=axes[0].transAxes,
            ha="center",
            va="top",
            fontsize=9,
        )

        sns.boxplot(
            data=plot_df,
            x="Tier",
            y=pheno_col,
            order=TIER_ORDER,
            hue="Tier",
            palette=TIER_COLORS,
            legend=False,
            ax=axes[1],
            showfliers=False,
        )
        sns.stripplot(
            data=plot_df,
            x="Tier",
            y=pheno_col,
            order=TIER_ORDER,
            color="black",
            alpha=0.06,
            size=2,
            ax=axes[1],
        )
        axes[1].set_xlabel("")
        axes[1].set_ylabel("")
        axes[1].tick_params(axis="x", rotation=20)
        axes[1].set_xticks(range(len(TIER_ORDER)))
        axes[1].set_xticklabels(label_with_n(plot_df, "Tier", TIER_ORDER))
        axes[1].text(
            -0.05, 1.05, "B", transform=axes[1].transAxes,
            fontsize=14, fontweight="bold", va="bottom", ha="right",
        )

        sns.boxplot(
            data=plot_df,
            x="Site_bin",
            y=pheno_col,
            order=SITE_BIN_LABELS,
            hue="Site_bin",
            palette=SITE_BIN_COLORS,
            legend=False,
            ax=axes[2],
            showfliers=False,
        )
        sns.stripplot(
            data=plot_df,
            x="Site_bin",
            y=pheno_col,
            order=SITE_BIN_LABELS,
            color="black",
            alpha=0.06,
            size=2,
            ax=axes[2],
        )
        axes[2].set_xlabel("")
        axes[2].set_ylabel("")
        axes[2].tick_params(axis="x", rotation=20)
        axes[2].set_xticks(range(len(SITE_BIN_LABELS)))
        axes[2].set_xticklabels(label_with_n(plot_df, "Site_bin", SITE_BIN_LABELS))
        axes[2].text(
            -0.05, 1.05, "C", transform=axes[2].transAxes,
            fontsize=14, fontweight="bold", va="bottom", ha="right",
        )

        groups = [
            plot_df.loc[plot_df["Site_bin"] == cat, pheno_col]
            for cat in SITE_BIN_LABELS
        ]
        kw_stat, kw_pval = kruskal(*groups)
        kw_text = "p < 0.001" if kw_pval < 0.001 else f"p = {kw_pval:.3g}"
        axes[2].text(
            0.5,
            0.98,
            f"Kruskal-Wallis {kw_text}",
            transform=axes[2].transAxes,
            ha="center",
            va="top",
            fontsize=9,
        )

        fig.suptitle(short_name)
        fig.tight_layout()

        out_name = f"boxplot_{short_name.replace(' ', '_')}.png"
        fig.savefig(out_name, dpi=300)
        plt.close(fig)
        print(f"Saved {out_name}")


if __name__ == "__main__":
    main()
