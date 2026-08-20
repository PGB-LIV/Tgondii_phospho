"""
Build phosphoprotein quality tiers (Gold/Silver/Bronze) for T. gondii ME49 strain
and relate them to subcellular localisation (MCMC-LOPIT TAGM allocation).

Produces:
  - Stacked proportion bar chart: composition of G/S/B/N within each location
  - Stacked proportion bar chart: composition of Phosphoprotein/Not within each location
  - Heatmap of G/S/B/N proportions by location (alternative visual)
  - Dot plot of phosphoprotein proportion per location, ranked
  - Grouped bar chart of G/S/B/N proportions for the top 5 compartments
    by phosphoprotein rate
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

SITE_FILE = "G2S1B_0.05_protein_pos_all_prot_mapping_without_contam.csv"
LOPIT_FILE = "MCMC_LOPIT_data.tsv"
MIN_PROTEINS_PER_LOCATION = 30

TIER_COLORS = {
    "Gold": "#D4AF37",
    "Silver": "#A8A9AD",
    "Bronze": "#AD6F3B",
    "Unmodified": "#6E91C7",
}
TIER_ORDER = ["Gold", "Silver", "Bronze", "Unmodified"]

BINARY_COLORS = {
    "Phosphoprotein": "#7B3F9E",
    "Not phosphoprotein": "#6E91C7",
}
BINARY_ORDER = ["Phosphoprotein", "Not phosphoprotein"]


def load_protein_tiers(strain_prefix):
    df = pd.read_csv(SITE_FILE, usecols=["Protein", "PTM_FLR_category"])
    df = df[df["Protein"].str.startswith(strain_prefix)].copy()
    df["Gene ID"] = df["Protein"].str.split("-", n=1).str[0]

    tier_rank = {"Gold": 3, "Silver": 2, "Bronze": 1}
    df["rank"] = df["PTM_FLR_category"].map(tier_rank)

    best = df.groupby("Gene ID")["rank"].max().reset_index()
    rank_to_tier = {3: "Gold", 2: "Silver", 1: "Bronze"}
    best["Tier"] = best["rank"].map(rank_to_tier)
    return best[["Gene ID", "Tier"]]


def main():
    tiers = load_protein_tiers("TGME49_")

    lopit = pd.read_csv(
        LOPIT_FILE, sep="\t", usecols=["Accession", "tagm.mcmc.allocation"]
    )
    lopit = lopit.rename(
        columns={"Accession": "Gene ID", "tagm.mcmc.allocation": "Location"}
    )

    merged = lopit.merge(tiers, on="Gene ID", how="left")
    merged["Tier"] = merged["Tier"].fillna("Unmodified")
    merged["Phospho_status"] = merged["Tier"].apply(
        lambda t: "Not phosphoprotein" if t == "Unmodified" else "Phosphoprotein"
    )

    loc_counts = merged["Location"].value_counts()
    kept_locations = loc_counts[loc_counts >= MIN_PROTEINS_PER_LOCATION].index.tolist()
    dropped = loc_counts[loc_counts < MIN_PROTEINS_PER_LOCATION]
    if len(dropped):
        print(f"Dropping locations with < {MIN_PROTEINS_PER_LOCATION} proteins:")
        print(dropped)
        print()

    plot_df = merged[merged["Location"].isin(kept_locations)].copy()
    plot_df.to_csv("ME49_phosphoprotein_localisation_merged.csv", index=False)

    # order locations by total protein count, descending
    loc_order = (
        plot_df["Location"].value_counts().sort_values(ascending=False).index.tolist()
    )

    tier_props = (
        plot_df.groupby("Location")["Tier"]
        .value_counts(normalize=True)
        .unstack(fill_value=0)
        .reindex(columns=TIER_ORDER)
        .reindex(loc_order)
    )
    binary_props = (
        plot_df.groupby("Location")["Phospho_status"]
        .value_counts(normalize=True)
        .unstack(fill_value=0)
        .reindex(columns=BINARY_ORDER)
        .reindex(loc_order)
    )
    loc_n = plot_df["Location"].value_counts().reindex(loc_order)
    loc_labels = [f"{loc} (n={loc_n[loc]})" for loc in loc_order]

    plot_stacked_bar(
        tier_props, loc_labels, TIER_COLORS, TIER_ORDER,
        "Proportion of proteins", "GT1 phosphoprotein tier composition by ME49 subcellular location",
        "stacked_tier_by_location.png",
    )
    plot_stacked_bar(
        binary_props, loc_labels, BINARY_COLORS, BINARY_ORDER,
        "Proportion of proteins", "Phosphoprotein composition by ME49 subcellular location",
        "stacked_binary_by_location.png",
    )
    plot_heatmap(tier_props, loc_labels)
    plot_dotplot(binary_props, loc_labels)
    plot_tier_location_bar(plot_df, top_n=5)
    plot_tier_location_bar(plot_df, top_n=None)


def plot_stacked_bar(props, loc_labels, colors, order, ylabel, title, out_name):
    fig, ax = plt.subplots(figsize=(8, max(5, 0.35 * len(props))))
    left = np.zeros(len(props))
    y_pos = np.arange(len(props))
    for cat in order:
        ax.barh(y_pos, props[cat], left=left, color=colors[cat], label=cat)
        left += props[cat].to_numpy()
    ax.set_yticks(y_pos)
    ax.set_yticklabels(loc_labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(ylabel)
    ax.set_xlim(0, 1)
    ax.set_title(title, fontsize=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.06 - 0.002 * len(props)), ncol=len(order), frameon=False)
    fig.tight_layout()
    fig.savefig(out_name, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_name}")


def plot_heatmap(tier_props, loc_labels):
    fig, ax = plt.subplots(figsize=(5, max(5, 0.35 * len(tier_props))))
    sns.heatmap(
        tier_props.to_numpy(),
        annot=True,
        fmt=".2f",
        cmap="viridis",
        yticklabels=loc_labels,
        xticklabels=tier_props.columns,
        cbar_kws={"label": "Proportion of proteins"},
        ax=ax,
    )
    ax.set_title("G/S/B/N proportion by location", fontsize=10)
    fig.tight_layout()
    fig.savefig("heatmap_tier_by_location.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved heatmap_tier_by_location.png")


def plot_dotplot(binary_props, loc_labels):
    overall_rate = binary_props["Phosphoprotein"].mean()
    order_idx = binary_props["Phosphoprotein"].to_numpy().argsort()[::-1]
    sorted_props = binary_props.iloc[order_idx]
    sorted_labels = [loc_labels[i] for i in order_idx]

    fig, ax = plt.subplots(figsize=(6, max(5, 0.35 * len(sorted_props))))
    y_pos = np.arange(len(sorted_props))
    ax.scatter(
        sorted_props["Phosphoprotein"], y_pos, color=BINARY_COLORS["Phosphoprotein"], s=60, zorder=3
    )
    ax.axvline(overall_rate, color="gray", linestyle="--", linewidth=1, label=f"Overall rate ({overall_rate:.2f})")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sorted_labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Proportion phosphoprotein")
    ax.set_xlim(0, 1)
    ax.legend(loc="lower right", frameon=False, fontsize=8)
    ax.set_title("Phosphoprotein enrichment by subcellular location", fontsize=10)
    fig.tight_layout()
    fig.savefig("dotplot_phosphoprotein_rate_by_location.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved dotplot_phosphoprotein_rate_by_location.png")


def plot_tier_location_bar(plot_df, top_n=None):
    # counts of proteins per (Tier, Location)
    counts = (
        plot_df.groupby(["Tier", "Location"]).size().unstack(fill_value=0)
        .reindex(index=TIER_ORDER, fill_value=0)
    )
    tier_totals = counts.sum(axis=1)

    # locations ordered by raw Gold protein count, optionally truncated to top_n
    locs = counts.loc["Gold"].sort_values(ascending=False).index.tolist()
    if top_n is not None:
        locs = locs[:top_n]

    # within each tier, proportion of that tier's proteins found in each location
    tier_loc_props = counts[locs].div(tier_totals, axis=0)

    fig, ax = plt.subplots(figsize=(max(9, 0.6 * len(locs)), 6))
    x = np.arange(len(locs))
    bar_width = 0.2
    for i, tier in enumerate(TIER_ORDER):
        ax.bar(
            x + (i - 1.5) * bar_width,
            tier_loc_props.loc[tier],
            width=bar_width,
            color=TIER_COLORS[tier],
            label=f"{tier} (n={tier_totals[tier]})",
            edgecolor="white",
            linewidth=0.5,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(locs, fontsize=9, rotation=35, ha="right")
    ax.set_ylabel("Proportion of tier's proteins in compartment")
    title = "Compartments ranked by Gold phosphoprotein count"
    if top_n is not None:
        title = f"Top {top_n} {title.lower()}"
    ax.set_title(title, fontsize=11)
    ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.32))
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    out_name = "bar_all_compartments_GSBN.png" if top_n is None else "bar_top5_compartments_GSBN.png"
    fig.savefig(out_name, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_name}")


if __name__ == "__main__":
    main()
