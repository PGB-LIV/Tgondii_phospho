"""
Cross-check orthogroup-based site mapping against direct peptidoform evidence.

Some phosphosites that are clearly the *same* site (same peptidoform observed
across proteins from multiple strains, with matching protein modification
positions) end up assigned to *different* (Group_ID, Aligned_pos) entries in
phosphosites_crossgenome_sites_v2.tsv. This happens when the ME49 gene model
has been split/merged relative to the other strains and the syntenic ortholog
mapping (AllXXX_genes_from_other_strains.tsv) was not updated to match.

This script:
  1. Builds a (Protein, Protein_pos) -> (Group_ID, Aligned_pos) lookup from
     phosphosites_crossgenome_sites_v2.tsv.
  2. Scans all_datasets_merged_Site_Peptidoform_centric_UniProt_tgondii.tsv for
     passing, non-decoy peptidoform rows whose "Proteins" / "Protein
     Modification Positions" span >=2 of our 4 strains.
  3. For each such row, looks up the (Group_ID, Aligned_pos) for each matched
     protein/position. If more than one distinct (Group_ID, Aligned_pos) is
     implicated by the SAME peptidoform, flags this as an ambiguous split.
  4. Writes:
       - ambiguous_orthogroup_splits.tsv  (peptidoform-level evidence)
       - ambiguous_site_pairs_summary.tsv (aggregated group/position pairs)
       - phosphosites_crossgenome_sites_v3.tsv (v2 + Ambiguous_split_with column)
"""
import csv
import os
import re
from collections import defaultdict

BASE        = os.path.dirname(os.path.abspath(__file__))
CROSS_FILE  = f"{BASE}/phosphosites_crossgenome_sites_v2.tsv"
PEP_FILE    = f"{BASE}/all_datasets_merged_Site_Peptidoform_centric_UniProt_tgondii.tsv"
OUT_EVIDENCE = f"{BASE}/ambiguous_orthogroup_splits.tsv"
OUT_SUMMARY  = f"{BASE}/ambiguous_site_pairs_summary.tsv"
OUT_FLAGGED  = f"{BASE}/phosphosites_crossgenome_sites_v3.tsv"

STRAINS = ["ME49", "GT1", "RH88", "VEG"]

SITE_RE = re.compile(r"^(.+):(-?\d+)\(([A-Za-z]+)\)$")


def strain_of(protein_id):
    if protein_id.startswith("TGME49_"):
        return "ME49"
    if protein_id.startswith("TGGT1_"):
        return "GT1"
    if protein_id.startswith("TGVEG_"):
        return "VEG"
    if re.match(r"^KAF\d+\.\d+$", protein_id):
        return "RH88"
    return None


# ══════════════════════════════════════════════════════════════════════════
# Step 1 - build (Protein, Protein_pos) -> (Group_ID, Aligned_pos) lookup
# ══════════════════════════════════════════════════════════════════════════
print("Loading crossgenome site table...")

site_lookup = {}   # (protein_id, protein_pos) -> (group_id, aligned_pos)
cross_rows = []

with open(CROSS_FILE, encoding="utf-8", newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        cross_rows.append(row)
        gid = row["Group_ID"]
        ap = row["Aligned_pos"]
        for s in STRAINS:
            sites = row[f"{s}_sites"]
            if not sites:
                continue
            for entry in sites.split(";"):
                m = SITE_RE.match(entry)
                if not m:
                    continue
                prot, pos, _res = m.group(1), int(m.group(2)), m.group(3)
                site_lookup[(prot, pos)] = (gid, ap)

print(f"  {len(site_lookup):,} (protein, pos) entries indexed from {len(cross_rows):,} crossgenome rows")


# ══════════════════════════════════════════════════════════════════════════
# Step 2 - scan peptidoform file for cross-strain agreement
# ══════════════════════════════════════════════════════════════════════════
print("\nScanning peptidoform-centric evidence file...")

evidence_rows = []          # rows for ambiguous_orthogroup_splits.tsv
pair_counts = defaultdict(lambda: {"n_peptidoforms": 0, "examples": []})

n_total = n_pass = n_multi_strain = n_ambiguous = 0

with open(PEP_FILE, encoding="utf-8", newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        n_total += 1
        if row["Site Passes Threshold [0.05]"] != "1":
            continue
        if row["Decoy Modification Site"] != "0":
            continue
        n_pass += 1

        proteins = row["Proteins"].split(":")
        positions = row["Protein Modification Positions"].split(":")

        matched = []   # (strain, protein, pos, group_id, aligned_pos)
        for prot, pos_str in zip(proteins, positions):
            strain = strain_of(prot)
            if strain is None:
                continue
            try:
                pos = int(pos_str)
            except ValueError:
                continue
            key = (prot, pos)
            if key in site_lookup:
                gid, ap = site_lookup[key]
                matched.append((strain, prot, pos, gid, ap))

        strains_seen = set(m[0] for m in matched)
        if len(strains_seen) < 2:
            continue
        n_multi_strain += 1

        distinct_targets = set((m[3], m[4]) for m in matched)
        if len(distinct_targets) <= 1:
            continue

        n_ambiguous += 1
        pep_id = row["Peptidoform Site ID"]
        psm_count = row.get("PSM Count Passing Threshold [0.05]", "")

        for strain, prot, pos, gid, ap in matched:
            evidence_rows.append({
                "Peptidoform_Site_ID": pep_id,
                "Strain": strain,
                "Protein": prot,
                "Protein_pos": pos,
                "Group_ID": gid,
                "Aligned_pos": ap,
                "PSM_count": psm_count,
                "Source_dataset": row.get("Source Dataset Identifier", ""),
            })

        # Aggregate pairwise across all distinct targets seen together
        targets = sorted(distinct_targets)
        for i in range(len(targets)):
            for j in range(i + 1, len(targets)):
                key = (targets[i], targets[j])
                pc = pair_counts[key]
                pc["n_peptidoforms"] += 1
                if len(pc["examples"]) < 3:
                    pc["examples"].append(pep_id)

print(f"  Total peptidoform rows:                 {n_total:,}")
print(f"  Passing 0.05 / non-decoy:                {n_pass:,}")
print(f"  Spanning >=2 of our strains (matched):   {n_multi_strain:,}")
print(f"  Ambiguous (>1 distinct group/pos hit):   {n_ambiguous:,}")


# ══════════════════════════════════════════════════════════════════════════
# Step 3 - write peptidoform-level evidence
# ══════════════════════════════════════════════════════════════════════════
print(f"\nWriting {OUT_EVIDENCE} ...")
EV_FIELDS = ["Peptidoform_Site_ID", "Strain", "Protein", "Protein_pos",
             "Group_ID", "Aligned_pos", "PSM_count", "Source_dataset"]
with open(OUT_EVIDENCE, "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=EV_FIELDS, delimiter="\t")
    writer.writeheader()
    writer.writerows(evidence_rows)
print(f"  {len(evidence_rows):,} rows written")


# ══════════════════════════════════════════════════════════════════════════
# Step 4 - write aggregated pair summary
# ══════════════════════════════════════════════════════════════════════════
print(f"\nWriting {OUT_SUMMARY} ...")
SUM_FIELDS = ["Group_ID_A", "Aligned_pos_A", "Group_ID_B", "Aligned_pos_B",
              "n_supporting_peptidoforms", "example_peptidoform_ids"]
with open(OUT_SUMMARY, "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=SUM_FIELDS, delimiter="\t")
    writer.writeheader()
    for (a, b), info in sorted(pair_counts.items(),
                                key=lambda kv: -kv[1]["n_peptidoforms"]):
        writer.writerow({
            "Group_ID_A": a[0],
            "Aligned_pos_A": a[1],
            "Group_ID_B": b[0],
            "Aligned_pos_B": b[1],
            "n_supporting_peptidoforms": info["n_peptidoforms"],
            "example_peptidoform_ids": ";".join(info["examples"]),
        })
print(f"  {len(pair_counts):,} distinct (group,pos) pairs flagged")


# ══════════════════════════════════════════════════════════════════════════
# Step 5 - flag the crossgenome table itself
# ══════════════════════════════════════════════════════════════════════════
print(f"\nWriting {OUT_FLAGGED} ...")

flag_map = defaultdict(set)   # (group_id, aligned_pos) -> set of partner (group_id, aligned_pos)
for (a, b) in pair_counts:
    flag_map[a].add(b)
    flag_map[b].add(a)

FIELDS = list(cross_rows[0].keys()) + ["Ambiguous_split_with"]
with open(OUT_FLAGGED, "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
    writer.writeheader()
    for row in cross_rows:
        key = (row["Group_ID"], row["Aligned_pos"])
        partners = flag_map.get(key)
        row["Ambiguous_split_with"] = (
            ";".join(f"{g}:{p}" for g, p in sorted(partners)) if partners else ""
        )
        writer.writerow(row)

n_flagged_rows = sum(1 for k in flag_map)
print(f"  {n_flagged_rows:,} crossgenome rows flagged as ambiguous")
print("\nDone.")
