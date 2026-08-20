"""
Build the final cross-genome phosphosite table.

One row per (orthogroup, aligned_position).
For each strain: a count column (0/1/>1 for complex groups) and a detail
column listing protein_id:original_pos(residue) entries.

Three route into aligned positions:
  1. Same-length groups  → aligned_pos is already the protein pos (no MSA needed)
  2. MSA groups with .afa file → map protein pos to alignment column
  3. MSA groups with missing .afa → kept as "alignment_missing"
"""
import csv
import os
from collections import defaultdict

BASE       = os.path.dirname(os.path.abspath(__file__))
ANNOT_FILE = f"{BASE}/phosphosites_orthogroup_annotated.tsv"
ALIGN_DIR  = f"{BASE}/orthogroup_alignments"
OUT_FILE   = f"{BASE}/phosphosites_crossgenome_sites_v2.tsv"

STRAINS = ["ME49", "GT1", "RH88", "VEG"]


# ══════════════════════════════════════════════════════════════════════════════
# Helper: parse a MUSCLE .afa file -> {gene_id: {seq_pos: align_col}}
# Header format: >gene_id|protein_id|strain|len=NNN
# ══════════════════════════════════════════════════════════════════════════════
def parse_afa(path):
    pos_maps = {}   # gene_id -> {1-based seq_pos: 1-based alignment column}
    current_gene = None
    parts = []

    def _flush():
        if current_gene is None:
            return
        aligned = "".join(parts)
        pm = {}
        seq_pos = 0
        for col_idx, ch in enumerate(aligned, start=1):
            if ch != "-":
                seq_pos += 1
                pm[seq_pos] = col_idx
        pos_maps[current_gene] = pm

    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                _flush()
                current_gene = line[1:].split("|")[0]
                parts = []
            else:
                parts.append(line)
    _flush()
    return pos_maps


# ══════════════════════════════════════════════════════════════════════════════
# Step 1 – Read annotated phosphosites; split by route
# ══════════════════════════════════════════════════════════════════════════════
print("Reading annotated phosphosite file...")

# site_records: list of dicts with keys we will aggregate on
site_records = []   # {group_id, aligned_pos (int or sentinel), strain, protein_id, protein_pos (int), residue}

msa_needed = defaultdict(list)   # group_id -> rows needing MSA lookup
n_direct = n_msa = n_missing_group = 0

with open(ANNOT_FILE, encoding="utf-8", newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        ap = row["Aligned_pos"]
        gid = row["Group_ID"]

        if ap == "no_orthogroup":
            n_missing_group += 1
            continue

        rec = {
            "group_id":   gid,
            "strain":     row["Strain"],
            "protein_id": row["Protein"],
            "gene_id":    row["Gene_ID"],
            "protein_pos": int(row["Protein_pos"]),
            "residue":    row["PTM_residue"],
            "quality":    row["PTM_FLR_category"],
        }

        if ap == "needs_MSA":
            msa_needed[gid].append(rec)
            n_msa += 1
        else:
            rec["aligned_pos"] = int(ap)
            site_records.append(rec)
            n_direct += 1

print(f"  Direct (same-length):  {n_direct:,}")
print(f"  Needs MSA:             {n_msa:,} rows across {len(msa_needed):,} groups")
print(f"  No orthogroup:         {n_missing_group:,}")


# ══════════════════════════════════════════════════════════════════════════════
# Step 2 – Parse alignment files and map positions
# ══════════════════════════════════════════════════════════════════════════════
print(f"\nParsing alignment files in {ALIGN_DIR}/ ...")

n_mapped = n_align_missing = n_pos_missing = 0

for group_id, rows in msa_needed.items():
    afa_path = os.path.join(ALIGN_DIR, f"{group_id}.afa")

    if not os.path.exists(afa_path) or os.path.getsize(afa_path) == 0:
        # Alignment job didn't complete or produced empty output — flag for re-run
        for rec in rows:
            rec["aligned_pos"] = "alignment_missing"
            site_records.append(rec)
        n_align_missing += len(rows)
        continue

    pos_maps = parse_afa(afa_path)

    for rec in rows:
        gene_id    = rec["gene_id"]
        prot_pos   = rec["protein_pos"]
        pm = pos_maps.get(gene_id)

        if pm is None:
            # Gene not found in alignment (shouldn't happen, flag it)
            rec["aligned_pos"] = "gene_not_in_alignment"
            n_pos_missing += 1
        elif prot_pos not in pm:
            # Position out of range of alignment (shouldn't happen)
            rec["aligned_pos"] = "pos_out_of_range"
            n_pos_missing += 1
        else:
            rec["aligned_pos"] = pm[prot_pos]
            n_mapped += 1

        site_records.append(rec)

print(f"  MSA-mapped:            {n_mapped:,}")
print(f"  Alignment file missing: {n_align_missing:,} rows")
print(f"  Position lookup failed: {n_pos_missing:,} rows")


# ══════════════════════════════════════════════════════════════════════════════
# Step 3 – Aggregate by (group_id, aligned_pos)
# ══════════════════════════════════════════════════════════════════════════════
print("\nAggregating by (orthogroup, aligned_pos)...")

# site_map: (group_id, aligned_pos) -> {strain: [(protein_id, protein_pos, residue, quality)]}
site_map = defaultdict(lambda: {s: [] for s in STRAINS})

for rec in site_records:
    key = (rec["group_id"], rec["aligned_pos"])
    strain = rec["strain"]
    site_map[key][strain].append(
        (rec["protein_id"], rec["protein_pos"], rec["residue"], rec["quality"])
    )

print(f"  Unique (group, aligned_pos) entries: {len(site_map):,}")


# ══════════════════════════════════════════════════════════════════════════════
# Step 4 – Write output
# ══════════════════════════════════════════════════════════════════════════════
def format_sites(hits):
    """[('protId', pos, res, qual), ...] -> 'protId:pos(res)'  semicolon-joined"""
    return ";".join(f"{pid}:{pos}({res})" for pid, pos, res, _ in hits)

def best_quality(hits):
    order = {"Gold": 0, "Silver": 1, "Bronze": 2, "": 3}
    return min((h[3] for h in hits), key=lambda q: order.get(q, 99), default="")

FIELDS = (
    ["Group_ID", "Aligned_pos", "n_strains_with_site"]
    + [f"n_{s}" for s in STRAINS]
    + [f"{s}_sites" for s in STRAINS]
    + [f"{s}_best_FLR_cat" for s in STRAINS]
)

print(f"Writing {OUT_FILE} ...")

# Sort: numeric aligned_pos first, then string sentinels; within group keep genomic order
def sort_key(item):
    gid, apos = item[0]
    return (gid, 0 if isinstance(apos, int) else 1, apos if isinstance(apos, int) else 0)

with open(OUT_FILE, "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
    writer.writeheader()

    for (group_id, aligned_pos), strain_hits in sorted(site_map.items(), key=sort_key):
        n_strains = sum(1 for s in STRAINS if strain_hits[s])
        row = {
            "Group_ID":          group_id,
            "Aligned_pos":       aligned_pos,
            "n_strains_with_site": n_strains,
        }
        for s in STRAINS:
            hits = strain_hits[s]
            row[f"n_{s}"]           = len(hits)
            row[f"{s}_sites"]       = format_sites(hits)
            row[f"{s}_best_FLR_cat"] = best_quality(hits)
        writer.writerow(row)

total_rows = len(site_map)
sentinels  = sum(1 for (_, ap) in site_map if not isinstance(ap, int))
print(f"  Done. {total_rows:,} rows written ({sentinels} with unresolved positions).")


# ══════════════════════════════════════════════════════════════════════════════
# Step 5 – Quick summary
# ══════════════════════════════════════════════════════════════════════════════
from collections import Counter
strain_counts = Counter()
for (_, ap), hits in site_map.items():
    if not isinstance(ap, int):
        continue
    for s in STRAINS:
        if hits[s]:
            strain_counts[s] += 1

multi_strain = sum(
    1 for (_, ap), hits in site_map.items()
    if isinstance(ap, int) and sum(bool(hits[s]) for s in STRAINS) > 1
)
all_four = sum(
    1 for (_, ap), hits in site_map.items()
    if isinstance(ap, int) and all(hits[s] for s in STRAINS)
)

print("\n-- Site presence summary (resolved positions only) --")
for s in STRAINS:
    print(f"  {s}: {strain_counts[s]:,} aligned positions with detected sites")
print(f"  Found in >1 strain at same aligned position: {multi_strain:,}")
print(f"  Found in all 4 strains:                      {all_four:,}")
