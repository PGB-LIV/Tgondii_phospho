"""
Build syntenic ortholog group summary from four Toxoplasma strain files.

Strategy: connected components on the pairwise syntenic links in
'Input Ortholog(s)', NOT grouping by the broad OG family ID.
This correctly separates paralogs that share an OG but have distinct
syntenic counterparts (e.g. MIC17A/B/C all in OG6_134570 but with
different RH88 partners).
"""
import csv
import os
from collections import defaultdict, deque

BASE = os.path.dirname(os.path.abspath(__file__))

FILES = {
    "ME49": f"{BASE}/AllME49_genes_from_other_strains.tsv",
    "GT1":  f"{BASE}/AllGT1_genes_from_other_strains.tsv",
    "RH88": f"{BASE}/AllRH88_genes_from_other_strains.tsv",
    "VEG":  f"{BASE}/AllVEG_genes_from_other_strains.tsv",
}

STRAIN_PREFIXES = {
    "TGME49": "ME49",
    "TGGT1":  "GT1",
    "TGRH88": "RH88",
    "TGVEG":  "VEG",
}

def get_strain(gene_id):
    for prefix, strain in STRAIN_PREFIXES.items():
        if gene_id.startswith(prefix):
            return strain
    return None   # human/mouse/other — ignore


# ── Step 1: read all files ───────────────────────────────────────────────────
# gene_info: gene_id -> {strain, og, description}
# adj:       undirected adjacency list (Toxo genes only)

gene_info = {}
adj = defaultdict(set)

for strain, path in FILES.items():
    seen = set()
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_id = row["Gene ID"].strip()
            og      = row["Ortholog Group"].strip()
            desc    = row["Product Description"].strip()
            ortho   = row["Input Ortholog(s)"].strip()

            if gene_id not in gene_info:
                gene_info[gene_id] = {"strain": strain, "og": og, "description": desc}

            if gene_id in seen:
                continue
            seen.add(gene_id)

            for partner in ortho.split(","):
                partner = partner.strip()
                if partner and get_strain(partner) is not None:
                    adj[gene_id].add(partner)
                    adj[partner].add(gene_id)   # ensure reverse edge too


# ── Step 2: connected components (BFS) ──────────────────────────────────────
# Universe = all genes that appear in any file OR are named as an ortholog
all_nodes = set(gene_info.keys()) | set(adj.keys())

visited    = set()
components = []

for start in sorted(all_nodes):
    if start in visited:
        continue
    component = []
    queue = deque([start])
    visited.add(start)
    while queue:
        node = queue.popleft()
        component.append(node)
        for nbr in adj.get(node, ()):
            if nbr not in visited:
                visited.add(nbr)
                queue.append(nbr)
    components.append(component)


# ── Step 3: assemble output rows ─────────────────────────────────────────────
def best_desc(descs):
    skip = {"hypothetical protein", "unspecified product", ""}
    informative = sorted(d for d in descs if d.lower() not in skip)
    return informative[0] if informative else (sorted(descs)[0] if descs else "")


rows = []
for component in components:
    members = {"ME49": [], "GT1": [], "RH88": [], "VEG": []}
    ogs   = set()
    descs = set()

    for gid in component:
        info   = gene_info.get(gid)
        strain = (info["strain"] if info else None) or get_strain(gid)
        if strain in members:
            members[strain].append(gid)
        if info:
            if info["og"]:
                ogs.add(info["og"])
            if info["description"]:
                descs.add(info["description"])

    for s in members:
        members[s].sort()

    me49, gt1, rh88, veg = members["ME49"], members["GT1"], members["RH88"], members["VEG"]

    rep      = (me49 or gt1 or rh88 or veg)[0]
    group_id = f"Tg_ortho_{rep}"

    counts   = [len(me49), len(gt1), len(rh88), len(veg)]
    n_strains = sum(1 for c in counts if c > 0)

    rows.append({
        "Group_ID":          group_id,
        "Ortholog_Group":    ";".join(sorted(ogs)),
        "ME49":              ";".join(me49),
        "GT1":               ";".join(gt1),
        "RH88":              ";".join(rh88),
        "VEG":               ";".join(veg),
        "n_ME49":            len(me49),
        "n_GT1":             len(gt1),
        "n_RH88":            len(rh88),
        "n_VEG":             len(veg),
        "n_strains_present": n_strains,
        "Complex_group":     "Yes" if max(counts) > 1 else "No",
        "Product_Description": best_desc(descs),
    })

rows.sort(key=lambda r: r["Group_ID"])


# ── Step 4: write TSV ────────────────────────────────────────────────────────
out_path = f"{BASE}/syntenic_orthogroup_summary_v2.tsv"
FIELDS   = [
    "Group_ID", "Ortholog_Group",
    "ME49", "GT1", "RH88", "VEG",
    "n_ME49", "n_GT1", "n_RH88", "n_VEG",
    "n_strains_present", "Complex_group",
    "Product_Description",
]

with open(out_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)


# ── Step 5: summary stats ────────────────────────────────────────────────────
total        = len(rows)
complex_grps = sum(1 for r in rows if r["Complex_group"] == "Yes")
all_four     = sum(1 for r in rows if r["n_strains_present"] == 4)
missing_me49 = sum(1 for r in rows if r["n_ME49"] == 0)
only_one     = sum(1 for r in rows if r["n_strains_present"] == 1)
multi_og     = sum(1 for r in rows if ";" in r["Ortholog_Group"])

print(f"Output: {out_path}")
print(f"  Total syntenic groups     : {total:,}")
print(f"  Present in all 4 strains  : {all_four:,}  ({100*all_four/total:.1f}%)")
print(f"  Missing from ME49         : {missing_me49:,}")
print(f"  Present in only 1 strain  : {only_one:,}")
print(f"  Complex (any strain >=2)  : {complex_grps:,}")
print(f"  Groups spanning >1 OG     : {multi_og:,}")
