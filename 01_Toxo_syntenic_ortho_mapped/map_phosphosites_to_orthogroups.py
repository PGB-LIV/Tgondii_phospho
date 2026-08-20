"""
Map phosphosites to aligned orthogroup positions.

For orthogroups where every member protein has the same length:
  aligned_pos = protein_pos  (no indels assumed)

For orthogroups where lengths differ:
  write a multi-sequence FASTA to orthogroup_fastas_for_muscle/
  so MUSCLE can produce the alignment for a subsequent mapping step.
"""
import csv
import os
import re
import sys
from collections import defaultdict, Counter

if sys.stdout.encoding is not None and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

BASE        = os.path.dirname(os.path.abspath(__file__))
FASTA_DIR   = f"{BASE}/protein_fasta_files"
PHOSPHO_IN  = f"{BASE}/G2S1B_0.05_protein_pos_all_prot_mapping_without_contam.csv"
ORTHO_FILE  = f"{BASE}/syntenic_orthogroup_summary_v2.tsv"
OUT_ANNOT   = f"{BASE}/phosphosites_orthogroup_annotated.tsv"
MUSCLE_DIR  = f"{BASE}/orthogroup_fastas_for_muscle"

FASTA_FILES = {
    "ME49": f"{FASTA_DIR}/ToxoDB-68_TgondiiME49_AnnotatedProteins.fasta",
    "GT1":  f"{FASTA_DIR}/ToxoDB-68_TgondiiGT1_AnnotatedProteins.fasta",
    "RH88": f"{FASTA_DIR}/ToxoDB-68_TgondiiRH88_AnnotatedProteins.fasta",
    "VEG":  f"{FASTA_DIR}/ToxoDB-68_TgondiiVEG_AnnotatedProteins.fasta",
}

TOXO_STRAINS = {"ME49", "GT1", "RH88", "VEG"}


# ══════════════════════════════════════════════════════════════════════════════
# 1. Parse FASTA files
#    protein_db:       protein_id -> {gene_id, strain, length, seq}
#    gene_to_proteins: gene_id    -> [protein_id, ...]   (ordered)
# ══════════════════════════════════════════════════════════════════════════════
print("Parsing FASTA files...")

protein_db       = {}                       # protein_id -> dict
gene_to_proteins = defaultdict(list)        # gene_id    -> [protein_id]

def parse_fasta(path, strain):
    prot_id = info = None
    seq_parts = []

    def _flush():
        if prot_id is None:
            return
        seq = "".join(seq_parts)
        info["seq"]    = seq
        info["length"] = len(seq)
        protein_db[prot_id] = info
        gene_to_proteins[info["gene_id"]].append(prot_id)

    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                _flush()
                header  = line[1:]
                prot_id = header.split("|")[0].strip()
                m       = re.search(r"\bgene=(\S+)", header)
                gene_id = m.group(1) if m else prot_id
                info      = {"gene_id": gene_id, "strain": strain}
                seq_parts = []
            else:
                seq_parts.append(line)
    _flush()

for strain, path in FASTA_FILES.items():
    parse_fasta(path, strain)
    n_genes = len(set(protein_db[p]["gene_id"] for p in protein_db
                      if protein_db[p]["strain"] == strain))
    print(f"  {strain}: {n_genes:,} genes, "
          f"{sum(1 for p in protein_db if protein_db[p]['strain']==strain):,} proteins")

# Detect multi-isoform genes
multi_iso = {g: ps for g, ps in gene_to_proteins.items() if len(ps) > 1}
print(f"  Genes with >1 protein isoform: {len(multi_iso)}")


# ══════════════════════════════════════════════════════════════════════════════
# 2. Load orthogroup summary
#    gene_to_group:   gene_id  -> Group_ID
#    group_members:   Group_ID -> [(gene_id, strain)]   (all members)
# ══════════════════════════════════════════════════════════════════════════════
print("\nLoading orthogroup summary...")

gene_to_group  = {}                    # gene_id  -> Group_ID
group_members  = defaultdict(list)     # Group_ID -> [(gene_id, strain)]

with open(ORTHO_FILE, encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        group_id = row["Group_ID"]
        for strain in TOXO_STRAINS:
            for gid in row[strain].split(";"):
                gid = gid.strip()
                if gid:
                    gene_to_group[gid] = group_id
                    group_members[group_id].append((gid, strain))

print(f"  {len(group_members):,} orthogroups, {len(gene_to_group):,} gene->group mappings")


# ══════════════════════════════════════════════════════════════════════════════
# 3. Read phosphosite file — keep Toxo rows only
#    For each Toxo protein: resolve gene_id and group_id
# ══════════════════════════════════════════════════════════════════════════════
print("\nReading phosphosite file...")

phospho_rows = []          # all Toxo phosphosite rows with annotations added
phospho_header = None

# Track which protein_id is used per gene (to pick the best representative
# for FASTA construction — most-used protein isoform per gene)
gene_phospho_prot_counts = defaultdict(Counter)  # gene_id -> Counter(protein_id)

with open(PHOSPHO_IN, encoding="utf-8", newline="") as fh:
    reader = csv.DictReader(fh)
    phospho_header = reader.fieldnames

    for row in reader:
        prot_id = row["Protein"]

        # Skip non-Toxo (human, mouse sp| tr| etc.)
        if prot_id not in protein_db:
            continue

        info    = protein_db[prot_id]
        strain  = info["strain"]
        gene_id = info["gene_id"]
        group_id = gene_to_group.get(gene_id, "")

        gene_phospho_prot_counts[gene_id][prot_id] += 1

        phospho_rows.append({
            **row,
            "_strain":    strain,
            "_gene_id":   gene_id,
            "_group_id":  group_id,
            "_prot_len":  info["length"],
        })

print(f"  Toxo phosphosite rows: {len(phospho_rows):,}")
print(f"  Rows with no orthogroup: "
      f"{sum(1 for r in phospho_rows if not r['_group_id']):,}")


# ══════════════════════════════════════════════════════════════════════════════
# 4. For each orthogroup, pick the representative protein per gene:
#      (a) the most-used phosphosite protein if the gene has phosphosites
#      (b) else the first protein in gene_to_proteins
#    Then collect all representative lengths to decide same-length / needs-MSA.
# ══════════════════════════════════════════════════════════════════════════════
print("\nDetermining representative proteins and length homogeneity...")

def rep_protein(gene_id):
    """Best single representative protein for a gene."""
    if gene_id in gene_phospho_prot_counts:
        return gene_phospho_prot_counts[gene_id].most_common(1)[0][0]
    prots = gene_to_proteins.get(gene_id, [])
    return prots[0] if prots else None

# group_id -> list of (gene_id, strain, prot_id, length, seq)
# Only include genes whose representative protein is in protein_db
group_rep_proteins = {}

for group_id, members in group_members.items():
    entries = []
    for gene_id, strain in members:
        pid = rep_protein(gene_id)
        if pid is None or pid not in protein_db:
            continue
        info = protein_db[pid]
        entries.append((gene_id, strain, pid, info["length"], info["seq"]))
    group_rep_proteins[group_id] = entries

def group_is_same_length(group_id):
    entries = group_rep_proteins.get(group_id, [])
    if not entries:
        return True   # no proteins to compare — treat as same
    lengths = {length for _, _, _, length, _ in entries}
    return len(lengths) == 1

# Pre-compute for all groups that have phosphosites
groups_with_phosphosites = set(r["_group_id"] for r in phospho_rows if r["_group_id"])

same_length_groups  = set()
diff_length_groups  = set()
for gid in groups_with_phosphosites:
    if group_is_same_length(gid):
        same_length_groups.add(gid)
    else:
        diff_length_groups.add(gid)

print(f"  Groups with phosphosites:   {len(groups_with_phosphosites):,}")
print(f"    Same-length (direct map): {len(same_length_groups):,}")
print(f"    Diff-length (needs MSA):  {len(diff_length_groups):,}")


# ══════════════════════════════════════════════════════════════════════════════
# 5. Write per-orthogroup FASTA files for diff-length groups
# ══════════════════════════════════════════════════════════════════════════════
print(f"\nWriting FASTA files to {MUSCLE_DIR}/...")
os.makedirs(MUSCLE_DIR, exist_ok=True)

for group_id in sorted(diff_length_groups):
    entries = group_rep_proteins.get(group_id, [])
    if not entries:
        continue
    fasta_path = os.path.join(MUSCLE_DIR, f"{group_id}.fasta")
    with open(fasta_path, "w") as fh:
        for gene_id, strain, prot_id, length, seq in entries:
            fh.write(f">{gene_id}|{prot_id}|{strain}|len={length}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")

print(f"  Written {len(diff_length_groups):,} FASTA files")


# ══════════════════════════════════════════════════════════════════════════════
# 6. Write annotated phosphosite table
# ══════════════════════════════════════════════════════════════════════════════
print(f"\nWriting annotated table to {OUT_ANNOT}...")

NEW_COLS = ["Strain", "Gene_ID", "Group_ID", "Protein_length",
            "All_group_same_length", "Aligned_pos"]

with open(OUT_ANNOT, "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=NEW_COLS + phospho_header,
        delimiter="\t",
        extrasaction="ignore",
    )
    writer.writeheader()

    for row in phospho_rows:
        group_id = row["_group_id"]
        same_len = group_id in same_length_groups if group_id else None

        if group_id and same_len:
            aligned_pos = row["Protein_pos"]   # same length → direct mapping
        elif group_id and not same_len:
            aligned_pos = "needs_MSA"
        else:
            aligned_pos = "no_orthogroup"

        out = {
            "Strain":               row["_strain"],
            "Gene_ID":              row["_gene_id"],
            "Group_ID":             group_id,
            "Protein_length":       row["_prot_len"],
            "All_group_same_length": "Yes" if same_len else ("No" if same_len is False else ""),
            "Aligned_pos":          aligned_pos,
        }
        out.update({k: row[k] for k in phospho_header})
        writer.writerow(out)

print(f"  Done. {len(phospho_rows):,} rows written.")


# ══════════════════════════════════════════════════════════════════════════════
# 7. Summary statistics
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Summary ────────────────────────────────────────────────")
total = len(phospho_rows)
n_mapped   = sum(1 for r in phospho_rows if r["_group_id"] in same_length_groups)
n_msa      = sum(1 for r in phospho_rows if r["_group_id"] in diff_length_groups)
n_nogroup  = sum(1 for r in phospho_rows if not r["_group_id"])

print(f"  Total Toxo phosphosite rows   : {total:,}")
print(f"  Directly aligned (same length): {n_mapped:,}  ({100*n_mapped/total:.1f}%)")
print(f"  Flagged for MSA (diff length) : {n_msa:,}  ({100*n_msa/total:.1f}%)")
print(f"  No orthogroup assigned        : {n_nogroup:,}  ({100*n_nogroup/total:.1f}%)")
print(f"  Multi-isoform genes flagged   : {len(multi_iso):,}")
