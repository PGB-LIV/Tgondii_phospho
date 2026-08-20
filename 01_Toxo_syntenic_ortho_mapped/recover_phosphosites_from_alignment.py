"""
Recover "missed" phosphosites that are real biological sites but were not
detected by proteomics in a given strain (commonly ME49, the reference, which
is studied far less than GT1).

Rationale: at a given (Group_ID, Aligned_pos) where >=1 strain has a detected
phosphosite on S/T, sequence divergence near the site can mean the tryptic
peptide differs enough that the site simply isn't observed by MS in another
strain -- even though the S/T residue (and its defining motif context, e.g.
the +1 residue: SP, SD, etc.) is conserved in the alignment.

For every strain lacking a detected site at a given aligned position, recover
it if:
  1. the strain's residue at that aligned column is S or T, AND
  2. the strain's residue at the *next* aligned column (first non-gap column
     after Aligned_pos) matches the +1 residue seen in >=1 strain that DOES
     have a detected site there.

Reuses the (Group_ID, Aligned_pos) aggregation already computed by
build_aligned_site_table.py (read straight from phosphosites_crossgenome_sites_v3.tsv)
and pulls residues from the per-orthogroup alignment files -- both the
original MSA-derived ones (orthogroup_alignments/, diff-length groups) and
the newly generated ones for same-length groups
(orthogroup_alignments_same_length/), which are colinear FASTAs renamed .afa.
"""
import csv
import os
from collections import defaultdict

BASE        = os.path.dirname(os.path.abspath(__file__))
SITES_FILE  = f"{BASE}/phosphosites_crossgenome_sites_v3.tsv"
ORTHO_FILE  = f"{BASE}/syntenic_orthogroup_summary_v2.tsv"
ALIGN_DIRS  = [f"{BASE}/orthogroup_alignments", f"{BASE}/orthogroup_alignments_same_length"]
OUT_FILE    = f"{BASE}/phosphosites_crossgenome_sites_v4_recovered.tsv"

STRAINS = ["ME49", "GT1", "RH88", "VEG"]
STP_RESIDUES = {"S", "T"}


# ══════════════════════════════════════════════════════════════════════════════
# 1. Load orthogroup membership: Group_ID -> {strain: [gene_id, ...]}
#    (needed to know which gene_id in the alignment belongs to which strain,
#    since a strain can have 0, 1, or >1 members in an orthogroup)
# ══════════════════════════════════════════════════════════════════════════════
print("Loading orthogroup membership...")

group_strain_genes = defaultdict(lambda: defaultdict(list))  # group_id -> strain -> [gene_id]

with open(ORTHO_FILE, encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        group_id = row["Group_ID"]
        for strain in STRAINS:
            for gid in row[strain].split(";"):
                gid = gid.strip()
                if gid:
                    group_strain_genes[group_id][strain].append(gid)

print(f"  {len(group_strain_genes):,} orthogroups")


# ══════════════════════════════════════════════════════════════════════════════
# 2. Locate alignment file per group (search both directories)
# ══════════════════════════════════════════════════════════════════════════════
align_path_for_group = {}
for d in ALIGN_DIRS:
    if not os.path.isdir(d):
        continue
    for fname in os.listdir(d):
        if fname.endswith(".afa"):
            align_path_for_group[fname[:-4]] = os.path.join(d, fname)

print(f"  {len(align_path_for_group):,} alignment files found across {len(ALIGN_DIRS)} directories")


# ══════════════════════════════════════════════════════════════════════════════
# 3. Parse an .afa -> {gene_id: aligned_seq_string}  (cache one group at a time)
# ══════════════════════════════════════════════════════════════════════════════
def parse_afa_seqs(path):
    seqs = {}
    current_gene = None
    parts = []

    def _flush():
        if current_gene is not None:
            seqs[current_gene] = "".join(parts)

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
    return seqs


def residue_at(aligned_seq, col_1based):
    """1-based alignment column -> residue char, or '' if out of range/gap-only tail."""
    idx = col_1based - 1
    if idx < 0 or idx >= len(aligned_seq):
        return ""
    return aligned_seq[idx]


def next_residue_after(aligned_seq, col_1based):
    """First non-gap residue strictly after the given 1-based column."""
    idx = col_1based  # 0-based index of the next column
    while idx < len(aligned_seq):
        if aligned_seq[idx] != "-":
            return aligned_seq[idx]
        idx += 1
    return ""


# ══════════════════════════════════════════════════════════════════════════════
# 4. Walk the existing cross-genome site table and attempt recovery for each
#    strain missing a site at rows with a numeric Aligned_pos.
# ══════════════════════════════════════════════════════════════════════════════
print(f"\nReading {SITES_FILE} ...")

OUT_FIELDS = [
    "Group_ID", "Aligned_pos", "n_strains_with_site_orig",
    "n_strains_with_site_recovered",
] + [f"n_{s}_orig" for s in STRAINS] + [f"{s}_recovered" for s in STRAINS] + [
    "n_newly_recovered", "Recovered_residue", "Recovered_next_residue",
]

n_rows = 0
n_candidate_rows = 0   # rows with >=1 site and >=1 strain missing, numeric Aligned_pos
n_recovered_events = 0
n_alignment_missing = 0
n_no_genes_in_alignment = 0

recovery_rows = []

with open(SITES_FILE, encoding="utf-8", newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        n_rows += 1
        gid = row["Group_ID"]
        ap_raw = row["Aligned_pos"]

        try:
            aligned_pos = int(ap_raw)
        except ValueError:
            continue  # skip sentinel rows (alignment_missing, pos_out_of_range, etc.)

        has_site = {s: int(row[f"n_{s}"]) > 0 for s in STRAINS}
        missing_strains = [s for s in STRAINS if not has_site[s]]
        present_strains = [s for s in STRAINS if has_site[s]]

        if not missing_strains or not present_strains:
            continue  # nothing to recover, or no site at all at this position

        n_candidate_rows += 1

        afa_path = align_path_for_group.get(gid)
        if afa_path is None:
            n_alignment_missing += 1
            continue

        seqs = parse_afa_seqs(afa_path)

        # Determine the +1 motif residue(s) seen among strains with a detected site
        present_next_residues = set()
        present_residues = set()
        for s in present_strains:
            for gene_id in group_strain_genes[gid].get(s, []):
                seq = seqs.get(gene_id)
                if seq is None:
                    continue
                present_residues.add(residue_at(seq, aligned_pos))
                present_next_residues.add(next_residue_after(seq, aligned_pos))
        present_residues.discard("")
        present_next_residues.discard("")

        if not present_residues:
            n_no_genes_in_alignment += 1
            continue

        recovered = {}
        for s in missing_strains:
            recovered_here = False
            res_seen = next_seen = ""
            for gene_id in group_strain_genes[gid].get(s, []):
                seq = seqs.get(gene_id)
                if seq is None:
                    continue
                res = residue_at(seq, aligned_pos)
                nxt = next_residue_after(seq, aligned_pos)
                if res in STP_RESIDUES and nxt in present_next_residues:
                    recovered_here = True
                    res_seen, next_seen = res, nxt
                    break
            recovered[s] = (recovered_here, res_seen, next_seen)

        n_new = sum(1 for s in missing_strains if recovered[s][0])
        if n_new == 0:
            continue

        n_recovered_events += n_new

        out_row = {
            "Group_ID": gid,
            "Aligned_pos": aligned_pos,
            "n_strains_with_site_orig": row["n_strains_with_site"],
            "n_strains_with_site_recovered": int(row["n_strains_with_site"]) + n_new,
            "n_newly_recovered": n_new,
        }
        rec_residues = set()
        rec_next = set()
        for s in STRAINS:
            out_row[f"n_{s}_orig"] = row[f"n_{s}"]
            if s in missing_strains and recovered[s][0]:
                out_row[f"{s}_recovered"] = "Yes"
                rec_residues.add(recovered[s][1])
                rec_next.add(recovered[s][2])
            else:
                out_row[f"{s}_recovered"] = ""
        out_row["Recovered_residue"] = "/".join(sorted(rec_residues))
        out_row["Recovered_next_residue"] = "/".join(sorted(rec_next))
        recovery_rows.append(out_row)

print(f"  Total rows in input table:                 {n_rows:,}")
print(f"  Candidate rows (site + missing strain):    {n_candidate_rows:,}")
print(f"  Alignment file missing for group:          {n_alignment_missing:,}")
print(f"  No member genes found in alignment:        {n_no_genes_in_alignment:,}")
print(f"  Rows with >=1 newly recovered site:         {len(recovery_rows):,}")
print(f"  Total newly recovered (group,pos,strain):   {n_recovered_events:,}")


# ══════════════════════════════════════════════════════════════════════════════
# 5. Write output
# ══════════════════════════════════════════════════════════════════════════════
print(f"\nWriting {OUT_FILE} ...")
with open(OUT_FILE, "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=OUT_FIELDS, delimiter="\t")
    writer.writeheader()
    for r in recovery_rows:
        writer.writerow(r)

print(f"  Done. {len(recovery_rows):,} rows written.")


# ══════════════════════════════════════════════════════════════════════════════
# 6. Per-strain recovery summary
# ══════════════════════════════════════════════════════════════════════════════
print("\n-- Per-strain newly recovered sites --")
strain_recovered_counts = defaultdict(int)
for r in recovery_rows:
    for s in STRAINS:
        if r[f"{s}_recovered"] == "Yes":
            strain_recovered_counts[s] += 1

for s in STRAINS:
    print(f"  {s}: {strain_recovered_counts[s]:,} newly recovered sites")
