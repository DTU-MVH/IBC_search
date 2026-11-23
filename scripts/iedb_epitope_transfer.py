"""Epitope-vs-epitope retrieval: do IEDB epitope motifs find OTHER epitopes?

The question this answers, which the count experiments did not
--------------------------------------------------------------
Every earlier run scored a ball by HOW MANY sites it matched. That is a
statement about commonness. It cannot distinguish "this geometry is rare"
from "this geometry matches the right things". This asks the second question:

    do motifs drawn from a curated IEDB epitope retrieve EPITOPE residues in
    other, unrelated antigens more often than motifs drawn from the same
    protein's non-epitope surface do?

Design
------
QUERIES  the 70 IEDB-3D curated antigens, scanned with ibc-search scan at
         5.5 A. A motif is *epitope* when >= POS_COV of its labelled residues
         are curated epitope contacts, *non_epitope* when none are, and is
         discarded when it straddles the boundary. Both classes come from the
         same protein, so fold, size and composition cancel.

TARGET   the sabdab2 antigen index: 5,951 antibody-antigen complexes with
         per-residue epitope labels and SASA. The label on the target side is
         structural contact, independent of the IEDB curation on the query
         side, so query labels and target labels do not share a source.

SCORE    of the matched residues that are EXPOSED in the hit antigen (a buried
         residue cannot be an antibody epitope, so counting it only dilutes),
         what fraction are labelled epitope. Compared against the base rate of
         epitope among exposed residues in the index.

Leakage control
---------------
A hit is discarded when the target antigen
  * comes from the same PDB entry as the query,
  * belongs to an ag_cluster that any same-PDB target belongs to, or
  * shares a 9-mer of sequence with the query antigen.
Without this an epitope motif retrieves a near-copy of its own antigen and
finds its own epitope, which is a guaranteed and meaningless positive.
"""
from __future__ import annotations
import csv, os, subprocess, sys, time, statistics as st
from collections import defaultdict
from pathlib import Path

REPO = Path("/Users/mvh_mac/documents_local_only/IBC_search")
IEDB = REPO / "investigation_databases/reveived_Ag_target_annotated_proteins"
SAB = Path.home() / "documents_local_only/master_thesis_data_files/sabdab2_antigen_index"
INDEX = SAB / "index/sabdab2_antigen"
WORK = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/iedb_transfer")
RADIUS, TOP, POS_COV, KMER = 5.5, 50, 0.5, 9
WORK.mkdir(parents=True, exist_ok=True)

THREE2ONE = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}


def tsv(p):
    with Path(p).open() as fh:
        yield from csv.DictReader(fh, delimiter="\t")


# ---------------------------------------------------------------- phase 0
print("loading labels", flush=True)
sab_struct = {r["instance"]: r for r in tsv(SAB / "structures.tsv")}

sab_epi: dict[tuple[str, str, int], int] = {}
sab_seq: dict[str, list[str]] = defaultdict(list)
for r in tsv(SAB / "epitope_residues.tsv"):
    try:
        sab_epi[(r["instance"], r["chain"], int(r["res_id"]))] = int(r["is_epitope"])
    except ValueError:
        continue
    sab_seq[r["instance"]].append(THREE2ONE.get(r["res_name"], "X"))

sab_exp: dict[tuple[str, str, int], int] = {}
n_exposed_tot = n_exposed_epi = 0
for r in tsv(SAB / "residue_sasa.tsv"):
    try:
        k = (r["instance"], r["chain"], int(r["res_id"]))
    except ValueError:
        continue
    e = int(r["is_exposed"])
    sab_exp[k] = e
    if e:
        n_exposed_tot += 1
        n_exposed_epi += sab_epi.get(k, 0)
BASE = n_exposed_epi / n_exposed_tot
print(f"  sabdab: {len(sab_struct)} antigens, {len(sab_epi):,} residues; "
      f"base epitope rate among EXPOSED residues = {BASE:.4f}", flush=True)

# query-side labels
iedb_epi: dict[tuple[str, str, int], int] = {}
iedb_seq: dict[str, list[str]] = defaultdict(list)
for r in tsv(IEDB / "residue_labels.tsv"):
    try:
        iedb_epi[(r["instance"], r["chain"], int(r["res_id"]))] = int(r["is_annotated"])
    except ValueError:
        continue
    iedb_seq[r["instance"]].append(THREE2ONE.get(r["res_name"], "X"))
iedb_struct = {r["instance"]: r for r in tsv(IEDB / "structures.tsv")}
print(f"  iedb  : {len(iedb_struct)} antigens, {len(iedb_epi):,} labelled residues, "
      f"{sum(iedb_epi.values()):,} epitope", flush=True)

# sequence k-mer sets, for the homology screen
def kmers(seq: list[str]) -> set[str]:
    s = "".join(seq)
    return {s[i:i + KMER] for i in range(len(s) - KMER + 1)}

sab_kmers = {i: kmers(s) for i, s in sab_seq.items()}
pdb_to_inst = defaultdict(list)
for i, r in sab_struct.items():
    pdb_to_inst[r["pdb_id"]].append(i)
inst_cluster = {i: {c for c in (r["ag_cluster"] or "").split("/") if c}
                for i, r in sab_struct.items()}

# ---------------------------------------------------------------- phase 1
rows_out = []
insts = sorted(iedb_struct)
t0 = time.time()
for n, inst in enumerate(insts, 1):
    cif = IEDB / "structures" / f"{inst}.cif"
    out = WORK / inst
    env = {**os.environ, "PYTHONPATH": str(REPO / "src")}
    if not (out / "folddisco_queries.tsv").is_file():
        sc = subprocess.run([str(REPO / ".venv/bin/python"), "-m", "ibc_search",
                             "scan", str(cif), "--radius", str(RADIUS),
                             "-o", str(out), "--log-level", "ERROR"],
                            capture_output=True, text=True, env=env)
        if sc.returncode != 0 or not (out / "folddisco_queries.tsv").is_file():
            print(f"  [{n}/70] {inst}: scan failed", flush=True); continue
    hits = out / "hits.tsv"
    if not hits.is_file():
        with hits.open("w") as fh:
            q = subprocess.run(["folddisco", "query", "-i", str(INDEX), "-q",
                                str(out / "folddisco_queries.tsv"), "-t", "1",
                                "--top", str(TOP), "--header", "--per-match",
                                "--format-output",
                                "tid,idf,rmsd,query_residues,matching_residues"],
                               stdout=fh, stderr=subprocess.PIPE, text=True)
        if q.returncode != 0:
            print(f"  [{n}/70] {inst}: query failed", flush=True); continue

    # motif -> residues
    motif_res = defaultdict(set)
    for r in tsv(out / "patch_members.tsv"):
        motif_res[r["patch_id"]].add((r["auth_chain_id"], int(r["auth_seq_id"])))
    qmap = {r["query_string"]: motif_res.get(r["patch_id"], set())
            for r in tsv(out / "folddisco_queries.index.tsv")}

    # leakage set for this antigen
    own_pdb = f"pdb_0000{iedb_struct[inst]['pdb']}"
    bad_clusters = set()
    for i in pdb_to_inst.get(own_pdb, []):
        bad_clusters |= inst_cluster[i]
    qk = kmers(iedb_seq[inst])
    excluded = set()
    for i, r in sab_struct.items():
        if r["pdb_id"] == own_pdb or (bad_clusters & inst_cluster[i]) \
           or (qk & sab_kmers.get(i, set())):
            excluded.add(i)

    agg = defaultdict(lambda: [0, 0, 0])      # hits, exposed matched, epitope
    n_rows = n_drop = 0
    for line in hits.open():
        f = line.rstrip("\n").split("\t")
        if len(f) != 5 or f[0] == "tid":
            continue
        n_rows += 1
        hinst = f[0].split("/")[-1].removesuffix(".cif")
        if hinst not in sab_struct or hinst in excluded:
            n_drop += 1
            continue
        a = agg[f[3]]
        a[0] += 1
        for t in f[4].split(","):
            if not t or t == "_":
                continue
            try:
                key = (hinst, t[0], int(t[1:]))
            except ValueError:
                continue
            if not sab_exp.get(key, 0):
                continue
            a[1] += 1
            a[2] += sab_epi.get(key, 0)

    for qs, residues in qmap.items():
        known = [iedb_epi[(inst, c, r)] for c, r in residues
                 if (inst, c, r) in iedb_epi]
        if not known:
            continue
        cov = sum(known) / len(known)
        cls = "epitope" if cov >= POS_COV else ("non_epitope" if cov == 0 else "straddle")
        h, ex, ep = agg.get(qs, [0, 0, 0])
        rows_out.append({"antigen": inst, "query": qs, "size": len(residues),
                         "cls": cls, "cov": cov, "n_hits": h,
                         "n_exposed": ex, "n_epitope": ep})
    if n % 10 == 0 or n == len(insts):
        print(f"  [{n}/{len(insts)}] {time.time()-t0:.0f}s  "
              f"(last: {n_rows} rows, {n_drop} dropped for leakage)", flush=True)

with (WORK / "transfer.tsv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows_out[0]), delimiter="\t",
                       lineterminator="\n")
    w.writeheader(); w.writerows(rows_out)
print(f"\nwrote {WORK/'transfer.tsv'}  ({len(rows_out)} motifs)")
