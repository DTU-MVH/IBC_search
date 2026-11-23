"""Statistics for the epitope-vs-epitope retrieval test."""
import csv, sys, statistics as st
from collections import defaultdict
from pathlib import Path

W = Path(sys.argv[1]); BASE = float(sys.argv[2])
rows = list(csv.DictReader((W / "transfer.tsv").open(), delimiter="\t"))
for r in rows:
    for k in ("size", "n_hits", "n_exposed", "n_epitope"):
        r[k] = int(r[k])
    r["cov"] = float(r["cov"])

print("=" * 74)
print("MOTIF CLASSES")
print("=" * 74)
for c in ("epitope", "non_epitope", "straddle"):
    g = [r for r in rows if r["cls"] == c]
    live = [r for r in g if r["n_exposed"] > 0]
    print(f"  {c:12s} {len(g):6d} motifs   {len(live):6d} with >=1 exposed matched residue")

print()
print("=" * 74)
print("POOLED TRANSFER RATE   (epitope share of matched exposed residues)")
print("=" * 74)
print(f"  {'class':12s} {'matched exposed':>16s} {'of which epitope':>18s} {'rate':>8s} {'vs base':>9s}")
pooled = {}
for c in ("epitope", "non_epitope"):
    g = [r for r in rows if r["cls"] == c]
    ex = sum(r["n_exposed"] for r in g); ep = sum(r["n_epitope"] for r in g)
    pooled[c] = ep / ex if ex else 0.0
    print(f"  {c:12s} {ex:16,d} {ep:18,d} {pooled[c]:8.4f} {pooled[c]/BASE:8.2f}x")
print(f"  {'base rate':12s} {'':16s} {'':18s} {BASE:8.4f} {1.0:8.2f}x")
print(f"\n  epitope / non_epitope ratio: {pooled['epitope']/pooled['non_epitope']:.3f}")

# per-antigen paired comparison
print()
print("=" * 74)
print("PAIRED WITHIN ANTIGEN  (each protein contributes both classes)")
print("=" * 74)
by = defaultdict(lambda: {"epitope": [0, 0], "non_epitope": [0, 0]})
for r in rows:
    if r["cls"] in ("epitope", "non_epitope"):
        by[r["antigen"]][r["cls"]][0] += r["n_exposed"]
        by[r["antigen"]][r["cls"]][1] += r["n_epitope"]
diffs = []
for a, d in by.items():
    e, ne = d["epitope"], d["non_epitope"]
    if e[0] == 0 or ne[0] == 0:
        continue
    diffs.append((a, e[1] / e[0], ne[1] / ne[0]))
up = sum(1 for _, x, y in diffs if x > y)
dn = sum(1 for _, x, y in diffs if x < y)
tie = len(diffs) - up - dn
print(f"  antigens with both classes populated : {len(diffs)}")
print(f"  epitope rate HIGHER than non-epitope : {up}")
print(f"  epitope rate LOWER                   : {dn}")
print(f"  equal                                : {tie}")
if diffs:
    d = [x - y for _, x, y in diffs]
    print(f"  mean paired difference               : {st.mean(d):+.4f}"
          f"  (median {st.median(d):+.4f}, SD {st.stdev(d) if len(d)>1 else 0:.4f})")
    try:
        from scipy.stats import wilcoxon, binomtest
        print(f"  Wilcoxon signed-rank p               : {wilcoxon(d).pvalue:.4g}")
        print(f"  sign test p                          : {binomtest(up, up+dn).pvalue:.4g}")
    except Exception as exc:
        print(f"  (scipy unavailable: {exc})")

# motif-level AUC
print()
print("=" * 74)
print("MOTIF-LEVEL AUC  (can a motif's transfer rate tell the classes apart?)")
print("=" * 74)
a = [r["n_epitope"] / r["n_exposed"] for r in rows
     if r["cls"] == "epitope" and r["n_exposed"] > 0]
b = [r["n_epitope"] / r["n_exposed"] for r in rows
     if r["cls"] == "non_epitope" and r["n_exposed"] > 0]
if a and b:
    m = sorted([(v, 1) for v in a] + [(v, 0) for v in b])
    i, n, cum = 0, len(m), 0.0
    while i < n:
        j = i
        while j < n and m[j][0] == m[i][0]:
            j += 1
        avg = (i + j + 1) / 2.0
        cum += sum(avg for k in range(i, j) if m[k][1] == 1)
        i = j
    na, nb = len(a), len(b)
    auc = (cum - na * (na + 1) / 2.0) / (na * nb)
    print(f"  epitope motifs     n={na:5d}  mean rate {st.mean(a):.4f}  median {st.median(a):.4f}")
    print(f"  non-epitope motifs n={nb:5d}  mean rate {st.mean(b):.4f}  median {st.median(b):.4f}")
    print(f"  AUC = {auc:.4f}   (null 0.5)")

print()
print("=" * 74)
print("RETRIEVAL VOLUME")
print("=" * 74)
for c in ("epitope", "non_epitope"):
    g = [r for r in rows if r["cls"] == c]
    h = [r["n_hits"] for r in g]
    print(f"  {c:12s} hits kept after leakage filter: total {sum(h):,d}  "
          f"mean/motif {st.mean(h):.1f}  median {st.median(h):.0f}")
