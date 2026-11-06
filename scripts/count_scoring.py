"""Score each residue by how many DISTINCT hits it took part in.

The method, in one sentence
---------------------------
A residue's score is the number of distinct reference sites at which that
query residue explicitly matched: count the hits, deduplicate them, and stop.

What this is an alternative to
------------------------------
``residue_scores.py`` gives a residue the **maximum IDF z-score** of any motif
containing it. That is a rarity statistic, standardised against a background
of structures, and it answers "how unusual is the best motif this residue
belongs to". This module answers a different and much blunter question:

    how many distinct places in the reference database look like a motif this
    residue is part of?

The answer is a count of matched sites, so it is reported as a count. **It is
not z-scored and must not be.** A z-score is a position within a fitted null
distribution; standardising a count against a background would make the output
a rarity statistic again, which is the thing this module exists not to be. The
count can be compared between residues of one structure directly, because the
reference database is the same for all of them.

What "distinct" means
---------------------
A hit is a (reference structure, matched residue set) pair -- one matched
*site*. Two hit rows naming the same ``tid`` and the same matched residues are
the same site found twice, so they count once for a residue.

This matters because motifs overlap heavily. A residue sits in many motifs,
those motifs differ by one or two residues, and they match the same place in
the same reference structure over and over. Counting hit rows would mostly
count that redundancy, and would score a residue by how many motifs happened
to contain it rather than by how much the database actually matched.

Deduplication is per residue, not global: the same site found through two
residues' motifs counts once for each of them, because it is one piece of
evidence about each.

Matched tokens are compared as a **set**, so token order cannot split one site
into two. FoldDisco's ``_`` placeholders for unmatched positions are dropped
before comparison (``parse_hit_residues``), so a site is identified by the
residues that actually matched.

What a hit is attributed to
---------------------------
By default only query positions with a non-placeholder target receive credit.
Rows whose query and target fields have different lengths are rejected: this
FoldDisco build sometimes omits unmatched positions, so guessing or skipping
these rows would produce misleading counts. Whole-motif attribution remains
available only through the explicit ``--attribution motif`` compatibility mode.
Matched attribution requires complete positional residue fields; superposition
columns alone are not decoded by this scorer.

Normalising by search effort
----------------------------
``dumb_search`` places one ball per residue, so a residue in the middle of a
chain falls inside more balls than one near a terminus and is searched more
often. The raw count therefore mixes "the database matches this residue's
neighbourhood often" with "this residue was searched often".

``hits_per_motif`` divides the count by ``n_motifs``, the number of emitted
queries the residue belonged to, which removes the second part. Both columns
are written and both are reported: the raw count is still what the 3-D
figures colour on, and the ratio is the effort-corrected view.

The ratio is NOT the mean hit count of the residue's motifs. Its numerator is
deduplicated across those motifs, so motifs that overlap and match the same
reference site are counted once, and the ratio lies below that mean by an
amount that grows with how redundant the motifs are. It is also bounded above
by FoldDisco's ``--top`` cap, so on a censored run it compresses toward that
cap just as the raw count does.

Zero is a result; absent is not
-------------------------------
A residue whose motifs were searched and matched nothing scores **0**: it was
assessed and the database had nothing like it. A residue in no emitted query
is **absent** from the output, because it was never assessed, and writing 0
for it would hide the difference. ``residue_scores.py`` makes the same
distinction; here the two cases are far more common, since any motif too large
or too small for the FoldDisco contract emits no query.

What a count does not mean
--------------------------
A hit establishes that two residue constellations are geometrically similar
under a hash-based tolerance. A high count means many reference structures
contain a similar constellation -- that the motif is **common**, which is
close to the opposite of what a high IDF z-score means. It is not evidence of
antibody binding, antigenicity, or anti-drug-antibody reactivity, and it is
not calibrated against anything. It is a count of matches.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ibc_search.folddisco_contract import (  # noqa: E402
    ContractError, parse_hit_residues, parse_query_token,
)

ResidueKey = tuple[str, int]
SiteKey = tuple[str, tuple[str, ...]]


def load_motif_residues(results_dir: Path) -> dict[str, list[ResidueKey]]:
    """patch_id -> its residues, from the scan."""
    out: dict[str, list[ResidueKey]] = defaultdict(list)
    path = results_dir / "patch_members.tsv"
    if not path.is_file():
        return out
    for r in csv.DictReader(path.open(), delimiter="\t"):
        out[r["patch_id"]].append((r["auth_chain_id"], int(r["auth_seq_id"])))
    return out


def load_queried_patches(results_dir: Path) -> set[str]:
    """The patches that were actually searched.

    Read from the query index rather than inferred from the hits, so that a
    motif that was searched and matched nothing is distinguishable from one
    that was never emitted as a query at all.
    """
    path = results_dir / "folddisco_queries.index.tsv"
    if not path.is_file():
        return set()
    return {
        r["patch_id"]
        for r in csv.DictReader(path.open(), delimiter="\t")
        if r.get("patch_id")
    }


def site_key(tid: str, matching_residues: str) -> SiteKey:
    """Identify the matched site a hit row names.

    The matched tokens are sorted, so a site is a *set* of reference residues
    and token order cannot make one site look like two. Unmatched-position
    placeholders are dropped, so the key names what matched and nothing else.
    """
    parsed = parse_hit_residues("", matching_residues)
    return tid, tuple(sorted(parsed.matched_tokens))


def split_patch_ids(patch_id: str) -> list[str]:
    """The patch ids a hit row names, splitting the pipe-joined form.

    ``folddisco_search.sh`` reports ``id1|id2`` when several patches share a
    query string, because FoldDisco's output cannot tell them apart -- they
    ARE the same residue set, so the hit belongs to every one of them. The
    join is deliberate there and this is the matching read.

    Reading the field whole instead silently loses those rows, since the
    joined string matches no single queried patch. MEASURED on the 5.5 A
    background run: 0.43% of hit rows, in most structures.
    """
    return [p for p in patch_id.split("|") if p]


def load_motif_sites(results_dir: Path) -> dict[str, set[SiteKey]]:
    """patch_id -> the distinct sites its hits matched."""
    out: dict[str, set[SiteKey]] = defaultdict(set)
    path = results_dir / "folddisco_hits.tsv"
    if not path.is_file():
        return out
    for r in csv.DictReader(path.open(), delimiter="\t"):
        patch_id = r.get("patch_id")
        if not patch_id:
            continue
        key = site_key(r.get("tid", ""), r.get("matching_residues", ""))
        if not key[1]:
            # Every position was a placeholder: nothing actually matched, so
            # there is no site to count.
            continue
        for one in split_patch_ids(patch_id):
            out[one].add(key)
    return out


def residue_hit_counts(results_dir: Path, *, attribution: str = "matched",
                       report: dict[str, int] | None = None,
                       ) -> dict[ResidueKey, int]:
    """Distinct matched sites per residue.

    Keys are only the residues that were assessed -- those in at least one
    emitted query. By default only explicitly matched query positions receive
    credit. ``motif`` retains the historical whole-motif attribution
    explicitly.

    Ambiguous rows are SKIPPED, per row, and counted
    -----------------------------------------------
    Where FoldDisco returns fewer target tokens than query tokens there is no
    positional correspondence to recover, so such a row cannot say which
    residue matched and is dropped from the ``matched`` tally.

    It is dropped alone. An earlier version raised instead, on the reasoning
    that a partial table is worse than none -- but MEASURED on the 187
    structure background run that discarded 43 whole structures (23%) over a
    few hundred rows each, while 144 structures had not one ambiguous row.
    The loss was also not random: it fell on structures whose motifs match
    partially, which is plausibly correlated with whatever is being measured.
    Dropping 0.5-6% of one structure's rows distorts its counts far less than
    deleting the structure distorts the set.

    The skipping is therefore **counted, not silent**: pass ``report`` and it
    receives ``rows_total``, ``rows_ambiguous_skipped`` and
    ``rows_no_match``, which ``main`` writes into the summary so the cost is
    visible wherever the numbers are read.
    """
    members = load_motif_residues(results_dir)
    queried = load_queried_patches(results_dir)
    if attribution not in {"matched", "motif"}:
        raise ValueError(f"unknown attribution: {attribution!r}")
    tally = {"rows_total": 0, "rows_ambiguous_skipped": 0, "rows_no_match": 0}
    per_residue: dict[ResidueKey, set[SiteKey]] = {}
    for patch_id in queried:
        for key in members.get(patch_id, ()):
            per_residue.setdefault(key, set())
    if attribution == "motif":
        sites = load_motif_sites(results_dir)
        for patch_id in queried:
            for key in members.get(patch_id, ()):
                per_residue[key].update(sites.get(patch_id, ()))
        if report is not None:
            report.update(tally)
        return {key: len(found) for key, found in per_residue.items()}

    path = results_dir / "folddisco_hits.tsv"
    if not path.is_file():
        raise ValueError(f"missing search output: {path}; cannot infer zero hits")
    with path.open() as fh:
        for line, row in enumerate(csv.DictReader(fh, delimiter="\t"), start=2):
            patch_id = row.get("patch_id", "")
            # A row may name several patches that share one query string.
            # They are the same residue set, so the hit is evidence about
            # each of them and every one must be queried and consistent.
            row_patches = split_patch_ids(patch_id)
            if not row_patches:
                raise ValueError(f"{path}:{line}: row names no patch")
            unknown = [p for p in row_patches if p not in queried]
            if unknown:
                raise ValueError(
                    f"{path}:{line}: unknown or unqueried patch "
                    f"{(unknown[0] if len(unknown) == 1 else unknown)!r}")
            tally["rows_total"] += 1
            query = row.get("query_residues", "")
            target = row.get("matching_residues", "")
            parsed = parse_hit_residues(query, target)
            prefix = f"{path}:{line}: patch {patch_id}: "
            try:
                keys = [parse_query_token(t) for t in parsed.query_tokens]
            except ContractError as exc:
                raise ValueError(prefix + str(exc)) from exc
            if not keys or len(set(keys)) != len(keys):
                raise ValueError(prefix + "query residues do not identify the motif members uniquely")
            # Patches joined on one row must agree on their member set; if
            # they did not, the query string could not have been shared, and
            # crediting them together would be wrong.
            for one in row_patches:
                if set(keys) != set(members.get(one, ())):
                    raise ValueError(
                        prefix + "query residues do not identify the motif "
                        f"members uniquely (patch {one})")
            if not parsed.positionally_alignable:
                # No positional correspondence exists for this row, so it
                # cannot name a matched residue. Skip it and count it; see
                # the docstring for why this is per row and not per
                # structure.
                tally["rows_ambiguous_skipped"] += 1
                continue
            tokens = [t for t in target.strip().split(",") if t]
            for token in parsed.matched_tokens:
                try:
                    parse_query_token(token)
                except ContractError as exc:
                    raise ValueError(prefix + str(exc)) from exc
            if not parsed.matched_tokens:
                tally["rows_no_match"] += 1
                continue
            if not row.get("tid", "").strip():
                raise ValueError(prefix + "matched site has no target structure ID")
            site = site_key(row["tid"], target)
            for key, token in zip(keys, tokens):
                if token != "_":
                    per_residue[key].add(site)
    if report is not None:
        report.update(tally)
    return {key: len(found) for key, found in per_residue.items()}


def residue_motif_counts(results_dir: Path) -> dict[ResidueKey, int]:
    """How many searched motifs each residue belongs to.

    Reported alongside the hit count as context -- whether a high count came
    from a well-covered residue or a thin one -- and used as the DENOMINATOR
    of ``hits_per_motif``, the normalised score.

    What the normalisation does and does not fix
    -------------------------------------------
    A residue near the middle of a chain sits in more balls than one near a
    terminus, so under ``dumb_search`` the raw count partly measures how
    many times the residue was searched at all. Dividing by that removes the
    part of the signal that is pure search effort.

    It does not turn the result into the mean hit count of the residue's
    motifs. The numerator is deduplicated across those motifs, so motifs
    that overlap and match the same reference site contribute once; the
    ratio therefore sits BELOW that mean, by an amount that grows with how
    redundant the residue's motifs are. Read it as distinct evidence per
    unit of search effort, not as a per-motif average.
    """
    members = load_motif_residues(results_dir)
    queried = load_queried_patches(results_dir)
    out: dict[ResidueKey, int] = defaultdict(int)
    for patch_id in queried:
        for key in members.get(patch_id, ()):
            out[key] += 1
    return dict(out)


def auc(pos: list[float], neg: list[float]) -> float | None:
    """Probability a positive outranks a negative (ties count half).

    Ties are the normal case for counts rather than an edge case, which is
    why the half-credit rule matters more here than it does for z-scores.
    """
    if not pos or not neg:
        return None
    merged = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    ranks: dict[int, float] = {}
    i = 0
    while i < len(merged):
        j = i
        while j < len(merged) and merged[j][0] == merged[i][0]:
            j += 1
        r = (i + j + 1) / 2
        for k in range(i, j):
            ranks[k] = r
        i = j
    rp = sum(ranks[k] for k, (_, g) in enumerate(merged) if g == 1)
    return (rp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path, help="directory holding panel.tsv")
    ap.add_argument("--epitopes", type=Path, default=None)
    ap.add_argument("--attribution", choices=("matched", "motif"), default="matched",
                    help="matched: credit only matched query positions, rejecting ambiguous mappings; "
                         "motif: historical credit to every motif member [matched]")
    ap.add_argument("--threshold", type=int, default=1, metavar="N",
                    help="distinct hits at or above which a residue is "
                         "called; 1 means 'matched anything at all' [1]")
    ap.add_argument("-o", "--out", type=Path, default=None)
    args = ap.parse_args(argv)

    # As in residue_scores.py: no epitope file is a legitimate state. A
    # designed binder predicted as a monomer has no antibody in the file, so
    # no residue can be labelled by contact, and counts are then PREDICTIONS
    # WITH NO GROUND TRUTH.
    epi_path = args.epitopes or next(iter(sorted(
        args.run_dir.glob("epitopes_*.tsv"))), None)
    epitope: dict[tuple[str, str, int], int] = {}
    if epi_path is not None:
        for r in csv.DictReader(epi_path.open(), delimiter="\t"):
            epitope[(r["instance"], r["chain"], int(r["res_id"]))] = \
                int(r["is_epitope"])
    labelled = bool(epitope)

    out_path = args.out or (args.run_dir / "residue_hit_counts.tsv")
    # Score the entire panel before opening output: a structure that fails
    # outright must not leave behind a plausible but partial score table.
    # An AMBIGUOUS ROW is no longer such a failure -- it is skipped and
    # counted per structure, and the totals land in the summary.
    panel = list(csv.DictReader((args.run_dir / "panel.tsv").open(), delimiter="\t"))
    panel_reports: dict[str, dict[str, int]] = {}
    try:
        panel_counts = {}
        for row in panel:
            rep: dict[str, int] = {}
            panel_counts[row["instance"]] = residue_hit_counts(
                args.run_dir / row["instance"],
                attribution=args.attribution, report=rep)
            panel_reports[row["instance"]] = rep
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    skipped_total = sum(r.get("rows_ambiguous_skipped", 0)
                        for r in panel_reports.values())
    rows_total = sum(r.get("rows_total", 0) for r in panel_reports.values())
    structures_affected = sum(1 for r in panel_reports.values()
                              if r.get("rows_ambiguous_skipped", 0))
    if skipped_total:
        print(f"skipped {skipped_total:,} of {rows_total:,} hit rows "
              f"({100 * skipped_total / rows_total:.3f}%) with no positional "
              f"residue mapping, across {structures_affected} structures",
              file=sys.stderr)
    per_structure = []
    all_rates: list[float] = []
    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["instance", "chain", "res_id", "is_epitope",
                    "n_distinct_hits", "n_motifs", "hits_per_motif",
                    "called_at_threshold"])
        for row in panel:
            inst = row["instance"]
            inst_dir = args.run_dir / inst
            counts = panel_counts[inst]
            if not counts:
                continue
            motifs = residue_motif_counts(inst_dir)
            pos: list[float] = []
            neg: list[float] = []
            rates: list[float] = []
            tp = fp = fn = tn = 0
            for (chain, res_id), n in sorted(counts.items()):
                lab = epitope.get((inst, chain, res_id))
                if lab is None and labelled:
                    continue  # residue not in the labelled antigen chains
                called = int(n >= args.threshold)
                n_motifs = motifs.get((chain, res_id), 0)
                # Exposure-normalised count: distinct sites per motif the
                # residue was actually searched in. Blank, not zero, where
                # the denominator is unknown -- see rate_note.
                rate = round(n / n_motifs, 4) if n_motifs else ""
                w.writerow([inst, chain, res_id,
                            "" if lab is None else lab, n,
                            n_motifs, rate, called])
                if rate != "":
                    rates.append(rate)
                    all_rates.append(rate)
                (pos if lab else neg).append(float(n))
                if lab and called: tp += 1
                elif lab: fn += 1
                elif called: fp += 1
                else: tn += 1
            if not labelled:
                per_structure.append({
                    "instance": inst, "n_scored": len(neg),
                    "n_called": fp, "call_rate": round(fp / len(neg), 4),
                    "mean_hits": round(st.mean(neg), 4),
                    "median_hits": st.median(neg),
                    "max_hits": int(max(neg)),
                    "n_zero_hit_residues": sum(1 for v in neg if v == 0),
                    "mean_hits_per_motif": (round(st.mean(rates), 4)
                                            if rates else None),
                    "median_hits_per_motif": (round(st.median(rates), 4)
                                              if rates else None),
                    "max_hits_per_motif": (round(max(rates), 4)
                                           if rates else None),
                })
                continue
            a = auc(pos, neg)
            if a is not None:
                per_structure.append({
                    "instance": inst, "n_scored": len(pos) + len(neg),
                    "n_epitope": len(pos), "auc": a,
                    "epitope_mean_hits": st.mean(pos),
                    "other_mean_hits": st.mean(neg),
                    "mean_hits_per_motif": (round(st.mean(rates), 4)
                                            if rates else None),
                    "tp": tp, "fp": fp, "fn": fn, "tn": tn,
                })

    if not per_structure:
        print("error: nothing counted", file=sys.stderr)
        return 2

    note_no_z = ("Counts of distinct matched sites. NOT z-scored: this is a "
                 "count, not a position in a fitted null distribution. A high "
                 "count means the motif is COMMON in the reference database.")

    if not labelled:
        scored = sum(p["n_scored"] for p in per_structure)
        called = sum(p["n_called"] for p in per_structure)
        hit = [p for p in per_structure if p["n_called"] > 0]
        summary = {
            "metric": "n_distinct_hits",
            "attribution": args.attribution,
            "rows_ambiguous_skipped": skipped_total,
            "rows_scored": rows_total,
            "structures_with_skipped_rows": structures_affected,
            "ambiguous_note": (
                "rows where FoldDisco returned fewer target tokens than "
                "query tokens carry no positional correspondence, so they "
                "name no matched residue and are skipped. Skipped per ROW, "
                "not per structure: aborting the structure would discard "
                "thousands of good rows over a few hundred bad ones."),
            "z_scored": False,
            "metric_note": note_no_z,
            "ground_truth_epitopes": False,
            "note": ("No antibody in these structures: every call is a "
                     "prediction with nothing to score it against. No AUC, "
                     "precision or recall is defined here."),
            "structures": len(per_structure),
            "residues_scored": scored,
            "threshold": args.threshold,
            "residues_called": called,
            "call_rate": round(called / scored, 4) if scored else None,
            "structures_with_a_call": len(hit),
            "structures_with_no_call": len(per_structure) - len(hit),
            "mean_hits": round(st.mean([p["mean_hits"] for p in per_structure]), 4),
            "rate_note": ("hits_per_motif = n_distinct_hits / n_motifs: the distinct matched sites a residue collected per motif it was searched in, so a residue in many motifs is not scored higher for that reason alone. The numerator is deduplicated ACROSS a residue's motifs, so this is not the mean hits of its motifs taken separately -- overlapping motifs that match the same site contribute once, which pulls the ratio below that mean. It is bounded above by the search --top cap. Blank where n_motifs is unknown; never zero, since zero is a real rate."),
            "residues_with_a_rate": len(all_rates),
            "pooled_mean_hits_per_motif": (
                round(st.mean(all_rates), 4) if all_rates else None),
            "pooled_median_hits_per_motif": (
                round(st.median(all_rates), 4) if all_rates else None),
            "pooled_sd_hits_per_motif": (
                round(st.pstdev(all_rates), 4) if len(all_rates) > 1 else None),
            "pooled_max_hits_per_motif": (
                round(max(all_rates), 4) if all_rates else None),
            "residues_with_zero_hits": sum(
                p["n_zero_hit_residues"] for p in per_structure),
            "median_calls_per_structure": st.median(
                [p["n_called"] for p in per_structure]),
            "max_calls_in_one_structure": max(
                p["n_called"] for p in per_structure),
        }
        (args.run_dir / "residue_hit_count_summary.json").write_text(json.dumps(
            {"summary": summary, "per_structure": per_structure}, indent=2))
        print(f"wrote {out_path}")
        for k, v in summary.items():
            print(f"  {k:<28} {v}")
        return 0

    aucs = [p["auc"] for p in per_structure]
    tp = sum(p["tp"] for p in per_structure); fp = sum(p["fp"] for p in per_structure)
    fn = sum(p["fn"] for p in per_structure); tn = sum(p["tn"] for p in per_structure)
    base = (tp + fn) / (tp + fp + fn + tn)
    prec = tp / (tp + fp) if tp + fp else float("nan")
    summary = {
        "metric": "n_distinct_hits",
        "attribution": args.attribution,
        "z_scored": False,
        "metric_note": note_no_z,
        "structures": len(per_structure),
        "residues_scored": tp + fp + fn + tn,
        "epitope_residues": tp + fn,
        "epitope_base_rate": round(base, 4),
        "threshold": args.threshold,
        "auc_mean": round(st.mean(aucs), 4),
        "auc_median": round(st.median(aucs), 4),
        "auc_above_0.5": sum(1 for a in aucs if a > 0.5),
        "auc_below_0.5": sum(1 for a in aucs if a < 0.5),
        "precision_at_threshold": round(prec, 4),
        "recall_at_threshold": round(tp / (tp + fn), 4) if tp + fn else None,
        "enrichment_over_base_rate": round(prec / base, 3) if base else None,
        "epitope_mean_hits": round(
            st.mean([p["epitope_mean_hits"] for p in per_structure]), 4),
        "other_mean_hits": round(
            st.mean([p["other_mean_hits"] for p in per_structure]), 4),
        "rate_note": ("hits_per_motif = n_distinct_hits / n_motifs: the distinct matched sites a residue collected per motif it was searched in, so a residue in many motifs is not scored higher for that reason alone. The numerator is deduplicated ACROSS a residue's motifs, so this is not the mean hits of its motifs taken separately -- overlapping motifs that match the same site contribute once, which pulls the ratio below that mean. It is bounded above by the search --top cap. Blank where n_motifs is unknown; never zero, since zero is a real rate."),
        "residues_with_a_rate": len(all_rates),
        "pooled_mean_hits_per_motif": (
            round(st.mean(all_rates), 4) if all_rates else None),
        "pooled_median_hits_per_motif": (
            round(st.median(all_rates), 4) if all_rates else None),
        "pooled_sd_hits_per_motif": (
            round(st.pstdev(all_rates), 4) if len(all_rates) > 1 else None),
    }
    (args.run_dir / "residue_hit_count_summary.json").write_text(
        json.dumps({"summary": summary, "per_structure": per_structure}, indent=2)
    )

    print(f"wrote {out_path}")
    for k, v in summary.items():
        print(f"  {k:<28} {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
