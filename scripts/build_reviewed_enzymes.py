"""Build reviewed UniProt enzyme cohorts with restricted catalytic/cofactor labels.

Outputs are separate from the older ligand-pocket enzyme dataset. Models are
AlphaFold DB full-length monomers with verified canonical sequence/numbering.
No proximity-derived labels, motif transfer, or invented catalytic triads.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import time
import threading
import urllib.error
import urllib.parse
import urllib.request

from biotite.structure.io.pdb import PDBFile

ROOT = Path(__file__).resolve().parents[1]
QUERIES = {
    "nitrogenases": "reviewed:true AND (ec:1.18.6.1 OR ec:1.19.6.1 OR ec:1.18.6.2 OR ec:1.19.6.2 OR ec:1.18.6.3 OR ec:1.19.6.3)",
    "serine_proteases": "reviewed:true AND keyword:KW-0720",
    "metalloproteases": "reviewed:true AND keyword:KW-0482",
}
NITROGENASE_LITERATURE = [
    "https://pubmed.ncbi.nlm.nih.gov/2500438/",
    "https://pubmed.ncbi.nlm.nih.gov/2287275/",
    "https://doi.org/10.1371/journal.pone.0072751",
]
AA3 = dict(zip("ARNDCQEGHILKMFPSTWYV", [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]))
REQUEST_LOCK = threading.Lock()
LAST_REQUEST = 0.0


def save_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_text(json.dumps(data, indent=2) + "\n")
    tmp.replace(path)


def fetch(url, path, retries=8):
    global LAST_REQUEST
    if path.exists():
        return path.read_bytes(), {}
    for attempt in range(retries):
        try:
            with REQUEST_LOCK:
                time.sleep(max(0., 1.1 - (time.monotonic() - LAST_REQUEST)))
                LAST_REQUEST = time.monotonic()
            req = urllib.request.Request(url, headers={"User-Agent": "IBC-reviewed-enzyme-builder/1.0"})
            with urllib.request.urlopen(req, timeout=90) as response:
                data = response.read(); headers = dict(response.headers)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            return data, {k.lower(): v for k, v in headers.items()}
        except urllib.error.HTTPError as exc:
            if exc.code in (400, 404):
                raise
            if attempt == retries - 1:
                raise
            if exc.code == 429:
                try:
                    delay = float(exc.headers.get("Retry-After", "60"))
                except ValueError:
                    delay = 60.
                print("Rate limited; waiting", delay, "seconds:", url, flush=True)
                time.sleep(max(60., delay))
        except (OSError, TimeoutError):
            if attempt == retries - 1:
                raise
        time.sleep(2 ** attempt)
    raise RuntimeError("retry exhausted")


def position(feature, sequence):
    loc = feature.get("location", {})
    start = loc.get("start", {}); end = loc.get("end", {})
    n = start.get("value")
    if (not isinstance(n, int) or n != end.get("value")
            or start.get("modifier") != "EXACT" or end.get("modifier") != "EXACT"
            or not 1 <= n <= len(sequence)):
        return None
    return n


def select_labels(entry, category):
    """Pure selector. Only exact single-residue UniProt features are eligible."""
    if entry.get("entryType") != "UniProtKB reviewed (Swiss-Prot)":
        return [], "not_reviewed"
    sequence = entry["sequence"]["value"]
    active = [f for f in entry.get("features", []) if f["type"] == "Active site"]
    selected = []
    if category == "serine_proteases":
        # Deliberately reject multi-domain proteins, dyads and incomplete triads.
        if len(active) != 3 or any("charge relay system" not in f.get("description", "").lower() for f in active):
            return [], "not_exactly_one_explicit_charge_relay_triad"
        pos = [position(f, sequence) for f in active]
        if None in pos or len(set(pos)) != 3:
            return [], "non_exact_triad"
        aa = sorted(sequence[p-1] for p in pos)
        if aa not in [sorted("HDS"), sorted("HES")]:
            return [], "noncanonical_triad_composition"
        selected = [(f, "catalytic_triad") for f in active]
    elif category == "nitrogenases":
        for f in entry.get("features", []):
            if f["type"] != "Binding site":
                continue
            ligand = f.get("ligand", {}).get("name", "").lower()
            if "cluster" not in ligand:
                continue
            if any(x in ligand for x in ("7fe-mo", "7fe-v", "8fe-9s")):
                role = "catalytic_cofactor_ligand"
            elif "8fe-7s" in ligand:
                role = "p_cluster_electron_transfer_ligand"
            elif "4fe-4s" in ligand:
                role = "fe_protein_electron_transfer_ligand"
            else:
                continue
            selected.append((f, role))
    elif category == "metalloproteases":
        # Require an actual catalytic residue, not merely a metal-binding protein.
        if not active:
            return [], "no_active_site_feature"
        selected = [(f, "catalytic_residue") for f in active]
        for f in entry.get("features", []):
            ligand = f.get("ligand", {})
            if f["type"] == "Binding site" and ligand.get("note", "").strip().lower() == "catalytic":
                selected.append((f, "catalytic_metal_ligand"))
    else:
        raise ValueError(category)
    labels = []
    for feature, role in selected:
        p = position(feature, sequence)
        if p is None:
            return [], "non_exact_selected_feature"
        if sequence[p-1] not in AA3:
            return [], "nonstandard_label_residue"
        labels.append({"position": p, "residue": sequence[p-1], "label_type": role,
                       "feature_type": feature["type"], "description": feature.get("description", ""),
                       "ligand": feature.get("ligand", {}), "evidence": feature.get("evidences", []),
                       "literature_context": NITROGENASE_LITERATURE if category == "nitrogenases" else []})
    if not labels:
        return [], "no_eligible_features"
    return sorted(labels, key=lambda x: (x["position"], x["label_type"])), "eligible"


def validate_model(text, sequence, labels):
    arr = PDBFile.read(io.StringIO(text)).get_structure(model=1)
    atoms = arr[~arr.hetero]
    chains = set(atoms.chain_id)
    if len(chains) != 1:
        raise ValueError("model must contain one chain")
    residues = {}
    for chain, n, ins, name in zip(atoms.chain_id, atoms.res_id, atoms.ins_code, atoms.res_name):
        if str(ins).strip():
            raise ValueError("insertion code in canonical model")
        key = int(n)
        if key in residues and residues[key] != str(name):
            raise ValueError("ambiguous residue identity")
        residues[key] = str(name)
    if set(residues) != set(range(1, len(sequence)+1)):
        raise ValueError("model numbering/coverage mismatch")
    for n, aa in enumerate(sequence, 1):
        if aa not in AA3 or residues[n] != AA3[aa]:
            raise ValueError("model sequence mismatch")
    for label in labels:
        if residues[label["position"]] != AA3[label["residue"]]:
            raise ValueError("label residue mismatch")
    return str(next(iter(chains)))


def model_for(task):
    entry, labels, out = task
    ac = entry["primaryAccession"]; sequence = entry["sequence"]["value"]
    try:
        raw, _ = fetch("https://alphafold.ebi.ac.uk/api/prediction/" + ac,
                       out / "cache/alphafold" / (ac + ".json"))
        models = json.loads(raw)
        matching = [m for m in models if m.get("uniprotAccession") == ac
                    and m.get("sequence") == sequence and m.get("uniprotSequence", sequence) == sequence
                    and m.get("sequenceStart") == 1 and m.get("sequenceEnd") == len(sequence)
                    and not m.get("isComplex", False) and m.get("pdbUrl")]
        if not matching:
            return None, "no_exact_full_length_afdb_model"
        model = sorted(matching, key=lambda x: x.get("entryId", x.get("modelEntityId", "")))[0]
        path = out / "cache/models" / (ac + ".pdb")
        raw, _ = fetch(model["pdbUrl"], path)
        chain = validate_model(raw.decode(), sequence, labels)
        return {"entry": entry, "labels": labels, "chain": chain, "model": model,
                "cached_structure": str(path), "sha256": hashlib.sha256(raw).hexdigest()}, None
    except (OSError, ValueError, KeyError, TypeError) as exc:
        if isinstance(exc, urllib.error.HTTPError) and exc.code == 429:
            raise RuntimeError("persistent rate limit; resume later rather than bias selection") from exc
        return None, type(exc).__name__ + ": " + str(exc)[:200]


def table(path, header, data):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header); w.writerows(data)


def build(args):
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    if (out / "build_manifest.json").exists():
        raise ValueError("finished dataset exists; choose a new output directory")
    manifest = {"built_utc": datetime.now(timezone.utc).isoformat(), "requested_per_class": args.count,
                "reviewed_only": True, "queries": QUERIES, "classes": {},
                "structure_policy": "AlphaFold DB canonical full-length monomers; exact sequence and residue numbering verified; predicted, no cofactors/assembly",
                "selection": "accession ascending; first eligible available exact models; no sequence clustering",
                "label_policy": "UniProt exact single-residue features only; no proximity, ATP, substrate-pocket, or structural-metal labels",
                "evidence_policy": "reviewed is NOT experimentally verified: preserve evidence codes including inference and absent evidence",
                "nitrogenase_scope": "individual reviewed EC-assigned subunits, NOT 100 complete enzyme complexes; catalytic and electron-transfer cofactors distinguished",
                "literature_context": NITROGENASE_LITERATURE,
                "limitations": ["Not a family-independent benchmark", "No direct literature verification for every homologue", "AF monomers omit catalytic metals, cofactors and inter-subunit geometry"]}
    all_panel = []
    for category, query in QUERIES.items():
        directory = out / category; directory.mkdir(exist_ok=True)
        (directory / "structures").mkdir(exist_ok=True)
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({
            "query": query, "format": "json", "size": 500, "sort": "accession asc"})
        selected = []; reasons = Counter(); page = 0; page_provenance = []
        while url and len(selected) < args.count:
            cache = out / "cache/uniprot" / f"{category}_{page}.json"
            header_path = cache.with_suffix(".headers.json")
            raw, headers = fetch(url, cache)
            if headers:
                save_json(header_path, headers)
            else:
                headers = json.loads(header_path.read_text())
            page_provenance.append({"url": url, "sha256": hashlib.sha256(raw).hexdigest(),
                                    "release": headers.get("x-uniprot-release"),
                                    "total_results": headers.get("x-total-results")})
            entries = json.loads(raw)["results"]
            eligible = []
            for entry in entries:
                labels, reason = select_labels(entry, category)
                if labels:
                    eligible.append((entry, labels, out))
                else:
                    reasons[reason] += 1
            print(category, "page", page, "eligible", len(eligible), flush=True)
            # Bounded batches prevent speculative download of all 500 proteins.
            with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                for start in range(0, len(eligible), args.jobs):
                    batch = eligible[start:start+args.jobs]
                    for result, reason in pool.map(model_for, batch):
                        if result is None:
                            reasons[reason] += 1
                        elif len(selected) < args.count:
                            selected.append(result)
                    print(category, "selected", len(selected), "of", args.count, flush=True)
                    if len(selected) >= args.count:
                        break
            link = headers.get("link", "")
            nxt = re.search(r'<([^>]+)>; rel="next"', link)
            url = nxt.group(1) if nxt else None
            page += 1
        annotation_rows = []; structure_rows = []; details = []; fasta = []
        for item in selected:
            entry = item["entry"]; ac = entry["primaryAccession"]; chain = item["chain"]
            dest = directory / "structures" / (ac + ".pdb")
            dest.write_bytes(Path(item["cached_structure"]).read_bytes())
            name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
            organism = entry["organism"]["scientificName"]
            structure_rows.append([ac, ac, str(dest), category, name, organism,
                                   len(entry["sequence"]["value"]), len({l['position'] for l in item['labels']}),
                                   "AlphaFold DB predicted monomer", item["model"]["pdbUrl"], item["sha256"]])
            all_panel.append([ac, str(dest), category])
            fasta.append(f">{ac} {name}\n{entry['sequence']['value']}\n")
            for label in item["labels"]:
                annotation_rows.append([ac, ac, chain, label["position"], AA3[label["residue"]], 1,
                                        label["label_type"], label["feature_type"], label["ligand"].get("name", ""),
                                        label["description"], json.dumps(label["evidence"], separators=(",", ":"))])
            details.append({"accession": ac, "uniprot_entry": entry, "labels": item["labels"],
                            "alphafold_metadata": item["model"], "structure_sha256": item["sha256"]})
        table(directory / "structures.tsv", ["instance", "uniprot_accession", "structure_file", "category", "protein_name", "organism", "sequence_length", "n_annotated", "structure_source", "structure_url", "sha256"], structure_rows)
        table(directory / "panel.tsv", ["instance", "structure_file", "category"], [[r[0], r[2], category] for r in structure_rows])
        table(directory / "annotated_residues.tsv", ["instance", "uniprot_accession", "chain", "res_id", "res_name", "is_annotated", "label_type", "feature_type", "site_ligand", "description", "evidence_json"], annotation_rows)
        (directory / "sequences.fasta").write_text("".join(fasta))
        save_json(directory / "provenance.json", details)
        manifest["classes"][category] = {"selected": len(selected), "annotation_rows": len(annotation_rows),
                                       "label_types": dict(Counter(r[6] for r in annotation_rows)),
                                       "excluded_or_unavailable": dict(reasons), "pages": page_provenance}
        save_json(out / "build_progress.json", manifest)
    table(out / "panel.tsv", ["instance", "structure_file", "category"], all_panel)
    manifest["complete"] = all(v["selected"] == args.count for v in manifest["classes"].values())
    save_json(out / "build_manifest.json", manifest)
    print(json.dumps(manifest["classes"], indent=2), flush=True)
    if not manifest["complete"]:
        raise RuntimeError("quota shortfall: see build_manifest.json; no padding performed")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=ROOT / "investigation_databases/enzymes_reviewed_catalytic")
    p.add_argument("--count", type=int, default=100)
    p.add_argument("--jobs", type=int, default=6)
    args = p.parse_args()
    if args.count < 1 or args.jobs < 1:
        p.error("count and jobs must be positive")
    build(args)


if __name__ == "__main__":
    main()