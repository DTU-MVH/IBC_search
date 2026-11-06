"""Command-line interface.

Commands
--------
``scan``
    Parse one structure (local file or PDB entry ID), compute exposure,
    extract discontinuous surface motifs, and write tables, the FoldDisco
    query file and the run manifest.

``show-config``
    Print the fully resolved configuration, so a run's parameters can be
    inspected before committing to it.

``folddisco-info``
    Report what is known about the installed FoldDisco build and the verified
    query contract. Runs no search.

``dumb-search``
    The same pipeline as ``scan``, but motifs come from
    :func:`ibc_search.dumb_search.run_dumb_search`: one ball per residue,
    walked from terminus to terminus, with none of the ball search's
    thresholds. A baseline, written to the same output format.

``list-databases``
    Print the reference databases in the registry, saying for each whether
    its index is present and whether IDF parameters have been fitted against
    it. The two are reported separately because they gate different things.

``folddisco-search``
    Run FoldDisco over a finished scan's queries and join the hits back to
    their patches. This is a thin delegation to
    ``scripts/folddisco_search.sh``: the search itself is still performed by
    the FoldDisco binary, not reimplemented here, so the invocation has one
    definition that is equally usable from a plain shell.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

from . import __version__
from .ball_search_algo import run_ball_search
from .config import (
    default_config_path,
    describe_databases,
    load_config,
    resolve_database,
)
from .dumb_search import run_dumb_search
from .few_balls_search import run_few_balls_search
from .errors import IbcSearchError
from .exposure import compute_exposure
from .logging_utils import setup_logging
from .manifest import build_manifest, folddisco_status
from .outputs import write_manifest, write_tables
from .structure_io import parse_structure, resolve_input


def _parse_overrides(pairs: list[str] | None) -> dict[str, str]:
    overrides: dict[str, str] = {}
    for pair in pairs or []:
        if "=" not in pair:
            raise IbcSearchError(
                f"--set expects 'section.key=value', got {pair!r}"
            )
        key, value = pair.split("=", 1)
        overrides[key.strip()] = value.strip()
    return overrides


def _add_structure_arguments(scan: argparse.ArgumentParser) -> None:
    """Add the options shared by every structure-scanning command.

    Shared rather than duplicated so that `scan` and `dumb-search` cannot
    drift apart on input handling, which is the part of the two commands that
    is genuinely the same.
    """
    scan.add_argument(
        "input",
        help=(
            "path to a .pdb/.ent/.cif/.mmcif file (optionally .gz), or a "
            "PDB entry ID in either form: 4-character (4CHA) or extended "
            "(pdb_00004cha)"
        ),
    )
    scan.add_argument(
        "-o", "--output-dir", type=Path, default=None,
        help="output directory (default: [run].output_dir from the config)",
    )
    scan.add_argument(
        "-c", "--config", type=Path, default=None,
        help="user config file, merged over the shipped defaults",
    )
    scan.add_argument(
        "--chain", action="append", default=None, metavar="ID",
        help="author chain ID to scan; repeatable. Default: all polymer chains",
    )
    scan.add_argument(
        "--model", type=int, default=None,
        help="deposited model number to scan (default: 1)",
    )
    scan.add_argument(
        "--radius", action="append", type=float, default=None, metavar="A",
        help=(
            "ball radius in Angstroms; repeatable. Overrides the configured "
            "radius set"
        ),
    )
    scan.add_argument(
        "--center-mode", choices=("cb", "sidechain_centroid"), default=None,
        help="ball-centre representation (default: cb)",
    )
    scan.add_argument(
        "--workers", type=int, default=None,
        help="worker threads for the ball queries (default 1; output is "
             "identical regardless)",
    )
    scan.add_argument(
        "--set", action="append", default=None, metavar="section.key=value",
        help="override any config key; repeatable",
    )
    scan.add_argument("--log-level", default=None, help="DEBUG/INFO/WARNING/ERROR")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ibc-search",
        description=(
            "Extract candidate discontinuous structural motifs from a protein "
            "surface and prepare them as FoldDisco queries. Patches are "
            "candidate motifs, NOT epitopes, and establish nothing about "
            "antigenicity or immunogenicity."
        ),
    )
    parser.add_argument("--version", action="version", version=f"ibc-search {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    scan = sub.add_parser(
        "scan",
        help="scan one structure from a file path or a PDB entry ID (4-character or extended)",
    )
    _add_structure_arguments(scan)
    scan.add_argument(
        "--seed-rsa", type=float, default=None, metavar="FRAC",
        help=(
            "RSA threshold for seed eligibility. An implementation choice, "
            "not an established epitope criterion"
        ),
    )

    dumb = sub.add_parser(
        "dumb-search",
        help="one ball per residue, terminus to terminus, with no thresholds",
        description=(
            "The dumb_search method. Walk the scanned chains from the "
            "N-terminal residue to the C-terminal residue and place one ball "
            "at every residue on the way; each ball's contents are a motif. "
            "Ball geometry (radii, centre mode, nearest-heavy-atom "
            "membership) is identical to `scan`; what is absent is every "
            "filter -- no RSA threshold for seeds or members, no surface "
            "connectivity check, no deduplication, no size limits. It exists "
            "as a baseline for telling results apart from the effect of "
            "`scan`'s thresholds. Output format is identical, and a motif is "
            "no more an epitope here than there."
        ),
    )
    _add_structure_arguments(dumb)
    dumb.add_argument(
        "--seed-residues", type=Path, default=None, metavar="TSV",
        help="place balls ONLY on the residues this file names, instead of "
             "on every residue. A TSV with chain/res_id (or "
             "auth_chain_id/auth_seq_id); a pdb_id/instance column filters "
             "it to the structure being scanned, and an is_annotated column "
             "is honoured. Membership is unaffected: a ball still collects "
             "every residue within the radius, named or not")
    dumb.add_argument(
        "--seed-rounds", type=int, default=None, metavar="N",
        help="rounds of seeding, with --seed-residues. 1 (default) places "
             "balls on the named residues only. 2 also places a ball on "
             "every residue that falls inside one of those balls, so the "
             "walk covers the annotated SITE and each of its residues is "
             "scored from its own centre. The added residues are exactly "
             "the members round 1 already collected")
    dumb.add_argument(
        "--seed-structure-id", default=None, metavar="ID",
        help="the ID to match the seed table's pdb_id/instance column "
             "against, when the structure file is not named after the "
             "entry (e.g. a `11KP-assembly1.cif.gz` whose annotations are "
             "filed under `11KP`). Defaults to the file stem")

    few = sub.add_parser(
        "few-balls-search",
        help="one ball per annotated residue plus an equal, dispersed "
             "background set on the same protein",
        description=(
            "The few_balls_search method. Place one ball on every annotated "
            "residue, then place the SAME NUMBER of further balls elsewhere "
            "on the same protein, chosen by farthest-point sampling so they "
            "are spread out in three dimensions and stratified by solvent "
            "accessibility so they span buried and exposed positions in the "
            "protein's own proportions. A floor of --min-background balls "
            "applies however few annotations there are. Ball geometry and "
            "membership are identical to `scan` and `dumb-search`; no member "
            "is filtered and nothing is deduplicated. The point is an "
            "internal control: both classes come from one structure, one "
            "radius and one search, so a difference between them cannot be "
            "explained by fold, motif count or search effort. Every motif is "
            "labelled `annotated` or `background` in the seed_class column, "
            "which downstream analysis must read rather than infer."
        ),
    )
    _add_structure_arguments(few)
    few.add_argument(
        "--seed-residues", type=Path, default=None, metavar="TSV",
        help="REQUIRED. The annotated residues, as for `dumb-search`. This "
             "method is defined by an annotated set, so an absent list is an "
             "error rather than a fallback to every residue")
    few.add_argument(
        "--seed-structure-id", default=None, metavar="ID",
        help="the ID to match the seed table's pdb_id/instance column "
             "against, when the structure file is not named after the entry "
             "(e.g. `11KP-assembly1.cif.gz` filed under `11KP`)")
    few.add_argument(
        "--min-background", type=int, default=None, metavar="N",
        help="floor on the background count, however few annotated residues "
             "there are (default: [few_balls_search].min_background_balls). "
             "Where it binds, the two classes are deliberately unequal and "
             "the manifest says so")
    few.add_argument(
        "--background-mode", choices=("protein", "surface"), default=None,
        help="where background balls may sit. `protein` (default) samples "
             "the whole protein, stratified by burial. `surface` restricts "
             "them to residues at or above --background-rsa, which is what "
             "a surface annotation (an interface, an epitope) needs: "
             "otherwise the control differs from the annotated class in "
             "burial as well as in being annotated")
    few.add_argument(
        "--background-rsa", type=float, default=None, metavar="FRAC",
        help="the buried/exposed boundary for stratifying the background "
             "(default: [few_balls_search].background_rsa_threshold)")
    few.add_argument(
        "--no-background-stratify", action="store_true",
        help="take whatever farthest-point sampling returns instead of "
             "matching the protein's buried/exposed composition. NOT "
             "advised: farthest-point sampling favours the hull, so the "
             "background becomes a surface sample and burial is confounded "
             "with annotation")

    show = sub.add_parser("show-config", help="print the resolved configuration")
    show.add_argument("-c", "--config", type=Path, default=None)
    show.add_argument("--set", action="append", default=None, metavar="section.key=value")

    sub.add_parser(
        "list-databases",
        help="print the configured reference databases and their readiness",
        description=(
            "List [folddisco.databases]. A database pairs an index with the "
            "IDF parameters fitted against it, and the two are reported "
            "separately: without the index there is no search at all, while "
            "without fitted parameters the hits can be counted but not "
            "z-scored. Pass a name to `folddisco-search -d` or "
            "`scripts/run_panel.py -d`."
        ),
    )

    sub.add_parser(
        "folddisco-info",
        help="report the installed FoldDisco build and the verified query contract",
    )

    search = sub.add_parser(
        "folddisco-search",
        help=(
            "run FoldDisco over a finished scan's queries (delegates to "
            "scripts/folddisco_search.sh)"
        ),
        description=(
            "Hand a results directory to the real FoldDisco binary and join "
            "the hits back to the patches that produced them. This does not "
            "implement a search in Python: it invokes the documented tool via "
            "scripts/folddisco_search.sh, so there is exactly one "
            "implementation of the invocation. A hit is a geometric motif "
            "match and establishes nothing about antigenicity."
        ),
    )
    search.add_argument(
        "results_dir",
        help="a directory written by `ibc-search scan` (holds folddisco_queries.tsv)",
    )
    search.add_argument("-c", "--config", type=Path, default=None)
    search.add_argument(
        "-d", "--database", default=None,
        help=(
            "named reference database from [folddisco.databases] "
            "(default: [folddisco].default_database). Selecting by name keeps "
            "the index and its fitted IDF z-score parameters together"
        ),
    )
    search.add_argument(
        "-i", "--index", default=None,
        help="explicit index path, overriding --database",
    )
    search.add_argument(
        "-t", "--threads", type=int, default=None,
        help=(
            "FoldDisco threads (default: [folddisco].search_threads = 1). "
            "Above 1 this build tears output rows; see the script's --help"
        ),
    )
    search.add_argument(
        "-o", "--output", default=None,
        help="hits file (default: <results_dir>/folddisco_hits.tsv)",
    )
    search.add_argument(
        "--top", type=int, default=None,
        help=(
            "keep the top N structures per query by IDF, before matching "
            "(default: [folddisco].search_top). A retrieval cutoff, and what "
            "makes a whole-PDB search terminate"
        ),
    )
    search.add_argument(
        "--no-top", action="store_true",
        help="remove the cap; expect a whole-PDB search to run for hours",
    )
    search.add_argument(
        "--num-residue", type=int, default=None,
        help=(
            "skip candidate structures over N residues at the match stage "
            "(default: [folddisco].search_num_residue). The count is summed "
            "over models, so this is what excludes the multi-model NMR "
            "ensembles that dominate the runtime. It changes which "
            "candidates come back"
        ),
    )
    search.add_argument(
        "--no-num-residue", action="store_true",
        help="no candidate size limit, leaving FoldDisco's own 50000",
    )
    search.add_argument(
        "--strict", action="store_true",
        help="fail if FoldDisco emitted any torn row",
    )
    search.add_argument(
        "--keep-raw", action="store_true",
        help="keep the unjoined FoldDisco output",
    )
    search.add_argument(
        "extra", nargs="*", default=[],
        help="after `--`, arguments passed straight to `folddisco query`",
    )
    return parser


def _apply_cli_to_config(args: argparse.Namespace, config: dict[str, Any]) -> None:
    """Fold explicit CLI flags into the resolved configuration.

    Applied after the config file so a flag always wins, and written into the
    same resolved mapping that the manifest records, so there is exactly one
    description of what the run used.
    """
    if getattr(args, "chain", None):
        config["input"]["chains"] = list(args.chain)
    if getattr(args, "model", None) is not None:
        config["input"]["model"] = int(args.model)
    if getattr(args, "radius", None):
        config["ball_search"]["radii_angstrom"] = [float(r) for r in args.radius]
    if getattr(args, "center_mode", None):
        config["ball_search"]["center_mode"] = args.center_mode
    if getattr(args, "seed_residues", None) is not None:
        config["ball_search"]["seed_residues_path"] = str(args.seed_residues)
    if getattr(args, "seed_rounds", None) is not None:
        config["ball_search"]["seed_rounds"] = int(args.seed_rounds)
    if getattr(args, "seed_structure_id", None):
        config["ball_search"]["seed_structure_id"] = str(args.seed_structure_id)
    if getattr(args, "min_background", None) is not None:
        config["few_balls_search"]["min_background_balls"] = int(
            args.min_background)
    if getattr(args, "background_mode", None):
        config["few_balls_search"]["background_mode"] = str(
            args.background_mode)
    if getattr(args, "background_rsa", None) is not None:
        config["few_balls_search"]["background_rsa_threshold"] = float(
            args.background_rsa)
    if getattr(args, "no_background_stratify", False):
        config["few_balls_search"]["background_stratify_by_burial"] = False
    if getattr(args, "seed_rsa", None) is not None:
        config["ball_search"]["seed_rsa_threshold"] = float(args.seed_rsa)
    if getattr(args, "workers", None) is not None:
        config["performance"]["workers"] = int(args.workers)
    if getattr(args, "log_level", None):
        config["run"]["log_level"] = args.log_level
    if getattr(args, "output_dir", None) is not None:
        config["run"]["output_dir"] = str(args.output_dir)


def cmd_scan(
    args: argparse.Namespace,
    argv: list[str],
    method: Callable[..., Any] = run_ball_search,
) -> int:
    """Run one structure end to end with ``method`` as the extraction step.

    ``scan`` and ``dumb-search`` differ in that one argument and nothing else:
    input resolution, exposure, output tables and the manifest are shared, so
    the two methods are compared on their extraction and not on incidental
    differences in how their runs were set up.
    """
    config = load_config(args.config, _parse_overrides(args.set))
    _apply_cli_to_config(args, config)
    from .config import validate_config

    validate_config(config)

    logger = setup_logging(config["run"]["log_level"])
    timings: dict[str, float] = {}

    t0 = time.perf_counter()
    source = resolve_input(args.input, config)
    timings["resolve_input_s"] = time.perf_counter() - t0
    logger.info(
        "input: %s (%s, %s)%s",
        source.path,
        source.file_format,
        f"sha256={source.sha256[:12]}",
        " [cached]" if source.from_cache else "",
    )

    t0 = time.perf_counter()
    parsed = parse_structure(source, config)
    timings["parse_structure_s"] = time.perf_counter() - t0
    logger.info(
        "parsed %d residues / %d heavy atoms; chains scanned: %s",
        len(parsed.residues),
        len(parsed.coords),
        ",".join(parsed.scanned_chains),
    )
    for warning in parsed.warnings:
        logger.warning("%s", warning)

    t0 = time.perf_counter()
    exposure = compute_exposure(parsed, config)
    timings["compute_exposure_s"] = time.perf_counter() - t0
    for warning in exposure.warnings:
        logger.warning("%s", warning)

    t0 = time.perf_counter()
    result = method(parsed, exposure, config)
    timings["ball_search_total_s"] = time.perf_counter() - t0
    for warning in result.warnings:
        logger.warning("%s", warning)
    logger.info(
        "extracted %d patches (%d discontinuous, %d with a FoldDisco query) "
        "from %d seeds",
        result.stats.get("n_patches_total", 0),
        result.stats.get("n_patches_discontinuous", 0),
        result.stats.get("n_patches_query_emitted", 0),
        result.stats.get("n_seeds", 0),
    )

    out_dir = Path(config["run"]["output_dir"])
    t0 = time.perf_counter()
    # write_tables can append warnings of its own -- notably that the
    # FoldDisco query format cannot carry a model number -- so anything added
    # here is logged too, rather than only reaching the manifest.
    warnings_before = len(result.warnings)
    paths = write_tables(out_dir, parsed, exposure, result, config)
    timings["write_tables_s"] = time.perf_counter() - t0
    for warning in result.warnings[warnings_before:]:
        logger.warning("%s", warning)

    manifest = build_manifest(
        command=argv,
        config=config,
        config_source=str(args.config or default_config_path()),
        parsed=parsed,
        exposure=exposure,
        result=result,
        timings=timings,
        output_paths=paths,
    )
    paths["manifest"] = write_manifest(out_dir, manifest)
    logger.info("wrote %d output files to %s", len(paths), out_dir)
    for name, path in paths.items():
        logger.info("  %-22s %s", name, path)
    logger.info(
        "REMINDER: extracted patches are candidate structural motifs, not "
        "epitopes; nothing here establishes antigenicity or immunogenicity"
    )
    return 0


def cmd_show_config(args: argparse.Namespace) -> int:
    config = load_config(args.config, _parse_overrides(args.set))
    print(json.dumps(config, indent=2))
    return 0


def cmd_list_databases() -> int:
    """Print the registry, marking the default and anything not ready."""
    rows = describe_databases(load_config())
    if not rows:
        print("no [folddisco.databases] are configured", file=sys.stderr)
        return 2
    for row in rows:
        mark = "*" if row["is_default"] else " "
        state = "ready" if row["index_present"] else "INDEX MISSING"
        print(f" {mark} {row['name']:<18}{state}")
        print(f"   {'':18}index   {row['index_path']}")
        if row["idf_params_present"]:
            fitted = f"fitted at --top {row['fitted_top']}"
            # The radius is only recorded for databases whose fit noted it;
            # inventing a default here would misreport the fit.
            if row["fitted_radius_angstrom"] is not None:
                fitted += f", radius {row['fitted_radius_angstrom']} A"
            print(f"   {'':18}idf     {row['idf_params_path'].name} ({fitted})")
        elif row["idf_params_path"] is None:
            print(f"   {'':18}idf     none fitted -- hits can be COUNTED but "
                  "not z-scored")
        else:
            print(f"   {'':18}idf     configured but MISSING at "
                  f"{row['idf_params_path']}")
    print("\n * = [folddisco].default_database")
    print("Select one with `folddisco-search -d NAME` or "
          "`scripts/run_panel.py -d NAME`.")
    return 0


def cmd_folddisco_info() -> int:
    from .folddisco_contract import (
        FOLDDISCO_VERIFIED_COMMIT,
        MIN_QUERY_RESIDUES,
        UNMATCHED_TOKEN,
    )

    config = load_config()
    status = folddisco_status(config["folddisco"]["executable"])
    info = {
        "folddisco": status,
        "verified_query_contract": {
            "token_grammar": "<single-character chain><non-negative integer>",
            "example_query": "B57,B102,C195",
            "range_syntax": "B57-60 expands to B57,B58,B59,B60",
            "query_file_format": "<structure_path>\\t<comma-separated tokens>",
            "min_query_residues": MIN_QUERY_RESIDUES,
            "unmatched_position_token": UNMATCHED_TOKEN,
            "cannot_express": [
                "insertion codes (B57A rejected)",
                "multi-character chain IDs (AAA57 rejected)",
                "negative residue numbers (B-5 rejected; '-' is the range separator)",
            ],
            "verified_against_commit": FOLDDISCO_VERIFIED_COMMIT,
        },
        "note": "This component runs no FoldDisco search.",
    }
    print(json.dumps(info, indent=2))
    return 0


def _search_script() -> Path:
    """Locate scripts/folddisco_search.sh.

    The script lives beside the package rather than inside it, because it is
    also meant to be run on its own from a shell. ``$IBC_FOLDDISCO_SCRIPT``
    overrides the lookup for an installed (non-source-tree) layout.
    """
    override = os.environ.get("IBC_FOLDDISCO_SCRIPT")
    if override:
        path = Path(override)
        if not path.is_file():
            raise IbcSearchError(
                f"IBC_FOLDDISCO_SCRIPT points at {path}, which is not a file"
            )
        return path
    path = Path(__file__).resolve().parents[2] / "scripts" / "folddisco_search.sh"
    if not path.is_file():
        raise IbcSearchError(
            f"could not find the search script at {path}. Set "
            "IBC_FOLDDISCO_SCRIPT to its location, or run it directly."
        )
    return path


def cmd_folddisco_search(args: argparse.Namespace) -> int:
    """Delegate to the shell script. No search is implemented in Python.

    Config supplies the defaults; an explicit flag always wins, matching how
    ``scan`` resolves its own options.
    """
    config = load_config(args.config)
    fd = config["folddisco"]

    # An explicit --index wins; otherwise resolve the named database, falling
    # back to the legacy single-index key so older configs keep working.
    index = args.index or ""
    if not index:
        try:
            entry = resolve_database(config, args.database)
            index = str(entry["index_path"])
            print(f"database: {entry['name']}", file=sys.stderr)
        except IbcSearchError:
            if args.database:
                raise
            index = fd.get("reference_index_path") or ""
    threads = args.threads if args.threads is not None else fd.get("search_threads", 1)

    cmd: list[str] = [str(_search_script()), "search", args.results_dir]
    if index:
        cmd += ["-i", str(index)]
    cmd += ["-t", str(threads)]
    if args.output:
        cmd += ["-o", args.output]
    if fd.get("search_format_output"):
        cmd += ["--format-output", str(fd["search_format_output"])]
    if args.no_top:
        cmd.append("--no-top")
    else:
        top = args.top if args.top is not None else fd.get("search_top")
        if top:
            cmd += ["--top", str(top)]
    if args.no_num_residue:
        cmd.append("--no-num-residue")
    else:
        numres = (args.num_residue if args.num_residue is not None
                  else fd.get("search_num_residue"))
        if numres:
            cmd += ["--num-residue", str(numres)]
    if args.strict:
        cmd.append("--strict")
    if args.keep_raw:
        cmd.append("--keep-raw")

    # argparse.REMAINDER keeps the literal `--`; the script expects it.
    extra = list(args.extra or [])
    if extra and extra[0] != "--":
        extra = ["--", *extra]
    cmd += extra

    completed = subprocess.run(cmd, check=False)  # noqa: S603 - argv, no shell
    return completed.returncode


def main(argv: list[str] | None = None) -> int:
    raw = list(sys.argv[1:] if argv is None else argv)
    full_command = ["ibc-search", *raw]
    parser = build_parser()
    args = parser.parse_args(raw)
    try:
        if args.command == "scan":
            return cmd_scan(args, full_command)
        if args.command == "dumb-search":
            return cmd_scan(args, full_command, method=run_dumb_search)
        if args.command == "few-balls-search":
            return cmd_scan(args, full_command, method=run_few_balls_search)
        if args.command == "show-config":
            return cmd_show_config(args)
        if args.command == "list-databases":
            return cmd_list_databases()
        if args.command == "folddisco-info":
            return cmd_folddisco_info()
        if args.command == "folddisco-search":
            return cmd_folddisco_search(args)
    except IbcSearchError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    parser.error(f"unknown command {args.command!r}")
    return 2


if __name__ == "__main__":
    sys.exit(main())
