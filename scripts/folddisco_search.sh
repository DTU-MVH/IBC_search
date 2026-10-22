#!/usr/bin/env bash
#
# Run FoldDisco over the queries emitted by `ibc-search scan`.
#
# This is the handoff step. `ibc-search scan` extracts candidate surface
# motifs and writes `folddisco_queries.tsv`; this script feeds that file to the
# real FoldDisco binary against a local index and joins the hits back to the
# patch that produced them.
#
#   status   report what is installed, extracted and runnable
#   extract  unpack the downloaded PDB index archive (slow, ~100 GB)
#   search   run FoldDisco over a results directory
#
# Nothing here establishes antigenicity. A hit is a geometric motif match
# against a structure database, and nothing more.

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults. Every one is overridable by flag or environment variable.
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# The downloaded-database directory, one level above the repo.
DB_DIR="${IBC_FOLDDISCO_DB_DIR:-$(cd -- "${REPO_ROOT}/.." && pwd)/folddisco_databases}"
ARCHIVE="${IBC_FOLDDISCO_ARCHIVE:-${DB_DIR}/pdb_folddisco.tar.lz4}"

# The index table itself. FoldDisco also reads sidecars (.lookup/.offset/.type)
# sitting next to it; it reports those itself if they are missing.
INDEX="${IBC_FOLDDISCO_INDEX:-${DB_DIR}/index/pdb_folddisco}"

FOLDDISCO="${IBC_FOLDDISCO_BIN:-folddisco}"

# The build the query contract in docs/FOLDDISCO_CONTRACT.md was verified
# against. A mismatch does not stop a search; it prints a loud warning,
# because every rule on that page would need re-checking.
VERIFIED_COMMIT="2a756d9522d5bab8f41bc73ec72ec2327e31cd97"

# Deliberately 1, not the core count. This build does not serialise its writes
# to stdout: with -t > 1, lines from concurrent threads interleave mid-line and
# those rows are unrecoverable. Measured over repeated -t 10 runs of the same
# 219 queries, 4-8 of ~1225 hit rows were torn each time, varying run to run,
# and a few survived with a plausible field count but a corrupted query string.
# At -t 1 the output was clean every time. Correctness is the default; raise it
# with -t and the script will warn and quarantine what it can detect.
THREADS="${IBC_FOLDDISCO_THREADS:-1}"

# Keep only the top N structures per query, ranked by IDF, BEFORE the
# per-match stage. This is not a cosmetic limit: it is what makes a search
# over the whole PDB terminate. Measured against the 230655-structure index,
# one 3-residue query:
#
#   --top 10     0.18 s        --top 1000    4.19 s
#   --top 100    0.41 s        (uncapped)    > 15 min, no output
#
# The hash lookup itself is fast (0.18 s); the cost is the match stage, which
# loads and superposes every candidate structure. A short motif matches a
# large fraction of the PDB, so uncapped the work is effectively unbounded.
#
# THIS IS A RETRIEVAL CUTOFF AND IT CHANGES WHAT YOU SEE. A true match ranked
# below N is not reported. Raise it, or pass --no-top, when recall matters
# more than runtime.
TOP="${IBC_FOLDDISCO_TOP:-1000}"

# Per-match columns. `query_residues` is the reliable side of the round trip
# (it echoes exactly what was sent); `matching_residues` is NOT reliably
# positionally aligned to it -- see docs/FOLDDISCO_CONTRACT.md.
FORMAT_OUTPUT="tid,idf,rmsd,node_count,query_residues,matching_residues"

die() { printf 'error: %s\n' "$*" >&2; exit 2; }
note() { printf '%s\n' "$*" >&2; }

usage() {
    sed -n '3,16p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    cat >&2 <<'EOF'

usage:
  folddisco_search.sh status
  folddisco_search.sh extract [--force]
  folddisco_search.sh search <results_dir> [options] [-- <extra folddisco args>]
  folddisco_search.sh search-ensemble <ensemble_dir> [-j N] [search options]
  folddisco_search.sh join <results_dir> [-o PATH] [--raw PATH] [--strict]

search options:
  -i, --index PATH       index table to search   [$IBC_FOLDDISCO_INDEX]
  -t, --threads N        FoldDisco threads [1 -- see the note below]
  -o, --output PATH      hits file  [<results_dir>/folddisco_hits.tsv]
      --format-output S  comma-separated per-match columns
      --top N            keep the top N structures per query by IDF [1000]
      --no-top           no cap -- see the warning below
      --rmsd X           maximum RMSD cutoff
      --keep-raw         keep the unjoined FoldDisco output alongside the hits
      --strict           fail if FoldDisco emitted any torn (interleaved) row
      --no-exclude-self  keep hits to the scanned entry itself (excluded by default)
      --exclude-id ID    also drop hits to this PDB entry; repeatable

Anything after `--` is passed through to `folddisco query` untouched, so the
filtering flags this script does not wrap (--sampling-ratio, --covered-node,
--tm-score, ...) remain reachable.

On --top: against the whole-PDB index a single 3-residue query took 0.18 s at
--top 10, 0.41 s at --top 100, 4.2 s at --top 1000, and did not finish in 15
minutes uncapped. The hash lookup is fast; the match stage loads and
superposes every candidate, and a short motif matches much of the PDB. The
cap is therefore on by default. It is a retrieval cutoff: a real match ranked
below N is not reported.

On threads: the installed build does not write stdout line-atomically. At
-t 10 roughly 0.5% of hit rows came back torn, differing run to run; at -t 1
none ever did. The default is therefore 1. Torn rows are detected and
quarantined rather than parsed, but they are still lost, so prefer -t 1 for
anything whose numbers you intend to report.
EOF
}

# ---------------------------------------------------------------------------
# status
# ---------------------------------------------------------------------------

cmd_status() {
    printf 'folddisco binary\n'
    if command -v "${FOLDDISCO}" >/dev/null 2>&1; then
        local resolved installed
        resolved="$(command -v "${FOLDDISCO}")"
        installed="$("${FOLDDISCO}" version 2>/dev/null | tr -d '[:space:]')"
        printf '  path      %s\n  version   %s\n' "${resolved}" "${installed:-<unreported>}"
        if [[ "${installed}" == "${VERIFIED_COMMIT}" ]]; then
            printf '  contract  matches the verified build\n'
        else
            printf '  contract  MISMATCH -- verified against %s; re-check docs/FOLDDISCO_CONTRACT.md\n' \
                "${VERIFIED_COMMIT}"
        fi
    else
        printf '  NOT FOUND on PATH (looked for %s)\n' "${FOLDDISCO}"
    fi

    printf '\nindex\n  expected  %s\n' "${INDEX}"
    if [[ -f "${INDEX}" ]]; then
        printf '  status    extracted (%s)\n' "$(du -h "${INDEX}" | cut -f1)"
    else
        printf '  status    NOT EXTRACTED\n'
        if [[ -f "${ARCHIVE}" ]]; then
            printf '  archive   %s (%s) -- run: %s extract\n' \
                "${ARCHIVE}" "$(du -h "${ARCHIVE}" | cut -f1)" "$0"
        else
            printf '  archive   not found at %s\n' "${ARCHIVE}"
        fi
    fi

    printf '\nstructure files\n'
    if [[ -d "$(dirname -- "${INDEX}")/pdb" ]]; then
        printf '  %s\n' "$(dirname -- "${INDEX}")/pdb"
    else
        printf '  not extracted (the archive carries index/pdb/**/pdb*.ent.gz)\n'
    fi
}

# ---------------------------------------------------------------------------
# extract
# ---------------------------------------------------------------------------

cmd_extract() {
    local force=0
    [[ "${1:-}" == "--force" ]] && force=1

    [[ -f "${ARCHIVE}" ]] || die "archive not found: ${ARCHIVE}"
    command -v lz4 >/dev/null 2>&1 || die "lz4 not found (brew install lz4)"

    if [[ -f "${INDEX}" && ${force} -eq 0 ]]; then
        die "index already present at ${INDEX} (use --force to re-extract)"
    fi

    # Measured, not estimated: the tar stream is 65.87 GiB logical across
    # 230659 files, which is 66.31 GiB once rounded to 4 KiB blocks. lz4 buys
    # almost nothing here (the .ent.gz files are already compressed and the
    # index is high-entropy), so the unpack is barely larger than the archive.
    # 80 GiB of headroom covers it; the archive itself is kept, so a disk
    # holding both needs ~132 GiB.
    local avail_gb
    avail_gb="$(df -g "${DB_DIR}" | awk 'NR==2 {print $4}')"
    printf 'extracting %s\n  into      %s\n  needs     ~67 GiB (230659 files)\n  free disk %s GiB\n' \
        "${ARCHIVE}" "${DB_DIR}" "${avail_gb}" >&2
    if [[ "${avail_gb}" -lt 80 ]]; then
        die "only ${avail_gb} GiB free; the unpack needs ~67 GiB and this script wants 80 GiB of headroom"
    fi
    # Decompressing the whole archive measured 19 s; the cost is writing 66 GiB
    # and creating 230k files, so expect minutes rather than hours.
    note 'Expect a few minutes: ~66 GiB written across 230k files.'

    if command -v pv >/dev/null 2>&1; then
        pv "${ARCHIVE}" | lz4 -dc | tar -xf - -C "${DB_DIR}"
    else
        note 'install pv for a progress bar; running without one'
        lz4 -dc "${ARCHIVE}" | tar -xf - -C "${DB_DIR}"
    fi

    note "done; index should now be at ${INDEX}"
    cmd_status
}


# ---------------------------------------------------------------------------
# Self-hit exclusion
#
# The index is the whole PDB, so scanning a deposited entry and searching it
# retrieves that entry matching itself. Those rows are real but circular, and
# they dominate: on 2N74 the entry was its own top hit for 56 of 87 motifs and
# 29% of all hit rows were self-matches.
#
# The entry ID comes from the scan's own manifest (input.structure_id) and is
# used ONLY when it looks like a PDB entry ID, matching the same grammar the
# package uses. For a structure that is not a deposited entry -- a design, a
# local model -- there is nothing to exclude and nothing is guessed.
#
# NOTE what this does NOT do: it removes the entry itself, not homologues or
# other depositions of the same protein. Those are equally uninformative and
# cannot be identified from an entry ID alone. Pass --exclude-id for each one
# you know of.
# ---------------------------------------------------------------------------

detect_self_id() {
    local results_dir="$1" manifest="$1/run_manifest.json"
    [[ -f "${manifest}" ]] || return 0
    python3 - "${manifest}" <<'PY'
import json, re, sys
try:
    sid = (json.load(open(sys.argv[1])).get("input") or {}).get("structure_id") or ""
except Exception:
    sys.exit(0)
# Same grammar as ibc_search.structure_io.PDB_ID_RE.
print(sid.lower() if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", sid) else "")
PY
}

# ---------------------------------------------------------------------------
# join: attach patch identity to FoldDisco output
#
# Separated from `search` so a raw file can be joined on its own -- notably
# when a long search was interrupted and the completed query blocks are still
# worth keeping.
# ---------------------------------------------------------------------------

do_join() {
    local index_file="$1" raw="$2" output="$3" strict="${4:-0}" exclude_ids="${5:-}"
    local torn_count_file
    torn_count_file="$(mktemp -t folddisco_torn)"

    # The join key is the query string itself, which is safe to rely on:
    # `query_residues` is an exact echo of the tokens that were sent (unlike
    # `matching_residues`, see the contract doc). A handful of patches can
    # share one query string -- the same member set extracted at two radii --
    # and those are genuinely indistinguishable in FoldDisco's output, so they
    # are reported as `id1|id2` rather than silently resolved to one.
    awk -F'\t' -v OFS='\t' \
        -v quarantine="${output%.tsv}.malformed.tsv" \
        -v tornflag="${torn_count_file}" \
        -v idxfile="${index_file}" \
        -v exclude_ids="${exclude_ids}" \
        -v selfhits="${output%.tsv}.selfhits.tsv" '
        # pass 1: query string -> patch identity, from the scan index file
        NR == FNR {
            if (FNR > 1) {
                k = $7
                if (k in patch) {
                    patch[k] = patch[k] "|" $2
                    if (radius[k] != $4) radius[k] = radius[k] "|" $4
                    if (scale[k]  != $5) scale[k]  = scale[k]  "|" $5
                    ambiguous[k] = 1
                } else {
                    patch[k] = $2; radius[k] = $4; scale[k] = $5
                }
                total_queries++
            }
            next
        }
        # pass 2: the FoldDisco output, whose blocks may be in any order
        header == "" {
            header = $0
            n_blocks = 1
            ncol = NF
            for (i = 1; i <= NF; i++) {
                if ($i == "query_residues") qcol = i
                if ($i == "tid") tidcol = i
            }
            n_excl = split(exclude_ids, excl, " ")
            if (n_excl && !tidcol) {
                print "WARNING: --format-output has no tid column, so self-hits " \
                      "cannot be excluded" > "/dev/stderr"
                n_excl = 0
            }
            if (!qcol) {
                print "error: --format-output must include query_residues; " \
                      "it is the key that joins a hit back to its patch" > "/dev/stderr"
                exit 3
            }
            print "patch_id", "radius_angstrom", "scale", $0
            next
        }
        $0 == header { n_blocks++; next } # a later query block repeats it
        # This build does not write stdout line-atomically when -t > 1:
        # lines from concurrent threads interleave mid-line, yielding rows
        # with the wrong field count. Measured on 5 identical -t 10 runs over
        # 1225 hit rows: 2, 0, 0, 2 and 4 torn rows; at -t 1, never any.
        # Such a row is unparseable, so it is quarantined, never emitted.
        NF != ncol {
            n_torn++
            print > quarantine
            next
        }
        # Drop rows belonging to the scanned entry itself. They are kept in a side file
        # rather than discarded silently, so the exclusion can be audited.
        n_excl && tidcol {
            tl = tolower($tidcol)
            for (e = 1; e <= n_excl; e++) {
                if (index(tl, "pdb" excl[e] ".") || tl == excl[e]) {
                    n_self++
                    print > selfhits
                    next
                }
            }
        }
        {
            k = $qcol
            if (k in patch) {
                seen[k] = 1
                if (k in ambiguous) n_amb++
                print patch[k], radius[k], scale[k], $0
            } else {
                n_orphan++
                print "UNKNOWN", "", "", $0
            }
        }
        END {
            hit = 0
            for (k in seen) hit++
            # n_blocks is how many queries FoldDisco actually emitted a header
            # for. On a complete run it equals the query count; on an
            # interrupted one it is smaller, and the difference is queries that
            # never ran -- which must not be reported as "found nothing".
            printf "queries searched: %d of %d, with hits: %d, with none: %d\n", \
                n_blocks, total_queries, hit, n_blocks - hit > "/dev/stderr"
            if (n_blocks < total_queries)
                printf "WARNING: PARTIAL -- %d of %d queries were never searched " \
                       "(output truncated/interrupted), not searched-and-empty\n", \
                       total_queries - n_blocks, total_queries > "/dev/stderr"
            if (n_amb)
                printf "note: %d rows match a query string shared by >1 patch; " \
                       "their patch_id lists every candidate\n", n_amb > "/dev/stderr"
            if (n_orphan)
                printf "WARNING: %d hit rows carry a query string absent from %s\n", \
                       n_orphan, idxfile > "/dev/stderr"
            if (n_self)
                printf "excluded %d self-hit rows (%s); kept in %s\n", \
                       n_self, exclude_ids, selfhits > "/dev/stderr"
            if (n_torn) {
                printf "WARNING: %d rows were torn by non-atomic " \
                       "multithreaded writes and were DROPPED\n", n_torn > "/dev/stderr"
                printf "         quarantined in %s; re-run with -t 1 for a clean stream\n", \
                       quarantine > "/dev/stderr"
            }
            print n_torn + 0 > tornflag
        }
    ' "${index_file}" "${raw}" > "${output}"

    local n_torn
    n_torn="$(cat "${torn_count_file}")"
    rm -f "${torn_count_file}"
    if [[ "${n_torn}" -eq 0 ]]; then
        rm -f "${output%.tsv}.malformed.tsv"
    elif [[ ${strict} -eq 1 ]]; then
        die "${n_torn} torn output rows and --strict is set; re-run with -t 1"
    fi
}

cmd_join() {
    local results_dir="" output="" raw="" strict=0 exclude_self=1 extra_ids=""
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -o|--output) output="$2"; shift 2 ;;
            --raw)       raw="$2"; shift 2 ;;
            --strict)    strict=1; shift ;;
            --no-exclude-self) exclude_self=0; shift ;;
            --exclude-id) extra_ids="${extra_ids} $(printf %s "$2" | tr "A-Z" "a-z")"; shift 2 ;;
            -h|--help)   usage; exit 0 ;;
            -*)          die "unknown option $1" ;;
            *)           [[ -n "${results_dir}" ]] && die "unexpected argument $1"
                         results_dir="$1"; shift ;;
        esac
    done

    [[ -n "${results_dir}" ]] || { usage; die "a results directory is required"; }
    [[ -d "${results_dir}" ]] || die "not a directory: ${results_dir}"

    [[ -n "${output}" ]] || output="${results_dir}/folddisco_hits.tsv"
    [[ -n "${raw}" ]]    || raw="${output%.tsv}.raw.tsv"
    local index_file="${results_dir}/folddisco_queries.index.tsv"

    [[ -f "${raw}" ]] || die "no raw FoldDisco output at ${raw}"
    [[ -f "${index_file}" ]] || die "no folddisco_queries.index.tsv in ${results_dir}"

    note "joining ${raw}"
    local exclude_ids=""
    if [[ ${exclude_self} -eq 1 ]]; then
        exclude_ids="$(detect_self_id "${results_dir}")"
        if [[ -n "${exclude_ids}" ]]; then
            note "  self    excluding the scanned entry itself (${exclude_ids})"
        fi
    fi
    exclude_ids="$(printf %s "${exclude_ids}${extra_ids}" | tr -s " ")"
    do_join "${index_file}" "${raw}" "${output}" "${strict}" "${exclude_ids}"
    local n_hits
    n_hits="$(( $(grep -c '' "${output}") - 1 ))"
    note "wrote ${n_hits} hit rows to ${output}"
}

# ---------------------------------------------------------------------------
# search-ensemble: one search per member, run concurrently
#
# Parallelism is by PROCESS, not by thread, and that is the whole point. The
# installed build does not write stdout line-atomically, so -t > 1 tears rows
# (see the note on threads). Separate processes each own their own output
# file and cannot interleave, so this gets the throughput of many cores while
# every member still runs at the -t 1 that was measured to be reproducible.
#
# It only parallelises when there is actually an ensemble to parallelise.
# With fewer than two members it says so and runs the single search serially,
# rather than paying setup cost and hiding output behind log files for no
# reason.
# ---------------------------------------------------------------------------

cmd_search_ensemble() {
    local ensemble_dir="" jobs=0
    local -a fwd=()
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -j|--jobs)   jobs="$2"; shift 2 ;;
            -h|--help)   usage; exit 0 ;;
            --)          shift; fwd+=("$@"); break ;;
            -*)          fwd+=("$1"); shift ;;
            *)           if [[ -z "${ensemble_dir}" ]]; then ensemble_dir="$1"; else fwd+=("$1"); fi; shift ;;
        esac
    done

    [[ -n "${ensemble_dir}" ]] || { usage; die "an ensemble directory is required"; }
    [[ -d "${ensemble_dir}" ]] || die "not a directory: ${ensemble_dir}"

    # A member is any run directory holding an unsearched query file.
    local -a members=()
    local d
    for d in "${ensemble_dir}"/run_baseline "${ensemble_dir}"/run_[0-9][0-9]; do
        [[ -f "${d}/folddisco_queries.tsv" ]] && members+=("${d}")
    done

    if [[ ${#members[@]} -eq 0 ]]; then
        die "no run_* directories with a folddisco_queries.tsv in ${ensemble_dir}
  Scan the members first (ibc-search scan ... -o ${ensemble_dir}/run_NN)."
    fi

    if [[ ${#members[@]} -lt 2 ]]; then
        note "only ${#members[@]} member found: not an ensemble, running serially"
        cmd_search "${members[0]}" ${fwd[@]+"${fwd[@]}"}
        return $?
    fi

    # Each process is ~0.7 GB RSS and memory-maps the shared index, so the
    # limit is cores rather than RAM; half the cores leaves the machine usable.
    if [[ "${jobs}" -le 0 ]]; then
        local cores
        cores="$( (sysctl -n hw.physicalcpu 2>/dev/null || nproc 2>/dev/null || echo 4) )"
        jobs=$(( cores / 2 )); [[ ${jobs} -lt 1 ]] && jobs=1
        [[ ${jobs} -gt 8 ]] && jobs=8
    fi

    note "ensemble: ${#members[@]} members, ${jobs} concurrent searches (each -t 1)"

    local started=0 i=0
    for d in "${members[@]}"; do
        ( cmd_search "${d}" ${fwd[@]+"${fwd[@]}"} ) > "${d}/search.log" 2>&1 &
        started=$(( started + 1 ))
        i=$(( i + 1 ))
        if [[ $(( i % jobs )) -eq 0 ]]; then wait; fi
    done
    wait

    local failed=0
    note ""
    note "member                         queries  hit rows  self-excl"
    for d in "${members[@]}"; do
        local q h x
        q="$(grep -o "queries searched: [0-9]* of [0-9]*" "${d}/search.log" | tail -1 | awk "{print \$3\"/\"\$5}" || true)"
        h="$(grep -o "wrote [0-9]* hit rows" "${d}/search.log" | tail -1 | awk "{print \$2}" || true)"
        x="$(grep -o "excluded [0-9]* self-hit" "${d}/search.log" | tail -1 | awk "{print \$2}" || true)"
        if [[ -z "${h}" ]]; then
            failed=$(( failed + 1 ))
            note "$(printf "%-30s %-8s %s" "$(basename "${d}")" "${q:-?}" "FAILED - see ${d}/search.log")"
        else
            note "$(printf "%-30s %-8s %-9s %s" "$(basename "${d}")" "${q:-?}" "${h}" "${x:-0}")"
        fi
        grep -h "WARNING" "${d}/search.log" | sed "s/^/    /" >&2 || true
    done

    if [[ ${failed} -gt 0 ]]; then
        die "${failed} of ${#members[@]} member searches failed"
    fi
    note ""
    note "all ${started} members searched"
}

# ---------------------------------------------------------------------------
# search
# ---------------------------------------------------------------------------

cmd_search() {
    local results_dir="" output="" keep_raw=0 strict=0
    local exclude_self=1 extra_ids=""
    local -a passthrough=()

    while [[ $# -gt 0 ]]; do
        case "$1" in
            -i|--index)         INDEX="$2"; shift 2 ;;
            -t|--threads)       THREADS="$2"; shift 2 ;;
            -o|--output)        output="$2"; shift 2 ;;
            --format-output)    FORMAT_OUTPUT="$2"; shift 2 ;;
            --top)              TOP="$2"; shift 2 ;;
            --no-top)           TOP=""; shift ;;
            --rmsd)             passthrough+=(--rmsd "$2"); shift 2 ;;
            --keep-raw)         keep_raw=1; shift ;;
            --no-exclude-self)  exclude_self=0; shift ;;
            --exclude-id)       extra_ids="${extra_ids} $(printf %s "$2" | tr "A-Z" "a-z")"; shift 2 ;;
            --strict)           strict=1; shift ;;
            -h|--help)          usage; exit 0 ;;
            --)                 shift; passthrough+=("$@"); break ;;
            -*)                 die "unknown option $1 (pass FoldDisco flags after --)" ;;
            *)
                [[ -n "${results_dir}" ]] && die "unexpected argument $1"
                results_dir="$1"; shift ;;
        esac
    done

    [[ -n "${results_dir}" ]] || { usage; die "a results directory is required"; }
    [[ -d "${results_dir}" ]] || die "not a directory: ${results_dir}"

    local query_file="${results_dir}/folddisco_queries.tsv"
    local index_file="${results_dir}/folddisco_queries.index.tsv"
    [[ -f "${query_file}" ]] || die "no folddisco_queries.tsv in ${results_dir} -- run \`ibc-search scan\` first"
    [[ -s "${query_file}" ]] || die "${query_file} is empty: that scan emitted no expressible queries"
    [[ -f "${index_file}" ]] || die "no folddisco_queries.index.tsv in ${results_dir}; cannot join hits to patches"

    command -v "${FOLDDISCO}" >/dev/null 2>&1 || die "folddisco not found on PATH (looked for ${FOLDDISCO})"
    if [[ ! -f "${INDEX}" ]]; then
        die "index not found: ${INDEX}
  The downloaded database is a single .tar.lz4 archive until it is unpacked.
  Run:  $0 extract        (or point --index at another index)"
    fi

    local installed
    installed="$("${FOLDDISCO}" version 2>/dev/null | tr -d '[:space:]')"
    if [[ -n "${installed}" && "${installed}" != "${VERIFIED_COMMIT}" ]]; then
        note "WARNING: installed FoldDisco is ${installed}, but the query contract was"
        note "         verified against ${VERIFIED_COMMIT}. Hits are still produced, but"
        note "         docs/FOLDDISCO_CONTRACT.md may no longer describe this build."
    fi

    [[ -n "${output}" ]] || output="${results_dir}/folddisco_hits.tsv"
    local raw="${output%.tsv}.raw.tsv"
    local n_queries
    n_queries="$(grep -c '' "${query_file}")"

    note "searching ${n_queries} queries from ${query_file}"
    note "  index   ${INDEX}"
    note "  threads ${THREADS}"
    if [[ "${THREADS}" -gt 1 ]]; then
        note "WARNING: -t ${THREADS} > 1. This build interleaves concurrent writes, so"
        note "         some output rows will be torn and dropped (~0.5% observed)."
        note "         Use -t 1 for a result you can reproduce exactly."
    fi

    # --header is load-bearing, not cosmetic: it names the columns, and the
    # join below locates `query_residues` by name so --format-output stays
    # free. Note what the join does NOT do: it never relies on output order.
    # FoldDisco prints one header block per query, but with -t > 1 those
    # blocks arrive in COMPLETION order, not file order (verified: at -t 10
    # the second block emitted was query line 36). Counting header blocks to
    # recover a query's line number is therefore only valid single-threaded,
    # and silently wrong otherwise.
    local -a top_arg=()
    if [[ -n "${TOP}" ]]; then
        top_arg=(--top "${TOP}")
        note "  top     ${TOP} structures per query (by IDF, before matching)"
    else
        note "WARNING: --no-top. Against a whole-PDB index an uncapped per-match"
        note "         search did not finish in 15 minutes for a single 3-residue"
        note "         query. Expect this to run for a very long time."
    fi

    "${FOLDDISCO}" query \
        -i "${INDEX}" \
        -q "${query_file}" \
        -t "${THREADS}" \
        --header \
        --per-match \
        ${top_arg[@]+"${top_arg[@]}"} \
        --format-output "${FORMAT_OUTPUT}" \
        ${passthrough[@]+"${passthrough[@]}"} \
        > "${raw}"

    local exclude_ids=""
    if [[ ${exclude_self} -eq 1 ]]; then
        exclude_ids="$(detect_self_id "${results_dir}")"
        if [[ -n "${exclude_ids}" ]]; then
            note "  self    excluding the scanned entry itself (${exclude_ids})"
        fi
    fi
    exclude_ids="$(printf %s "${exclude_ids}${extra_ids}" | tr -s " ")"
    do_join "${index_file}" "${raw}" "${output}" "${strict}" "${exclude_ids}"

    local n_hits
    n_hits="$(( $(grep -c '' "${output}") - 1 ))"
    note "wrote ${n_hits} hit rows to ${output}"

    if [[ ${keep_raw} -eq 1 ]]; then
        note "raw FoldDisco output kept at ${raw}"
    else
        rm -f "${raw}"
    fi

    note 'REMINDER: a hit is a geometric motif match. It is not an epitope and'
    note '          establishes nothing about antibody binding or immunogenicity.'
}

# ---------------------------------------------------------------------------

case "${1:-}" in
    status)        shift; cmd_status "$@" ;;
    extract)       shift; cmd_extract "$@" ;;
    search)        shift; cmd_search "$@" ;;
    search-ensemble) shift; cmd_search_ensemble "$@" ;;
    join)          shift; cmd_join "$@" ;;
    -h|--help|"")  usage; exit 0 ;;
    *)             usage; die "unknown command '$1'" ;;
esac
