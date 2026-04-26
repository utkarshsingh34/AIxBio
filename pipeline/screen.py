"""
screen.py

Screen a corpus FASTA against a hazardous BLAST DB. Produces TWO output CSVs:

  1. {basename}_fragments.csv  (NEW)
       - One row per fragment. The "raw measurement": for each fragment, what
         was BLAST's best hit (min e-value, max bitscore, n_hits)?
       - This is the basis for the per-fragment detection-rate analysis
         (Figure 1: BLAST sensitivity profile across fragment lengths).

  2. {basename}_orders.csv
       - One row per order. "Aggregation applied": min e-value across the
         order's fragments, count of fragments with any hit, etc.
       - This is the basis for the per-order ROC analysis (Figure 2:
         deployment scenario with any-flag aggregation).

Conceptually: fragments.csv is the fundamental sensitivity measurement;
orders.csv is what happens when you bundle fragments into realistic orders
and use any-flag aggregation as the screener decision rule.

USAGE (from your pipeline folder, with blast-env active):
    python screen.py --corpus ../data/corpora/T50/honest.fasta --db hazardous_db
    python screen.py --corpus ../data/corpora/T50/evasion_pure.fasta --db hazardous_db

Outputs (next to the input corpus):
    {basename}_hits.tsv        - raw BLAST tabular output
    {basename}_fragments.csv   - per-fragment best-hit summary (NEW)
    {basename}_orders.csv      - per-order any-flag aggregation
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# Fragment ID format from build_corpus.py:
#   order{ID}__frag{IDX}__role{honest|hazard|filler}__parent{NAME}__{START}_{END}
ID_RE = re.compile(
    r"^order(?P<order>\d+)__frag(?P<frag>\d+)__role(?P<role>\w+)"
    r"__parent(?P<parent>.+)__(?P<start>\d+)_(?P<end>\d+)$"
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--corpus", required=True, help="Input corpus FASTA")
    p.add_argument("--db", required=True, help="BLAST DB prefix (e.g. hazardous_db)")
    p.add_argument("--task", default="blastn-short",
                   choices=["blastn", "blastn-short", "megablast", "dc-megablast"],
                   help="blastn task preset (default: %(default)s)")
    p.add_argument("--evalue", type=float, default=10.0,
                   help="BLAST -evalue REPORTING threshold (permissive on purpose; "
                        "evaluate.py decides the screening cutoff). Default: %(default)s")
    p.add_argument("--word-size", type=int, default=None,
                   help="Override word_size (default: task default)")
    p.add_argument("--max-target-seqs", type=int, default=5,
                   help="Max hits per query (default: %(default)s)")
    p.add_argument("--threads", type=int, default=4,
                   help="BLAST threads (default: %(default)s)")
    return p.parse_args()


def parse_fragment_id(qseqid):
    m = ID_RE.match(qseqid)
    if not m:
        return None
    d = m.groupdict()
    return {
        "order_id": int(d["order"]),
        "frag_idx": int(d["frag"]),
        "role": d["role"],
        "parent": d["parent"],
        "frag_start": int(d["start"]),
        "frag_end": int(d["end"]),
        "frag_len": int(d["end"]) - int(d["start"]),
    }


def main():
    args = parse_args()
    corpus = Path(args.corpus).resolve()
    if not corpus.exists():
        sys.exit(f"Corpus FASTA not found: {corpus}")

    base = corpus.stem
    out_dir = corpus.parent
    hits_path = out_dir / f"{base}_hits.tsv"
    frags_path = out_dir / f"{base}_fragments.csv"
    orders_path = out_dir / f"{base}_orders.csv"

    # ----- Run BLAST -----
    cmd = [
        "blastn",
        "-task", args.task,
        "-query", str(corpus),
        "-db", args.db,
        "-evalue", str(args.evalue),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend "
                   "sstart send evalue bitscore",
        "-out", str(hits_path),
        "-max_target_seqs", str(args.max_target_seqs),
        "-num_threads", str(args.threads),
    ]
    if args.word_size is not None:
        cmd += ["-word_size", str(args.word_size)]

    print(f"[blast] {' '.join(cmd)}")
    res = subprocess.run(cmd)
    if res.returncode != 0:
        sys.exit(f"blastn failed (code {res.returncode})")

    # ----- Read hits -----
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    if hits_path.stat().st_size == 0:
        hits = pd.DataFrame(columns=cols)
    else:
        hits = pd.read_csv(hits_path, sep="\t", names=cols)
    print(f"[blast] {len(hits):,} raw hits")

    # ----- Build per-fragment index from the corpus FASTA itself -----
    # We need every fragment, not just the ones with hits, so we can compute
    # n_frags and detect orders/fragments with zero hits.
    frag_records = []
    with open(corpus) as f:
        for line in f:
            if line.startswith(">"):
                qid = line[1:].strip().split()[0]
                meta = parse_fragment_id(qid)
                if meta:
                    meta["qseqid"] = qid
                    frag_records.append(meta)
    frags = pd.DataFrame(frag_records)
    if frags.empty:
        sys.exit(f"No parseable fragment IDs in {corpus}. "
                 f"Did this corpus come from build_corpus.py?")

    # ----- Per-fragment: best hit (min evalue, max bitscore) -----
    # This is the "raw measurement" — for each query fragment, what was
    # BLAST's best hit? Used directly for Figure 1 (sensitivity curve).
    if hits.empty:
        per_frag = frags.assign(min_evalue=np.nan, max_bitscore=np.nan, n_hits=0)
    else:
        agg = (hits.groupby("qseqid")
                   .agg(min_evalue=("evalue", "min"),
                        max_bitscore=("bitscore", "max"),
                        n_hits=("sseqid", "size"))
                   .reset_index())
        per_frag = frags.merge(agg, on="qseqid", how="left")
        per_frag["n_hits"] = per_frag["n_hits"].fillna(0).astype(int)

    # Ground truth at the fragment level: hazardous fragments are positives.
    per_frag["truth_hazard"] = per_frag["role"] == "hazard"

    # Save per-fragment CSV (NEW). Drop the redundant qseqid since fragment
    # identity is fully captured by the structured columns.
    frag_cols = ["order_id", "frag_idx", "role", "parent", "frag_start",
                 "frag_end", "frag_len", "n_hits", "min_evalue", "max_bitscore",
                 "truth_hazard"]
    per_frag[frag_cols].to_csv(frags_path, index=False)
    print(f"[done] {len(per_frag):,} fragments -> {frags_path.name}")
    print(f"        hazardous: {per_frag['truth_hazard'].sum()}, "
          f"with any hit: {(per_frag['n_hits'] > 0).sum()}")

    # ----- Per-order: any-flag aggregation -----
    # Order score = min evalue across all fragments (NaN-safe -> Inf if no hits).
    # Order is "flagged" iff at least one fragment had a hit (under the ROC
    # threshold sweep that evaluate.py performs).
    per_frag["min_evalue_filled"] = per_frag["min_evalue"].fillna(np.inf)

    orders = (per_frag.groupby("order_id")
                      .agg(n_frags=("frag_idx", "size"),
                           n_frags_with_hits=("n_hits", lambda x: (x > 0).sum()),
                           total_hits=("n_hits", "sum"),
                           order_min_evalue=("min_evalue_filled", "min"),
                           order_max_bitscore=("max_bitscore", "max"),
                           has_hazard_fragment=("role",
                               lambda x: (x == "hazard").any()))
                      .reset_index())

    # Ground-truth label: order is positive (i.e. evasion) iff it contains any
    # hazard fragment. Honest orders contain only "honest" fragments.
    orders["truth_evasion"] = orders["has_hazard_fragment"]

    orders.to_csv(orders_path, index=False)
    print(f"[done] {len(orders):,} orders -> {orders_path.name}")
    print(f"        positives (truth_evasion=True): {orders['truth_evasion'].sum()}")
    print(f"        flagged (any hit, evalue <= {args.evalue}): "
          f"{(orders['n_frags_with_hits'] > 0).sum()}")


if __name__ == "__main__":
    main()