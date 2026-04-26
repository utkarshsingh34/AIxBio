"""
build_corpus.py

Generate three order corpora at a given fragment-size T:
  - honest:        fragments randomly sampled from a non-hazardous source
  - evasion_pure:  ALL fragments from one hazardous parent per order
  - evasion_dilute: ONE hazardous fragment + benign filler per order

Each "order" is a small bundle of fragments (5-20 by default). Output is one
FASTA file per corpus, with fragment IDs that encode order_id + role + parent.

NEW IN v1.5: --mutation-rate applies random substitutions to HAZARDOUS fragments
(only) before writing them. Honest and filler fragments are never mutated. This
models an adversary who introduces variation to evade exact-match screening.
The mutation rate is per-base independent uniform random substitution to a
different nucleotide.

USAGE (from your pipeline folder, with blast-env active):
    # Baseline (no mutation)
    python build_corpus.py --T 50

    # 5% substitution rate on hazardous fragments
    python build_corpus.py --T 50 --mutation-rate 0.05

    # Different output directory (default appends _mut{rate}% to T{N})
    python build_corpus.py --T 50 --mutation-rate 0.05 --out-dir ../data/corpora/T50_mut5

The fragment IDs follow the format:
    order{ID}__frag{IDX}__role{honest|hazard|filler}__parent{NAME}__{START}_{END}
so the screening + analysis steps can recover order membership and ground truth
without consulting a separate metadata file.
"""

import argparse
import json
import random
import sys
from pathlib import Path

from Bio import SeqIO


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--T", type=int, required=True,
                   help="Fragment size in bp (all fragments will be this length)")
    p.add_argument("--n-orders", type=int, default=200,
                   help="Number of orders per corpus (default: %(default)s)")
    p.add_argument("--frags-per-order-min", type=int, default=5)
    p.add_argument("--frags-per-order-max", type=int, default=20)
    p.add_argument("--honest-source", default="../data/honest_source.fasta",
                   help="FASTA of non-hazardous sequences to sample honest fragments from")
    p.add_argument("--hazardous", default="../data/hazardous.fasta",
                   help="FASTA of hazardous parent sequences for evasion orders")
    p.add_argument("--out-dir", default=None,
                   help="Output directory (default: ../data/corpora/T{T}_mut{rate}%%)")
    p.add_argument("--mutation-rate", type=float, default=0.0,
                   help="Per-base random substitution rate for HAZARDOUS fragments only "
                        "(0.0 = no mutation, 0.05 = 5%%, etc.). Default: %(default)s")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def load_fasta(path):
    p = Path(path).resolve()
    if not p.exists():
        sys.exit(f"FASTA not found: {p}")
    records = list(SeqIO.parse(str(p), "fasta"))
    if not records:
        sys.exit(f"FASTA is empty: {p}")
    return records


def random_slice(seq_str, length, rng):
    if len(seq_str) < length:
        return None
    start = rng.randint(0, len(seq_str) - length)
    return start, start + length, seq_str[start:start + length]


def mutate(seq_str, rate, rng):
    """Apply random substitutions at the given per-base rate.
    Each mutated position gets a uniformly chosen DIFFERENT nucleotide."""
    if rate <= 0:
        return seq_str
    bases = "ACGT"
    out = []
    for ch in seq_str:
        if rng.random() < rate:
            cu = ch.upper()
            if cu in bases:
                # Pick a different base
                alts = [b for b in bases if b != cu]
                new = rng.choice(alts)
                # Preserve case
                out.append(new if ch.isupper() else new.lower())
            else:
                # Non-standard char (e.g. N) — leave alone rather than guess
                out.append(ch)
        else:
            out.append(ch)
    return "".join(out)


def make_id(order_id, frag_idx, role, parent_id, start, end):
    safe_parent = (parent_id.replace(" ", "_").replace("|", "_")
                            .replace(":", "_").replace(".", "_"))
    return (f"order{order_id:05d}__frag{frag_idx:02d}__role{role}"
            f"__parent{safe_parent}__{start}_{end}")


def build_honest(out_path, honest_records, n_orders, frags_min, frags_max, T, rng):
    n_frags = 0
    with open(out_path, "w") as f:
        for order_id in range(n_orders):
            n_frags_this_order = rng.randint(frags_min, frags_max)
            for frag_idx in range(n_frags_this_order):
                for _ in range(20):
                    rec = rng.choice(honest_records)
                    sl = random_slice(str(rec.seq), T, rng)
                    if sl is not None:
                        s, e, frag = sl
                        fid = make_id(order_id, frag_idx, "honest", rec.id, s, e)
                        f.write(f">{fid}\n{frag}\n")
                        n_frags += 1
                        break
    return n_frags


def build_evasion_pure(out_path, hazard_records, n_orders, frags_min, frags_max,
                       T, mutation_rate, rng):
    """Each order = N random T-bp slices from ONE hazardous parent, with mutations."""
    n_frags = 0
    skipped_orders = 0
    with open(out_path, "w") as f:
        for order_id in range(n_orders):
            eligible = [r for r in hazard_records if len(r.seq) >= T]
            if not eligible:
                skipped_orders += 1
                continue
            parent = rng.choice(eligible)
            n_frags_this_order = rng.randint(frags_min, frags_max)
            for frag_idx in range(n_frags_this_order):
                sl = random_slice(str(parent.seq), T, rng)
                if sl is None:
                    continue
                s, e, frag = sl
                frag = mutate(frag, mutation_rate, rng)  # ← mutation applied here
                fid = make_id(order_id, frag_idx, "hazard", parent.id, s, e)
                f.write(f">{fid}\n{frag}\n")
                n_frags += 1
    if skipped_orders:
        print(f"  [warn] skipped {skipped_orders} pure-evasion orders "
              f"(no hazardous parent >= {T} bp)")
    return n_frags


def build_evasion_dilute(out_path, hazard_records, honest_records, n_orders,
                         frags_min, frags_max, T, mutation_rate, rng):
    """Each order = 1 hazardous fragment (mutated) + (N-1) honest fragments (not mutated)."""
    n_frags = 0
    skipped_orders = 0
    with open(out_path, "w") as f:
        for order_id in range(n_orders):
            eligible = [r for r in hazard_records if len(r.seq) >= T]
            if not eligible:
                skipped_orders += 1
                continue
            n_frags_this_order = rng.randint(frags_min, frags_max)
            hazard_pos = rng.randint(0, n_frags_this_order - 1)

            for frag_idx in range(n_frags_this_order):
                if frag_idx == hazard_pos:
                    parent = rng.choice(eligible)
                    sl = random_slice(str(parent.seq), T, rng)
                    if sl is None:
                        continue
                    s, e, frag = sl
                    frag = mutate(frag, mutation_rate, rng)  # ← mutation applied here
                    fid = make_id(order_id, frag_idx, "hazard", parent.id, s, e)
                else:
                    for _ in range(20):
                        rec = rng.choice(honest_records)
                        sl = random_slice(str(rec.seq), T, rng)
                        if sl is not None:
                            break
                    if sl is None:
                        continue
                    s, e, frag = sl
                    fid = make_id(order_id, frag_idx, "filler", rec.id, s, e)
                f.write(f">{fid}\n{frag}\n")
                n_frags += 1
    if skipped_orders:
        print(f"  [warn] skipped {skipped_orders} dilute-evasion orders "
              f"(no hazardous parent >= {T} bp)")
    return n_frags


def main():
    args = parse_args()
    rng = random.Random(args.seed)

    # Default output dir reflects both T and mutation rate
    if args.out_dir:
        out_dir = Path(args.out_dir).resolve()
    else:
        mut_pct = int(round(args.mutation_rate * 100))
        out_dir = Path(f"../data/corpora/T{args.T}_mut{mut_pct}").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] honest source: {args.honest_source}")
    honest = load_fasta(args.honest_source)
    print(f"  {len(honest)} record(s), total {sum(len(r.seq) for r in honest):,} bp")

    print(f"[load] hazardous: {args.hazardous}")
    hazard = load_fasta(args.hazardous)
    print(f"  {len(hazard)} record(s), total {sum(len(r.seq) for r in hazard):,} bp")
    eligible_count = sum(1 for r in hazard if len(r.seq) >= args.T)
    print(f"  {eligible_count} hazardous record(s) >= {args.T} bp")

    print(f"[build] T={args.T}, n_orders={args.n_orders}, "
          f"frags_per_order={args.frags_per_order_min}-{args.frags_per_order_max}, "
          f"mutation_rate={args.mutation_rate}")
    counts = {}

    h_path = out_dir / "honest.fasta"
    counts["honest"] = build_honest(h_path, honest, args.n_orders,
                                    args.frags_per_order_min, args.frags_per_order_max,
                                    args.T, rng)
    print(f"  honest:          {counts['honest']:>6} fragments  -> {h_path.name}")

    ep_path = out_dir / "evasion_pure.fasta"
    counts["evasion_pure"] = build_evasion_pure(
        ep_path, hazard, args.n_orders,
        args.frags_per_order_min, args.frags_per_order_max,
        args.T, args.mutation_rate, rng)
    print(f"  evasion_pure:    {counts['evasion_pure']:>6} fragments  -> {ep_path.name}")

    ed_path = out_dir / "evasion_dilute.fasta"
    counts["evasion_dilute"] = build_evasion_dilute(
        ed_path, hazard, honest, args.n_orders,
        args.frags_per_order_min, args.frags_per_order_max,
        args.T, args.mutation_rate, rng)
    print(f"  evasion_dilute:  {counts['evasion_dilute']:>6} fragments  -> {ed_path.name}")

    meta = {
        "T": args.T,
        "mutation_rate": args.mutation_rate,
        "n_orders": args.n_orders,
        "frags_per_order_min": args.frags_per_order_min,
        "frags_per_order_max": args.frags_per_order_max,
        "seed": args.seed,
        "honest_source": str(Path(args.honest_source).resolve()),
        "hazardous": str(Path(args.hazardous).resolve()),
        "fragment_counts": counts,
    }
    with open(out_dir / "corpus_meta.json", "w") as f:
        json.dump(meta, f, indent=2)
    print(f"[done] Metadata at {out_dir / 'corpus_meta.json'}")


if __name__ == "__main__":
    main()