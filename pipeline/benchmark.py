"""
benchmark.py

Measures how the pipeline scales with corpus size. Runs build_corpus.py +
screen.py at several --n-orders values, records wall-clock time per stage,
and reports a scaling fit.

USAGE (from pipeline\ folder, blast-env active):
    python benchmark.py
    python benchmark.py --T 50 --sizes 50,100,200,500,1000
    python benchmark.py --db hazardous_db --sizes 100,200,500

Output: prints a table to stdout, writes benchmark_results.csv next to this
script. Scaling exponent close to 1.0 = linear; > 1.0 = super-linear (bad);
< 1.0 = sub-linear (good, usually means startup overhead dominates).
"""

import argparse
import shutil
import subprocess
import time
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--T", type=int, default=50,
                   help="Fragment size to benchmark at (default: %(default)s)")
    p.add_argument("--sizes", default="50,100,200,500,1000",
                   help="Comma-separated n_orders values to test")
    p.add_argument("--db", default="hazardous_db", help="BLAST DB prefix")
    p.add_argument("--corpora", choices=["honest", "evasion_pure", "evasion_dilute", "all"],
                   default="all", help="Which corpora to time")
    p.add_argument("--out", default="benchmark_results.csv")
    return p.parse_args()


def time_cmd(cmd, cwd=None):
    """Run cmd, return wall-clock seconds (and crash if it fails)."""
    t0 = time.perf_counter()
    res = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    if res.returncode != 0:
        print(res.stdout)
        print(res.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return elapsed


def fit_power_law(sizes, times):
    """Fit times = a * sizes^b. Return (a, b). b is the scaling exponent."""
    import math
    logs = [math.log(s) for s in sizes]
    logt = [math.log(t) for t in times]
    n = len(sizes)
    mean_logs = sum(logs) / n
    mean_logt = sum(logt) / n
    num = sum((logs[i] - mean_logs) * (logt[i] - mean_logt) for i in range(n))
    den = sum((logs[i] - mean_logs) ** 2 for i in range(n))
    b = num / den if den > 0 else float("nan")
    a = math.exp(mean_logt - b * mean_logs)
    return a, b


def main():
    args = parse_args()
    here = Path(__file__).resolve().parent
    sizes = [int(x) for x in args.sizes.split(",") if x.strip()]
    corpora = (["honest", "evasion_pure", "evasion_dilute"]
               if args.corpora == "all" else [args.corpora])

    # Use a sandbox so we don't trample the real corpora
    sandbox = here / "_benchmark_sandbox"
    if sandbox.exists():
        shutil.rmtree(sandbox)
    sandbox.mkdir()

    rows = []
    for n in sizes:
        out_dir = sandbox / f"T{args.T}_n{n}"
        print(f"\n--- n_orders={n} ---")

        # Build
        build_cmd = [
            "python", "build_corpus.py",
            "--T", str(args.T),
            "--n-orders", str(n),
            "--out-dir", str(out_dir),
        ]
        t_build = time_cmd(build_cmd, cwd=str(here))
        print(f"  build: {t_build:.2f}s")

        for corpus in corpora:
            corpus_path = out_dir / f"{corpus}.fasta"
            screen_cmd = [
                "python", "screen.py",
                "--corpus", str(corpus_path),
                "--db", args.db,
            ]
            t_screen = time_cmd(screen_cmd, cwd=str(here))
            print(f"  screen[{corpus}]: {t_screen:.2f}s")

            # How many fragments did we actually screen?
            n_frags = sum(1 for line in open(corpus_path) if line.startswith(">"))

            rows.append({
                "n_orders": n,
                "corpus": corpus,
                "n_fragments": n_frags,
                "build_seconds": t_build,
                "screen_seconds": t_screen,
                "total_seconds": t_build + t_screen,
            })

    df = pd.DataFrame(rows)
    df.to_csv(here / args.out, index=False)

    print("\n=== Per-corpus scaling fits (screen time vs n_orders) ===")
    print("Exponent ~1.0 = linear; >1.0 = super-linear (bad)")
    for corpus in corpora:
        sub = df[df["corpus"] == corpus]
        if len(sub) >= 2:
            a, b = fit_power_law(sub["n_orders"].tolist(),
                                 sub["screen_seconds"].tolist())
            print(f"  {corpus:>15}: exponent = {b:.2f}  "
                  f"(predicted screen time at n=10000: {a * 10000 ** b:.1f}s)")

    print("\n=== Raw measurements ===")
    print(df.to_string(index=False))

    # Cleanup
    shutil.rmtree(sandbox)
    print(f"\n[done] CSV: {here / args.out}")


if __name__ == "__main__":
    main()