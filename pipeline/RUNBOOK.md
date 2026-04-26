# Pipeline runbook

Four scripts that together produce the headline ROC figure for your project. Each
is a separate concern, communicates via files on disk, and can be re-run independently.

## File map

```
scripts/
  fetch_data.py     # Pulls E. coli K-12 reference from NCBI (honest source)
  build_corpus.py   # Generates honest + evasion order corpora at one T value
  screen.py         # Runs BLAST on a corpus, applies any-flag aggregation
  evaluate.py       # Combines all T values, computes ROC, plots the headline figure

data/
  honest_source.fasta             # written by fetch_data.py
  hazardous.fasta                 # YOU PROVIDE (use placeholder for dry run)
  corpora/
    T20/
      honest.fasta                # written by build_corpus.py
      evasion_pure.fasta
      evasion_dilute.fasta
      corpus_meta.json
      honest_orders.csv           # written by screen.py
      evasion_pure_orders.csv
      evasion_dilute_orders.csv
      *_hits.tsv                  # raw BLAST output (kept for debugging)
    T30/  T50/  T75/  T100/  T150/  T200/

results/
  roc_data.csv      # written by evaluate.py
  summary.csv
  roc_panels.png    # the headline figure
```

## End-to-end run

All commands are run from the `scripts/` folder with `(blast-env)` active.

### Step 0 — One-time setup

Place your hazardous FASTA at `../data/hazardous.fasta`. For a dry run, copy
the placeholder:

```
mkdir ..\data
copy hazardous_PLACEHOLDER.fasta ..\data\hazardous.fasta
```

### Step 1 — Fetch the honest source genome

```
python fetch_data.py
```

Downloads E. coli K-12 to `../data/honest_source.fasta` (~5 MB). Idempotent.

### Step 2 — Build the hazardous BLAST DB

```
makeblastdb -in ..\data\hazardous.fasta -dbtype nucl -out hazardous_db -title "Hazardous reference"
```

This writes `hazardous_db.{nhr,nin,nsq,...}` into the current `scripts/` folder.

### Step 3 — Generate corpora at every T value

Run once per T. Loop in PowerShell or copy-paste:

```
python build_corpus.py --T 20
python build_corpus.py --T 30
python build_corpus.py --T 50
python build_corpus.py --T 75
python build_corpus.py --T 100
python build_corpus.py --T 150
python build_corpus.py --T 200
```

Each run writes three FASTAs (honest, evasion_pure, evasion_dilute) and a
metadata JSON to `../data/corpora/T{T}/`. With defaults (200 orders × 5-20
fragments), each call produces ~2,500 fragments per file and runs in seconds.

### Step 4 — Screen each corpus against the hazardous DB

Three calls per T (one per corpus file). For T=50:

```
python screen.py --corpus ..\data\corpora\T50\honest.fasta --db hazardous_db
python screen.py --corpus ..\data\corpora\T50\evasion_pure.fasta --db hazardous_db
python screen.py --corpus ..\data\corpora\T50\evasion_dilute.fasta --db hazardous_db
```

Repeat for each T. Each call writes `*_hits.tsv` and `*_orders.csv` next to
the input. BLAST runs in seconds for small T, up to ~1 minute for T=200 with
the full evasion corpus.

(There's a one-liner PowerShell loop at the bottom of this file.)

### Step 5 — Evaluate and plot

```
python evaluate.py
```

Auto-discovers all `T*/` subdirs, computes per-T ROC for both evasion types,
prints a summary table, writes `roc_panels.png` to `../results/`.

## Key knobs you might want to turn

### In screen.py
- `--task blastn-short` is the default (best for small fragments). Switch to
  `blastn` to see how the default preset performs.
- `--evalue 10` is permissive on purpose - we want to see all potential hits and
  do thresholding ourselves in evaluate.py. Tightening this here just truncates
  your ROC curve on the left.
- `--word-size 7` would push detection of very short fragments harder.

### In build_corpus.py
- `--n-orders 500` for tighter confidence intervals (slower).
- `--frags-per-order-min 1 --frags-per-order-max 1` to test per-fragment behavior
  (any-flag becomes trivial; useful as a sanity check).

### In evaluate.py
- `--fpr-targets 0.001,0.01,0.05` to see TPR at very low FPR (the operating
  region a real screener cares about).

## PowerShell convenience: screen all T values at once

```powershell
20,30,50,75,100,150,200 | ForEach-Object {
    $T = $_
    Write-Host "=== T=$T ==="
    python screen.py --corpus "..\data\corpora\T$T\honest.fasta" --db hazardous_db
    python screen.py --corpus "..\data\corpora\T$T\evasion_pure.fasta" --db hazardous_db
    python screen.py --corpus "..\data\corpora\T$T\evasion_dilute.fasta" --db hazardous_db
}
```

## What to look at in the output

- **`results/roc_panels.png`** — your headline figure. One panel per T, each with
  two ROC curves (pure and dilute evasion vs honest baseline). The diagonal is
  random-chance; the further northwest a curve lies, the better.
- **`results/summary.csv`** — AUC + TPR-at-fixed-FPR for each (T, evasion type).
  This is the table you'd paste into a paper or slide.
- **`results/roc_data.csv`** — the raw (threshold, FPR, TPR) points if you want
  to make a custom plot.

## Sanity-check expectations (with placeholder data)

Don't read too much into the dry-run numbers - the placeholder hazardous file
is full of repeats and won't behave like real toxin genes. But you should at
least see:

- `summary.csv` has 14 rows (7 T values × 2 evasion types)
- `roc_panels.png` has 7 panels
- AUCs for `evasion_pure` are higher than `evasion_dilute` at every T (more
  hazard fragments = more chances to be caught under any-flag)
- Both AUCs degrade as T shrinks

If those qualitative patterns hold, the pipeline is wired correctly. Real
hazardous sequences will give you the real numbers.
