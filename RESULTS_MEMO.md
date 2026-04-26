# Results memo — BLAST sensitivity benchmark

> Read this before reading the figures or CSVs. It tells you what the experiment measured, what the headline findings are, and what's safe vs unsafe to claim in the writeup.

## TL;DR (one paragraph)

Per-fragment BLAST detection of hazardous DNA degrades smoothly along two axes: shorter fragments and higher adversarial mutation. At the OSTP-mandated 50 bp screening floor and a standard biosecurity cutoff (e ≤ 1e-9), BLAST detects >99% of unmutated hazardous fragments but only **16% of fragments with 20% mutation**. Per-order detection under any-flag aggregation is robust against pure evasion (>99% AUC everywhere) but degrades sharply against diluted evasion at small T with high mutation — at T=50 and 20% mutation, dilute evasion succeeds 40% of the time at 1% FPR. The headline finding is not "find the optimal threshold" (no single threshold works) but rather **the policy threshold needs to be set with explicit reference to assumed adversary capability**, and the gap between unmutated and 20%-mutated adversaries is large enough that the OSTP threshold is appropriate against the former but inadequate against the latter.

## What was measured

For each combination of fragment length T ∈ {20, 30, 50, 75, 100, 150, 200} bp and adversarial mutation rate M ∈ {0, 1, 5, 10, 15, 20}%, we generated three 1000-order corpora and ran BLAST against a curated hazardous DB. Per-fragment results give us BLAST's raw sensitivity profile (Figure 1). Per-order results show what happens when fragments are bundled into realistic orders and aggregated with any-flag (Figure 2).

**Two distinct evasion threat models:**
- **Pure evasion** — every fragment in the order is hazardous (mutated). Tests "max signal, can BLAST find it among many shots?" Under any-flag, any single detected fragment flags the order.
- **Dilute evasion** — one mutated hazardous fragment + 4–19 unmutated benign fillers. Tests "single shot, can BLAST find one needle in benign hay?" The hazardous fragment must be detectable on its own.

The dilute case is the more demanding test of the screener because it removes the redundancy that any-flag aggregation exploits.

## The three findings, in order of strength

### Finding 1 — BLAST per-fragment sensitivity has a clean degradation profile

**The shape:** BLAST sensitivity at any fixed e-value cutoff degrades monotonically as either fragment length shrinks or mutation rate rises. At the loose cutoff (e ≤ 1e-3) BLAST detects unmutated 20 bp fragments nearly perfectly; at the strict cutoff (e ≤ 1e-30) it requires 75 bp even unmutated.

**The numbers (from `fragment_sensitivity.csv`):** Holding the cutoff at e ≤ 1e-9 (typical biosecurity operating point):

| Fragment length | 0% mutation | 5% | 10% | 20% |
|---|---|---|---|---|
| 30 bp | 1.000 | 0.600 | 0.223 | 0.014 |
| 50 bp | 1.000 | 0.993 | 0.822 | 0.161 |
| 75 bp | 1.000 | 1.000 | 0.984 | 0.408 |
| 100 bp | 1.000 | 1.000 | 0.999 | 0.617 |
| 150 bp | 1.000 | 1.000 | 1.000 | 0.818 |
| 200 bp | 1.000 | 1.000 | 1.000 | 0.924 |

**Interpretation:** Mutation matters as much as length. A fragment that's 50 bp at 20% mutation is harder for BLAST to detect (TPR=0.16) than a fragment that's 30 bp at 5% mutation (TPR=0.60). Said differently, **a 4× length advantage can be erased by a 4× mutation increase**.

**Honest fragment FPR is essentially zero across all configurations** (see right column of Figure 1). Mutation doesn't apply to honest fragments, and E. coli K-12 has no meaningful homology to the toxin DB. So when BLAST flags something, it's a real positive — the TPR/FPR gap is the entire story.

### Finding 2 — Any-flag aggregation rescues "pure" attacks but not "dilute" attacks

**The shape:** Per-order ROC under any-flag aggregation stays at AUC=1.00 for nearly every (T, mutation) cell of the pure-evasion grid. Dilute-evasion AUC degrades sharply at small T with high mutation.

**The numbers (from `summary.csv`):** TPR at 1% FPR, dilute evasion only:

| T | 0% mut | 5% | 10% | 15% | 20% |
|---|---|---|---|---|---|
| 20 | 1.000 | 0.789 | 0.490 | 0.241 | 0.123 |
| 30 | 1.000 | 0.966 | 0.804 | 0.513 | 0.284 |
| 50 | 1.000 | 1.000 | 0.984 | 0.850 | 0.596 |
| 75 | 1.000 | 1.000 | 1.000 | 0.957 | 0.774 |
| 100 | 1.000 | 1.000 | 1.000 | 0.991 | 0.892 |

The pure-evasion analog of this table is essentially all 1.000s. **The pure-vs-dilute gap is the operationally meaningful finding.**

**Interpretation:** Any-flag aggregation gives the screener N independent chances to detect the hazardous content. When all N fragments carry hazardous signal (pure), even imperfect per-fragment detection compounds to near-certain order-level detection. When only one fragment carries signal (dilute), per-fragment detection is the only line of defense, and small/mutated fragments slip through. This means **the dilute attack is the realistic threat model for any-flag screeners**, and the OSTP 50 bp threshold catches only 60% of dilute attacks at 1% FPR with 20% mutation.

### Finding 3 — There is no single optimal threshold; there is a frontier

The headline framing for the writeup. The data does not support a claim that "T = X bp is the right threshold." It supports a more nuanced claim:

> The right policy threshold depends on the assumed adversary. Against a naive adversary (no mutation), even T = 20 bp catches everything. Against a realistic adversary (10–20% mutation), T must be ≥ 100 bp to maintain high catch rates against dilute evasion. The OSTP-mandated 50 bp threshold sits in the middle: appropriate against weak adversaries, inadequate against strong ones. This shifts the policy question from "what threshold is best" to "what adversary capability are we defending against, and is the chosen threshold consistent with that assumption?"

This is the strongest framing because it engages with the policy debate as it actually exists rather than pretending there is a magic number.

## What's safe to claim vs not

**Safe to claim:**
- BLAST per-fragment sensitivity degrades smoothly with length and mutation
- Honest FPR is essentially zero against this hazardous DB
- Any-flag aggregation makes pure evasion easy to catch
- Dilute evasion is the realistic attack mode for any-flag screeners
- The OSTP 50 bp threshold is well-chosen against unmutated adversaries
- The OSTP 50 bp threshold catches only ~60% of dilute attacks at 1% FPR when adversary mutation is 20%

**Unsafe to claim:**
- "BLAST is bad" — Figure 1 shows BLAST works extremely well in most regimes
- "Threshold X is optimal" — data argues against a single optimum
- "Commec would fail similarly" — commec adds taxonomic + HMM layers we haven't tested
- "Screening is broken" — the failure mode is specific (small T, high mutation, dilute) not general

## Caveats the writeup should disclose

1. **Honest source diversity.** Single bacterial genome (E. coli K-12). Real customer orders span more taxa. FPR is plausibly higher in reality, which would worsen all numbers.
2. **Mutation model.** Uniform-random substitution. Real adversaries can do better with synonymous substitutions (preserve protein function while disrupting nucleotide alignment). Our 20% rate is probably easier on BLAST than a 5% codon-optimized adversary would be.
3. **Hazardous DB size and composition.** 19 records, all ≤10 kb, treated as nucleotide. Production screeners use thousands of records, often via translated-protein search. Numbers are illustrative of the methodology, not directly transferable to commec or commercial screeners.
4. **No cross-order experiment.** Adversaries can split fragments across N orders to defeat any per-order screener. Quantifying this is v2.
5. **BLAST only.** Commec comparison is the natural next study.

## Files to read

- **`fig1_fragment_sensitivity.png`** — the main visual finding
- **`fig2_order_roc_panels.png`** — the deployment-scenario visual; large grid, focus on T=20–50 columns at high mutation
- **`summary.csv`** — the table behind Figure 2; filter to interesting cells
- **`fragment_sensitivity.csv`** — the table behind Figure 1
- **`roc_data.csv`** — raw threshold sweep data if anyone wants to make custom plots
- **`hazardous_manifest.tsv`** — what's actually in the DB (audit trail)
