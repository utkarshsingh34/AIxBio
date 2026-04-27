# Results

This directory contains two complete result sets from the project. Both are referenced in the writeup.

## `v1_5/` — Primary results

The results reported in the main figures and tables of the writeup. Configuration:
- Honest source: E. coli K-12 MG1655 only (4.6 Mb)
- Hazardous DB: 19 records, ≤10 kb cap (~250 kb total)
- 1000 orders per (T, mutation) cell

These numbers are slightly more conservative than the v2 sweep and have the cleanest interpretability (toxin-CDS-scale records only).

## `v2_robustness/` — Robustness check (appendix)

Re-run with expanded data, included as an appendix robustness check. Configuration:
- Honest source: 10 sources spanning bacterial, viral, eukaryotic, plant, synthetic origins (~100 Mb)
- Hazardous DB: 40 records, ≤50 kb cap (~700 kb-1 MB total)
- 1000 orders per (T, mutation) cell

Qualitative findings preserved across both configurations. Absolute TPR numbers are slightly lower in v2 due to two confounded effects: (1) a larger BLAST DB increases e-value penalties for the same alignment, and (2) the 50 kb cap admits records with more non-coding content. Disentangling these effects is left to future work.

## File layout in each subdirectory

- `fig1_fragment_sensitivity.png` — Per-fragment TPR/FPR vs T × mutation
- `fig2_order_roc_panels.png` — Per-order ROC grid
- `fig3_evasion_heatmap.png` — Pure vs dilute heatmap (headline visual)
- `fig4_operating_points.png` — Catch rates at OSTP threshold across adversary types
- `fragment_sensitivity.csv` — Tidy data behind Figure 1
- `roc_data.csv` — Threshold sweep data behind Figure 2
- `summary.csv` — AUC + TPR @ fixed FPR per (T, mutation, evasion type)