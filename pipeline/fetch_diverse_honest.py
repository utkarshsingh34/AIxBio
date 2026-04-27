"""
fetch_diverse_honest.py

Replaces the single-source E. coli honest FASTA with a multi-source consolidated
FASTA covering several homology classes (bacterial, eukaryotic, viral, plant).
This makes false-positive rate measurements meaningful: a real synthesis customer
might order from any of these taxa, so honest fragments should reflect that
diversity rather than being drawn from one bacterial genome.

USAGE (from pipeline folder, with blast-env active):
    python fetch_diverse_honest.py
    python fetch_diverse_honest.py --email you@example.com   # NCBI prefers a real email
    python fetch_diverse_honest.py --out ../data/honest_source.fasta

Output:
    ../data/honest_sources/{tag}.fasta    one file per source (cached, skipped if already there)
    ../data/honest_source.fasta           consolidated FASTA used by build_corpus.py
    ../data/honest_manifest.tsv           per-record audit trail

The script is idempotent: per-source downloads are cached. Delete a per-source
file to force re-download. Delete the consolidated FASTA to force re-merge.
"""

import argparse
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO


# Curated source list. Picked to span homology classes that matter for FPR:
#   - 3 bacteria from different phyla (different operon structures, distant
#     homology to virulence genes)
#   - 1 archaeon (more diverse coding patterns)
#   - 1 yeast (eukaryotic, broad customer-order relevance)
#   - 1 plant (plant biology orders)
#   - 1 mammalian transcript subset (medical / research customer orders)
#   - 2 phages (mobile elements, sometimes share content with toxin-bearing prophages)
#   - 1 generic synthetic biology vector (typical lab customer order)
#
# All accessions verified as RefSeq or curated.
SOURCES = [
    # tag                    accession         description
    ("ecoli_k12",           "NC_000913.3",    "E. coli K-12 MG1655 reference genome"),
    ("bsubtilis_168",       "NC_000964.3",    "B. subtilis 168 reference genome"),
    ("paeruginosa_pao1",    "NC_002516.2",    "P. aeruginosa PAO1 reference genome"),
    ("mtuberculosis_h37rv", "NC_000962.3",    "M. tuberculosis H37Rv reference genome"),
    ("scerevisiae_chr1",    "NC_001133.9",    "S. cerevisiae chromosome 1"),
    ("athaliana_chr1",      "NC_003070.9",    "A. thaliana chromosome 1"),
    ("phage_lambda",        "NC_001416.1",    "Enterobacteria phage lambda"),
    ("phage_t4",            "AF158101.6",     "Enterobacteria phage T4"),
    ("pUC19_vector",        "L09137.2",       "pUC19 cloning vector"),
    ("human_chr22",         "NC_000022.11",   "Human chromosome 22 (compact)"),
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--email", default="hackathon@example.com",
                   help="Email for NCBI Entrez (they ask for one; not validated)")
    p.add_argument("--out", default="../data/honest_source.fasta",
                   help="Consolidated output FASTA path (default: %(default)s)")
    p.add_argument("--sources-dir", default="../data/honest_sources",
                   help="Per-source download directory (default: %(default)s)")
    p.add_argument("--manifest", default="../data/honest_manifest.tsv",
                   help="Output TSV manifest (default: %(default)s)")
    p.add_argument("--max-mb-per-source", type=float, default=None,
                   help="Optional cap on per-source size in MB. If set, large records "
                        "are truncated to this many MB. Useful for keeping the consolidated "
                        "FASTA from being dominated by one huge source like human chr22.")
    return p.parse_args()


def fetch_one(tag, accession, description, sources_dir, email):
    """Download a single source. Idempotent — skips if already cached."""
    out_path = sources_dir / f"{tag}.fasta"
    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"  [cached] {tag}: {out_path.stat().st_size:,} bytes")
        return out_path

    Entrez.email = email
    print(f"  [fetch]  {tag} ({accession}) ...", end="", flush=True)
    t0 = time.time()
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession,
                               rettype="fasta", retmode="text")
        rec = SeqIO.read(handle, "fasta")
        handle.close()
        # Rewrite the description to include our tag so downstream tools can trace origin
        rec.id = f"{tag}__{rec.id}"
        rec.description = f"{tag} | {description} | {rec.description}"
        SeqIO.write([rec], str(out_path), "fasta")
        elapsed = time.time() - t0
        print(f" {len(rec.seq):,} bp in {elapsed:.1f}s")
        return out_path
    except Exception as e:
        print(f" FAILED ({e})")
        return None


def main():
    args = parse_args()
    out_path = Path(args.out).resolve()
    sources_dir = Path(args.sources_dir).resolve()
    manifest_path = Path(args.manifest).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sources_dir.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    max_bytes = int(args.max_mb_per_source * 1_000_000) if args.max_mb_per_source else None

    print(f"[fetch] Downloading {len(SOURCES)} honest source(s) to {sources_dir}")
    fetched = []
    for tag, accession, description in SOURCES:
        # NCBI asks for ~3 second pause between requests when not using API key
        path = fetch_one(tag, accession, description, sources_dir, args.email)
        if path:
            fetched.append((tag, accession, description, path))
        time.sleep(0.4)  # be nice to NCBI

    if not fetched:
        sys.exit("No sources downloaded successfully. Check your network / NCBI status.")

    print(f"\n[merge] Consolidating into {out_path}")
    manifest_rows = []
    total_bp = 0
    n_records = 0
    with open(out_path, "w") as out_fh:
        for tag, accession, description, path in fetched:
            for rec in SeqIO.parse(str(path), "fasta"):
                seq_str = str(rec.seq)
                if max_bytes and len(seq_str) > max_bytes:
                    print(f"  [trim]  {tag}: {len(seq_str):,} bp -> {max_bytes:,} bp")
                    seq_str = seq_str[:max_bytes]
                out_fh.write(f">{rec.id} {rec.description}\n{seq_str}\n")
                manifest_rows.append({
                    "source_tag": tag,
                    "accession": accession,
                    "description": description,
                    "record_id": rec.id,
                    "length_bp": len(seq_str),
                })
                n_records += 1
                total_bp += len(seq_str)

    with open(manifest_path, "w") as mf:
        mf.write("source_tag\taccession\tdescription\trecord_id\tlength_bp\n")
        for row in manifest_rows:
            mf.write(f"{row['source_tag']}\t{row['accession']}\t{row['description']}"
                     f"\t{row['record_id']}\t{row['length_bp']}\n")

    print(f"\n[done] {n_records} record(s), {total_bp:,} bp total -> {out_path}")
    print(f"       per-source manifest: {manifest_path}")
    print(f"\n=== Per-source breakdown ===")
    by_source = {}
    for row in manifest_rows:
        by_source.setdefault(row["source_tag"], 0)
        by_source[row["source_tag"]] += row["length_bp"]
    for tag, bp in sorted(by_source.items(), key=lambda x: -x[1]):
        pct = 100 * bp / total_bp if total_bp else 0
        print(f"  {tag:25} {bp:>12,} bp  ({pct:5.1f}%)")


if __name__ == "__main__":
    main()
