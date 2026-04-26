"""
fetch_data.py

Pulls the honest-source genome (E. coli K-12 MG1655, NC_000913.3) from NCBI
and writes it to disk. Run once before building corpora.

USAGE (from your scripts\ folder, with blast-env active):
    python fetch_data.py
    python fetch_data.py --email you@example.com   # NCBI prefers you set this
    python fetch_data.py --out ../data/honest_source.fasta

The script is idempotent: if the output file already exists and is non-empty,
it prints a notice and exits. Delete the file to force re-download.
"""

import argparse
import sys
from pathlib import Path

from Bio import Entrez, SeqIO


HONEST_ACCESSION = "NC_000913.3"  # E. coli K-12 MG1655 reference genome


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--email", default="hackathon@example.com",
                   help="Email for NCBI Entrez (they ask for one; not validated)")
    p.add_argument("--out", default="../data/honest_source.fasta",
                   help="Output FASTA path (default: %(default)s)")
    p.add_argument("--accession", default=HONEST_ACCESSION,
                   help="NCBI accession to fetch (default: %(default)s)")
    return p.parse_args()


def main():
    args = parse_args()
    out = Path(args.out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    if out.exists() and out.stat().st_size > 0:
        print(f"[skip] {out} already exists ({out.stat().st_size:,} bytes). "
              f"Delete to re-download.")
        return

    Entrez.email = args.email
    print(f"[fetch] Downloading {args.accession} from NCBI ...")
    try:
        handle = Entrez.efetch(db="nucleotide", id=args.accession,
                               rettype="fasta", retmode="text")
        rec = SeqIO.read(handle, "fasta")
        handle.close()
    except Exception as e:
        sys.exit(f"NCBI fetch failed: {e}")

    SeqIO.write([rec], str(out), "fasta")
    print(f"[done] Wrote {out} ({len(rec.seq):,} bp, id={rec.id})")


if __name__ == "__main__":
    main()
