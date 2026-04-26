"""
aggregate_hazardous.py

Walk a folder of .gb (GenBank) and .fasta files, extract every nucleotide
sequence record, and write them all to one consolidated FASTA ready for
makeblastdb. Renames record IDs to avoid collisions across source files
(important when files share generic names like "sequence.gb" / "sequence.fasta"
from NCBI default downloads).

USAGE (from pipeline folder, with blast-env active):
    python aggregate_hazardous.py
    python aggregate_hazardous.py --in-dir ../data/hazardous_data --out ../data/hazardous.fasta
    python aggregate_hazardous.py --include-protein   # also include protein records (off by default)

Default behavior:
- Reads every .gb / .gbk / .genbank / .fasta / .fa / .fna file under --in-dir
- Treats EVERY file independently (no FASTA/GenBank pairing logic; files like
  "sequence.gb" and "sequence (1).fasta" are independent records, not pairs)
- Skips records whose sequence content is missing/undefined (common in partial
  GenBank exports that have annotations but no ORIGIN block)
- Keeps only nucleotide records (skips proteins by default)
- Renames IDs to {filename_stem}__{original_id} to prevent collisions
- Writes one consolidated FASTA + a manifest TSV listing what came from where
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError


GENBANK_EXTS = {".gb", ".gbk", ".genbank"}
FASTA_EXTS = {".fasta", ".fa", ".fna", ".ffn"}
NUCL_ALPHABET = set("ATGCNUatgcnu-")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-dir", default="../data/hazardous_data",
                   help="Directory containing .gb / .fasta files (default: %(default)s)")
    p.add_argument("--out", default="../data/hazardous.fasta",
                   help="Output consolidated FASTA (default: %(default)s)")
    p.add_argument("--manifest", default="../data/hazardous_manifest.tsv",
                   help="Output TSV listing each record's source (default: %(default)s)")
    p.add_argument("--include-protein", action="store_true",
                   help="Also include protein records (default: nucleotide only)")
    p.add_argument("--min-length", type=int, default=20,
                   help="Skip records shorter than this many bp (default: %(default)s)")
    return p.parse_args()


def classify_seq(seq_str: str) -> str:
    """Return 'nucl' if sequence consists only of nucleotide characters, else 'prot'."""
    if not seq_str:
        return "empty"
    letters = set(seq_str) - set(" \t\r\n")
    return "nucl" if letters <= NUCL_ALPHABET else "prot"


def safe_id_part(text: str) -> str:
    """Make a string safe for FASTA IDs (no whitespace, no pipes, no colons)."""
    return (text.replace(" ", "_").replace("|", "_")
                .replace(":", "_").replace(",", "_")
                .replace(";", "_").replace("\t", "_")
                .replace("(", "").replace(")", ""))


def iter_input_files(in_dir: Path):
    """Yield (path, format) tuples for every parseable file under in_dir."""
    for p in sorted(in_dir.rglob("*")):
        if not p.is_file():
            continue
        ext = p.suffix.lower()
        if ext in GENBANK_EXTS:
            yield p, "genbank"
        elif ext in FASTA_EXTS:
            yield p, "fasta"


def main():
    args = parse_args()
    in_dir = Path(args.in_dir).resolve()
    out_path = Path(args.out).resolve()
    manifest_path = Path(args.manifest).resolve()

    if not in_dir.is_dir():
        sys.exit(f"Input directory not found: {in_dir}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    files = list(iter_input_files(in_dir))
    if not files:
        sys.exit(f"No .gb / .fasta files found under {in_dir}")

    print(f"[scan] {len(files)} input file(s) under {in_dir}")
    for p, fmt in files:
        print(f"  {fmt:>8}  {p.relative_to(in_dir)}")

    n_kept = 0
    n_skipped_short = 0
    n_skipped_protein = 0
    n_skipped_empty = 0
    n_skipped_undefined = 0
    seen_ids = Counter()
    manifest_rows = []

    with open(out_path, "w") as out_fh:
        for path, fmt in files:
            stem = safe_id_part(path.stem)
            try:
                records = list(SeqIO.parse(str(path), fmt))
            except Exception as e:
                print(f"  [warn] could not parse {path.name}: {e}")
                continue

            for rec in records:
                # Defensive: GenBank records can have undefined sequences (records
                # with annotation only and no ORIGIN block). Catch and skip cleanly.
                try:
                    seq_str = str(rec.seq)
                except UndefinedSequenceError:
                    n_skipped_undefined += 1
                    print(f"  [skip] {path.name} record {rec.id}: "
                          f"sequence undefined (no ORIGIN block?)")
                    continue
                except Exception as e:
                    n_skipped_undefined += 1
                    print(f"  [skip] {path.name} record {rec.id}: "
                          f"could not extract sequence ({e})")
                    continue

                kind = classify_seq(seq_str)
                if kind == "empty":
                    n_skipped_empty += 1
                    continue
                if kind == "prot" and not args.include_protein:
                    n_skipped_protein += 1
                    continue
                if len(seq_str) < args.min_length:
                    n_skipped_short += 1
                    continue

                # Build a collision-free, BLAST-safe ID. With NCBI default
                # filenames like "sequence (3).gb", the parens cause issues so
                # safe_id_part strips them; we still get a unique ID via the
                # stem (sequence_3) plus the original record ID.
                orig_id = safe_id_part(rec.id or "noid")
                new_id = f"{stem}__{orig_id}"
                if seen_ids[new_id] > 0:
                    new_id = f"{new_id}__{seen_ids[new_id] + 1}"
                seen_ids[new_id] += 1

                description = (rec.description or "").replace("\n", " ").strip()
                if description and len(description) > 150:
                    description = description[:147] + "..."

                header = new_id if not description else f"{new_id} {description}"
                out_fh.write(f">{header}\n{seq_str}\n")
                n_kept += 1

                manifest_rows.append({
                    "new_id": new_id,
                    "source_file": str(path.relative_to(in_dir)),
                    "format": fmt,
                    "kind": kind,
                    "length_bp": len(seq_str),
                    "original_id": orig_id,
                })

    with open(manifest_path, "w") as mf:
        mf.write("new_id\tsource_file\tformat\tkind\tlength_bp\toriginal_id\n")
        for row in manifest_rows:
            mf.write(f"{row['new_id']}\t{row['source_file']}\t{row['format']}\t"
                     f"{row['kind']}\t{row['length_bp']}\t{row['original_id']}\n")

    print(f"\n[done] {n_kept} records written to {out_path}")
    print(f"       manifest: {manifest_path}")
    if n_skipped_undefined:
        print(f"       skipped {n_skipped_undefined} record(s) with undefined sequence")
    if n_skipped_protein:
        print(f"       skipped {n_skipped_protein} protein record(s) "
              f"(use --include-protein to keep them)")
    if n_skipped_short:
        print(f"       skipped {n_skipped_short} record(s) shorter than {args.min_length} bp")
    if n_skipped_empty:
        print(f"       skipped {n_skipped_empty} empty record(s)")

    if manifest_rows:
        lengths = [r["length_bp"] for r in manifest_rows]
        print(f"\n=== Output summary ===")
        print(f"  records: {len(lengths)}")
        print(f"  total bp: {sum(lengths):,}")
        print(f"  shortest: {min(lengths):,} bp")
        print(f"  longest:  {max(lengths):,} bp")
        print(f"  median:   {sorted(lengths)[len(lengths) // 2]:,} bp")


if __name__ == "__main__":
    main()