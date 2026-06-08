#!/usr/bin/env python3
"""
Extract variants from a consensus VCF (no VEP annotation required).

Produces the same TSV schema as db_vep_vcf_to_variants_all.py so the
output is compatible with the BigQuery variants table. Annotation columns
(gene_symbol, mane_select, hgvsc, variant_impact, clinvar_clnsig,
clinvar_alleleid) are left empty — they are populated downstream by the
database annotation layer.

Output columns:
    sample_id, assay_type, gene_symbol, mane_select, chrom, pos, ref, alt,
    variant_id, gt, alt_allele_count, is_hom_alt, hgvsc, variant_impact,
    clinvar_clnsig, clinvar_alleleid
"""

import argparse
import gzip
import sys
from typing import Dict, Optional, Tuple


def open_maybe_gz(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def gt_to_alt_count(gt: str) -> Tuple[Optional[int], int]:
    if gt in (".", "./.", ".|."):
        return None, 0
    sep = "|" if "|" in gt else "/"
    alleles = gt.split(sep)
    if any(a == "." for a in alleles):
        return None, 0
    alt_count = sum(1 for a in alleles if a != "0")
    is_hom_alt = 1 if len(alleles) == 2 and alleles[0] != "0" and alleles[0] == alleles[1] else 0
    return alt_count, is_hom_alt


HEADER = (
    "sample_id\tassay_type\tgene_symbol\tmane_select\tchrom\tpos\tref\talt\t"
    "variant_id\tgt\talt_allele_count\tis_hom_alt\thgvsc\tvariant_impact\t"
    "clinvar_clnsig\tclinvar_alleleid\n"
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True, help="Sample ID")
    ap.add_argument("--assay",  required=True, help="Assay type")
    ap.add_argument("--vcf",    required=True, help="Input VCF (plain or .gz)")
    ap.add_argument("--out",    required=False, help="Output TSV path (default: <sample>_<assay>_variants.tsv)")
    args = ap.parse_args()

    out_path = args.out or f"{args.sample}_{args.assay}_variants.tsv"

    n_written = 0
    with open_maybe_gz(args.vcf) as fh, open(out_path, "w") as out:
        out.write(HEADER)

        for line in fh:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue

            chrom, pos, _id, ref, alt_str, _qual, flt, _info, fmt = cols[:9]
            sample_field = cols[9]

            # Keep PASS and unfiltered variants only
            if flt not in ("PASS", ".", ""):
                continue

            # Strip chr prefix for consistency
            if chrom.startswith("chr"):
                chrom = chrom[3:]

            # Parse genotype
            fmt_keys = fmt.split(":")
            fmt_vals = sample_field.split(":")
            fmt_map: Dict[str, str] = {
                k: (fmt_vals[i] if i < len(fmt_vals) else "")
                for i, k in enumerate(fmt_keys)
            }
            gt = fmt_map.get("GT", "./.")
            alt_count, is_hom_alt = gt_to_alt_count(gt)

            # Skip non-variant calls
            if alt_count == 0:
                continue

            # One row per ALT allele
            for alt in alt_str.split(","):
                if alt in (".", "*"):
                    continue

                variant_id = f"{chrom}-{pos}-{ref}-{alt}"
                alt_count_str = "" if alt_count is None else str(alt_count)

                out.write(
                    f"{args.sample}\t{args.assay}\t"   # sample_id, assay_type
                    f"\t"                               # gene_symbol  (no VEP)
                    f"\t"                               # mane_select  (no VEP)
                    f"{chrom}\t{pos}\t{ref}\t{alt}\t"
                    f"{variant_id}\t{gt}\t"
                    f"{alt_count_str}\t{is_hom_alt}\t"
                    f"\t"                               # hgvsc         (no VEP)
                    f"\t"                               # variant_impact (no VEP)
                    f"\t"                               # clinvar_clnsig (no VEP)
                    f"\n"                               # clinvar_alleleid (no VEP)
                )
                n_written += 1

    print(f"Wrote {n_written} variant row(s) to {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
