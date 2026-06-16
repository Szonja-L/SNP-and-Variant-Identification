#!/usr/bin/env python3
"""
variant_analysis_pipeline.py
=============================
Standalone Python script that mirrors the analytical steps performed
in Galaxy, but runs entirely on the synthetic results data.

Use this script to reproduce the variant analysis, generate summary
statistics, and verify figure outputs WITHOUT needing a Galaxy server.

Pipeline steps reproduced here
--------------------------------
1.  Load and validate variants from annotated VCF (or CSV)
2.  Quality filter   (QUAL >= 200, DP >= 10, AF >= 0.90)
3.  Functional classification (SnpEff ANN field parsing)
4.  Resistance gene annotation  (curated database)
5.  Summary statistics + export
6.  All figures via generate_figures.py

Usage
-----
    cd snp-variant-identification/
    python scripts/variant_analysis_pipeline.py

    # With a real VCF from Galaxy:
    python scripts/variant_analysis_pipeline.py --vcf path/to/real.vcf

Requirements
------------
    pip install pandas numpy matplotlib seaborn
"""

import argparse
import os
import sys
import json
import textwrap
from datetime import datetime

import numpy as np
import pandas as pd

BASE_DIR    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR    = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
FIG_DIR     = os.path.join(BASE_DIR, "figures")

# ─────────────────────────────────────────────────────────────────────────────
# 1.  RESISTANCE DATABASE  (curated from CARD / ResFinder literature)
# ─────────────────────────────────────────────────────────────────────────────

RESISTANCE_DB = {
    # Gene → {mutation: (drug_class, mechanism, MIC_fold)}
    "gyrA": {
        "D87N": ("Fluoroquinolone", "Altered target (DNA gyrase QRDR)", 32),
        "A90V": ("Fluoroquinolone", "Altered target (DNA gyrase QRDR)", 16),
        "S83L": ("Fluoroquinolone", "Altered target (DNA gyrase QRDR)", 8),
    },
    "rpoB": {
        "S531L": ("Rifampicin", "Altered target (RNA polymerase β-subunit)", 256),
        "H526Y": ("Rifampicin", "Altered target (RNA polymerase β-subunit)", 128),
        "D516V": ("Rifampicin", "Altered target (RNA polymerase β-subunit)", 64),
    },
    "marR": {
        "G103S": ("Multidrug (MarA-mediated)", "Regulatory derepression → AcrAB-TolC upregulation", 8),
        "Y137H": ("Multidrug (MarA-mediated)", "Regulatory derepression → AcrAB-TolC upregulation", 4),
    },
    "acrA": {
        "I355V": ("Multidrug (efflux)", "Efflux pump structural change (MFP subunit)", 4),
    },
    "tolC": {
        "L457P": ("Multidrug (efflux)", "Altered outer-membrane channel specificity", 6),
    },
    "ompC": {
        "K163*": ("Broad spectrum", "Loss-of-function porin (reduced outer-membrane permeability)", 4),
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# 2.  LOAD VARIANTS
# ─────────────────────────────────────────────────────────────────────────────

def load_vcf_to_df(vcf_path: str) -> pd.DataFrame:
    """Parse a FreeBayes/SnpEff VCF into a tidy DataFrame."""
    records = []
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, pos, vid, ref, alt, qual, flt, info, fmt = parts[:9]
            sample = parts[9] if len(parts) > 9 else "."
            # Parse INFO
            info_d = {}
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_d[k] = v
            # Parse ANN (SnpEff)
            ann_parts = info_d.get("ANN", "").split("|")
            aa_change = ann_parts[0] if ann_parts else ""
            effect    = ann_parts[1] if len(ann_parts) > 1 else ""
            func_cls  = ann_parts[2] if len(ann_parts) > 2 else ""
            gene      = info_d.get("GENE", ann_parts[3] if len(ann_parts) > 3 else "")
            records.append({
                "CHROM":            chrom,
                "POS":              int(pos),
                "REF":              ref,
                "ALT":              alt,
                "QUAL":             float(qual) if qual != "." else np.nan,
                "FILTER":           flt,
                "DP":               int(info_d.get("DP", 0)),
                "AF":               float(info_d.get("AF", 0)),
                "GENE":             gene,
                "AA_CHANGE":        aa_change,
                "EFFECT":           effect,
                "FUNCTIONAL_CLASS": func_cls,
            })
    return pd.DataFrame(records)


def load_synthetic_csv() -> pd.DataFrame:
    path = os.path.join(RESULTS_DIR, "variants_annotated.csv")
    return pd.read_csv(path)


# ─────────────────────────────────────────────────────────────────────────────
# 3.  QUALITY FILTER
# ─────────────────────────────────────────────────────────────────────────────

def quality_filter(df: pd.DataFrame,
                   min_qual: float = 200,
                   min_dp: int = 10,
                   min_af: float = 0.90) -> pd.DataFrame:
    """Apply hard quality filters."""
    before = len(df)
    df = df[
        (df["QUAL"] >= min_qual) &
        (df["DP"]   >= min_dp)   &
        (df["AF"]   >= min_af)
    ].copy()
    print(f"  Quality filter: {before} → {len(df)} variants retained "
          f"(QUAL≥{min_qual}, DP≥{min_dp}, AF≥{min_af})")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 4.  ANNOTATE WITH RESISTANCE DATABASE
# ─────────────────────────────────────────────────────────────────────────────

def annotate_resistance(df: pd.DataFrame) -> pd.DataFrame:
    """Cross-reference variants with curated resistance DB."""
    drug_classes, mechanisms, mic_folds, is_resist = [], [], [], []
    for _, row in df.iterrows():
        gene = row.get("GENE", "")
        aa   = row.get("AA_CHANGE", "")
        # Extract short mutation code from HGVS-like string e.g. "p.Asp87Asn (D87N)" → "D87N"
        mut_code = ""
        if "(" in aa and ")" in aa:
            mut_code = aa[aa.index("(") + 1 : aa.index(")")]

        db_entry = RESISTANCE_DB.get(gene, {}).get(mut_code, None)
        if db_entry:
            drug_classes.append(db_entry[0])
            mechanisms.append(db_entry[1])
            mic_folds.append(db_entry[2])
            is_resist.append(True)
        else:
            drug_classes.append("N/A")
            mechanisms.append("N/A")
            mic_folds.append(None)
            is_resist.append(False)
    df = df.copy()
    df["DRUG_CLASS"]  = drug_classes
    df["MECHANISM"]   = mechanisms
    df["MIC_FOLD"]    = mic_folds
    df["IS_RESISTANCE"] = is_resist
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 5.  SUMMARY STATISTICS
# ─────────────────────────────────────────────────────────────────────────────

def compute_summary(df: pd.DataFrame, align_stats: dict, qc_stats: dict) -> dict:
    fc_counts = df["FUNCTIONAL_CLASS"].value_counts().to_dict()
    resist_df = df[df["IS_RESISTANCE"] == True]

    summary = {
        "run_timestamp":         datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "sample_id":             qc_stats.get("sample_id", "EC_UTI_2023_001"),
        "reference":             align_stats.get("reference", "NC_000913.3"),
        "total_variants":        len(df),
        "nonsynonymous":         fc_counts.get("NONSYNONYMOUS", 0),
        "synonymous":            fc_counts.get("SYNONYMOUS", 0),
        "modifier":              fc_counts.get("MODIFIER", 0),
        "nonsense":              fc_counts.get("NONSENSE", 0),
        "resistance_variants":   len(resist_df),
        "resistance_genes":      resist_df["GENE"].nunique(),
        "drug_classes_affected": resist_df["DRUG_CLASS"].nunique(),
        "mean_qual":             round(df["QUAL"].mean(), 1),
        "mean_depth":            round(df["DP"].mean(), 1),
        "mean_vaf":              round(df["AF"].mean(), 3),
        "mapping_rate_pct":      align_stats.get("mapping_rate_pct"),
        "mean_coverage_x":       align_stats.get("mean_coverage_x"),
        "raw_reads":             qc_stats.get("raw_total_reads"),
        "trimmed_reads":         qc_stats.get("trim_total_reads"),
        "pct_q30_raw":           qc_stats.get("raw_pct_q30"),
        "pct_q30_trimmed":       qc_stats.get("trim_pct_q30"),
    }
    return summary


# ─────────────────────────────────────────────────────────────────────────────
# 6.  PRETTY PRINT REPORT
# ─────────────────────────────────────────────────────────────────────────────

def print_report(df: pd.DataFrame, summary: dict):
    SEP = "─" * 72

    print(f"\n{SEP}")
    print("  SNP & VARIANT ANALYSIS REPORT")
    print(f"  Sample : {summary['sample_id']}")
    print(f"  Ref    : {summary['reference']}")
    print(f"  Date   : {summary['run_timestamp']}")
    print(SEP)

    print("\n  ▶  QC & ALIGNMENT STATISTICS")
    print(f"     Raw reads          : {summary['raw_reads']:>12,}")
    print(f"     Trimmed reads      : {summary['trimmed_reads']:>12,}")
    print(f"     %Q30 (raw)         : {summary['pct_q30_raw']:>11.1f}%")
    print(f"     %Q30 (trimmed)     : {summary['pct_q30_trimmed']:>11.1f}%")
    print(f"     Mapping rate       : {summary['mapping_rate_pct']:>11.2f}%")
    print(f"     Mean coverage      : {summary['mean_coverage_x']:>10.1f}×")

    print("\n  ▶  VARIANT SUMMARY")
    print(f"     Total variants     : {summary['total_variants']:>12}")
    print(f"     Nonsynonymous      : {summary['nonsynonymous']:>12}")
    print(f"     Synonymous         : {summary['synonymous']:>12}")
    print(f"     Modifier           : {summary['modifier']:>12}")
    print(f"     Nonsense (stop)    : {summary['nonsense']:>12}")
    print(f"     Mean QUAL          : {summary['mean_qual']:>12.1f}")
    print(f"     Mean read depth    : {summary['mean_depth']:>12.1f}×")
    print(f"     Mean VAF           : {summary['mean_vaf']:>12.1%}")

    print("\n  ▶  RESISTANCE FINDINGS")
    print(f"     Resistance variants: {summary['resistance_variants']:>12}")
    print(f"     Resistance genes   : {summary['resistance_genes']:>12}")
    print(f"     Drug classes hit   : {summary['drug_classes_affected']:>12}")

    resist_df = df[df["IS_RESISTANCE"] == True].sort_values("QUAL", ascending=False)
    if len(resist_df):
        print()
        header = f"  {'GENE':<8} {'MUTATION':<14} {'DRUG CLASS':<30} {'QUAL':>6} {'VAF':>7}"
        print(header)
        print(f"  {'─'*72}")
        for _, row in resist_df.iterrows():
            mut = row["AA_CHANGE"].split("(")[-1].rstrip(")") if "(" in row["AA_CHANGE"] else row["AA_CHANGE"]
            print(f"  {row['GENE']:<8} {mut:<14} {row['DRUG_CLASS']:<30} {row['QUAL']:>6.0f} {row['AF']:>7.1%}")

    print(f"\n{SEP}\n")


# ─────────────────────────────────────────────────────────────────────────────
# 7.  MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Bacterial SNP & Variant Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Examples
        --------
          # Run with synthetic demonstration data:
          python scripts/variant_analysis_pipeline.py

          # Run with a real VCF from Galaxy:
          python scripts/variant_analysis_pipeline.py --vcf results/my_variants.vcf

          # Adjust quality thresholds:
          python scripts/variant_analysis_pipeline.py --min-qual 300 --min-af 0.95
        """),
    )
    parser.add_argument("--vcf",      default=None,  help="Path to annotated VCF file (optional)")
    parser.add_argument("--min-qual", default=200.0, type=float, help="Minimum QUAL score (default: 200)")
    parser.add_argument("--min-dp",   default=10,    type=int,   help="Minimum read depth (default: 10)")
    parser.add_argument("--min-af",   default=0.90,  type=float, help="Minimum allele frequency (default: 0.90)")
    parser.add_argument("--figures",  action="store_true",       help="Regenerate all figures")
    args = parser.parse_args()

    print("\n" + "═" * 72)
    print("  E. coli SNP & VARIANT IDENTIFICATION PIPELINE")
    print("  Portfolio project — NGS Quality Control & Variant Calling")
    print("═" * 72)

    # Load alignment & QC stats
    with open(os.path.join(DATA_DIR, "alignment_stats.json")) as f:
        align_stats = json.load(f)
    with open(os.path.join(DATA_DIR, "qc_summary.json")) as f:
        qc_stats = json.load(f)

    print("\n[1/5] Loading variants…")
    if args.vcf:
        print(f"      Source: {args.vcf}")
        df = load_vcf_to_df(args.vcf)
    else:
        print("      Source: synthetic data (results/variants_annotated.csv)")
        df = load_synthetic_csv()
    print(f"      Loaded {len(df)} raw variant calls.")

    print("\n[2/5] Applying quality filters…")
    df_filt = quality_filter(df, args.min_qual, args.min_dp, args.min_af)

    print("\n[3/5] Annotating with resistance database…")
    df_annot = annotate_resistance(df_filt)
    n_resist = df_annot["IS_RESISTANCE"].sum()
    print(f"      Found {n_resist} variants matching known resistance mutations.")

    print("\n[4/5] Computing summary statistics…")
    summary = compute_summary(df_annot, align_stats, qc_stats)

    # Write outputs
    out_path = os.path.join(RESULTS_DIR, "pipeline_results.csv")
    df_annot.to_csv(out_path, index=False)
    print(f"      Results written to: {out_path}")

    summary_path = os.path.join(RESULTS_DIR, "pipeline_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"      Summary written to: {summary_path}")

    print("\n[5/5] Printing report…")
    print_report(df_annot, summary)

    if args.figures:
        print("Regenerating figures…")
        import subprocess
        subprocess.run([sys.executable,
                        os.path.join(BASE_DIR, "scripts", "generate_figures.py")])

    print("Pipeline complete.\n")


if __name__ == "__main__":
    main()
