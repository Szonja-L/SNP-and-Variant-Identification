"""
generate_synthetic_data.py
==========================
Generates biologically realistic synthetic data for the E. coli SNP/Variant
Identification portfolio project.

Simulates outputs from a Galaxy pipeline:
  FastQC → Trimmomatic → BWA-MEM → SAMtools → FreeBayes → SnpEff

Biological scenario
-------------------
Sample   : E. coli clinical isolate (UTI patient, NL hospital, 2023)
Reference: E. coli K-12 MG1655 (NCBI: NC_000913.3 / U00096.3)
Platform : Illumina NovaSeq 6000, 150bp paired-end
Goal     : Identify mutations conferring antibiotic resistance
"""

import numpy as np
import pandas as pd
import json
import os

np.random.seed(42)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

# ─────────────────────────────────────────────────────────────────────────────
# 1.  QC METRICS  (FastQC-style before and after Trimmomatic)
# ─────────────────────────────────────────────────────────────────────────────

def generate_per_base_quality():
    """Simulate per-base sequence quality for 150bp reads (Phred scores)."""
    positions = np.arange(1, 151)
    
    # Raw reads: quality drops at the ends (realistic Illumina degradation)
    raw_mean = np.where(
        positions < 10,
        28 + positions * 0.6,            # ramp up from ~28 at pos 1
        np.where(
            positions < 130,
            36 - (positions - 10) * 0.015,  # slight decline in middle
            36 - (positions - 10) * 0.015 - (positions - 130) * 0.25  # 3′ drop
        )
    )
    raw_q1  = raw_mean - np.random.uniform(3, 5, 150)
    raw_q3  = raw_mean + np.random.uniform(1, 3, 150)
    raw_q10 = raw_mean - np.random.uniform(6, 9, 150)
    raw_q90 = raw_mean + np.random.uniform(2, 4, 150)
    raw_q10 = np.clip(raw_q10, 2, 40)
    raw_q90 = np.clip(raw_q90, 10, 41)

    # Trimmed reads: cleaner, no terminal drops
    trim_mean = np.where(
        positions <= 137,
        np.where(positions < 10,
                 32 + positions * 0.4,
                 37.5 - (positions - 10) * 0.008),
        37.5 - (137 - 10) * 0.008 - (positions - 137) * 0.3
    )
    trim_q1  = trim_mean - np.random.uniform(1.5, 3, 150)
    trim_q3  = trim_mean + np.random.uniform(0.5, 1.5, 150)
    trim_q10 = trim_mean - np.random.uniform(3, 5, 150)
    trim_q90 = trim_mean + np.random.uniform(1, 2.5, 150)
    trim_q10 = np.clip(trim_q10, 20, 41)
    trim_q90 = np.clip(trim_q90, 30, 41)

    df = pd.DataFrame({
        "position": positions,
        "raw_mean":  np.clip(raw_mean,  2, 40),
        "raw_q10":   np.clip(raw_q10,  2, 40),
        "raw_q25":   np.clip(raw_q1,   2, 40),
        "raw_q75":   np.clip(raw_q3,   2, 40),
        "raw_q90":   np.clip(raw_q90,  2, 40),
        "trim_mean": np.clip(trim_mean, 2, 41),
        "trim_q10":  np.clip(trim_q10, 2, 41),
        "trim_q25":  np.clip(trim_q1,  2, 41),
        "trim_q75":  np.clip(trim_q3,  2, 41),
        "trim_q90":  np.clip(trim_q90, 2, 41),
    })
    return df


def generate_gc_content():
    """Simulate per-sequence GC content distribution (E. coli ~50.7% GC)."""
    gc_values = np.arange(0, 101)
    # Theoretical normal centred at 50.7
    theoretical = 100 * np.exp(-0.5 * ((gc_values - 50.7) / 6.5) ** 2)
    theoretical /= theoretical.sum() / 100
    # Actual with slight right skew
    actual = 100 * np.exp(-0.5 * ((gc_values - 50.7) / 7.2) ** 2)
    actual[45:55] *= np.random.uniform(0.97, 1.04, 10)
    actual /= actual.sum() / 100
    return pd.DataFrame({"gc_pct": gc_values, "theoretical": theoretical, "actual": actual})


def generate_qc_summary():
    """Top-line QC numbers for the report."""
    stats = {
        "sample_id":               "EC_UTI_2023_001",
        "sequencer":               "Illumina NovaSeq 6000",
        "read_length":             150,
        "raw_total_reads":         3_142_088,
        "raw_r1_reads":            1_571_044,
        "raw_r2_reads":            1_571_044,
        "raw_pct_q30":             87.4,
        "raw_pct_gc":              50.7,
        "raw_mean_quality":        35.8,
        "trim_total_reads":        3_056_442,
        "trim_surviving_pct":      97.3,
        "trim_pct_q30":            93.1,
        "trim_mean_quality":       37.2,
        "adapters_removed_reads":  85_646,
        "low_quality_removed":     10_221,
        "trimmomatic_version":     "0.39",
        "fastqc_version":          "0.12.1",
    }
    return stats


# ─────────────────────────────────────────────────────────────────────────────
# 2.  ALIGNMENT STATS  (BWA-MEM + SAMtools flagstat / coverage)
# ─────────────────────────────────────────────────────────────────────────────

def generate_coverage_data():
    """Simulate per-window read depth across the 4.64 Mbp E. coli genome."""
    genome_size = 4_641_652          # NC_000913.3
    window       = 5_000             # 5 kb windows
    n_windows    = genome_size // window

    positions = np.arange(n_windows) * window + window // 2

    # Base coverage ~48.7x; local variation driven by GC content & mappability
    base_cov  = 48.7
    gc_bias   = 3.0 * np.sin(2 * np.pi * positions / 200_000) + \
                1.5 * np.cos(2 * np.pi * positions / 750_000)
    noise     = np.random.normal(0, 4.0, n_windows)
    coverage  = np.clip(base_cov + gc_bias + noise, 0, 120)

    # Add a low-coverage region (rRNA operon – common in assemblies)
    rrna_start = int(4_166_659 / window)
    coverage[rrna_start:rrna_start + 10] *= 0.35

    # Highlight gene regions (approximate positions of key resistance genes)
    key_genes = {
        "gyrA":  (2_337_072, "#FF6B6B"),
        "rpoB":  (4_185_061, "#4ECDC4"),
        "marR":  (1_623_580, "#FFE66D"),
        "acrA":  (487_397,   "#A8DADC"),
        "tolC":  (3_177_130, "#C77DFF"),
    }

    df = pd.DataFrame({"position": positions, "coverage": coverage})
    return df, key_genes


def generate_alignment_stats():
    """SAMtools flagstat-equivalent statistics."""
    stats = {
        "total_reads":           3_056_442,
        "mapped_reads":          3_012_609,
        "mapping_rate_pct":      98.57,
        "properly_paired_pct":   97.21,
        "singletons":            11_482,
        "mean_coverage_x":       48.7,
        "median_coverage_x":     47.3,
        "pct_genome_gt_0x":      99.87,
        "pct_genome_gt_20x":     98.94,
        "pct_genome_gt_50x":     44.12,
        "insert_size_mean_bp":   372,
        "insert_size_sd_bp":     48,
        "bwa_mem_version":       "0.7.17",
        "samtools_version":      "1.18",
        "reference":             "NC_000913.3 (E. coli K-12 MG1655, 4,641,652 bp)",
    }
    return stats


# ─────────────────────────────────────────────────────────────────────────────
# 3.  VARIANTS  (FreeBayes VCF + SnpEff annotation)
# ─────────────────────────────────────────────────────────────────────────────

VARIANTS = [
    # CHROM      POS        REF  ALT  QUAL  DP  AF   GENE  AA_CHANGE   EFFECT              FUNCTIONAL_CLASS   CLINICAL_NOTE
    ("NC_000913.3", 2_337_074, "G",  "A",  998, 52, 0.981, "gyrA",  "p.Asp87Asn (D87N)", "MISSENSE_VARIANT",   "NONSYNONYMOUS", "Fluoroquinolone resistance; quinolone-resistance-determining region (QRDR)"),
    ("NC_000913.3", 4_185_063, "C",  "T",  967, 61, 1.000, "rpoB",  "p.Ser531Leu (S531L)","MISSENSE_VARIANT",  "NONSYNONYMOUS", "Rifampicin resistance; RNA polymerase beta-subunit hotspot"),
    ("NC_000913.3", 1_623_580, "G",  "A",  823, 44, 0.955, "marR",  "p.Gly103Ser (G103S)", "MISSENSE_VARIANT", "NONSYNONYMOUS", "Activates marA/marB; upregulates AcrAB-TolC efflux pump (MDR)"),
    ("NC_000913.3",   487_404, "A",  "G",  791, 58, 0.966, "acrA",  "p.Ile355Val (I355V)", "MISSENSE_VARIANT", "NONSYNONYMOUS", "Efflux pump component; may enhance MDR efflux capacity"),
    ("NC_000913.3", 3_177_140, "T",  "C",  612, 39, 0.923, "tolC",  "p.Leu457Pro (L457P)", "MISSENSE_VARIANT", "NONSYNONYMOUS", "Outer-membrane channel of AcrAB-TolC; altered specificity"),
    ("NC_000913.3", 2_337_102, "C",  "T",  743, 55, 0.982, "gyrA",  "p.Ala90Val (A90V)",   "MISSENSE_VARIANT", "NONSYNONYMOUS", "Secondary QRDR mutation; synergistic fluoroquinolone resistance"),
    ("NC_000913.3",   892_114, "A",  "T",  541, 47, 0.936, "ompC",  "p.Lys163* (K163*)",   "STOP_GAINED",      "NONSENSE",      "Truncation of OmpC porin; reduces outer-membrane permeability"),
    ("NC_000913.3", 1_304_018, "G",  "C",  489, 43, 0.907, "dnaQ",  "p.Pro72Arg (P72R)",   "MISSENSE_VARIANT", "NONSYNONYMOUS", "DNA polymerase III proofreading subunit; may increase mutation rate"),
    # Silent / synonymous
    ("NC_000913.3",   251_012, "C",  "T",  667, 51, 0.961, "recA",  "p.Leu276Leu (L276L)", "SYNONYMOUS_VARIANT","SYNONYMOUS",   "Silent; recA codon wobble position"),
    ("NC_000913.3", 3_892_445, "G",  "A",  589, 48, 0.979, "rpsL",  "p.Ala157Ala (A157A)", "SYNONYMOUS_VARIANT","SYNONYMOUS",   "Silent; 30S ribosomal protein S12"),
    # Intergenic
    ("NC_000913.3",   104_839, "T",  "A",  312, 38, 0.868, "intergenic", "n/a",             "INTERGENIC_VARIANT","MODIFIER",     "Intergenic; possible regulatory effect on nearby ompT"),
    ("NC_000913.3", 2_198_067, "G",  "T",  287, 34, 0.853, "intergenic", "n/a",             "INTERGENIC_VARIANT","MODIFIER",     "Intergenic; upstream of fimA fimbrial operon"),
    # UTR
    ("NC_000913.3", 4_038_912, "C",  "A",  401, 42, 0.905, "rpsA",  "5'UTR",               "5_PRIME_UTR_VARIANT","MODIFIER",    "5′-UTR of rpsA; ribosomal protein S1"),
]


def build_variants_df():
    cols = ["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AF",
            "GENE", "AA_CHANGE", "EFFECT", "FUNCTIONAL_CLASS", "CLINICAL_NOTE"]
    df = pd.DataFrame(VARIANTS, columns=cols)
    df["FILTER"] = "PASS"
    df["VAF_PCT"] = (df["AF"] * 100).round(1)
    return df


def write_vcf(df):
    """Write a properly formatted VCF 4.2 file."""
    header = """##fileformat=VCFv4.2
##fileDate=20231015
##source=FreeBayes v1.3.6
##reference=NC_000913.3
##contig=<ID=NC_000913.3,length=4641652,species="Escherichia coli K-12 MG1655">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
##INFO=<ID=ANN,Number=.,Type=String,Description="SnpEff annotation (AA_CHANGE|EFFECT|FUNCTIONAL_CLASS)">
##FILTER=<ID=PASS,Description="Variant passes all quality filters">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele frequency">
##SnpEffVersion="5.1d (build 2022-09-09)"
##SnpEffCmd="SnpEff -v Escherichia_coli_k_12"
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEC_UTI_2023_001
"""
    lines = [header]
    for _, row in df.iterrows():
        ann = f"{row['AA_CHANGE']}|{row['EFFECT']}|{row['FUNCTIONAL_CLASS']}|{row['GENE']}"
        info = f"DP={row['DP']};AF={row['AF']:.3f};GENE={row['GENE']};ANN={ann}"
        fmt  = f"GT:DP:AF\t1/1:{row['DP']}:{row['AF']:.3f}"
        lines.append(
            f"{row['CHROM']}\t{row['POS']}\t.\t{row['REF']}\t{row['ALT']}\t"
            f"{row['QUAL']}\t{row['FILTER']}\t{info}\t{fmt}\n"
        )
    return "".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# 4.  SAMPLE FASTQ  (tiny synthetic snippet for format demonstration)
# ─────────────────────────────────────────────────────────────────────────────

FASTQ_SNIPPET = """\
@EC_UTI_2023_001.1 1/1
ATGCGATCGTAGCTAGCTAGCGATCGATCGATGCGATCGATGCTAGCTAGCTAGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@EC_UTI_2023_001.2 2/1
GCTAGCTAGCGATCGATCGATGCGATCGATGCTAGCTAGCTAGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@EC_UTI_2023_001.3 3/1
CGATCGATCGATGCTAGCTAGCTAGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
+
IIIIIIHHHHHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIH
"""


# ─────────────────────────────────────────────────────────────────────────────
# 5.  WRITE ALL FILES
# ─────────────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # QC data
    qc_df = generate_per_base_quality()
    qc_df.to_csv(os.path.join(DATA_DIR, "per_base_quality.csv"), index=False)

    gc_df = generate_gc_content()
    gc_df.to_csv(os.path.join(DATA_DIR, "gc_content.csv"), index=False)

    qc_summary = generate_qc_summary()
    with open(os.path.join(DATA_DIR, "qc_summary.json"), "w") as f:
        json.dump(qc_summary, f, indent=2)

    # Alignment data
    cov_df, key_genes = generate_coverage_data()
    cov_df.to_csv(os.path.join(DATA_DIR, "genome_coverage.csv"), index=False)

    align_stats = generate_alignment_stats()
    with open(os.path.join(DATA_DIR, "alignment_stats.json"), "w") as f:
        json.dump(align_stats, f, indent=2)

    # Variants
    var_df = build_variants_df()
    var_df.to_csv(os.path.join(RESULTS_DIR, "variants_annotated.csv"), index=False)

    vcf_text = write_vcf(var_df)
    with open(os.path.join(RESULTS_DIR, "variants.vcf"), "w") as f:
        f.write(vcf_text)

    # Resistance genes only
    resist_df = var_df[var_df["FUNCTIONAL_CLASS"] == "NONSYNONYMOUS"].copy()
    resist_df.to_csv(os.path.join(RESULTS_DIR, "resistance_variants.csv"), index=False)

    # Sample FASTQ
    with open(os.path.join(DATA_DIR, "EC_UTI_2023_001_R1_example.fastq"), "w") as f:
        f.write(FASTQ_SNIPPET)

    # Key genes JSON (for other scripts)
    with open(os.path.join(DATA_DIR, "key_gene_positions.json"), "w") as f:
        json.dump({k: {"pos": v[0], "color": v[1]} for k, v in key_genes.items()}, f, indent=2)

    print(f"✓ QC data          → {DATA_DIR}/per_base_quality.csv")
    print(f"✓ GC content       → {DATA_DIR}/gc_content.csv")
    print(f"✓ QC summary       → {DATA_DIR}/qc_summary.json")
    print(f"✓ Coverage data    → {DATA_DIR}/genome_coverage.csv")
    print(f"✓ Alignment stats  → {DATA_DIR}/alignment_stats.json")
    print(f"✓ Annotated VCF    → {RESULTS_DIR}/variants.vcf")
    print(f"✓ Variants CSV     → {RESULTS_DIR}/variants_annotated.csv")
    print(f"✓ Resistance CSV   → {RESULTS_DIR}/resistance_variants.csv")
    print(f"✓ Sample FASTQ     → {DATA_DIR}/EC_UTI_2023_001_R1_example.fastq")
    print(f"\nTotal variants: {len(var_df)}")
    print(f"  Nonsynonymous : {(var_df.FUNCTIONAL_CLASS == 'NONSYNONYMOUS').sum()}")
    print(f"  Synonymous    : {(var_df.FUNCTIONAL_CLASS == 'SYNONYMOUS').sum()}")
    print(f"  Modifier      : {(var_df.FUNCTIONAL_CLASS == 'MODIFIER').sum()}")


if __name__ == "__main__":
    main()
