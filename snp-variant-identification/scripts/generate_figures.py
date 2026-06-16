"""
generate_figures.py
===================
Produces all publication-quality figures for the SNP & Variant
Identification portfolio project.

Figures generated
-----------------
1. per_base_quality.png    – FastQC-style before/after Trimmomatic
2. gc_content.png          – Per-sequence GC content distribution
3. genome_coverage.png     – Read depth across the 4.64 Mbp E. coli genome
4. variant_vaf.png         – Variant allele frequency distribution
5. variant_annotation.png  – Stacked bar: variant effects by category
6. resistance_summary.png  – Clinical resistance gene panel
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import os, json

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR    = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
FIG_DIR     = os.path.join(BASE_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# ── Colour palette (dark terminal / sequencing aesthetic) ──────────────────
BG          = "#0D1117"
PANEL_BG    = "#161B22"
GRID_COL    = "#21262D"
TEXT_MAIN   = "#E6EDF3"
TEXT_DIM    = "#7D8590"
CYAN        = "#39D0D8"
LIME        = "#7EE787"
AMBER       = "#E3B341"
CORAL       = "#F78166"
VIOLET      = "#C778DD"
BLUE        = "#58A6FF"

plt.rcParams.update({
    "figure.facecolor":  BG,
    "axes.facecolor":    PANEL_BG,
    "axes.edgecolor":    GRID_COL,
    "axes.labelcolor":   TEXT_MAIN,
    "axes.titlecolor":   TEXT_MAIN,
    "xtick.color":       TEXT_DIM,
    "ytick.color":       TEXT_DIM,
    "text.color":        TEXT_MAIN,
    "grid.color":        GRID_COL,
    "grid.linewidth":    0.6,
    "legend.facecolor":  PANEL_BG,
    "legend.edgecolor":  GRID_COL,
    "font.family":       "DejaVu Sans",
    "font.size":         10,
})


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 – Per-base sequence quality (before & after Trimmomatic)
# ─────────────────────────────────────────────────────────────────────────────

def fig_per_base_quality():
    df = pd.read_csv(os.path.join(DATA_DIR, "per_base_quality.csv"))
    pos = df["position"].values

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), facecolor=BG)
    fig.suptitle("Per-Base Sequence Quality — Before and After Adapter Trimming",
                 fontsize=13, fontweight="bold", color=TEXT_MAIN, y=1.01)

    for ax, prefix, title, accent in [
        (axes[0], "raw",  "Raw Reads  (pre-Trimmomatic)", CORAL),
        (axes[1], "trim", "Trimmed Reads  (post-Trimmomatic)", LIME),
    ]:
        mean = df[f"{prefix}_mean"].values
        q10  = df[f"{prefix}_q10"].values
        q25  = df[f"{prefix}_q25"].values
        q75  = df[f"{prefix}_q75"].values
        q90  = df[f"{prefix}_q90"].values

        # Background quality bands
        ax.axhspan(0,  20, color="#3D0000", alpha=0.35, zorder=0)
        ax.axhspan(20, 28, color="#3D2800", alpha=0.35, zorder=0)
        ax.axhspan(28, 41, color="#00270D", alpha=0.35, zorder=0)

        # Percentile bands
        ax.fill_between(pos, q10,  q90,  color=accent, alpha=0.12, label="Q10–Q90")
        ax.fill_between(pos, q25,  q75,  color=accent, alpha=0.25, label="Q25–Q75 (IQR)")
        ax.plot(pos, mean, color=accent, linewidth=1.8, label="Median quality", zorder=5)

        # Q30 threshold line
        ax.axhline(30, color=AMBER, linewidth=0.9, linestyle="--", alpha=0.8, zorder=4)
        ax.text(2, 30.6, "Q30", color=AMBER, fontsize=8, alpha=0.9)

        ax.set_xlim(1, 150);  ax.set_ylim(0, 41)
        ax.set_xlabel("Position in Read (bp)", labelpad=6)
        ax.set_ylabel("Phred Quality Score", labelpad=6)
        ax.set_title(title, fontsize=11, pad=8, color=TEXT_MAIN)
        ax.grid(True, axis="y", alpha=0.5)

        # Annotations
        pct_q30 = (mean >= 30).mean() * 100
        ax.text(0.97, 0.05, f"%Q30 ≥ {pct_q30:.0f}%",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=9, color=TEXT_DIM)

        legend = ax.legend(fontsize=8, loc="lower left")
        for text in legend.get_texts():
            text.set_color(TEXT_MAIN)

    plt.tight_layout()
    out = os.path.join(FIG_DIR, "per_base_quality.png")
    plt.savefig(out, dpi=180, bbox_inches="tight", facecolor=BG)
    plt.close()
    print(f"✓ {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 – GC content distribution
# ─────────────────────────────────────────────────────────────────────────────

def fig_gc_content():
    df = pd.read_csv(os.path.join(DATA_DIR, "gc_content.csv"))
    fig, ax = plt.subplots(figsize=(8, 4.5), facecolor=BG)

    ax.fill_between(df["gc_pct"], df["theoretical"], color=CYAN, alpha=0.18, label="Theoretical (E. coli K-12)")
    ax.plot(df["gc_pct"], df["theoretical"], color=CYAN, linewidth=1.4, linestyle="--")
    ax.fill_between(df["gc_pct"], df["actual"], color=LIME, alpha=0.25, label="Observed (EC_UTI_2023_001)")
    ax.plot(df["gc_pct"], df["actual"], color=LIME, linewidth=2.0)

    ax.axvline(50.7, color=AMBER, linewidth=0.9, linestyle=":", alpha=0.85)
    ax.text(51.5, ax.get_ylim()[1] * 0.92, "50.7% GC\n(K-12 expected)",
            color=AMBER, fontsize=8, va="top")

    ax.set_xlabel("GC Content (%)", labelpad=6)
    ax.set_ylabel("% of Reads", labelpad=6)
    ax.set_title("Per-Sequence GC Content Distribution", fontsize=12, fontweight="bold", pad=10)
    ax.set_xlim(20, 80);  ax.grid(True, alpha=0.5)
    legend = ax.legend(fontsize=9)
    for text in legend.get_texts():
        text.set_color(TEXT_MAIN)

    plt.tight_layout()
    out = os.path.join(FIG_DIR, "gc_content.png")
    plt.savefig(out, dpi=180, bbox_inches="tight", facecolor=BG)
    plt.close()
    print(f"✓ {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 – Genome-wide read depth (E. coli 4.64 Mb)
# ─────────────────────────────────────────────────────────────────────────────

def fig_genome_coverage():
    cov_df   = pd.read_csv(os.path.join(DATA_DIR, "genome_coverage.csv"))
    with open(os.path.join(DATA_DIR, "key_gene_positions.json")) as f:
        key_genes = json.load(f)

    gene_colors = {
        "gyrA": CORAL,
        "rpoB": CYAN,
        "marR": AMBER,
        "acrA": BLUE,
        "tolC": VIOLET,
    }

    fig, ax = plt.subplots(figsize=(14, 4.5), facecolor=BG)

    pos_mb = cov_df["position"].values / 1e6

    # Coverage area
    ax.fill_between(pos_mb, cov_df["coverage"], alpha=0.35, color=LIME)
    ax.plot(pos_mb, cov_df["coverage"], linewidth=0.7, color=LIME, alpha=0.8)

    # Reference lines
    mean_cov = cov_df["coverage"].mean()
    ax.axhline(mean_cov, color=CYAN, linewidth=1.1, linestyle="--", alpha=0.8, zorder=5)
    ax.text(0.01, mean_cov + 1.5, f"Mean {mean_cov:.1f}×", color=CYAN, fontsize=8)
    ax.axhline(20, color=AMBER, linewidth=0.8, linestyle=":", alpha=0.7)
    ax.text(0.01, 21.5, "20× minimum", color=AMBER, fontsize=7.5, alpha=0.8)

    # Gene markers
    for gene, info in key_genes.items():
        gpos = info["pos"] / 1e6
        col  = gene_colors.get(gene, "#FFFFFF")
        ax.axvline(gpos, color=col, linewidth=1.4, alpha=0.85, zorder=8)
        ax.text(gpos + 0.015, ax.get_ylim()[1] * 0.88 if gene not in ("marR",) else ax.get_ylim()[1] * 0.76,
                gene, color=col, fontsize=8, fontweight="bold", rotation=90, va="top")

    ax.set_xlabel("Genome Position (Mb)", labelpad=6)
    ax.set_ylabel("Read Depth (×)", labelpad=6)
    ax.set_title("Genome-Wide Read Depth — E. coli K-12 MG1655 (NC_000913.3)",
                 fontsize=12, fontweight="bold", pad=10)
    ax.set_xlim(0, 4.64);  ax.set_ylim(0, 90)
    ax.grid(True, axis="y", alpha=0.4)

    legend_elements = [Line2D([0], [0], color=gene_colors[g], linewidth=1.5, label=g)
                       for g in gene_colors]
    legend = ax.legend(handles=legend_elements, fontsize=8, loc="upper right",
                       title="Key Resistance Genes", ncol=5)
    legend.get_title().set_color(TEXT_DIM)
    for text in legend.get_texts():
        text.set_color(TEXT_MAIN)

    plt.tight_layout()
    out = os.path.join(FIG_DIR, "genome_coverage.png")
    plt.savefig(out, dpi=180, bbox_inches="tight", facecolor=BG)
    plt.close()
    print(f"✓ {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4 – Variant allele frequency distribution
# ─────────────────────────────────────────────────────────────────────────────

def fig_variant_vaf():
    df = pd.read_csv(os.path.join(RESULTS_DIR, "variants_annotated.csv"))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), facecolor=BG)
    fig.suptitle("Variant Allele Frequency (VAF) Profile",
                 fontsize=13, fontweight="bold", color=TEXT_MAIN)

    # Left: VAF histogram
    ax = axes[0]
    bins = np.linspace(0.80, 1.02, 14)
    ax.hist(df["AF"], bins=bins, color=CYAN, edgecolor=BG, linewidth=0.5, alpha=0.9, zorder=4)
    ax.axvline(0.95, color=AMBER, linewidth=1.1, linestyle="--", alpha=0.85)
    ax.text(0.955, ax.get_ylim()[1] * 0.6, "95% VAF\nthreshold", color=AMBER, fontsize=8, va="center")
    ax.set_xlabel("Variant Allele Frequency", labelpad=6)
    ax.set_ylabel("Number of Variants", labelpad=6)
    ax.set_title("VAF Distribution\n(clonal bacterial population)", fontsize=10, pad=8)
    ax.set_xlim(0.80, 1.02);  ax.grid(True, axis="y", alpha=0.5)

    # Right: scatter – VAF vs QUAL coloured by functional class
    ax2 = axes[1]
    fc_colors = {
        "NONSYNONYMOUS": CORAL,
        "SYNONYMOUS":    LIME,
        "MODIFIER":      AMBER,
        "NONSENSE":      VIOLET,
    }
    for fc, grp in df.groupby("FUNCTIONAL_CLASS"):
        col = fc_colors.get(fc, BLUE)
        ax2.scatter(grp["AF"], grp["QUAL"],
                    color=col, s=80, edgecolors=BG, linewidths=0.6,
                    label=fc.capitalize().replace("_", " "), zorder=5, alpha=0.92)
        # Label key genes
        for _, row in grp.iterrows():
            if row["GENE"] in ("gyrA", "rpoB", "marR", "ompC"):
                ax2.annotate(row["GENE"],
                             (row["AF"], row["QUAL"]),
                             xytext=(5, 3), textcoords="offset points",
                             fontsize=7.5, color=col, fontweight="bold")

    ax2.set_xlabel("Variant Allele Frequency", labelpad=6)
    ax2.set_ylabel("Variant Quality Score (QUAL)", labelpad=6)
    ax2.set_title("VAF vs. Quality Score\nColoured by Functional Class", fontsize=10, pad=8)
    ax2.set_xlim(0.82, 1.02);  ax2.grid(True, alpha=0.4)
    legend = ax2.legend(fontsize=8, loc="lower left")
    for text in legend.get_texts():
        text.set_color(TEXT_MAIN)

    plt.tight_layout()
    out = os.path.join(FIG_DIR, "variant_vaf.png")
    plt.savefig(out, dpi=180, bbox_inches="tight", facecolor=BG)
    plt.close()
    print(f"✓ {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 5 – Variant annotation breakdown
# ─────────────────────────────────────────────────────────────────────────────

def fig_variant_annotation():
    df = pd.read_csv(os.path.join(RESULTS_DIR, "variants_annotated.csv"))

    effect_map = {
        "MISSENSE_VARIANT":   "Missense",
        "SYNONYMOUS_VARIANT": "Synonymous",
        "INTERGENIC_VARIANT": "Intergenic",
        "STOP_GAINED":        "Nonsense (stop)",
        "5_PRIME_UTR_VARIANT":"5′ UTR",
    }
    df["EFFECT_LABEL"] = df["EFFECT"].map(effect_map).fillna(df["EFFECT"])
    eff_counts = df["EFFECT_LABEL"].value_counts()

    palette = [CORAL, LIME, AMBER, VIOLET, BLUE, CYAN]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), facecolor=BG)
    fig.suptitle("Variant Annotation Summary (SnpEff)",
                 fontsize=13, fontweight="bold", color=TEXT_MAIN)

    # Pie chart
    ax = axes[0]
    wedges, texts, autotexts = ax.pie(
        eff_counts.values,
        labels=eff_counts.index,
        autopct="%1.0f%%",
        startangle=140,
        colors=palette[:len(eff_counts)],
        wedgeprops=dict(edgecolor=BG, linewidth=1.8),
        textprops=dict(color=TEXT_MAIN, fontsize=9),
    )
    for at in autotexts:
        at.set_color(BG)
        at.set_fontweight("bold")
        at.set_fontsize(9)
    ax.set_title("Variant Effects Distribution", fontsize=11, pad=12)

    # Horizontal bar – count per gene (nonsynonymous only)
    ax2 = axes[1]
    ns_df = df[df["FUNCTIONAL_CLASS"].isin(["NONSYNONYMOUS", "NONSENSE"])].copy()
    gene_counts = ns_df["GENE"].value_counts()
    bars = ax2.barh(gene_counts.index, gene_counts.values,
                    color=[CORAL, VIOLET, AMBER, CYAN, BLUE][:len(gene_counts)],
                    edgecolor=BG, linewidth=0.8, alpha=0.92)

    # Annotate bars with VAF
    for bar, gene in zip(bars, gene_counts.index):
        rows = ns_df[ns_df["GENE"] == gene]
        mean_vaf = rows["AF"].mean()
        ax2.text(bar.get_width() + 0.03, bar.get_y() + bar.get_height() / 2,
                 f"VAF {mean_vaf:.1%}", va="center", fontsize=8, color=TEXT_DIM)

    ax2.set_xlabel("Number of Nonsynonymous Variants", labelpad=6)
    ax2.set_title("Nonsynonymous Variants per Gene", fontsize=11, pad=10)
    ax2.set_xlim(0, gene_counts.max() + 0.7)
    ax2.grid(True, axis="x", alpha=0.4)
    ax2.invert_yaxis()

    plt.tight_layout()
    out = os.path.join(FIG_DIR, "variant_annotation.png")
    plt.savefig(out, dpi=180, bbox_inches="tight", facecolor=BG)
    plt.close()
    print(f"✓ {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 6 – Antibiotic resistance gene panel (key clinical findings)
# ─────────────────────────────────────────────────────────────────────────────

def fig_resistance_summary():
    resist_data = [
        # Gene     Mutation    Drug Class              MIC_fold  VAF
        ("gyrA",  "D87N",     "Fluoroquinolone",       32,      0.981),
        ("gyrA",  "A90V",     "Fluoroquinolone",       16,      0.982),
        ("rpoB",  "S531L",    "Rifampicin",            256,     1.000),
        ("marR",  "G103S",    "Multidrug (AcrAB-TolC)",  8,    0.955),
        ("acrA",  "I355V",    "Multidrug (efflux)",      4,     0.966),
        ("ompC",  "K163*",    "Broad spectrum",           4,    0.923),
        ("tolC",  "L457P",    "Multidrug (efflux)",       6,    0.923),
    ]
    genes   = [r[0] for r in resist_data]
    muts    = [r[1] for r in resist_data]
    drugs   = [r[2] for r in resist_data]
    mics    = [r[3] for r in resist_data]
    vafs    = [r[4] for r in resist_data]

    drug_colors = {
        "Fluoroquinolone":           CORAL,
        "Rifampicin":                VIOLET,
        "Multidrug (AcrAB-TolC)":   AMBER,
        "Multidrug (efflux)":        CYAN,
        "Broad spectrum":            LIME,
    }
    bar_colors = [drug_colors.get(d, BLUE) for d in drugs]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), facecolor=BG)
    fig.suptitle("Antibiotic Resistance Variant Panel  —  EC_UTI_2023_001",
                 fontsize=13, fontweight="bold", color=TEXT_MAIN)

    # Left: MIC fold-change bar chart (log scale)
    ax = axes[0]
    labels = [f"{g} {m}" for g, m in zip(genes, muts)]
    bars = ax.barh(labels, np.log2(mics),
                   color=bar_colors, edgecolor=BG, linewidth=0.8, alpha=0.9)

    ax.set_xlabel("MIC Fold-Change vs Susceptible Strain (log₂)", labelpad=6)
    ax.set_title("Predicted MIC Impact by Mutation", fontsize=11, pad=8)
    ax.set_xlim(0, 10)
    ax.set_xticks([0, 2, 4, 6, 8])
    ax.set_xticklabels(["1×", "4×", "16×", "64×", "256×"])
    ax.axvline(4, color=AMBER, linewidth=0.9, linestyle="--", alpha=0.7)
    ax.text(4.05, -0.4, "Clinical\nresistance\nbreakpoint", color=AMBER,
            fontsize=7, va="top")
    ax.grid(True, axis="x", alpha=0.4)
    ax.invert_yaxis()

    for bar, val in zip(bars, mics):
        ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height() / 2,
                f"{val}×", va="center", fontsize=8.5, color=TEXT_DIM)

    # Right: VAF lollipop
    ax2 = axes[1]
    y    = np.arange(len(labels))
    ax2.hlines(y, 0.80, vafs, color=TEXT_DIM, linewidth=1.2, alpha=0.7)
    ax2.scatter(vafs, y, color=bar_colors, s=110, zorder=5, edgecolors=BG, linewidths=0.8)
    ax2.axvline(0.95, color=LIME, linewidth=0.9, linestyle="--", alpha=0.8)
    ax2.text(0.952, len(labels) - 0.3, "95%", color=LIME, fontsize=8)

    for v, yi in zip(vafs, y):
        ax2.text(v + 0.004, yi, f"{v:.1%}", va="center", fontsize=8, color=TEXT_DIM)

    ax2.set_yticks(y)
    ax2.set_yticklabels(labels, fontsize=9)
    ax2.set_xlabel("Variant Allele Frequency (VAF)", labelpad=6)
    ax2.set_title("Clonality — VAF per Resistance Variant", fontsize=11, pad=8)
    ax2.set_xlim(0.79, 1.04)
    ax2.grid(True, axis="x", alpha=0.4)
    ax2.invert_yaxis()

    # Drug class legend
    legend_patches = [mpatches.Patch(color=c, label=d) for d, c in drug_colors.items()]
    legend = ax2.legend(handles=legend_patches, fontsize=7.5, loc="lower right",
                        title="Drug Class")
    legend.get_title().set_color(TEXT_DIM)
    for text in legend.get_texts():
        text.set_color(TEXT_MAIN)

    plt.tight_layout()
    out = os.path.join(FIG_DIR, "resistance_summary.png")
    plt.savefig(out, dpi=180, bbox_inches="tight", facecolor=BG)
    plt.close()
    print(f"✓ {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("Generating figures…")
    fig_per_base_quality()
    fig_gc_content()
    fig_genome_coverage()
    fig_variant_vaf()
    fig_variant_annotation()
    fig_resistance_summary()
    print("\nAll figures written to:", FIG_DIR)
