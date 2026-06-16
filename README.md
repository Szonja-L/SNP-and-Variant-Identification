# 🧬 SNP & Variant Identification in *E. coli* — NGS Pipeline Portfolio Project

[![Python](https://img.shields.io/badge/Python-3.9+-3776AB?logo=python&logoColor=white)](https://python.org)
[![Galaxy](https://img.shields.io/badge/Galaxy-usegalaxy.eu-2D9BF0?logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAxMDAgMTAwIj48Y2lyY2xlIGN4PSI1MCIgY3k9IjUwIiByPSI0NSIgZmlsbD0id2hpdGUiLz48L3N2Zz4=)](https://usegalaxy.eu)
[![FreeBayes](https://img.shields.io/badge/FreeBayes-1.3.6-FF6B6B)](https://github.com/freebayes/freebayes)
[![BWA](https://img.shields.io/badge/BWA--MEM-0.7.17-4ECDC4)](https://github.com/lh3/bwa)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Portfolio project demonstrating:** Next-Generation Sequencing (NGS) quality control, reference genome alignment, Bayesian variant calling, functional annotation, and clinical interpretation of antibiotic resistance mutations.

---

## 🔬 Biological Scenario

A *Escherichia coli* clinical isolate (UTI patient, Netherlands, 2023) was whole-genome sequenced on an **Illumina NovaSeq 6000** (150 bp paired-end). The reads were aligned to the *E. coli* K-12 MG1655 reference genome (NC_000913.3) and variants were called to identify mutations that could explain the isolate's **multidrug-resistant (MDR)** phenotype.

**Key questions:**
- What SNPs distinguish this clinical strain from the K-12 laboratory reference?
- Which mutations are in known resistance genes?
- What drug classes are predicted to be affected?

---

## 🚀 Quick Start

```bash
# Clone the repository
git clone https://github.com/szfekete/snp-variant-identification.git
cd snp-variant-identification

# Install dependencies
pip install numpy pandas matplotlib seaborn scipy

# Generate synthetic demonstration data
python scripts/generate_synthetic_data.py

# Generate all figures
python scripts/generate_figures.py

# Run the full analysis pipeline
python scripts/variant_analysis_pipeline.py

# Optional: regenerate figures alongside analysis
python scripts/variant_analysis_pipeline.py --figures
```

---

## 📂 Repository Structure

```
snp-variant-identification/
│
├── 📓 notebooks/
│   └── snp_analysis_pipeline.ipynb     # Full Jupyter narrative
│
├── 🐍 scripts/
│   ├── generate_synthetic_data.py      # Creates realistic NGS data
│   ├── generate_figures.py             # All publication-quality figures
│   └── variant_analysis_pipeline.py    # Main analysis pipeline
│
├── 📊 data/
│   ├── per_base_quality.csv            # FastQC-equivalent QC metrics
│   ├── gc_content.csv                  # Per-sequence GC distribution
│   ├── genome_coverage.csv             # Per-5kb read depth
│   ├── alignment_stats.json            # BWA-MEM + SAMtools flagstat
│   ├── qc_summary.json                 # Trimmomatic QC summary
│   └── EC_UTI_2023_001_R1_example.fastq # Sample FASTQ (format demo)
│
├── 📈 results/
│   ├── variants.vcf                    # FreeBayes/SnpEff annotated VCF
│   ├── variants_annotated.csv          # Full variant table
│   ├── resistance_variants.csv         # Resistance-associated variants
│   └── pipeline_summary.json           # Computed summary statistics
│
├── 🖼️  figures/
│   ├── per_base_quality.png            # FastQC before/after Trimmomatic
│   ├── gc_content.png                  # GC content distribution
│   ├── genome_coverage.png             # Genome-wide read depth
│   ├── variant_vaf.png                 # VAF distribution + VAF vs QUAL
│   ├── variant_annotation.png          # SnpEff functional breakdown
│   └── resistance_summary.png          # Resistance gene panel
│
├── 🌌 galaxy_workflow/
│   └── ecoli_variant_calling.ga        # Importable Galaxy workflow
│
└── 📄 report/
    └── snp_portfolio_report.html       # Interactive HTML report
```

---

## 🔧 Galaxy Workflow

The complete pipeline is available as an **importable Galaxy workflow** (`galaxy_workflow/ecoli_variant_calling.ga`).

To use it:
1. Go to [usegalaxy.eu](https://usegalaxy.eu) and create a free account
2. Navigate to **Workflow → Import** and upload the `.ga` file
3. Upload your FASTQ files and reference genome, then run

### Pipeline Steps

| # | Tool | Version | Purpose |
|---|------|---------|---------|
| 1 | **FastQC** | 0.12.1 | Raw read quality assessment |
| 2 | **Trimmomatic** | 0.39 | Adapter removal & quality trimming (PE) |
| 2b | **FastQC** | 0.12.1 | Post-trim QC verification |
| 3 | **BWA Index** | 0.7.17 | Index reference genome |
| 4 | **BWA-MEM** | 0.7.17 | Paired-end alignment to reference |
| 5 | **SAMtools Sort** | 1.18 | Coordinate sort BAM |
| 5b | **SAMtools Markdup** | 1.18 | Mark PCR duplicates |
| 5c | **SAMtools Flagstat** | 1.18 | Alignment statistics |
| 5d | **BEDtools genomeCoverage** | 2.31 | Per-base read depth |
| 6 | **FreeBayes** | 1.3.6 | Bayesian variant calling (ploidy=1) |
| 7 | **VCFfilter** | 1.0 | Quality filter (QUAL≥200, AF≥0.90) |
| 8 | **SnpEff** | 5.1d | Functional annotation |

---

## 📊 Results Overview

### Quality Control

| Metric | Raw | After Trimmomatic |
|--------|-----|-------------------|
| Total reads | 3,142,088 | 3,056,442 |
| % Reads surviving | — | 97.3% |
| Mean Phred quality | 35.8 | 37.2 |
| **% Q30** | **87.4%** | **93.1%** |
| GC content | 50.7% | 50.7% |

### Alignment Statistics (BWA-MEM → NC_000913.3)

| Metric | Value |
|--------|-------|
| Mapping rate | **98.57%** |
| Properly paired | 97.21% |
| Mean coverage depth | **48.7×** |
| Genome covered ≥ 20× | 98.94% |
| Mean insert size | 372 bp |

### Variant Summary

| Category | Count |
|----------|-------|
| Total high-quality variants | **11** |
| Nonsynonymous (missense) | 7 |
| Nonsense (stop-gained) | 1 |
| Synonymous | 2 |
| Modifier (intergenic/UTR) | 1 |
| **Resistance variants** | **7** |
| Genes with resistance mutations | 6 |
| Drug classes affected | 5 |

---

## 🛡️ Key Resistance Findings

| Gene | Mutation | Drug Class | Mechanism | MIC Fold-↑ | VAF |
|------|----------|-----------|-----------|-----------|-----|
| `gyrA` | **D87N** | Fluoroquinolone | Altered target — QRDR | 32× | 98.1% |
| `gyrA` | **A90V** | Fluoroquinolone | Altered target — QRDR (additive) | 16× | 98.2% |
| `rpoB` | **S531L** | Rifampicin | Altered RNA polymerase β | 256× | 100.0% |
| `marR` | **G103S** | Multidrug (MarA) | MDR regulator → AcrAB-TolC efflux | 8× | 95.5% |
| `acrA` | **I355V** | Multidrug efflux | Efflux pump MFP structural change | 4× | 96.6% |
| `tolC` | **L457P** | Multidrug efflux | Efflux channel altered specificity | 6× | 92.3% |
| `ompC` | **K163\*** | Broad spectrum | Porin truncation → reduced permeability | 4× | 93.6% |

> ⚠️ **MDR alert:** The double `gyrA` mutation (D87N + A90V) at positions 87 and 90 of the quinolone-resistance-determining region (QRDR) predicts **high-level fluoroquinolone resistance**. Combined with marR-mediated efflux upregulation and ompC porin loss, this isolate displays a multi-mechanism MDR genotype.

---

## 🖼️ Figures

<details>
<summary>Click to expand figure descriptions</summary>

**Figure 1 — Per-Base Sequence Quality**
FastQC-style before/after plot showing the characteristic Illumina 3′ quality drop (raw) corrected by Trimmomatic (trimmed). %Q30 improves from 87.4% → 93.1%.

**Figure 2 — GC Content Distribution**
Per-sequence GC content closely matches the theoretical *E. coli* K-12 distribution (50.7% GC). No contamination peaks detected.

**Figure 3 — Genome-Wide Coverage Depth**
Highly uniform 48.7× coverage across the 4.64 Mb chromosome. The rRNA operon trough (~4.17 Mb) is expected due to repetitive sequences. Resistance gene positions marked with vertical lines.

**Figure 4 — Variant Allele Frequency (VAF)**
All 11 variants cluster at VAF > 88%, confirming a clonal bacterial population. High-confidence resistance variants show VAF > 95% with QUAL scores of 700–1000.

**Figure 5 — Variant Annotation (SnpEff)**
Breakdown of variant effects: predominantly missense (63.6%), one nonsense mutation in ompC, two synonymous changes. `gyrA` carries 2 nonsynonymous variants — both in the QRDR.

**Figure 6 — Resistance Gene Panel**
MIC fold-change predictions and VAF lollipop for all 7 resistance variants. `rpoB` S531L shows the highest individual MIC impact (256×). All variants are at high allele frequency.
</details>

---

## 🧪 Using a Real VCF from Galaxy

Once you've run the Galaxy workflow on your own data, you can plug the SnpEff-annotated VCF into the analysis script:

```bash
python scripts/variant_analysis_pipeline.py \
  --vcf results/my_snpeff_annotated.vcf \
  --min-qual 200 \
  --min-dp 10 \
  --min-af 0.90 \
  --figures
```

---

## 📚 Key Concepts Demonstrated

| Concept | Where demonstrated |
|---------|-------------------|
| Illumina sequencing quality metrics | FastQC plots, Phred score explanation |
| Adapter trimming rationale | Trimmomatic parameter justification |
| Seed-and-extend alignment | BWA-MEM section of notebook |
| PCR duplicate marking | SAMtools markdup step |
| Haplotype-based variant calling | FreeBayes ploidy=1 explanation |
| VCF format & INFO fields | VCF inspection cell in notebook |
| Functional effect prediction | SnpEff ANN annotation parsing |
| Clinical resistance interpretation | Resistance profile section |

---

## 🔗 References & Databases

- **NCBI Reference:** [NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3) — *E. coli* K-12 MG1655
- **CARD Database:** McArthur et al. (2023) *Nucleic Acids Research*
- **FreeBayes:** Garrison & Marth (2012) *arXiv:1207.3907*
- **SnpEff:** Cingolani et al. (2012) *Fly* 6(2):80–92
- **gyrA QRDR mutations:** Jacoby (2005) *Clin Infect Dis* 41:S120–S126
- **rpoB S531L:** Goldstein (2014) *Clin Infect Dis* 59:S93–S106
- **MarA/MarR circuit:** Martin & Rosner (2002) *Mol Microbiol* 44:1–7
- **Galaxy platform:** Afgan et al. (2022) *Nucleic Acids Research* 50:W345–W351

---

## 📋 Skills Demonstrated

- **Bioinformatics platforms:** Galaxy (workflow design, tool parameterisation, history management)
- **NGS QC:** FastQC interpretation, Trimmomatic PE trimming, quality metric benchmarks
- **Alignment:** BWA-MEM, read group tagging, SAMtools sort/markdup/flagstat, coverage analysis
- **Variant calling:** FreeBayes haploid mode, VCF format, quality filtering strategies
- **Annotation:** SnpEff functional classification, HGVS notation, effect prediction
- **Resistance genomics:** CARD/ResFinder cross-referencing, QRDR mutation biology, efflux mechanisms
- **Python data science:** pandas, matplotlib/seaborn, numpy, argparse CLI design

---

*Szonja Fekete · Bioinformatics Portfolio · 2024*
