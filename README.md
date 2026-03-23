
# HCR-prober v1.9.6: Toolkit for HCR Probe Design

## Overview

HCR-prober is a powerful, flexible, and user-friendly command-line pipeline for designing DNA probes for third-generation *in situ* Hybridization Chain Reaction (HCR v3.0). It is a professional-grade bioinformatics platform that automates the tedious and critical process of selecting dozens of probe pairs that tile across a target transcript. It provides robust filtering options to ensure high specificity and optimal performance in a wide range of model and non-model organisms.

This toolkit is designed for molecular biologists and bioinformaticians who need to generate reliable, high-quality probe sets for single or multiple gene targets, including complex scenarios involving gene isoforms. This version fixes a critical bug and improves the clarity of the output visualizations.

## Key Features

- **Robust Probe Design Engine:** Implements a core filtering cascade based on GC content, melting temperature (Tm), homopolymer runs, and inter-probe spacing.
- **Optimal Probe Spacing Algorithm:** Uses a dynamic programming algorithm to select the **absolute maximum number of non-overlapping probes**, ensuring the densest possible coverage across the transcript.
- **Advanced Specificity Screening:** Leverages the power of local NCBI BLAST+ to screen probes against entire transcriptomes, ensuring they only bind to their intended target.
- **Advanced Positive Selection Filtering:** Control how probes are selected based on their BLAST hits against a reference transcriptome using the `--positive-selection-strategy` flag:
  - `any-strong-hit` (Default): Keeps any probe with at least one high-quality BLAST hit.
  - `best-coverage`: A robust algorithm that identifies the transcript with the best coverage (most unique probe hits) and highest average quality (bitscore), then selects probes unique to that transcript. Highly recommended for avoiding false positives from paralogs.
  - `specific-id`: Keeps only probes that *uniquely* hit a user-specified transcript ID.
- **Automated Isoform Analysis with Hybrid BLAST Strategy:** A powerful `isoform-split` command that automatically designs probes for common and unique regions of gene isoforms. It uses smart, context-aware defaults (`--common-strategy any-strong-hit` and `--unique-strategy best-coverage`) to maximize probe counts for shared regions while ensuring specificity for unique regions.
- **Batch Processing:** Built-in support for designing probes for all transcripts in a FASTA file at once, and for swapping HCR amplifiers on entire directories of existing probe sets.
- **Comprehensive & Transparent Reporting:** For every run, the tool generates a detailed summary file with run parameters and a filtering funnel.
- **NEW: Highlighted Probes in BLAST Report:** The detailed BLAST report now includes a `SELECTED` column, marking probes in the final set with an asterisk (`*`) for easy traceability.
- **IMPROVED: Clearer Probe Map Visualization:** The output SVG probe map renders all probes on a single side of the transcript for an intuitive view of tiling and density.
- **Extensible Amplifier System:** Add new HCR amplifiers to the system by creating a simple JSON file—no code changes required.

---

## 1. Prerequisites

The only external dependency you need to install yourself is **NCBI BLAST+**. The probe design pipeline relies on `makeblastdb` and `blastn` to perform specificity screening. This must be installed and available in your system's PATH.

Other Python dependencies installed automatically by HCR-prober:
```
'biopython>=1.79', 'pandas>=1.3.0', 'openpyxl>=3.0.0',
'numpy>=1.20.0', 'matplotlib>=3.3.0', 'PyYAML>=5.4.0', 'loguru>=0.5.3'
```

**On Linux (Ubuntu/Debian):**
```bash
sudo apt-get update && sudo apt-get install ncbi-blast+
```

**On MacOS (using Homebrew):**
```bash
brew install blast
```

You can verify the installation by typing `blastn -version` in your terminal.

---

## 2. Installation

You can install HCR-prober from the provided `.zip` file or by cloning the Git repository.

### Method A: From a ZIP file

1.  Unzip the project file (`HCR-prober-project-v1.9.6.zip`).
2.  Navigate into the newly created directory.
    ```bash
    cd HCR-prober-project-v1.9.6
    ```
3.  Install the package using `pip`. The `.` tells pip to install from the current directory.
    ```bash
    pip install .
    ```

### Method B: From a Git Repository

1.  Clone the repository to your machine.
    ```bash
    git clone https://github.com/invertome/hcr-prober.git
    ```
2.  Navigate into the repository directory.
    ```bash
    cd hcr-prober
    ```
3.  Install the package.
    ```bash
    pip install .
    ```

After installation, the command `hcr-prober` will be available everywhere in your system. You can verify this by running:
```bash
hcr-prober --help
```

---

## 3. The HCR-prober Pipeline: How It Works

Understanding the pipeline's logic is key to effective troubleshooting. When you run the `design` or `isoform-split` command, the following process occurs:

1.  **Antisense Conversion:** The input sense-strand mRNA sequence is converted to its reverse complement. All subsequent steps work on this antisense strand.
2.  **Candidate Generation (Sliding Window):** The script moves a 52bp window across the antisense sequence, generating thousands of potential probe pair candidates.
3.  **Thermo/Sequence Filtering (`prober.py`):** Each candidate is passed through a series of filters:
    - *Homopolymer Filter:* Rejects candidates with long runs of a single nucleotide.
    - *GC Content Filter:* Rejects candidates whose overall GC% is outside the specified range.
    - *GC Balance Filter:* Rejects candidates where the two 25bp arms have a large difference in GC content.
    - *Melting Temperature (Tm) Filter:* Rejects candidates where either arm's Tm is outside the specified range.
4.  **Specificity Screening (`blast_wrapper.py`):**
    - The target-binding region (52bp) of every passing probe is queried against a BLAST database built from your reference transcriptome (`--blast-ref`).
    - Results are filtered by your specified cutoffs (`--min-bitscore`, `--max-evalue`).
5.  **Optimal Spacing & Selection:** A dynamic programming algorithm selects the **maximum possible number of non-overlapping probes** from the specific candidates, honoring the `--min-probe-distance` parameter.
6.  **Formatting:** Passing candidates are assigned a unique ID (e.g., `MyGene_pair_1`) and combined with the specified HCR amplifier sequences.
7.  **Final Subsampling:** If more probes pass all filters than requested with `--max-probes`, the script selects a final, evenly-spaced subset.
8.  **Output Generation (`file_io.py`):** The final probe set is written to user-friendly output files, including an Excel order sheet, the improved SVG probe map, and a detailed summary report.

---

## 4. Usage and Detailed Parameter Guide

The tool has three main commands: `design`, `isoform-split`, and `swap`.

### `hcr-prober design`
This is the primary command for designing probes against one or more transcripts.

#### **Core Arguments**
- `-i, --input`: (Required) Path to the input FASTA file containing your target transcript(s).
- `-o, --output-dir`: Directory where all output files will be saved. Defaults to `hcr_prober_output`.
- `--amplifier`: (Required) One or more HCR amplifier IDs to use (e.g., `--amplifier B1 B2`).
- `--gene-name`: (Optional) If your input file contains many transcripts, use this to process only one specific transcript by its FASTA header ID.
- `--pool-name`: (Optional) A custom name for the probe pool.

#### **Filtering Arguments**
- `--skip-5prime`: Number of nucleotides to ignore at the 5' end of the transcript. Defaults to `100`.
- `--max-probes`: The maximum number of probe pairs to include in the final set. Defaults to `33`.
- `--min-probe-distance`: Minimum distance (in nt) between adjacent probes. Defaults to `2`.
- `--min-gc`, `--max-gc`: Set a required GC content range for a probe to be kept. **Disabled by default.**
- `--min-tm`, `--max-tm`: Set a required melting temperature (Tm) range. **Disabled by default.**
- `--max-homopolymer`: The maximum allowed run of a single base. Defaults to `4`.

#### **BLAST Specificity Arguments**
- `--blast-ref`: **(Highly Recommended)** Path to your reference transcriptome in FASTA format.
- `--min-bitscore`: The minimum BLAST bitscore for a hit to be considered 'strong'. Defaults to `75`.
- `--max-evalue`: The maximum BLAST e-value for a hit to be considered 'strong'. Defaults to `1e-10`.

---

### `hcr-prober isoform-split`
This command automatically designs probes for common and unique regions of gene isoforms.

#### **Core Arguments**
- `-i, --input`: (Required) Path to a FASTA file containing **all isoforms** for the gene(s) of interest.
- `-o, --output-dir`: Directory to save the output. Defaults to `hcr_isoform_output`.
- `--gene-prefix`: (Required) The common prefix in the FASTA headers that identifies a gene group (e.g., for `MyGene_iso1`, the prefix is `MyGene`).
- `--delimiter`: The character that separates the gene prefix from the isoform identifier. Defaults to `_`.

#### **Hybrid Strategy Arguments**
- `--common-strategy`: BLAST strategy for the common probe set. Defaults to `any-strong-hit`.
- `--unique-strategy`: BLAST strategy for the isoform-specific sets. Defaults to `best-coverage`.

*This command uses all the same Filtering and BLAST arguments as the `design` command.*

---

### `hcr-prober swap`
This utility swaps the HCR amplifiers on an existing probe set order form.

#### **Core Arguments**
- `--input-probes`: (Required) Path to a single `.xlsx` order file or an entire directory containing them.
- `--output-dir`: Directory to save the newly swapped probe files. Defaults to `swapped_probes`.
- `--new-amplifier`: (Required) The ID of the new amplifier to swap to (e.g., `B3`).

---

## 5. Use Cases and Examples

### Example 1: Maximizing Probe Count
**Goal:** Get the densest possible probe set for a gene, `MyGene`, by allowing probes to be immediately adjacent.
```bash
hcr-prober design \
    -i MyGene.fasta \
    -o MyGene_probes/ \
    --amplifier B1 \
    --blast-ref transcriptome.fasta \
    --min-gc 35 --max-gc 65 \
    --min-probe-distance 0
```
**Explanation:**
By setting `--min-probe-distance 0` and relaxing GC constraints, you allow the optimal spacing algorithm to pack the maximum number of probes onto the transcript. A value of 0 is not recommended.

### Example 2: Optimal, One-Shot Isoform Probe Design
**Goal:** Design probe sets for an Opsin gene with two isoforms: one set that detects all isoforms and specific sets for each unique isoform. Sequence IDs are in the format of Opsin1_iso1, Opsin1_iso2, Opsin2_iso1, etc.
```bash
hcr-prober isoform-split \
    -i all_opsin_isoforms.fasta \
    -o opsin_probes/ \
    --gene-prefix Opsin \
    --amplifier B3 \
    --blast-ref transcriptome.fasta
```
**Explanation:**
The script automatically uses the best BLAST strategy for each phase: `any-strong-hit` for common regions (to maximize shared probes) and `best-coverage` for unique regions (to ensure high specificity). This gives you the best possible result in a single command.

### Example 3: Designing Probes with Strict Thermo-Filtering
**Goal:** Design probes for the *Clock* gene in the land crab *Gecarcinus lateralis*, ensuring they will work well at a specific hybridization temperature by controlling the Tm.
```bash
hcr-prober design \
    -i Gecarcinus_clock.fasta \
    -o crab_probes/ \
    --amplifier B4 \
    --min-tm 60 --max-tm 68 \
    --min-gc 45 --max-gc 65 \
    --blast-ref Gecarcinus_transcriptome.fasta
```
**Explanation:**
Here, we explicitly opt-in to thermodynamic filtering. Constraining the melting temperature (`--min-tm`, `--max-tm`) increases the likelihood that all probes will behave similarly in the experiment.
