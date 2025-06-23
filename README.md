# HCR-prober v1.5.0: Toolkit for HCR Probe Design

## Overview

HCR-prober is a powerful, flexible, and user-friendly command-line pipeline for designing DNA probes for third-generation *in situ* Hybridization Chain Reaction (HCR v3.0). It automates the tedious and critical process of selecting dozens of probe pairs that tile across a target transcript, while providing robust filtering options to ensure high specificity and optimal performance in a wide range of model and non-model organisms.

This toolkit is designed for molecular biologists and bioinformaticians who need to generate reliable, high-quality probe sets for single or multiple gene targets, including complex scenarios involving gene isoforms.

## Key Features

- **Robust Probe Design Engine:** Implements a core filtering cascade based on GC content, melting temperature (Tm), homopolymer runs, and inter-probe spacing.
- **Advanced Specificity Screening:** Leverages the power of local NCBI BLAST+ to screen probes against entire transcriptomes, ensuring they only bind to their intended target.
 - **Advanced Positive Selection Filtering:** You can now control how probes are selected based on their BLAST hits against a reference transcriptome using the `--positive-selection-strategy` flag:
   - `any-strong-hit` (Default): The original behavior. Keeps any probe with at least one high-quality BLAST hit.
   - `best-coverage`: A new, robust algorithm that identifies the transcript with the best coverage (most unique probe hits) and highest average quality (bitscore), then selects probes that are unique to that transcript. This is highly recommended for avoiding false positives from paralogs or repetitive elements.
   - `specific-id`: Keeps only probes that *uniquely* hit a user-specified transcript ID (provided with `--target-transcript-id`).
- **Automated Isoform Analysis:** Includes a powerful `isoform-split` command that automatically designs two sets of probes: one to detect all versions of a gene, and other sets to detect each unique isoform individually.
- **Batch Processing:** Built-in support for designing probes for all transcripts in a FASTA file at once, and for swapping HCR amplifiers on entire directories of existing probe sets.
- **Comprehensive & Transparent Reporting:** For every run, successful or not, the tool generates a detailed summary file with run parameters, a filtering funnel showing how many probes passed each step, and exhaustive BLAST hit tables for easy troubleshooting.
- **Extensible Amplifier System:** Add new HCR amplifiers to the system by creating a simple JSON file—no code changes required.

---

## 1. Prerequisites

The only external dependency you need to install yourself is **NCBI BLAST+**. The probe design pipeline relies on `makeblastdb` and `blastn` to perform specificity screening. This must be installed and available in your system's PATH.

Other dependencies installed with HCR-prober:         
      `  'biopython>=1.79', 'pandas>=1.3.0', 'openpyxl>=3.0.0',
        'numpy>=1.20.0', 'matplotlib>=3.3.0', 'PyYAML>=5.4.0', 'loguru>=0.5.3'`

**On Linux (Ubuntu/Debian):**
```bash
sudo apt-get update && sudo apt-get install ncbi-blast+
```

**On MacOS (using Homebrew):**
```bash
brew install blast
```

You can verify the installation by typing `blastn -version` in your terminal. All other Python dependencies will be handled automatically by the installation script.

---

## 2. Installation

You can install HCR-prober from the provided `.zip` file or by cloning a Git repository.

### Method A: From the ZIP file

1.  Unzip the project file provided.
    ```bash
    unzip HCR-prober-project-v1.4.0.zip
    ```
2.  Navigate into the newly created directory.
    ```bash
    cd HCR-prober-project-v1.4.0
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
    pip install ./
    ```

After installation, the command `hcr-prober` will be available everywhere in your system. You can verify this by running:
```bash
hcr-prober --help
```

---

## 3. The HCR-prober Pipeline: How It Works

Understanding the pipeline's logic is key to effective troubleshooting. When you run the `design` or `isoform-split` command, the following process occurs for each target sequence:

1.  **Antisense Conversion:** The input sense-strand mRNA sequence is converted to its reverse complement. All subsequent steps work on this antisense strand.
2.  **Candidate Generation (Sliding Window):** The script moves a 52bp window across the antisense sequence, generating thousands of potential probe pair candidates.
3.  **Thermo/Sequence Filtering (`prober.py`):** Each candidate is passed through a series of filters:
    - *Homopolymer Filter:* Rejects candidates with long runs of a single nucleotide (e.g., 'AAAAA').
    - *GC Content Filter:* Rejects candidates whose overall GC% is outside the range you specify with `--min-gc` and `--max-gc` (disabled by default).
    - *GC Balance Filter:* Rejects candidates where the two 25bp arms have a large difference in GC content.
    - *Melting Temperature (Tm) Filter:* Rejects candidates where either arm's Tm is outside the range you specify with `--min-tm` and `--max-tm` (disabled by default).
4.  **Overlap Filtering:** To ensure probes are tiled across the transcript, the script selects a set of non-overlapping candidates from the filtered pool.
5.  **Formatting:** The passing candidates are formatted into their final structure, assigned a unique ID (e.g., `MyGene_pair_1`), and combined with the specified HCR amplifier sequences.
6.  **Specificity Screening (`blast_wrapper.py`):**
    - If a reference transcriptome (`--blast-ref`) is provided, the target-binding region (52bp) of every passing probe is queried against the BLAST database.
    - The results are filtered by your specified cutoffs (`--min-bitscore`, `--max-evalue`).
    - Only probes that pass this rigorous check are kept.
7.  **Final Subsampling:** If more probes pass all filters than the number requested with `--max-probes`, the script selects a final, evenly-spaced subset.
8.  **Output Generation (`file_io.py`):** The final set of probes is written to the various output files, including the Excel order sheet and the detailed summary report.

---

## 4. Usage and Detailed Parameter Guide

The tool has three main commands: `design`, `isoform-split`, and `swap`.

### `hcr-prober design`
This is the primary command for designing probes against one or more transcripts.

#### **Core Arguments**
- `-i, --input`: (Required) Path to the input FASTA file containing your target transcript(s).
- `-o, --output-dir`: Directory where all output files will be saved. Defaults to `hcr_prober_output`.
- `--amplifier`: (Required) One or more HCR amplifier IDs to use (e.g., `--amplifier B1 B2`). The script will run the full pipeline for each amplifier specified.
- `--gene-name`: (Optional) If your input file contains many transcripts, use this to process only one specific transcript by its FASTA header ID.
- `--pool-name`: (Optional) A custom name for the probe pool. If not provided, a name is generated automatically.

#### **Filtering Arguments**
- `--skip-5prime`: Number of nucleotides to ignore at the 5' end of the transcript. Defaults to `100`.
- `--max-probes`: The maximum number of probe pairs to include in the final set. Defaults to `33`.
- `--min-gc`, `--max-gc`: Set a required GC content range for a probe to be kept. **Disabled by default.**
- `--min-tm`, `--max-tm`: Set a required melting temperature (Tm) range. **Disabled by default.**
- `--max-homopolymer`: The maximum allowed run of a single base (e.g., 'AAAA'). Defaults to `4`.

#### **BLAST Specificity Arguments**
- `--blast-ref`: **(Highly Recommended)** Path to your reference transcriptome in FASTA format. This is used for the positive screen to ensure probes hit a target.
- `--min-bitscore`: The minimum BLAST bitscore for a hit to be considered 'strong'. Defaults to `75`.
- `--max-evalue`: The maximum BLAST e-value for a hit to be considered 'strong'. Defaults to `1e-10`.

---

### `hcr-prober isoform-split`
This command automatically designs probes for common and unique regions of gene isoforms.

#### **Core Arguments**
- `-i, --input`: (Required) Path to a FASTA file containing **all isoforms** for the gene(s) of interest.
- `-o, --output-dir`: Directory to save the output. Defaults to `hcr_isoform_output`.
- `--gene-prefix`: (Required) The common prefix in the FASTA headers that identifies a gene group (e.g., for `MyGene_iso1` and `MyGene_iso2`, the prefix is `MyGene`).
- `--delimiter`: The character that separates the gene prefix from the isoform identifier. Defaults to `_`.
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

### Example 1: Basic Probe Design for a Sensory Gene
**Goal:** Design a set of HCR probes for the TRPA1 channel in the nudibranch *Berghia stephanieae*.

**Command:**
```bash
hcr-prober design \
    -i berghia_trpa1.fasta \
    -o berghia_probes/ \
    --amplifier B1 B2 \
    --blast-ref berghia_transcriptome.fasta \
    --min-bitscore 50
```
**Explanation:**
- We provide the single gene sequence in `berghia_trpa1.fasta`.
- We request probes using both B1 and B2 amplifiers.
- We use the full *Berghia* transcriptome to confirm our probes have strong hits.
- We lower the `--min-bitscore` to `50` as a slightly more permissive starting point for a non-traditional model.

### Example 2: Designing Probes with Strict Thermo-Filtering
**Goal:** Design probes for the *Clock* gene in the land crab *Gecarcinus lateralis*, ensuring they will work well at a specific hybridization temperature by controlling the Tm.

**Command:**
```bash
hcr-prober design \
    -i gecarcinus_clock.fasta \
    -o crab_probes/ \
    --amplifier B3 \
    --min-tm 60 \
    --max-tm 68 \
    --min-gc 45 \
    --max-gc 65 \
    --blast-ref gecarcinus_transcriptome.fasta
```
**Explanation:**
- Here, we explicitly opt-in to thermodynamic filtering.
- `--min-tm 60` and `--max-tm 68` constrain the probes to a narrow melting temperature range, increasing the likelihood they will all behave similarly in the experiment.
- `--min-gc 45` and `--max-gc 65` further ensure the probes have good sequence characteristics.

### Example 3: Automated Isoform-Specific Probe Design
**Goal:** Design two probe sets for an Opsin gene in the cave amphipod *Niphargus hrabei*. One set should detect all opsin isoforms, and another should be specific to just one of them.

**Input File (`niphargus_opsins.fa`):**
```fasta
>Nhrab_Opsin_1 ...
ATGCGTACG...
>Nhrab_Opsin_2 ...
ATGCGTTCG...
>AnotherGene_1 ...
ATGGGGCCC...
```

**Command:**
```bash
hcr-prober isoform-split \
    -i niphargus_opsins.fa \
    -o niphargus_probes/ \
    --gene-prefix Nhrab_Opsin \
    --amplifier B4 \
    --blast-ref niphargus_transcriptome.fasta
```
**Explanation:**
- The `isoform-split` command is used.
- `--gene-prefix Nhrab_Opsin` tells the script to group `Nhrab_Opsin_1` and `Nhrab_Opsin_2` together for analysis.
- The script will create a directory structure like this:
  ```
  niphargus_probes/
  └── Nhrab_Opsin/
      ├── common_probes/      (Probes that hit both isoforms)
      └── isoform_specific_probes/ (Probes for each unique part)
  ```
