# HCR-prober v1.9.6

This is a professional-grade bioinformatics platform for designing and troubleshooting Hybridization Chain Reaction (HCR) v3.0 probe sets. This version fixes a critical bug and improves the clarity of the output visualizations.

## Key Features in This Version

- **NEW: Highlighted Probes in BLAST Report:** The detailed BLAST report in the output summary file now includes a `SELECTED` column. Hits that correspond to probes in the final, filtered set are marked with an asterisk (`*`), making it much easier to trace which candidates passed all filtering steps.
- **IMPROVED: Clearer Probe Map Visualization:** The output SVG probe map renders all probes on a single side of the transcript for an intuitive view of tiling and density.
- **Optimal Probe Spacing Algorithm:** The pipeline uses a dynamic programming algorithm to select the **absolute maximum number of non-overlapping probes**, ensuring the densest possible coverage.
- **Hybrid BLAST Strategy for Isoform Analysis:** The `isoform-split` command uses smart, context-aware defaults (`--common-strategy any-strong-hit` and `--unique-strategy best-coverage`) to maximize probe counts for shared regions while ensuring specificity for unique regions.

---

## 1. Prerequisites

The only external dependency you need to install yourself is **NCBI BLAST+**. This must be installed and available in your system's PATH.

**On Linux (Ubuntu/Debian):**
```bash
sudo apt-get update && sudo apt-get install ncbi-blast+
```

**On MacOS (using Homebrew):**
```bash
brew install blast
```

---

## 2. Installation

1.  Unzip the project file provided (`HCR-prober-project-v1.9.5.zip`).
2.  Navigate into the newly created directory.
    ```bash
    cd HCR-prober-project-v1.9.5
    ```
3.  Install the package using `pip`.
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


---

## 3. The HCR-prober Pipeline: How It Works

1.  **Candidate Generation:** The script generates thousands of potential probe candidates.
2.  **Thermo/Sequence Filtering:** Candidates are filtered for GC content, Tm, homopolymers, etc.
3.  **Specificity Screening (BLAST):** Candidates are BLASTed against a reference database.
4.  **Optimal Spacing & Selection:** A dynamic programming algorithm selects the maximum possible number of non-overlapping probes from the specific candidates.
5.  **Output Generation:** The final, dense probe set is written to user-friendly output files, including the improved SVG probe map.

---

## 4. Usage and Detailed Parameter Guide

### `hcr-prober design`
The primary command for designing probes.

#### **Core Arguments**
- `-i, --input`: (Required) Path to the input FASTA file.
- `-o, --output-dir`: Directory for output files.
- `--amplifier`: (Required) One or more HCR amplifier IDs.
- `--min-probe-distance`: Minimum distance (in nt) between adjacent probes. Defaults to `2`.

### `hcr-prober isoform-split`
Designs probes for common and unique regions of gene isoforms.

#### **Hybrid Strategy Arguments**
- `--common-strategy`: BLAST strategy for the common probe set. Defaults to `any-strong-hit`.
- `--unique-strategy`: BLAST strategy for the isoform-specific sets. Defaults to `best-coverage`.

---

## 5. Use Cases and Examples

### Example 1: Maximizing Probe Count
**Goal:** Get the densest possible probe set for a gene, `MyGene`.
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
By setting `--min-probe-distance 0` and relaxing GC constraints, you allow the optimal spacing algorithm to pack the maximum number of probes onto the transcript.

### Example 2: Optimal, One-Shot Isoform Probe Design
**Goal:** Design probes for an Opsin gene with two isoforms.
```bash
hcr-prober isoform-split \
    -i all_opsin_isoforms.fasta \
    -o opsin_probes/ \
    --gene-prefix Opsin \
    --amplifier B3 \
    --blast-ref transcriptome.fasta
```
**Explanation:**
The script automatically uses the best BLAST strategy for each phase: `any-strong-hit` for common regions and `best-coverage` for unique regions, giving you the best possible result in a single command.
