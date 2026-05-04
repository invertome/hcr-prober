# HCR-prober v1.12.0

A command-line pipeline for designing DNA probes for third-generation *in situ* Hybridization Chain Reaction (HCR v3.0). Given an mRNA target, hcr-prober produces an order-ready set of split-initiator probe pairs that tile the transcript, with optional BLAST-based specificity screening, secondary-structure filtering, and a fully reproducible audit trail.

Designed for molecular biologists working with model and non-model organisms, including complex isoform families.

---

## Table of contents

1. [Quickstart](#quickstart)
2. [Installation](#installation)
3. [How the pipeline works](#how-the-pipeline-works)
4. [Commands](#commands)
5. [Hybridisation buffer conditions](#hybridisation-buffer-conditions)
6. [Output files](#output-files)
7. [Examples](#examples)
8. [Reproducibility](#reproducibility)
9. [Amplifier sources](#amplifier-sources)
10. [What's new in v1.12.0](#whats-new-in-v1120)

---

## Quickstart

```bash
# install
pip install -e .

# design probes for one mRNA, screening against the host transcriptome
hcr-prober design \
    -i my_gene.fasta \
    --amplifier B1 \
    --blast-ref host_transcriptome.fasta \
    -o my_gene_probes/

# design isoform-aware probes (common-region + per-isoform unique sets)
hcr-prober isoform-split \
    -i all_isoforms.fasta \
    --gene-prefix MyGene \
    --amplifier B3 \
    --blast-ref host_transcriptome.fasta

# swap amplifier on an existing order sheet
hcr-prober swap \
    --input-probes my_gene_probes/MyGene/B1/MyGene_B1_order.xlsx \
    --new-amplifier B5
```

Defaults reflect HCR hybridisation conditions (5×SSC + 50 % formamide at 37 °C). For PCR conditions, pass `--buffer-preset pcr`.

---

## Installation

### Prerequisites

- Python ≥ 3.8
- **NCBI BLAST+ ≥ 2.10** (`blastn`, `makeblastdb`) on your PATH

```bash
# Linux (Debian/Ubuntu)
sudo apt-get install ncbi-blast+

# macOS
brew install blast

# verify
blastn -version
```

### Install hcr-prober

```bash
git clone https://github.com/invertome/hcr-prober.git
cd hcr-prober
pip install -e .
hcr-prober --help
```

Python dependencies (installed automatically): `biopython>=1.80`, `pandas>=1.3.0`, `openpyxl>=3.0.0`, `numpy>=1.20.0`, `matplotlib>=3.3.0`, `PyYAML>=5.4.0`, `loguru>=0.5.3`, `primer3-py>=2.0.0`.

---

## How the pipeline works

For each input transcript, the `design` command runs the following filtering funnel. The exact counts after each step are reported in the run summary so you can see where candidates are lost.

| # | Step | What it does |
|---|------|--------------|
| 1 | Antisense conversion | Reverse-complement the input mRNA. All windows are sampled from the antisense strand so each probe arm is the complementary sequence to the mRNA target. |
| 2 | Sliding window | Generate all 52-nt windows (= 25-nt arm + 2-nt gap + 25-nt arm). |
| 3 | 5'-skip | Discard windows in the first `--skip-5prime` nt of the transcript (default 50). |
| 4 | ACGT-only filter | Reject any window containing N or other IUPAC ambiguity codes (R, Y, S, W, K, M…). |
| 5 | Region/sequence masking | Optional: mask user-specified intervals (`--mask-regions`) or sequences (`--mask-sequences`). |
| 6 | Homopolymer filter | Reject windows containing runs of ≥ `--max-homopolymer + 1` of the same base (default rejects runs ≥ 5). |
| 7 | Per-arm GC filter | Both 25-nt arms must independently fall inside `[--min-gc, --max-gc]` AND differ by no more than `--max-gc-diff`. |
| 8 | Tm filter | Both arms must have a (formamide-corrected) Tm inside `[--min-tm, --max-tm]`. |
| 9 | Secondary-structure filter | primer3-py rejects arms that form stable hairpins or homodimers, or pairs that form stable heterodimers, at the configured ionic conditions and 37 °C. |
| 10 | Specificity (BLAST) | Each candidate's joined 25+NN+25 sequence is BLASTed against `--blast-ref`. The 2-nt `NN` gap is intentional: HCR signal requires *both* arms to bind in close proximity, so a query that scores as a 52-mer enforces co-localisation. The query is searched on the reverse-complement strand (`-strand minus`) since probes are antisense to mRNA. |
| 11 | Selection strategy | One of `any-strong-hit`, `best-coverage`, `specific-id` (see below). |
| 12 | Negative BLAST | Optional: reject probes with strong off-target hits. The negative reference is auto-derived from `--blast-ref` minus the target gene and any IDs supplied via `--target-transcript-id`. |
| 13 | Optimal spacing | A weighted-interval-scheduling DP picks the **maximum** non-overlapping subset honouring `--min-probe-distance`. |
| 14 | Optional Tm uniformity | If `--max-tm-sigma N` is set, iteratively drop the most extreme-Tm probe until σ(Tm) ≤ N. |
| 15 | Subsampling | If more than `--max-probes` survive, pick a quasi-uniform subset. |
| 16 | Amplifier finalisation | Per amplifier, append the up/dn split-initiator handles + spacers; produce one `_pair_N` ID per pair. |

---

## Commands

The CLI exposes three subcommands: `design`, `isoform-split`, and `swap`.

### `hcr-prober design`

Design probes against one or more transcripts.

#### Core
| Flag | Default | Notes |
|------|---------|-------|
| `-i / --input` | — (required) | Input FASTA; one or more transcripts. |
| `-o / --output-dir` | `hcr_prober_output` | All outputs are written here. |
| `--amplifier` | — (required) | One or more HCR amplifier IDs (e.g. `--amplifier B1 B2 B3`). |
| `--gene-name` | — | Process a single transcript ID from the input. |
| `--pool-name` | — | Override the default pool name in the order sheet. |

#### Process control
| Flag | Default | Notes |
|------|---------|-------|
| `--seed N` | `0` | RNG seed; output is deterministic given the seed. |
| `--threads N` | `1` | Pass `-num_threads N` to `blastn`. |
| `--dry-run` | off | Run filters and produce the audit funnel; skip BLAST entirely. |
| `--verbose` / `--quiet` | off | Set log level to DEBUG / WARNING. |
| `--db-path DIR` | tempdir | Persistent location for BLAST databases (cached by content-hashed key, so two refs with the same basename do not collide). |
| `--force` | off | Ignore cached BLAST DBs. |

#### Sliding-window / spacing
| Flag | Default | Notes |
|------|---------|-------|
| `--window-size N` | `52` | Total probe-pair footprint. Almost never changed. |
| `--probe-len N` | `25` | Length of each arm. |
| `--spacer-len N` | `2` | Gap between the two arms in the joined query. |
| `--skip-5prime N` | `50` | nt to exclude from the 5′ end of the transcript. |
| `--max-probes N` | `33` | Cap on the final probe count. |
| `--min-probe-distance N` | `2` | Minimum gap between adjacent probe footprints. |
| `--mask-regions S` | — | Comma-separated 1-based intervals to exclude (e.g. `'1-100,500-650'`). |
| `--mask-sequences PATH` | — | FASTA of sequences to exclude (e.g. repeats). |

#### Sequence / thermo filters
| Flag | Default | Notes |
|------|---------|-------|
| `--min-gc` / `--max-gc` | `40` / `60` (%) | **Per-arm** GC range. |
| `--max-gc-diff` | `15` (%) | Maximum GC % difference between the two arms. |
| `--min-tm` / `--max-tm` | `40` / `55` (°C) | Per-arm Tm range, **after formamide correction**. |
| `--max-homopolymer N` | `4` | Reject windows with runs ≥ N+1 of any base. |
| `--max-tm-sigma N` | off | Optional: enforce σ(Tm) ≤ N across the final pool. |

#### Hybridisation buffer
| Flag | Default | Notes |
|------|---------|-------|
| `--buffer-preset {hcr-5xssc, pcr}` | `hcr-5xssc` | Convenience: sets Na+, Mg²⁺, formamide together. |
| `--na-conc` | preset (825 / 50 mM) | Sodium concentration, mM. 5×SSC ≈ 825 mM Na⁺. |
| `--mg-conc` | preset (0 mM) | Magnesium concentration, mM. |
| `--dntp-conc` | `0` mM | dNTP concentration. |
| `--dna-conc` | `25` nM | Probe oligo concentration. |
| `--formamide-pct` | preset (50 / 0 %) | Formamide percentage. Tm is reduced by 0.65 °C per percent. |

#### Secondary structure
| Flag | Default | Notes |
|------|---------|-------|
| `--max-hairpin-dg` | `-3.0` (kcal/mol) | Reject arm if hairpin ΔG is more negative. |
| `--max-homodimer-dg` | `-5.0` (kcal/mol) | Reject arm if homodimer ΔG is more negative. |
| `--max-heterodimer-dg` | `-5.0` (kcal/mol) | Reject pair if cross-arm heterodimer ΔG is more negative. |

All three values are evaluated under the same Na⁺/Mg²⁺/formamide conditions as the Tm filter.

#### BLAST specificity
| Flag | Default | Notes |
|------|---------|-------|
| `--blast-ref PATH` | — | Reference transcriptome FASTA. Without this, no specificity screen is run. |
| `--positive-selection-strategy {any-strong-hit, best-coverage, specific-id}` | `any-strong-hit` | See below. |
| `--target-transcript-id ID [ID ...]` | — | One or more target IDs. Required for `specific-id`; also added to the auto-negative-screen exclusion list (so paralogs/antisense contigs of the same gene are not treated as off-targets). |
| `--min-bitscore` | `75` | Bitscore cutoff for "strong" hits. |
| `--max-evalue` | `1e-10` | E-value cutoff for "strong" hits. |
| `--blast-extra-args 'flags…'` | — | Extra flags appended to every `blastn` invocation. |
| `--blast-negative-ref PATH` | auto-derived | Optional explicit negative reference. |
| `--negative-bitscore` / `--negative-evalue` | inherits | Override thresholds for the negative screen only. |

#### Positive-selection strategies

- `any-strong-hit` — keep any probe with at least one BLAST hit passing the bitscore/e-value cutoffs. Fast and permissive; appropriate when the reference is curated and small.
- `best-coverage` — rank reference transcripts by (probe-breadth, mean bitscore), pick the best, and keep probes that hit only that transcript. Good default for non-model organisms with paralogs.
- `specific-id` — keep only probes that hit at least one of the supplied `--target-transcript-id`s and nothing else.

### `hcr-prober isoform-split`

Design probes for common and unique regions of an isoform group in one command. All `design` flags apply.

| Flag | Default | Notes |
|------|---------|-------|
| `-i / --input` | — (required) | FASTA containing all isoforms for one or more genes. |
| `--gene-prefix P [P ...]` | — (required) | Gene-prefix(es) shared across each gene's isoform IDs. |
| `--delimiter S` | `_` | Character separating the prefix from the isoform identifier (e.g. `MyGene_iso1` → prefix `MyGene`). |
| `--common-strategy` | `any-strong-hit` | Selection strategy for common-region probes. |
| `--unique-strategy` | `best-coverage` | Selection strategy for isoform-unique probes. |

The negative screen for the common-region job auto-excludes **all** sister isoforms (so the screen does not reject probes whose entire purpose is to bind every isoform). Unique-region jobs exclude every isoform except the one being targeted.

### `hcr-prober swap`

Replace the amplifier handles on an existing order sheet without re-running the design pipeline.

| Flag | Default | Notes |
|------|---------|-------|
| `--input-probes PATH` | — (required) | A single `.xlsx` order sheet OR a directory containing one or more. |
| `--output-dir DIR` | `swapped_probes` | Output directory. |
| `--new-amplifier ID` | — (required) | New amplifier (e.g. `B5`). |

Sequences are validated by length (handle + spacer + 25-nt arm) before stripping; rows whose layout does not match are passed through unchanged with a warning. Round-trips are deterministic given the same `--seed`.

---

## Hybridisation buffer conditions

HCR v3.0 hybridises in **5×SSC + 50 % formamide at 37 °C**. Reported Tm and ΔG values reflect those conditions by default, so the `--min-tm`/`--max-tm` range you set is the Tm your probes will actually see in the assay.

If you want to tune in PCR conditions (50 mM Na⁺, no formamide) — for example to reproduce values from a primer-design tool — pass `--buffer-preset pcr`. You can also override individual values; explicit flags win over the preset.

The Tm correction uses the standard empirical coefficient `Tm_actual = Tm_NN(salt-corrected) − 0.65 × formamide_pct` (Casey & Davidson 1977; McConaughy et al. 1969).

---

## Output files

For each (gene, amplifier) pair, hcr-prober writes the following into `<output-dir>/<gene>/<amplifier>/`:

| File | Contents |
|------|----------|
| `<gene>_<amp>_order.xlsx` | Two columns ("Pool name", "Sequence"). Ready to paste into an oligo-order form. |
| `<gene>_<amp>_probes.fasta` | Each pair as two FASTA records (`<id>_A` for the dn probe, `<id>_B` for the up probe). |
| `<gene>_<amp>_probe_map.svg` | Visual map: probes drawn as coloured rectangles above the transcript track, one distinct HSL hue per probe. Hover for tooltips. |
| `<gene>_<amp>_summary.txt` | Run provenance (command, versions, host, timestamp, seed), input parameters, the full filtering funnel, and a detailed BLAST report (probes in the final set marked with `*`). |
| `<gene>_<amp>_details.csv` | Per-probe table: pair_id, pair_num, start/end position on sense, both arms, GC, Tm, hairpin/homodimer/heterodimer ΔG. |

The SVG, summary, and CSV together let you reproduce or audit any run from the file system alone.

---

## Examples

### 1. Default HCR run with BLAST screen

```bash
hcr-prober design \
    -i MyGene.fasta \
    -o probes/MyGene/ \
    --amplifier B1 \
    --blast-ref host_transcriptome.fasta \
    --threads 4
```

### 2. Multiple amplifiers in one shot

```bash
hcr-prober design \
    -i MyGene.fasta \
    --amplifier B1 B2 B3 \
    --blast-ref host_transcriptome.fasta \
    -o probes/
```

Three amplifier-specific output folders are produced under `probes/MyGene/`.

### 3. Isoform-aware design

For a gene `Opsin1` with two isoforms `Opsin1_A` and `Opsin1_B`:

```bash
hcr-prober isoform-split \
    -i all_opsins.fasta \
    --gene-prefix Opsin1 \
    --amplifier B3 \
    --blast-ref host_transcriptome.fasta \
    -o probes/opsins/
```

You get one `common_probes/Opsin1_common/B3/` folder (probes that bind both isoforms) and per-isoform `isoform_specific_probes/Opsin1_A/B3/` etc. (probes unique to each isoform).

### 4. Tight Tm uniformity for multiplex experiments

```bash
hcr-prober design \
    -i MyGene.fasta \
    --amplifier B1 \
    --blast-ref host_transcriptome.fasta \
    --min-tm 42 --max-tm 50 \
    --max-tm-sigma 1.0 \
    -o probes/
```

After spacing and Tm-window filtering, probes whose Tm is most extreme are dropped one-by-one until σ(Tm) ≤ 1.0 °C across the final pool.

### 5. Running under PCR conditions

```bash
hcr-prober design \
    -i MyGene.fasta \
    --amplifier B1 \
    --buffer-preset pcr \
    --min-tm 60 --max-tm 70 \
    -o probes/
```

### 6. Filter-tuning without paying BLAST cost

```bash
hcr-prober design \
    -i MyGene.fasta \
    --amplifier B1 \
    --blast-ref host_transcriptome.fasta \
    --dry-run \
    --min-gc 35 --max-gc 65 \
    --min-tm 38 --max-tm 52
```

The audit funnel still runs; the BLAST stack is bypassed entirely.

### 7. Excluding paralogs and antisense contigs from the negative screen

If your reference contains the target's antisense contig and a paralog, mark them all as targets so the auto-negative-screen does not reject every probe:

```bash
hcr-prober design \
    -i MyGene_sense.fasta \
    --amplifier B1 \
    --blast-ref host_transcriptome.fasta \
    --target-transcript-id MyGene_sense MyGene_AS MyGene_paralog \
    -o probes/
```

### 8. Swap an existing order to a different amplifier

```bash
hcr-prober swap \
    --input-probes probes/MyGene/B1/MyGene_B1_order.xlsx \
    --new-amplifier B5 \
    --output-dir probes_b5/
```

Or pass a directory to swap every `.xlsx` inside it.

---

## Reproducibility

`--seed N` (default `0`) seeds Python's `random` and `numpy.random` at the top of `main()`. With the seed fixed, IUPAC spacer resolution, subsampling, and any other randomised step produce byte-identical output across runs.

The summary file embeds every piece of context needed to reproduce the run: the exact command, hcr-prober / BLAST+ / primer3 versions, host, UTC timestamp, and seed.

---

## Amplifier sources

The package ships with two JSON files defining HCR v3.0 split-initiator handles:

- `hcr_prober/config/amplifiers/HCR_v3_Choi_2018.json` — **B1, B2, B3, B4, B5**, verified against Choi, Beck, Pierce 2014 *ACS Nano* 8(5):4284-4294 SI page 43 (S6); architecture as published in Choi et al. 2018 *Development* 145:dev165753.
- `hcr_prober/config/amplifiers/HCR_v3_Wang_2020.json` — **B7, B9**, attributed to Wang et al. 2020. **These have not yet been independently cross-checked against the Wang 2020 supplement** in the current release; treat as plausible-but-unverified until you confirm against the source.

Each amplifier entry carries a `_source` field with its citation. Adding a new amplifier is a matter of dropping in another JSON file with `up`, `dn`, `upspc`, `dnspc`, and `_source` — no code change required.

---

## What's new in v1.12.0

This release is a comprehensive audit of correctness, scientific defaults, robustness, and reproducibility. The full per-task plan lives at `docs/plans/2026-05-03-hcr-prober-audit-fixes.md`.

**Correctness fixes (output-affecting)**
- Per-arm GC range filter (was applied to the joined 52-mer; arms with extreme GC could pass).
- Reject windows containing N or other IUPAC ambiguity codes before downstream Tm/primer3/BLAST.
- IUPAC spacer resolution (B7/B9 use `ww` = A/T) is now done once at amplifier load under the seeded RNG; previously each finalize call re-randomised, breaking reproducibility and swap round-trips.
- Negative BLAST screen for `isoform-split` now propagates the full sister-isoform set into the exclusion list (previously every common-region probe was rejected).
- Negative BLAST DB path was being passed as a directory rather than the DB-name prefix; the screen failed silently. Fixed.
- BLAST `-strand minus` (was `-strand plus`, which produced 0 on-target hits because probes are antisense to the mRNA reference).

**Scientific defaults aligned with HCR**
- Default `--na-conc` is now 825 mM (5×SSC) instead of 50 mM (PCR).
- New `--formamide-pct` flag (default 50 %) corrects Tm by −0.65 °C per percent.
- New `--buffer-preset {hcr-5xssc, pcr}` convenience flag.
- primer3 secondary-structure calls receive HCR ionic conditions instead of primer3's PCR-like defaults.
- `--skip-5prime` default reduced from 100 to 50.

**New flags**
- `--seed`, `--threads`, `--dry-run`, `--verbose`, `--quiet`, `--max-tm-sigma`, `--target-transcript-id` (now accepts multiple values).

**Robustness**
- BLAST DB cache key now includes a SHA-256 of the absolute reference path, preventing same-basename collisions in a shared `--db-path`.
- `read_fasta` raises a clear `ValueError` on empty / malformed input (was silently returning `{}`).
- swap matches probes by exact `len(handle)+len(spacer)+25` length and by handle position (5′ for `up`, 3′ for `dn`); rows that don't match are passed through unchanged.
- `check_dependencies` parses and logs the BLAST+ version; warns below 2.10.

**Provenance & reproducibility**
- Summary file opens with a Provenance block: command, hcr-prober / BLAST+ / primer3 versions, host, UTC timestamp, seed.
- Details CSV adds `end_pos_on_sense` alongside `start_pos_on_sense`.
- Glob-driven file iteration is sorted, so amplifier-load order and directory-mode swap are stable across machines.

**Test coverage**
- 38 → 104 tests covering correctness fixes, edge cases, end-to-end BLAST integration, isoform-split with real BLAST, swap workflow, and provenance metadata.

**Algorithmic / UX**
- `subsample_probes` no longer pins both endpoints (was `np.linspace`, which always included index 0 and n−1); now picks bucket midpoints.
- SVG probe map uses one distinct HSL hue per probe (was a 6-colour recycled palette).

**Deferred to future releases** — Probe `@dataclass` refactor and full type-hint sweep (Phase 6.1, 6.2); soft-penalty DP weight (Phase 7.2); MAFFT-based MSA path for many-isoform genes (Phase 7.3); independent verification of B7/B9 against the Wang 2020 supplement (Phase 1.6, blocked on supplement access).

---

## License

See `LICENSE` for terms.
