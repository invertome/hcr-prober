# HCR-prober v1.13.3

A command-line pipeline for designing DNA probes for third-generation *in situ* Hybridization Chain Reaction (HCR v3.0). Given an mRNA target, hcr-prober produces an order-ready set of split-initiator probe pairs that tile the transcript, with optional BLAST-based specificity screening, secondary-structure filtering, and a fully reproducible audit trail.

Designed for molecular biologists working with model and non-model organisms, including complex isoform families.

---

## Table of contents

1. [Quickstart](#quickstart)
2. [Installation](#installation)
3. [How the pipeline works](#how-the-pipeline-works)
4. [Commands](#commands)
5. [Hybridization buffer conditions](#hybridisation-buffer-conditions)
6. [Output files](#output-files)
7. [Examples](#examples)
8. [Reproducibility](#reproducibility)
9. [Amplifier sources](#amplifier-sources)
10. [What's new in v1.13.3](#whats-new-in-v1133)

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

By default hcr-prober makes no Tm-related adjustments and the Tm filter is off — pass `--min-tm` / `--max-tm` to opt into a Tm window, and pass `--buffer-preset {hcr-5xssct-urea, hcr-5xssc, pcr}` to apply a denaturant correction (urea-HCR, formamide HCR, or PCR-style).

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
| 8 | Tm filter | **Off by default** — pass `--min-tm` / `--max-tm` to opt in. When enabled, both arms must have a (denaturant-corrected) Tm inside the window. |
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
| `--min-tm` / `--max-tm` | off / off (°C) | Per-arm Tm range, **after denaturant correction**. Off by default; pass either bound to enable filtering. Typical urea-HCR window: `--min-tm 60 --max-tm 75`. |
| `--max-homopolymer N` | `4` | Reject windows with runs ≥ N+1 of any base. |
| `--max-tm-sigma N` | off | Optional: enforce σ(Tm) ≤ N across the final pool. |

#### Hybridization buffer
| Flag | Default | Notes |
|------|---------|-------|
| `--buffer-preset {hcr-5xssct-urea, hcr-5xssc, pcr}` | unset | Optional convenience preset. When unset, no denaturant correction is applied. |
| `--na-conc` | `825` mM | Sodium concentration. 5×SSC(T) ≈ 825 mM Na⁺. `pcr` preset sets to 50. |
| `--mg-conc` | `0` mM | Magnesium concentration. |
| `--dntp-conc` | `0` mM | dNTP concentration. |
| `--dna-conc` | `25` nM | Probe oligo concentration. |
| `--formamide-pct` | `0` % | Formamide percent. `hcr-5xssc` preset sets to 50. Tm is reduced by 0.65 °C per percent. |
| `--urea-M` | `0` M | Urea molarity. `hcr-5xssct-urea` preset sets to 4. Tm is reduced by 2.25 °C per molar. |

#### Secondary structure
| Flag | Default | Notes |
|------|---------|-------|
| `--max-hairpin-dg` | `-3.0` (kcal/mol) | Reject arm if hairpin ΔG is more negative. |
| `--max-homodimer-dg` | `-5.0` (kcal/mol) | Reject arm if homodimer ΔG is more negative. |
| `--max-heterodimer-dg` | `-5.0` (kcal/mol) | Reject pair if cross-arm heterodimer ΔG is more negative. |

All three values are evaluated under the same Na⁺/Mg²⁺/formamide/urea conditions as the Tm filter.

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

## Hybridization buffer conditions

By default hcr-prober makes no buffer-related Tm adjustments. Reported Tm is the **salt-corrected nearest-neighbor Tm at 5×SSCT** (Na⁺ ≈ 825 mM, no Mg²⁺, no denaturant). The Tm-window filter is also off by default — HCR's split-initiator architecture is tolerant enough that GC and homopolymer constraints alone are usually sufficient.

You can opt into denaturant corrections by passing `--buffer-preset`:

| `--buffer-preset` | Na⁺ | Mg²⁺ | Formamide | Urea | When to use |
|---|---|---|---|---|---|
| *(unset, default)* | 825 mM | 0 mM | 0 % | 0 M | No Tm adjustments. Reported Tm = salt-corrected NN. |
| `hcr-5xssct-urea` | 825 mM | 0 mM | 0 % | 4 M | Urea HCR (Choi et al. 2018, *Development* 145:dev165753). |
| `hcr-5xssc` | 825 mM | 0 mM | 50 % | 0 M | Original formamide HCR (Choi 2014, 2018). |
| `pcr` | 50 mM | 0 mM | 0 % | 0 M | DNA primers/probes for PCR-style assays. |

Individual flags (`--na-conc`, `--mg-conc`, `--formamide-pct`, `--urea-M`) win over the preset, so you can mix and match. To opt into a Tm filter, pass `--min-tm` / `--max-tm`; with the urea preset a typical window is 60–75 °C, with the formamide preset 40–55 °C, with no preset 70–85 °C.

The Tm correction (when a denaturant is present) is

```
Tm_actual = Tm_NN(salt-corrected) − 0.65 × formamide_pct − 2.25 × urea_M
```

(Casey & Davidson 1977; McConaughy et al. 1969 for formamide; Hutton 1977 for urea. The two corrections are independent and additive.)

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

Each example is self-contained: paste the command, point it at your own FASTA files, and inspect the output files described in [Output files](#output-files). If you're new to HCR probe design, work through Example 1 in full before skimming the rest — every later example builds on the same mental model.

You'll typically need two FASTA files:

- **Target FASTA** (e.g. `MyGene.fasta`) — one or more mRNA records, sense strand. If the file has more than one record, hcr-prober processes all of them; use `--gene-name` to pick a single record.
- **Reference transcriptome** (e.g. `host_transcriptome.fasta`) — every transcript in the expression background of your organism. hcr-prober uses this both to confirm each probe hits its intended target (the positive screen) and to reject probes that *also* hit something else (the auto-derived negative screen).

A reference transcriptome is optional but strongly recommended: without `--blast-ref`, no specificity check is performed and probes might bind paralogs or unrelated transcripts.

### 1. Default HCR run with BLAST screen

The most common scenario: you have one mRNA target, you want a tilable set of probe pairs that work under standard urea-HCR conditions (5×SSCT + 4 M urea + 37 °C), and you have a transcriptome to screen against.

```bash
hcr-prober design \
    -i MyGene.fasta \
    -o probes/MyGene/ \
    --amplifier B1 \
    --blast-ref host_transcriptome.fasta \
    --threads 4
```

What each flag does:

- `-i MyGene.fasta` — the target. The first record (or the one matching `--gene-name`) becomes the design target.
- `-o probes/MyGene/` — the top-level output directory. hcr-prober creates `probes/MyGene/<gene>/B1/` underneath it and writes the five output files there (`order.xlsx`, `probes.fasta`, `probe_map.svg`, `summary.txt`, `details.csv`).
- `--amplifier B1` — uses the B1 split-initiator handles from the Choi 2018 amplifier set. See [Amplifier sources](#amplifier-sources) for the full list and provenance.
- `--blast-ref host_transcriptome.fasta` — runs blastn against this file to (a) confirm each probe pair hits the intended target with bitscore ≥ 75 and (b) reject any probe that *also* hits a non-target transcript at the same bitscore threshold.
- `--threads 4` — passes `-num_threads 4` to blastn. BLAST is the slow step; this scales nearly linearly with available cores.

When the run finishes, open `<gene>_B1_summary.txt` first. It shows the filtering funnel (how many candidate windows survived each filter), the exact command, software versions, host, and seed needed to reproduce the run. `<gene>_B1_probe_map.svg` shows where the surviving probes landed along the transcript; `<gene>_B1_order.xlsx` is what you upload to your oligo vendor.

### 2. Multiple amplifiers in one shot

In multiplexed HCR, each gene gets a different amplifier so each gets a different fluorophore. You can design probe sets for multiple amplifiers against the same target in a single command — useful when you haven't decided which channel to use yet, or when you want to keep your options open.

```bash
hcr-prober design \
    -i MyGene.fasta \
    --amplifier B1 B2 B3 \
    --blast-ref host_transcriptome.fasta \
    -o probes/
```

You'll get three sibling output directories: `probes/<gene>/B1/`, `probes/<gene>/B2/`, `probes/<gene>/B3/`. The probe target sequences are identical across the three; only the appended HCR initiator handles differ. (If you only realize later that you needed a different amplifier, see Example 8 — swapping is cheaper than re-running.)

### 3. Isoform-aware design

Many genes have several mRNA isoforms (alternative splice forms). Depending on your biology you may want to:

- visualize *all* isoforms with one set of probes — target the regions they share ("common probes"),
- distinguish individual isoforms — target the regions unique to each ("isoform-specific probes"),
- or both.

`isoform-split` runs a pairwise alignment across the input FASTA, partitions each transcript into common and unique regions, and runs the full design pipeline on each region category separately.

For a gene `Opsin1` with two isoforms `Opsin1_A` and `Opsin1_B`:

```bash
hcr-prober isoform-split \
    -i all_opsins.fasta \
    --gene-prefix Opsin1 \
    --amplifier B3 \
    --blast-ref host_transcriptome.fasta \
    -o probes/opsins/
```

`--gene-prefix Opsin1` tells the tool to group every record whose ID starts with `Opsin1_` into one isoform set. Pass multiple prefixes (`--gene-prefix Opsin1 Opsin2`) to process several genes in one run. The output layout is:

```
probes/opsins/
├── common_probes/Opsin1_common/B3/         # binds both isoforms
└── isoform_specific_probes/
    ├── Opsin1_A/B3/                        # binds only Opsin1_A
    └── Opsin1_B/B3/                        # binds only Opsin1_B
```

If a region is too short for any valid 52 nt window (for example, a unique stretch only 30 nt long), no probe is reported for that region — check the relevant `*_summary.txt` to see which slots came up empty and why.

### 4. Tight Tm uniformity for multiplex experiments

When you pool probes from many genes into one hybridization, every probe needs to anneal at the same temperature. The Tm filter is off by default; pass `--min-tm` / `--max-tm` to opt into a window. If your protocol assumes a single hybridization Tm across genes, you also want a tighter spread.

```bash
# urea HCR (default preset): typical window is 60-75 C
hcr-prober design \
    -i MyGene.fasta \
    --amplifier B1 \
    --blast-ref host_transcriptome.fasta \
    --min-tm 64 --max-tm 70 \
    --max-tm-sigma 1.0 \
    -o probes/
```

`--min-tm` / `--max-tm` define the per-probe Tm window (under whichever buffer preset is active). `--max-tm-sigma 1.0` adds a *post-selection* step: after the normal pipeline has chosen its set, hcr-prober drops the probe with the most extreme Tm, recomputes σ(Tm), and repeats until either σ ≤ 1.0 °C or fewer than four probes remain. The summary file logs how many probes were peeled and the final σ.

For formamide-HCR users, the equivalent window is ~23 °C lower (the formamide correction subtracts 32.5 °C vs urea's 9 °C):

```bash
hcr-prober design \
    -i MyGene.fasta \
    --amplifier B1 \
    --buffer-preset hcr-5xssc \
    --blast-ref host_transcriptome.fasta \
    --min-tm 42 --max-tm 50 \
    --max-tm-sigma 1.0 \
    -o probes/
```

Trade-off: tighter σ usually means fewer probes. If `--max-probes 33` was already hit before this filter, expect to lose several here.

### 5. Running under PCR conditions

By default hcr-prober reports thermo values for HCR hyb buffer (high salt, high formamide). If you're using this tool to design DNA primers or probes for a PCR-style assay (50 mM Na⁺, no formamide, higher annealing temperatures), switch buffer presets:

```bash
hcr-prober design \
    -i MyGene.fasta \
    --amplifier B1 \
    --buffer-preset pcr \
    --min-tm 60 --max-tm 70 \
    -o probes/
```

`--buffer-preset pcr` sets `--na-conc 50` and `--formamide-pct 0` automatically. The reported Tm values then reflect the PCR buffer rather than HCR hyb buffer, so your `--min-tm` / `--max-tm` window should be raised accordingly (PCR Tms run roughly 20 °C higher than the corresponding HCR Tms). If you need finer control, override `--na-conc`, `--mg-conc`, `--dntp-conc`, `--dna-conc`, or `--formamide-pct` directly — explicit values always win over the preset.

### 6. Filter-tuning without paying BLAST cost

BLAST is the slow step. While you're still tuning GC% and Tm windows for an unfamiliar transcript, you don't need it — you just want to see how many candidate windows survive your sequence/thermo filters.

```bash
hcr-prober design \
    -i MyGene.fasta \
    --amplifier B1 \
    --blast-ref host_transcriptome.fasta \
    --dry-run \
    --min-gc 35 --max-gc 65 \
    --min-tm 38 --max-tm 52
```

`--dry-run` runs every filter (length, IUPAC, homopolymer, GC, Tm, secondary structure) and reports the full funnel in `*_summary.txt`, but skips BLAST entirely. Iterate quickly on the GC/Tm windows until the funnel shows a healthy candidate count, then drop `--dry-run` for the real run.

### 7. Excluding paralogs and antisense contigs from the negative screen

The auto-derived negative screen treats every transcript in `--blast-ref` that *isn't* explicitly named as a potential off-target. That bites you in two common cases:

- the reference transcriptome contains both the sense and antisense assembly of your target locus,
- the reference contains the target's paralogs or alternative isoforms.

Without help, every probe will hit those near-identical contigs and get rejected as non-specific — and the run will report zero surviving probes for no obvious reason.

```bash
hcr-prober design \
    -i MyGene_sense.fasta \
    --amplifier B1 \
    --blast-ref host_transcriptome.fasta \
    --target-transcript-id MyGene_sense MyGene_AS MyGene_paralog \
    -o probes/
```

`--target-transcript-id` accepts multiple IDs. Each one is added to the negative-screen exclusion set, so probes that hit any of them are kept rather than discarded. (The first ID is also used by the `specific-id` positive-selection strategy if you opt into it via `--positive-selection-strategy specific-id`.)

If you're not sure which IDs to pass, run the design once without this flag and look at the negative-screen rejection table in `*_summary.txt` — the most-hit reference contigs are usually the culprits.

### 8. Swap an existing order to a different amplifier

You already designed and validated probes against your target with B1. Now you want the same probes in B5 — for instance, to combine them with another B1-labeled gene in a multiplex panel. There's no need to re-run `design`; the binding sequence is the same, only the appended initiator handles change.

```bash
hcr-prober swap \
    --input-probes probes/MyGene/B1/MyGene_B1_order.xlsx \
    --new-amplifier B5 \
    --output-dir probes_b5/
```

`swap` reads the order spreadsheet, strips the B1 handles + spacer from each oligo, prepends the matching B5 handles + spacer, and writes a new spreadsheet. Pass a directory to `--input-probes` to swap every `.xlsx` inside it in one go. The swap is round-trip-safe: B1 → B5 → B1 reproduces the original sequences byte-for-byte (locked in by a regression test), so this is always a non-destructive operation.

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

## What's new in v1.13.3

Bug fix: hcr-prober failed to import at all on Python 3.8–3.11 with `SyntaxError: f-string: expecting '}'`. Six f-strings in `file_io.py`, `isoform_analyzer.py`, and `blast_wrapper.py` reused the same quote character for the outer string and an inner expression — a syntax that PEP 701 only legalised in Python 3.12+. setup.py declares `python_requires>=3.8`, so the code must remain parseable on 3.8–3.11. Each offending f-string was rewritten by switching the outer quote (or pre-computing the inner expression). New static-analysis test `test_python_compat.py` scans every shipped `.py` file for the same-quote-nested pattern and fails if any are reintroduced.

## What's new in v1.13.2

Bug fix: `hcr-prober design -h` and `hcr-prober isoform-split -h` errored out with `ValueError: unsupported format character` because a literal `%` in the `--formamide-pct` help string collided with argparse's printf-style help-string interpolation. The `%` is now properly escaped (`%%`) so subcommand help renders cleanly. New regression test in `test_integration.py` invokes `-h` for every subcommand.

## What's new in v1.13.1

`--buffer-preset` is now **unset by default**. Without an explicit preset, hcr-prober makes no denaturant-related Tm adjustments and reports the salt-corrected NN Tm at 5×SSCT. Users opt into urea or formamide corrections by passing `--buffer-preset hcr-5xssct-urea` (or `hcr-5xssc`). Individual flags (`--na-conc`, `--formamide-pct`, `--urea-M`, `--mg-conc`) gain explicit defaults (825 / 0 / 0 / 0) so they work standalone.

Why: keeps default behaviour minimal and predictable. The tool no longer assumes which buffer you use, matching the philosophy of simpler community probe-design tools.

## What's new in v1.13.0

Default behaviour now reflects contemporary urea-HCR practice (Choi et al. 2018) rather than the older formamide protocol, and the Tm filter is opt-in instead of always-on. Both changes align hcr-prober with how most HCR labs actually run experiments today.

**Default changes (behaviour-affecting)**
- New default `--buffer-preset hcr-5xssct-urea` (5×SSCT + 4 M urea, 37 °C). Reported Tm and ΔG values now reflect urea-HCR conditions out of the box.
- `--min-tm` and `--max-tm` default to **off**. Pass either bound to opt into Tm-window filtering. Other filters (GC, homopolymer, secondary structure, BLAST) are unchanged.

**New flags & options**
- `--urea-M` (default from preset) — urea molarity, with empirical Tm correction of −2.25 °C / M (Hutton 1977).
- New buffer preset `hcr-5xssct-urea`. The previous default `hcr-5xssc` (50 % formamide) and `pcr` remain as alternates.

**Migration**

| You want                                          | Add to your command                                                |
|---------------------------------------------------|--------------------------------------------------------------------|
| v1.12.0 reporting **and** filtering               | `--buffer-preset hcr-5xssc --min-tm 40 --max-tm 55`                |
| v1.12.0 filter only, with v1.13.0 urea reporting  | `--min-tm 60 --max-tm 75` (urea-corrected equivalent of 40–55 °C)  |
| Stay with formamide reporting, no Tm filter       | `--buffer-preset hcr-5xssc`                                        |

The 23.5 °C offset between the two windows comes directly from the difference in denaturant correction (50 × 0.65 = 32.5 °C for 50 % formamide, vs. 4 × 2.25 = 9 °C for 4 M urea).

**Documentation**
- Hyb-conditions section rewritten around urea-HCR with formamide as alternate.
- Tm correction now: `Tm_actual = Tm_NN(salt-corrected) − 0.65 × formamide_pct − 2.25 × urea_M` (corrections are independent and additive).

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

**Algorithmic / UX**
- `subsample_probes` no longer pins both endpoints (was `np.linspace`, which always included index 0 and n−1); now picks bucket midpoints.
- SVG probe map uses one distinct HSL hue per probe (was a 6-colour recycled palette).

---

## License

See `LICENSE` for terms.
