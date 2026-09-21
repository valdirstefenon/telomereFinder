# telomereFinder

**A Universal GUI tool for detecting and organizing telomeric sequences in any *de novo* assembled genome**

[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](#license)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](#system-requirements)
[![Version](https://img.shields.io/badge/version-1.0-orange.svg)](#version-history)
[![GUI: Tkinter](https://img.shields.io/badge/GUI-Tkinter-orange.svg)](#installation-guide)

**telomereFinder** is a Python-based tool that identifies telomeric repeats in genome assemblies and relocates missassembled telomeric sequences to the correct positions at the scaffold ends. Misassemblies frequently place telomeric repeats in the middle of a scaffold instead of at its terminus; this tool detects those events and produces a cleaned, reorganized multi-FASTA.

The tool is organism-agnostic and ships with curated presets for the major eukaryotic telomere types:

- Plants (`CCCTAAA` / `TTTAGGG`)
- Vertebrates (`TTAGGG` / `CCCTAA`)
- Insects (`TTAGG` / `CCTAA`)
- Nematodes (`TTAGGC` / `GCCTAA`)
- Ciliates (`TTGGGG` / `CCCCAA`)
- Fungi (*Candida*, *Neurospora*)
- Trypanosomes (`TTAGGG` / `CCCTAA`)
- Custom — user-defined motifs

The pipeline combines:

- Canonical and fuzzy motif detection for telomeric repeats
- Automatic generation of fuzzy variants from any canonical motif
- Automatic deduplication of overlapping regex hits
- Terminal vs. internal classification based on a user-defined margin
- Orientation-aware relocation — forward repeats to the start, reverse repeats to the end
- Motif suggestion by scanning scaffold ends
- Detailed reporting with per-scaffold statistics and a full movement log

This makes **telomereFinder** suitable both for assembly quality control and for preprocessing genomes before downstream analyses (e.g., chromosome-level scaffolding, telomere-to-telomere assembly validation).

---

## Table of Contents

- [Overview](#overview)
- [What telomereFinder Does](#what-telomerefinder-does)
- [New in Version 1.0](#new-in-version-10)
- [Pipeline Workflow](#pipeline-workflow)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Organism Presets](#organism-presets)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Using the Graphical Interface](#using-the-graphical-interface)
- [Understanding the Results](#understanding-the-results)
- [Understanding the Dry Run Mode](#understanding-the-dry-run-mode)
- [Understanding Deduplication](#understanding-deduplication)
- [Motif Suggestion](#motif-suggestion)
- [Troubleshooting](#troubleshooting)
- [FAQ](#frequently-asked-questions)
- [Citation](#citation)
- [License](#license)
- [Version History](#version-history)

---

## Overview

Genome assemblies produced by long-read or hybrid pipelines frequently contain telomeric repeats at internal positions of a scaffold. This usually indicates:

- Two scaffolds were incorrectly joined at a telomere.
- A single chromosome was split, and the two halves were merged with the wrong orientation.
- A haplotig was collapsed into the primary assembly.

**telomereFinder** scans every sequence of a multi-FASTA, identifies these internal telomeric repeats, and relocates them to the appropriate scaffold terminus based on their orientation:

- Forward motifs → moved to the start of the scaffold.
- Reverse motifs → moved to the end of the scaffold.

The result is a reorganized multi-FASTA where internal telomeres are no longer embedded in the middle of sequences, plus a detailed textual report documenting every relocation.

The tool works with any organism whose telomeres are tandem repeats: the user selects a preset or enters custom motifs, and the pipeline automatically generates fuzzy variants and runs the analysis.

---

## What telomereFinder Does

### Core Functions

#### Telomere Detection

- Scans each scaffold for canonical and fuzzy telomeric motifs.
- Accepts any tandem-repeat motif in ACGTN alphabet.
- Automatically generates fuzzy variants by allowing IUPAC substitutions at positions 2 and 5 of the motif.
- Supports a configurable minimum number of consecutive repeat units.

#### Internal vs. Terminal Classification

- Classifies every detected telomeric region as terminal (within `margin` bp of a scaffold end) or internal (ITS).
- Terminal telomeres are left untouched.
- Internal telomeres are flagged for relocation.

#### Orientation-Aware Relocation

- Forward ITSs are moved to the 5-prime start of the scaffold.
- Reverse ITSs are moved to the 3-prime end of the scaffold.
- Multiple ITSs on the same scaffold are handled simultaneously, with the body of the sequence preserved in its original order.

#### Deduplication of Overlapping Hits

- Multiple regex patterns (canonical + fuzzy) can hit the same physical telomere multiple times.
- Overlapping hits are merged, keeping the longest region per locus.

#### Motif Suggestion

- Optional scan of scaffold ends to suggest the most likely organism preset.
- Ranks all presets by the number of tandem arrays found in terminal regions.
- Updates the GUI automatically with the best match.

#### Reporting

- Per-scaffold summary (length, terminal telomeres, ITS count, bases relocated).
- Global statistics (scaffolds affected, total bases relocated).
- Detailed table of every ITS relocation with its original position, orientation, destination, and matching pattern.
- Record of the organism preset and motif set used.

#### Graphical Interface

- Modern orange-themed Tkinter GUI.
- Threaded execution — the interface remains responsive.
- Real-time log console.
- Progress bar.
- Dry-run mode for inspection without writing files.

---

## New in Version 1.0

- First public release of **telomereFinder**.
- Full support for any organism via presets or custom motifs.
- Automatic generation of fuzzy variants from user-provided motifs.
- Overlap deduplication.
- Motif suggestion based on terminal-region scanning.
- Graphical user interface with dry-run mode.
- Detailed text report with global and per-scaffold statistics.

---

## Pipeline Workflow

```text
+-----------------------------------------------------------------------------+
|                              INPUT FILE                                     |
|                          multi-FASTA genome                                 |
+-------------------------------+---------------------------------------------+
                                |
                                v
+-----------------------------------------------------------------------------+
|                     ORGANISM PRESET SELECTION                               |
|  - User picks a preset (Plants, Vertebrates, Insects, ...)                  |
|  - Or enters custom forward / reverse motifs                                |
|  - Optional: Suggest Motifs scans scaffold ends to auto-pick a preset       |
+-----------------------------------------------------------------------------+
                                |
                                v
+-----------------------------------------------------------------------------+
|                          FASTA PARSING                                      |
|  - Read all sequences                                                       |
|  - Uppercase nucleotides                                                    |
|  - Preserve original order                                                  |
+-----------------------------------------------------------------------------+
                                |
                                v
+-----------------------------------------------------------------------------+
|                       TELOMERE DETECTION                                    |
|  - Canonical patterns from selected motif set                               |
|  - Fuzzy variants generated automatically                                   |
|  - Minimum repeat threshold                                                 |
|  - Minimum ITS length filter                                                |
+-----------------------------------------------------------------------------+
                                |
                                v
+-----------------------------------------------------------------------------+
|                        DEDUPLICATION                                        |
|  - Merge overlapping hits                                                   |
|  - Keep longest region per locus                                            |
+-----------------------------------------------------------------------------+
                                |
                                v
+-----------------------------------------------------------------------------+
|                    TERMINAL vs. INTERNAL CLASSIFICATION                     |
|  - Region within margin bp of 5-prime end  -> terminal (start)              |
|  - Region within margin bp of 3-prime end  -> terminal (end)                |
|  - Otherwise                               -> internal (ITS)                |
+-----------------------------------------------------------------------------+
                                |
                                v
+-----------------------------------------------------------------------------+
|                       SCAFFOLD REORGANIZATION                               |
|  - Remove all internal ITSs from the body                                   |
|  - Prepend forward ITSs (ascending position order)                          |
|  - Append reverse ITSs (ascending position order)                           |
+-----------------------------------------------------------------------------+
                                |
                                v
+-----------------------------------------------------------------------------+
|                        REPORTING AND OUTPUT                                 |
+-----------------------------------------------------------------------------+
|  - Reorganized multi-FASTA                                                  |
|  - Text report with per-scaffold statistics                                 |
|  - Detailed ITS movement table                                              |
+-----------------------------------------------------------------------------+
```

---

## System Requirements

### Minimum Hardware

| Component | Recommendation |
| --- | --- |
| CPU | 2+ cores recommended |
| RAM | 4 GB minimum; 8 GB+ recommended for large genomes |
| Storage | Enough space for the input and output FASTA files |

### Software Dependencies

| Requirement | Details |
| --- | --- |
| Operating System | Linux, macOS, or Windows |
| Python | 3.8 or higher |
| Tkinter | Bundled with most Python distributions (see installation notes for minimal Linux distros) |

### Required Python Modules

All modules used are part of the Python standard library:

| Module | Purpose |
| --- | --- |
| re | Regex-based motif detection |
| tkinter | Graphical user interface |
| threading | Background execution |
| queue | GUI thread communication |
| csv | Report generation |
| os, sys, datetime, collections | Utilities |

No `pip install` is required.

---

## Installation Guide

### Step 1 — Verify Python version

```bash
python3 --version
```

Python 3.8 or higher is required.

### Step 2 — Ensure Tkinter is available

On most systems, Tkinter is bundled with Python. Verify with:

```bash
python3 -c "import tkinter; print('Tkinter OK')"
```

If the import fails on a minimal Linux distribution:

```bash
# Debian / Ubuntu
sudo apt-get install python3-tk

# Fedora / RHEL
sudo dnf install python3-tkinter

# Arch
sudo pacman -S tk
```

### Step 3 — Clone the repository

```bash
git clone https://github.com/<your-username>/telomereFinder.git
cd telomereFinder
```

### Step 4 — (Optional) Create a Conda environment

```bash
conda create -n telomerefinder python=3.10 -y
conda activate telomerefinder
conda install -c conda-forge tk -y
```

### Step 5 — Launch the application

```bash
python3 telomereFinder.py
```

---

## Organism Presets

**telomereFinder** ships with curated presets for the major eukaryotic telomere types. Selecting a preset automatically fills the forward and reverse motif fields, which remain read-only unless Custom is selected.

| Preset | Forward motifs | Reverse motifs | Notes |
| --- | --- | --- | --- |
| Plants | `CCCTAAA`, `CCCTAA` | `TTTAGGG`, `TTAGGG` | Arabidopsis-type plant telomeres |
| Vertebrates | `TTAGGG` | `CCCTAA` | Mammals, birds, fish, reptiles, amphibians |
| Insects (non-Drosophila) | `TTAGG` | `CCTAA` | Excludes Drosophila, which uses retrotransposons |
| Nematodes | `TTAGGC` | `GCCTAA` | Includes Caenorhabditis elegans |
| Ciliates | `TTGGGG` | `CCCCAA` | Includes Tetrahymena |
| Fungi (Candida) | `TGTGGG` | `CCCACA` | Candida-type fungal telomeres |
| Fungi (Neurospora) | `TTAGGG` | `CCCTAA` | Neurospora-type fungal telomeres |
| Trypanosomes | `TTAGGG` | `CCCTAA` | Trypanosoma telomeres |
| Custom | user-defined | user-defined | Any ACGTN motif; fuzzy variants generated automatically |

**Note on overlap between presets.** Some organisms share the same telomeric motif (e.g., Vertebrates, Trypanosomes, and Neurospora all use `TTAGGG` / `CCCTAA`). In those cases, the choice of preset only affects labeling in the report — the detection behavior is identical.

**Note on Drosophila and yeasts.** Drosophila uses retrotransposons (HeT-A, TART, TAHRE) rather than tandem repeats, and Saccharomyces cerevisiae uses an irregular TG(1-3) repeat. Neither is supported by **telomereFinder** without custom logic. Use the Custom preset only when your organism uses a regular tandem-repeat motif.

### Fuzzy Variants

When the **Use fuzzy variants** checkbox is enabled, **telomereFinder** automatically generates variant patterns from each canonical motif by allowing IUPAC substitutions at positions 2 and 5 (0-indexed 1 and 4) of the motif. For example, from `CCCTAAA` the tool generates:

- `C[CT]CTAAA`
- `CCCT[AT]AA`

This increases sensitivity to degenerate telomeres at the cost of a slightly higher false-positive rate.

---

## Input Files

### Mandatory Input

| File | Format | Description |
| --- | --- | --- |
| Input FASTA | `.fasta`, `.fa`, `.fas`, `.fna` | Multi-FASTA genome assembly to analyze and reorganize |

### Optional Output Paths

| Field | Description |
| --- | --- |
| Output FASTA | Destination for the reorganized assembly. Defaults to `<input>_reorganized.fasta`. |
| Report file | Destination for the textual report. Defaults to `<input>_telomere_report.txt`. |

### File Preparation

- Sequences must be in standard FASTA format.
- Each record must start with `>`.
- Only the first whitespace-delimited token after `>` is used as the scaffold identifier.
- Sequences are internally uppercased; case in the input file is not relevant.

---

## Output Files

### 1. Reorganized FASTA

| Property | Description |
| --- | --- |
| Format | Multi-FASTA, 60 characters per line |
| Header | Original scaffold ID, optionally annotated with relocation counts |

Scaffolds that were modified receive an annotated header:

```text
>scaffold_1 [REORGANIZED F:2 R:1]
```

- `F:n` — number of forward ITSs moved to the start.
- `R:n` — number of reverse ITSs moved to the end.

Scaffolds without internal telomeres keep their original header unchanged.

### 2. Text Report

The report is a plain text file with the following sections:

| Section | Contents |
| --- | --- |
| Header | Timestamp, input path, output path, elapsed time |
| Organism and Motifs | Preset name, forward motifs, reverse motifs, fuzzy status |
| Detection Parameters | Min repeats, min ITS length, margin, max ITS per scaffold |
| Global Summary | Total scaffolds analyzed, scaffolds with terminal telomeres, scaffolds with internal telomeres, total ITSs relocated, total bases relocated |
| Per-Scaffold Summary | Table with length, terminal telomere counts, forward/reverse ITS counts, and bases moved per scaffold |
| Detailed ITS Movements | Full table of every relocation with scaffold ID, start, end, length, orientation, destination, and matching pattern |

Example excerpt:

```text
GLOBAL SUMMARY
--------------------------------------------------------------------------------
Total scaffolds analyzed          : 1,247
Scaffolds with terminal telomeres : 812 start / 790 end
Scaffolds with internal telomeres : 23
Total ITSs relocated              : 41
  - Forward moved to START        : 27
  - Reverse moved to END          : 14
Total bases relocated             : 18,432 bp
```

---

## Using the Graphical Interface

### Launching the GUI

```bash
python3 telomereFinder.py
```

### Main Window Layout

#### 1. Files Panel

| Field | Description |
| --- | --- |
| Input FASTA | Browse to select the multi-FASTA to analyze |
| Output FASTA | Browse to choose the destination of the reorganized assembly |
| Report file | Browse to choose the destination of the textual report |

When you select an input file, the output and report paths are auto-filled based on the input filename.

#### 2. Organism Preset Panel

| Control | Description |
| --- | --- |
| Organism group | Dropdown with the curated presets plus Custom |
| Suggest Motifs | Scans scaffold ends and auto-selects the best-matching preset |
| Description | Contextual explanation of the selected preset |
| Forward motifs | Read-only for presets; editable for Custom. Comma-separated. |
| Reverse motifs | Read-only for presets; editable for Custom. Comma-separated. |
| Use fuzzy variants | Auto-generates fuzzy patterns from each canonical motif |

#### 3. Detection Parameters Panel

| Parameter | Default | Description |
| --- | --- | --- |
| Min repeats | `3` | Minimum number of consecutive telomeric repeat units required to call a telomere |
| Min ITS length (bp) | `50` | Minimum length in base pairs for an internal telomeric region to be considered an ITS |
| Terminal margin (bp) | `5000` | Distance from each scaffold end within which a telomere is classified as terminal |
| Max ITS per scaffold | empty = unlimited | Optional cap on the number of ITSs relocated per scaffold |

#### 4. Actions

| Button | Description |
| --- | --- |
| Reorganize | Runs the full pipeline and writes the reorganized FASTA and the report |
| Dry Run | Runs detection only; writes no files |
| Clear Log | Clears the on-screen log console |

#### 5. Progress Bar

Shows the fraction of scaffolds analyzed.

#### 6. Analysis Log

Real-time console output with color-coded messages:

- Green — successful operations
- Red — errors
- Orange bold — section headers
- Dim brown — secondary notes

---

## Understanding the Results

### Classification Rules

A detected region is classified as:

- Start-terminal — if its start coordinate is at most `margin` bp from the 5-prime end.
- End-terminal — if its end coordinate is at most `margin` bp from the 3-prime end.
- Internal (ITS) — otherwise.

Only internal regions are relocated.

### Relocation Rules

- Forward ITSs are inserted, in ascending order of their original position, at the beginning of the scaffold.
- Reverse ITSs are appended, in ascending order of their original position, at the end of the scaffold.
- The sequence between and around the ITSs (the body) preserves its original order and orientation.

### Interpreting the Report

The per-scaffold summary table is the fastest way to spot problematic scaffolds:

- A scaffold with many forward and reverse ITSs likely represents a misjoin between two chromosome arms.
- A scaffold with a single ITS may indicate a split chromosome or a collapsed haplotig.
- A scaffold with no ITSs but with terminal telomeres at both ends is likely a complete chromosome.

---

## Understanding the Dry Run Mode

The **Dry Run** button executes the entire detection pipeline — parsing, motif scanning, deduplication, classification — but writes no files. It is intended as a safety inspection step before applying any modification.

### What Dry Run does

| Action | Dry Run | Reorganize |
| --- | --- | --- |
| Read the input FASTA | yes | yes |
| Detect terminal telomeres | yes | yes |
| Detect internal ITSs | yes | yes |
| Modify sequences | no | yes |
| Write output FASTA | no | yes |
| Write report | no | yes |
| Print movements to log | yes | yes |

### When to use it

- Calibrating parameters before committing to a run.
- Estimating impact — how many scaffolds will be affected and how many bases will move.
- Sanity-checking that no spurious ITS is being captured (e.g., low-complexity regions caught by fuzzy patterns).
- Reviewing the movement list on screen before accepting the result.

The recommended workflow is to always run **Dry Run first**, adjust parameters if needed, then click **Reorganize**.

---

## Understanding Deduplication

The detection engine scans each scaffold with multiple overlapping regex patterns (one canonical plus one or more fuzzy variants per motif). Because a single physical telomere can match several patterns simultaneously, the raw hit list contains multiple entries for the same region.

The deduplication step:

1. Sorts hits by start position, then by length (longest first).
2. Merges any hit that overlaps with an already-accepted region.
3. Keeps the longest region per locus.

Without deduplication:

- Statistics would count the same telomere several times.
- The reorganization step would cut and reinsert the same region multiple times, corrupting the sequence.

Deduplication ensures that each physical telomere is treated as a single entity.

---

## Motif Suggestion

The **Suggest Motifs** button helps when you do not know the telomeric motif of your organism. It performs the following steps:

1. Reads the input FASTA.
2. Extracts the first and last 2000 bp of every scaffold (or the whole scaffold if shorter).
3. Counts occurrences of tandem arrays (at least 3 repeats) for each preset's motifs.
4. Ranks the presets by total count and selects the best match.
5. Updates the GUI dropdown and logs the top 3 candidates with their scores.

The suggestion is intended as a starting point, not a verdict. Terminal-region scans can be confounded by contamination, low-complexity sequence, or incomplete assemblies. Always review the log output and the dry-run results before proceeding.

**Tip.** If the suggestion is ambiguous (two presets with similar scores), inspect the raw counts in the log and check whether the two presets share the same motif (e.g., Vertebrates and Trypanosomes both use `TTAGGG`).

---

## Troubleshooting

| Issue | Solution |
| --- | --- |
| Tkinter not found | Install `python3-tk` (Debian/Ubuntu) or equivalent for your distribution. |
| No ITSs detected | Reduce the terminal margin (e.g., `1000`) or the minimum repeat count (e.g., `2`). |
| Too many ITSs detected | Increase the min ITS length or the min repeats value. |
| Fuzzy patterns capture non-telomeric regions | Increase min repeats; fuzzy motifs are intentionally permissive. |
| Suggestion returns no match | The genome may be too fragmented or contain no clear terminal arrays. Try a preset manually. |
| Invalid motif error | Use only A, C, G, T, N in custom motifs. |
| Analysis is slow on large genomes | Expected for genomes with thousands of scaffolds. Consider running on a subset first. |
| Out of memory | Extremely large genomes may require more RAM. The tool holds all sequences in memory. |
| Report contains no movements | Verify that the input genome actually has internal telomeres; run Dry Run with reduced margin to inspect. |

### Interpreting a Suspicious ITS

If a reported ITS looks like a low-complexity region rather than a real telomere:

- Check the pattern column in the report — fuzzy patterns (`*_fuzzy`) are more permissive.
- Increase min repeats (e.g., to 4 or 5) to require a longer tandem array.
- Increase min ITS length to exclude short spurious matches.

---

## Frequently Asked Questions

**Q: Does telomereFinder work for any organism?**

Yes, as long as the organism's telomeres are tandem repeats of a short motif. Presets are provided for plants, vertebrates, insects, nematodes, ciliates, some fungi, and trypanosomes. For other organisms, choose Custom and enter the motifs manually. Organisms whose telomeres are not tandem repeats (e.g., Drosophila, Saccharomyces) are not supported.

**Q: Does telomereFinder split scaffolds at internal telomeres?**

No. It relocates ITSs to the scaffold ends. If the biologically correct action is to split a misjoined scaffold, that must be done separately (e.g., with a scaffold-splitting tool). The output is a reorganized FASTA, not a split assembly.

**Q: What happens to the sequence between ITSs?**

It is preserved exactly as in the input, in its original order and orientation. Only the ITS blocks are moved.

**Q: Can two ITSs be relocated to the same end?**

Yes. Multiple forward ITSs are concatenated at the start in ascending order of their original position, and multiple reverse ITSs are concatenated at the end.

**Q: Does the tool modify scaffolds with no ITSs?**

No. Scaffolds without internal telomeres pass through unchanged and keep their original header.

**Q: Why do some scaffolds have both forward and reverse ITSs?**

This is typical of misjoined scaffolds where two chromosome arms of opposite orientation were merged. **telomereFinder** will move forward ITSs to the start and reverse ITSs to the end, producing a scaffold that is telomere-capped on both sides.

**Q: Can I run telomereFinder without the GUI?**

The current version is GUI-only. The core logic is organized in a `TelomereAnalyzer` class that can be imported and scripted if needed.

**Q: How long does the analysis take?**

For a typical genome with thousands of scaffolds of a few kb each, the analysis completes in seconds to a couple of minutes, depending on repeat count and margin settings.

**Q: How does Suggest Motifs work?**

It scans the first and last 2000 bp of every scaffold, counts tandem arrays for each preset's motifs, and picks the preset with the highest count. See the Motif Suggestion section for details and caveats.

**Q: What if my organism uses a telomere motif not in the preset list?**

Select Custom and enter the forward and reverse motifs manually (comma-separated, ACGTN only). Fuzzy variants are generated automatically when the checkbox is enabled.

**Q: Does telomereFinder require external software like HMMER or gffread?**

No. Only Python and Tkinter are required. **telomereFinder** works directly on FASTA sequences.

---

## Citation

If you use **telomereFinder** in your research, please cite:

> **telomereFinder v1.0: A universal GUI tool for detecting and relocating internal telomeric sequences in genome assemblies.**

---

## License

**telomereFinder** is distributed under the **MIT License**.

---

## Version History

| Version | Date | Changes |
| --- | --- | --- |
| v1.0 | Sep 2026 | Initial public release with multi-organism presets, fuzzy variant generation, and motif suggestion |

---

*telomereFinder v1.0 — Making telomere-aware assembly curation accessible to everyone.*
