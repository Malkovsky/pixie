# AGENTS.md - Academic Materials

## Overview

This directory contains academic materials related to the Pixie succinct data structures library. It serves as a centralized repository for published papers, presentations, technical reports, bibliography, and other scholarly resources that document the design, implementation, and evaluation of Pixie.

Materials are authored in [Quarto](https://quarto.org) for rich rendering (math, tables, code, citations). Research notes may use plain markdown.

## Directory Structure

```
academic/
├── _quarto.yml              # Book project (reports, presentations, notes)
├── index.qmd                # Project landing page
├── _extensions/             # Shared Quarto extensions (ACM journal format)
│   └── quarto-journals/acm/ # ACM acm-pdf extension
├── papers/                  # Research papers (standalone Quarto projects)
│   └── <paper-name>/        # Each paper in its own directory
│       ├── _quarto.yml      # Paper's own Quarto project (type: default)
│       ├── _extensions -> ../../_extensions  # Symlink to shared extensions
│       ├── paper.qmd        # Paper source (ACM acm-pdf format)
│       └── references.bib   # Paper-specific bibliography
├── presentations/           # Conference talks, seminar slides
├── reports/                 # Technical reports, design documents
├── bibliography/
│   └── references.bib       # Master bibliography
├── notes/                   # Research notes (plain markdown OK)
└── figures/                 # Shared figures, diagrams, plots
```

## Architecture

### Two Quarto Project Types

The `academic/` directory contains **two kinds** of Quarto projects:

1. **Book project** (`academic/_quarto.yml`, `type: book`): Covers reports, presentations, and notes. Rendered from the `academic/` root.
2. **Standalone paper projects** (`papers/<name>/_quarto.yml`, `type: default`): Each paper is its own Quarto project. Rendered from the paper's directory.

Papers are standalone because Quarto book projects cannot resolve journal extensions (like `quarto-journals/acm`). Each paper directory has its own `_quarto.yml` and symlinks `_extensions/` to the shared extensions at `academic/_extensions/`.

### ACM Extension

The `acm-pdf` format is NOT built-in to Quarto. It requires the `quarto-journals/acm` extension, installed at `academic/_extensions/quarto-journals/acm/`.

**Installation** (one-time, already done):
```bash
cd academic/
quarto add quarto-journals/acm
```

To install for a new paper directory:
```bash
cd papers/<new-paper>/
ln -s ../../_extensions _extensions
```

## Rendering

### Papers (ACM format)

Each paper is rendered from its own directory:

```bash
cd papers/sample-paper/
quarto render                # Produces paper.pdf + paper.tex
quarto preview               # Live preview with reload
```

### Book (reports, presentations, notes)

```bash
cd academic/
quarto render                # Produces _book/ directory
quarto preview               # Live preview
```

## Conventions

### Paper Template (ACM Format)

Each paper lives in its own subdirectory under `papers/` with:

- `_quarto.yml` -- `project: type: default` (standalone, not book)
- `_extensions/` -- Symlink to `../../_extensions`
- `paper.qmd` -- Main source file using ACM `acm-pdf` format
- `references.bib` -- Paper-specific BibTeX

The paper's `paper.qmd` must include:

1. `author:` (not `authors:`) with nested `affiliation:` objects
2. `acm-metadata:` block with: `final`, `acmart-options`, `copyright-year`, `acm-year`, `copyright`, `doi`, `conference-*`, `price`, `isbn`, `ccs`, `keywords`, `abstract`
3. `format: acm-pdf:` with `keep-tex: true` and `biblio-style: ACM-Reference-Format`
4. `bibliography: references.bib`

To create a new paper: copy `papers/sample-paper/` directory and modify. Ensure the `_extensions` symlink is preserved.

### File Naming

- Use descriptive, lowercase filenames with hyphens as separators
- Include year prefix for chronological sorting: `2026-wavelet-tree-construction/`
- BibTeX citation keys: `AuthorYear` format (e.g., `Navarro2016`)

### BibTeX Entries

- Maintain a master bibliography file: `bibliography/references.bib`
- Use standard BibTeX entry types (`@article`, `@inproceedings`, `@book`, etc.)
- Include DOI and URL fields when available
- Do NOT use `@` inside BibTeX comments — BibTeX misparses them

## Key References

The following works are foundational to Pixie's design:

- **SPIDER** (Williams & Kurlin, 2023): Reference for BitVector implementation
- **pasta-toolbox**: Reference implementation for comparison benchmarks
- **sdsl-lite** (Gog et al., 2014): Predecessor succinct data structure library
- **Navarro (2016)**: "Compact Data Structures" -- comprehensive textbook
- **Jacobson (1989)**: Original PhD thesis on succinct graph representations

## Guidelines for AI Agents

- Do not modify compiled PDFs or binary files directly
- When adding new papers or references, update `bibliography/references.bib` accordingly
- Preserve existing citation keys and BibTeX formatting
- New papers should copy `papers/sample-paper/` as a template
- Each paper gets its own directory with `_quarto.yml`, `_extensions` symlink, `paper.qmd`, and `references.bib`
- Use `@AuthorYear` citation syntax in `.qmd` files
- LaTeX code blocks should use `{cpp}` or `{latex}` language tags
- Do not commit large binary files (>10MB) without explicit approval; consider using Git LFS
- Render papers from their own directory, NOT from the `academic/` book project root
- Do NOT put `@` symbols in BibTeX comments — BibTeX tries to parse them as entries
