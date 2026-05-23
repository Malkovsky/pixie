# AGENTS.md - Academic Materials

## Overview

This directory contains academic materials related to the Pixie succinct data structures library. It serves as a centralized repository for published papers, presentations, technical reports, bibliography, and other scholarly resources that document the design, implementation, and evaluation of Pixie.

Materials are authored in [Quarto](https://quarto.org) for rich rendering (math, tables, code, citations). Research notes may use plain markdown.

## Directory Structure

```
academic/
├── _quarto.yml              # Book project (reports, presentations, notes)
├── index.qmd                # Project landing page
├── scripts/
│   └── init.sh              # Scaffolds new papers (supports acm/ieee styles)
├── papers/                  # Research papers (standalone Quarto projects)
│   └── <paper-name>/        # Each paper in its own directory
│       ├── _quarto.yml      # Paper's own Quarto project (type: default)
│       ├── _extensions/     # Quarto extensions (not in git, installed by init.sh)
│       ├── paper.qmd        # Paper source (acm-pdf or ieee-pdf format)
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

Papers are standalone because Quarto book projects cannot resolve journal extensions (like `quarto-journals/acm`). Each paper directory has its own `_quarto.yml` and `_extensions/` (installed locally, not tracked in git).

### Journal Extensions

Journal formatting (ACM, IEEE, etc.) is NOT built-in to Quarto. It requires per-paper extensions, installed in `_extensions/`.

Extensions are **not tracked in git** (they are in `.gitignore`). Use the init script to scaffold a new paper with the extension pre-installed:

```bash
cd academic/
./scripts/init.sh <paper-name>                      # ACM (default)
./scripts/init.sh --style ieee <paper-name>         # IEEE
```

The script copies `papers/sample-paper/` for shared files (_quarto.yml, references.bib, .gitignore), then replaces paper.qmd with the style-specific template from `scripts/templates/` and installs the matching extension.

Supported styles and their extensions:
| Style | Extension | Format |
|-------|-----------|--------|
| `acm` | `quarto-journals/acm` | `acm-pdf` |
| `ieee` | `dfolio/quarto-ieee` | `ieee-pdf` |

To manually install an extension into an existing paper directory:
```bash
cd papers/<paper-name>/
quarto add quarto-journals/acm --no-prompt
quarto add dfolio/quarto-ieee --no-prompt
```

## Rendering

### Papers

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

### Paper Templates

Each paper lives in its own subdirectory under `papers/` with:

- `_quarto.yml` -- `project: type: default` (standalone, not book)
- `_extensions/` -- Quarto extensions (installed by init script, not in git)
- `paper.qmd` -- Main source file (format depends on style: `acm-pdf` or `ieee-pdf`)
- `references.bib` -- Paper-specific BibTeX

Style templates live in `scripts/templates/<style>.qmd`. To add a new journal style, create a new template file there and add the extension mapping to `init.sh`.

To create a new paper: run `./scripts/init.sh [--style <style>] <paper-name>` from the `academic/` directory.

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
- New papers should be created with `./scripts/init.sh [--style <style>] <paper-name>` (from `academic/`)
- Each paper gets its own directory with `_quarto.yml`, `_extensions/`, `paper.qmd`, and `references.bib`
- Use `@AuthorYear` citation syntax in `.qmd` files
- LaTeX code blocks should use `{cpp}` or `{latex}` language tags
- Do not commit large binary files (>10MB) without explicit approval; consider using Git LFS
- Render papers from their own directory, NOT from the `academic/` book project root
- Do NOT put `@` symbols in BibTeX comments — BibTeX tries to parse them as entries
