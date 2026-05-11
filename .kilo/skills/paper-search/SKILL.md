---
name: paper-search
description: "Search for academic papers across Semantic Scholar, arXiv, and CrossRef APIs. Returns unified results with title, authors, year, abstract, DOI, venue, and citation counts. Integrates with Zotero MCP tools for adding found papers to a Zotero library and generating BibTeX entries. Use when the user asks to find papers, search for related work, look up a DOI, or discover references on a topic."
---

# Paper Search

Search external academic APIs for papers. Provides a unified interface across Semantic Scholar, arXiv, and CrossRef with optional Zotero integration.

## Workflow

### 1. Search for Papers

Run the search script from the skill's `scripts/` directory:

```bash
python3 scripts/search_papers.py --query "topic" --source semantic_scholar --limit 10 --format compact
```

Available sources:
- `semantic_scholar` — Default. Best for comprehensive search with citation counts.
- `arxiv` — Best for preprints and recent unpublished work.
- `crossref` — Best for published works and DOI-based metadata.
- `all` — Query all three sources (slower, results combined).

Output formats:
- `json` — Full JSON output (default). Good for programmatic use.
- `compact` — Human-readable summary with title, authors, year, venue, citations, and truncated abstract.

### 2. DOI Lookup

Look up a specific paper by DOI:

```bash
python3 scripts/search_papers.py --doi "10.1145/1234567.1234568" --format compact
```

### 3. Download PDFs

Download open-access PDFs directly from search results:

```bash
python3 scripts/search_papers.py --query "wavelet tree" --source arxiv --limit 3 --download ~/papers
```

- arXiv papers always have PDFs available.
- Semantic Scholar provides `openAccessPdf` URLs when available.
- CrossRef may provide PDF links via publisher APIs.

The `--download` flag adds a `downloaded_path` field to each result in JSON output.

### 4. Add to Zotero

**Option A: Via DOI/URL (metadata only)**

After finding relevant papers, add them to Zotero using the Zotero MCP tools:

- `zotero_add_by_doi` — Preferred when DOI is available. Fetches full metadata from CrossRef.
- `zotero_add_by_url` — Use for arXiv papers or when only a URL is available.

**Option B: Via downloaded PDF (metadata + attachment)**

Download the PDF first, then add to Zotero with the PDF file:

```bash
# Step 1: Download PDFs and get paths in JSON
python3 scripts/search_papers.py --doi "10.1007/978-3-540-73420-8_13" --download ~/papers --format json

# Step 2: Use zotero_zotero_add_from_file with the downloaded_path
```

The agent should call `zotero_zotero_add_from_file` with the `downloaded_path` from the JSON output. This attaches the PDF to the Zotero item and attempts DOI-based metadata extraction.

**Option C: Download + Zotero in one step**

Use `--zotero` to download PDFs with paths formatted for easy Zotero import:

```bash
python3 scripts/search_papers.py -q "succinct data structures" -s arxiv -n 3 --zotero --download ~/papers
```

After adding papers, update the semantic search database:

```
zotero_update_search_database
```

### 5. Generate BibTeX

For papers already in Zotero, use `zotero_get_item_metadata` with `format: "bibtex"` to get BibTeX entries. Alternatively, use `zotero_fetch` for full metadata.

For papers NOT in Zotero, BibTeX can be constructed from the search results' JSON fields (`authors`, `year`, `title`, `venue`, `doi`).

## Guidance

- Start with `semantic_scholar` for general queries — it has the broadest coverage and citation data.
- Use `arxiv` when looking for very recent work or preprints in CS/ML/physics.
- Use `crossref` for DOI lookups or when Semantic Scholar returns no results.
- When using `--source all`, results may contain duplicates (same paper from different sources). Deduplicate by DOI or title similarity.
- Citation counts are approximate and may differ across sources.
- arXiv results return the arXiv ID (e.g., `2301.12345`) which can be used with `zotero_add_by_url` via `https://arxiv.org/abs/2301.12345`.

## API Quirks

- **arXiv `atom:id` is NOT a DOI** — it contains an arXiv URL like `http://arxiv.org/abs/2301.12345`. Store the extracted ID in `arxiv_id` only; set `doi` to `None` for arXiv results. Writing the arXiv URL into `doi` produces invalid DOI metadata downstream (e.g., Zotero import).
- **CrossRef `select` must include `link`** — the `link` field is needed for `pdf_url` extraction. If omitted from `select`, the API won't return link metadata and `pdf_url` will silently be empty for all CrossRef results.
