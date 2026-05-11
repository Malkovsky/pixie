#!/usr/bin/env python3
"""Search external APIs for academic papers.

Sources: Semantic Scholar, arXiv, CrossRef.
Outputs unified JSON to stdout.

Usage:
    python3 search_papers.py --query "wavelet tree succinct" --source semantic_scholar --limit 10
    python3 search_papers.py --query "succinct data structures" --source arxiv --limit 5
    python3 search_papers.py --doi "10.1145/123" --source crossref
    python3 search_papers.py --query "rank select" --source all --limit 5
    python3 search_papers.py --query "wavelet tree" --source arxiv --limit 1 --download ~/papers
    python3 search_papers.py --doi "10.1007/978-3-540-73420-8_13" --download ~/papers --zotero
"""

import argparse
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


def _get(url: str, headers: dict[str, str] | None = None, timeout: int = 30,
         retries: int = 2) -> dict:
    for attempt in range(retries + 1):
        req = urllib.request.Request(url, headers=headers or {})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 429 and attempt < retries:
                wait = 2 ** attempt
                print(f"Rate limited, retrying in {wait}s...", file=sys.stderr)
                time.sleep(wait)
                continue
            print(f"HTTP {e.code}: {e.reason} for {url}", file=sys.stderr)
            return {}
        except urllib.error.URLError as e:
            print(f"URL error: {e.reason} for {url}", file=sys.stderr)
            return {}
    return {}


def search_semantic_scholar(query: str, limit: int = 10) -> list[dict[str, Any]]:
    """Search Semantic Scholar API."""
    params = urllib.parse.urlencode({
        "query": query,
        "limit": limit,
        "fields": "title,authors,year,abstract,externalIds,venue,publicationDate,citationCount,url,openAccessPdf",
    })
    url = f"https://api.semanticscholar.org/graph/v1/paper/search?{params}"
    data = _get(url, headers={"Accept": "application/json"})
    results = []
    for paper in data.get("data", []):
        ext_ids = paper.get("externalIds") or {}
        pdf_info = paper.get("openAccessPdf") or {}
        results.append({
            "source": "semantic_scholar",
            "title": paper.get("title", ""),
            "authors": [a.get("name", "") for a in paper.get("authors", [])],
            "year": paper.get("year"),
            "abstract": paper.get("abstract", ""),
            "doi": ext_ids.get("DOI"),
            "arxiv_id": ext_ids.get("ArXiv"),
            "venue": paper.get("venue", ""),
            "citation_count": paper.get("citationCount"),
            "url": paper.get("url", ""),
            "pdf_url": pdf_info.get("url"),
        })
    return results


def search_arxiv(query: str, limit: int = 10) -> list[dict[str, Any]]:
    """Search arXiv API."""
    words = query.split()
    if len(words) == 1:
        search_term = f"all:{query}"
    elif len(words) == 2:
        # Phrase search for 2-word queries
        search_term = f'all:"{query}"'
    else:
        # Use OR of phrase and individual terms for 3+ words
        # This catches exact phrase matches AND papers with all terms
        phrase = f'all:"{query}"'
        and_terms = " AND ".join(f"all:{w}" for w in words)
        search_term = f"({phrase}) OR ({and_terms})"
    params = urllib.parse.urlencode({
        "search_query": search_term,
        "start": 0,
        "max_results": limit,
    })
    url = f"http://export.arxiv.org/api/query?{params}"
    req = urllib.request.Request(url)
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            xml_data = resp.read().decode()
    except (urllib.error.URLError, urllib.error.HTTPError) as e:
        print(f"arXiv API error: {e}", file=sys.stderr)
        return []

    import xml.etree.ElementTree as ET
    root = ET.fromstring(xml_data)
    ns = {"atom": "http://www.w3.org/2005/Atom"}
    results = []
    for entry in root.findall("atom:entry", ns):
        title = entry.findtext("atom:title", "", ns).strip().replace("\n", " ")
        abstract = entry.findtext("atom:summary", "", ns).strip().replace("\n", " ")
        authors = [a.findtext("atom:name", "", ns) for a in entry.findall("atom:author", ns)]
        published = entry.findtext("atom:published", "", ns)
        year = int(published[:4]) if published else None
        arxiv_id = ""
        for link in entry.findall("atom:link", ns):
            href = link.get("href", "")
            if "arxiv.org/abs/" in href:
                arxiv_id = href.split("/abs/")[-1]
                break
        results.append({
            "source": "arxiv",
            "title": title,
            "authors": authors,
            "year": year,
            "abstract": abstract,
            "doi": None,
            "arxiv_id": arxiv_id,
            "venue": "arXiv",
            "citation_count": None,
            "url": f"https://arxiv.org/abs/{arxiv_id}" if arxiv_id else "",
            "pdf_url": f"https://arxiv.org/pdf/{arxiv_id}" if arxiv_id else None,
        })
    return results


def search_crossref(query: str, limit: int = 10) -> list[dict[str, Any]]:
    """Search CrossRef API."""
    params = urllib.parse.urlencode({
        "query": query,
        "rows": limit,
        "select": "DOI,title,author,published-print,abstract,container-title,is-referenced-by-count,URL,type,link",
    })
    url = f"https://api.crossref.org/works?{params}"
    data = _get(url, headers={"Accept": "application/json"})
    results = []
    for item in data.get("message", {}).get("items", []):
        title_list = item.get("title", [])
        title = title_list[0] if title_list else ""
        authors = []
        for a in item.get("author", []):
            name = f"{a.get('given', '')} {a.get('family', '')}".strip()
            if name:
                authors.append(name)
        pub_date = item.get("published-print", {}).get("date-parts", [[None]])
        year = pub_date[0][0] if pub_date and pub_date[0] else None
        venue_list = item.get("container-title", [])
        venue = venue_list[0] if venue_list else ""
        pdf_url = None
        for link in item.get("link", []):
            if "pdf" in link.get("content-type", ""):
                pdf_url = link.get("URL")
                break
        results.append({
            "source": "crossref",
            "title": title,
            "authors": authors,
            "year": year,
            "abstract": item.get("abstract", ""),
            "doi": item.get("DOI"),
            "arxiv_id": None,
            "venue": venue,
            "citation_count": item.get("is-referenced-by-count"),
            "url": item.get("URL", ""),
            "pdf_url": pdf_url,
        })
    return results


def lookup_doi(doi: str) -> dict[str, Any] | None:
    """Look up a single paper by DOI via CrossRef."""
    url = f"https://api.crossref.org/works/{urllib.parse.quote(doi, safe='')}"
    data = _get(url)
    item = data.get("message")
    if not item:
        return None
    title_list = item.get("title", [])
    title = title_list[0] if title_list else ""
    authors = []
    for a in item.get("author", []):
        name = f"{a.get('given', '')} {a.get('family', '')}".strip()
        if name:
            authors.append(name)
    pub_date = item.get("published-print", {}).get("date-parts", [[None]])
    year = pub_date[0][0] if pub_date and pub_date[0] else None
    venue_list = item.get("container-title", [])
    venue = venue_list[0] if venue_list else ""
    pdf_url = None
    for link in item.get("link", []):
        if "pdf" in link.get("content-type", ""):
            pdf_url = link.get("URL")
            break
    return {
        "source": "crossref",
        "title": title,
        "authors": authors,
        "year": year,
        "abstract": item.get("abstract", ""),
        "doi": item.get("DOI"),
        "arxiv_id": None,
        "venue": venue,
        "citation_count": item.get("is-referenced-by-count"),
        "url": item.get("URL", ""),
        "pdf_url": pdf_url,
    }


SOURCES = {
    "semantic_scholar": search_semantic_scholar,
    "arxiv": search_arxiv,
    "crossref": search_crossref,
}


def _sanitize_filename(title: str) -> str:
    """Generate a clean filename from paper title."""
    name = re.sub(r'[^\w\s-]', '', title.lower())
    name = re.sub(r'[\s]+', '_', name.strip())
    return name[:80]


def download_pdf(url: str, output_dir: str, paper: dict[str, Any]) -> str | None:
    """Download a PDF and return the local path."""
    filename = _sanitize_filename(paper.get("title", "paper")) + ".pdf"
    output_path = Path(output_dir) / filename
    output_path.parent.mkdir(parents=True, exist_ok=True)

    req = urllib.request.Request(url, headers={
        "User-Agent": "Mozilla/5.0 (academic paper-search script)"
    })
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            content_type = resp.headers.get("Content-Type", "")
            if "pdf" not in content_type and "octet-stream" not in content_type:
                print(f"Warning: unexpected content type '{content_type}' for {url}",
                      file=sys.stderr)
            with open(output_path, "wb") as f:
                f.write(resp.read())
        print(f"Downloaded: {output_path}", file=sys.stderr)
        return str(output_path)
    except (urllib.error.URLError, urllib.error.HTTPError) as e:
        print(f"Download failed for {url}: {e}", file=sys.stderr)
        return None


def main():
    parser = argparse.ArgumentParser(description="Search for academic papers")
    parser.add_argument("--query", "-q", help="Search query")
    parser.add_argument("--doi", help="Look up a specific DOI")
    parser.add_argument("--source", "-s", default="semantic_scholar",
                        choices=["semantic_scholar", "arxiv", "crossref", "all"],
                        help="API source (default: semantic_scholar)")
    parser.add_argument("--limit", "-n", type=int, default=10,
                        help="Max results per source (default: 10)")
    parser.add_argument("--format", "-f", default="json",
                        choices=["json", "compact"],
                        help="Output format (default: json)")
    parser.add_argument("--download", "-d", metavar="DIR",
                        help="Download PDFs to DIR (requires pdf_url in results)")
    parser.add_argument("--zotero", "-z", action="store_true",
                        help="Download PDFs and output paths for Zotero import (implies --download)")
    args = parser.parse_args()

    if not args.query and not args.doi:
        parser.error("Either --query or --doi is required")

    if args.zotero and not args.download:
        args.download = "."

    results = []
    if args.doi:
        paper = lookup_doi(args.doi)
        if paper:
            results.append(paper)
    elif args.source == "all":
        for name, func in SOURCES.items():
            try:
                results.extend(func(args.query, args.limit))
            except Exception as e:
                print(f"Error searching {name}: {e}", file=sys.stderr)
            time.sleep(1)
    else:
        results = SOURCES[args.source](args.query, args.limit)

    if args.download:
        for r in results:
            pdf_url = r.get("pdf_url")
            if pdf_url:
                path = download_pdf(pdf_url, args.download, r)
                r["downloaded_path"] = path
            else:
                r["downloaded_path"] = None

    if args.format == "json":
        print(json.dumps(results, indent=2))
    else:
        for i, r in enumerate(results, 1):
            authors = ", ".join(r["authors"][:3])
            if len(r["authors"]) > 3:
                authors += " et al."
            doi_str = f"  DOI: {r['doi']}" if r.get("doi") else ""
            arxiv_str = f"  arXiv: {r['arxiv_id']}" if r.get("arxiv_id") else ""
            cite_str = f"  Citations: {r['citation_count']}" if r.get("citation_count") else ""
            pdf_str = f"  PDF: {r['pdf_url']}" if r.get("pdf_url") else "  PDF: N/A"
            dl_str = ""
            if r.get("downloaded_path"):
                dl_str = f"  Downloaded: {r['downloaded_path']}"
            print(f"[{i}] {r['title']}")
            print(f"    {authors} ({r.get('year', '?')}) — {r.get('venue', '')}")
            print(f"    {r.get('url', '')}{doi_str}{arxiv_str}{cite_str}")
            print(f"    {pdf_str}{dl_str}")
            if r.get("abstract"):
                abstract = r["abstract"][:200]
                if len(r["abstract"]) > 200:
                    abstract += "..."
                print(f"    {abstract}")
            print()


if __name__ == "__main__":
    main()
