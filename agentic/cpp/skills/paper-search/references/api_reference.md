# External Paper Search APIs

## Semantic Scholar

- **Base URL**: `https://api.semanticscholar.org/graph/v1/`
- **Rate limit**: 1 req/sec without API key, 10 req/sec with key
- **No auth required** for basic usage
- **Fields**: title, authors, year, abstract, externalIds (DOI, ArXiv), venue, citationCount
- **Best for**: Comprehensive academic search with citation counts

## arXiv

- **Base URL**: `http://export.arxiv.org/api/query`
- **Rate limit**: Be nice, ~3 sec between requests
- **No auth required**
- **Returns**: XML (Atom feed)
- **Best for**: Preprints, recent work not yet published

## CrossRef

- **Base URL**: `https://api.crossref.org/`
- **Rate limit**: 50 req/sec with polite pool (include `mailto` header)
- **No auth required**
- **Best for**: DOI lookup, published works, metadata enrichment

## Zotero Integration

After finding papers via external search, use Zotero MCP tools:

1. `zotero_add_by_doi` — Add paper by DOI (fetches metadata from CrossRef)
2. `zotero_add_by_url` — Add paper by URL (arXiv, DOI URLs)
3. `zotero_update_search_database` — Update semantic search index after adding
