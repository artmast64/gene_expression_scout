# Step 1 Functions

import time

import requests
import pandas as pd
from goatools.obo_parser import GODag


def get_go_terms_for_gene(gene_names: list[str], taxon_id: int, category_map: dict) -> list[dict]:
    """
    Query UniProt for multiple genes in a single API call and return
    their GO annotations. Filters to reviewed (Swiss-Prot) entries.
    Returns only the first UniProt hit per gene.
    """
    # Combine all gene names into a single OR query
    gene_names_lower = {g.lower(): g for g in gene_names}  # lowercase → original
    gene_queries = {f"gene:{g}" for g in gene_names_lower.values()}  # add gene_exact key to all entries
    gene_clause = f"({' OR '.join(gene_queries)})"

    # Build query componentes dynamically
    query_parts = [gene_clause]
    query_parts.append(f"organism_id:{taxon_id}")
    #query_parts.append(f"reviewed:true")

    # Combine parts with AND
    full_query = " AND ".join(query_parts)

    uniprot_search_url = "https://rest.uniprot.org/uniprotkb/search"

    params = {
        "query": full_query,
        "fields": "gene_names,go_id,go",   # request GO IDs + full GO annotation field
        "format": "json",
        "size": len(gene_names) * 3,  # allow headroom for multiple hits per gene
    }

    # Call UniProt API
    response = requests.get(uniprot_search_url, params=params)
    response.raise_for_status()
    results = response.json().get("results", [])

    if not results:
        print(f"  [!] No UniProt entry found for query: {full_query}")
        return []
    
    seen_genes = set()
    go_records = []

    for entry in results:
        accession = entry.get("primaryAccession", "N/A")

        # Collect all name variants from this entry: primary, synonyms, locus names, ORF names
        entry_gene_names = []
        for g in entry.get("genes", []):
            if "geneName" in g:
                entry_gene_names.append(g["geneName"]["value"])
            for synonym in g.get("synonyms", []):
                entry_gene_names.append(synonym["value"])
            for locus in g.get("orderedLocusNames", []):
                entry_gene_names.append(locus["value"])
            for orf in g.get("orfNames", []):
                entry_gene_names.append(orf["value"])

        # Find which queried gene name this entry matches
        matched_gene = next(
            (gene_names_lower[n.lower()] for n in entry_gene_names if n.lower() in gene_names_lower),
            None
        )

        # Skip if we've already processed this gene
        if matched_gene is None or matched_gene.lower() in seen_genes:
            continue
        seen_genes.add(matched_gene.lower())
        print(f"  [{matched_gene}] → UniProt accession: {accession}")

        for xref in entry.get("uniProtKBCrossReferences", []):
            if xref.get("database") != "GO":
                continue
            go_id = xref.get("id")
            properties = {p["key"]: p["value"] for p in xref.get("properties", [])}
            raw_term = properties.get("GoTerm", "")
            category_code = raw_term[0] if raw_term else ""
            term_name = raw_term[2:] if len(raw_term) > 2 else raw_term
            category_label = category_map.get(category_code, "Unknown")

            go_records.append({
                "gene":        matched_gene,
                "uniprot_id":  accession,
                "go_id":       go_id,
                "go_term":     term_name,
                "go_category": category_label,
            })

    for g in gene_names:
        if g.lower() not in seen_genes:
            print(f"  [!] No UniProt entry found for gene: {g}")

    return go_records


def filter_obsolete_terms(go_ids, batch_size):
    """
    Splits a list of GO IDs into active and obsolete sets using the
    isObsolete field from QuickGO's ontology endpoint.
    """
    active, obsolete = [], []

    # Retrieve batches
    for i in range(0, len(go_ids), batch_size):
        batch = go_ids[i:i + batch_size]
        url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{','.join(batch)}"
        r = requests.get(url, headers={"Accept": "application/json"})
        r.raise_for_status()
        results = r.json()["results"]

        for term in results:
            if term.get("isObsolete"):
                obsolete.append(term["id"])
            else:
                active.append(term["id"])

    return active, obsolete


def fetch_all_go_terms(gene_list: list[str], taxon_id: int, category_map: dict) -> pd.DataFrame:
    """
    Query UniProt for GO annotations of all provided genes.
    """
    all_records = []

    print(f"Fetching GO terms for genes: {gene_list}")
    records = get_go_terms_for_gene(gene_list, taxon_id, category_map)

    # Add gene records to master records list
    all_records.extend(records)
    time.sleep(0.5)

    return all_records


def get_go_depths(go_ids: list[str], godag) -> dict[str, int]:
    """
    Returns level and depth for each GO ID using GOATOOLS.
    - level: shortest path from root
    - depth: longest path from root (recommended for GO specificity filtering)
    """
    results = {}
    for go_id in go_ids:
        term = godag.get(go_id)
        if term is None or term.is_obsolete:
            print(f"  [!] GO term not found or obsolete: {go_id}")
            continue
        results[go_id] = {
            "level": term.level,
            "depth": term.depth,
        }
    return results

def filter_records_by_depth(records: list[dict], godag, min_depth: int = 0) -> list[dict]:
    """
    Given a flat list of GO annotation records (each with a "go_id" key),
    fetches depths for all unique GO IDs in one batched pass,
    then filters out records below min_depth.
    """
    # Collect unique GO IDs across all records
    unique_go_ids = list({r["go_id"] for r in records if r.get("go_id")})
    print(f"\nFetching depths for {len(unique_go_ids)} unique GO terms...")

    depth_info = get_go_depths(unique_go_ids, godag)

    # Filter records
    filtered = []
    removed = []
    for record in records:
        go_id = record.get("go_id")
        info = depth_info.get(go_id)
        if info is None:
            continue
        if info["depth"] >= min_depth:
            filtered.append({**record, "go_depth": info["depth"], "go_level": info["level"]})
        else:
            removed.append((go_id, record.get("go_term", ""), info["depth"]))

    if removed:
        print(f"Removed {len(removed)} records below depth {min_depth}:")
        for go_id, term, depth in sorted(set(removed), key=lambda x: x[2]):
            print(f"  {go_id} ({term}): depth {depth}")

    print(f"Kept {len(filtered)}/{len(records)} records at depth >= {min_depth}")
    return filtered