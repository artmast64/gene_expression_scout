# Step 2 Functions

import time
from io import StringIO

import requests
import pandas as pd


def get_domain(taxon_id):
    """
    Retrieves the NCBI Taxon ID's domain using the NCBI API.

    Queries the NCBI taxonomy database via HTTP requests to retrieve the 
    lineage data for a specific ID in JSON format.
    """
    # NCBI EFetch endpoint for taxonomy
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{taxon_id}/dataset_report"
    # Define query parameters (asking for JSON format)
    
    print(f"Fetching domain type for taxon {taxon_id}...")
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an error for bad status codes (4xx, 5xx)
        data = response.json()
        
        # Extract the domain ID from the nested JSON structure
        domain_id = data["reports"][0]["taxonomy"]["classification"]["domain"]["id"]

        print(f"Taxon {taxon_id} has domain ID {domain_id}.")
        return domain_id
    except requests.exceptions.RequestException as e:
        return f"Network or API Error: {e}"
    except ValueError:
        return "Error parsing JSON response or Taxon ID not found."
    

def fetch_go_annotations(go_ids, taxon_id, original_genes, domain_id, aspect=None):
    """
    Fetch gene annotations for multiple GO terms in a single QuickGO
    downloadSearch request.

    Parameters:
    go_ids : list[str]
        e.g. ["GO:0006915", "GO:0008219", "GO:0007049", ...]
    taxon_id : int
        NCBI taxonomy ID
    aspect : str or None
        Restrict to "molecular_function", "biological_process",
        or "cellular_component". None = no restriction.
    """
    base_url = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch"

    params = {
        "goId": ",".join(go_ids),
        "taxonId": taxon_id,
        "goUsage": "exact",   # restrict to the exact GO terms only, no descendants
        "geneProductType": "protein",   # filter out non-genes/non-proteins
        "includeFields": "goName",
        "selectedFields": "geneProductId,symbol,qualifier,goId,goName,goAspect,evidenceCode,goEvidence,reference,withFrom,taxonId,assignedBy,extensions,date",
    }
    if aspect:
        params["aspect"] = aspect

    headers = {"Accept": "text/tsv"}

    print("Downloading GO terms...")
    response = requests.get(base_url, params=params, headers=headers, timeout=120)

    if response.status_code != 200:
        print("Request URL:", response.url)
        print("Response body:", response.text)
        response.raise_for_status()

    df = pd.read_csv(StringIO(response.text), sep="\t")

    records = []

    for _, row in df.iterrows():
        symbol = row.get("SYMBOL") or row.get("GENE PRODUCT ID", "")
        gene_name = symbol.split("_")[0] if "_" in symbol else symbol

        if gene_name.lower() in original_genes:
            continue

        # If bacteria, skip if first character is uppercase — likely a protein, not a gene
        if domain_id == 2:  # 2 = "Bacteria"
            if not gene_name or not gene_name[0].islower():
                continue

        records.append({
            "source_go_id":  row.get("GO TERM", ""),
            "go_term":       row.get("GO NAME", ""),
            "go_category":   row.get("GO ASPECT", ""),
            "new_gene":      gene_name,
            "uniprot_id":    row.get("GENE PRODUCT ID", ""),
            "evidence_code": row.get("GO EVIDENCE CODE", ""),
        })
    return records


def get_genes_for_go_term(go_ids: list[str] | str, original_genes: set[str], taxon_id: int, domain_id: int, batch_size: int = 5) -> list[dict]:
    """
    Query QuickGO for all genes annotated with any of the given
    GO terms. Accepts either a single GO ID string or a list of GO IDs,
    batching requests to avoid URL length limits.
    Excludes genes already in the original seed list.
    """
    if isinstance(go_ids, str):
        go_ids = [go_ids]

    batches = [go_ids[i:i + batch_size]
               for i in range(0, len(go_ids), batch_size)]

    records = []

    for batch in batches:
        print(f"Batch: {batch}")
        #batch_records = fetch_annotation_batch(batch, original_genes, domain_id, max_pages)
        batch_records = fetch_go_annotations(batch, taxon_id, original_genes, domain_id)

        # If the batch failed with a 500, fall back to one GO ID at a time
        if batch_records is None:
            print(f"  [!] Batch failed, retrying individually: {batch}")
            for go_id in batch:
                print(f"Trying: {go_id}")
                #single_records = fetch_annotation_batch([go_id], original_genes, domain_id, max_pages)
                single_records = fetch_go_annotations(batch, taxon_id, original_genes, domain_id)
                if single_records is None:
                    print(f"  [!] Skipping {go_id} — server error on individual query")
                else:
                    records.extend(single_records)
                time.sleep(0.5)
        else:
            print("  Batch downloaded successfully")
            records.extend(batch_records)

        time.sleep(1.0) # <-- This was 0.5

    return records


def get_genes_for_go_term_set(
    go_ids: list[str],
    original_genes: set[str],
    taxon_id: int,
    domain_id: int,
    batch_size: int = 5,
    max_pages: int = 5,
) -> set[str]:
    """
    Retrieves genes annotated with ALL GO terms in the given set,
    by fetching genes for each GO term individually and intersecting
    the results.
    """
    if isinstance(go_ids, str):
        go_ids = [go_ids]

    go_id_set = set(go_ids)
    records = get_genes_for_go_term(go_ids, original_genes, taxon_id, domain_id, batch_size)

    # Build gene → {go_ids seen} map from records
    gene_go_hits: dict[str, set[str]] = {}
    for record in records:
        gene      = record["new_gene"]
        go_id_hit = record["source_go_id"]
        if gene not in gene_go_hits:
            gene_go_hits[gene] = set()
        gene_go_hits[gene].add(go_id_hit)

    # Keep only genes annotated with every GO term in the set
    return {
        gene for gene, hits in gene_go_hits.items()
        if go_id_set.issubset(hits)
    }


def expand_genes_from_go_terms(go_df: pd.DataFrame, taxon_id: int, go_term_grouping: str, domain_id: int, batch_size: int = 5) -> pd.DataFrame:
    """
    Takes the GO term DataFrame from Step 1 and retrieves all additional
    E. coli genes associated with each unique GO term.

    GO term grouping:
      none — for each unique GO term, retrieve all genes annotated
                   with it
      all — for each original gene, retrieve only genes annotated
                   with the exact same combination of GO terms as that gene
      categories — for each original gene, retrieve only genes annotated
                   with the same exact combinations of GO term category as that gene
    """
    # Exclude original genes from results
    original_genes = set(go_df["gene"].str.lower().unique())
    print(f"\noriginal_genes: {original_genes}")

    all_records    = []

    if go_term_grouping == "none":
        unique_go_ids = go_df["go_id"].dropna().unique().tolist()
        print(f"\nGo term grouping: {go_term_grouping}")
        print(f"Querying {len(unique_go_ids)} unique GO terms...\n")
        records = get_genes_for_go_term(unique_go_ids, original_genes, taxon_id, domain_id, batch_size)
        all_records.extend(records)
    
    elif go_term_grouping == "all" or go_term_grouping == "categories":
        # Group GO terms by original gene
        if go_term_grouping == "all":
            gene_go_groups = (
                go_df.groupby("gene")["go_id"]
                .apply(lambda x: sorted(x.dropna().unique().tolist()))
                .to_dict()
            )
        if go_term_grouping == "categories":
            gene_go_groups = (
                go_df.groupby(["gene", "go_category"])["go_id"]
                .apply(lambda x: sorted(x.dropna().unique().tolist()))
                .to_dict()
            )

        # Ensure every value is a list
        gene_go_groups = {
            gene: ([gene_go_ids] if isinstance(gene_go_ids, str) else gene_go_ids)
            for gene, gene_go_ids in gene_go_groups.items()
        }

        print(f"\nGo term grouping: {go_term_grouping}")
        print(f"Querying {len(gene_go_groups)} GO term groups...\n")

        # Convert to list to avoid iterator exhaustion issues,
        # and rename loop variable to avoid shadowing
        for gene, gene_go_ids in list(gene_go_groups.items()):
            print(f"  Finding genes sharing all {len(gene_go_ids)} GO terms "
                f"of {gene}: {gene_go_ids}")
            try:
                matched_genes = get_genes_for_go_term_set(gene_go_ids, original_genes, taxon_id, domain_id)
                print(f"    → {len(matched_genes)} genes share the exact GO term set")

                if go_term_grouping == "all":
                    for new_gene in matched_genes:
                        all_records.append({
                            "source_gene":       gene,
                            "source_go_ids":     "|".join(gene_go_ids),
                            "new_gene":          new_gene,
                            "go_term_hit_count": len(gene_go_ids),
                        })
                if go_term_grouping == "categories":
                    for new_gene in matched_genes:
                        all_records.append({
                            "source_gene":       gene[0],
                            "go_category":       gene[1],
                            "source_go_ids":     "|".join(gene_go_ids),
                            "new_gene":          new_gene,
                            "go_term_hit_count": len(gene_go_ids),
                        })
            except requests.HTTPError as e:
                print(f"    [!] HTTP error for {gene}: {e}")
            time.sleep(0.5)

    if not all_records:
        print("No new genes found.")
        return pd.DataFrame()

    df = pd.DataFrame(all_records)

    # Deduplicate: one row per (new_gene, source_go_id) pair
    df = df.drop_duplicates(subset=["new_gene",
                                     "source_go_id" if go_term_grouping == "none"
                                     else "source_gene"])

    # Summary: how many GO terms link to each new gene (a relevance signal)
    if go_term_grouping == "none":
        df["go_term_hit_count"] = df.groupby("new_gene")["source_go_id"].transform("count")

    df = df.sort_values(["go_term_hit_count", "new_gene"], ascending=[False, True])
    return df