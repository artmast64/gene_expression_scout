# Gene Expression Scout

from datetime import datetime
from collections import Counter
import json
from pprint import pprint

import pandas as pd
from goatools.obo_parser import GODag
import GEOparse
import numpy as np

import retrieve_go_terms
import retrieve_genes_from_go_terms
import retrieve_experiment_keywords
import retrieve_experiments_and_expression_levels

# Default settings
#go_term_grouping = "none" # none, categories, all
#batch_size = 50
#min_depth = 4
#max_series_return = 50

#model_name = "gemini-2.5-flash"
#api_key = "AQ.Ab8RN6J4ptcINnUhIf3jd_f5JpcpkjF6isW0tYKGFi0dz_YVcg"

def get_params():
    # GFF files are downloaded per-experiment based on the strain detected in GEO metadata.
    # They are cached locally by assembly accession to avoid repeat downloads.
    gff_cache_dir = "gff_cache"

    # Register your email with NCBI (required by their usage policy)
    ncbi_email = "bjohnsonhill@oakland.edu" # <-- TEMPORARY email !!!
    #ncbi_email = input("\n[!] NCBI's API requires an email address\nPlease input your email: ")

    #skip_steps_flag = 1

    # Map UniProt GO category codes to human-readable labels
    category_map = {
        "F": "Molecular Function",
        "P": "Biological Process",
        "C": "Cellular Component",
    }

    return gff_cache_dir, ncbi_email, category_map


def print_settings(species, condition, gene_list, model_name, go_term_grouping, batch_size, min_depth, max_series_return, warn_file_size_mb, drop_unmatched_genes):
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("\n--- Settings: ---")
    print(f"  Species name: {species}")
    print(f"  Factor: {condition}")
    print(f"  List of genes: {gene_list}")
    print(f"  LLM model name: {model_name}")
    print(f"  Gene Retrieval Method: {go_term_grouping}")
    print(f"  Batch Size: {batch_size}")
    print(f"  Minimum GO Term Depth: {min_depth}")
    print(f"  Maximum GEO series results per keyword: {max_series_return}")
    print(f"  Warning size for large supplementary files: {warn_file_size_mb} MB")
    print(f"  Drop unmatched genes from expression matrix: {drop_unmatched_genes}")
    print()


def identify_taxon(query: str, ncbi_email: str):
    """
    Identifies a NCBI Taxon ID for a given search query string.
    
    1. Looks up 'query' in `species_df` matching against scientific or common names.
    2. If not found, queries NCBI's Taxonomy API via E-utilities.
    
    Returns:
        str: Taxon ID (e.g., '562') if found/resolved, else None.
    """
    # --- Lookup table of common species ---
    # Replaced 562 with 83333 for E. coli
    species_df = pd.DataFrame({
        "Common Name": ["Mouse", "Human", "Rat", "Thale cress", "Zebrafish", "Fruit fly", "Baker's yeast", "Roundworm", "E. coli", "Cattle"],
        "Scientific Name": ["Mus musculus", "Homo sapiens", "Rattus norvegicus", "Arabidopsis thaliana", "Danio rerio", "Drosophilia melanogaster",
                            "Saccharomyces cerevisiae", "Caenorhabditis elegans", "Escherichia coli", "Bos taurus"],
        "Taxon ID": [10090, 9606, 10116, 3702, 7955, 7227, 4932, 6239, 562, 9913],
    })

    # --- Check local species_df for matches ---
    import re

    if not query or not isinstance(query, str):
        return None

    raw_query = query.strip()
    clean_query = re.sub(r'[^a-zA-Z0-9\s]', '', raw_query).lower()
    clean_query = " ".join(clean_query.split())  # collapse extra spaces

    for _, row in species_df.iterrows():
        com_name = re.sub(r'[^a-zA-Z0-9\s]', '', str(row["Common Name"])).lower()
        sci_name = re.sub(r'[^a-zA-Z0-9\s]', '', str(row["Scientific Name"])).lower()
        tax_id = str(row["Taxon ID"]).strip()

        if clean_query in (com_name, sci_name):
            print(f"[Taxon Lookup] Found local match in species_df: '{query}' -> TaxID {tax_id}")
            return tax_id
        
    # --- Fallback to NCBI Taxonomy API ---
    import requests

    print(f"[Taxon Lookup] '{query}' not found locally. Querying NCBI Taxonomy API...")
    
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "taxonomy",
        "term": query,
        "retmode": "json",
        "tool": "streamlit_transcriptomics_app",
        "email": ncbi_email
    }

    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        
        # Sanitize potential control characters from NCBI response
        clean_text = re.sub(r'[\x00-\x1F\x7F]', '', response.text)
        data = json.loads(clean_text)
        
        id_list = data.get("esearchresult", {}).get("idlist", [])
        
        if id_list:
            ncbi_tax_id = str(id_list[0])
            print(f"[Taxon Lookup] NCBI API resolved '{query}' -> TaxID {ncbi_tax_id}")
            return ncbi_tax_id
        else:
            print(f"[Taxon Lookup] Warning: NCBI could not find a Taxon ID for '{query}'.")
            return None

    except Exception as e:
        print(f"[Taxon Lookup] Error contacting NCBI Taxonomy API: {e}")
        return None
    

def sanitize_matrix_to_floats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Scans the entire body of the expression matrix, strips leading/trailing
    apostrophes or quotes, removes invisible whitespace/commas, forces all
    values to native float64 types, and rounds to a specified precision to
    prevent Excel string conversions.
    """
    def clean_cell(val):
        if pd.isna(val):
            return np.nan
        
        # If value is already a float or int, return directly
        if isinstance(val, (int, float)):
            return float(val)
            
        # Convert to string to clean characters
        val_str = str(val).strip()
        
        # Strip literal leading/trailing apostrophes, single quotes, or double quotes
        val_str = val_str.lstrip("'").rstrip("'").lstrip('"').rstrip('"').strip()
        
        # Strip non-breaking spaces (\xa0) and commas used in formatting (e.g. 1,234.56)
        val_str = val_str.replace('\xa0', '').replace(',', '')
        
        # Handle empty strings or string representations of missing values
        if val_str == "" or val_str.lower() in ["n/a", "na", "nan", "null", "none", "-"]:
            return np.nan
            
        try:
            return float(val_str)
        except ValueError:
            return np.nan

    # Apply element-wise cleaning to the entire numeric body of the DataFrame
    cleaned_df = df.map(clean_cell)

    # Round float values to eliminate extreme floating-point artifacts (Excel threshold)
    # 8 digits is reasonably precise (Maximum digits Excel can handle is 15)
    cleaned_df = cleaned_df.round(8)
    
    # Cast entire DataFrame explicitly to float64
    return cleaned_df.astype("float64")


def finalize_expression_matrix(combined_df: pd.DataFrame, all_genes: list, drop_unmatched_genes: bool = False) -> pd.DataFrame:
    """
    Ensures all genes (original or generated) are included in the final output dataframe,
    and sorts the gene list alphabetically.
    """
    # Ensure all input/added genes are present as rows (index)
    all_genes_sorted = sorted(list(set(all_genes)))
    
    # Re-index the DataFrame to include all genes
    # Genes missing from the expression matrix will have NaN values across samples
    final_df = combined_df.reindex(index=all_genes_sorted)
    
    # Assign clean gene index name
    final_df.index.name = "Gene_Symbol"

    # Drop rows that contain NaN across all sample columns
    # Check that user has drop_unmatched_genes set to True
    if drop_unmatched_genes:
        unmapped_count = final_df.isna().all(axis=1).sum()
        if unmapped_count > 0:
            print(f"[Matrix Post-Processing] Dropping {unmapped_count} genes with no expression data.")
            final_df = final_df.dropna(how="all", axis=0)
    
    return final_df


# --- Run ---
def run_retrieve_go_terms(gene_list, taxon_id, category_map, min_depth, batch_size):
    print("--------------------------------------------")
    print("--- Step 1: Retrieve Gene Ontology Terms ---")
    print("--------------------------------------------")

    all_records = retrieve_go_terms.fetch_all_go_terms(gene_list, taxon_id, category_map)
    
    # Filter out GO terms by depth
    godag = GODag("go-basic.obo") # Open OBO file (Download once before running: http://purl.obolibrary.org/obo/go/go-basic.obo)
    all_records = retrieve_go_terms.filter_records_by_depth(all_records, godag, min_depth)

    df = pd.DataFrame(all_records)

    # Check that genes were returned
    if df.empty:
        print()
        raise KeyError("[!] No GO terms were found.")
    
    # Check and remove obsolete genes
    print("\nChecking for obsolete GO terms...")
    go_ids = df["go_id"].tolist()
    active_terms, obsolete_terms = retrieve_go_terms.filter_obsolete_terms(go_ids, batch_size)
    df = df[~df["go_id"].isin(obsolete_terms)]
    print(f"  {len(obsolete_terms)} obsolete terms removed: {obsolete_terms}")

    # Find genes that were and were not matched
    matched_genes = set(df["gene"].unique())
    unmatched_genes = set(gene_list) - matched_genes

    print(f"\nEntries found for genes: {matched_genes}")
    print(f"No entries found for genes: {unmatched_genes}")
    print(f"Total GO annotations retrieved: {len(df)}")
    print(f"Unique GO annotations retrieved: {len(df["go_id"].unique().tolist())}")

    # Calculate how many of each type of GO term was found
    go_type_distribution = df["go_category"].value_counts()
    print("\nGO annotations by category:")
    print(go_type_distribution)

    # Print out first 20 GO terms
    print("\nFirst 20 GO terms:")
    print(df.head(20).to_string(index=False)) # <-- Print out all go terms instead? !!!

    return df, godag

# --- Run ---
def run_retrieve_genes_from_go_terms(go_df, taxon_id, godag, category_map, go_term_grouping, batch_size):
    print("--------------------------------------------")
    print("--- Step 2: Retrieve Genes from GO Terms ---")
    print("--------------------------------------------")

    # Check if the taxon is bacteria (for detecting upper case/lower case gene symbol stuff)
    domain_id = retrieve_genes_from_go_terms.get_domain(taxon_id)

    expanded_df = retrieve_genes_from_go_terms.expand_genes_from_go_terms(go_df, taxon_id, go_term_grouping, domain_id, batch_size)

    print(f"\nTotal new gene GO associations found: {len(expanded_df)}\n")

    # Count number of genes from each GO term
    term_counts = Counter(expanded_df["source_go_id"])

    if "godag" not in globals() or not godag:
        godag = GODag("go-basic.obo")
    print("Top 10 GO terms by new genes added:")
    for go_id, count in term_counts.most_common(10):
        term = godag.get(go_id)
        depth = term.depth if term else "?"
        name = term.name if term else "?"
        print(f"  {go_id} (depth {depth}) ({name}): {count} genes")

    # Calculate how many of each type of GO term was found
    if go_term_grouping in ["none", "categories"]: # <-- "all" will group across GO term categories
        expanded_df["go_category"] = expanded_df["go_category"].replace(category_map)
        type_gene_distribution = expanded_df["go_category"].value_counts()
        print("\nAdded associations by GO category:")
        print(type_gene_distribution)

    print("\nFirst 20 expanded gene rows:")
    print(expanded_df.head(20).to_string(index=False))

    # Unique new genes discovered
    new_genes = expanded_df["new_gene"].unique()
    print(f"\nUnique new genes discovered: {len(new_genes)}")

    return expanded_df

# --- Run ---
def run_retrieve_experiment_keywords(condition, model_name, api_key):
    print("--------------------------------------------")
    print("--- Step 3: Retrieve Experiment Keywords ---")
    print("--------------------------------------------")

    condition_list = retrieve_experiment_keywords.retrieve_exp_keys(condition, model_name, api_key)

    return condition_list

# --- Run ---
def run_retrieve_experiments(gene_list, expanded_df, species, condition_list, ncbi_email, max_series_return):
    print("------------------------------------------------------------------------")
    print("--- Step 4: Retrieve Experiments and Gene Expression Levels from GEO ---")
    print("------------------------------------------------------------------------")

    original_genes = list(set(gene_list))
    new_genes = expanded_df["new_gene"].unique().tolist()
    all_genes = list(set(original_genes + new_genes))

    print(f"\nTotal genes to query: {len(all_genes)}")
    print(f"Searching GEO for {len(condition_list)} conditions...")
    
    experiments = retrieve_experiments_and_expression_levels.run_geo_expression_pipeline(species, condition_list, ncbi_email, max_series_return)
    print(f"\nRetrieved {len(experiments)} experiments")

    # Remove duplicate experiments
    print("Removing duplicate experiments...")
    seen_accessions = []
    filtered_experiments = []
    for exp in experiments:
        acc = exp["accession"]

        if acc in seen_accessions:
            continue
        seen_accessions.append(acc)
        filtered_experiments.append(exp)

    print(f"Unique experiments to process: {len(filtered_experiments)}\n")

    return filtered_experiments, all_genes


def run_retrieve_expression_levels(experiments, all_genes, species, warn_file_size_mb, drop_unmatched_genes, gff_cache_dir, ncbi_email,
                                   approved_files=None, rejected_files=None):
    #series_metadata = retrieve_experiments_and_expression_levels.get_geo_series_metadata(seen_accessions, ncbi_email)

    # Flags for large files
    if approved_files is None:
        approved_files = set()
    if rejected_files is None:
        rejected_files = set()
    
    all_results = {}

    for exp in experiments:
        acc = exp["accession"]
        title = exp.get("title", "")
        print(f"Processing {acc}: {title}")

        try:
            # Resolve the strain-specific GFF for this experiment so that locus
            # tags in the expression data are mapped using the correct assembly,
            # not a cross-strain approximation.
            gse = GEOparse.get_GEO(geo=acc, destdir="./geo_cache/", silent=True)
            #gse = series_metadata[acc]

            # Resolve GFF symbol map
            id_to_symbol = retrieve_experiments_and_expression_levels.get_id_to_symbol_map_for_gse(gse, species, gff_cache_dir, ncbi_email)
            if not id_to_symbol:
                print(f"[{acc}] Warning: no strain-specific GFF resolved — "
                    "locus tag normalisation will be skipped.")
    
            # Extract expression levels (includes file size logic)
            # Pass user approvals/rejections into get_expression
            df, pending_prompt = retrieve_experiments_and_expression_levels.get_expression(
                geo_accession=acc,
                genes_of_interest=all_genes,
                id_to_symbol=id_to_symbol,
                warn_file_size_mb=warn_file_size_mb,
                approved_files=approved_files,
                rejected_files=rejected_files,
                series_title=title
            )

            # If a file was hit that needs approval, stop immediately and return information back to Streamlit
            if pending_prompt:
                print(f"[{acc}] Execution paused: Large file pending approval.")
                return pd.DataFrame(), pending_prompt

            # Format matrix if data was successfully fetched
            if df is not None and not df.empty:
                # Force entire DataFrame to float64 to purge any lingering string types
                df = df.apply(pd.to_numeric, errors='coerce')
                df = retrieve_experiments_and_expression_levels.generate_formatted_matrix(df, gse, acc)
                if not df.empty:
                    all_results[acc] = df
                    print(f"[{acc}] Success: {df.shape[0]} genes x {df.shape[1]} samples\n")
            else:
                print(f"[{acc}] Skipped or empty DataFrame.\n")

        except Exception as e:
            print(f"[{acc}] Failed: {e}\n")

    # Combine all collected matrices
    if all_results:
        # Concatenate on axis=1 using the gene names as the joining index anchor
        combined_df = pd.concat(all_results.values(), axis=1)

        # Clean and force all matrix cells to native numeric float64
        combined_df = sanitize_matrix_to_floats(combined_df)

        # Make sure all genes (original or generated) are included and sorted alphabetically
        combined_df = finalize_expression_matrix(combined_df, all_genes, drop_unmatched_genes)

        print(f"Returning combined matrix dataframe: {combined_df.shape}")
        return combined_df, None
    
    return pd.DataFrame(), None