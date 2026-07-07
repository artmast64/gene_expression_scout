# Gene Expression Scout

from collections import Counter
from pprint import pprint

import pandas as pd
from goatools.obo_parser import GODag
import GEOparse

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


def print_settings(species, condition, gene_list, go_term_grouping, batch_size, min_depth, max_series_return):
    print("\n--- Settings: ---")
    print(f"  Species name: {species}")
    print(f"  Factor: {condition}")
    print(f"  List of genes: {gene_list}")
    print(f"  Gene Retrieval Method: {go_term_grouping}")
    print(f"  Batch Size: {batch_size}")
    print(f"  Minimum GO Term Depth: {min_depth}")
    print(f"  Maximum GEO series results per keyword: {max_series_return}")
    print()


def identify_taxon(species: str):
    # Lookup table of common species
    # Replaced 562 with 83333 for E. coli
    species_df = pd.DataFrame({
        "Common Name": ["Mouse", "Human", "Rat", "Thale cress", "Zebrafish", "Fruit fly", "Baker's yeast", "Roundworm", "E. coli", "Cattle"],
        "Scientific Name": ["Mus musculus", "Homo sapiens", "Rattus norvegicus", "Arabidopsis thaliana", "Danio rerio", "Drosophilia melanogaster",
                            "Saccharomyces cerevisiae", "Caenorhabditis elegans", "Escherichia coli", "Bos taurus"],
        "Taxon ID": [10090, 9606, 10116, 3702, 7955, 7227, 4932, 6239, 83333, 9913],
    })

    # Find taxon ID from species name (try common name first, then scientific name)
    try:
        taxon_id = species_df.loc[species_df["Common Name"] == species, "Taxon ID"].values[0]
    except:
        try:
            taxon_id = species_df.loc[species_df["Scientific Name"] == species, "Taxon ID"].values[0]
        except:
            print(f"Taxon ID not found for species '{species}'")
            taxon_id = input("Type in the taxon ID for the species: ")
    print(f"Taxon ID for {species}: {taxon_id}")
    return taxon_id


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

    # Save to CSV <-- for testing !!!
    df.to_csv("output_go_terms.csv", index=False)
    print("\nSaved to output_go_terms.csv")

    return df, godag

#go_df, godag = run_retrieve_go_terms(gene_list, taxon_id, category_map, min_depth, batch_size)


#go_df = pd.read_csv("output_go_terms.csv") # <-- load saved .csv for testing !!!
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
    # print(sorted(new_genes))

    # Save to .csv <-- for testing !!!
    expanded_df.to_csv("output_expanded_genes.csv", index=False)
    print("\nSaved to output_expanded_genes.csv")

    return expanded_df

#expanded_df = run_retrieve_genes_from_go_terms(go_df, taxon_id, godag, category_map, go_term_grouping, batch_size)


# --- Run ---
def run_retrieve_experiment_keywords(condition, model_name, api_key):
    print("--------------------------------------------")
    print("--- Step 3: Retrieve Experiment Keywords ---")
    print("--------------------------------------------")

    condition_list = retrieve_experiment_keywords.retrieve_exp_keys(condition, model_name, api_key)

    return condition_list

#condition_list = run_retrieve_experiment_keywords(condition, model_name, api_key)


# df = pd.read_csv("output_go_terms.csv") # <-- load saved .csv for testing !!!
# expanded_df = pd.read_csv("output_expanded_genes.csv") # <-- load saved .csv for testing !!!
# --- Run ---
def run_retrieve_experiments(go_df, expanded_df, species, condition_list, ncbi_email, max_series_return):
    print("------------------------------------------------------------------------")
    print("--- Step 4: Retrieve Experiments and Gene Expression Levels from GEO ---")
    print("------------------------------------------------------------------------")

    original_genes = go_df["gene"].unique().tolist()
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


def run_retrieve_expression_levels(experiments, all_genes, species, gff_cache_dir, ncbi_email):
    #series_metadata = retrieve_experiments_and_expression_levels.get_geo_series_metadata(seen_accessions, ncbi_email)
    all_results = {}

    for exp in experiments:
        acc = exp["accession"]
        print(f"Processing {acc}: {exp['title']}")
        try:
            # Resolve the strain-specific GFF for this experiment so that locus
            # tags in the expression data are mapped using the correct assembly,
            # not a cross-strain approximation.
            gse = GEOparse.get_GEO(geo=acc, destdir="./geo_cache/", silent=True)
            #gse = series_metadata[acc]
            id_to_symbol = retrieve_experiments_and_expression_levels.get_id_to_symbol_map_for_gse(gse, species, gff_cache_dir, ncbi_email)
            if not id_to_symbol:
                print(f"[{acc}] Warning: no strain-specific GFF resolved — "
                    "locus tag normalisation will be skipped.")
    
            df = retrieve_experiments_and_expression_levels.get_expression(acc, all_genes, id_to_symbol)
            if not df.empty:
                all_results[acc] = df
                print(f"[{acc}] Success: {df.shape[0]} genes x {df.shape[1]} samples\n")
            else:
                print(f"[{acc}] Returned empty DataFrame — check supplemental files.\n")
        except Exception as e:
            print(f"[{acc}] Failed: {e}\n")

    # Combine across studies (genes as rows, multi-index columns by study+sample)
    if all_results:
        combined_df = pd.concat(all_results, axis=1)
        print("First 5 rows of combined dataset:")
        print(combined_df.head())

        # Export to .csv file
        combined_df.to_csv("output_expression_matrix.csv")
        print(f"\nSaved combined matrix: {combined_df.shape}")
        print("Saved to output_expression_matrix.csv")

        return combined_df

#experiments, all_genes = run_retrieve_experiments(all_genes, species, condition_list, ncbi_email, max_series_return)
#combined_df = run_retrieve_expression_levels(experiments, all_genes, species, gff_cache_dir, ncbi_email)