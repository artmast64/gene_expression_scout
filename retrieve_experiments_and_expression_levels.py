# Step 4 Functions

import os
import io
import time
from datetime import datetime
import json
import re
import gzip
import zipfile
import tarfile
import math

import requests
import pandas as pd
import GEOparse

# Runtime global cache to prevent repeating the eutils/Datasets API pipeline for the same RefSeq/assembly ID
# processed_refseq_acc_set = set()
global_gff_mapping_cache: dict[str, dict[str, str]] = {}


def best_assembly_report(reports: list[dict]) -> dict | None:
    """
    From a list of NCBI Datasets assembly reports, returns the one with the
    highest priority: reference genome > representative genome > other.
    """
    if not reports:
        return None
 
    def priority(r):
        level = r.get("assembly_info", {}).get("refseq_category", "")
        if level == "reference genome":
            return 0
        if level == "representative genome":
            return 1
        return 2
 
    return sorted(reports, key=priority)[0]


def nc_accession_to_assembly_accession(nc_accession: str, ncbi_email: str) -> str | None:
    """
    Converts a RefSeq sequence accession (NC_ / NZ_) to a GCF_ assembly
    accession using the eutils esummary endpoint.
 
    The nuccore esummary record contains an 'AssemblyAcc' field (e.g.
    'GCF_000017765.1') which can then be passed directly to the NCBI
    Datasets genome/accession endpoint.
    """
    # Convert Nuccore string ID to an internal Assembly Database UID
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    params = {
        "dbfrom":  "nuccore",
        "db":      "assembly",
        "id":      nc_accession,
        "retmode": "json",
        "email":   ncbi_email,
    }
    print(f"[GFF] Looking up assembly accession for sequence '{nc_accession}' via eutils...")
    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        
        # Sanitize dangerous control characters from raw text before parsing JSON
        clean_text = re.sub(r'[\x00-\x1F\x7F]', '', resp.text)
        data = json.loads(clean_text)
    except Exception as e:
        print(f"[GFF] eutils nuccore-assembly link failed: {e}")
        return None

    # Extract internal NCBI Assembly ID
    try:
        linksets = data["linksets"][0]
        linksetdbs = linksets["linksetdbs"][0]
        assembly_uid = linksetdbs["links"][0]
        print(f"[GFF] Found internal NCBI Assembly ID: {assembly_uid}")
    except (KeyError, IndexError):
        print(f"Could not link nucleotide {nc_accession} to an Assembly record.")
        return None
    
    # Query Assembly database using the unique numeric UID
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "assembly",
        "id": assembly_uid,
        "retmode": "json",
        "email": ncbi_email,
    }

    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        print(f"[GFF] eutils assembly failed: {e}")
        return None

    # Extract the literal 'assemblyaccession' field containing the target GCF string
    try:
        # Entrez organizes summaries under the string representation of the target UID
        uid_string_key = str(assembly_uid)
        assembly_info = data["result"][uid_string_key]
        
        assembly_acc = assembly_info["assemblyaccession"]
    except KeyError:
        raise KeyError("The esummary JSON payload did not match the expected Entrez layout structure.")
 
    print(f"[GFF] Sequence '{nc_accession}' → assembly '{assembly_acc}'")
    return assembly_acc


def resolve_gff_acc_from_refseq_accession(accession: str, ncbi_email: str) -> str | None:
    """
    Given a RefSeq sequence accession (e.g. NC_012967.1) or assembly
    accession (e.g. GCF_000005845.2), returns the GFF assembly accession
    for the corresponding assembly.
 
    GCF_ / GCA_ assembly accessions are passed directly to the NCBI
    Datasets API. NC_ / NZ_ sequence accessions are first resolved to
    a GCF_ assembly accession via eutils esummary, then looked up the
    same way. The Datasets API has no sequence_accession endpoint.
    """
    print(f"[GFF] Resolving assembly for accession '{accession}'...")
 
    # Resolve NC_ / NZ_ to a GCF_ assembly accession first
    if not accession.startswith(("GCF_", "GCA_")):
        accession = nc_accession_to_assembly_accession(accession, ncbi_email)
        if not accession:
            return None
    
    return accession


def resolve_gff_acc_from_strain_name(strain: str, species: str = "Escherichia coli") -> str | None:
    """
    Fallback: queries NCBI Datasets by species + strain name and returns the
    GFF accession for the best matching assembly.
    Used only when no sequence/assembly accession could be extracted from GEO.
    """
    search_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{requests.utils.quote(species)}/dataset_report"
    params = {
        "filters.assembly_source": "refseq",
        "filters.assembly_level":  "complete_genome",
        "filters.search_text":     strain,
        "page_size":               10,
    }
 
    print(f"[GFF] Querying NCBI Datasets for '{species}' strain '{strain}'...")
    try:
        resp = requests.get(search_url, params=params, timeout=30)
        resp.raise_for_status()
    except requests.RequestException as e:
        print(f"[GFF] NCBI Datasets query failed: {e}")
        return None
 
    report = best_assembly_report(resp.json().get("reports", []))
    if not report:
        print(f"[GFF] No assemblies found for strain '{strain}'.")
        return None
 
    accession = report.get('accession')
    print(f"[GFF] Resolved assembly accession: {report.get('accession')}")
    return accession


def extract_gff_hints_from_gsm(gsm) -> tuple[str | None, str | None]:
    """
    Inspects a single GSM's metadata and returns (refseq_accession, strain_name),
    either of which may be None if not found.
 
    Checks in priority order:
      1. data_processing — often contains the reference accession used for
         alignment (e.g. "Assembly: NC_012967.1"), which is the most precise
         identifier and avoids any strain-name ambiguity.
      2. source_name — frequently contains the strain name (e.g. "REL606",
         "MG1655") when no explicit 'strain:' characteristic is present.
      3. characteristics_ch1 — structured key:value pairs; checked for an
         explicit 'strain:' entry.
    """
    refseq_accession = None
    strain_name      = None

    # Regex for RefSeq sequence accessions (NC_, NZ_) and assembly accessions (GCF_, GCA_)
    refseq_accession_re = re.compile(
        r'\b((?:GCF|GCA)_\d+\.\d+|N[CZ]_\d+\.\d+)\b'
    )
    # Regex for common E. coli strain names used in free text
    strain_name_re = re.compile(
        r'\b(K-?12|MG1655|W3110|BL21(?:\s*\(DE3\))?|DH5\w*|DH10B|'
        r'JM109|XL1-?Blue|Top10|HB101|C600|MC4100|AB1157|REL606)\b',
        re.IGNORECASE,
    )
 
    # 1. data_processing: look for a RefSeq / assembly accession
    for entry in gsm.metadata.get("data_processing", []):
        match = refseq_accession_re.search(entry)
        if match:
            refseq_accession = match.group(1)
            break   # accession found — no need to check further entries
 
    # 2. source_name: treat the whole value as a potential strain name
    if not strain_name:
        for entry in gsm.metadata.get("source_name", []):
            m = strain_name_re.search(entry)
            if m:
                strain_name = m.group(0)
                break
 
    # 3. characteristics_ch1: look for an explicit 'strain:' key
    if not strain_name:
        for char in gsm.metadata.get("characteristics_ch1", []):
            if char.lower().startswith("strain:"):
                strain_name = char.split(":", 1)[1].strip()
                break
 
    return refseq_accession, strain_name


def download_gff(assembly_accession: str, gff_cache_dir: str) -> str:
    """
    Downloads a GFF file from NCBI (if not already cached) and returns
    the local file path. Checks gff_cache_dir first to prevent redownloads.
    """
    os.makedirs(gff_cache_dir, exist_ok=True)
    cache_path = os.path.join(gff_cache_dir, f"{assembly_accession}_gff.zip")

    if os.path.exists(cache_path):
        print(f"[GFF CACHE HIT] Using locally cached archive file: {cache_path}")
        return cache_path
    
    gff_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{assembly_accession}/download"
    params = {
        "include_annotation_type": [
            "GENOME_FASTA", "GENOME_GFF", "RNA_FASTA", "CDS_FASTA", "PROT_FASTA", "SEQUENCE_REPORT"
        ],
        "hydrated": "FULLY_HYDRATED",
    }

    print(f"[GFF DOWNLOAD] Cached file not found in '{gff_cache_dir}'. Fetching from NCBI: {gff_url}")
    response = requests.get(gff_url, params=params, timeout=120)
    response.raise_for_status()

    with open(cache_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=1024 * 64):
            f.write(chunk)

    print(f"[GFF] Saved new file to cache: {cache_path}")
    return cache_path


def get_id_to_symbol_map_for_gse(gse, species: str, gff_cache_dir: str, ncbi_email: str) -> dict[str, str]:
    """
    Scans ALL GSMs in the series to identify every unique taxon/strain involved, 
    resolving and merging their GFF mappings to support multi-taxon experiments.
    Uses global in-memory cache if an assembly/RefSeq ID was already parsed.
    """
    merged_id_to_symbol = {}
    resolved_assemblies = set()
    reported_current_gse = set()

    print(f"[{gse.name}] Scanning samples for taxons/strains...")
    for gsm_name, gsm in gse.gsms.items():
        # Fallback to sample-specific organism metadata if available
        gsm_species = gsm.metadata.get("organism_ch1", [species])[0]

        # Safely pull hints first
        refseq_accession, strain_name = extract_gff_hints_from_gsm(gsm)

        # Determine a signature identifier to check for spam
        # Priority on RefSeq identifier, fallback to specific strain text
        identifier_signature = refseq_accession if refseq_accession else strain_name

        if not identifier_signature:
            continue

        # Check local experiment cache for repetition within the exact same GSE series
        if identifier_signature in reported_current_gse:
            continue
        reported_current_gse.add(identifier_signature)

        # Check global cache across separate experiments
        if identifier_signature in global_gff_mapping_cache:
            print(f"[{gse.name}] Reusing cached GFF mapping for identifier '{identifier_signature}'.")
            merged_id_to_symbol.update(global_gff_mapping_cache[identifier_signature])
            continue

        # Resolve assembly if brand new
        assembly_accession = None
        if refseq_accession:
            assembly_accession = resolve_gff_acc_from_refseq_accession(refseq_accession, ncbi_email)
        if not assembly_accession and strain_name:
            assembly_accession = resolve_gff_acc_from_strain_name(strain_name, gsm_species)
            
        if assembly_accession:
            # Check if assembly accession itself is in the global cache
            if assembly_accession in global_gff_mapping_cache:
                print(f"[{gse.name}] Reusing cached GFF mapping for assembly '{assembly_accession}'.")
                taxon_map = global_gff_mapping_cache[assembly_accession]
            else:
                try:
                    gff_path = download_gff(assembly_accession, gff_cache_dir)
                    taxon_map = build_id_to_symbol_map(gff_path)
                    # Store in global cache
                    global_gff_mapping_cache[assembly_accession] = taxon_map
                    print(f"[{gse.name}] Successfully merged GFF mappings for assembly: {assembly_accession} ({gsm_species})")
                except Exception as e:
                    print(f"[{gse.name}] Failed parsing GFF for assembly {assembly_accession}: {e}")
                    continue

            # Map signature to cached dictionary
            global_gff_mapping_cache[identifier_signature] = taxon_map
            merged_id_to_symbol.update(taxon_map)
            resolved_assemblies.add(assembly_accession)

    if not merged_id_to_symbol and not resolved_assemblies:
        print(f"[{gse.name}] Note: No new unique strain GFF mappings were needed or resolved.")
    return merged_id_to_symbol


def build_id_to_symbol_map(gff_path: str) -> dict[str, str]:
    """
    Parses an NCBI GFF file and returns a dictionary that maps every
    alternative identifier to the official gene symbol, including:
      - locus tags     (e.g. b0001, JW0001)
      - old locus tags (legacy aliases)
      - gene IDs       (e.g. 945571)
      - Dbxref entries (e.g. EcoGene:EG10998)

    Only 'gene' feature rows are processed; CDS / rRNA / etc. are skipped.
    """
    id_to_symbol: dict[str, str] = {}
    rows_parsed = 0

    with zipfile.ZipFile(gff_path, "r") as z:
        # Dynamically find the path to the internal .gff file (e.g., ncbi_dataset/data/.../genomic.gff)
        gff_internal_path = None
        for name in z.namelist():
            if name.endswith('.gff') or '/genomic.gff' in name:
                gff_internal_path = name
                break

        if gff_internal_path is None:
            raise FileNotFoundError(f"Could not find a .gff file inside the cached archive. Contents: {z.namelist()}")

        # Open that specific internal file in text mode ('r')
        with z.open(gff_internal_path, "r") as fh:
            # Zipfile reads bytes, so wrap it in TextIOWrapper to read lines as text
            text_file = io.TextIOWrapper(fh, encoding='utf-8')
            for line in text_file:
                #print(line.strip())

                if line.startswith("#"):
                    continue                         # skip comment / header lines

                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                feature_type = parts[2]
                if feature_type != "gene":           # only process gene rows
                    continue

                attributes = parts[8]
                attr_dict: dict[str, str] = {}
                for item in attributes.split(";"):
                    item = item.strip()
                    if "=" in item:
                        key, _, value = item.partition("=")
                        attr_dict[key.strip()] = value.strip()

                symbol = attr_dict.get("Name", "").strip()
                if not symbol:
                    continue                         # no gene symbol — skip

                rows_parsed += 1

                # Map the symbol to itself (handles exact-match lookups)
                id_to_symbol[symbol] = symbol
                id_to_symbol[symbol.lower()] = symbol  # case-insensitive fallback

                # Map locus tag (e.g. b0001, ECK12_RS00005) → symbol
                locus_tag = attr_dict.get("locus_tag", "").strip()
                if locus_tag:
                    id_to_symbol[locus_tag] = symbol
                    id_to_symbol[locus_tag.lower()] = symbol

                # Map old locus tags (comma-separated) → symbol
                old_locus = attr_dict.get("old_locus_tag", "").strip()
                for old_tag in old_locus.split(","):
                    old_tag = old_tag.strip()
                    if old_tag:
                        id_to_symbol[old_tag] = symbol
                        id_to_symbol[old_tag.lower()] = symbol

                # Map each Dbxref entry (e.g. "GeneID:945571", "EcoGene:EG10998") → symbol
                dbxrefs = attr_dict.get("Dbxref", "")
                for xref in dbxrefs.split(","):
                    xref = xref.strip()
                    if xref:
                        id_to_symbol[xref] = symbol
                        # Also map just the numeric/string part after the colon
                        if ":" in xref:
                            id_to_symbol[xref.split(":", 1)[1]] = symbol

    print(f"[GFF] Parsed {rows_parsed} gene rows → {len(id_to_symbol)} identifier mappings")
    return id_to_symbol


def normalize_index_to_symbols(
    df: pd.DataFrame,
    id_to_symbol: dict[str, str],
    geo_accession: str = "",
) -> pd.DataFrame:
    """
    Attempts to convert the DataFrame's index from locus tags / other IDs
    to gene symbols using the GFF-derived mapping.

    If fewer than 5 % of rows map successfully the index is likely already
    gene symbols (or an incompatible ID format), so the original index is
    kept and a warning is printed.

    Rows that cannot be mapped are dropped.
    Returns a new DataFrame with a gene-symbol index.
    """    
    original_index = df.index.tolist()
    mapped = []
    for idx in original_index:
        idx_str = str(idx).strip()
        
        # Direct match
        symbol = id_to_symbol.get(idx_str) or id_to_symbol.get(idx_str.lower())
        
        # Match without version dot suffix (e.g. "ECB_RS00005.1" -> "ECB_RS00005")
        if not symbol and "." in idx_str:
            no_dot = idx_str.split(".")[0]
            symbol = id_to_symbol.get(no_dot) or id_to_symbol.get(no_dot.lower())
            
        mapped.append(symbol)

    mapped_count = sum(1 for m in mapped if m is not None)
    tag = f"[{geo_accession}] " if geo_accession else ""

    if mapped_count == 0:
        print(f"{tag}[GFF] No index values mapped — index may already be gene symbols. Keeping original.")
        return df

    frac = mapped_count / len(original_index)
    if frac < 0.05:
        print(f"{tag}[GFF] Only {mapped_count}/{len(original_index)} rows mapped "
              f"({frac:.1%}) — keeping original index.")
        return df

    print(f"{tag}[GFF] Mapped {mapped_count}/{len(original_index)} rows "
          f"({frac:.1%}) to gene symbols.")

    df = df.copy()
    df.index = pd.Index(mapped, name="gene")
    df = df[df.index.notna()]          # drop unmapped rows (None index)
    df = df[~df.index.duplicated()]    # keep first occurrence of each symbol
    return df


# Find relevant experiments
def find_geo_experiments(species: str, condition: str, ncbi_email: str, retmax: int = 50) -> list[dict]:
    """
    Uses NCBI ESearch to find GEO Series (GSE) matching species + condition,
    then ESummary to fetch their metadata.

    ESearch query syntax:
      <condition>[ALL] AND <species>[ORGN] AND gse[ETYP]
    [ORGN] is the organism field — server-side filtering, not client-side.
    [ETYP] restricts to Series records only (not platforms or samples).
    """
    # ESearch: get matching GEO UIDs
    eutils_url  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    search_params = {
        "db":      "gds",
        "term":    f"{condition}[ALL] AND {species}[ORGN] AND gse[ETYP]",
        "sort":    "relevance", # default
        "retmax":  retmax, # was 200
        "retmode": "json",
        "email":   ncbi_email,
    }

    # Retry logic for HTTP 429 (too many requests)
    print(f"\nSearching GEO for condition '{condition}' and species '{species}'...")
    max_retries = 5
    for attempt in range(max_retries):
        search_resp = requests.get(f"{eutils_url}/esearch.fcgi", params=search_params, timeout=30)
        
        if search_resp.status_code == 429:
            retry_after = search_resp.headers.get("Retry-After")
            wait_time = 5  # default fallback wait time
            
            if retry_after:
                try:
                    wait_time = int(retry_after)
                except ValueError:
                    try:
                        target_time = datetime.strptime(retry_after, "%a, %d %b %Y %H:%M:%S GMT")
                        wait_time = max(1, int((target_time - datetime.utcnow()).total_seconds()))
                    except ValueError:
                        wait_time = 5

            print(f"Rate limited (429). Waiting {wait_time}s before retry (Attempt {attempt + 1}/{max_retries})...")
            time.sleep(wait_time)
        else:
            # Not a 429, break out of retry loop to handle response below
            break

    # Raise exception if request failed (e.g. 404, 500, or persistent 429 after all retries)
    search_resp.raise_for_status()
    search_data = search_resp.json()

    uid_list = search_data.get("esearchresult", {}).get("idlist", [])
    print(f"ESearch returned {len(uid_list)} GEO UIDs for "
          f"'{condition}' in {species} (maximum {retmax} scaled by keyword importance)")

    if not uid_list:
        return []

    time.sleep(0.4)

    # ESummary: fetch metadata for each UID
    summary_params = {
        "db":      "gds",
        "id":      ",".join(uid_list),
        "retmode": "json",
        "version": "2.0",
        "email":   ncbi_email,
    }

    print(f"Retrieving metadata from GEO results...")
    for attempt in range(max_retries):
        summary_resp = requests.get(f"{eutils_url}/esummary.fcgi", params=summary_params)
        if summary_resp.status_code == 429:
            time.sleep(3)
        else:
            break
    
    summary_resp.raise_for_status()
    summary_data = summary_resp.json()
    print("Metadata retrieved")

    experiments = []
    for uid, doc in summary_data.get("result", {}).items():
        if uid == "uids":
            continue  # skip the index key

        accession = doc.get("accession", "")
        if not accession.startswith("GSE"):
            continue  # skip non-Series records (e.g. platforms)

        # Extract sample condition labels from the 'variables' field
        variables = [v.get("description", "") for v in doc.get("variables", [])]

        experiments.append({
            "accession":   accession,
            "title":       doc.get("title", ""),
            "organism":    doc.get("taxon", ""),
            "exp_type":    doc.get("gdstype", ""),
            "n_samples":   doc.get("n_samples", 0),
            "summary":     doc.get("summary", ""),
            "variables":   variables,
            "geo_uid":     uid,
        })

    print(f"\nFound {len(experiments)} GSE Series:")
    for exp in experiments:
        print(f"  {exp['accession']}  |  {exp['title']}")
        print(f"    Type: {exp['exp_type']}  |  Samples: {exp['n_samples']}")
        print(f"    Variables: {exp['variables']}")

    return experiments


# Main pipeline
def run_geo_expression_pipeline(
    species:   str,
    condition_list: list[str],
    ncbi_email: str,
    max_series_return: int = 50,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    
    # Count conditions and create max return list
    max_return_list = []

    condition_num = 1
    for condition in condition_list:
        condition_max_return = math.floor(max_series_return - (0.05 * max_series_return * (condition_num - 1)))
        max_return_list.append(condition_max_return)
        condition_num += 1

    # Get experiments for all conditions in the expanded condition list
    experiments = []
    for condition, condition_max_return in zip(condition_list, max_return_list):
        condition_experiments = find_geo_experiments(species, condition, ncbi_email, condition_max_return)

        if not condition_experiments:
            print("No matching GEO experiments found.")
            # return pd.DataFrame(), pd.DataFrame()

        experiments.extend(condition_experiments)

    if not experiments:
        raise KeyError("No matching GEO experiments found.")
        # return pd.DataFrame(), pd.DataFrame()
    
    return experiments


def extract_nested_files(content_bytes: bytes, filename: str) -> list[tuple[str, bytes]]:
    """
    Recursively extracts files from compressed archives (.gz, .tar, .tar.gz, .zip).
    Returns a list of tuples containing (individual_filename, uncompressed_bytes).
    """
    extracted_files = []
    lower_name = filename.lower()

    try:
        if lower_name.endswith(".tar.gz") or lower_name.endswith(".tgz"):
            with tarfile.open(fileobj=io.BytesIO(content_bytes), mode="r:gz") as tar:
                for member in tar.getmembers():
                    if member.isfile() and not os.path.basename(member.name).startswith('.'):
                        print(f"  [Archive Contents] Found file in TAR.GZ: {member.name}")
                        f_bytes = tar.extractfile(member).read()
                        extracted_files.extend(extract_nested_files(f_bytes, member.name))
                        
        elif lower_name.endswith(".tar"):
            with tarfile.open(fileobj=io.BytesIO(content_bytes), mode="r:") as tar:
                for member in tar.getmembers():
                    if member.isfile() and not os.path.basename(member.name).startswith('.'):
                        print(f"  [Archive Contents] Found file in TAR: {member.name}")
                        f_bytes = tar.extractfile(member).read()
                        extracted_files.extend(extract_nested_files(f_bytes, member.name))
                        
        elif lower_name.endswith(".zip"):
            with zipfile.ZipFile(io.BytesIO(content_bytes)) as z:
                for name in z.namelist():
                    if not os.path.basename(name).startswith('.'):
                        print(f"  [Archive Contents] Found file in ZIP: {name}")
                        f_bytes = z.read(name)
                        extracted_files.extend(extract_nested_files(f_bytes, name))
                        
        elif lower_name.endswith(".gz") and not lower_name.endswith(".tar.gz"):
            # Decompress single file
            decompressed = gzip.decompress(content_bytes)
            # Strip the .gz extension to check the inner filename format
            inner_name = filename[:-3]
            print(f"  [Archive Contents] Decompressed single GZ wrapper for: {inner_name}")
            extracted_files.extend(extract_nested_files(decompressed, inner_name))
            
        else:
            # Base case: uncompressed or final file level
            extracted_files.append((filename, content_bytes))
            
    except Exception as e:
        print(f"  [Archive Extraction Warning] Failed to unpack archive tier for {filename}: {e}")
        
    return extracted_files


def get_rnaseq_expression_for_genes(
    geo_accession: str,
    genes_of_interest: list[str],
    id_to_symbol: dict[str, str],
    warn_file_size_gb: int = 1,
    approved_files=None,
    rejected_files=None,
    series_title="",
    supp_cache_dir: str = "supp_cache",
) -> tuple[pd.DataFrame, dict | None]:
    """
    For RNA-seq GSE accessions where GEOparse returns empty tables,
    fetches and iterates through all valid supplemental files, unpacking nested 
    compressions to aggregate gene expression metrics across all formats.
    """
    if approved_files is None:
        approved_files = set()
    if rejected_files is None:
        rejected_files = set()

    # Find supplemental file URLs from the GSE metadata
    gse = GEOparse.get_GEO(geo=geo_accession, destdir="./geo_cache/", silent=True)
    suppl_files = gse.metadata.get("supplementary_file", [])
    print(f"[{geo_accession}] Supplemental files:")
    for f in suppl_files:
        print(f"  {f}")

    # Identify non-raw valid data files
    # RNA-seq studies typically provide a single counts or TPM matrix
    skip_keywords = [".cel", ".idat", ".bam", ".fastq", ".fasta"]
    valid_extensions = (".gz", ".tar.gz", ".tar", ".tgz", ".zip", ".txt", ".csv", ".tsv", ".xls", ".xlsx")
    
    candidates = [
        f for f in suppl_files
        if not any(k in f.lower() for k in skip_keywords) and f.lower().endswith(valid_extensions)
    ]

    if not candidates:
        print(f"[{geo_accession}] No suitable expression supplemental files found.")
        return pd.DataFrame(), None
    
    combined_suppl_df = pd.DataFrame()

    # Create GSE-specific subdirectory inside the supplementary cache folder
    gse_supp_dir = os.path.join(supp_cache_dir, geo_accession)
    os.makedirs(gse_supp_dir, exist_ok=True)

    for file_url in candidates:
        # Check if user rejected this file earlier
        if file_url in rejected_files:
            print(f"[{geo_accession}] Skipping rejected file: {file_url}")
            continue

        filename = os.path.basename(file_url)
        download_url = file_url.replace("ftp://", "https://")
        local_filepath = os.path.join(gse_supp_dir, filename)

        content_bytes = None

        # 1. Check local cache first
        if os.path.exists(local_filepath):
            print(f"[{geo_accession}] [SUPP CACHE HIT] Loading locally cached supplementary file: {local_filepath}")
            try:
                with open(local_filepath, "rb") as f:
                    content_bytes = f.read()
            except Exception as e:
                print(f"[{geo_accession}] Failed reading cached file {local_filepath}: {e}")

        # 2. Download if not cached
        if content_bytes is None:
            # Check remote file size before downloading
            size_gb = get_remote_file_size_gb(file_url)

            if size_gb and size_gb > warn_file_size_gb:
                if file_url not in approved_files:
                    print(f"[{geo_accession}] Large file detected ({size_gb} GB). Requesting user approval...")
                    prompt_info = {
                        "series": f"{geo_accession}: {series_title}",
                        "filename": filename,
                        "url": file_url,
                        "size_gb": size_gb
                    }
                    # Pause and return prompt to caller
                    return pd.DataFrame(), prompt_info               
                
            print("  ↳")
            print(f"[{geo_accession}] Processing supplemental root file entry: {file_url}")
            print(f"[{geo_accession}] File size {file_url} ({size_gb} GB)")

            headers = {
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
            }

            try:
                response = requests.get(download_url, headers=headers, timeout=90)
                response.raise_for_status()
                content_bytes = response.content

                # Save file to local cache subdirectory
                with open(local_filepath, "wb") as f:
                    f.write(content_bytes)
                print(f"[{geo_accession}] Saved supplementary file to cache: {local_filepath}")

            except Exception as e:
                print(f"[{geo_accession}] Error downloading/processing {file_url}: {e}")
                continue

        # 3. Unpack and parse metrics
        try:
            # Unpack all potential nested data tables
            all_files = extract_nested_files(content_bytes, filename)
            
            for fname, fbytes in all_files:
                is_valid_format = fname.lower().endswith((".csv", ".tsv", ".txt", ".xls", ".xlsx"))

                if not is_valid_format:
                    print(f"  [-] IGNORED file (unsupported tabular extension): {fname}")
                    continue

                print(f"  [+] KEEPING & PARSING data table: {fname}")
                df = detect_header_rows_and_parse_tabular(fbytes, fname)

                if df.empty:
                    print(f"  [-] SKIPPED data table (empty dataframe parsed): {fname}")
                    continue

                # Disambiguate generic column names (e.g. "Count") using file names
                sample_label = os.path.splitext(fname)[0] # Extract filename without extension
                # If there's only 1 column and it has a non-unique generic name, rename it to sample_label
                generic_names = {"count", "counts", "readcount", "read_count", "val", "value", "expression", "fpkm", "tpm"}
                renamed_cols = {col: sample_label for col in df.columns if str(col).strip().lower() in generic_names}
                if renamed_cols:
                    df = df.rename(columns=renamed_cols)
                
                # Normalize index to gene symbols
                df = normalize_index_to_symbols(df, id_to_symbol, geo_accession)
                
                # Align rows with targets
                found_genes = df.index.intersection(genes_of_interest)
                df_filtered = df.loc[found_genes]
                
                if not df_filtered.empty:
                    # Combine columns into our running collection
                    # Force all values to float to guarantee no leftover string cells from .xls files
                    df_filtered = df_filtered.apply(pd.to_numeric, errors='coerce')
                    combined_suppl_df = pd.concat([combined_suppl_df, df_filtered], axis=1)
                    
        except Exception as e:
            print(f"[{geo_accession}] Error unpacking/processing contents of {filename}: {e}")

    if combined_suppl_df.empty:
        return pd.DataFrame(), None

    # Handle duplicate column names if identical sample files were processed twice
    combined_suppl_df = combined_suppl_df.loc[:, ~combined_suppl_df.columns.duplicated()]
    combined_suppl_df = combined_suppl_df.groupby(combined_suppl_df.index).mean()
    
    print(f"[{geo_accession}] Final aggregated matrix shape: {combined_suppl_df.shape}")
    return combined_suppl_df, None


def get_expression_for_genes(
    geo_accession: str,
    genes_of_interest: list[str],
    id_to_symbol: dict[str, str],
    probe_id_col: str = None,       # e.g. "ID_REF", "ID", "SPOT_ID"
    value_col: str = None,          # e.g. "VALUE", "LOG2_RATIO", "RMA_VALUE"
    gene_symbol_col: str = None,    # e.g. "Gene Symbol", "ORF", "GENE", "gene_id"
) -> pd.DataFrame:
    """
    Retrieves, processes, and filters gene expression data from NCBI GEO.

    This function downloads a GEO Series (GSE) dataset, constructs an expression
    matrix across all available samples (GSMs), maps platform probe IDs (GPL) to 
    gene identifiers, normalizes identifiers to standard symbols, and filters 
    the final matrix to a targeted list of genes.

    If columns for probe IDs, expression values, or gene symbols are not explicitly 
    provided, the function attempts to heuristically auto-detect them from the 
    dataset metadata.
    """
    
    gse = GEOparse.get_GEO(geo=geo_accession, destdir="./geo_cache/", silent=True)

    # --- Resolve GSM columns ---
    first_gsm = list(gse.gsms.values())[0]
    gsm_cols = list(first_gsm.table.columns)

    # Auto-detect probe ID column if not specified
    if probe_id_col is None:
        candidates = ["ID_REF", "ID", "SPOT_ID", "Reporter", "PROBE_ID"]
        probe_id_col = next((c for c in candidates if c in gsm_cols), gsm_cols[0])
        print(f"[{geo_accession}] Using probe ID column: '{probe_id_col}'")

    # Auto-detect value column if not specified
    if value_col is None:
        candidates = ["VALUE", "LOG2_RATIO", "RMA_VALUE", "SIGNAL", "INTENSITY"]
        value_col = next((c for c in candidates if c in gsm_cols), gsm_cols[1])
        print(f"[{geo_accession}] Using value column: '{value_col}'")

    # --- Build expression matrix ---
    try:
        expr_matrix = pd.DataFrame({
            gsm_name: gsm.table.set_index(probe_id_col)[value_col]
            for gsm_name, gsm in gse.gsms.items()
        })
        
        # Clean invisible spaces/commas and force SOFT table values to float64
        expr_matrix = expr_matrix.astype(str).replace(r'[\xa0\r\t,]', '', regex=True)
        expr_matrix = expr_matrix.apply(pd.to_numeric, errors='coerce')
        
    except KeyError as e:
        raise KeyError(
            f"Column not found in GSM table. "
            f"Available columns: {gsm_cols}. "
            f"Pass probe_id_col= and value_col= explicitly."
        ) from e

    # --- Resolve GPL columns ---
    if not gse.gpls:
        print(f"[{geo_accession}] No GPL found; returning probe-level matrix.")
        return expr_matrix

    platform_id = list(gse.gpls.keys())[0]
    gpl = gse.gpls[platform_id]
    gpl_cols = list(gpl.table.columns)

    # Auto-detect gene symbol column if not specified
    if gene_symbol_col is None:
        candidates = [
            "Gene Symbol", "GENE_SYMBOL", "ORF", "gene_assignment",
            "GENE", "gene_id", "Gene_Symbol", "Symbol",
            # E. coli-specific
            "Gene name", "gene_name", "GeneName", "Locus"
        ]
        gene_symbol_col = next((c for c in candidates if c in gpl_cols), None)
        if gene_symbol_col:
            print(f"[{geo_accession}] Using gene symbol column: '{gene_symbol_col}'")
        else:
            print(
                f"[{geo_accession}] No gene symbol column found. "
                f"Available GPL columns: {gpl_cols}. "
                f"Pass gene_symbol_col= explicitly."
            )
            # return expr_matrix  # return probe-level so you can still inspect
            # # ^^^ This was returning non-gene matches, so it is commented out
            # Ex: GSE49914 was returning its entire matrix, even though genes were listed as sequentially increasing digits rather than symbols

    # --- Map probes to genes ---
    # GPL index column is usually "ID" — find it defensively
    gpl_id_col = "ID" if "ID" in gpl_cols else gpl_cols[0]
    probe_to_gene = gpl.table.set_index(gpl_id_col)[gene_symbol_col]

    expr_matrix["gene"] = expr_matrix.index.map(probe_to_gene)
    expr_matrix = expr_matrix.dropna(subset=["gene"])

    # Some platforms put multiple gene names in one cell (e.g. "recA /// reca2")
    # Split and explode so each gene gets its own row, then average
    expr_matrix["gene"] = expr_matrix["gene"].str.split(r"\s*///\s*")
    expr_matrix = expr_matrix.explode("gene")
    expr_matrix["gene"] = expr_matrix["gene"].str.strip()

    expr_matrix = (
        expr_matrix
        .groupby("gene", sort=False)
        .mean(numeric_only=True)
    )

    # Normalise whatever gene identifiers the platform used to standard symbols
    expr_matrix = normalize_index_to_symbols(expr_matrix, id_to_symbol, geo_accession)

    # --- Filter to genes of interest ---
    found = expr_matrix.index.intersection(genes_of_interest)
    missing = set(genes_of_interest) - set(found)
    if missing:
        print(f"[{geo_accession}] Genes not found in platform: {missing}")

    return expr_matrix.loc[found]


def get_expression(
    geo_accession: str,
    genes_of_interest: list[str],
    id_to_symbol: dict[str, str],
    warn_file_size_gb: int = 1,
    approved_files=None,
    rejected_files=None,
    series_title="",
    probe_id_col: str = None,
    value_col: str = None,
    gene_symbol_col: str = None,
    supp_cache_dir: str = "supp_cache",
) -> pd.DataFrame:
    """
    Dispatcher: automatically routes to microarray or RNA-seq path
    depending on whether GEOparse finds table data in the SOFT file.

    id_to_symbol should be the strain-specific map produced by
    get_id_to_symbol_map_for_gse() for this experiment's GSE object.
    """
    # Download metadata
    # "brief" skips downloading the actual expression table
    gse = GEOparse.get_GEO(geo=geo_accession, how="brief", destdir="./geo_cache/", silent=True)

    # Check if the first GSM has table data
    first_gsm = list(gse.gsms.values())[0]
    has_table_data = not first_gsm.table.empty

    if has_table_data:
        print(f"[{geo_accession}] Microarray data detected — using SOFT table path.")
        return get_expression_for_genes(
            geo_accession, genes_of_interest,
            probe_id_col=probe_id_col,
            value_col=value_col,
            gene_symbol_col=gene_symbol_col,
            id_to_symbol=id_to_symbol
        )
    else:
        print(f"[{geo_accession}] Empty SOFT tables — using RNA-seq supplemental file path.")
        return get_rnaseq_expression_for_genes(
            geo_accession,
            genes_of_interest,
            id_to_symbol,
            warn_file_size_gb,
            approved_files,
            rejected_files,
            series_title,
            supp_cache_dir,
        )


def get_remote_file_size_gb(file_url: str) -> float | None:
    """
    Performs an HTTP HEAD request to fetch Content-Length without downloading the body.
    Returns size in Gigabytes (GB), or None if unavailable.
    """
    url = file_url.replace("ftp://", "https://")
    # Note: look into using API keys instead of modifying the header - it allows for more frequent requests to NCBI/GEO and may remove HTTP 403 errors
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
    }
    try:
        response = requests.head(url, headers=headers, allow_redirects=True, timeout=(5,30)) # 5 seconds to connect, 30 to return headers for large files
        content_length = response.headers.get("Content-Length")
        if content_length:
            size_bytes = int(content_length)
            return round(size_bytes / (1024 * 1024 * 1024), 2)  # Convert bytes to GB
    except Exception as e:
        print(f"[Size Check Warning] Couldn't fetch size for {url}: {e}")
    return None


def detect_header_rows_and_parse_tabular(file_bytes: bytes, filename: str) -> pd.DataFrame:
    """
    Dynamically identifies and parses metadata/header rows in tabular data.
    Strips comment lines, detects headers, drops known genomic metadata columns,
    and retains only numeric sample data columns.
    """
    fn_lower = filename.lower()
    header_rows_count = 0

    # Row array extraction based on file extension
    if fn_lower.endswith((".xls", ".xlsx")):
        try:
            # Read the full raw layout as a baseline matrix with no indexing columns or headers
            raw_excel_df = pd.read_excel(io.BytesIO(file_bytes), header=None)
            # Transform rows into structural tokens to feed into the scanner loop
            lines = raw_excel_df.astype(str).values.tolist()
        except Exception as e:
            print(f"  [Parser Error] Failed to read layout for Excel file {filename}: {e}")
            return pd.DataFrame()
    else:
        # Decode text content
        text_content = file_bytes.decode("utf-8", errors="ignore")
        raw_lines = text_content.splitlines()
        
        # Strip out metadata comment lines starting with '#' or '!'
        # '#' - metadata comments for featureCounts, VCFs, GFFs
        # '!' - metadata comments for GEO SOFT files
        clean_lines = [line for line in raw_lines if line.strip() and not line.strip().startswith(("#", "!"))]
        
        if not clean_lines:
            return pd.DataFrame()
        
        # Re-encode clean text without comments
        file_bytes = "\n".join(clean_lines).encode("utf-8")
        
        sep = "\t" if fn_lower.endswith((".tsv", ".txt")) else ","
        # Keep tokens as raw string cells for checking
        lines = [line.strip().split(sep) for line in clean_lines]

    if not lines:
        return pd.DataFrame()

    # Scan rows to find where numeric data starts
    for idx, tokens in enumerate(lines):
        numeric_signals = []
        # Check all cells starting from index 1 (ignoring column 0, which holds the gene identifier strings)
        for token in tokens[1:]:
            clean_token = str(token).strip()
            # Catch standard empty cells or cosmetic N/A labels
            if not clean_token or clean_token.lower() in ["n/a", "na", "nan", "null", "none", "."]:
                continue
            try:
                float(clean_token)
                numeric_signals.append(True)
            except ValueError:
                numeric_signals.append(False)
        
        # If row contains numeric gene expression values, the header block ends here
        if numeric_signals and any(numeric_signals) and all(x for x in numeric_signals if x is True):
            header_rows_count = idx
            break
            
    print(f"  [Parser] Detected {header_rows_count} header/metadata row(s) in {filename}")

    # Parsing final matrix with automatic multi-header assignment
    try:
        if fn_lower.endswith((".xls", ".xlsx")):
            header_arg = list(range(header_rows_count)) if header_rows_count > 1 else (0 if header_rows_count == 1 else None)
            df = pd.read_excel(io.BytesIO(file_bytes), header=header_arg, index_col=0)
        else:
            sep = "\t" if fn_lower.endswith((".tsv", ".txt")) else ","
            header_arg = list(range(header_rows_count)) if header_rows_count > 1 else (0 if header_rows_count == 1 else None)
            df = pd.read_csv(io.BytesIO(file_bytes), sep=sep, header=header_arg, index_col=0)
    except Exception as e:
        print(f"  [Parser Error] Failed to construct DataFrame from {filename}: {e}")
        return pd.DataFrame()
    
    # Drop known genomic annotation columns by name (case-insensitive match)
    genomic_metadata_cols = {
        "chr", "chromosome", "seqid", "contig",
        "start", "end", "stop", "strand", "length", 
        "gene_length", "transcript_length", "biotype", "description"
    }
    
    filtered_cols = []
    for col in df.columns:
        # Check standard string representation of column names
        col_name_clean = str(col).strip().lower()
        if col_name_clean not in genomic_metadata_cols:
            filtered_cols.append(col)
        else:
            print(f"  [-] Dropping genomic metadata column: '{col}'")

    df = df[filtered_cols]

    # Keep remaining columns that contain numeric expression data
    numeric_cols = []
    for col in df.columns:
        col_series = df[col]

        # 1. Thoroughly clean string values if column is object type
        if col_series.dtype == "object":
            col_series = (
                col_series.astype(str)
                # Remove non-breaking spaces (\xa0), carriage returns, and commas in numbers (e.g. 1,234.5)
                .str.replace(r'[\xa0\r\t]', '', regex=True)
                .str.replace(',', '', regex=False)
                .str.strip()
            )
            # Replace standard text representations of missing data
            col_series = col_series.replace(["", "n/a", "N/A", "NA", "nan", "NaN", "null", "None", "-"], pd.NA)

        # 2. Force conversion to float
        converted = pd.to_numeric(col_series, errors='coerce')

        # 3. Keep column if it contains valid numeric expression numbers
        if converted.notna().sum() > 0:
            df[col] = converted
            numeric_cols.append(col)

    if not numeric_cols:
        print(f"  [Parser Warning] No numeric expression columns remaining in {filename}")
        return pd.DataFrame()

    return df[numeric_cols]


def generate_formatted_matrix(df: pd.DataFrame, gse, geo_accession: str) -> pd.DataFrame:
    """
    Cleans raw data matrices and maps metadata to a MultiIndex column structure,
    preserving the gene index for downstream concatenation.
    """
    # 1. Clean out cosmetic or statistical summary columns consistently
    unwanted_keywords = ["average", "mean", "std", "avg", "unnamed"]
    valid_cols = [col for col in df.columns if not any(kw in str(col).lower() for kw in unwanted_keywords)]
    df = df.loc[:, valid_cols]
    
    # Drop rows or columns that are entirely unpopulated
    df = df.dropna(how="all", axis=0).dropna(how="all", axis=1)

    # 2. Compile custom metadata alignment rows
    geo_series_row  = []
    gse_name_row    = []
    geo_sample_row  = []
    gsm_name_row    = []

    for col in df.columns:
        # If the column header is a tuple (due to detecting multi-row headers), extract the active text string
        col_str = str(col[-1]).strip() if isinstance(col, tuple) else str(col).strip()
        
        gsm_obj = gse.gsms.get(col_str)
        
        # Fallback check: see if the column name can be tracked anywhere inside sample titles
        if not gsm_obj:
            for gsm_id, gsm in gse.gsms.items():
                if gsm_id in col_str or col_str in gsm.metadata.get("title", [""])[0]:
                    gsm_obj = gsm
                    break

        if gsm_obj:
            geo_series_row.append(geo_accession)
            gse_name_row.append(gse.metadata.get("title", [""])[0])
            geo_sample_row.append(gsm_obj.name)
            gsm_name_row.append(gsm_obj.metadata.get("title", [col_str])[0])
        else:
            geo_series_row.append(geo_accession)
            gse_name_row.append(gse.metadata.get("title", [""])[0])
            geo_sample_row.append(col_str)
            gsm_name_row.append(col_str)

    # 3. Apply MultiIndex columns, keeping Gene Symbols as the standard dataframe Index
    export_df = df.copy()
    export_df.columns = pd.MultiIndex.from_arrays(
        [geo_series_row, gse_name_row, geo_sample_row, gsm_name_row],
        names=["GSE Accession", "GSE Name", "GSM Accession", "GSM Name"]
    )
    
    return export_df