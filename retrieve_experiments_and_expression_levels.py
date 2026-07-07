# Step 4 Functions

import os
import io
import time
import re
import gzip
import zipfile
import tarfile
import math

import requests
import pandas as pd
import GEOparse


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
        data = resp.json()
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
            print(f"[GFF] RefSeq accession from data_processing: '{refseq_accession}'")
            break   # accession found — no need to check further entries
 
    # 2. source_name: treat the whole value as a potential strain name
    if not strain_name:
        for entry in gsm.metadata.get("source_name", []):
            m = strain_name_re.search(entry)
            if m:
                strain_name = m.group(0)
                print(f"[GFF] Strain from source_name: '{strain_name}'")
                break
 
    # 3. characteristics_ch1: look for an explicit 'strain:' key
    if not strain_name:
        for char in gsm.metadata.get("characteristics_ch1", []):
            if char.lower().startswith("strain:"):
                strain_name = char.split(":", 1)[1].strip()
                print(f"[GFF] Strain from characteristics_ch1: '{strain_name}'")
                break
 
    return refseq_accession, strain_name


def download_gff(assembly_accession: str, gff_cache_dir: str) -> str:
    """
    Downloads a GFF file from NCBI (if not already cached) and returns
    the local file path. Files are cached in GFF_CACHE_DIR by assembly
    accession so each assembly is only downloaded once per session.
    """
    os.makedirs(gff_cache_dir, exist_ok=True)
    cache_path = os.path.join(gff_cache_dir, f"{assembly_accession}_gff.zip")

    if os.path.exists(cache_path):
        print(f"[GFF] Using cached file: {cache_path}")
        return cache_path
    
    gff_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{assembly_accession}/download"
    params = {
        "include_annotation_type": ["GENOME_FASTA",
        "GENOME_GFF",
        "RNA_FASTA",
        "CDS_FASTA",
        "PROT_FASTA",
        "SEQUENCE_REPORT"
        ],
        "hydrated": "FULLY_HYDRATED",
    }

    print(f"[GFF] Downloading GFF from NCBI: {gff_url}")
    response = requests.get(gff_url, params=params, timeout=120)
    response.raise_for_status()

    with open(cache_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=1024 * 64):
            f.write(chunk)

    print(f"[GFF] Saved to: {cache_path}")
    return cache_path


def get_id_to_symbol_map_for_gse(gse, species: str, gff_cache_dir: str, ncbi_email: str) -> dict[str, str]:
    """
    High-level helper: resolves the correct strain-specific GFF for this
    experiment and returns the id→symbol mapping dict.
 
    Resolution priority (stops at first success):
      1. Sample the first GSM for a RefSeq/assembly accession in
         data_processing — the most precise identifier, directly tied to
         the genome the authors actually aligned to.
      2. Use a strain name from source_name or characteristics_ch1 to
         search NCBI Datasets — less precise but usually sufficient.
      3. Return an empty dict, causing normalize_index_to_symbols to leave
         the expression matrix index unchanged rather than mismapping genes.
    """
    # Sample the first GSM for hints
    first_gsm = next(iter(gse.gsms.values()), None)
    #first_gsm = next(iter(gse["samples"]), None)

    refseq_accession, strain_name = (None, None)
    if first_gsm:
        refseq_accession, strain_name = extract_gff_hints_from_gsm(first_gsm)
 
    # Priority 1: resolve directly from the reference accession
    assembly_accession = None
    if refseq_accession:
        assembly_accession = resolve_gff_acc_from_refseq_accession(refseq_accession, ncbi_email)
 
    # Priority 2: fall back to strain name search
    if not assembly_accession and strain_name:
        assembly_accession = resolve_gff_acc_from_strain_name(strain_name, species)
 
    if not assembly_accession:
        print("[GFF] Could not resolve a GFF for this experiment — "
              "locus tag normalisation will be skipped.")
        return {}

    gff_path = download_gff(assembly_accession, gff_cache_dir)
    return build_id_to_symbol_map(gff_path)


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
    mapped = [id_to_symbol.get(idx) or id_to_symbol.get(str(idx).lower()) for idx in original_index]

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

    # Call GEO API
    print(f"\nSearching GEO for condition '{condition}' and species '{species}'...")
    search_resp = requests.get(f"{eutils_url}/esearch.fcgi", params=search_params)
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

    # Call GEO API
    print(f"Retrieving metadata from GEO results...")
    summary_resp = requests.get(f"{eutils_url}/esummary.fcgi", params=summary_params)
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


def get_rnaseq_expression_for_genes(
    geo_accession: str,
    genes_of_interest: list[str],
    id_to_symbol: dict[str, str],
) -> pd.DataFrame:
    """
    For RNA-seq GSE accessions where GEOparse returns empty tables,
    fetch the series matrix or supplemental count file from GEO FTP.
    """
    gse = GEOparse.get_GEO(geo=geo_accession, destdir="./geo_cache/", silent=True)

    # Step 1: find supplemental file URLs from the GSE metadata
    suppl_files = gse.metadata.get("supplementary_file", [])
    print(f"[{geo_accession}] Supplemental files:")
    for f in suppl_files:
        print(f"  {f}")

    # Step 2: filter to likely count/expression matrix files
    # RNA-seq studies typically provide a single counts or TPM matrix
    #preferred_keywords = ["count", "tpm", "fpkm", "rpkm", "expression", "matrix"]
    skip_keywords = [".cel", ".idat", ".bam", ".fastq"] # ["raw", ".cel", ".idat", ".bam", ".fastq"]

    candidates = [
        f for f in suppl_files
        #if any(k in f.lower() for k in preferred_keywords)
        #and not any(k in f.lower() for k in skip_keywords)
        if not any(k in f.lower() for k in skip_keywords)
        and f.endswith((".gz", ".tar.gz", ".tar", ".txt", ".csv", ".tsv", ".xls", ".xlsx"))
    ]

    if not candidates:
        print(f"[{geo_accession}] No suitable supplemental file found automatically.")
        # print("Try passing a file URL directly via the suppl_url= parameter.")
        return pd.DataFrame()

    # Use the first candidate (inspect the printed list to pick a better one if needed)
    file_url = candidates[0]
    download_url = file_url.replace("ftp://", "https://")

    # Get file size using a HEAD request
    try:
        # headers=None or standard user-agent to ensure NCBI accepts the request
        head_response = requests.head(download_url, timeout=30, allow_redirects=True)
        head_response.raise_for_status()
        
        # Get content length in bytes (default to 0 if header is missing)
        content_length = int(head_response.headers.get('Content-Length', 0))
        file_size_mb = content_length / (1024 * 1024)
        
        if content_length > 0:
            print(f"[{geo_accession}] File size: {file_size_mb:.2f} MB")
            
            # Size threshold warning prompt
            if file_size_mb > 100:
                user_input = input(f"[WARNING]: This file is large ({file_size_mb:.2f} MB). Do you want to download it? (y/n): ")
                if user_input.strip().lower() != 'y':
                    print(f"[{geo_accession}] Download skipped.")
                    return pd.DataFrame()
        else:
            print(f"[{geo_accession}] Could not determine file size automatically from server.")

    except requests.exceptions.RequestException as e:
        print(f"[{geo_accession}] Warning: Could not verify file size ({e}).")


    print(f"[{geo_accession}] Downloading: {file_url}")
    response = requests.get( download_url, timeout=60)
    response.raise_for_status()

    # Step 3: read into a DataFrame — try common separators
    raw = io.BytesIO(response.content)

    #print(f">>>>>>>>>> file_url: {file_url}")
    if file_url.endswith(".xls.gz") or file_url.endswith(".xlsx.gz"):
        with gzip.GzipFile(fileobj=raw) as f:
            excel_bytes = io.BytesIO(f.read())
        df = pd.read_excel(excel_bytes, index_col=0)
    elif file_url.endswith(".tar.gz"):
        with tarfile.open(fileobj=raw, mode="r:gz") as tar:
            members = tar.getnames() # Get a list of all file names inside the tar archive
            data_file_name = next((m for m in members if not m.startswith('.') and # Look for the data file you want to read (skipping hidden system files)
                                (m.endswith('.xls') or m.endswith('.xlsx') or m.endswith('.tsv') or m.endswith('.csv'))), None)
            
            if data_file_name:
                extracted_file = tar.extractfile(data_file_name)
                
                if data_file_name.endswith(('.xls', '.xlsx')):
                    excel_bytes = io.BytesIO(extracted_file.read())
                    df = pd.read_excel(excel_bytes, index_col=0)
                    
                else:
                    sep = "\t" if data_file_name.endswith(".tsv") or data_file_name.endswith(".txt") else ","
                    df = pd.read_csv(extracted_file, sep=sep, index_col=0)
            else:
                raise FileNotFoundError("No valid data file found inside the TAR archive.")
    elif file_url.endswith(".gz"):
        with gzip.GzipFile(fileobj=raw) as f:
            sep = "\t" if file_url.endswith(".tsv.gz") else ","
            df = pd.read_csv(f, sep=sep, index_col=0)
    elif file_url.endswith(".tar"):
        with tarfile.open(fileobj=raw, mode="r:") as tar: # Change mode to "r:" for uncompressed TAR files
            members = tar.getnames() # Get a list of all file names inside the tar archive
            data_file_name = next((m for m in members if not m.startswith('.') and # Look for the data file you want to read (skipping hidden system files)
                                (m.endswith('.xls') or m.endswith('.xlsx') or m.endswith('.tsv') or m.endswith('.csv'))), None)
            if data_file_name:
                extracted_file = tar.extractfile(data_file_name)
                
                if data_file_name.endswith(('.xls', '.xlsx')):
                    excel_bytes = io.BytesIO(extracted_file.read())
                    df = pd.read_excel(excel_bytes, index_col=0)
                    
                else:
                    sep = "\t" if data_file_name.endswith(".tsv") or data_file_name.endswith(".txt") else ","
                    df = pd.read_csv(extracted_file, sep=sep, index_col=0)
            else:
                raise FileNotFoundError("No valid data file found inside the TAR archive.")
    else:
        sep = "\t" if file_url.endswith(".tsv") or data_file_name.endswith(".txt") else ","
        df = pd.read_csv(raw, sep=sep, index_col=0)

    print(f"[{geo_accession}] Matrix shape: {df.shape}")
    print(f"[{geo_accession}] First few index values: {list(df.index[:5])}")
    print(f"[{geo_accession}] Columns: {list(df.columns[:5])}")

    # Step 4: Normalize index to gene symbols using the GFF reference map
    df = normalize_index_to_symbols(df, id_to_symbol, geo_accession)

    # Step 5: filter to genes of interest
    # Index may be gene names, b-numbers, or Ensembl IDs — inspect the print above
    found = df.index.intersection(genes_of_interest)
    missing = set(genes_of_interest) - set(found)
    if missing:
        print(f"[{geo_accession}] Genes not found: {missing}")

    return df.loc[found]


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
    probe_id_col: str = None,
    value_col: str = None,
    gene_symbol_col: str = None,
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
        return get_rnaseq_expression_for_genes(geo_accession, genes_of_interest, id_to_symbol)
    
def get_geo_series_metadata(gse_ids: list[str], ncbi_email: str):
    """
    Retrieves GEO metadata for multiple GSE accessions using a single API call.
    """
    metadata_results = {}
    clean_gse_uids = []

    # Prepare to fetch from NCBI. Convert GSE string to numerical UID string
    for gse in gse_ids:
        # Example: "GSE11111" -> "200011111"
        if gse.startswith("GSE"):
            numeric_part = gse[3:]                # Extract numbers (e.g., "1234")
            padded_numeric = numeric_part.zfill(6)      # Pads out with zeros (e.g., "001234")
            uid = "200" + padded_numeric                # Result is always 9 digits: "200001234"
            clean_gse_uids.append(uid)

    # Fetch records in a single bulk request
    id_string = ",".join(clean_gse_uids)
    
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "gds",
        "id": id_string,
        "retmode": "json",
        "email": ncbi_email
    }

    print(f"Fetching metadata for {len(clean_gse_uids)} GEO series...")

    try:
        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()
        raw_data = response.json().get("result", {})
        
        # Align and cache retrieved records
        for uid in clean_gse_uids:
            if uid in raw_data:
                entry = raw_data[uid]
                
                metadata_results[entry.get("accession")] = entry
            else:
                print(f"Warning: Accession {uid} not found at NCBI.")
                
    except requests.exceptions.RequestException as e:
        print(f"Network error trying to fetch metadata: {e}")
    
    print()
    return metadata_results