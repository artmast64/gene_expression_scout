# Streamlit App

import io
import re
import zipfile
from contextlib import redirect_stdout

import streamlit as st
import pandas as pd

import main # Main script

# Set tab title and icon
st.set_page_config(page_title="Gene Expression Scout", page_icon="🧬", layout="wide")

# Initialize session state keys to avoid KeyError on first load
# Streamlit app variables
if "condition_input" not in st.session_state:
    st.session_state.condition_input = ""
if "species_input" not in st.session_state:
    st.session_state.species_input = ""
if "gene_list_input" not in st.session_state:
    st.session_state.gene_list_input = ""
if "processing_active" not in st.session_state:
    st.session_state.processing_active = False
if "zip_file_bytes" not in st.session_state:
    st.session_state.zip_file_bytes = None
if "expression_df_preview" not in st.session_state:
    st.session_state.expression_df_preview = None

# Pipeline data state
if "taxon_id" not in st.session_state:
    st.session_state.taxon_id = None
if "go_df" not in st.session_state:
    st.session_state.go_df = None
if "expanded_df" not in st.session_state:
    st.session_state.expanded_df = None
if "condition_list" not in st.session_state:
    st.session_state.condition_list = None
if "experiments" not in st.session_state:
    st.session_state.experiments = None
if "all_genes" not in st.session_state:
    st.session_state.all_genes = None

# File size permission state
if "approved_files" not in st.session_state:
    st.session_state.approved_files = set()
if "rejected_files" not in st.session_state:
    st.session_state.rejected_files = set()
if "pending_prompt" not in st.session_state:
    st.session_state.pending_prompt = None

# Log file state
if "pipeline_log" not in st.session_state:
    st.session_state.pipeline_log = ""


# Page content
st.title("Gene Expression Scout")
st.caption("An Agentic AI-powered transcriptomics research tool")
st.markdown(
    """
    Provide the system with a species name, biological factor and genes of interest, and the system retrieves relevant gene expression data.

    Please visit our [Github](https://github.com/artmast64/gene_expression_scout) for this project's source code and more information.

    Created by Brady Johnson-Hill at Oakland University.

    :red[NOTE: This project is still under development. Please let us know if you encounter anything unusual.]
    """
)

st.divider()

with st.sidebar:
    st.image("ges.png")
    st.markdown("""
    View the project source code on [Github](https://github.com/artmast64/gene_expression_scout).
    """)
    st.markdown('<a href="mailto:bjohnsonhill@oakland.edu">Contact us</a>', unsafe_allow_html=True)

# --- SETTINGS ---
st.subheader("Settings")

col1, col2 = st.columns(2) # Split screen into two columns

with col1:
    # LLM selection
    models_df = pd.DataFrame({
        "Available models": ["Gemini 2.5 Flash Lite", "Gemini 3.1 Flash Lite", "Claude Opus 4.8", "ChatGPT 3.5 Turbo"],
        "Model names": ["gemini-2.5-flash-lite", "gemini-3.1-flash-lite", "claude-opus-4-8", "gpt-3.5-turbo"]
    })
    model_mapping = dict(zip(models_df["Model names"], models_df["Available models"]))

    st.selectbox("LLM selection:",
                 options=models_df["Model names"],
                 format_func=lambda x: model_mapping.get(x), # Show "Available models" but pass "Model names"
                 key="model_name")
    
    # LLM API key
    st.text_input("LLM API Key", key="llm_api_key")

    st.markdown("""
    Need an API key? Follow these links for more info:
    - [Gemini](https://aistudio.google.com/app/api-keys) (Google)
    - [Claude](https://platform.claude.com/settings/keys) (Anthropic)
    - [ChatGPT](https://platform.openai.com/account/api-keys) (OpenAI)
    """)

with col2:
    # GO term grouping
    go_term_grouping_options = pd.DataFrame({
        "Options": ["none", "categories", "all"],
        "Description": ["For each unique GO term, retrieve all genes annotated with it",
                        "For each original gene, retrieve only genes annotated with the exact same combination of GO terms as that gene",
                        "For each original gene, retrieve only genes annotated with the same combination of GO terms for a GO category as that gene"]
    })
    st.radio(
        "GO term grouping",
        key = "go_term_grouping", # saved in st.session_state.go_term_grouping
        options = go_term_grouping_options["Options"]
    )

    # Batch size for API calls
    st.number_input("Batch size for API calls (default 50)", value=50, key="batch_size" )

    # Minimum GO term depth
    st.number_input("Minimum GO term depth (recommended 4+)", value=4, key="min_depth")

    # Maximum GEO series results per keyword
    st.number_input("Maximum GEO series results per keyword (default 50)", value=50, key="max_series_return")

    # Warning size for large supplementary files
    st.number_input("Warning size for large supplementary files (in MB, default 512)", value = 512, key="warn_file_size_mb")

    # Drop unmatched genes from expression matrix
    st.radio(
        "Drop unmatched genes from expression matrix (default false)",
        key = "drop_unmatched_genes",
        options = ["True", "False"],
        index = 1 # default = False
    )

st.divider()

# --- INPUTS ---
st.subheader("Inputs")

def run_program():
    # Clear log file
    st.session_state.pipeline_log = ""

    # Hard-clear output objects
    st.session_state.zip_file_bytes = None
    st.session_state.expression_df_preview = None
    
    # Clear pipeline cached dataframes
    st.session_state.taxon_id = None
    st.session_state.go_df = None
    st.session_state.godag = None
    st.session_state.expanded_df = None
    st.session_state.condition_list = None
    st.session_state.experiments = None
    st.session_state.all_genes = None

    # Reset file authorization state
    st.session_state.approved_files = set()
    st.session_state.rejected_files = set()
    st.session_state.pending_prompt = None
    
    st.session_state.processing_active = True # Sets a flag at the bottom of this script to run the program

# Define a callback function to safely reset values before widgets reload
def clear_form_callback():
    st.session_state.condition_input = ""
    st.session_state.species_input = ""
    st.session_state.gene_list_input = ""

with st.form(key="gene_list_form", clear_on_submit=False):
    
    st.text_input("Factor", key="condition_input")
    st.text_input("Species name", key="species_input")
    st.text_area("List of gene symbols (case-sensitive)", placeholder="Ex: BRCA1 TP53 TNF", key="gene_list_input")
    st.markdown("*Paste your list of gene symbols in any format. Spaces, commas, new lines, quotes, and special characters will be handled automatically.*")

    col1, col2, col3 = st.columns([1,1,10]) # Push buttons closer together
    with col1:
        submit_button = st.form_submit_button("Submit", type="primary", on_click=run_program) # type = theme color
    with col2:
        clear_button = st.form_submit_button("Clear", type="secondary", on_click=clear_form_callback)

st.divider()

# --- RUN THE PROGRAM ---
if st.session_state.processing_active:
    # Check if variables are defined
    if not st.session_state.species_input:
        st.error("Species name is required!")
        st.session_state.processing_active = False
    elif not st.session_state.condition_input:
        st.error("Factor is required!")
        st.session_state.processing_active = False
    elif not st.session_state.gene_list_input:
        st.error("List of gene symbols is required!")
        st.session_state.processing_active = False
    elif not st.session_state.llm_api_key:
        st.error("LLM API key! is required!")
        st.session_state.processing_active = False
    
    else:
        # Run the program
        with st.status("Running transcriptomics pipeline...", expanded=True) as status:

            # Open an in-memory string stream to catch print() statements
            log_stream = io.StringIO()

            with redirect_stdout(log_stream):

                status.write("Parsing inputs...")
                # Process inputs
                species = st.session_state.species_input
                condition = st.session_state.condition_input
                # Splits on commas, semicolons, pipes, whitespace (spaces, tabs, newlines), brackets, and quotes
                gene_list = [
                    gene.strip() 
                    for gene in re.split(r'[\s,;|\\/\[\]()"\']+', st.session_state.gene_list_input) 
                    if gene.strip()
                ]
                # gene_list = st.session_state.gene_list_input.replace(",", " ").split() # Remove commas, then seperate on any whitespace

                model_name = st.session_state.model_name
                api_key = st.session_state.llm_api_key

                # Pull defaults using .get() to prevent crashes during top-to-bottom initialization
                go_term_grouping = st.session_state.get("go_term_grouping", "none")
                batch_size = st.session_state.get("batch_size", 50)
                min_depth = st.session_state.get("min_depth", 4)
                max_series_return = st.session_state.get("max_series_return", 50)
                warn_file_size_mb = st.session_state.get("warn_file_size_mb", 512)
                drop_unmatched_genes = st.session_state.get("drop_unmatched_genes", False)

                # Get hardcoded variables
                gff_cache_dir, ncbi_email, category_map = main.get_params()

                # Print settings to console
                main.print_settings(species, condition, gene_list, model_name, go_term_grouping, batch_size,
                                    min_depth, max_series_return, warn_file_size_mb, drop_unmatched_genes)

                # Get taxon id
                # Check if it was cached (don't retrieve if pending prompt)
                if st.session_state.taxon_id is None and st.session_state.pending_prompt is None:
                    taxon_id = main.identify_taxon(species, ncbi_email)
                    st.session_state.taxon_id = taxon_id
                elif st.session_state.taxon_id is not None and st.session_state.pending_prompt is None:
                    taxon_id = st.session_state.taxon_id
                    print("\nRetrieved taxon ID from session state")

                # Get dataframes from session state (don't retrieve if pending prompt)
                # If the file size warning forced the app to rerun, the objects get deleted, but remain in the session state
                if st.session_state.go_df is not None and st.session_state.pending_prompt is None:
                    go_df = st.session_state.go_df
                    print("Retrieved GO terms dataframe from session state")
                if st.session_state.expanded_df is not None and st.session_state.pending_prompt is None:
                    expanded_df = st.session_state.expanded_df
                    print("Retrieved expanded genes dataframe from session state\n")
                # combined_df wouldn't be saved, so just set it to an empty dataframe
                combined_df = pd.DataFrame()

                # --- STEP 1: GO Terms ---
                if st.session_state.go_df is None:
                    status.write("Retrieving Gene Ontology terms from original genes...")
                    go_df, godag = main.run_retrieve_go_terms(gene_list, taxon_id, category_map, min_depth, batch_size)
                    st.session_state.go_df = go_df
                    st.session_state.godag = godag

                # # --- STEP 2: Expand Genes ---
                if st.session_state.go_df is not None and st.session_state.expanded_df is None:
                    status.write("Expanding gene list via GO terms...")
                    expanded_df = main.run_retrieve_genes_from_go_terms(go_df, taxon_id, godag, category_map, go_term_grouping, batch_size)
                    st.session_state.expanded_df = expanded_df

                # --- STEP 3: LLM Keywords ---
                if st.session_state.expanded_df is not None and st.session_state.condition_list is None:
                    status.write("Identifying condition keyword synonyms...")
                    condition_list = main.run_retrieve_experiment_keywords(condition, model_name, api_key)
                    st.session_state.condition_list = condition_list

                # --- STEP 4A: Search GEO Experiments ---
                if st.session_state.condition_list is not None and st.session_state.experiments is None:
                    status.write("Fetching Gene Expression Omnibus series metadata...")
                    experiments, all_genes = main.run_retrieve_experiments(gene_list, expanded_df, species, condition_list, ncbi_email, max_series_return)
                    st.session_state.experiments = experiments
                    st.session_state.all_genes = all_genes

                # --- STEP 4B: Fetch Expression Levels ---
                if st.session_state.experiments is not None:

                    # If there is a pending confirmation prompt, show the UI dialog and halt execution
                    if st.session_state.pending_prompt:
                        prompt = st.session_state.pending_prompt
                        st.warning(f"⚠️ **Large File Warning**: {prompt['series']}")
                        st.write(f"File **`{prompt['filename']}`** is **{prompt['size_mb']:.2f} MB**.")
                        
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            if st.button("✅ Download File", key="approve_btn"):
                                # Write to log file
                                print("\n--------------------------------------------")
                                print("--- User chose to download the file      ---")
                                print("--- Rerunning app.py with cached outputs ---")
                                print("--------------------------------------------\n")

                                # Track approval & clear pending state
                                st.session_state.approved_files.add(prompt['url'])
                                st.session_state.pending_prompt = None

                                # Append log file
                                st.session_state.pipeline_log += log_stream.getvalue()

                                # Make sure processing stays active and rerun
                                st.session_state.processing_active = True
                                st.rerun() # Re-runs app.py; jumps straight back to Step 4B
                                
                        with col2:
                            if st.button("❌ Skip File", key="reject_btn"):
                                # Write to log file
                                print("\n--------------------------------------------")
                                print("--- User chose to skip the file          ---")
                                print("--- Rerunning app.py with cached outputs ---")
                                print("--------------------------------------------\n")

                                # Track rejection & clear pending state
                                st.session_state.rejected_files.add(prompt['url'])
                                st.session_state.pending_prompt = None

                                # Append log file
                                st.session_state.pipeline_log += log_stream.getvalue()

                                # Make sure processing stays active and rerun
                                st.session_state.processing_active = True
                                st.rerun()

                    # Otherwise, proceed with expression matrix retrieval
                    else:
                        status.write("Downloading and processing GEO series expression data...")

                        # Pass approval/rejection state filters to expression downloader
                        combined_df, pending_prompt = main.run_retrieve_expression_levels(
                            experiments=st.session_state.experiments,
                            all_genes=st.session_state.all_genes,
                            species=species,
                            warn_file_size_mb=warn_file_size_mb,
                            drop_unmatched_genes=drop_unmatched_genes,
                            gff_cache_dir=gff_cache_dir,
                            ncbi_email=ncbi_email,
                            approved_files=st.session_state.approved_files,
                            rejected_files=st.session_state.rejected_files
                        )

                        # If the retrieval hit a large file limit, save to state and rerun
                        if pending_prompt:
                            st.session_state.pending_prompt = pending_prompt

                            # Ensure downstream outputs are cleared so nothing renders below the prompt
                            st.session_state.expression_df_preview = None
                            st.session_state.zip_file_bytes = None

                            # Append log gathered up to this point
                            st.session_state.pipeline_log += log_stream.getvalue()
                            
                            st.rerun()

                        # Process final expression data preview if available
                        if combined_df is not None and not combined_df.empty:
                            # Format for Streamlit UI Preview
                            # Flatten MultiIndex into clean string column names so PyArrow doesn't throw a mixed-type warning
                            preview_df = combined_df.copy()
                            preview_df.columns = [f"{col[0]} | {col[2]} ({col[3]})" if isinstance(col, tuple) else str(col) for col in preview_df.columns]
                            st.session_state.expression_df_preview = preview_df.reset_index()
                        else:
                            st.error("No expression matrix data could be collected.")

            # Append new print statements from this run to the persistent log session
            st.session_state.pipeline_log += log_stream.getvalue()

            # Only create ZIP and complete Status if no prompt is pending
            if not st.session_state.pending_prompt:
                # Create the ZIP Archive inside RAM
                status.write("Creating .zip file...")

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:

                    # Convert Pandas dataframes to CSV strings and pack directly into the ZIP archive
                    zip_file.writestr("gene_ontology_terms.csv", go_df.to_csv(index=False))
                    zip_file.writestr("expanded_genes.csv", expanded_df.to_csv(index=False))

                    # Write expression matrix to CSV
                    # Keep index=True (Gene_Symbol) and header=True (MultiIndex Metadata Rows)
                    # This exports multi-row headers directly without adding numerical row indices (0, 1, 2...)
                    zip_file.writestr("geo_expression_levels.csv", combined_df.to_csv(index=True, header=True).encode('utf-8'))
                    
                    # Write the collected terminal stdout string into a text log file inside the ZIP archive
                    zip_file.writestr("pipeline_console_log.txt", st.session_state.pipeline_log)

                # Extract raw byte contents amd commit to global session state
                st.session_state.zip_file_bytes = zip_buffer.getvalue()

                status.update(label="Analysis Pipeline Complete", state="complete", expanded=False)
                st.success("Data processed successfully! Download the data output below.")
            
                # Mark processing as finished so preview & download show up
                st.session_state.processing_active = False

# --- DOWNLOAD & PREVIEW TABLE---
# Only display outputs if a run finished successfully AND no file approval prompt is waiting
if not st.session_state.processing_active and st.session_state.pending_prompt is None:

    # Download Button
    # This button is placed outside of the form and is persistent after generation completes
    if st.session_state.zip_file_bytes is not None:
        st.subheader("📦 Export Files")
        st.download_button(
            label="Download Analysis ZIP Archive",
            data=st.session_state.zip_file_bytes, # Reads directly out of RAM memory cache
            file_name="gene_expression_scout_output.zip",
            mime="application/zip",
            type="primary"
        )

    # Preview output_expression_matrix.csv in the Streamlit UI
    # Hide the column names row so it looks clean
    if st.session_state.expression_df_preview is not None:
        st.subheader("Expression Table Preview")
        st.dataframe(st.session_state.expression_df_preview, hide_index=True)