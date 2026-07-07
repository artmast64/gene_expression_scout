# Streamlit App

import io
import zipfile
from contextlib import redirect_stdout

import streamlit as st
import pandas as pd

import main # Main script

# Set tab title and icon
st.set_page_config(page_title="Gene Expression Scout", page_icon="🧬", layout="wide")

# Initialize session state keys to avoid KeyError on first load
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

# Page content
st.title("Gene Expression Scout")
st.caption("An Agentic AI-powered transcriptomics research tool")
st.markdown(
    """
    Provide the system with a biological factor and genes of interest, and the system retrieves relevant gene expression data.

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
    st.text_input("LLM API Key", key="llm_api_key", type="password")

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

st.divider()

# --- INPUTS ---
st.subheader("Inputs")

def run_program():
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
    st.markdown("*Gene symbols can be seperated by spaces, commas, or new lines*")

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
        with st.status("Running transcriptomics pipeline...", expanded=True) as status:

            # Open an in-memory string stream to catch print() statements
            log_stream = io.StringIO()

            with redirect_stdout(log_stream):

                status.write("Parsing inputs...")
                # Process inputs
                species = st.session_state.species_input
                condition = st.session_state.condition_input
                gene_list = st.session_state.gene_list_input.replace(",", " ").split() # Remove commas, then seperate on any whitespace

                model_name = st.session_state.model_name
                api_key = st.session_state.llm_api_key

                # Pull defaults using .get() to prevent crashes during top-to-bottom initialization
                go_term_grouping = st.session_state.get("go_term_grouping", "none")
                batch_size = st.session_state.get("batch_size", 50)
                min_depth = st.session_state.get("min_depth", 4)
                max_series_return = st.session_state.get("max_series_return", 50)

                # Get hardcoded variables
                gff_cache_dir, ncbi_email, category_map = main.get_params()

                # Print settings to console
                main.print_settings(species, condition, gene_list, go_term_grouping, batch_size, min_depth, max_series_return)

                taxon_id = main.identify_taxon(species)

                status.write("Retrieving Gene Ontology terms from original genes...")
                go_df, godag = main.run_retrieve_go_terms(gene_list, taxon_id, category_map, min_depth, batch_size)

                status.write("Expanding gene list via GO terms...")
                expanded_df = main.run_retrieve_genes_from_go_terms(go_df, taxon_id, godag, category_map, go_term_grouping, batch_size)

                status.write("Identifying condition keyword synonyms...")
                condition_list = main.run_retrieve_experiment_keywords(condition, model_name, api_key)
                # Run LLM step 3 first?

                status.write("Fetching Gene Expression Omnibus series metadata...")
                experiments, all_genes = main.run_retrieve_experiments(go_df, expanded_df, species, condition_list, ncbi_email, max_series_return)

                status.write("Downloading and processing GEO series expression data...")
                combined_df = main.run_retrieve_expression_levels(experiments, all_genes, species, gff_cache_dir, ncbi_email)

            # Create the ZIP Archive inside RAM
            status.write("Creating .zip file...")
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
                
                # Convert Pandas dataframes to CSV strings and pack directly into the ZIP archive
                zip_file.writestr("gene_ontology_terms.csv", go_df.to_csv(index=False))
                zip_file.writestr("expanded_genes.csv", expanded_df.to_csv(index=False))
                zip_file.writestr("geo_expression_levels.csv", combined_df.to_csv(index=False))
                
                # Write the collected terminal stdout string into a text log file inside the ZIP archive
                zip_file.writestr("pipeline_console_log.txt", log_stream.getvalue())

            # Extract raw byte contents amd commit to global session state
            st.session_state.zip_file_bytes = zip_buffer.getvalue()

            status.update(label="Analysis Pipeline Complete", state="complete", expanded=False)
            st.success("Data processed successfully! Download the data output below.")
        
            st.session_state.processing_active = False

# --- DOWNLOAD ---
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