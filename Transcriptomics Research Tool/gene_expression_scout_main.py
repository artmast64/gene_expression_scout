# Transcriptomics Research Tool
# Version 1.0
# Created by Brady Johnson-Hill

from openai import OpenAI
import streamlit as st
import pandas as pd

## https://github.com/streamlit/llm-examples/blob/main/Chatbot.py

# Set tab title and icon
st.set_page_config(page_title="Gene Expression Scout", page_icon="🧬")

st.title("Gene Expression Scout")
st.caption("An Agentic AI-powered transcriptomics research tool")
st.markdown(
    """
    Provide the system with a query about biological process and genes of interest, and the system retrieves and evaluates relevant gene expression data.

    Please visit our [Github](https://github.com/artmast64/gene_expression_scout) for this project's source code and more information.

    Created by Brady Johnson-Hill at Oakland University.

    :red[NOTE: This project is still under development. Please let us know if you encounter anything unusual.]
    """
)

models_df = pd.DataFrame({
    "Available models": ["Gemini 3.1", "Gemini 2.5 Flash", "Claude", "ChatGPT 3.5 Turbo"]
})

llm_option = st.selectbox(
    f"LLM selection:",
    models_df["Available models"])

with st.sidebar:
    llm_api_key = st.text_input("LLM API Key", key="chatbot_api_key", type="password")

    st.markdown("""
    Need an API key? Follow these links for more info:
    - [Gemini](https://aistudio.google.com/app/api-keys) (Google)
    - [Claude](https://platform.claude.com/settings/keys) (Anthropic)
    - [ChatGPT](https://platform.openai.com/account/api-keys) (OpenAI)
    """)
    st.divider()
    st.markdown("""
    View the project source code on [Github](https://github.com/artmast64/gene_expression_scout).
    """)
    st.markdown('<a href="mailto:bjohnsonhill@oakland.edu">Contact us</a>', unsafe_allow_html=True)

st.divider()

if "messages" not in st.session_state:
    st.session_state["messages"] = [{"role": "assistant", "content": "How can I help you?"}]

for msg in st.session_state.messages:
    st.chat_message(msg["role"]).write(msg["content"])

if prompt := st.chat_input():
    if not llm_api_key:
        st.info("Please add an API key to continue.")
        st.stop()

    client = OpenAI(api_key=llm_api_key)
    st.session_state.messages.append({"role": "user", "content": prompt})
    st.chat_message("user").write(prompt)
    response = client.chat.completions.create(model="gpt-3.5-turbo", messages=st.session_state.messages)
    msg = response.choices[0].message.content
    st.session_state.messages.append({"role": "assistant", "content": msg})
    st.chat_message("assistant").write(msg)