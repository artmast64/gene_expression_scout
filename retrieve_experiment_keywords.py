# Step 3 Functions

import os
from typing import List

from pydantic import BaseModel, Field
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph.graph import StateGraph, START, END

def retrieve_exp_keys(keyword: str, model_name: str, api_key: str):
    # Define the desired output structure using Pydantic
    class GeoKeywords(BaseModel):
        expanded_keywords: List[str] = Field(
            description="A list of specific, highly relevant search terms, synonyms, and variations for GEO."
        )

    # Define the Graph State
    class SearchState(BaseModel):
        user_keyword: str
        search_terms: List[str] = []

    # Initialize Gemini with Structured Output
    os.environ["GOOGLE_API_KEY"] = api_key.strip()
    llm = ChatGoogleGenerativeAI(model=model_name,
                                google_api_key=api_key.strip(),
                                temperature=0.1)
    structured_llm = llm.with_structured_output(GeoKeywords)

    # Define the Agent Node
    def generate_geo_terms(state: SearchState):
        prompt = f"""
        You are an expert bioinformatics assistant specializing in transcriptomics data retrieval.
        The user wants to find data series in the NCBI Gene Expression Omnibus (GEO).
        
        Original Keyword: "{state.user_keyword}"
        
        Generate a comprehensive list of expanded search terms. Include:
        - Standard medical synonyms and acronyms
        - Common transcriptomics variations (e.g., cell line names, model organisms equivalents)

        The list of expanded search terms should include around 10 to 15 results.
        
        Keep the terms precise for NCBI GEO query optimization. Do not include organism names. Do not include boolean operators like AND/OR. 
        """
        # Print the prompt to the log file
        print(f"Model name: {model_name}")
        print("Prompt:")
        print(prompt)
        
        response = structured_llm.invoke(prompt)
        return {"search_terms": response.expanded_keywords}

    # Build and Compile the Graph
    builder = StateGraph(SearchState)
    builder.add_node("generator", generate_geo_terms)

    builder.add_edge(START, "generator")
    builder.add_edge("generator", END)

    geo_keyword_graph = builder.compile()

    # Execute the System
    initial_state = {"user_keyword": keyword}
    print("\nLLM responding...")
    result = geo_keyword_graph.invoke(initial_state)

    print(f"Original: {result["user_keyword"]}")
    print("Expanded Search Terms:")
    print(result["search_terms"])

    return result["search_terms"]