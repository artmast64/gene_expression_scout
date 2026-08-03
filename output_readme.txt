========================================================================
                     GENE EXPRESSION SCOUT - OUTPUT README
========================================================================

This .zip file contains the identified expression level results, raw
data, and metadata from the transcriptomics data retrieval pipeline.

------------------------------------------------------------------------
FILE CONTENTS
------------------------------------------------------------------------

1. geo_expression_levels.csv
   - The primary output matrix containing parsed and formatted gene expression 
     data collected across all matching GEO series.
   - Rows correspond to gene symbols (original input and GO-expanded genes).
   - Multi-level columns represent sample metadata (GSE Accession, GSE Title, 
     GSM Accession, GSM Sample Name).

2. expanded_genes.csv
   - Detailed breakdown of genes identified during the Gene Ontology (GO) term 
     expansion step.
   - Includes information on the source gene, GO IDs, and hit counts linking 
     the original gene list to the newly discovered genes.

3. expanded_genes_grouped.csv
   - An alternative version of expanded_genes.csv that groups together rows by
     new gene symbols.
   - Genes found in multiple GO terms have concatenated values.

3. gene_ontology_terms.csv
   - Full list of GO terms and annotations (Biological Process, Molecular 
     Function, Cellular Component) mapped to the initial seed gene list.
   - Includes info on filtering using specified GO term depth settings.

4. geo_series_links.txt
   - Text file listing all GEO Series (GSE) accessions identified and processed 
     during this execution run, including series titles and direct web links to NCBI.

5. pipeline_console_log.txt
   - Execution log recording console outputs, status checks, taxonomy lookups, 
     and pipeline settings used for this run.

6. raw_files/ (Directory)
   - Contains raw downloaded GEO series metadata (.SOFT files) and processed 
     supplementary tables organized into subfolders by GSE accession number.

------------------------------------------------------------------------
ABOUT & CITATION
------------------------------------------------------------------------
Gene Expression Scout was developed by Brady Johnson-Hill, Dr. Vijayan 
Sugumaran, and Dr. Fabia Battistuzzi at Oakland University.

Source code & documentation: https://github.com/artmast64/gene_expression_scout
Contact: bjohnsonhill@oakland.edu
========================================================================