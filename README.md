# Pubmed Analysis 
Program for comparing trends and applications of two search terms in Bioinformatics domain.
Generates 6 png files (8 plots) and 2 xlsx files containing keywords related to specified search terms.

## Usage:
    python3 src/pubmed.py

    or

    python src/pubmed.py

## Make sure to set relevant search parameters in **src/settings.py** file
Name | Default value | Description
------------ | ------------- | -------------
ENTREZ_EMAIL | 'example@student.uj.edu.pl' | your email (used to query Entrez)
SEARCH_TERM_1 | 'deep learning' | first term to query Entrez for
MAX_RETURNED_TERM_1 | 30 | maximum amount of returned records for first term
SEARCH_TERM_2 | 'artificial intelligence' | second term to query Entrez for
MAX_RETURNED_TERM_2 | 60 | maximum amount of returned records for second term
OUTPUT_PATH | 'output' | path to store .png plots and .xlxs keywords files in

**Tip:** add *[MeSH Terms]* to the end of the term to have Entrez find aliases for term and search for those as well - i.e. 'artificial intelligence[MeSH Terms]')
