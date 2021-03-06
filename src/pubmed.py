import sys, os
import time
import numpy as np
import logging
from dictionaries import create_maps, age_of_cited_work, create_article_year_by_topic_from_map, keywords_map
from settings import SEARCH_TERM_1, SEARCH_TERM_2, MAX_RETURNED_TERM_1, MAX_RETURNED_TERM_2, OUTPUT_PATH
from util import create_keywords_excel, create_unverified_context
from plotting import growth_rate, create_amount_of_article_citations_per_year_plot, create_avg_citations_age_per_year_plot, percentage_per_year_pie
""" 
Program for comparing trends and applications of two search terms in Bioinformatics terms.
Generates 6 png files (8 plots) and 2 xlxs files containing keywords related to search terms.

Usage:
    python3 src/pubmed.py
    or
    python src/pubmed.py

Make sure to set relevant search parameters in src/settings.py file
    ENTREZ_EMAIL - your email (used to query Entrez)
    SEARCH_TERM_1 - first term to query Entrez for 
                    (add [MeSH Terms] to the end of the term to have Entrez find aliases for term 
                    and search for those as well - i.e. 'artificial intelligence[MeSH Terms]') 
    MAX_RETURNED_TERM_1 - maximum amount of returned records for first term
    SEARCH_TERM_2 - second term to query Entrez for
    MAX_RETURNED_TERM_2 - maximum amount of returned records for second term
    OUTPUT_PATH - path to store .png plots and .xlxs keywords files in
"""

if __name__ == '__main__':
    """ Initializing logger """
    logging.basicConfig(
        filename='logs/pubmedAnalysis.log',
        format='%(asctime)s - [%(levelname)s] - %(message)s',
        level=logging.DEBUG)

    # =============== INIT ================
    """ Setting unverified context - do not verify HTTPS certificates """
    create_unverified_context()
    """ Creating output directory """
    if not os.path.isdir(OUTPUT_PATH):
        logging.info('Created ' + OUTPUT_PATH + ' directory')
        os.mkdir(OUTPUT_PATH)

    # =============== FIRST TERM ARTICLES AND CITATIONS ================
    logging.info('Searching for ' + SEARCH_TERM_1 + ' with max_returned=' +
                 str(MAX_RETURNED_TERM_1))
    ml_papers_map, ml_articles_links_map, ml_citations_map = create_maps(
        SEARCH_TERM_1, MAX_RETURNED_TERM_1)  # + ' NOT ' + SEARCH_TERM_2

    cit_ml_dict_sort = create_article_year_by_topic_from_map(ml_citations_map)
    ml_map = create_article_year_by_topic_from_map(ml_papers_map)
    ml_age_of_cited_work = age_of_cited_work(
        ml_papers_map, ml_articles_links_map, ml_citations_map)
    ml_keywords_map = keywords_map(ml_papers_map)

    # =============== SECOND TERM ARTICLES AND CITATIONS ================
    logging.info('Searching for ' + SEARCH_TERM_2 + ' with max_returned=' +
                 str(MAX_RETURNED_TERM_2))
    ai_papers_map, ai_articles_links_map, ai_citations_map = create_maps(
        SEARCH_TERM_2, MAX_RETURNED_TERM_2)  # + ' NOT ' + SEARCH_TERM_1

    cit_ai_dict_sort = create_article_year_by_topic_from_map(ai_citations_map)
    ai_map = create_article_year_by_topic_from_map(ai_papers_map)
    ai_age_of_cited_work = age_of_cited_work(
        ai_papers_map, ai_articles_links_map, ai_citations_map)
    ai_keywords_map = keywords_map(ai_papers_map)

    logging.info('Saving keywords to xslsx')
    create_keywords_excel(SEARCH_TERM_1 + '_keywords.xlsx', ml_keywords_map)
    create_keywords_excel(SEARCH_TERM_2 + '_keywords.xlsx', ai_keywords_map)

    # ====================== PLOTTING ======================
    create_amount_of_article_citations_per_year_plot(
        ai_map, ml_map, cit_ai_dict_sort, cit_ml_dict_sort)
    create_avg_citations_age_per_year_plot(ai_age_of_cited_work,
                                           ml_age_of_cited_work)
    percentage_per_year_pie(ai_map, SEARCH_TERM_2)
    percentage_per_year_pie(ml_map, SEARCH_TERM_1)
    growth_rate(ai_map, SEARCH_TERM_2)
    growth_rate(ml_map, SEARCH_TERM_1)
