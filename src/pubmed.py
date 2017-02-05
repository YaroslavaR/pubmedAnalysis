import sys, os
import time
from dictionaries import create_maps, age_of_cited_work, create_article_year_by_topic_from_map, keywords_map
from settings import SEARCH_TERM_1, SEARCH_TERM_2, MAX_RETURNED_TERM_1, MAX_RETURNED_TERM_2, OUTPUT_PATH
from util import create_keywords_excel, avg, avg_list
from plotting import growth_rate, create_amount_of_article_citations_per_year_plot, create_avg_citations_age_per_year_plot
import ssl
import numpy as np

if __name__ == '__main__':

    # =============== INIT ================
    try:
        _create_unverified_https_context = ssl._create_unverified_context
    except AttributeError:
        # Legacy Python that doesn't verify HTTPS certificates by default
        pass
    else:
        # Handle target environment that doesn't support HTTPS verification
        ssl._create_default_https_context = _create_unverified_https_context
    if not os.path.isdir(OUTPUT_PATH):
        os.mkdir(OUTPUT_PATH)

    # =============== FIRST TERM ARTICLES AND CITATIONS ================
    ml_papers_map, ml_articles_links_map, ml_citations_map = create_maps(
        SEARCH_TERM_1 + ' NOT ' + SEARCH_TERM_2, MAX_RETURNED_TERM_1)

    cit_ml_dict_sort = create_article_year_by_topic_from_map(ml_citations_map)
    ml_map = create_article_year_by_topic_from_map(ml_papers_map)
    ml_age_of_cited_work = age_of_cited_work(
        ml_papers_map, ml_articles_links_map, ml_citations_map)
    ml_keywords_map = keywords_map(ml_papers_map)

    # =============== SECOND TERM ARTICLES AND CITATIONS ================

    ai_papers_map, ai_articles_links_map, ai_citations_map = create_maps(
        SEARCH_TERM_2 + ' NOT ' + SEARCH_TERM_1, MAX_RETURNED_TERM_2)

    cit_ai_dict_sort = create_article_year_by_topic_from_map(ai_citations_map)
    ai_map = create_article_year_by_topic_from_map(ai_papers_map)
    ai_age_of_cited_work = age_of_cited_work(
        ai_papers_map, ai_articles_links_map, ai_citations_map)
    ai_keywords_map = keywords_map(ai_papers_map)

    create_keywords_excel('ml_keywords.xlsx', ml_keywords_map)
    create_keywords_excel('ai_keywords.xlsx', ai_keywords_map)

    # ====================== PLOTTING ======================

    create_amount_of_article_citations_per_year_plot(
        ai_map, ml_map, cit_ai_dict_sort, cit_ml_dict_sort)
    create_avg_citations_age_per_year_plot(cit_ai_dict_sort, cit_ml_dict_sort,
                                           ai_age_of_cited_work,
                                           ml_age_of_cited_work)
    growth_rate(ai_map, 'ai')
    growth_rate(ml_map, 'ml')
