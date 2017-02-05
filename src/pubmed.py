import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import sys
from Bio import Entrez
import time
from queries import lookup, get_details_by_id, get_citations_ids
from dictionaries import get_citations_ids_map, avg, avg_list, age_of_cited_work, \
    create_articles_map, create_book_year_by_topic_map, sort_map_by_keys, zscore, create_article_year_by_topic_from_map
from settings import ENTREZ_EMAIL, SEARCH_TERM_1, SEARCH_TERM_2, MAX_RETURNED_TERM_1, MAX_RETURNED_TERM_2
import ssl
import numpy as np


if __name__ == '__main__':
    Entrez.email = ENTREZ_EMAIL


    # =============== FIRST TERM ARTICLES AND CITATIONS ================
    start_time = time.time()
    try:
      _create_unverified_https_context = ssl._create_unverified_context
    except AttributeError:
      # Legacy Python that doesn't verify HTTPS certificates by default
      pass
    else:
      # Handle target environment that doesn't support HTTPS verification
      ssl._create_default_https_context = _create_unverified_https_context

    results = lookup(SEARCH_TERM_1 + ' NOT '+ SEARCH_TERM_2, MAX_RETURNED_TERM_1)
    elapsed_time = time.time() - start_time
    print('====== ML - MESH RESULTS =========')
    print(results)      
    print('====== ML - 160 - ELAPSED TIME =========')
    print(elapsed_time)
    ml_id_list = results['IdList']
    ml_webenv = results["WebEnv"]

    ml_papers = get_details_by_id(ml_id_list, ml_webenv)

    ml_papers_map = create_articles_map(ml_papers)
    print('====== KEYWORDS =========')
    print(ml_papers_map[str(ml_id_list[0])]['pub_keywords'])

    print('====== SIZE papers_map =========')
    print(sys.getsizeof(ml_papers_map))

    start_time = time.time()
    links_ml = get_citations_ids(ml_id_list, ml_webenv)
    ml_articles_links_map = get_citations_ids_map(ml_id_list, ml_webenv)
    citations_ml = get_details_by_id(links_ml, ml_webenv)
    ml_citations_map = create_articles_map(citations_ml)  
    elapsed_time = time.time() - start_time
    print('====== ML CITATIONS - 160 - ELAPSED TIME =========')
    print(ml_articles_links_map)
    print(elapsed_time)

    cit_ml_dict_sort = create_article_year_by_topic_from_map(ml_citations_map)
    ml_map = create_article_year_by_topic_from_map(ml_papers_map)
    ml_age_of_cited_work = age_of_cited_work(ml_papers_map, ml_articles_links_map, ml_citations_map)


    # =============== SECOND TERM ARTICLES AND CITATIONS ================
    start_time = time.time()
    print('====== AI - MESH LOOKUP STRING =========')
    lookup_string_ai = SEARCH_TERM_2 + ' NOT '+ SEARCH_TERM_1
    print(lookup_string_ai)
    ai_results = lookup(lookup_string_ai, MAX_RETURNED_TERM_2)
    print('====== AI - MESH RESULTS =========')
    print(ai_results)
    elapsed_time = time.time() - start_time      
    print('====== AI - 670 - ELAPSED TIME =========')
    print(elapsed_time)
    ai_id_list = ai_results['IdList']
    ai_webenv = ai_results["WebEnv"]

    ai_papers = get_details_by_id(ai_id_list, ai_webenv)
    ai_papers_map = create_articles_map(ai_papers)    

    start_time = time.time()
    links_ai = get_citations_ids(ai_id_list, ai_webenv)
    ai_articles_links_map = get_citations_ids_map(ai_id_list, ai_webenv)
    citations_ai = get_details_by_id(links_ai, ai_webenv)
    ai_citations_map = create_articles_map(citations_ai)
    elapsed_time = time.time() - start_time
    print('====== AI CITATIONS - 670 - ELAPSED TIME =========')
    print(ai_articles_links_map)
    print(elapsed_time)

    cit_ai_dict_sort = create_article_year_by_topic_from_map(ai_citations_map)
    ai_map = create_article_year_by_topic_from_map(ai_papers_map)

    print('=============== ARTICLESSS =================')
    ai_age_of_cited_work = age_of_cited_work(ai_papers_map, ai_articles_links_map, ai_citations_map)


    # ====================== PLOTTING ======================

    print('=============== PLOT KEYS =================')
    print(list(ai_map.keys()))
    print(list(cit_ai_dict_sort.keys()))
    print(list(ml_map.keys()))
    print(list(cit_ml_dict_sort.keys()))
    print('=============== PLOT VALUES =================')
    print(list(ai_map.values()))
    print(list(cit_ai_dict_sort.values()))
    print(list(ml_map.values()))
    print(list(cit_ml_dict_sort.values()))
    plt.figure(1, figsize=(10,5))
    plt.subplot(211).xaxis.set_major_formatter(tick.FormatStrFormatter('%0.0f'))
    plt.plot(list(ai_map.keys()),list(ai_map.values()), ':r', label=SEARCH_TERM_2)
    plt.plot(list(ml_map.keys()),list(ml_map.values()), '--c', label=SEARCH_TERM_1)
    plt.title("Amount of articles published about " +  SEARCH_TERM_2 + " and " + SEARCH_TERM_1 + " per year")
    plt.xlabel("Year")
    plt.ylabel("Amount of articles");

    plt.legend(loc=2)
    #plt.show()

    plt.subplot(212).xaxis.set_major_formatter(tick.FormatStrFormatter('%0.0f'))
    plt.plot(list(cit_ai_dict_sort.keys()),list(cit_ai_dict_sort.values()), ':r', label=SEARCH_TERM_2)
    plt.plot(list(cit_ml_dict_sort.keys()),list(cit_ml_dict_sort.values()), '--c', label=SEARCH_TERM_1)
    plt.title("Amount of article citations per year for " +  SEARCH_TERM_2 + " vs for " + SEARCH_TERM_1)
    plt.xlabel("Year")
    plt.ylabel("Amount of articles");

    plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)

    plt.legend(loc=2)
    plt.show()

    # myList[:] = [x / myInt for x in myList]
    plt.figure(2, figsize=(10,5))
    plt.subplot(211).xaxis.set_major_formatter(tick.FormatStrFormatter('%0.0f'))
    plt.plot(list(cit_ai_dict_sort.keys()),avg(list(cit_ai_dict_sort.values())),':r', label=SEARCH_TERM_2)
    plt.plot(list(cit_ml_dict_sort.keys()),avg(list(cit_ml_dict_sort.values())),'--c', label=SEARCH_TERM_1)
    plt.title("Average amount of article citations per year for " +  SEARCH_TERM_2 + " vs for " + SEARCH_TERM_1)
    plt.xlabel("Year")
    plt.ylabel("Average amount of articles");
    plt.legend(loc=2)

    plt.subplot(212).xaxis.set_major_formatter(tick.FormatStrFormatter('%0.0f'))
    print(list(ai_age_of_cited_work.keys()))
    print(avg_list(ai_age_of_cited_work.values()))
    ind = np.arange(len(ai_age_of_cited_work))
    plt.bar(ind-0.15, avg_list(ai_age_of_cited_work.values()),0.3, color='r', label=SEARCH_TERM_2)
    plt.bar(ind+0.15, avg_list(ml_age_of_cited_work.values()), 0.3, color='c', label=SEARCH_TERM_1)
    #, ':r', label='AI')
    plt.xticks(ind, list(ai_age_of_cited_work.keys()))
    plt.title("Average age cited works per year for " +  SEARCH_TERM_2 + " vs for " + SEARCH_TERM_1)
    plt.xlabel("Year")
    plt.ylabel("Age of articles");

    plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)

    plt.legend(loc=2)
    plt.show()

    
    
    growth_rate = np.exp(np.diff(np.log(list(ai_map.values())))) - 1
    print(growth_rate)
    x = np.array(list(ai_map.keys()))
    print(x)
    growth_rate = np.insert(growth_rate,0,0)
    print(growth_rate)
    print(len(growth_rate))
    print(len(x))

    plt.figure(3)
    plt.subplot(111).xaxis.set_major_formatter(tick.FormatStrFormatter('%0.0f'))
    plt.plot(x,growth_rate)
    plt.show()

# papers_map - articles map #'26466993': {'pub_type': 'Comparative Study', 'pub_title': 'Prediction of delayed graft function after kidney transplantation: comparison between logistic regression and machine learning methods.', 'pub_year': '2015', 'pub_lang': 'eng'}
# link_map - articleslinks map #'26125433': ['26125433', '15718595', '15847374', '15668886', '19540562', '12857667', '20164018', '22463804', '25927149', '10165328', '21965283', '17129792', '15217254', '3606341', '8526559', '22970841', '22063731', '15901437', '9640709', '15503664', '23171903', '10097629', '21284440', '24075599', '23171902', '23672084', '11747236', '25415180', '23171886', '22341607', '20129646', '15845147', '22925804', '17868776', '21937077', '19552938', '21489193', '19501988', '21366463', '9593340', '2403974', '12097174', '21309285', '23372984', '19552937', '11857500', '18035694', '8526555', '9308736', '24091812', '19596461', '21477260', '24840034', '8526556', '20081805', '23237667', '22925786', '23402499', '2488697', '9987448', '16445248', '22429209', '15875727', '8526554', '10840835', '8526558', '21903681', '9013826', '22558887', '20854068', '3606340', '8316495', '11710693', '22429212', '1297434', '25132434', '10674035', '2929514', '7988115', '22498580', '15893454', '18673501', '11273201', '12463845', '11831620', '27563491', '26958230', '27898974', '27762389', '27898975', '26551779', '26955694', '26993572', '28087372', '27441580', '26103652', '26632567', '23745154', '23115580', '27769367', '28080147', '26816261', '27194850']
# citations_map - linksdetails map
# find average age of cited work
# papers_map[link_map['pub_year']] : citations_map[link_map.values()['pub_year']]
# link_map.values() -> citations_map[link_map.values()['pub_year']]
#   year_map ={}
#   if not year in year_map:
#     year_map[year] = 1
#   else:
#     year_map[year] += 1




