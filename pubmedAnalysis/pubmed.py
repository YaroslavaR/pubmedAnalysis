import matplotlib.pyplot as plt
import sys
from Bio import Entrez
Entrez.email = "yaroslava.romanets@uj.edu.pl"
import time

from queries import lookup, get_details_by_id, get_citations_ids
from dictionaries import get_citations_ids_map, create_article_year_by_topic_map, \
    create_articles_map, create_book_year_by_topic_map, sort_map_by_keys, zscore


if __name__ == '__main__':
    start_time = time.time()
    results = lookup('machine learning',160)
    elapsed_time = time.time() - start_time      
    print '====== ML - 160 - ELAPSED TIME ========='
    print elapsed_time
    ml_id_list = results['IdList']
    ml_webenv = results["WebEnv"]
    #ml_query_key = results["QueryKey"] 
    ml_papers = get_details_by_id(ml_id_list, ml_webenv)#, ml_query_key)
    ml_papers_map = create_articles_map(ml_papers)
    #print '====== ML - PAPERS_MAP ========='
    #print papers_map
    print '====== SIZE papers_map ========='
    print(sys.getsizeof(ml_papers_map))

    ml_papers_year_dict = create_article_year_by_topic_map(ml_papers)
    #ml_books_year_dict = create_book_year_by_topic_map(papers)

    print '====== SIZE papers_year_dict, books_year_dict ========='
    print(sys.getsizeof(ml_papers_year_dict))#, sys.getsizeof(ml_books_year_dict))

    print '====== SIZE papers before ========='
    print(sys.getsizeof(ml_papers))

    ml_papers = []
    print '====== SIZE papers after ========='
    print(sys.getsizeof(ml_papers))

    ml_pdo = sort_map_by_keys(ml_papers_year_dict)
    #books_pdo = sort_map_by_keys(books_year_dict)
    
    start_time = time.time()
    ai_results = lookup('artificial intelligence',700)
    elapsed_time = time.time() - start_time      
    print '====== AI - 670 - ELAPSED TIME ========='
    print elapsed_time
    ai_id_list = ai_results['IdList']
    ai_webenv = ai_results["WebEnv"]
    #ai_query_key = ai_results["QueryKey"] 
    ai_papers = get_details_by_id(ai_id_list, ai_webenv)#, ai_query_key)

    ai_papers_year_dict = create_article_year_by_topic_map(ai_papers)
    #ai_books_year_dict = create_book_year_by_topic_map(ai_papers)

    ai_pdo = sort_map_by_keys(ai_papers_year_dict)
    #ai_books_pdo = sort_map_by_keys(ai_books_year_dict)


#print '========== ML PubmedArticle ======='
#for year in pdo:
 # print year
  #print pdo[year]

#print '========== ML PubmedBookArticle ======='

#for year in books_year_dict:
 # print year
 # print books_year_dict[year]


#print '========== AI PubmedArticle ======='
#for year in ai_pdo:
 # print year
  #print ai_pdo[year]

#print '========== AI PubmedBookArticle ======='
#for year in ai_books_pdo:
 # print year
  #print ai_books_pdo[year]

start_time = time.time()
#links_ai = get_citations_ids(ai_id_list)
links_ai = get_citations_ids(ai_id_list, ai_webenv)#, ai_query_key)
#ai_cit_webenv = links_ai["WebEnv"]
#ml_query_key = results["QueryKey"]
#citations_ai = get_citation_details_by_id(links_ai)
citations_ai = get_details_by_id(links_ai, ai_webenv)#, ai_query_key)
# citations_map = create_articles_map(citations_ai)
# print '====== CITATIONS_MAP ========='
# print citations_map
elapsed_time = time.time() - start_time
print '====== AI CITATIONS - 670 - ELAPSED TIME ========='
print elapsed_time

start_time = time.time()
#links_ml = get_citations_ids(ml_id_list)
links_ml = get_citations_ids(ml_id_list, ml_webenv)#, ml_query_key)
#citations_ml = get_citation_details_by_id(links_ml)
citations_ml = get_details_by_id(links_ml, ml_webenv)#, ml_query_key)
elapsed_time = time.time() - start_time
print '====== ML CITATIONS - 160 - ELAPSED TIME ========='
print elapsed_time
#print '============= LINKS ============'
#print links

# print '====== AI CITATIONS MAP - 670 - ELAPSED TIME ========='
# start_time = time.time()
# link_map = get_citations_ids_map(ai_id_list)
# elapsed_time = time.time() - start_time
# print elapsed_time
# print link_map

citations_ai_dict = create_article_year_by_topic_map(citations_ai)
cit_ai_dict_sort = sort_map_by_keys(citations_ai_dict)

citations_ml_dict = create_article_year_by_topic_map(citations_ml)
cit_ml_dict_sort = sort_map_by_keys(citations_ml_dict)

plt.figure(1)
plt.subplot(211)
plt.plot(ai_pdo.keys(), ai_pdo.values(), ':r', label='AI')
plt.plot(ml_pdo.keys(),ml_pdo.values(), '--c', label='ML')
plt.title("Amount of articles published about ML and AI per year")
plt.xlabel("Year")
plt.ylabel("Amount of articles");

plt.legend(loc=2)
#plt.show()

plt.subplot(212)
plt.plot(cit_ai_dict_sort.keys(),cit_ai_dict_sort.values(), ':r', label='AI')
plt.plot(cit_ml_dict_sort.keys(),cit_ml_dict_sort.values(), '--c', label='ML')
plt.title("Amount of article citations per year for AI vs for ML")
plt.xlabel("Year")
plt.ylabel("Amount of articles");

plt.legend(loc=2)
plt.show()

#plt.bar(range(len(pdo)), pdo.values(), align='center')
#plt.xticks(range(len(pdo)), pdo.keys())
#plt.bar(range(len(books_year_dict)), books_year_dict.values(), color='red', align='center')
#plt.xticks(range(len(books_year_dict)), books_year_dict.keys())

#plt.show()



import numpy as np

growth_rate = np.exp(np.diff(np.log(ai_pdo.values()))) - 1
print growth_rate
x = np.array(ai_pdo.keys())
print x
growth_rate = np.insert(growth_rate,0,0)
print growth_rate
print len(growth_rate)
print len(x)


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




#plt.plot(ai_books_pdo.keys(), ai_books_pdo.values(), ':r', label='AI')
#plt.plot(books_pdo.keys(),books_pdo.values(), '--c', label='ML')
#plt.title("Amount of books published about ML and AI per year")
#plt.xlabel("Year")
#plt.ylabel("Amount of books");

#plt.legend()
#plt.show()
