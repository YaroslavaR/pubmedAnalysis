from collections import OrderedDict
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
    results = lookup('machine learning',2)
    elapsed_time = time.time() - start_time      
    print '====== ML - 160 - ELAPSED TIME ========='
    print elapsed_time
    id_list = results['IdList'] 
    papers = get_details_by_id(id_list)
    papers_map = create_articles_map(papers)
    #print '====== ML - PAPERS_MAP ========='
    #print papers_map
    print '====== SIZE papers_map ========='
    print(sys.getsizeof(papers_map))

    papers_year_dict = create_article_year_by_topic_map(papers)
    books_year_dict = create_book_year_by_topic_map(papers)

    print '====== SIZE papers_year_dict, books_year_dict ========='
    print(sys.getsizeof(papers_year_dict), sys.getsizeof(books_year_dict))

    print '====== SIZE papers before ========='
    print(sys.getsizeof(papers))

    papers = []
    print '====== SIZE papers after ========='
    print(sys.getsizeof(papers))

    
    #print type(papers)
    #print papers['PubmedArticle'][0]
      #print paper['MedlineCitation']['Article']['ArticleTitle']
#     print paper['MedlineCitation']['Article']['ArticleDate']    
      #print paper['MedlineCitation']['DateCreated']['Year']
      #print paper['PubmedData']['History'][0]['Year']#['PubMedPubDate']
      #title = paper['MedlineCitation']['Article']['ArticleTitle']
      #print book['BookDocument']['Book']['BookTitle']
      #print book['PubmedBookData']['History'][0]['Year']
      #print '\n'
#      print book

    #pdo = OrderedDict(sorted(papers_year_dict.items()))
    pdo = sort_map_by_keys(papers_year_dict)
    books_pdo = sort_map_by_keys(books_year_dict)
    
    start_time = time.time()
    ai_results = lookup('artificial intelligence',2)
    elapsed_time = time.time() - start_time      
    print '====== AI - 670 - ELAPSED TIME ========='
    print elapsed_time
    ai_id_list = ai_results['IdList']
    ai_papers = get_details_by_id(ai_id_list)

    ai_papers_year_dict = create_article_year_by_topic_map(ai_papers)
    ai_books_year_dict = create_book_year_by_topic_map(ai_papers)

    ai_pdo = sort_map_by_keys(ai_papers_year_dict)
    ai_books_pdo = sort_map_by_keys(ai_books_year_dict)


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
links = get_citations_ids(ai_id_list)
citations = get_details_by_id(links)
citations_map = create_articles_map(citations)
print '====== CITATIONS_MAP ========='
print citations_map
elapsed_time = time.time() - start_time
print '====== AI CITATIONS - 670 - ELAPSED TIME ========='
print elapsed_time
start_time = time.time()
links_ml = get_citations_ids(id_list)
citations_ml = get_details_by_id(links_ml)
elapsed_time = time.time() - start_time
print '====== ML CITATIONS - 160 - ELAPSED TIME ========='
print elapsed_time
#print '============= LINKS ============'
#print links

print '====== AI CITATIONS MAP - 670 - ELAPSED TIME ========='
start_time = time.time()
link_map = get_citations_ids_map(ai_id_list)
elapsed_time = time.time() - start_time
print elapsed_time
print link_map

citations_dict = create_article_year_by_topic_map(citations)
cit_dict_sort = sort_map_by_keys(citations_dict)

citations_ml_dict = create_article_year_by_topic_map(citations_ml)
cit_ml_dict_sort = sort_map_by_keys(citations_ml_dict)


plt.plot(cit_ml_dict_sort.keys(),cit_ml_dict_sort.values(), '--r', label='ML')
plt.plot(cit_dict_sort.keys(),cit_dict_sort.values(), '--c', label='AI')
plt.title("Amount of article citations per year for AI vs for ML")
plt.xlabel("Year")
plt.ylabel("Amount of articles");

plt.legend()
plt.show()

#plt.bar(range(len(pdo)), pdo.values(), align='center')
#plt.xticks(range(len(pdo)), pdo.keys())
#plt.bar(range(len(books_year_dict)), books_year_dict.values(), color='red', align='center')
#plt.xticks(range(len(books_year_dict)), books_year_dict.keys())

#plt.show()

plt.plot(ai_pdo.keys(), ai_pdo.values(), ':r', label='AI')
plt.plot(pdo.keys(),pdo.values(), '--c', label='ML')
plt.title("Amount of articles published about ML and AI per year")
plt.xlabel("Year")
plt.ylabel("Amount of articles");

plt.legend()
plt.show()

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
#plt.plot(ai_books_pdo.keys(), ai_books_pdo.values(), ':r', label='AI')
#plt.plot(books_pdo.keys(),books_pdo.values(), '--c', label='ML')
#plt.title("Amount of books published about ML and AI per year")
#plt.xlabel("Year")
#plt.ylabel("Amount of books");

#plt.legend()
#plt.show()
