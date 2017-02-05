from collections import OrderedDict
from Bio import Entrez
Entrez.email = "yaroslava.romanets@uj.edu.pl"

def get_citations_ids_map(id_list, webenv):
  print('============== ID LIST ================')
  print(id_list)
  linked = {}
  for i in range(0, len(id_list)):
    handle = Entrez.elink(dbfrom="pubmed", id=id_list[i], linkname="pubmed_pubmed_refs", webenv=webenv)
    results = Entrez.read(handle)
    print('============== RESULTS ================')
    print(results)
    handle.close()
    if len(results[0]["LinkSetDb"]) != 0:
      linked[id_list[i]] = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
      print('============== LINKED ================')
      print(linked)
      print('============== ID ================')
      print(id_list[i])
  return linked

def create_article_year_by_topic_map(query_result):
    year_by_topic = {}
    for paper in query_result['PubmedArticle']:
      year = paper['PubmedData']['History'][0]['Year']#['PubMedPubDate']
      if not year in year_by_topic:
        year_by_topic[year] = 1
      else:
        year_by_topic[year] += 1
    return year_by_topic

def create_articles_map(query_result):
    articles_map = {}
    print(query_result)
    for paper in query_result['PubmedArticle']:
      attributes_map = {}
      attributes_map['pub_year'] = paper['PubmedData']['History'][0]['Year']
      attributes_map['pub_lang'] = paper['MedlineCitation']['Article']['Language'][0]
      attributes_map['pub_type'] = str(paper['MedlineCitation']['Article']['PublicationTypeList'][0])
      attributes_map['pub_title'] = paper['MedlineCitation']['Article']['ArticleTitle']
      attributes_map['pub_keywords'] = paper['MedlineCitation']['KeywordList']
      articles_map[str(paper['MedlineCitation']['PMID'])] = attributes_map
    print('=========    MAPMAPMAP     =========')
    print(articles_map)
    return articles_map

def create_article_year_by_topic_from_map(a_map):
    year_by_topic = {}
    for i in a_map:
      year = a_map[i]['pub_year']
      if not year in year_by_topic:
        year_by_topic[year] = 1
      else:
        year_by_topic[year] += 1
    print('=========    YEAR_BY_TOPIC_MAPMAPMAP     =========')
    print(sort_map_by_keys(year_by_topic))
    return sort_map_by_keys(year_by_topic)

def create_book_year_by_topic_map(query_result):
    year_by_topic = {}
    for book in query_result['PubmedBookArticle']:
      year = book['PubmedBookData']['History'][0]['Year']#['PubMedPubDate']
      if not year in year_by_topic:
        year_by_topic[year] = 1
      else:
        year_by_topic[year] += 1
    return year_by_topic

def sort_map_by_keys(dictionary):
    d = OrderedDict(sorted(dictionary.items()))
    return d

def avg(slist):  
    print('========== DICT ============')
    print(slist)
    # for i in slist:  
    #     total = total + i
    for i in range(0,len(slist)):  
        slist[i] = float(slist[i])/len(slist)
    return slist

def avg_list(slist):  
    print('========== DICT ============')
    print(slist)
    rlist = []
    # for i in slist:  
    #     total = total + i
    for elem in slist:
      rlist.append(sum(elem)/float(len(elem)))
    slist = rlist
    print(slist)
    return slist

def age_of_cited_work(art_map, art_cit_map, cit_map):
  age = {}
  for article in art_cit_map:
    age_list = []
    print('========== ARTICLE ============')
    print(article)
    for citation in art_cit_map[article]:
      print('========== CITATION ============')
      print(citation)
      print(art_map[article]['pub_year'])
      print(cit_map[citation]['pub_year'])
      age_list.append(int(art_map[article]['pub_year']) - int(cit_map[citation]['pub_year']))
      age[art_map[article]['pub_year']] = age_list
      print('========== AGE ============')
      print(age)
  return sort_map_by_keys(age)

#def split_keys_into_chunks(id_list):
 # for i in range(0, len(id_list), 70):

from math import sqrt

def zscore(obs, pop):
  # Size of population.
  number = float(len(pop))
  # Average population value.
  avg = sum(pop) / number
  # Standard deviation of population.
  std = sqrt(sum(((c - avg) ** 2) for c in pop) / number)
  # Zscore Calculation.
  return (obs - avg) / std

  # def get_citations_ids_map(id_list):
#   #ids = ''
#   #if (not isinstance(id_list, list)):
#   print '============== ID LIST ================'
#   print id_list
#   #ids = ','.join(id_list)
#   #print '============== IDS ================'
#   #print ids
#   #else:
#    # ids = id_list
#   linked = {}
#   for i in range(0, len(id_list)):
#     handle = Entrez.elink(dbfrom="pubmed", id=id_list[i], linkname="pubmed_pubmed")
#     results = Entrez.read(handle)
#     print '============== RESULTS ================'
#     print results
#     handle.close()
#     link_year = {}
#     for link in results[0]["LinkSetDb"][0]["Link"]:
#       link_year[link["Id"]] = get_details_by_id(link["Id"])['PubmedArticle'][0]['PubmedData']['History'][0]['Year']
#     linked[id_list[i]] = link_year

#     print '============== LINKED ================'
#     print linked
#     print '============== ID ================'
#     print id_list[i]
#     #print linked
#   return linked