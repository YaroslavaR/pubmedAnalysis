from collections import OrderedDict
from Bio import Entrez
from queries import lookup, get_details_by_id, get_citations_ids
from settings import ENTREZ_EMAIL
Entrez.email = ENTREZ_EMAIL

""" Stages:
        1. lookup(term, max_returned) - search for term in Entrez, return list of Ids (IdList)
        2. save IdList and WebEnv
        3. get_details_by_id(IdList, WebEnv) - search for detailed information for each id from IdList, 
                        return Entrez Parse object
        4. create_articles_map(get_details_by_id(IdList, WebEnv)) - create articles dictionary, return dictionary
        5. get_citations_ids(IdList, WebEnv) - get ids of articles referenced by articles from IdList, return IdList
        6. get_citations_ids_map(IdList, WebEnv) - create dictionary of articles mapped to articles they reference, 
                        return dictionary
        7. get_details_by_id(get_citation_ids(IdList, WebEnv), WebEnv) - search for detailed information 
                        for each id from IdList of citations, return Entrez Parse object
        8. create_articles_map(get_details_by_id(get_citation_ids(IdList, WebEnv), WebEnv)) - create citations dictionary, 
                        return dictionary
"""

"""         
        9. create_article_year_by_topic_from_map(create_articles_map(get_details_by_id(IdList, WebEnv)))
        10. create_article_year_by_topic_from_map(create_articles_map(get_details_by_id(get_citation_ids(IdList, WebEnv), WebEnv)))
        11. age_of_cited_work(create_articles_map(get_details_by_id(IdList, WebEnv)), get_citation_ids_map(IdList, WebEnv), create_articles_map(get_details_by_id(get_citation_ids(IdList, WebEnv), WebEnv)))
"""

def create_maps(term, max_returned):
  results = lookup(term, max_returned)
  id_list = results['IdList']
  webenv = results["WebEnv"]
  articles = get_details_by_id(id_list, webenv)
  articles_dict = create_articles_map(articles)
  cit_id_list = get_citations_ids(id_list, webenv)
  articles_links_dict = get_citations_ids_map(id_list, webenv)
  citations = get_details_by_id(cit_id_list, webenv)
  citations_dict = create_articles_map(citations)
  return articles_dict, articles_links_dict, citations_dict

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
      if citation in cit_map:
        print(cit_map[citation]['pub_year'])
        age_list.append(int(art_map[article]['pub_year']) - int(cit_map[citation]['pub_year']))
        age[art_map[article]['pub_year']] = age_list
        print('========== AGE ============')
        print(age)
  return sort_map_by_keys(age)

def keywords_map(articles_map):
  keywords_map = {}
  for article in articles_map:
    if len(articles_map[article]['pub_keywords']) != 0:
      keyword_list = articles_map[article]['pub_keywords'][0] # for each element in articles_map[article]['pub_keywords'] list get keyword
      print(keyword_list)
      for keyword in keyword_list:
        if not keyword.lower() in keywords_map:
          keywords_map[str(keyword).lower()] = 1
        else:
          keywords_map[str(keyword).lower()] += 1
  return keywords_map

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