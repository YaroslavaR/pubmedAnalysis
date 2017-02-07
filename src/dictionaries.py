import logging
from collections import OrderedDict
from Bio import Entrez
from queries import lookup, get_details_by_id, get_citations_ids
from util import create_unverified_context
from settings import ENTREZ_EMAIL
Entrez.email = ENTREZ_EMAIL


def create_maps(term, max_returned):
    """ 
    Execute Entrez queries and create dictionaries, return three dictionaries 

    Method setting in motion most of the program's logic. Returns three
    dictionaries: articles_dict - a dictionary containing detailed data about
    articles; citations_dict - a dictionary containing detailed data about
    citations; articles_links_dict - a dictionary representing the mappings of
    article to the citations it references.
    """
    results = lookup(term, max_returned)
    id_list = results['IdList']

    articles = get_details_by_id(id_list)
    articles_dict = create_articles_map(articles)
    cit_id_list = get_citations_ids(id_list)
    articles_links_dict = get_citations_ids_map(id_list)
    citations = get_details_by_id(cit_id_list)
    citations_dict = create_articles_map(citations)
    return articles_dict, articles_links_dict, citations_dict


def get_citations_ids_map(id_list):
    """ 
    Create dictionary of articles mapped to articles they reference, 
                        return dictionary 

    Returns a dictionary representing the mappings of
    article to the citations it references.
    """
    create_unverified_context()
    logging.debug('============== IN get_citations_ids_map: ================')
    logging.debug('============== ID LIST: ================')
    logging.debug(id_list)
    linked = {}
    for i in range(0, len(id_list)):
        handle = Entrez.elink(
            dbfrom="pubmed", id=id_list[i], linkname="pubmed_pubmed_refs")
        results = Entrez.read(handle)
        logging.debug('============== RESULTS: ================')
        logging.debug(results)
        handle.close()
        if len(results[0]["LinkSetDb"]) != 0:
            linked[id_list[i]] = [
                link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]
            ]
            logging.debug('============== LINKED ARTICLES: ================')
            logging.debug(linked)
            logging.debug('============== ARTICLE ID: ================')
            logging.debug(id_list[i])
    return linked


def create_articles_map(query_result):
    """ 
    Create articles dictionary, return dictionary 

    Returns a dictionary containing detailed data about query_result - 
    Entrez Parser object parsed to a regular Python dictionary.
    """
    logging.debug('============== IN create_articles_map: ================')
    articles_map = {}
    logging.debug(query_result)
    for paper in query_result['PubmedArticle']:
        attributes_map = {}
        attributes_map['pub_year'] = paper['PubmedData']['History'][0]['Year']
        attributes_map['pub_lang'] = paper['MedlineCitation']['Article'][
            'Language'][0]
        attributes_map['pub_type'] = str(paper['MedlineCitation']['Article'][
            'PublicationTypeList'][0])
        attributes_map['pub_title'] = paper['MedlineCitation']['Article'][
            'ArticleTitle']
        attributes_map['pub_keywords'] = paper['MedlineCitation'][
            'KeywordList']
        articles_map[str(paper['MedlineCitation']['PMID'])] = attributes_map
    logging.debug('=========    ARTICLES MAP:     =========')
    logging.debug(articles_map)
    return articles_map


def create_article_year_by_topic_from_map(a_map):
    """ 
    Create year by topic dictionary, return ordered dictionary 

    Returns a dictionary containing information about the amount
    of articles written in particular year.
    """
    logging.debug(
        '============== IN create_article_year_by_topic_from_map: ================'
    )
    year_by_topic = {}
    for i in a_map:
        year = a_map[i]['pub_year']
        if not year in year_by_topic:
            year_by_topic[year] = 1
        else:
            year_by_topic[year] += 1
    logging.debug('=========    YEAR BY TOPIC MAP:     =========')
    logging.debug(sort_map_by_keys(year_by_topic))
    return sort_map_by_keys(year_by_topic)


def sort_map_by_keys(dictionary):
    """ Sort dictionary by keys, return dictionary """
    logging.debug('============== IN sort_map_by_keys: ================')
    d = OrderedDict(sorted(dictionary.items()))
    return d


def age_of_cited_work(art_map, art_cit_map, cit_map):
    """ 
    Create age of cited work dictionary, return dictionary 

    Calculates the ages of works cited by given article with
    respect to the year the article was written - returns a 
    dictionary mapping the year of article creation to the ages
    of articles referenced by it.
    """
    logging.debug('============== IN age_of_cited_work: ================')
    age = {}
    for article in art_cit_map:
        age_list = []
        logging.debug('========== ARTICLE: ============')
        logging.debug(article)
        for citation in art_cit_map[article]:
            logging.debug('========== CITATION: ============')
            logging.debug(citation)
            logging.debug(art_map[article]['pub_year'])
            if citation in cit_map:
                logging.debug(cit_map[citation]['pub_year'])
                age_list.append(
                    int(art_map[article]['pub_year']) - int(cit_map[citation][
                        'pub_year']))
                age[art_map[article]['pub_year']] = age_list
                logging.debug('========== AGE OF CITED WORK: ============')
                logging.debug(age)
    return sort_map_by_keys(age)


def keywords_map(articles_map):
    """ Create a dictionary of keywords (keyword, number of occurences), return dictionary """
    logging.debug('============== IN keywords_map: ================')
    keywords_map = {}
    for article in articles_map:
        if len(articles_map[article]['pub_keywords']) != 0:
            keyword_list = articles_map[article]['pub_keywords'][
                0]  # for each element in articles_map[article]['pub_keywords'] list get keyword
            logging.debug(keyword_list)
            for keyword in keyword_list:
                if not keyword.lower() in keywords_map:
                    keywords_map[keyword.encode('ascii', 'ignore').lower()] = 1
                else:
                    keywords_map[keyword.encode('ascii', 'ignore').lower(
                    )] += 1
    return keywords_map
