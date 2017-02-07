import logging
import sys
from Bio import Entrez
from util import create_unverified_context
from settings import ENTREZ_EMAIL
Entrez.email = ENTREZ_EMAIL


def lookup(query, max_returned):
    """ Search for term in Entrez, return Entrez Parse object """
    create_unverified_context()
    handle = Entrez.esearch(
        db='pubmed',
        sort='relevance',
        retmax=max_returned,
        retmode='xml',
        term=query,
        usehistory="y")
    logging.debug('========== LOOKUP HANDLE: =======')
    results = Entrez.read(handle)
    logging.debug('=========== LOOKUP RESULTS: =======')
    logging.debug(results)
    return results


def get_details_by_id(id_list):
    """ Search for detailed information for each id from id_list, 
                        return Entrez Parse object """
    create_unverified_context()
    ids = ','.join(id_list)
    handle = Entrez.efetch(db='pubmed', retmode='xml', id=ids)
    try:
        results = Entrez.read(handle)
    except BaseException as e:
        logging.exception('Failed to create relevant dictionaries for given term: id_list=' + str(id_list) + \
        ' Exception: ' + str(e) + ' Please, try different search parameters - consider increasing max_returned values to 15.')
        sys.exit('Program executed with errors - see log for details.')
    return results


def get_citations_ids(id_list):
    """ Get ids of articles referenced by articles from id_list, return list """
    create_unverified_context()
    ids = ','.join(id_list)
    linked = []
    for i in range(0, len(ids), 200):
        handle = Entrez.elink(
            dbfrom="pubmed", id=ids[i:i + 200], linkname="pubmed_pubmed_refs")

        try:
            results = Entrez.read(handle)
        except BaseException as e:
            logging.exception('Failed to create relevant dictionaries for given term: id_list=' + str(id_list) + \
            ' Exception: ' + str(e) + ' Please, try different search parameters - consider increasing max_returned values to 15.')
            sys.exit('Program executed with errors - see log for details.')
        handle.close()
        logging.debug(
            '============== get_citations_ids RESULTS: ================')
        logging.debug(results)
        if len(results[0]["LinkSetDb"]) != 0:
            linked = linked + (
                [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]])
            logging.debug('====== REFERENCED ARTICLES: =========')
            logging.debug(linked)
    return linked
