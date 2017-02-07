import logging
import sys
from Bio import Entrez
from settings import ENTREZ_EMAIL
Entrez.email = ENTREZ_EMAIL


def lookup(query, max_returned):
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


def get_details_by_id(id_list, webenv):
    ids = ','.join(id_list)
    handle = Entrez.efetch(db='pubmed', retmode='xml', id=ids, webenv=webenv)
    try:
        results = Entrez.read(handle)
    except BaseException as e:
        logging.exception('Failed to create relevant dictionaries for given term: id_list=' + str(id_list) + \
        ' Exception: ' + str(e) + ' Please, try different search parameters - consider increasing max_returned values to 15.')
        sys.exit('Program executed with errors - see log for details.')
    return results


def get_citations_ids(id_list, webenv='0'):
    ids = ','.join(id_list)
    linked = []
    for i in range(0, len(ids), 200):
        handle = Entrez.elink(
        dbfrom="pubmed",
        id=ids[i:i + 200],
        linkname="pubmed_pubmed_refs"
        , webenv=webenv)  #, query_key=query_key)
        try:
            results = Entrez.read(handle)
        except BaseException as e:
            logging.exception('Failed to create relevant dictionaries for given term: id_list=' + str(id_list) + \
            ' Exception: ' + str(e) + ' Please, try different search parameters - consider increasing max_returned values to 15.')
            sys.exit('Program executed with errors - see log for details.')
        handle.close()
        logging.debug('============== get_citations_ids RESULTS: ================')
        logging.debug(results)
        if len(results[0]["LinkSetDb"]) != 0:
            linked = linked + (
                [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]])
            logging.debug('====== REFERENCED ARTICLES: =========')
            logging.debug(linked)
    return linked

    # results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", id=ids[i:i+200]))
    # #handle.close()
    # pmc_ids = [link["Id"] for link in results[0]["LinkSetDb"]]#[0]["Link"]]
    # print '============== pmc_ids ================'
    # print pmc_ids
    # results2 = Entrez.read(Entrez.elink(dbfrom="pmc", db="pubmed", LinkName="pmc_pubmed", id=",".join(pmc_ids)))
    # print '============== RESULTS2 ================'
    # print results2
