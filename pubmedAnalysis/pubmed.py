from Bio import Entrez
Entrez.email = "yaroslava.romanets@uj.edu.pl"

def lookup(query):
    handle = Entrez.esearch(db='pubmed', 
    sort='relevance', 
    retmax='1600',
    retmode='xml', 
    term=query)
    results = Entrez.read(handle)
    return results

def get_details_by_id(id_list):
    ids = ','.join(id_list)
    handle = Entrez.efetch(db='pubmed',
    retmode='xml',
    id=ids)
    results = Entrez.read(handle)
    return results

if __name__ == '__main__':
    results = lookup('machine learning')
    id_list = results['IdList'] 
    papers = get_details_by_id(id_list)
    print type(papers)
    print papers['PubmedArticle'][0]
    for paper in papers['PubmedArticle']:
      print paper['MedlineCitation']['Article']['ArticleTitle']
#     print paper['MedlineCitation']['Article']['ArticleDate']    
      print paper['MedlineCitation']['DateCreated']['Year']
      print paper['PubmedData']['History'][0]['Year']#['PubMedPubDate']
    print '========== PubmedBookArticle ======='
    for book in papers['PubmedBookArticle']:
      print book
