from collections import OrderedDict
import matplotlib.pyplot as plt
from Bio import Entrez
Entrez.email = "yaroslava.romanets@uj.edu.pl"

def lookup(query,max_returned):
    handle = Entrez.esearch(db='pubmed', 
    sort='relevance', 
    retmax=max_returned,
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

def create_article_year_by_topic_map(query_result):
    year_by_topic = {}
    for paper in query_result['PubmedArticle']:
      year = paper['PubmedData']['History'][0]['Year']#['PubMedPubDate']
      if not year in year_by_topic:
        year_by_topic[year] = 1
      else:
        year_by_topic[year] += 1
    return year_by_topic

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

if __name__ == '__main__':
    results = lookup('machine learning',1600)
    id_list = results['IdList'] 
    papers = get_details_by_id(id_list)
    
    papers_year_dict = create_article_year_by_topic_map(papers)
    books_year_dict = create_book_year_by_topic_map(papers)
    
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

    ai_results = lookup('artificial intelligence',70000)
    ai_id_list = ai_results['IdList']
    ai_papers = get_details_by_id(ai_id_list)

    ai_papers_year_dict = create_article_year_by_topic_map(ai_papers)
    ai_books_year_dict = create_book_year_by_topic_map(ai_papers)

    ai_pdo = sort_map_by_keys(ai_papers_year_dict)
    ai_books_pdo = sort_map_by_keys(ai_books_year_dict)


print '========== ML PubmedArticle ======='
for year in pdo:
  print year
  print pdo[year]

print '========== ML PubmedBookArticle ======='

for year in books_year_dict:
  print year
  print books_year_dict[year]


print '========== ML PubmedArticle ======='
for year in ai_pdo:
  print year
  print ai_pdo[year]

print '========== ML PubmedBookArticle ======='
for year in ai_books_pdo:
  print year
  print ai_books_pdo[year]
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


#plt.plot(ai_books_pdo.keys(), ai_books_pdo.values(), ':r', label='AI')
#plt.plot(books_pdo.keys(),books_pdo.values(), '--c', label='ML')
#plt.title("Amount of books published about ML and AI per year")
#plt.xlabel("Year")
#plt.ylabel("Amount of books");

#plt.legend()
#plt.show()
