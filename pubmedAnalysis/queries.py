from Bio import Entrez
Entrez.email = "yaroslava.romanets@uj.edu.pl"

def lookup(query,max_returned):
    handle = Entrez.esearch(db='pubmed', 
      sort='relevance', 
      retmax=max_returned,
      retmode='xml', 
      term=query)
    print '========== HANDLE ======='
    print type(handle)
    results = Entrez.read(handle)
    print '=========== RESULTS ======='
    print type(results)
    print results
    return results
       # print handle.readline().strip()

def get_details_by_id(id_list):
    #ids = '' 
    #if (not isinstance(id_list, list)):
    #print 'id_list'
    #print type(id_list)
    #if (isinstance(id_list,list)):
     # print id_list
    ids = ','.join(id_list)
    #print 'ids'
    #print type(ids)
    #else:
     # ids = id_list
    handle = Entrez.efetch(db='pubmed',
    retmode='xml',
    id=ids)
    results = Entrez.read(handle)
    return results

def get_citations_ids(id_list):
  #ids = ''
  #if (not isinstance(id_list, list)):
  ids = ','.join(id_list)
  #else:
   # ids = id_list
  linked = []
  for i in range(0, len(ids), 200):
    handle = Entrez.elink(dbfrom="pubmed", id=ids[i:i+200], linkname="pubmed_pubmed")
    results = Entrez.read(handle)
    handle.close()
    print '====== LINKS BEFORE ========='
    print linked
    linked = linked + ([link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]])
    print '====== LINKS AFTER ========='
    print linked
    #print linked
  return linked

