from Bio import Entrez
Entrez.email = "yaroslava.romanets@uj.edu.pl"

def lookup(query,max_returned):
    handle = Entrez.esearch(db='pubmed', 
      sort='relevance', 
      retmax=max_returned,
      retmode='xml', 
      term=query,
      usehistory="y")
    print '========== HANDLE ======='
    print type(handle)
    results = Entrez.read(handle)
    print '=========== RESULTS ======='
    print type(results)
    print results
    return results
       # print handle.readline().strip()

def get_details_by_id(id_list, webenv):#, query_key):
    ids = ','.join(id_list)
    handle = Entrez.efetch(db='pubmed',
    retmode='xml',
    id=ids,
    webenv=webenv)#,query_key=query_key)
    results = Entrez.read(handle)
    return results

#def get_citations_ids(id_list):
def get_citations_ids(id_list, webenv):#, query_key):
  ids = ','.join(id_list)
  linked = []
  for i in range(0, len(ids), 200):
    handle = Entrez.elink(dbfrom="pubmed", id=ids[i:i+200], linkname="pubmed_pubmed", webenv=webenv)#, query_key=query_key)
    results = Entrez.read(handle)
    handle.close()

    # results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", id=ids[i:i+200]))
    # print '============== RESULTS ================'
    # print results
    # #handle.close()
    # pmc_ids = [link["Id"] for link in results[0]["LinkSetDb"]]#[0]["Link"]]
    # print '============== pmc_ids ================'
    # print pmc_ids
    # results2 = Entrez.read(Entrez.elink(dbfrom="pmc", db="pubmed", LinkName="pmc_pubmed", id=",".join(pmc_ids)))
    # print '============== RESULTS2 ================'
    # print results2

    print '====== LINKS BEFORE ========='
    print linked
    linked = linked + ([link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]])
    print '====== LINKS AFTER ========='
    print linked
    #print linked
  return linked
  #results["WebEnv"], results["QueryKey"] 
  #{'linked':linked, 'webenv':webenv,'query_key':query_key }

