from queries import lookup, get_citations_ids, get_details_by_id
from settings import ENTREZ_EMAIL
from dictionaries import create_maps, create_articles_map
from Bio import Entrez
import unittest
import ssl, os
Entrez.email = ENTREZ_EMAIL


class FunctionalTest(unittest.TestCase):

    def test_lookup(self):
        """ Testing queries.py lookup """
        term = 'artificial intelligence'
        result = lookup(term, 10)
        self.assertTrue(isinstance(result, dict))
        self.assertIn('IdList', result)
        self.assertTrue(len(result['IdList']) == 10)

    def test_get_details_by_id(self):
        """ Testing queries.py get_details_by_id """
        article_id = ['21876725']
        result = get_details_by_id(article_id)
        self.assertTrue(isinstance(result, dict))
        self.assertIn('PubmedArticle', result)
        self.assertEqual(result['PubmedArticle'][0]['MedlineCitation']['PMID'],
                         '21876725')
        self.assertEqual(result['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle'], 
            'Valuing ecological systems and services.')

    def test_get_details_by_id_exception(self):
        """ Testing queries.py get_details_by_id """
        article_id = []
        with self.assertRaises(BaseException) as context:
            get_details_by_id(article_id)
            self.assertTrue(
                'Please, try different search parameters' in context.exception)

    def test_get_citations_ids(self):
        """ Testing queries.py get_citations_ids """
        article_id = ['21876725']
        citations_list = ['19906066', '18202288']
        self.assertEqual(get_citations_ids(article_id), citations_list)

    def test_get_citations_ids_exception(self):
        """ Testing queries.py get_citations_ids """
        article_id = []
        with self.assertRaises(BaseException) as context:
            get_citations_ids(article_id)
            self.assertTrue(
                'Supplied id parameter is empty' in context.exception)

    def test_create_articles_map(self):
        """ Testing dictionaries.py create_articles_map """
        article_id = ['21876725']
        result = get_details_by_id(article_id)
        testdict = create_articles_map(result)
        self.assertTrue(isinstance(testdict, dict))
        self.assertEqual(testdict['21876725']['pub_title'],
                         'Valuing ecological systems and services.')

    def test_create_maps(self):
        """ Testing dictionaries.py create_maps """
        art_dic, art_cit_dic, cit_dic = create_maps('artificial intelligence',
                                                    30)
        self.assertTrue(isinstance(art_dic, dict))
        self.assertTrue(isinstance(art_cit_dic, dict))
        self.assertTrue(isinstance(cit_dic, dict))
        self.assertIn(list(art_cit_dic.keys())[0], art_dic)
        self.assertIn(list(art_cit_dic.values())[0][0], cit_dic)

if __name__ == '__main__':

    unittest.main()