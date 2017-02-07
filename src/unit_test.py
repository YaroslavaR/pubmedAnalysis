from util import percentage, avg_list, create_keywords_excel
from queries import lookup, get_citations_ids, get_details_by_id
from dictionaries import age_of_cited_work, sort_map_by_keys, create_article_year_by_topic_from_map, keywords_map, \
    get_citations_ids_map, create_articles_map, create_maps
from settings import ENTREZ_EMAIL, OUTPUT_PATH
from collections import OrderedDict
from Bio.Entrez.Parser import StringElement, ListElement
from Bio import Entrez
import ssl, os
import unittest
Entrez.email = ENTREZ_EMAIL


class MyTest(unittest.TestCase):
    def test_lookup(self):
        """ Testing queries.py lookup """
        term = 'artificial intelligence'
        result = lookup(term, 10)
        self.assertTrue(isinstance(result, dict))
        self.assertIn('IdList', result)
        self.assertTrue(len(result['IdList']) > 5)

    def test_get_details_by_id(self):
        """ Testing queries.py get_details_by_id """
        article_id = ['21876725']
        result = get_details_by_id(article_id)
        self.assertTrue(isinstance(result, dict))
        self.assertIn('PubmedArticle', result)
        self.assertEqual(result['PubmedArticle'][0]['MedlineCitation']['PMID'],
                         '21876725')

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

    def test_get_citations_ids_map(self):
        """ Testing dictionaries.py get_citations_ids_map """
        article_id = ['21876725']
        article_citations_dict = {'21876725': ['19906066', '18202288']}
        self.assertEqual(
            get_citations_ids_map(article_id), article_citations_dict)

    def test_create_articles_map(self):
        """ Testing dictionaries.py create_articles_map """
        article_id = ['21876725']
        result = get_details_by_id(article_id)
        testdict = create_articles_map(result)
        self.assertTrue(isinstance(testdict, dict))
        self.assertEqual(testdict['21876725']['pub_title'],
                         'Valuing ecological systems and services.')

    def test_age_of_cited_work(self):
        """ Testing dictionaries.py age_of_cited_work """
        art_map = {
            '25153475': {
                'pub_year': '2014',
                'pub_lang': 'eng',
                'pub_type': 'Journal Article',
                'pub_title':
                'Artificial intelligence in public health prevention of legionelosis in drinking water systems.',
                'pub_keywords': []
            }
        }
        art_cit_map = {'25153475': ['20122385', '18276473', '7818640']}
        cit_map = {
            '20122385': {
                'pub_year': '2010',
                'pub_lang': 'eng',
                'pub_type': 'Journal Article',
                'pub_title':
                'Preliminary report: outbreak of Legionnaires disease in the cities of Ulm and Neu-Ulm in Germany, December 2009 - January 2010.',
                'pub_keywords': []
            },
            '18276473': {
                'pub_year': '1992',
                'pub_lang': 'eng',
                'pub_type': 'Journal Article',
                'pub_title':
                'Neural networks designed on approximate reasoning architecture and their applications.',
                'pub_keywords': []
            },
            '7818640': {
                'pub_year': '1994',
                'pub_lang': 'eng',
                'pub_type': 'Journal Article',
                'pub_title':
                'A massive outbreak in Milwaukee of cryptosporidium infection transmitted through the public water supply.',
                'pub_keywords': []
            }
        }
        self.assertEqual(
            age_of_cited_work(art_map, art_cit_map, cit_map),
            OrderedDict({
                '2014': [4, 22, 20]
            }))

    def test_sort_map_by_keys(self):
        """ Testing dictionaries.py sort_map_by_keys """
        testdict = {'2014': 1, '1994': 10, '1995': 8, '2017': 2, '1986': 1}
        self.assertEqual(
            sort_map_by_keys(testdict),
            OrderedDict([('1986', 1), ('1994', 10), ('1995', 8), ('2014', 1),
                         ('2017', 2)]))

    def test_create_article_year_by_topic_from_map(self):
        """ Testing dictionaries.py create_article_year_by_topic_from_map """
        testdict = {
            '23110900': {
                'pub_year': '2012',
                'pub_lang': 'eng',
                'pub_type': 'Comment',
                'pub_title': 'In silico plant biology comes of age.',
                'pub_keywords': []
            },
            '23110897': {
                'pub_year': '2012',
                'pub_lang': 'eng',
                'pub_type': 'Journal Article',
                'pub_title':
                'Multiscale systems analysis of root growth and development: modeling beyond the network and cellular scales.',
                'pub_keywords': []
            },
            '23110896': {
                'pub_year': '2012',
                'pub_lang': 'eng',
                'pub_type': 'Journal Article',
                'pub_title':
                'Modeling regulatory networks to understand plant development: small is beautiful.',
                'pub_keywords': []
            },
            '22353237': {
                'pub_year': '2011',
                'pub_lang': 'eng',
                'pub_type': 'Journal Article',
                'pub_title':
                'An ovary transcriptome for all maturational stages of the striped bass (Morone saxatilis), a highly advanced perciform fish.',
                'pub_keywords': []
            }
        }
        self.assertEqual(
            create_article_year_by_topic_from_map(testdict),
            OrderedDict([('2011', 1), ('2012', 3)]))

    def test_keywords_map(self):
        """ Testing dictionaries.py keywords_map """
        testdict = {
            '26265491': {
                'pub_year': '2015',
                'pub_lang': 'eng',
                'pub_type': 'Journal Article',
                'pub_title':
                'Thirty years of artificial intelligence in medicine (AIME) conferences: A review of research themes.',
                'pub_keywords': [
                    ListElement([
                        StringElement('Artificial Intelligence in Medicine'),
                        StringElement('History of science'),
                        StringElement('Literature review')
                    ])
                ]
            },
            '23743214': {
                'pub_year': '2012',
                'pub_lang': 'eng',
                'pub_type': 'Journal Article',
                'pub_title':
                'Patient behavior and the benefits of artificial intelligence: the perils of "dangerous" literacy and illusory patient empowerment.',
                'pub_keywords': [
                    ListElement([
                        StringElement('Health literacy'),
                        StringElement('literature review'),
                        StringElement('Patient education'),
                        StringElement('Patient empowerment')
                    ])
                ]
            },
            '26021669': {
                'pub_year': '2015',
                'pub_lang': 'eng',
                'pub_type': 'Journal Article',
                'pub_title':
                'Applying artificial intelligence technology to support decision-making in nursing: A case study in Taiwan.',
                'pub_keywords':
                [ListElement([StringElement('nursing informatics')])]
            }
        }
        self.assertEqual(
            keywords_map(testdict), {
                'artificial intelligence in medicine': 1,
                'history of science': 1,
                'literature review': 2,
                'health literacy': 1,
                'patient education': 1,
                'patient empowerment': 1,
                'nursing informatics': 1
            })

    def test_create_maps(self):
        """ Testing dictionaries.py create_maps """
        art_dic, art_cit_dic, cit_dic = create_maps('artificial intelligence',
                                                    30)
        self.assertTrue(isinstance(art_dic, dict))
        self.assertTrue(isinstance(art_cit_dic, dict))
        self.assertTrue(isinstance(cit_dic, dict))
        self.assertIn(list(art_cit_dic.keys())[0], art_dic)
        self.assertIn(list(art_cit_dic.values())[0][0], cit_dic)

    def test_avg_list(self):
        """ Testing util.py avg_list """
        testlist = [[2, 3, 4, 2, 3], [2, 3, 4, 2, 3]]
        self.assertEqual(avg_list(testlist), [2.8, 2.8])

    def test_percentage(self):
        """ Testing util.py percentage """
        testlist = [2, 3, 4, 2, 3]
        self.assertEqual(
            percentage(testlist), [14.29, 21.43, 28.57, 14.29, 21.43])

    @staticmethod
    def test_create_keywords_excel(self):
        """ Testing util.py create_keywords_excel """
        testdict = {
            'artificial intelligence in medicine': 1,
            'history of science': 1,
            'literature review': 2
        }
        create_keywords_excel('output', testdict)
        assert os.path.exists(OUTPUT_PATH + '/output') == 1


if __name__ == '__main__':

    unittest.main()
