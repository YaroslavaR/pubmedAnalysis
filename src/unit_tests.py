from util import percentage, avg_list
from queries import get_citations_ids
from dictionaries import age_of_cited_work
from collections import OrderedDict
import ssl
import unittest


class MyTest(unittest.TestCase):
    def test_avg_list(self):
        testlist = [[2, 3, 4, 2, 3], [2, 3, 4, 2, 3]]
        self.assertEqual(avg_list(testlist), [2.8, 2.8])

    def test_percentage(self):
        testlist = [2, 3, 4, 2, 3]
        self.assertEqual(
            percentage(testlist), [14.29, 21.43, 28.57, 14.29, 21.43])

    # def test_get_citations_ids(self):
    # 	article_id = '21876725'
    # 	citations_list = ['19906066', '18202288']
    # 	self.assertEqual(get_citations_ids(article_id), citations_list)
    def test_age_of_cited_work(self):
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


if __name__ == '__main__':
    try:
        _create_unverified_https_context = ssl._create_unverified_context
    except AttributeError:
        # Legacy Python that doesn't verify HTTPS certificates by default
        pass
    else:
        # Handle target environment that doesn't support HTTPS verification
        ssl._create_default_https_context = _create_unverified_https_context

    unittest.main()
