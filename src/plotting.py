import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
from util import avg, avg_list
from settings import SEARCH_TERM_1, SEARCH_TERM_2, MAX_RETURNED_TERM_1, MAX_RETURNED_TERM_2


def create_amount_of_article_citations_per_year_plot(
        term1_map, term2_map, term1_cit_map, term2_cit_map):
    plt.figure(1, figsize=(10, 5))
    plt.subplot(211).xaxis.set_major_formatter(
        tick.FormatStrFormatter('%0.0f'))
    plt.plot(
        list(term1_map.keys()),
        list(term1_map.values()),
        ':r',
        label=SEARCH_TERM_2)
    plt.plot(
        list(term2_map.keys()),
        list(term2_map.values()),
        '--c',
        label=SEARCH_TERM_1)
    plt.title("Amount of articles published about " + SEARCH_TERM_2 + " and " +
              SEARCH_TERM_1 + " per year")
    plt.xlabel("Year")
    plt.ylabel("Amount of articles")

    plt.legend(loc=2)

    plt.subplot(212).xaxis.set_major_formatter(
        tick.FormatStrFormatter('%0.0f'))
    plt.plot(
        list(term1_cit_map.keys()),
        list(term1_cit_map.values()),
        ':r',
        label=SEARCH_TERM_2)
    plt.plot(
        list(term2_cit_map.keys()),
        list(term2_cit_map.values()),
        '--c',
        label=SEARCH_TERM_1)
    plt.title("Amount of article citations per year for " + SEARCH_TERM_2 +
              " vs for " + SEARCH_TERM_1)
    plt.xlabel("Year")
    plt.ylabel("Amount of articles")

    plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)

    plt.legend(loc=2)
    plt.savefig('output/amount_of_article_citations_per_year.png')
    plt.show()


def create_avg_citations_age_per_year_plot(term1_cit_map, term2_cit_map,
                                           term1_age_of_cited_work,
                                           term2_age_of_cited_work):
    plt.figure(2, figsize=(10, 5))
    plt.subplot(211).xaxis.set_major_formatter(
        tick.FormatStrFormatter('%0.0f'))
    plt.plot(
        list(term1_cit_map.keys()),
        avg(list(term1_cit_map.values())),
        ':r',
        label=SEARCH_TERM_2)
    plt.plot(
        list(term2_cit_map.keys()),
        avg(list(term2_cit_map.values())),
        '--c',
        label=SEARCH_TERM_1)
    plt.title("Average amount of article citations per year for " +
              SEARCH_TERM_2 + " vs for " + SEARCH_TERM_1)
    plt.xlabel("Year")
    plt.ylabel("Average amount of articles")
    plt.legend(loc=2)

    plt.subplot(212).xaxis.set_major_formatter(
        tick.FormatStrFormatter('%0.0f'))
    print(list(term1_age_of_cited_work.keys()))
    print(avg_list(term1_age_of_cited_work.values()))
    ind_ai = np.arange(len(term1_age_of_cited_work))
    ind_ml = np.arange(len(term2_age_of_cited_work))
    plt.bar(ind_ai - 0.15,
            avg_list(term1_age_of_cited_work.values()),
            0.3,
            color='r',
            label=SEARCH_TERM_2)
    plt.bar(ind_ml + 0.15,
            avg_list(term2_age_of_cited_work.values()),
            0.3,
            color='c',
            label=SEARCH_TERM_1)
    plt.xticks(ind_ai, list(term1_age_of_cited_work.keys()))
    plt.xticks(ind_ml, list(term2_age_of_cited_work.keys()))

    plt.title("Average age cited works per year for " + SEARCH_TERM_2 +
              " vs for " + SEARCH_TERM_1)
    plt.xlabel("Year")
    plt.ylabel("Age of articles")

    plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)

    plt.legend(loc=2)
    plt.savefig('output/avg_citations_age_per_year.png')
    plt.show()


def growth_rate(term_map, term):
    growth_rate = np.exp(np.diff(np.log(list(term_map.values())))) - 1
    print(growth_rate)
    x = np.array(list(term_map.keys()))
    print(x)
    growth_rate = np.insert(growth_rate, 0, 0)
    print(growth_rate)
    print(len(growth_rate))
    print(len(x))

    plt.figure(3)
    plt.subplot(111).xaxis.set_major_formatter(
        tick.FormatStrFormatter('%0.0f'))
    plt.plot(x, growth_rate)
    plt.savefig('output/growth_' + term + '.png')
    plt.show()
