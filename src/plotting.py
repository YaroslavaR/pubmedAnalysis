import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
import logging
from util import percentage, avg_list
from settings import SEARCH_TERM_1, SEARCH_TERM_2, MAX_RETURNED_TERM_1, MAX_RETURNED_TERM_2


def create_amount_of_article_citations_per_year_plot(
        term1_map, term2_map, term1_cit_map, term2_cit_map):
    """ Create plot showing amount of articles and their citations per year """
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


def create_avg_citations_age_per_year_plot(term1_age_of_cited_work,
                                           term2_age_of_cited_work):
    """ Create bar charts showing average age of cited works per year """
    plt.figure(2, figsize=(10, 5))
    plt.subplot(211).xaxis.set_major_formatter(
        tick.FormatStrFormatter('%0.0f'))
    logging.debug(list(term1_age_of_cited_work.keys()))
    logging.debug(avg_list(term1_age_of_cited_work.values()))

    ind_ml = np.arange(len(term1_age_of_cited_work))
    plt.bar(ind_ml,
            avg_list(term1_age_of_cited_work.values()),
            0.3,
            color='r',
            label=SEARCH_TERM_1)
    plt.xticks(ind_ml, list(term1_age_of_cited_work.keys()))

    plt.title("Average age of cited works per year for " + SEARCH_TERM_1)
    plt.xlabel("Year")
    plt.ylabel("Age of articles")

    plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)

    plt.legend(loc=2)

    plt.subplot(212).xaxis.set_major_formatter(
        tick.FormatStrFormatter('%0.0f'))

    logging.debug(list(term2_age_of_cited_work.keys()))
    logging.debug(avg_list(term2_age_of_cited_work.values()))

    ind_ml = np.arange(len(term2_age_of_cited_work))
    plt.bar(ind_ml,
            avg_list(term2_age_of_cited_work.values()),
            0.3,
            color='c',
            label=SEARCH_TERM_2)
    plt.xticks(ind_ml, list(term2_age_of_cited_work.keys()))

    plt.title("Average age of cited works per year for " + SEARCH_TERM_2)
    plt.xlabel("Year")
    plt.ylabel("Age of articles")

    plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)

    plt.legend(loc=2)
    plt.savefig('output/avg_citations_age_per_year.png')
    plt.show()


def percentage_per_year_pie(term_map, term):
    """ Create a pie chart showing percentage of articles per year """
    labels = list(term_map.keys())
    sizes = percentage(list(term_map.values()))

    fig1, ax1 = plt.subplots(figsize=(10, 7.5))
    ax1.pie(sizes,
            labels=labels,
            autopct='%1.1f%%',
            shadow=True,
            startangle=90,
            labeldistance=1.2)
    ax1.axis(
        'equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.title("Percentage of articles per year for " + term, y=1.08)

    plt.savefig('output/percentage_per_year_pie_' + term + '.png')
    plt.show()


def growth_rate(term_map, term):
    """ Create a plot showing yearly growth trend """
    growth_rate = np.exp(np.diff(np.log(list(term_map.values())))) - 1
    logging.debug(growth_rate)
    x = np.array(list(term_map.keys()))
    logging.debug(x)
    growth_rate = np.insert(growth_rate, 0, 0)
    logging.debug(growth_rate)
    logging.debug(len(growth_rate))
    logging.debug(len(x))

    plt.figure(3)
    plt.subplot(111).xaxis.set_major_formatter(
        tick.FormatStrFormatter('%0.0f'))
    plt.plot(list(term_map.keys()), growth_rate)
    plt.title("Growth trend for " + term + " per year")
    plt.savefig('output/growth_' + term + '.png')
    plt.show()
