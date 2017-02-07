import logging
import xlsxwriter
import ssl
from settings import OUTPUT_PATH


def create_keywords_excel(doc, keywords_map):
    """ Write keywords_map to excel doc """
    workbook = xlsxwriter.Workbook(OUTPUT_PATH + '/' + doc)
    worksheet = workbook.add_worksheet()
    row = 0
    col = 0
    for key in keywords_map.keys():
        worksheet.write(row, col, key)
        worksheet.write(row, col + 1, keywords_map[key])
        row += 1
    workbook.close()


def percentage(slist):
    """ Calculate percentage of each item in list, return list """
    logging.debug('========== DICT ============')
    logging.debug(slist)
    rlist = []
    for i in range(0, len(slist)):
        rlist.append(round(float(slist[i]) / sum(slist) * 100, 2))
    return rlist


def avg_list(slist):
    """ Calculate average of each item in list, return list """
    logging.debug('========== DICT ============')
    logging.debug(slist)
    rlist = []
    for elem in slist:
        rlist.append(sum(elem) / float(len(elem)))
    slist = rlist
    logging.debug(slist)
    return slist


def create_unverified_context():
    try:
        _create_unverified_https_context = ssl._create_unverified_context
    except AttributeError:
        # Legacy Python that doesn't verify HTTPS certificates by default
        logging.info(
            'Using Python in v. < 3.0 - using unverified context by default')
    else:
        # Handle target environment that doesn't support HTTPS verification
        ssl._create_default_https_context = _create_unverified_https_context
        logging.debug('Using Python in v. > 3.0 - setting unverified context')
