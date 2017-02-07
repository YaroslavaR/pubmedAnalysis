import logging
import xlsxwriter
from settings import OUTPUT_PATH


def create_keywords_excel(doc, keywords_map):
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
    logging.debug('========== DICT ============')
    logging.debug(slist)
    rlist = []
    for i in range(0, len(slist)):
        rlist.append(round(((float(slist[i]) / sum(slist)) * 100), 2))
    return rlist


def avg_list(slist):
    logging.debug('========== DICT ============')
    logging.debug(slist)
    rlist = []
    for elem in slist:
        rlist.append(sum(elem) / float(len(elem)))
    slist = rlist
    logging.debug(slist)
    return slist
