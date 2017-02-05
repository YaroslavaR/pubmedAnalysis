import xlsxwriter


def create_keywords_excel(doc, keywords_map):
    workbook = xlsxwriter.Workbook(doc)
    worksheet = workbook.add_worksheet()
    row = 0
    col = 0
    for key in keywords_map.keys():
        worksheet.write(row, col, key)
        worksheet.write(row, col + 1, keywords_map[key])
        row += 1
    workbook.close()

def avg(slist):  
    print('========== DICT ============')
    print(slist)
    for i in range(0,len(slist)):  
        slist[i] = float(slist[i])/len(slist)
    return slist

def avg_list(slist):  
    print('========== DICT ============')
    print(slist)
    rlist = []
    for elem in slist:
      rlist.append(sum(elem)/float(len(elem)))
    slist = rlist
    print(slist)
    return slist