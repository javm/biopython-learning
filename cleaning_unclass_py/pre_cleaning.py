#!/usr/bin/env python
############################ Open and write files ##############################

# input
data = open('unclassified_unclean.txt', 'r')
ids = open('uniprot_ids.txt', 'r')

# output
not_found = open('unclassified_set.txt', 'w')
found = open('cleaning_set.txt', 'w')

# variables
lista_2 = {}
A1 = []
A2 = []

######################### Procesamiento de datos  ##############################

# data
for line1 in data:
    lines_d = line1.strip()
    columns_d = lines_d.split()
    A1.append(columns_d)                               # guarda en A1 (arreglo)

# k101_1570       len=558 1118_g;36-550;+;2       ATCTTTGCCATAACCGCCTTTCTGCCCAAG
# IFAITAFLPKPYAYWVFILGMVATQILYMVPSISVLRFERFVPRLGHMAERFALLTLIVLGEGFFKLVVTLSEKGIYK
# 1-171[+]        .       .       .       .       .       .       .       .
# .       .       False

#ids
for line2 in ids:
    lines_i = line2.strip()
    columns_i =  lines_i.split()
    A2.append(columns_i)                               # guarda en A2 (arreglo)

# k101_1090469

####################### Comparaciones de archivos ##############################

# ids
for i in A2:
    lista_2[i[0]] = i[0]

# data
for line in A1:
    element = line[0]
    value = element in lista_2                          # busca id en dic (ids)
    if value == True:
        word = ''                                               #word en un str
        for columns in line:      #columns = cada elemento del renglon de linea
            word = word+columns+'\t'
        found.write(word+'\n')
    if value == False:
        word = ''
        for columns in line:
            word = word+columns+'\t'
        not_found.write(word+'\n')

############################### Close files ####################################

data.close()
ids.close()
not_found.close()
found.close()
