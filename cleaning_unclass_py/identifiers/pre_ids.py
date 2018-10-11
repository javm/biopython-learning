#!/usr/bin/env python

############################ Open and write files ##############################
# input
data = open('unclassified_contig_unclean.txt', 'r')
ids = open('unclassified_uniprot.txt', 'r')
# files for the output
not_found = open('unclass_pre_ids.txt', 'w')
found = open('uniprot_pre_ids.txt', 'w')
# arreglos y variables
lista_2 = {}
A1 = []
A2 = []

######################### Procesamiento de datos  ##############################

# data
for line1 in data:
    lines_d = line1.strip()
    columns_d = lines_d.split()
    A1.append(columns_d)                               # guarda en A1 (arreglo)

# k101_1570	.
# k101_1570	ROXA_ECOLI^ROXA_ECOLI^Q:109-481,H:3-370^48.53%ID^E:4e-122^.^

#ids
for line2 in ids:
    lines_i = line2.strip()
    columns_i =  lines_i.split()
    A2.append(columns_i)                               # guarda en A2 (arreglo)
#print (A2[0])

# ROXA_ECOLI^ROXA_ECOLI^Q:109-481,H:3-370^48.53%ID^E:4e-122^.^

####################### Comparaciones de archivos ##############################

# ids
for i in A2:                        #dic con ids
        lista_2[i[0]] = i[0]

# data

for line in A1:
    element = line[1]
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
