#!/usr/bin/env python
import re
############################ Open and write files ##############################

uniprot_db = open('uniprot_test.txt', 'r')
uniprot_lines = open('uniprot_lines.txt','w')

A1 = []
A2 = []
element = ''
######################### Procesamiento de datos  ##############################

# data
for line1 in uniprot_db:
    line_parts = line1.strip().split()
    if(len(line_parts) > 2){
        
    }
    A1.append(lines_db)                              # guarda en A1 (arreglo)
# print(A1[0])

for i in range(0, len(A1)):
    element = element+' '+A1[i]
    #print (element)
    if i < len(A1)-1:
        # print ('entra')
        temp = A1[i+1]
        # print (temp)
        if re.match('\w=', temp):
            print ('entra')
            A2.append(element)
            element = ''

for i in A2:
    print (i)




    #print(A1[i])

# if "blah" not in somestring:
#     continue

# if element[0] == re.findall(r"^\w\=", element)
#    uniprot_lines.write()
# uniprot_lines.close()

# for element in A1:
#     m = re.match("(^\w)\W(d\w+)", element)
#     if m:
#         print(m.groups())
