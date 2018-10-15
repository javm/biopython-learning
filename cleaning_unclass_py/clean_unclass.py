#!/usr/bin/env python
import re
############################ Open and write files ##############################

uniprot_db = open("uniprot.txt", "r")
uniprot_lines = open('uniprot_lines.txt','w')
#
# A1 = []
# A2 = []
# element = ''

######################### Procesamiento de datos  ##############################

# data
# for line1 in uniprot_db:
#     lines_db = line1.strip()
#     A1.append(lines_db)                              # guarda en A1 (arreglo)
# print(A1[0])

# data
uniprot_db_lines = uniprot_db.readlines()
#print uniprot_db_lines

l = uniprot_db_lines[0].strip()
one_line = [];

for i in range(1 , len(uniprot_db_lines)):
     line1 = uniprot_db_lines[i]
     line1 = line1.strip()

     line_parts = line1.split()

     if(len(line_parts) > 2):
         one_line.append(l)
         l = line1
     else:
         l = l + ' ' + line1
for line in one_line:
    uniprot_lines.write(line+'\n');


    #A1.append(lines_db)
#
#
# for i in range(0, len(A1)):
#     element = element+' '+A1[i]
#     #print (element)
#     if i < len(A1)-1:
#         # print ('entra')
#         temp = A1[i+1]
#         # print (temp)
#         if re.match('\w=', temp):
#             print ('entra')
#             A2.append(element)
#             element = ''
#
# for i in A2:
#     print (i)


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
