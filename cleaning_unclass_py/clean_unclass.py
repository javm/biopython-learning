#!/usr/bin/env python
import re
############################ Open and write files ##############################

uniprot_db = open("uniprot.txt", "r")
unclassified = open("unclassified_set.txt", "r")
uniprot_lines = open('uniprot_lines.txt','w')

# Putting uniprot in one line

uniprot_db_lines = uniprot_db.readlines()
#print uniprot_db_lines

l = uniprot_db_lines[0].strip()
one_line = [];

for i in range(1 , len(uniprot_db_lines)):
     line1 = uniprot_db_lines[i]
     line1 = line1.strip()

     #if(re.match('^[A-Z0-9]+\s[A-z{1}\(s|\t)+[0-9]+\:',line1)):
     if(re.search('[0-9]+\:',line1)):
         one_line.append(l)
         l = line1
     else:
         l = l + ';' + line1
for line in one_line:
    uniprot_lines.write(line+'\n')

uniprot_lines.write(line1+'\n')


# Reading uniprot lines
one_line_uniprot = open("uniprot_lines.txt", "r")
olu_lines = one_line_uniprot.readlines()
uniprot_dict = {}
for l in olu_lines:
    parts = l.strip().split()
    domain = ''
    dl = parts[1]
    if(dl == 'V'):
        domain = 'Virus'
    elif (dl == 'E'):
        domain = 'Eukaryota'
    elif (dl == 'A'):
        domain = 'Archaea'
    elif (dl == 'B'):
        domain = 'Bacteria'
    elif (dl == 'X'):
        domain = 'X'
    # {ABCDF: {domain: 'Virus', otra: ''} }
    # h['ABCDF']

    uniprot_dict[parts[0]] = {
        'domain': domain,
        'number': parts[2][:-1],
        'classification': parts[3]
    }

#print uniprot_dict


def uniprot_info(h, ui):
     s = ''
     d = h[ui]
     classification = d['classification]
 
     s = '^'+d['domain']+';\t'


######################### Procesamiento de datos  ##############################

unclassified_lines = unclassified.readlines()
unclassified_ids = {}

for line in unclassified_lines:
     parts = line.split()
     cols = parts[6].split('_')
     if (len(cols) > 1):
          uniprot_id = cols[1]
          unclassified_ids[parts[0]] = uniprot_id.split('^')[0]
          
          #subs(parts[6], 'False', uniprot_extract(uniprot_id))
          #out.write(parts)
     #out.write(line)
print unclassified_ids
