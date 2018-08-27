import re

summarized_data = open("summarized.txt", 'r')
eukaryota = open("eukaryota.txt", 'w')
bacteria = open("bacteria.txt", 'w')
archaea = open("archaea.txt", 'w')
virus = open("virus.txt", 'w')
unclassified = open("unclassified_group.txt", 'w')

for l in summarized_data:
    cols = l.split()
    group_data = re.search('\^(Eukaryota|Bacteria|Archaea|Virus)', \
        cols[-1]);
    if(group_data):
        g = group_data.group()
    else:
        g = ''
    if(g == '^Eukaryota'):
        eukaryota.write(l + "\n")
    elif (g == '^Bacteria'):
        bacteria.write(l + "\n")
    elif (g == '^Archaea'):
        archaea.write(l + "\n")
    elif (g == '^Virus'):
        virus.write(l + "\n")
    else:
        unclassified.write(l + "\n")
eukaryota.close()
bacteria.close()
archaea.close()
virus.close()
unclassified.close()
