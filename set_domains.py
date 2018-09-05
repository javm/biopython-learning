import re

#Comes from coverage
#summarized_data = open("sum_coverage.txt", 'r')
summarized_data = open("summarized.txt", 'r')
eukaryota = open("eukaryota_set.txt", 'w')
bacteria = open("bacteria_set.txt", 'w')
archaea = open("archaea_set.txt", 'w')
virus = open("virus_set.txt", 'w')
unclassified = open("unclassified_set.txt", 'w')

for l in summarized_data:
    cols = l.split()
    cols_annotations = re.split('\t(\d+_g;\d+-d+;(\+|\-)\d\|*)+\t',l)
    print "cols_annotations[0]: "+ cols_annotations[0]
    print "\n cols_annotations[-1]: "+ cols_annotations[-1]
    group_data = re.search('\^(Eukaryota|Bacteria|Archaea|Virus)', l);

    if(group_data):
        g = group_data.group()
        print l
        print "contig: "+cols[0] +"group: "+g
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
