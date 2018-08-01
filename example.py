#/usr/bin/python

# PART 1: OPEN FILES

contigs = open('contigs.txt', 'r')
exons = open('exons.txt', 'r')
nucleotides = open('nuc.txt', 'r')
proteins = open('prot.txt', 'r')
anotation = open('anotation.txt', 'r')
data = open('data.txt','w')
data_array = []

# PART 1.1: OPEN BY LINES

contigs_lines = contigs.readlines()
exons_lines = exons.readlines()
nucleotides_lines = nucleotides.readlines()
proteins_lines = proteins.readlines()
anotation = anotation.readlines()

# PART 1.2: REMOVE 1ST CHARACTER OF THE FILE CONTIGS.TXT

for i in range(0, len(contigs_lines) - 1, 2):
    contig = contigs_lines[i].lstrip('>')
    contig = contig.strip().split()
    row = {'id': contig[0], 'len': contig[3].split('=')[1]}
    sequence = contigs_lines[i+1].strip()
    row['seq'] = sequence
    print(row)

#lista_1 = {}
#lista_2 = {}
#repetidos = []
#no_repetidos = []

#i=0
#while i<len(A1)-1:
#        lista_1[A1[1]
#for element in A1:
#        b_names[element]=A1[i]
#        i=i+1
