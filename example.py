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

# {k_101: {len: 123, seq: 'fsdfsdf'}, k_104: {}}
# {k_101: [{}, {}, {}]}

exon_hash = {}
for i in range(0, (len(exons_lines) -1)):
    exon_values = exons_lines[i].strip().split()
    id = exon_values[0];
    if(exon_hash[id]):
        #
    else:
        exon_hash[id] = {
            intervalo_a: exon_values[3],
            intervalo_b: exon_values[4],
            direction: exon_values[6],
        }
