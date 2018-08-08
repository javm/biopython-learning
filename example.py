#!/usr/bin/env python3
import pprint

# PART 1: OPEN FILES

contigs = open('contigs.txt', 'r')
exons = open('exons.txt', 'r')
nucleotides = open('nuc.txt', 'r')
proteins = open('prot.txt', 'r')
anotation = open('anotation.txt', 'r')

data_contigs = []
#File to write
data = open('data.txt','w')

# PART 1.1: OPEN BY LINES

contigs_lines = contigs.readlines()
exons_lines = exons.readlines()
nucleotides_lines = nucleotides.readlines()
proteins_lines = proteins.readlines()
anotation = anotation.readlines()

# PART 1.2: REMOVE 1ST CHARACTER OF THE FILE CONTIGS.TXT

for i in range(0, len(contigs_lines), 2):
    contig = contigs_lines[i].lstrip('>')
    contig = contig.strip().split()
    row = {'id': contig[0], 'len': contig[3].split('=')[1]}
    sequence = contigs_lines[i+1].strip()
    row['seq'] = sequence
    print(row)
    data_contigs.append(row)
# {k_101: {len: 123, seq: 'fsdfsdf'}, k_104: {}}
# {k_101: [{}, {}, {}]}

exon_hash = {}
for i in range(0, (len(exons_lines))):
    exon_values = exons_lines[i].strip().split()
    id = exon_values[0];
    gen_id = exon_values[9][1:-2]

    data_exon = {
            'intervalo_a': exon_values[3],
            'intervalo_b': exon_values[4],
            'direction': exon_values[6],
            'num': exon_values[7],
            'gen_id': gen_id
        }

    if (exon_hash.has_key(id)):
        # El arreglo no es vacio
        exon_hash[id].append(data_exon)
    else:
        exon_hash[id] = [data_exon]

for i in range(0, (len(data_contigs))):
    contig_id = data_contigs[i]['id']
    # Este es un arreglo
    contig_exons = exon_hash[contig_id]
    # TODO: crear un funcion para formato
    intervals = []
    for j in range(0, (len(contig_exons))):
        current = contig_exons[j]
        interval_data = [ current['gen_id'],
            ("-".join( [current['intervalo_a'],current['intervalo_b']])),
            current['direction'],
            current['num']
        ]
        interval = ";".join(interval_data)
        intervals.append(interval)

    s_intervals = "|".join(intervals)
    sequence = data_contigs[i]['seq']
    line = contig_id+"\tlen="+data_contigs[i]['len']+"\t"+s_intervals+" \t"+sequence+"\n"
    print(line)
    data.write(line)
data.close()
