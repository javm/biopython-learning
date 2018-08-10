#!/usr/bin/#!/usr/bin/env python3

# -----------------------------------------------------------
# Marisol
# UNAM - IBt
# This program takes as input NGS output form a WGS, in order
# to make a .txt summary that includes: assembly, annotation,
# nucleotides,proteins and genes predictions.
# -----------------------------------------------------------

# PART 1: OPEN FILES

# Test files
contigs = open('contigs.txt', 'r')
exons = open('exons.txt', 'r')
nucleotides = open('nuc.txt', 'r')
proteins = open('prot.txt', 'r')
annotation = open('annotation.txt', 'r')
# Files to write
data = open('summarized.txt','w')
data_contigs = []
# PART 1.1: OPEN BY LINES
contigs_lines = contigs.readlines()
exons_lines = exons.readlines()
nucleotides_lines = nucleotides.readlines()
proteins_lines = proteins.readlines()
annotation = annotation.readlines()
gen_sequences = {}
# Reading gen sequences
def read_gen_sequence(gen_sequences_out, sequence_name, lines):
    for i in range(0, len(lines), 2):

        gen_id = lines[i].lstrip('>')
        gen_id = gen_id.strip()
        if(not gen_sequences.has_key(gen_id)):
            gen_sequences_out[gen_id] = {}
        gen_sequences_out[gen_id][sequence_name] = ''
        sequence = lines[i+1].strip()
        gen_sequences_out[gen_id][sequence_name] = sequence
    #print(gen_sequences_out)
    #return gen_sequences_out

#gent_sequences =
read_gen_sequence(gen_sequences, 'prot', proteins_lines)
#gen_sequences =
read_gen_sequence(gen_sequences, 'nuc', nucleotides_lines)

# {'1_g: {'nuc': 'AGTC...', 'prot': 'GTCAA...'}}
print(gen_sequences)

def add_sequences(intervals, gen_sequences, gen_id):
    if(gen_sequences.has_key(gen_id)):
        intervals.append("\t".join([
            gen_sequences[current_gen_id]['nuc'],
            gen_sequences[current_gen_id]['prot']
        ]))


# PART 2:

for i in range(0, len(contigs_lines), 2):
    contig = contigs_lines[i].lstrip('>')
    contig = contig.strip().split()
    row = {'id': contig[0], 'len': contig[3].split('=')[1]}
    sequence = contigs_lines[i+1].strip()
    row['seq'] = sequence
    #print(row)
    data_contigs.append(row)

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
    # Crear un funcion para formato
    intervals = []

    for j in range(0, (len(contig_exons))):
        current = contig_exons[j]
        current_gen_id = current['gen_id']

        interval_data = [ current_gen_id,
            ("-".join( [current['intervalo_a'],
            current['intervalo_b']])),
            current['direction'],
            current['num']
        ]
        interval = ";".join(interval_data)
        intervals.append(interval)

        if( (len(contig_exons) - 1) == j):
            add_sequences(intervals, gen_sequences, current_gen_id)
        else:
            next_gen_id = contig_exons[j+1]['gen_id']
            if(not (current_gen_id == next_gen_id)):
                add_sequences(intervals, gen_sequences, current_gen_id)
    s_intervals = "|".join(intervals)

    sequence = data_contigs[i]['seq']
    line = contig_id+"\tlen="+data_contigs[i]['len']+"\t"+s_intervals+"\n"
    print(line)
    data.write(line)
data.close()
