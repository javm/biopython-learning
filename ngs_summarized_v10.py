#!/usr/bin/env python3
#------------------------------- Importing modules -----------------------------

import re
import defings as ngs

#----------------------------- Open and Write files ----------------------------

contigs = open('contigs.txt', 'r')
exons = open('exons.txt', 'r')
nucleotides = open('nuc.txt', 'r')
proteins = open('prot.txt', 'r')
annotation = open('annotation_out.txt', 'r')
data = open('summarized.txt','w')
unclassified = open('unclassified.txt', 'w')

contigs_lines = contigs.readlines()
exons_lines = exons.readlines()
nucleotides_lines = nucleotides.readlines()
proteins_lines = proteins.readlines()
annotation_lines = annotation.readlines()

sequences_dic = {}

#-------------------------- Module: Reading annotation -------------------------

#ngs.read_annotation()

def read_annotation(annotation_lines):
    data_annotation = []
    for i in range(0, len(annotation_lines)):
        annotation_dic = {}
        data_full = re.split('\d+\-\d+\[(\+|\-)\]', annotation_lines[i].strip())
        data = annotation_lines[i].strip().split()
        group_data = re.search('\^(Eukaryota|Bacteria|Archaea|Virus)',\
        data_full[-1])
        if group_data:
            g = group_data.group()                                  #classified
        else:
            g = False                                             #unclassified
            unclassified.write(annotation_lines[i])
        annotation_dic['gen_id'] = data[0]                        #key <gen_id>
        annotation_dic['global'] = {'interval_prot': data[6],     #key <global>
         'recname': data[7],                                     #key <RecName>
         'full': data_full[-1],                                     #key <Full>
         'classified': g                                  #key <Classification>
         }
        data_annotation.append(annotation_dic)
    for a in data_annotation:
        print "%s\n"%(a)                                               #print_1
    return data_annotation

#{'gen_id': '13_g', 'global': {'interval_prot': '1-151[+]',
# 'recname': 'CQSS_VIBCB^CQSS_VIBCB^Q:4-151,H:3-150^60.81%ID^E:1e-57^RecName:',
# 'classified': '^Bacteria',
# 'full': '\tCQSS_VIBCB^CQSS_VIBCB^Q:4-151,H:3-150^60.81%ID^E:1e-57^RecName:
#  Full=CAI-1 autoinducer sensor kinase/phosphatase CqsS;^Bacteria; Proteoba.'}}

#--------------------------- Module: Reading sequences -------------------------

#ngs.read_sequeces()

def read_sequences(sequences_out, sequences_name, lines):
    for i in range(0, len(lines), 2):
        gen_id = lines[i].lstrip('>')
        gen_id = gen_id.strip()
        if(not sequences_dic.has_key(gen_id)):
            sequences_out[gen_id] = {}
        sequences_out[gen_id][sequences_name] = ''
        sequence = lines[i+1].strip()
        sequences_out[gen_id][sequences_name] = sequence

#-------------------------------------------------------------------------------

read_sequences(sequences_dic, 'nuc', nucleotides_lines)
read_sequences(sequences_dic, 'prot', proteins_lines)

#{'1_g: {'nuc': 'AGTC...', 'prot': 'AKLVWY...'}}

#----------------------------- Module: Get sequences ---------------------------

#ngs.get_sequences()

def get_sequences(gen_id, sequences_dic, get_annotation):
    sequences = ""
    if(sequences_dic.has_key(gen_id)):
        sequences =  "\t".join([
            sequences_dic[current_gen_id]['nuc'],
            sequences_dic[current_gen_id]['prot']             #sequences_output
        ])
    if(get_annotation.has_key(gen_id)):
        get_a = get_annotation[gen_id]['global']
        sequences = "\t".join([sequences,
            get_a['interval_prot'],
            get_a['recname'],
            get_a['full'],
            str(get_a['classified'])                         #annotation_output
            ])
    return sequences

# OUTPUT: NUCLEOTIDES_ATCG AA_YWVAVL ANNOTATION_XXX

#----------------------------- Module: Read contigs ----------------------------

def read_contigs():
    data_contigs = []
    for i in range(0, len(contigs_lines), 2):
        contig = contigs_lines[i].lstrip('>')
        contig = contig.strip().split()
        row = {'id': contig[0], 'len': contig[3].split('=')[1]}
        data_contigs.append(row)                      #read contigs information
    return data_contigs

#------------------------------ Module: Read exons -----------------------------

def read_exons():
    exon_hash = {}
    for i in range(0, (len(exons_lines))):
        exons_values = exons_lines[i].strip().split()
        contigs_id = exons_values[0]
        gen_id = exons_values[9][1:-2]                                #relation
        exon_dic = {
            'intervalo_a': exons_values[3],
            'intervalo_b': exons_values[4],
            'direction': exons_values[6],
            'num': exons_values[7],
            'gen_id': gen_id
        }
        if (exon_hash.has_key(contigs_id)):
            # El arreglo no es vacio
            exon_hash[contigs_id].append(exon_dic)
        else:
            exon_hash[contigs_id] = [exon_dic]
    return exon_hash

#------------------------------ Module: Read exons -----------------------------

#annotation global
def read_annotation_global():
    annotations_list = read_annotation(annotation_lines)
    annotation_global = {}
    #print annotations_list[1]
    for i in range(len(annotations_list)):
        a = annotations_list[i]
        annotation_global[a['gen_id']] = {'global': a['global']}
    return annotation_global

data_contigs = read_contigs()
exon_hash = read_exons()
annotation_global = read_annotation_global()

for i in range(0, (len(data_contigs))):
    contigs_id = data_contigs[i]['id']
    # Este es un arreglo
    contig_exons = exon_hash[contigs_id]
    # Crear un funcion para formato
    intervals = []
    s_intervals = ""

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

        sequences = ""
        if( (len(contig_exons) - 1) == j):
            sequences = get_sequences(current_gen_id, sequences_dic, annotation_global)
            s_intervals = s_intervals + "\t".join(["|".join(intervals), sequences])
        else:
            next_gen_id = contig_exons[j+1]['gen_id']
            if(not (current_gen_id == next_gen_id)):
                sequences = get_sequences(current_gen_id, sequences_dic, annotation_global)
                s_intervals = "\t".join(["|".join(intervals), sequences]) + "\t"
                intervals = []

    annotation_line_data = ""
    line = contigs_id+"\tlen="+data_contigs[i]['len']+"\t"+s_intervals+annotation_line_data+"\n"
    #print(line)
    data.write(line)
data.close()
unclassified.close()
