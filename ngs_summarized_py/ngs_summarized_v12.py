#!/usr/bin/python2.7
# This Python file uses the following encoding: utf-8

#------------------------------- Importing modules -----------------------------

import re
#import fungs as fungs
import os, sys

#----------------------------- Open and Write files ----------------------------

contigs = open('contigs.txt', 'r')
exons = open('exons.txt', 'r')
nucleotides = open('nuc.txt', 'r')
proteins = open('prot.txt', 'r')
annotation = open('annotation_out.txt', 'r')

data = open('summarized.txt','w')                      # output: summarized.txt
not_found =  open('not_found.txt', 'w')                 # output: not_found.txt

contigs_lines = contigs.readlines()
exons_lines = exons.readlines()
nucleotides_lines = nucleotides.readlines()
proteins_lines = proteins.readlines()
annotation_lines = annotation.readlines()

sequences_dic = {}

#-------------------------- Module: Reading annotation -------------------------

# It reads the annotation data and set if we have a group domain

def read_annotation(annotation_lines):
    data_annotation = []
    for i in range(0, len(annotation_lines)):
        annotation_dic = {}
        data_full = re.split('\d+\-\d+\[(\+|\-)\]', annotation_lines[i].strip())
        data = annotation_lines[i].strip().split()
        group_data = re.search('\^(Eukaryota|Bacteria|Archaea|Virus)', data_full[-1])
        if group_data:
            g = group_data.group()                                  #classified
        else:
            g = False                                             #unclassified
        annotation_dic['gen_id'] = data[0]                        #key <gen_id>
        annotation_dic['global'] = {'interval_prot': data[6],     #key <global>
         'recname': data[7],                                     #key <RecName>
         'full': data_full[-1],                                     #key <Full>
         'classified': g                                  #key <Classification>
         }
        data_annotation.append(annotation_dic)
#    for a in data_annotation:
#        print "%s\n"%(a)                                               #print_1
    return data_annotation

# {'gen_id': '13_g', 'global': {'interval_prot': '1-151[+]',
# 'recname': 'CQSS_VIBCB^CQSS_VIBCB^Q:4-151,H:3-150^60.81%ID^E:1e-57^RecName:',
# 'classified': '^Bacteria',
# 'full': '\tCQSS_VIBCB^CQSS_VIBCB^Q:4-151,H:3-150^60.81%ID^E:1e-57^RecName:
#  Full=CAI-1 autoinducer sensor kinase/phosphatase CqsS;^Bacteria; Proteoba.'}}

#--------------------------- Module: Reading sequences -------------------------

#fungs.read_sequeces()

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

# Classify group domains

eukaryota = open("eukaryota_set.txt", 'w')           #output: eukaryota_set.txt
eukaryota_fa = open("eukaryota_set.fa", 'w')           #output: eukaryota_set.txt
bacteria = open("bacteria_set.txt", 'w')              #output: bacteria_set.txt
bacteria_fa = open("bacteria_set.fa", 'w')
archaea = open("archaea_set.txt", 'w')
archaea_fa = open("archaea_set.fa", 'w')                 #output: archaea_set.txt
virus = open("virus_set.txt", 'w')                       #output: virus_set.txt
virus_fa = open("virus_set.fa", 'w')
unclassified = open("unclassified_set.txt", 'w')         # unclassified_set.txt
unclassified_fa = open("unclassified_set.fa", 'w')


def write_domains(str_intervals, ga, data):
    #ga = global_annotation[gen_id]['global']
    contig_id = data['id']
    len = data['len']

    l =  "\t".join([contigs_id, 'len='+len, str_intervals, ga['full'], str(ga['classified'])])
    g = str(ga['classified'])

    if(g == '^Eukaryota'):
        eukaryota.write(l + "\n")
        write_FASTA(eukaryota_fa, data)
    elif (g == '^Bacteria'):
        bacteria.write(l + "\n")
        write_FASTA(bacteria_fa, data)
    elif (g == '^Archaea'):
        archaea.write(l + "\n")
        write_FASTA(archaea_fa, data)
    elif (g == '^Virus'):
        virus.write(l + "\n")
        write_FASTA(virus_fa, data)
    else:
        unclassified.write(l + "\n")
        write_FASTA(unclassified_fa, data)
        return False
    return True

# Get the sequence by gen_id and set the classified data. (domain group)
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
        seq = contigs_lines[i+1].lstrip('\t')
        row = {
            'id': contig[0],
            'flag': contig[1],
            'multi': contig[2],
            'len': contig[3].split('=')[1],
            'seq': seq
            }
        data_contigs.append(row)                      #read contigs information
    return data_contigs

#[{'id': 'k101_1', 'len': '366'}

#------------------------------ Module: Read exons -----------------------------

def read_exons():
    exon_dic = {}
    for i in range(0, (len(exons_lines))):
        exons_values = exons_lines[i].strip().split()
        contigs_id = exons_values[0]
        gen_id = exons_values[9][1:-2]
        exon_len = {
            'intervalo_a': exons_values[3],
            'intervalo_b': exons_values[4],
            'direction': exons_values[6],
            'num': exons_values[7],
            'gen_id': gen_id
        }                                                             #relation
        if (exon_dic.has_key(contigs_id)):              #El arreglo no es vacío
            exon_dic[contigs_id].append(exon_len)
        else:
            exon_dic[contigs_id] = [exon_len]
    return exon_dic

# 1_g;55-321;+;0

#------------------------- Module: read_global_annotation ----------------------

#fungs.read_global_annotation()

# Converts the list to a hash
def read_global_annotation():
    annotations_list = read_annotation(annotation_lines)
    global_annotation = {}
    for i in range(len(annotations_list)):
        a = annotations_list[i]
        global_annotation[a['gen_id']] = {'global': a['global']}
    return global_annotation         #Add the rest of the annotation (for sets)

#-------------------------------------------------------------------------------

def coverage(len, sum):
    if(sum == 0):
        return 0
    else:
        r = (float(sum)/float(len))*100
        return round(r,2)


data_contigs = read_contigs()
exon_dic = read_exons()
global_annotation = read_global_annotation()

#------------------------------ Module: write_FASTA ----------------------------

def write_FASTA(f, data):
    contigs_id = data['id']
    contig_len = data['len']
    contig_multi = data['multi']
    contig_flag = data['flag']
    contig_seq = data['seq']
    f.write(" ".join([
        ">"+contigs_id,
        contig_flag,
        contig_multi,
        "len="+contig_len,'\n']))
    f.write(contig_seq)

# Iterando sobre los contigs (write_summarized)
for i in range(0, (len(data_contigs))):
    contigs_id = data_contigs[i]['id']
    contig_len = data_contigs[i]['len']
    #contig_multi = data_contigs[i]['multi']
    #contig_flag = data_contigs[i]['flag']
    #contig_seq = data_contigs[i]['seq']
    if (not exon_dic.has_key(contigs_id)):              #El arreglo no es vacío
        write_FASTA(not_found, data_contigs[i])
        # not_found.write(" ".join([
        #     ">"+contigs_id,
        #     contig_flag,
        #     contig_multi,
        #     "len="+contig_len,'\n']))
        # not_found.write(contig_seq)
        continue
    # Relación entre contigs y exons
    contig_exons = exon_dic[contigs_id]
    intervals = []
    str_intervals = ""
    sum_len_4_coverage = 0
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
        count = False
        if( (len(contig_exons) - 1) == j):
            sequences = get_sequences(current_gen_id, sequences_dic,\
            global_annotation)

            str_rest = "\t".join(["|".join(intervals), sequences])
            str_intervals = str_intervals + str_rest
            # Test write_domains
            count = write_domains(
                str_rest,
                global_annotation[current_gen_id]['global'],
                data_contigs[i]
            );

            #str_intervals = ""

        else:
            next_gen_id = contig_exons[j+1]['gen_id']
            if(not (current_gen_id == next_gen_id)):
                sequences = get_sequences(current_gen_id, sequences_dic,\
                 global_annotation)
                str_intervals = "\t".join(["|".join(intervals), sequences])+'\t'
                #str_intervals = str_intervals + "\t".join(["|".join(intervals), sequences]) + "\t"

                # Test write_domains
                count = write_domains(
                    str_intervals,
                    global_annotation[current_gen_id]['global'],
                    data_contigs[i]
                );

                intervals = []
        # coverage calculation
        if(count):
            sum_len_4_coverage = sum_len_4_coverage + (int(current['intervalo_b']) - int(current['intervalo_a']))
    cover = coverage(contig_len, sum_len_4_coverage)
    annotation_line_data = ""
    line = contigs_id+"\tlen="+contig_len+"\t"+str(cover)+"%\t"+str_intervals+\
    annotation_line_data+"\n"
    data.write(line)

not_found.close()
data.close()
unclassified.close()

eukaryota.close()
eukaryota_fa.close()
bacteria.close()
bacteria_fa.close()
archaea.close()
archaea_fa.close()
virus.close()
virus_fa.close()
unclassified.close()
unclassified_fa.close()
