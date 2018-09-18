#!/usr/bin/python2.7
# This Python file uses the following encoding: utf-8

#------------------------------- Importing modules -----------------------------

import re
#import fungs as fungs
import os, sys
import read_functions as rf

#----------------------------- Open and Write files ----------------------------

contigs = open('contigs.txt', 'r')
exons = open('exons.txt', 'r')
nucleotides = open('nuc.txt', 'r')
proteins = open('prot.txt', 'r')
annotation = open('annotation_out.txt', 'r')
data = open('summarized.txt','w')                      # output: summarized.txt
#unclassified = open('unclassified.txt', 'w')         # output: unclassified.txt
not_found =  open('not_found.txt', 'w')                 # output: not_found.txt

contigs_lines = contigs.readlines()
exons_lines = exons.readlines()
nucleotides_lines = nucleotides.readlines()
proteins_lines = proteins.readlines()
annotation_lines = annotation.readlines()

sequences_dic = {}

#{'gen_id': '13_g', 'global': {'interval_prot': '1-151[+]',
# 'recname': 'CQSS_VIBCB^CQSS_VIBCB^Q:4-151,H:3-150^60.81%ID^E:1e-57^RecName:',
# 'classified': '^Bacteria',
# 'full': '\tCQSS_VIBCB^CQSS_VIBCB^Q:4-151,H:3-150^60.81%ID^E:1e-57^RecName:
#  Full=CAI-1 autoinducer sensor kinase/phosphatase CqsS;^Bacteria; Proteoba.'}}

#-------------------------------------------------------------------------------

sequences_dic = rf.read_sequences(sequences_dic, 'nuc', nucleotides_lines)
sequences_dic = rf.read_sequences(sequences_dic, 'prot', proteins_lines)

#{'1_g: {'nuc': 'AGTC...', 'prot': 'AKLVWY...'}}

#----------------------------- Module: Get sequences ---------------------------

#fungs.get_sequences()

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

#-------------------------------------------------------------------------------

def coverage(len, sum):
    print len, sum
    if(sum == 0):
        return 0
    else:
        r = (float(sum)/float(len))*100
        return round(r,2)


data_contigs = rf.read_contigs(contigs_lines)
exon_dic = rf.read_exons(exons_lines)
global_annotation = rf.read_global_annotation(annotation_lines)

#Write FASTA
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

    if (not exon_dic.has_key(contigs_id)):              #El arreglo no es vacío
        write_FASTA(not_found, data_contigs[i])
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

        else:
            next_gen_id = contig_exons[j+1]['gen_id']
            if(not (current_gen_id == next_gen_id)):
                sequences = get_sequences(current_gen_id, sequences_dic,\
                 global_annotation)
                str_intervals = "\t".join(["|".join(intervals), sequences])+'\t'

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
