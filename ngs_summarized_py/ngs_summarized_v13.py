#!/usr/bin/python2.7
# This Python file uses the following encoding: utf-8

#------------------------------- Importing modules -----------------------------

from collections import OrderedDict
import re
import os, sys

#----------------------------- Opens and write files ---------------------------

contigs = open('contigs.txt', 'r')
exons = open('exons.txt', 'r')
nucleotides = open('nuc.txt', 'r')
proteins = open('prot.txt', 'r')
annotation = open('annotation_out.txt', 'r')

data = open('summarized_set.txt','w')                  # output: summarized.txt
not_found =  open('not_found_set.txt', 'w')             # output: not_found.txt

#------------------------------------ Open sets --------------------------------

eukaryota = open("eukaryota_set.txt", 'w')           #output: eukaryota_set.txt
eukaryota_fa = open("eukaryota_set.fa", 'w')         #output: eukaryota_set.txt
bacteria = open("bacteria_set.txt", 'w')              #output: bacteria_set.txt
bacteria_fa = open("bacteria_set.fa", 'w')
archaea = open("archaea_set.txt", 'w')
archaea_fa = open("archaea_set.fa", 'w')               #output: archaea_set.txt
virus = open("virus_set.txt", 'w')                       #output: virus_set.txt
virus_fa = open("virus_set.fa", 'w')
unclassified = open("unclassified_set.txt", 'w')         # unclassified_set.txt
unclassified_fa = open("unclassified_set.fa", 'w')

#------------------------- Returns everything into lines -----------------------
print "Step_01:readlines_start"

contigs_lines = contigs.readlines()
exons_lines = exons.readlines()
nucleotides_lines = nucleotides.readlines()
proteins_lines = proteins.readlines()
annotation_lines = annotation.readlines()

print "Step_01:readlines_end"
#-------------------------- Module: Reading annotation -------------------------

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
    return data_annotation

#--------------------------- Module: Reading sequences -------------------------

sequences_dic = {}

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

print "Step_02:read_sequences_start"
read_sequences(sequences_dic, 'nuc', nucleotides_lines)
read_sequences(sequences_dic, 'prot', proteins_lines)
print "Step_02:read_sequences_end"

#----------------------------- Module: Get sequences ---------------------------

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
            #get_a['recname'], este no va porque full ya lo tiene
            get_a['full'],
            str(get_a['classified'])                         #annotation_output
            ])
    return sequences

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

#------------------------------ Module: Read exons -----------------------------

def read_exons():
    #exon_dic = {}
    exon_dic = OrderedDict()
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

#------------------------- Module: Read global annotations ---------------------

def read_global_annotation():
    annotations_list = read_annotation(annotation_lines)
    global_annotation = {}
    for i in range(len(annotations_list)):
        a = annotations_list[i]
        global_annotation[a['gen_id']] = {'global': a['global']}
    return global_annotation         #Add the rest of the annotation (for sets)

#----------------------------- Coverage calculation ----------------------------

def coverage(len, sum):
    if(sum == 0):
        return 0
    else:
        r = (float(sum)/float(len))*100
        return round(r,2)

print "Step_03:read_contigs_start"
data_contigs = read_contigs()
print "Step_03:read_contigs_end"
print "Step_04:read_exons_start"
exon_dic = read_exons()
print "Step_04:read_exons_end"
print "Step_05:read_global_annotation_start"
global_annotation = read_global_annotation()
print "Step_05:read_global_annotation_end"

#------------------------------ Module: Write fastas ---------------------------

def write_fasta(f, data):
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

#------------------------------- Module: Write sets ----------------------------

def write_domains(str_intervals, ga, data):
    contig_id = data['id']
    len = data['len']
    l =  "\t".join([contigs_id, 'len='+len, str_intervals])
    g = str(ga['classified'])

    if(g == '^Eukaryota'):
        eukaryota.write(l + "\n")
        write_fasta(eukaryota_fa, data)
    elif (g == '^Bacteria'):
        bacteria.write(l + "\n")
        write_fasta(bacteria_fa, data)
    elif (g == '^Archaea'):
        archaea.write(l + "\n")
        write_fasta(archaea_fa, data)
    elif (g == '^Virus'):
        virus.write(l + "\n")
        write_fasta(virus_fa, data)
    else:
        unclassified.write(l + "\n")
        write_fasta(unclassified_fa, data)
        return False
    return True

#------------------------------- > ID evaluation -------------------------------

def maxper(num):
    datos = annotation_lines[num].strip().split()
    count = 0
    A = []
    a = 0
    for element in datos:
        Z = re.search('\%'+'ID'+'\^',element)
        if bool(Z) == True:
            a = count
            A.append(a)
        count += 1
    if len(A)>0:
        C = []
        for b in A:
            B = re.split('\:',datos[b])
            C.append(B)
        D = []
        for element in C:
            D.append(element[2])
        E = []
        for element in D:
            x = re.split('\^',element)
            E.append(x[1])
        F=[]
        for element in E:
            x = re.split('\%',element)
            F.append(float(x[0]))
        value = max(F)
    else:
        value=0.
    return value

#------------------------------- > ID evaluation -------------------------------

print "Step_06:new_start"
New=[]
count = 0
for element in read_exons():
    New.append([])
    New[count].append(element)
    for gene in read_exons()[element]:
        if gene['gen_id'] not in New[count]:
            New[count].append(gene['gen_id'])
    count += 1
print "Step_06:new_end"

print "Step_07:for_start"
counter = 0

for i in range(0, (len(data_contigs))):                                        #
    contigs_id = data_contigs[i]['id']
    contig_len = data_contigs[i]['len']
    if (not exon_dic.has_key(contigs_id)):
        write_fasta(not_found, data_contigs[i])
        continue
    contig_exons = exon_dic[contigs_id]
    intervals = []
    str_intervals = ""
    sum_len_4_coverage = 0

    count2 = 1
    Parr = []
    while count2 < len(New[counter]):
        v = int(re.split('\_',New[counter][count2])[0])-1
        Parr.append(maxper(v))
        count2 += 1
    grt = max(Parr)

    count3 = 1
    if grt != 0.:
        while count3 < len(New[counter]):
            v = int(re.split('\_',New[counter][count3])[0])-1
            datos = annotation_lines[v].strip().split()
            for element in datos:
                Z1 = re.search(str(grt)+'\%'+'ID',element)
                Z2 = re.search(str(int(grt))+'\%'+'ID',element)
                if bool(Z1) == True or bool(Z2) == True:
                    ident = annotation_lines[v].strip().split()[0]
            count3 +=1
    else:
        ident='no'
    counter += 1
    for j in range(0, (len(contig_exons))):
        current = contig_exons[j]
        current_gen_id = current['gen_id']
        # if current_gen_id == ident:
        if current_gen_id == ident or ident == "no":
            interval_data = [ current_gen_id,
                ("-".join( [current['intervalo_a'],
                current['intervalo_b']])),
                current['direction'],
                current['num']]
            interval = ";".join(interval_data)
            intervals.append(interval)
            sequences = ""
            count = False
            if ((len(contig_exons) - 1) == j):
                sequences = get_sequences(current_gen_id, sequences_dic, global_annotation)
                str_rest = "\t".join(["|".join(intervals), sequences])
                str_intervals = str_intervals + str_rest
                count = write_domains(
                    str_rest,
                    global_annotation[current_gen_id]['global'],
                    data_contigs[i]
                    );
            else:
                next_gen_id = contig_exons[j+1]['gen_id']
                if(not (current_gen_id == next_gen_id)):
                    sequences = get_sequences(current_gen_id, sequences_dic, global_annotation)
                    str_intervals = "\t".join(["|".join(intervals), sequences])+'\t'
                    count = write_domains(
                        str_intervals,
                        global_annotation[current_gen_id]['global'],
                        data_contigs[i]
                    );
                    intervals = []
            if(count):
                sum_len_4_coverage = sum_len_4_coverage + (int(current['intervalo_b']) - int(current['intervalo_a']))
        # elif ident == 'no':
        #     if( (len(contig_exons) - 1) == j):
        #         sequences = get_sequences(current_gen_id, sequences_dic, global_annotation)
        #         str_rest = "\t".join(["|".join(intervals), sequences])
        #         str_intervals = str_intervals + str_rest
        #     else:
        #         next_gen_id = contig_exons[j+1]['gen_id']
        #         if(not (current_gen_id == next_gen_id)):
        #             sequences = get_sequences(current_gen_id, sequences_dic, global_annotation)
        #             str_intervals = "\t".join(["|".join(intervals), sequences])+'\t'
    cover = coverage(contig_len, sum_len_4_coverage)
    annotation_line_data = ""
    line = contigs_id+"\tlen="+contig_len+"\t"+str(cover)+"%\t"+str_intervals+annotation_line_data+"\n"
    data.write(line)
print "Step_07:for_end"

#--------------------------------- Close files ---------------------------------

not_found.close()
data.close()

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

#----------------------------------- Cluster -----------------------------------

# contigs = open('final.contigs.PE.fa', 'r')
# exons = open('genemark_reduced_cds.gtf', 'r')
# nucleotides = open('nuc_seq.mod.fna', 'r')
# proteins = open('prot_seq.mod.faa', 'r')
# annotation = open('annotation_out.xls', 'r')
#
# data = open('summarized_set.txt','w')                  # output: summarized.txt
# not_found =  open('not_found_set.txt', 'w')             # output: not_found.txt

#----------------------------------- Pruebas -----------------------------------

# contigs = open('contigs.txt', 'r')
# exons = open('exons.txt', 'r')
# nucleotides = open('nuc.txt', 'r')
# proteins = open('prot.txt', 'r')
# annotation = open('annotation_out.txt', 'r')
#
# data = open('summarized_set.txt','w')                  # output: summarized.txt
# not_found =  open('not_found_set.txt', 'w')             # output: not_found.txt
