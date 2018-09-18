import re
import os, sys

#-------------------------- Module: Reading annotation -------------------------
#ngs.read_annotation()
# It reads the annotation data and set if we have a group domain
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
            #unclassified.write(annotation_lines[i])
        annotation_dic['gen_id'] = data[0]                        #key <gen_id>
        annotation_dic['global'] = {'interval_prot': data[6],     #key <global>
         'recname': data[7],                                     #key <RecName>
         'full': data_full[-1],                                     #key <Full>
         'classified': g                                  #key <Classification>
         }
        data_annotation.append(annotation_dic)
    return data_annotation

#------------------------- Module: read_global_annotation ----------------------

#fungs.read_global_annotation()

# Converts the list to a hash
def read_global_annotation(annotation_lines):
    annotations_list = read_annotation(annotation_lines)
    global_annotation = {}
    for i in range(len(annotations_list)):
        a = annotations_list[i]
        global_annotation[a['gen_id']] = {'global': a['global']}
    return global_annotation         #Add the rest of the annotation (for sets)




#--------------------------- Module: Reading sequences -------------------------

#fungs.read_sequeces()

def read_sequences(sequences_out, sequences_name, lines):
    for i in range(0, len(lines), 2):
        gen_id = lines[i].lstrip('>')
        gen_id = gen_id.strip()
        if(not sequences_out.has_key(gen_id)):
            sequences_out[gen_id] = {}
        sequences_out[gen_id][sequences_name] = ''
        sequence = lines[i+1].strip()
        sequences_out[gen_id][sequences_name] = sequence
    return sequences_out


#----------------------------- Module: Read contigs ----------------------------

#fungs.read_contigs()

def read_contigs(contigs_lines):
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

#fungs.read_exons()

def read_exons(exons_lines):
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
        }
        if (exon_dic.has_key(contigs_id)):
            exon_dic[contigs_id].append(exon_len)
        else:
            exon_dic[contigs_id] = [exon_len]
    return exon_dic
