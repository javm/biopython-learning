# coding: utf-8
#!/usr/bin/python2.7
#------------------------------- Importing modules -----------------------------
import os, sys
import read_functions as rf
#------------------------------- Open data files -------------------------------

contigs = open('contigs.txt', 'r')
exons = open('exons.txt', 'r')
nucleotides = open('nuc.txt', 'r')
proteins = open('prot.txt', 'r')
annotation = open('annotation_out.txt', 'r')
sequences_dic = {}

#------------------------------------ Write sets -------------------------------

data = open('summarized_set.txt','w')                  # output: summarized.txt
not_found = open('not_found_set.fa', 'w')               # output: not_found.txt

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

#--------------------------------- Debugger No.1 -------------------------------
print "Step_01:readlines_start"
#------------------------- Returns everything into lines -----------------------

contigs_lines = contigs.readlines()
exons_lines = exons.readlines()
nucleotides_lines = nucleotides.readlines()
proteins_lines = proteins.readlines()
annotation_lines = annotation.readlines()

#--------------------------------- Debugger No.1 -------------------------------
print "Step_01:readlines_end"
#-------------------------- Library: read_functions.py -------------------------

#--------------------------------- Debugger No.2 -------------------------------
print "Step_02:read_sequences_start"
#-------------------------------------------------------------------------------

sequences_dic = rf.read_sequences(sequences_dic, 'nuc', nucleotides_lines)
sequences_dic = rf.read_sequences(sequences_dic, 'prot', proteins_lines)

#--------------------------------- Debugger No.2 -------------------------------
print "Step_02:read_sequences_end"
#----------------------------- Function: Write fastas --------------------------

def write_fasta(f, data):
    contig_id = data['id']
    contig_len = data['len']
    contig_multi = data['multi']
    contig_flag = data['flag']
    contig_seq = data['seq']
    f.write(" ".join([">" + contig_id, contig_flag, contig_multi, "len=" + contig_len, '\n']))
    f.write(contig_seq)

#---------------------------- Function: Write domains --------------------------

def write_domains(str_intervals, ga, data):
    contig_id = data['id']
    contig_len = data['len']
    g = str(ga['classified'])
    l =  "\t".join([contigs_id, 'len=' + contig_len, str_intervals, ga['full'], g])
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

#---------------------------- Function: Get sequences --------------------------
#     Get the sequence by gen_id and set the classified data (domain group)
#-------------------------------------------------------------------------------

def get_sequences(gen_id, sequences_dic, get_annotation):
    sequences = ""
    if(sequences_dic.has_key(gen_id)):
        sequences = "\t".join([sequences_dic[current_gen_id]['nuc'], sequences_dic[current_gen_id]['prot']])
    if(get_annotation.has_key(gen_id)):
        get_a = get_annotation[gen_id]['global']
        sequences = "\t".join([sequences,
            get_a['interval_prot'],
            get_a['recname'],
            get_a['full'],
            str(get_a['classified'])                         #annotation_output
            ])
    return sequences

#----------------------------- Coverage calculation ----------------------------

def coverage(leng, summ):
    if(summ == 0):
        return 0
    else:
        r = (float(summ)/float(leng))*100
        return round(r,2)

def Percent_gen(gen_id):
    per_str = global_annotation[gen_id]['global']['percentage']
    return per_str

#--------------------------------- Debugger No.3 -------------------------------
print "Step_03:read_contigs_start"
#-------------------------------------------------------------------------------
data_contigs = rf.read_contigs(contigs_lines)
#--------------------------------- Debugger No.3 -------------------------------
print "Step_03:read_contigs_end"
#-------------------------------------------------------------------------------
#--------------------------------- Debugger No.4 -------------------------------
print "Step_04:read_exons_start"
#-------------------------------------------------------------------------------
exon_dic = rf.read_exons(exons_lines)
#--------------------------------- Debugger No.4 -------------------------------
print "Step_04:read_exons_end"
#-------------------------------------------------------------------------------
#--------------------------------- Debugger No.5 -------------------------------
print "Step_05:read_global_annotation_start"
#-------------------------------------------------------------------------------
global_annotation = rf.read_global_annotation(annotation_lines)
#--------------------------------- Debugger No.5 -------------------------------
print "Step_05:read_global_annotation_end"

#------------------------------- > ID evaluation -------------------------------

for i in range(0, (len(data_contigs))):
    contigs_id = data_contigs[i]['id']#data_contigs = rf.read_contigs(contigs.readlines())
    contig_len = data_contigs[i]['len']
    if (not exon_dic.has_key(contigs_id)):#exon_dic = rf.read_exons(exons.readlines())
        write_fasta(not_found, data_contigs[i])
        continue
    contig_exons = exon_dic[contigs_id]
    intervals = []
    str_intervals = ""
    str_intervals_2 = ""
    sum_len_4_coverage = 0

    Percents = []
    for j in range(0, (len(contig_exons))):
        current = contig_exons[j]
        current_gen_id = current['gen_id']
        Percents.append(float(Percent_gen(current_gen_id)))

    max_per = max(Percents)
    max_per = str(max_per)

    for j in range(0, (len(contig_exons))):
        current = contig_exons[j]
        current_gen_id = current['gen_id']
        if Percent_gen(current_gen_id) == max_per:
            max_gen = current_gen_id

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

        if current_gen_id == max_gen:
            sequences = get_sequences(current_gen_id, sequences_dic, global_annotation)
            str_rest = "\t".join(["|".join(intervals), sequences])
            str_intervals = str_intervals + str_rest

        count = False
        #if((len(contig_exons) - 1) == j):
        if(j == (len(contig_exons) - 1)):
            sequences = get_sequences(current_gen_id, sequences_dic, global_annotation)
            #sequences_dic tiene gen y secuencias de nucleótidos y aminoácidos (nuc y prot)
            #global_annotation tiene gen_id, global: interval_prot, recname, full, classified, y porcentaje
            #get_sequences regresa la unión:'nuc','prot','global','interval_prot','recname','full','classified'
            str_rest = "\t".join(["|".join(intervals), sequences])
            str_intervals = str_intervals + str_rest
            #l =  "\t".join([contigs_id, 'len=' + contig_len, str_intervals, ga['full'], g])
            #count = write_domains(str_rest, global_annotation[current_gen_id]['global'], data_contigs[i])#;
            if current_gen_id == max_gen:
                str_rest_2 = "\t".join(["|".join(intervals), sequences])
                str_intervals_2 = str_intervals_2 + str_rest_2
                count = write_domains(str_rest_2, global_annotation[current_gen_id]['global'], data_contigs[i])#;
        else:
            next_gen_id = contig_exons[j+1]['gen_id']
            if(not (current_gen_id == next_gen_id)):
                sequences = get_sequences(current_gen_id, sequences_dic, global_annotation)
                str_intervals = "\t".join(["|".join(intervals), sequences])+'\t'
                if current_gen_id == max_gen:
                    str_intervals_2 = "\t".join(["|".join(intervals), sequences])+'\t'
                    count = write_domains(str_intervals_2,global_annotation[current_gen_id]['global'],data_contigs[i])#;
                #count = write_domains(str_intervals,global_annotation[current_gen_id]['global'],data_contigs[i])#;
                intervals = []
        if(count):
            sum_len_4_coverage = sum_len_4_coverage + (int(current['intervalo_b']) - int(current['intervalo_a']))
    cover = coverage(contig_len, sum_len_4_coverage)
    #annotation_line_data = ""
    line = contigs_id + "\tlen=" + contig_len + "\t" + str(cover) + "%\t" + str_intervals + "\n" #+annotation_line_data+"\n"
    data.write(line)

#--------------------------------- Close files ---------------------------------

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

#----------------------------------- Cluster -----------------------------------
# contigs = open('final.contigs.PE.fa', 'r')
# exons = open('genemark_reduced_cds.gtf', 'r')
# nucleotides = open('nuc_seq.mod.fna', 'r')
# proteins = open('prot_seq.mod.faa', 'r')
# annotation = open('annotation_out.xls', 'r')
#
# data = open('summarized_set.txt','w')                # output: summarized.txt
# not_found =  open('not_found_set.fa', 'w')            # output: not_found.txt
#----------------------------------- Pruebas -----------------------------------
# contigs = open('contigs.txt', 'r')
# exons = open('exons.txt', 'r')
# nucleotides = open('nuc.txt', 'r')
# proteins = open('prot.txt', 'r')
# annotation = open('annotation_out.txt', 'r')
#
# data = open('summarized_set.txt','w')                # output: summarized.txt
# not_found =  open('not_found_set.fa', 'w')            # output: not_found.txt
#-------------------------------------------------------------------------------
