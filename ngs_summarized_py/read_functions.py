#------------------------------- Importing modules -----------------------------
import re
#---------------------------- Module 1: Read contigs --------------------------#
# It reads the contigs data and make a dictionary with the header of each seq  #
#------------------------------------------------------------------------------#
# >k101_5711 flag=1 multi=4.0000 len=603
# CAAGGAAGGAACGGGAAGAGCATAGGTACAACAGCCCAATCAGCACACTCACTGACTGTGCCAATACCGTTGCAACAG
# >k101_14001 flag=1 multi=4.0000 len=711
# GCTGCCAGGCGGCGATGATAGCGGCCAGCGGGTCGGTAACACGCGCGGTGGGGCAGACGGATTGGATCGCGTCCGCCA
#------------------------------------------------------------------------------#
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
        data_contigs.append(row)
    return data_contigs
#----------------------------- Module 2: Read exons ---------------------------#
#It reads the prediction/exons data and make an association between contigs ids#
#                                 and genes ids                                #
#------------------------------------------------------------------------------#
# k101_113765     GeneMark.hmm    CDS     135     452     .       +       0
# gene_id "1_g";  transcript_id "1_t";
# k101_110238     GeneMark.hmm    CDS     29      581     .       +       2
# gene_id "2_g";  transcript_id "2_t";
#------------------------------------------------------------------------------#
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
#------------------------- Module 3: Reading annotation -----------------------#
#        It reads the annotation data and set if we have a group domain        #
#------------------------------------------------------------------------------#
# 1_g     comp1_c0_seq1   .       .       .       cds.comp1_c0_seq1|m.1
# 1-106[+]        .       .       .       .       .       .       .       .
# 7_g     comp7_c0_seq1   .       .       .       cds.comp7_c0_seq1|m.7
# 1-320[+]        MFD_RICFE^MFD_RICFE^Q:18-308,H:28-313^31.62%ID^E:2e-40^RecName:
# Full=Transcription-repair-coupling factor {ECO:0000255|HAMAP-Rule:MF_00969};
# ^Bacteria; Proteobacteria; Alphaproteobacteria; Rickettsiales; Rickettsiaceae;
# Rickettsieae; Rickettsia; spotted fever   group   .       .       .
# COG1197^transcriptioN-repair coupling factor GO:0005737^cellular_component^
# cytoplasm`GO:0005524^molecular_function^ATP binding`GO:0003684^molecular_
# function^damaged DNA binding`GO:0004386^molecular_function^helicase activity`
# GO:0006355^biological_process^regulation of transcription, DNA-templated`
# GO:0000716^biological_process^transcription-coupled nucleotide-excision repair,
# DNA damage recognition     .       .       .
#------------------------------------------------------------------------------#
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
#-------------------------- Module 4: read_id_percentage ----------------------#
#                        Extract the percentage of identity                    #
#------------------------------------------------------------------------------#
# UNC89_CAEEL^UNC89_CAEEL^Q:4-185,H:5185-5367^32.43%ID^E:9e-20^RecName:
# DNLJ_NOVAD^DNLJ_NOVAD^Q:325-499,H:8-186^76.54%ID^E:3e-85^RecName:
#------------------------------------------------------------------------------#
# def read_id_percentage(annotation_lines):
#     annotations_list = read_annotation(annotation_lines)
#     annotation_ids = {}
#     for i in range(len(annotation_list)):
#         get_recname = annotations_list['global']
#         recname = get_recname['recname']
#     return annotation_ids
#------------------------ Module 5: read_global_annotation --------------------#
#                           Converts the list to a hash                        #
#------------------------------------------------------------------------------#
def read_global_annotation(annotation_lines):
    annotations_list = read_annotation(annotation_lines)
    global_annotation = {}
    for i in range(len(annotations_list)):
        a = annotations_list[i]
        global_annotation[a['gen_id']] = {'global': a['global']}
    return global_annotation
#-------------------------- Module 6: Reading sequences -----------------------#
#    It makes an association between gen ids and the corresponding sequence    #
#------------------------------------------------------------------------------#
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
