#!/usr/bin/python2.7
# ---------------------------------------------------------------
# Marisol
# UNAM - IBt
# This program takes as input a text file with a percentage (%)
# value asociated to one of its columns, then gets only the
# floating value of this percentage and plots the frequency of
# events between 0%-100% separated by 10 bins.
#
#               Dependencies: < hisfunc.py >
# ---------------------------------------------------------------

#################################################################
#                             Option 3                          #
#################################################################

import hisfunc as hfunc

hfunc.plot_histogram(hfunc.plot_coverage_histogram('eukaryota'), 10, "Histograma de cobertura sobre anotaciones con ID", "histogram_cover_euk.png")
hfunc.plot_histogram(hfunc.plot_coverage_histogram('bacteria'), 10, "Histograma de cobertura sobre anotaciones con ID", "histogram_cover_bac.png")
hfunc.plot_histogram(hfunc.plot_coverage_histogram('archaea'), 10, "Histograma de cobertura sobre anotaciones con ID", "histogram_cover_arc.png")
hfunc.plot_histogram(hfunc.plot_coverage_histogram('virus'), 10, "Histograma de cobertura sobre anotaciones con ID", "histogram_cover_vir.png")

#################################################################
# Option 1: from a .json list with the values only
#################################################################

#import json
#import hisfunc as hfunc

#with open("list.json","r") as myfile:
#    list = json.load(myfile)

#hfunc.plot_histogram(list, 10, "Histograma", "histogram.png")

#################################################################
# Option 2: from a .txt file with percentage values in the column
# with position 2
#################################################################

#import re
#import hisfunc as hfunc

# list_percentage = []
# with open("/Users/solouli/Desktop/Git_hub/Ngs_scripts/ngs_histogram_cover_py/example_set.txt", "r") as set:
#     for i in set:
#         if re.search('^k',i):
#              list_percentage.append(round(float(i.split()[2].split('%')[0])))
#              #print(list_percentage)
#
# hfunc.plot_histogram(list_percentage, 10, "Histograma de cobertura sobre anotaciones", "histogram_cover.png")

#################################################################
# Option 2.1: from a .txt file with percentage values in the column
# with position 2 on TLAHUICA.
#################################################################

#import re
#import hisfunc as hfunc

#Eukaryota

#list_percentage_euk = []
#with open("/scratch01/mnavarro/SUMMARIZED/eukaryota_set.txt", "r") as set:

#    for i in set:
#        if re.search('^k',i):
#             list_percentage_euk.append(round(float(i.split()[2].split('%')[0])))
#hfunc.plot_histogram(list_percentage_euk, 10, "Histograma de cobertura sobre anotaciones en: EUKARYOTA", "histogram_cover_euk.png")

#Bacteria

#list_percentage_bac = []
#with open("/scratch01/mnavarro/SUMMARIZED/bacteria_set.txt", "r") as set:
#    for i in set:
#        if re.search('^k',i):
#             list_percentage_bac.append(round(float(i.split()[2].split('%')[0])))
#hfunc.plot_histogram(list_percentage_bac, 10, "Histograma de cobertura sobre anotaciones en: BACTERIA", "histogram_cover_bac.png")

# Archaea

#list_percentage_arc = []
#with open("/scratch01/mnavarro/SUMMARIZED/archaea_set.txt", "r") as set:
#    for i in set:
#        if re.search('^k',i):
#             list_percentage_arc.append(round(float(i.split()[2].split('%')[0])))
#hfunc.plot_histogram(list_percentage_arc, 10, "Histograma de cobertura sobre anotaciones en: ARCHAEA", "histogram_cover_arc.png")

# Virus

#list_percentage_vir = []
#with open("/scratch01/mnavarro/SUMMARIZED/virus_set.txt", "r") as set:
#    for i in set:
#        if re.search('^k',i):
#             list_percentage_vir.append(round(float(i.split()[2].split('%')[0])))
#hfunc.plot_histogram(list_percentage_vir, 10, "Histograma de cobertura sobre anotaciones en: VIRUS", "histogram_cover_vir.png")

#################################################################
# Option 3: from a .txt file with percentage values in the column
# with position 2 for the global matrix <summarized>, to each set
# of taxonomic groups on TLAHUICA.
#################################################################

#import hisfunc as hfunc

#def plot_coverage_histogram(set_name):
#    coverage = {}
#    summarized = open('./summarized.txt')
#    summ_lines = summarized.readlines()
#    for l in summ_lines:
#        el = l.split()
#        contig_id = el[0]
#        percentage = round(float(el[2][:-1]))
#        coverage[contig_id] = percentage
#    print coverage

#    list_percentage = []
#    with open('./' + set_name + '_set.txt', 'r') as set:
#        for i in set:
#            contig_id = i.split()[0]
#            list_percentage.append(coverage[contig_id])
    #hfunc.plot_histogram(list_percentage, 10, "Histograma de cobertura sobre anotaciones con ID", "histogram_cover.png")

#plot_coverage_histogram('bacteria')

#################################################################
# Option 3.1: from a .txt file with percentage values in the column
# with position 2 for the global matrix <summarized>, to whole sets
# of taxonomic groups on TLAHUICA.
#################################################################

#import hisfunc as hfunc

#hfunc.plot_histogram(hfunc.plot_coverage_histogram('bacteria'), 10, "Histograma de cobertura sobre anotaciones con ID", "histogram_cover_bac.png")
#hfunc.plot_histogram(hfunc.plot_coverage_histogram('eukaryota'), 10, "Histograma de cobertura sobre anotaciones con ID", "histogram_cover_euk.png")
