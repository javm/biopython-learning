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

import re
import hisfunc as hfunc

list_percentage = []
with open("/Users/solouli/Desktop/Git_hub/Ngs_scripts/ngs_histogram_cover_py/example_set.txt", "r") as set:
    for i in set:
        if re.search('^k',i):
             list_percentage.append(round(float(i.split()[2].split('%')[0])))
             #print(list_percentage)

hfunc.plot_histogram(list_percentage, 10, "Histograma de cobertura sobre anotaciones", "histogram_cover.png")

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
