# --------------------------------------------------------------
# Marisol
# UNAM - IBt
# This function file, contains several defined functions for
# plotting a histogram. It is called < hisfunc.py > and should
# be abbreviate as < hfunc.py >.
# --------------------------------------------------------------

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math
from collections import Counter

#################################################################
# Module p1
#################################################################

def bucketize(point, bucket_size):
      if bucket_size < 10 and bucket_size >= 4:
          return round(( math.floor( ( bucket_size * (point / bucket_size)  )*10)/10)*2)/2
      elif bucket_size < 4 and bucket_size >= 2:
          return round(( math.floor( ( bucket_size * (point / bucket_size)  )*10)/10)*3)/3
      elif bucket_size < 2:
          return  math.floor( ( bucket_size * (point / bucket_size)  )*10)/10
      else:
          return math.floor(bucket_size) * math.floor( point / bucket_size )

#################################################################
# Module p2
#################################################################

def make_histogram(points, bucket_size):
         return Counter(bucketize(float(x), bucket_size) for x in points)

#################################################################
# Module p3
#################################################################

def plot_histogram(points, bucket_size, title="", file_name=''):
    histogram = make_histogram(points, bucket_size)
    plt.bar(histogram.keys(), histogram.values(), width=bucket_size, color='#17becf')
    plt.title(title)
    #plt.savefig('/Users/solouli/Desktop/Git_hub/Git_scratch/biopython-learning/ngs_summarized_py/%s' % file_name)
    plt.savefig('./%s' % file_name)
    #plt.show()

#################################################################
# Module p4
#################################################################

def plot_coverage_histogram(set_name):
    coverage = {}
    summarized = open('./summarized.txt')
    summ_lines = summarized.readlines()
    for l in summ_lines:
        el = l.split()
        contig_id = el[0]
        percentage = round(float(el[2][:-1]))
        coverage[contig_id] = percentage
    print coverage

    list_percentage = []
    with open('./' + set_name + '_set.txt', 'r') as set:
        for i in set:
            contig_id = i.split()[0]
            list_percentage.append(coverage[contig_id])
    return list_percentage
