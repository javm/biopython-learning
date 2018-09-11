# --------------------------------------------------------------
# Marisol
# UNAM - IBt
# This file, contains several defined functions for plotting
# a histogram. It is called < hisfunc.py > and must be
# abbreviate as < hfunc.py >.
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
# Module p3 - computer
#################################################################

def plot_histogram(points, bucket_size, title="", file_name=''):
    histogram = make_histogram(points, bucket_size)
    plt.bar(histogram.keys(), histogram.values(), width=bucket_size, color='red')
    plt.title(title)
    plt.savefig('/Users/solouli/Desktop/Git_hub/Ngs_scripts/ngs_histogram_cover_py/%s' % file_name)
    plt.show()

#################################################################
# Module p3 - server tlahuica
#################################################################

#def plot_histogram(points, bucket_size, title="", file_name=''):
#    histogram = make_histogram(points, bucket_size)
#    plt.bar(histogram.keys(), histogram.values(), width=bucket_size, color='red')
#    plt.title(title)
#    plt.savefig('/scratch01/mnavarro/SUMMARIZED/%s' % file_name)
#    plt.show()
