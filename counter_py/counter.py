#!/usr/bin/env python

import re

eucaryota = open("euka_reduced.genes", "r")
eucaryota_out = open("euka_count.genes", "w")

eucaryota_lines = eucaryota.readlines()


for i in range (0, len(eucaryota_lines)-1, 2):
    k = eucaryota_lines[i].strip()
    parts = k.split('_')
    kid = format(int(parts[1]), '07')

    whole_id = parts[0]+'_'+kid
    nuc = eucaryota_lines[i+1].strip()
    eucaryota_out.write(whole_id+'\n')
    eucaryota_out.write(nuc+'\n')
