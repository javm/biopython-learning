import re

summarized = open('summarized.txt','r')
coverage = open('sum_coverage.txt','w')

A1 = []
for line in summarized:
    lines = line.strip()
    column = lines.split()
    A1.append(column)

A2 = []
for element in A1:
    A2.append(element[1])

total_len = []
for element in A2:
    num = element.replace('len=', '')
    num = float(num)
    total_len.append(num)

interval_len = []
for raw in A1:
    s = 0
    for col in raw:
        value = re.search('^\d+'+'_g;',col)
        if bool(value) == True:
            A = re.split('\|',col)
            for coor in A:
                B = re.split(';',coor)
                C = re.split('-',B[1])
                n1 = float(C[0])
                n2 = float(C[1])
                n3 = n2-n1
                s = s + n3
    interval_len.append(s)

count = 0
percen_arr = []
while count < len(total_len):
    percen = interval_len[count]/total_len[count]*100
    percen_arr.append(round(percen,2))
    count = count + 1

count = 0
while count < len(A1):
    a = str(percen_arr[count])
    b = a + '%'
    A1[count].insert(2,b)
    count = count + 1

for raw in A1:
    for col in raw:
        coverage.write(col + '\t')
    coverage.write('\n')
    
coverage.close()
