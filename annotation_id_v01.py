#!/usr/bin/env python

annotation = open ('annotation_out.txt','a')
A1 = []

counter = 0

for i in headers:
        line = i.strip()                                            #read lines
        column = line.split()                                     #read columns
        if counter != 0:
            A1.append(column)
        counter = counter + 1

#print(A1[0][0])

length_lines = len(A1)

#print(length_lines)

counter = 1

while counter <= length_lines:
    num = str(counter)
    id = num + '_g'
    A1[counter - 1][0] = id
    counter = counter + 1

for element in A1:
    for i in element:
        annotation.write(i + '\t')
    annotation.write('\n')
annotation.close()
