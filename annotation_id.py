headers = open ('trinotate_annotation_report.xls', 'r')
annotation = open ('annotation_out.xls','w')
A1 = []

headers_lines = headers.readlines()
for i in range(1, len(headers_lines)):
        line = headers_lines[i].strip()
        columns = line.split()
        A1.append(columns)

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
    annotation.write("\t".join(element)+'\n')

annotation.close()
