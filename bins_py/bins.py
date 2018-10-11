#!/usr/bin/python2.7

# USAGE: python script.py parent_fasta.fasta <how many records per file>

import sys
import Bio

def write_file(input_file,split_number):
    #get file_counter and base name of fasta_file
    parent_file_base_name = input_file(".")[0]
    counter = 1

    #our first file name
    file = parent_file_base_name + "_" + str(counter) + ".fasta"

    #carries all of our records to be written
    joiner = []
    #enumerate huge fasta
    for num,record in enumerate(Bio.SeqIO.parse(input_file, "fasta"),start=1):
        #append records to our list holder
        joiner.append(">" + record.id + "\n" + str(record.seq))

        #if we have reached the maximum numbers to be in that file, write to a file, and then clear
        #record holder
        if num % split_number == 0:
            joiner.append("")
            with open(file,'w') as f:
                f.write("\n".join(joiner))

            #change file name,clear record holder, and change the file count
            counter += 1
            file = parent_file_base_name + "_" + str(counter) + ".fasta"
            joiner = []
      if joiner:
        joiner.append("")
        with open(file,'w') as f:
          f.write("\n".join(joiner))

if __name__ == "__main__":
    input_file = sys.argv[1]
    split_number = sys.argv[2]
    write_file(input_file,split_number)
    print "fasta_splitter.py is finished."
