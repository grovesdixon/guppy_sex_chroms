#!/usr/bin/env python

#Get_files
'''
'''

from sys import argv

debug = False

#set up usage information
usage = "fasta_sequence_characters.py -fa $stickleGenome > chromLengths.txt\nsubset_bams_by_chrom_size.py chromLengths.txt 4 *sorted.bam > subsetBams"
try:
    argv[1]
except IndexError:
    print('\nUsage:')
    exit(usage)




infileName = argv[1]
nSubs = int(argv[2])
outfile_name = argv[3]
bams = argv[4:]

if debug: #print to see if list is correct
    print FileList


with open(infileName, 'r') as infile:
	with open(outfile_name, 'w') as out:
		for line in infile:
			line=line.strip("\n").split()
			chrom = line[0]
			length = int(line[1])
			step=length/(nSubs)
			stepList = range(0, nSubs)
			sections = [x*step for x in stepList]
			sections.append(length)
			select_bams = []
			for b in bams:
				if "_chr{}_".format(chrom) in b:
					select_bams.append(b)
			print('---------')
			print(chrom)
			print('{} bam files'.format(len(select_bams)))
			print(select_bams)
			for b in select_bams:
				for i in range(len(sections)-1):
					outFile = b.split(".bam")[0] + "_{}_sub{}.bam".format(chrom, i+1) 
					out.write("samtools view -b {} -o {} {}:{}-{}\n".format(b, outFile, chrom, sections[i], (sections[i+1]-1)))



