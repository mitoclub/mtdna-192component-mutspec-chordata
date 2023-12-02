import time
import timeit
import os

PATH_TO_DATA = 'data/abasic_sites/'

path_to_abasic_data = PATH_TO_DATA + '41467_2022_33594_MOESM26_ESM.csv'
genfilepath = PATH_TO_DATA + 'mm10_ChrM_oneline.fasta'
genfile = open(genfilepath,'r')

genfile.readline()
g = genfile.readline().strip()


def comp(nuc):
	if nuc == "A":
		return "T"
	elif nuc == "T":
		return "A"
	elif nuc == "G":
		return "C"
	elif nuc == "C":
		return "G"


nucs = ['A','T','G','C']
triplets = []
for i in nucs:
	for j in nucs:
		for k in nucs:
			triplets.append(i+j+k)


path3 = PATH_TO_DATA + "AbasicSitesMtDNAcontext_start_heavy.csv"
reps3 = open(path3, "w")

reps3.write("triplet;orig;count;deepcount;avg\n")


for tri in triplets:
	file = open(path_to_abasic_data,'r')
	file.readline()

	count = 0
	deepcount = 0
	while True:
		read = file.readline().strip()
		if read == "":
			break

		row = read.split(';')
		if row[4] == '-':
			#print(g[int(row[1])-1:int(row[1])+2] + " " + )
			if row[0] == "15423": # Exclude D-Loop
				print("The End")
				break

			if g[int(row[1])-1:int(row[1])+2] == tri:
				# print("match")
				count += 1
				deepcount += int(row[3])
	
	reps3.write("%s;%s;%s;%s;%s\n" % (comp(tri[2]).lower() + comp(tri[1]) + comp(tri[0]).lower(), tri[0].lower() + tri[1] + tri[2].lower(),count, deepcount, deepcount/count))



path3 = PATH_TO_DATA + "AbasicSitesMtDNAcontext_start_light.csv"
reps3 = open(path3, "w")
reps3.write("triplet;count;deepcount;avg\n")

for tri in triplets:

	file = open(path_to_abasic_data,'r')
	file.readline()

	count = 0
	deepcount = 0
	while True:
		read = file.readline().strip()
		if read == "":
			break

		row = read.split(';')
		if row[4] == '+':
			#print(g[int(row[1])-1:int(row[1])+2] + " " + )
			if row[0] == "15423": # Exclude D-Loop
				print("The End")
				break

			if g[int(row[1])-1:int(row[1])+2] == tri:
				# print("match")
				count += 1
				deepcount += int(row[3])

	print(tri)
	reps3.write("%s;%s;%s;%s\n" % (tri[0].lower() + tri[1] + tri[2].lower(), count, deepcount, deepcount/count))
