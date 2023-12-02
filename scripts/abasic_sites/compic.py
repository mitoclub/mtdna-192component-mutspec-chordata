
PATH_TO_DATA = 'data/abasic_sites/'


first = open(PATH_TO_DATA + 'AbasicSitesMtDNAcontext_start_heavy.csv','r')
first.readline().strip()

reps = open(PATH_TO_DATA + "AbasicSitesMtDNAcontext_HLcompare.csv", "w")
reps.write("tripletH;origtriH;countH;deepcountH;avgH;tripletL;countL;deepcountL;avgL\n")

while True:
	read1 = first.readline().strip()
	if read1 == "":
		break

	row1 = read1.split(';')

	second = open(PATH_TO_DATA + 'AbasicSitesMtDNAcontext_start_light.csv','r')
	second.readline().strip()

	while True:
		read2 = second.readline().strip()
		if read2 == "":
			break

		row2 = read2.split(';')
		if row2[0] == row1[0]:	
			reps.write("%s;%s\n" % (read1,read2))
