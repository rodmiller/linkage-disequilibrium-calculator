#
#UP TO STEP 4 IN EMAIL..
#
'''
Y1. Take a file single file and read the first sequence (note each file has a header).

Y2. Take that sequence and count the numbers of G, C, A and T (note that dashes are gaps - these are either insertion or deletion mutations).

Y3. Count the number of each amino acid in the sequence (hint: use a dictionary)

Y4. Now read in all the sequences in a file (putting them into a list) and count the number of sites which show variation across the strains. E.g. if we had the folliowing sequences

seq1   TCTTGA
seq2   TCTTGC
seq3   TATTGA

we would find that two sites were variable or "polymorphic" across the strains.

5. Calculate the linkage disequilibrium (LD) between all pairs of polymorphic sites. To calculate LD we will use the method of Hill and Robertson, r^2.

r^2 = D^2 / (P(AB) P(aB) P(Ab) P(ab))

where D = P(AB) - P(A)P(B)

and P(AB)..etc are the frequencies of the haplotypes (chromosomes) carrying the A and B alleles, and P(A)...etc are the frequencies of the A allele at the A locus. I can describe this in more detail if you need.

5. Finally we have to divide the polymorphisms as to whether they are synonymous or non-synonymous (non-synonymous mutations change an amino acid, synonymous ones do not). We will want to look at patterns of LD for synonymous and non-synonymnous separately.
'''

import os #For looping through files in a dir and determining system os for file to load

#A dictionary for mapping nucleotide codons to their respective amino acids. 
aminoAcidCodons = {
	'Ile':['ATT','ATC','ATA'],
	'Leu':['CTT','CTC','CTA','CTG','TTA','TTG'],
	'Val':['GTT','GTC','GTA','GTG'],
	'Phe':['TTT','TTC'],
	'Met':['ATG'],
	'Cys':['TGT','TGC'],
	'Ala':['GCT','GCC','GCA','GCG'],
	'Gly':['GGT','GGC','GGA','GGG'],
	'Pro':['CCT','CCC','CCA','CCG'],
	'Thr':['ACT','ACC','ACA','ACG'],
	'Ser':['TCT','TCC','TCA','TCG','AGT','AGC'],
	'Tyr':['TAT','TAC'],
	'Trp':['TGG'],
	'Gln':['CAA','CAG'],
	'Asn':['AAT','AAC'],
	'His':['CAT','CAC'],
	'Glu':['GAA','GAG'],
	'Asp':['GAT','GAC'],
	'Lys':['AAA','AAG'],
	'Arg':['CGT','CGC','CGA','CGG','AGA','AGG'],
	'STOP':['TAA','TAG','TGA']
	}


def openFile(dir):
	'''Takes the path to a file and opens it, returning the file object'''
	#FOR DEBUG USING FILE ortholog_000000.nt_ali.Bifidobacteriumanimalissubsplactis.fasta
	f = open(dir)
	return f
	
def getSequences(fileObject):
	'''Takes the file object and dismantles it into its sequences, stripping the headers
	then returning a list of the sequences in the file'''
	#print fileObject
	seqNuc = ''
	seqHeader = ''
	sequences = []
	for line in fileObject:
		line = line.strip() #Removes leading and trailing whitespace
		#print str(i) + ' ' + line
		if line[0] == '>': #If the line is a header it starts with '>'
			if len(seqHeader) > 0: #If true then this variable has already been set from the previous sequence in the file
				sequences.append(seqNuc) #Appends the stored sequence to the list
				#print '%s --> %s' % (seqHeader,seqNuc)
				seqNuc = '' #Resets the sequence
				#seqHeader = '' #Resets the header
			else:
				seqHeader = line #If the line starts with > but the variable hasn't been set then its the first header
			#print seqHeader
		else: #If the line isn't a header its a line of nucleotide sequence, append it to the current working string
			seqNuc = seqNuc + line
	return sequences #Returns a list of sequences

def iterateSequences(sequences):
	'''Takes a list of sequences from a file, iterates through them counting nucleotides
	 and the amino acids in each. It returns a tuple in the format (nucleotides,amino acids)'''
	for sequence in sequences: #Iterates through each sequence
		countedNucs = countNucs(sequence) #Counts the nucleotides, returns a dictionary
		countedAminoAcids = countAminoAcids(sequence) #Counts the amino acids, returns a dictionary
		return countedNucs,countedAminoAcids #Returns a tuple of the two above results
		

def countNucs(seqNucs):
	'''Takes a sequence of nucleotides, iterates though the sequence and counts the
	occurrences of each nucleotide. Returns a tuple in the format
	 (counted nucleotides, point mutations) The counted nucleotides is itself a dictionary
	 with keys for each nucleotide type counted. Point mutations is a list of locations of
	 insertion or deletion mutations in the DNA sequence'''
	countedNucs = {}
	pointMutations = []
	countA = 0
	countT = 0
	countC = 0
	countG = 0
	i = 0 #For indexing point mutations '-'
	for nuc in seqNucs:
		if nuc == 'A':
			countA = countA + 1
		elif nuc == 'T':
			countT = countT + 1
		elif nuc == 'C':
			countC = countC + 1
		elif nuc == 'G':
			countG = countG + 1
		elif nuc == '-': #A known insertion or deletion mutation
			pointMutations.append(i) #Append the index of the mutation to the list
		else:
			print 'STRANGE NUCLEOTIDE --> ' + nuc
		i = i + 1
	countedNucs['a'] = countA #Add to the returned dictionary for each nucleotide type in lowercase
	countedNucs['t'] = countT
	countedNucs['c'] = countC
	countedNucs['g'] = countG
	return countedNucs,pointMutations #Returns a tuple

def findPointMutations(seqNucs):
	'''Takes a sequence of nucleotides, iterates through them and finds only the insertion
	or deletion mutations. Returns a list of their locations in the nucleotide sequence'''
	pointMutations = []
	i = 0 #For indexing point mutations
	for nuc in seqNucs:
		if nuc == '-': #A known insertion or deletion mutation
			pointMutations.append(i)
		i = i + 1
	return pointMutations #Returns a list of point mutation locations in the sequence
	
def countAminoAcids(seqCodons):
	'''Takes a list of codons. It maps each codon against the dictionary of amino acids 
	and their codons to find the amino acid. It then counts	the occurances and returns
	 a dictionary in the format amino acid three letter code:number counted'''
	countAminoAcid = {
		'Ile':0,'Leu':0,'Val':0,'Phe':0,'Met':0,'Cys':0,'Ala':0,'Gly':0,'Pro':0,'Thr':0,
		'Ser':0,'Tyr':0,'Trp':0,'Gln':0,'Asn':0,'His':0,'Glu':0,'Asp':0,'Lys':0,'Arg':0,
		'STOP':0
		} #Dictionary for counting each type of AA
	
	for codon in seqCodons: #For each codon in the list of codons
		for aminoAcid in aminoAcidCodons: #For each of the amino acids in the codon dictionary
			if codon in aminoAcidCodons[aminoAcid]: #If the codon represents an AA
				countAminoAcid[aminoAcid] = countAminoAcid[aminoAcid] + 1 #Add one to the respective amino acid
				#print codon + ' --> ' + aminoAcid
	#print '----------------------------------'
	return countAminoAcid #Returns a dictionary

def countVariation(sequences):
	'''Takes a list of sequences, iterates through them one by one and compares them to 
	all of the other sequences in turn and finds any differences. It then adds the location
	of these polymorphisms to a list and returns that.'''
	testsequences = ['TCTTGA','TCTTGC','TATTGA'] #Just a test sequence to see if it is working. Correct result is two polymorphisms at indexes 1 and 5.
	polymorphismIndexes = []
	for baseSequence in sequences: #For each sequence in the list to compare the others to
		for sequence in sequences: #Compares each sequence to the base sequence chosen above
			for x in range(len(sequence)): #For each character in the sequence, using x to allow the same character index to be compared among the base and test sequence
				#print sequence[x]
				if sequence[x] != baseSequence[x]: #Not equal == polymorphism
					#print 'Polymorphism at index %s between sequence %s and sequence %s' % (str(x), sequences.index(sequence), sequences.index(baseSequence))
					if x not in polymorphismIndexes: #If the polymorphism location is already noted then no point adding it again
						polymorphismIndexes.append(x)
	return polymorphismIndexes #Returns a list of the indexes of the polymorphisms

def getAllCodons(sequences):
	'''Takes a list of sequences, iterates through them and passes each one to getCodons().
	Then returns a list of a list of codons (e.g. [[AGG, GAG...],[TCG, AGC...]])
	'''
	allCodons = [] #Initiate list to append to
	for sequence in sequences: #For each sequence in the supplied list of sequences
		codons = getCodons(sequence) #Passes each sequence to getCodons() which returns a list of codons
		allCodons.append(codons) #Append the list of codons to the allCodons list
	return allCodons #Returns a list of a list of codons
	
def getAllAminoAcids(allCodons):
	'''Takes a list of list of codons (e.g. [[AGG, GAG...],[TCG, AGC...]]),
	passes each list to getAminoAcids() and appends it to a list.
	Returns a list of list of Amino Acids
	'''
	allAminoAcids = [] #Initiate list to append to
	for codons in allCodons: #For each list of codons in the list of list of codons
		aminoAcids = getAminoAcids(codons) #Passes each codon list to getAminoAcids() which returns a list of Amino Acids
		allAminoAcids.append(aminoAcids) #Append the list of amino acids to the allAminoAcids list
	return allAminoAcids #Returns a list of a list of Amino Acids

def getCodons(seqNuc):
	'''Takes a sequence string, iterates through it and splits the string into a list of
	codons. ORF starts from 1st character'''
	i = 0 #For counting through the codons
	seqCodons = [] #Initiate the list to appended to
	while i <= len(seqNuc): #Avoid Str index out of range errors
		codon = seqNuc[i:i+3] #Take a range of chars from i to i + 3. One codon
		#print codon
		seqCodons.append(codon) #Append the codon to the end of the list
		i = i + 3 #increment by 3 to move the ORF to the next codon
	return seqCodons #Returns a list of codons
	
def getAminoAcids(seqCodons):
	'''Takes a list of codons and maps them against their respective amino acids.
	Returns a list of Amino Acids'''
	aminoAcids = []
	for codon in seqCodons: #For each codon in the list of codons
		for aminoAcid in aminoAcidCodons: #For each of the amino acids in the codon dictionary
			if codon in aminoAcidCodons[aminoAcid]: #If the codon represents an AA
				aminoAcids.append(aminoAcid) #Append the corresponding AA to the list
	return aminoAcids #Returns a list of Amino Acid 3-character names

def calcLinkageDisequilibrium(sequence):
	'''5. Calculate the linkage disequilibrium (LD) between all pairs of polymorphic sites. To calculate LD we will use the method of Hill and Robertson, r^2.

	r^2 = D^2 / (P(AB) P(aB) P(Ab) P(ab))

	where D = P(AB) - P(A)P(B)

	and P(AB)..etc are the frequencies of the haplotypes (chromosomes) carrying the A and B alleles, and P(A)...etc are the frequencies of the A allele at the A locus. I can describe this in more detail if you need.
	'''

def getNonSynonymousPolymorphisms(allSeqCodons):
	'''Takes a list of a list of Amino Acids (i.e. [[Gln, Asp...],[Asp, Glu...]]).
	Then iterates through the lists to compare each sequence to another. 
	Finds polymorphisms only if non-synonymous and returns a dictionary in the format
	{polymorphic AA index: (1stAA,2ndAA)} NB value is a tuple.
	'''	
	nonSynPolymorphisms = {}
	for baseSeqCodons in allSeqCodons: #The first sequence to compare the second to
		for testSeqCodons in allSeqCodons: #The second sequence to compare the first to
			if baseSeqCodons != testSeqCodons: #Only if the sequences are different
				for x in range(len(baseSeqCodons)): #Prevent str index out of range errors
					print '%s - %s' % (len(baseSeqCodons), len(testSeqCodons))
					if baseSeqCodons[x] != testSeqCodons[x]: #Only print if the polymorphism is non-synonymous
						print '%s: %s -> %s' % (x, baseSeqCodons[x], testSeqCodons[x]) #Prints out the index: 1stAA->2ndAA
						changeAminoAcid = (baseSeqCodons[x], testSeqCodons[x])
						print changeAminoAcid
						#nonSynPolymorphisms[x: changeAminoAcid] #Appends to the dictionary
	return nonSynPolymorphisms #Returns a dictionary
		
def iterateFiles(dir):
	'''Takes the path to a directory and iterates over the files in it, passing each one
	to getSequences() and countVariation(). It then returns a string for the number of
	polymorphisms and their location within the file.'''
	files = os.listdir(dir) #Lists the files in the supplied dir. 
	for file in files: #Iterates through each file in the list of files
		#print dir + file
		fileObject = openFile(dir + file) #Opens the file chosen using openFile()
		sequences = getSequences(fileObject) #Breaks the file down into sequences
		#i = 0
		for sequence in sequences:
			if len(sequences[0]) != len(sequence):
				print 'AHH Different sequence lengths in %s' % file
			else:
				print 'Same Lengths in %s' % file
			print '%s -> %s' % (file, len(sequence))
		#for sequence in sequences:
		#	i = i + 1
		#	countedNucs = countNucs(sequence)
		#	print '%s - sequence %s --> %s' % (file[9:15], str(i), str(countedNucs[0]))
		#variation = countVariation(sequences) #Works out the variation for each sequence
		#allCodons = getAllCodons(sequences)
		#print allCodons
		#allAminoAcids = getAllAminoAcids(allCodons)
		#print allAminoAcids
		#nonSynPoly = getNonSynonymousPolymorphisms(allAminoAcids)
		#if nonSynPoly:
		#	print nonSynPoly 
		#iterateSequences(sequences)
		#print '-------------------------'
#		if len(variation)==1: #If only one polymorphism
#			print '%s --> %s polymorphism. Location is: %s' % (file[9:15], len(variation), str(variation))
##		elif len(variation)==0: #If no polymorphisms. Comment out to not list all the uninteresting non-polymorphic files
##			print '%s --> No Polymorphisms' % (file)
#		elif len(variation) > 1: #If more than one polymorphism
#			print '%s --> %s polymorphisms. Locations are: %s' % (file[9:15], len(variation), str(variation))
		
def main():
	if os.name == 'nt': #If on windows, open this tester file
		fileObject = openFile('F:\USER FILES\Dropbox\Dropbox\Biomedicine\Yr 3\FYP\Bifidobacterium animalis lactis\ortholog_000000.nt_ali.Bifidobacteriumanimalissubsplactis.fasta')
	elif os.name == 'posix': #If on mac, open this tester file
		fileObject = openFile('/Users/robert/Dropbox/Biomedicine/Yr 3/FYP/Bifidobacterium animalis lactis/ortholog_000000.nt_ali.Bifidobacteriumanimalissubsplactis.fasta')
	#sequences = getSequences(fileObject)
		
	print '-------------------------'
	#print fileObject
	#print countVariation(sequences)
	#print seqNuc
	#print seqHeader
	iterateFiles('/Users/robert/Dropbox/Biomedicine/Yr 3/FYP/Bifidobacterium animalis lactis/') #Calculates for all the files in the bifidobacterium file
	
if __name__ == '__main__':
	main()
	