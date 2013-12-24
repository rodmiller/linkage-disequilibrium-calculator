#!/usr/bin/python
#coding: utf-8
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
import csv #For writing to csv files in order to make a graph in calc
import math #For doing math.floor
import matplotlib.pyplot as plt #For plotting the graphs
import sys #For getting command line arguments
import time #For sleeping the program for a few secs for emphasis
import pynotify #For notifying user every time a bacteria finishes

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
	'STOP':['TAA','TAG','TGA'],
	'NONE':['---']
	}


def openFile(dir):
	'''Takes the path to a file and opens it, returning the file object'''
	#FOR DEBUG USING FILE ortholog_000000.nt_ali.Bifidobacteriumanimalissubsplactis.fasta
	#print dir[-5:]
	if dir[-5:] == 'fasta':
		f = open(dir)
		#print dir
		return f
	else:
		return
	
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
	countU = 0
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
		elif nuc == 'U':
			countU = countU + 1
		elif nuc == '-': #A known insertion or deletion mutation
			pointMutations.append(i) #Append the index of the mutation to the list
		else:
			print 'STRANGE NUCLEOTIDE --> ' + nuc
		i = i + 1
	countedNucs['a'] = countA #Add to the returned dictionary for each nucleotide type in lowercase
	countedNucs['t'] = countT
	countedNucs['c'] = countC
	countedNucs['g'] = countG
	if countU > 0:
		countedNucs['u'] = countU
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
	countAminoAcid = {}
	for aminoAcidCodon in AminoAcidCodons:
		countAminoAcid[aminoAcidCodon] =  0 #Dictionary for counting each type of AA
	
	for codon in seqCodons: #For each codon in the list of codons
		for aminoAcid in aminoAcidCodons: #For each of the amino acids in the codon dictionary
			if codon in aminoAcidCodons[aminoAcid]: #If the codon represents an AA
				countAminoAcid[aminoAcid] = countAminoAcid[aminoAcid] + 1 #Add one to the respective amino acid
				#print codon + ' --> ' + aminoAcid
	#print '----------------------------------'
	return countAminoAcid #Returns a dictionary

def getVariation(sequences):
	'''Takes a list of sequences, iterates through them one by one and compares them to 
	all of the other sequences in turn and finds any differences. It then adds the location
	of these polymorphisms to a list and returns that.'''
	testsequences = ['TCTTGA','TCTTGC','TATTGA'] #Just a test sequence to see if it is working. Correct result is two polymorphisms at indexes 1 and 5.
	polymorphismIndexes = []
	for baseSequence in sequences: #For each sequence in the list to compare the others to
		for sequence in sequences: #Compares each sequence to the base sequence chosen above
			for x in range(len(sequence)): #For each character in the sequence, using x to allow the same character index to be compared among the base and test sequence
				#print sequence[x]
				if sequence[x] != baseSequence[x] and sequence[x] != '-' and baseSequence[x] != '-': #Not equal == polymorphism
					#print 'Polymorphism at index %s between sequence %s and sequence %s' % (str(x), sequences.index(sequence), sequences.index(baseSequence))
					if x not in polymorphismIndexes: #If the polymorphism location is already noted then no point adding it again
						polymorphismIndexes.append(x) #Adds it to a list of the polymorphism indexes in the file
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
	codons. ORF starts from 0th character'''
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

def calcLinkageDisequilibrium(variation, sequences, file):
	'''Takes a list of polymorphism locations in the sequence data and the sequence data,
	iterates through them and compares them pairwise to one another.
	It assigns alleles to nucleotides at the sites, with the 1st nucleotide coming across
	being bigA and the second being	littleA. Then it calculates the LD between the sites using r^2.
	
	r^2 = D^2 / (P(AB) P(aB) P(Ab) P(ab))

	where D = P(AB) - P(A)P(B)
	NEEDS REWRITE ASAP
	'''
	#testSeq = ['ATAACA', 'AAAAGA', 'ATAACA', 'AAAACA']
	#testVariation = [1, 4]
	#for debug only...
	#variation = testVariation
	#sequences = testSeq
	
	for i in range(len(variation)): #Use i to be sure we only compare polymorphisms occurring after the current one in the list
		basePolyList = [] #Initialise list to be used to compare the nucleotides at the first polymorphic site
		for sequence in sequences: #For each sequence in the list of sequences
			#print sequence[variation[i]]
			if sequence[variation[i]] not in basePolyList and sequence[variation[i]] != '-': #Prevent duplicates of the same nucleotide being added and prevent point mutations being added
				basePolyList.append(sequence[variation[i]]) #Append the polymorphic nucleotide to the list
				#print sequence[variation[i]] + ' Added!'
		#print '------- Base Poly is: ' + str(basePolyList)
		#print 'Len of Base Poly is: ' + str(len(basePolyList))
		#if len(basePolyList) != 2: #Only interested if there are 2 possible nucleotides, if more than 2 then its too complicated if less than then its not polymorphic!
		#	print 'Too many polymorphisms at this site to work with'
		#	break #Break back to next polymorphism index in file
		#else:
		#print basePolyList
		bigA = basePolyList[0] #First one to come across is bigA
		for basePoly in basePolyList:
			if basePoly != bigA:
				littleA = basePoly

		#littleA = basePolyList[1] #Second one is littleA
		#print bigA + ' ' + littleA
		for x in range(i, len(variation)): #Use x here to iterate through the list, using i as the start in the range to prevent comparing preceding polymorphisms
			#print x
			testPolyList = [] #Initialise list to be used to compare the nucleotides at the second polymorphic site
			for sequence in sequences: #For each sequence in the list of sequences
				if sequence[variation[x]] not in testPolyList and sequence[variation[i]] != '-': #Prevent duplicates of the same nucleotide being added and prevent point mutations being added
					testPolyList.append(sequence[variation[i]]) #Append the polymorphic nucleotide to the list
			if len(testPolyList) != 2: #Only interested if there are 2 possible nucleotides, if more than 2 then its too complicated if less than then its not polymorphic!
				#print 'Too many polymorphisms at this site to work with'
				break #Break back to next polymorphism index in file
			else:
				if basePolyList != testPolyList: #Only print if bigA and bigB are different or littleA and littleB are different
					print '----------- File is: ' + file + ' and index is: ' + str(variation[i])
					print 'Base Poly is: ' + str(basePolyList)
					print 'Test Poly is: ' + str(testPolyList)
				bigB = testPolyList[0] #First one to come across is bigB
				littleB = testPolyList[1] #Second one is littleB
			#Finished defining the nucleotides now moving to define the haplotype frequencies
			if bigA == bigB:
				pAB = 1.0
			if bigA == littleB:
				pAb = 1.0
			if littleA == bigB:
				paB = 1.0
			if littleA == littleB:
				pab = 1.0
			print 'A = %s' % bigA
			print 'B = %s' % bigB
			print 'b = %s' % littleB
			print 'a = %s' % littleA
			print basePolyList

def newCalcLinkageDisequilibrium(variation, sequences):
	'''Takes a list of the Non Synonymous Polymorphisms and a list of the Synonymous Polymorphisms'''
	variationNucleotides = {} #Initialise the dictionary used to collate the possible nucleotides at each site. In the format {polymorphism index: [nucleotide 1, nucleotide 2], ...}
	allPossibleNucleotides = {} #Initialise the dictionary used to list all the nucleotides at each site without ommitting duplicates
	#sequences = ['CA', 'CG', 'TG', 'CG', 'TG', 'TA']
	#variation = [0, 1]
	distanceAndLdSyn = {}
	distanceAndLdNonSyn = {}
	for x in range(len(variation)): #Prevent str index out of range errors in variation and maintain same index throughout
		variationNucleotides[variation[x]] = [] #Initialise the list for each polymorphism index encountered
		for sequence in sequences: #Then go through each sequence in turn
			
			#if sequence[variation[x]] not in variationNucleotides[variation[x]]: #Only add it to the list if the nucleotide isn't already there
			variationNucleotides[variation[x]].append(sequence[variation[x]]) # Add it to the list
		for variationNucleotide in variationNucleotides: #Iterate each nucleotide
			currentNuc = [] #A list for the current Nucleotides
			if variationNucleotide not in currentNuc:
				currentNuc.append(variationNucleotides[variationNucleotide])
			if len(currentNuc) > 2: #Remove complex cases w/ more than 2 nucleotides at each site
				del(variationNucleotide) # If len >2 delete the entry
	#print variationNucleotides #Print out the dictionary
	#Now finished creating the dictionary of the polymorphisms. Now iterate through it and do pairwise comparisons
	allSeqCodons = getAllCodons(sequences) #Get the output of allSeqCodons
	allSeqAA = getAllAminoAcids(allSeqCodons) #Get the output of allSeqAA
	synPoly = getSynonymousPolymorphisms(allSeqCodons) #Get the output of Synpoly
	nonSynPoly = getNonSynonymousPolymorphisms(allSeqAA) #Get the output of nonSynPoly
	for basePolyIndex in variationNucleotides:
		for testPolyIndex in variationNucleotides:
			ldEquationR2 = None
			bigA = None
			del(bigA)
			bigB = None
			del(bigB)
			littleA = None
			del(littleA)
			littleB = None
			del(littleB)
			if testPolyIndex == basePolyIndex:
				break
			
			
			distanceBetween = abs(abs(testPolyIndex) - abs(basePolyIndex))
			#if distanceBetween > 1000:
			#	
			##	print testPolyIndex
			#	print basePolyIndex
			#	print '---' + str(distanceBetween)
			for baseNuc in variationNucleotides[basePolyIndex]:
				try:
					bigA
				except NameError:
					#print 'Setting BigA to %s' % baseNuc
					bigA = baseNuc
				else:
					if baseNuc != bigA:
						littleA = baseNuc
				#print '-----------bigA Chosen %s, baseNuc is %s. From %s' % (bigA, baseNuc, variationNucleotides[basePolyIndex])
				for testNuc in variationNucleotides[testPolyIndex]:
					try:
						bigB
					except:
						bigB = testNuc
					else:
						#print bigB
						if testNuc != bigB:
							littleB = testNuc
					#print 'bigB Chosen %s, testNuc is %s. From %s' % (bigB, testNuc, variationNucleotides[testPolyIndex])	
			#Now work out some probabilities
			basePolyIndexTotal = len(variationNucleotides[basePolyIndex])
			testPolyIndexTotal = len(variationNucleotides[testPolyIndex])
			#print '%s %s' % (basePolyIndexTotal, testPolyIndexTotal)
			#print variationNucleotides
			bigACount = variationNucleotides[basePolyIndex].count(bigA)
			littleACount = variationNucleotides[basePolyIndex].count(littleA)
			bigBCount = variationNucleotides[testPolyIndex].count(bigB)
			littleBCount = variationNucleotides[testPolyIndex].count(littleB)
			#print '%s %s/%s %s %s/%s' % (bigACount, littleACount,basePolyIndexTotal, bigBCount, littleBCount, testPolyIndexTotal)
			pBigA = float(bigACount)/float(basePolyIndexTotal)
			pLittleA = float(littleACount)/float(basePolyIndexTotal)
			pBigB = float(bigBCount)/float(testPolyIndexTotal)
			pLittleB = float(littleBCount)/float(testPolyIndexTotal)
			if pBigA + pLittleA != 1 or pBigB + pLittleB != 1: #If they don't equal one then we're working with something which has more than one possible pLittleX allele, ignore these.
				break
			#Now moving onto the more complex 2 part probabilities
			#Finding pAB
			countAB = 0
			for i in range(basePolyIndexTotal):
				baseNuc = variationNucleotides[basePolyIndex][i]
				if baseNuc == bigA:
					testNuc = variationNucleotides[testPolyIndex][i]
					if testNuc == bigB:
						countAB = countAB + 1
			pAB = float(countAB)/float(basePolyIndexTotal)
			#print 'pAB is ' + str(pAB)
			#Finding pAb
			countAb = 0
			for i in range(basePolyIndexTotal):
				baseNuc = variationNucleotides[basePolyIndex][i]
				if baseNuc == bigA:
					testNuc = variationNucleotides[testPolyIndex][i]
					#print 'BaseNuc is %s, testNuc is %s' % (baseNuc, testNuc)
					if testNuc == littleB:
						countAb = countAb + 1
			#print 'countAb is ' + str(countAb)
			pAb = float(countAb)/float(basePolyIndexTotal)
			#print 'pAb is ' + str(pAb)
			#Finding paB
			countaB = 0
			for i in range(basePolyIndexTotal):
				baseNuc = variationNucleotides[basePolyIndex][i]
				if baseNuc == littleA:
					testNuc = variationNucleotides[testPolyIndex][i]
					if testNuc == bigB:
						countaB = countaB + 1
			paB = float(countaB)/float(basePolyIndexTotal)
			#print 'paB is ' + str(paB)
			#Finding pab
			countab = 0
			for i in range(basePolyIndexTotal):
				baseNuc = variationNucleotides[basePolyIndex][i]
				if baseNuc == littleA:
					testNuc = variationNucleotides[testPolyIndex][i]
					if testNuc == littleB:
						countab = countab + 1
			pab = float(countab)/float(basePolyIndexTotal)
			#print 'pab is ' + str(pab)

			#Work out if its synonymous or nonsynonymous
			
			synOrNonSyn =  chooseSynOrNonSynPoly(basePolyIndex, sequences, synPoly, nonSynPoly)
			#print synOrNonSyn
			if synOrNonSyn == 'syn':
				#Now do the calculation
				ldEquationD = float(pAB - (pBigA * pBigB))
				ldEquationR2Bottom = float(pBigA * pLittleA * pBigB * pLittleB)
				if ldEquationR2Bottom != 0:
					ldEquationR2 = (ldEquationD**2.0) / ldEquationR2Bottom
					if ldEquationR2 > 1.1:
						print '----\nld = %s; \n pA=%s; \n pa=%s; \n pB=%s; \n pb=%s; \n pAB=%s; \n pAb=%s; \n paB=%s; \n pab=%s; \n d=%s; \n bottom=%s;' % (ldEquationR2, pBigA, pLittleA, pBigB, pLittleB, pAB, pAb, paB, pab, ldEquationD, ldEquationR2Bottom)
						print '%s %s/%s %s %s/%s' % (bigACount, littleACount,basePolyIndexTotal, bigBCount, littleBCount, testPolyIndexTotal)
						print 'base: %s test:%s' % (variationNucleotides[basePolyIndex], variationNucleotides[testPolyIndex])
					if ldEquationR2 != 1.0:
						try:
							distanceAndLdSyn[float(abs(distanceBetween))]
						except:
							distanceAndLdSyn[float(abs(distanceBetween))] = [float(ldEquationR2)]
						else:
							distanceAndLdSyn[float(abs(distanceBetween))].append(float(ldEquationR2))
					
			elif synOrNonSyn == 'nonSyn':
				#Now do the calculation
				ldEquationD = float(pAB - (pBigA * pBigB))
				ldEquationR2Bottom = float(pBigA * pLittleA * pBigB * pLittleB)
				if ldEquationR2Bottom != 0:
					ldEquationR2 = (ldEquationD**2.0) / ldEquationR2Bottom
					if ldEquationR2 > 1.1:
						print '----\nld = %s; \n pA=%s; \n pa=%s; \n pB=%s; \n pb=%s; \n pAB=%s; \n pAb=%s; \n paB=%s; \n pab=%s; \n d=%s; \n bottom=%s;' % (ldEquationR2, pBigA, pLittleA, pBigB, pLittleB, pAB, pAb, paB, pab, ldEquationD, ldEquationR2Bottom)
						print '%s %s/%s %s %s/%s' % (bigACount, littleACount,basePolyIndexTotal, bigBCount, littleBCount, testPolyIndexTotal)
						print 'base: %s test:%s' % (variationNucleotides[basePolyIndex], variationNucleotides[testPolyIndex])
					if ldEquationR2 != 1.0:
						try:
							distanceAndLdNonSyn[float(abs(distanceBetween))]
						except:
							distanceAndLdNonSyn[float(abs(distanceBetween))] = [float(ldEquationR2)]
						else:
							distanceAndLdNonSyn[float(abs(distanceBetween))].append(float(ldEquationR2))
	return (distanceAndLdSyn, distanceAndLdNonSyn)
			
def chooseSynOrNonSynPoly(variationIndex, sequences, synPolymorphisms, nonSynPolymorphisms):
	#seqAllCodons = getAllCodons(sequences)
	#seqAllAA = getAllAminoAcids(seqAllCodons)
	variationCodonIndex = math.floor(float(variationIndex)/3)
	#print '-------------------'
	#print 'SynPoly = ' + str(synPolymorphisms)
	#print 'NonSynPoly = ' + str(nonSynPolymorphisms)
	#print variationCodonIndex
	if int(variationCodonIndex) in synPolymorphisms:
		#print 'Its a Synonymous Poly'
		return 'syn'
	elif int(variationCodonIndex) in nonSynPolymorphisms:
		#print 'Its a Non Synonymous Poly'
		return 'nonSyn'
	else:
		#print 'Couldnt find the poly in the list :('
		return None

	


def getNonSynonymousPolymorphisms(allSeqAA):
	'''Takes a list of a list of Amino Acids (i.e. [[Gln, Asp...],[Asp, Glu...]]).
	Then iterates through the lists to compare each sequence to another. 
	Finds polymorphisms only if non-synonymous and returns a dictionary in the format
	{polymorphic AA index: (1stAA,2ndAA)} NB value is a tuple.
	'''	
	nonSynPolymorphisms = {}
	for baseSeqAA in allSeqAA: #The first sequence to compare the second to
		for testSeqAA in allSeqAA: #The second sequence to compare the first to
			if baseSeqAA != testSeqAA: #Only if the Amino Acids are different
				for x in range(min(len(baseSeqAA), len(testSeqAA))): #Prevent str index out of range errors
					#print '%s - %s' % (len(baseSeqAA), len(testSeqAA))
					if baseSeqAA[x] != testSeqAA[x]: #Only print if the polymorphism is non-synonymous
						#print '%s: %s -> %s' % (x, baseSeqAA[x], testSeqAA[x]) #Prints out the index: 1stAA->2ndAA
						changeAminoAcid = (baseSeqAA[x], testSeqAA[x]) #The tuple to append to the dictionary (1stAA, 2ndAA)
						#print changeAminoAcid
						nonSynPolymorphisms[x] = changeAminoAcid #Appends to the dictionary
	return nonSynPolymorphisms #Returns a dictionary

def getSynonymousPolymorphisms(allSeqCodons):
	'''Takes a list of a list of codons (i.e. [[AGC, TAG...],[TGC, CGT...]])
	then iterates through the list comparing each list of codons to another. 
	It finds if the codons are different then compares them to the amino acid dictionary
	to ensure they code for the same AA and returns only Synonymous Polymorphisms.
	Returns a dictionary in the format {polymorphic codon index: (1st Codon, 2nd Codon)}
	'''
	synPolymorphisms = {} #Initialise the dictionary
	for baseSeqCodons in allSeqCodons: #The first codons to compare the others to
		for testSeqCodons in allSeqCodons: #The codons to compare the first to
			for x in range(len(baseSeqCodons)): #Iterate through the list of codons using x to maintain the same index
				if baseSeqCodons[x] != testSeqCodons[x]: #Only if the codons are different
					#print baseSeqCodons[x] + ' ' + testSeqCodons[x]
					for aminoAcid in aminoAcidCodons: #Iterate through the dictionary to find the AA the codon codes for
						#print aminoAcid
						if baseSeqCodons[x] in aminoAcidCodons[aminoAcid]: #If True we have found the AA for this codon
							#print 'Found an AA for base codon'
							#print aminoAcid
							baseSeqCodonAA = aminoAcid
							if testSeqCodons[x] in aminoAcidCodons[aminoAcid]: #The codon represents an amino acid and the two codons represent the same one
								#print 'Found the same codon, its ' + aminoAcid
								synPolymorphisms[x] = (baseSeqCodons[x], testSeqCodons[x]) #Append it to the dictionary
	return synPolymorphisms	

def getLinkageDisequilibrium(synPolymorphisms, nonSynPolymorphisms, allSeqCodons):
	'''Takes a list of synonymous polymorhisms and a list of non-synonymous
	polymorhisms and the sequence data. Then iterates through each list in
	turn, finding the exact polymorphic nucleotide and then adding it to a 
	list for the sequence'''
	for seqCodons in allSeqCodons: #Get the individual sequence of codons in the file
		for synPolyIndex in synPolymorhisms: #Get each index for the polymorphisms 
			print 'Does nothing'
		print 'Does Nothing'
	return


		
def categorisePolymorphisms(allSeqCodons):
	'''Takes a list of a list of codons (i.e. [[AGC, TAG...],[TGC, CGT...]])
	then iterates through comparing each list of codons to another.
	If the codons are different it compares them in the Amino Acid dictionary. 
	If they code for the same AA they are Synonymous and are added to the synPolymorphisms
	dictionary in the format {polymorphic codon index: (Amino Acid, 1st Codon, 2nd Codon)}
	If they code for different AAs they are Non Synonymous and are added to the
	nonSynPolymorphisms dictionary in the format {polymorphic codon index: (1st AA, 2nd AA)}
	Returns a tuple in the format (synPolymorphisms, nonSynPolymorphisms)
	'''
	#testData = [['ATT', 'ATT', 'ATT', 'ATT'], ['ATC', 'ATG', 'ATT', 'ATT']]
	#allSeqCodons = testData
	synPolymorphisms = {} #Initialise the Synonymous Polymorphisms dictionary
	nonSynPolymorphisms = {} #Initialise the Non Synonymous Polymorphisms dictionary
	for baseSeqCodons in allSeqCodons: #The first codons to compare the others to
		for testSeqCodons in allSeqCodons: #Codons to compare the first to
			#print 'Comparing %s to %s' % (baseSeqCodons, testSeqCodons)
			for x in range(len(baseSeqCodons)): #Iterate through using x to maintain the same index
				if baseSeqCodons[x] != testSeqCodons[x]: #The codons are different therefore polymorphism
					#print 'Difference found at index %s between %s and %s' % (x, baseSeqCodons[x], testSeqCodons[x])
					for aminoAcid in aminoAcidCodons: #Iterate through the Amino Acid dictionary
						if baseSeqCodons[x] in aminoAcidCodons[aminoAcid]: #If True we have found the AA for this codon
							baseSeqCodonAA = aminoAcid
							#print 'BASE: Found an Amino Acid for %s. Its %s' % (baseSeqCodons[x], baseSeqCodonAA)
							if testSeqCodons[x] in aminoAcidCodons[aminoAcid]: #if True the codon represents an amino acid and the two codons represent the same one therefore synonymous
								#print 'TEST: Its the Same Amino Acid in the test seq. Its %s' % aminoAcid
								synPolymorphisms[x] = (aminoAcid, baseSeqCodons[x], testSeqCodons[x]) #Append to the synonymous polymorphisms dictionary
							else: #If above statement is False, it could be Non Synonymous
								#print 'TEST: Its a different Amino Acid in the test seq.'
								for aminoAcid in aminoAcidCodons:
									if testSeqCodons[x] in aminoAcidCodons[aminoAcid]:
										testSeqCodonAA = aminoAcid
										nonSynPolymorphisms[x] = (baseSeqCodonAA, testSeqCodonAA)
	return (synPolymorphisms, nonSynPolymorphisms)

def iterateFolders(dir, ext):
	'''Takes the path to the sequence data folder and finds folders in this folder for each bacterium.
	Then calls iterateFiles on each one of those folders'''
	files = os.listdir(dir)
	doneBacteria = []
	valuesSyn = {}
	valuesNonSyn = {}
	totalNumBacteria = 0
	sequenceDataDir = "/home/robert/FYP/Sequence Data/goodBacteria/"
	osWalk = os.walk(sequenceDataDir)
	global allBacteria
	allBacteria = [x[0].rsplit('/', 1)[1] for x in osWalk]
	allBacteria.pop(0)
	print allBacteria
	
	for file in files:
		if '-a' in sys.argv:
			#Doing average graphs...
			#print file[:7]
			if file[:7] == 'Average':
				print file[8:-4]
				if file[8:-4] not in doneBacteria:
					doneBacteria.append(file[8:-4])
		else:
			if file[-3:] == ext or file[-3:] == 'pdf' or file[-3:] == 'png':
				#pdfPath = '"' + dir + file + '"'
				#dropboxUpload(pdfPath, '"/Biomedicine/Yr 3/FYP/Results for Adam/"')
				
				if file[:-4] not in doneBacteria:
					doneBacteria.append(file[:-4])
				
			else:
				totalNumBacteria = totalNumBacteria + 1
	#print 'Already done: ' + str(doneBacteria)
	if '-forcegraphs' in sys.argv:
		numToDo = len(allBacteria)
		doneBacteria = []
	else:
		numToDo = len(allBacteria) - len(doneBacteria)
		print doneBacteria
		print '-----------------' + str(numToDo)
		
	
	print 'Going to do %s bacteria. Done %s and total is %s' % (str(numToDo), len(doneBacteria), totalNumBacteria)
	#print doneBacteria
	currentNumber = 0
	for file in files:
		if file not in doneBacteria and file[-3:] != ext and file[-3:] != 'pdf' and file[-3:] != 'png' or '-forcegraphs' in sys.argv:
			
			if file[-3:] != ext and file[-3:] != 'pdf' and file[-3:] != 'png':
				currentNumber = currentNumber + 1
				print '----- ' + file.upper() + ' On file %s of %s ----- ' % (currentNumber, numToDo)
				#n = pynotify.Notification("Now doing %s \n On file %s of %s" % (file, currentNumber, numToDo))
				#n.show()
				dataPath = dir + file + '/'
				bacteriaFiles = os.listdir(dataPath)
				if '-a' in sys.argv:
					if 'averageSynDistance.txt' in bacteriaFiles:
						print 'Loading average data from file...'
						#n = pynotify.Notification("%s - Loading average data" % file)
						#n.show()
						averageData = loadAverageData(file)
						drawAverageGraph(averageData[0], averageData[1], averageData[2], averageData[3], file, ext)
					else:
						print 'Nothing calculated for this bacteria, calculating now...'
						iterateFiles(dataPath, file, ext)
				else:
					if 'distanceValuesNonSyn.txt' in bacteriaFiles and 'distanceValuesSyn.txt' in bacteriaFiles and 'ldValuesNonSyn.txt' in bacteriaFiles and 'ldValuesSyn.txt' in bacteriaFiles and '-forcecalc' not in sys.argv:
						print 'LD and distances already calculated, loading now...'
						f = open(dataPath + 'distanceValuesNonSyn.txt')
						distanceValuesNonSyn = eval(f.read())
						f.close()
						f = open(dataPath + 'distanceValuesSyn.txt')
						distanceValuesSyn = eval(f.read())
						f.close()
						f = open(dataPath + 'ldValuesNonSyn.txt')
						ldValuesNonSyn = eval(f.read())
						f.close()
						f = open(dataPath + 'ldValuesSyn.txt')
						ldValuesSyn = eval(f.read())
						f.close()
						print ' Done.'
						drawGraph(distanceValuesSyn, ldValuesSyn, distanceValuesNonSyn, ldValuesNonSyn, file, ext)
					else:
						#n = pynotify.Notification("%s - Calculating LD" % file)
						#n.show()
						iterateFiles(dataPath, file, ext)		
	return

def iterateFiles(dir, bacteriaName, ext):
	'''Takes the path to a directory and iterates over the files in it, passing each one
	to getSequences() and getVariation(). It then returns a string for the number of
	polymorphisms and their location within the file.'''
	files = os.listdir(dir) #Lists the files in the supplied dir. 
	#print files
	#os.remove('./syn.csv')
	#os.remove('./nonSyn.csv')
	#f = open('./syn.csv', 'w')
	#f.close()
	#f = open('./syn.csv', 'a')
	#f2 = open('./nonSyn.csv', 'w')
	#f2.close()
	#f2 = open('./nonSyn.csv', 'a')
	distanceValuesSyn = []
	ldValuesSyn = []
	distanceValuesNonSyn = []
	ldValuesNonSyn = []
	valuesSyn = {}
	valuesNonSyn = {}
	for file in files: #Iterates through each file in the list of files
		#print file
		fileObject = openFile(dir + file) #Opens the file chosen using openFile()
		if fileObject:
			#print '--------- %s' % file
			sequences = getSequences(fileObject) #Breaks the file down into sequences
			variation = getVariation(sequences) #Works out the variation for each sequence
			#print variation
			#allCodons = getAllCodons(sequences)
			#print allCodons
			#allAminoAcids = getAllAminoAcids(allCodons)
			#print allAminoAcids
			#nonSynPoly = getNonSynonymousPolymorphisms(allAminoAcids)
			#synPoly = getSynonymousPolymorphisms(allCodons)
			#if variation and len(variation) != len(nonSynPoly)+len(synPoly):
			#	print '-------------------------'
			#	print '%s => %s variation and %s polys' % (file, len(variation), str(len(nonSynPoly)+len(synPoly)))
			#	print 'Variation: %s' % variation
			#	print 'nonSynPoly: %s' % nonSynPoly
			#	print 'synPoly: %s' % synPoly
			#print file
			distanceAndLd = newCalcLinkageDisequilibrium(variation, sequences)
			print distanceAndLd
			for k,v in distanceAndLd[0].iteritems():
				for value in v:
					#print 'Distance Type: %s, LD Type: %s' % (k, v)
					#print value
					#if '-a' in sys.argv:
					#	try:
					#		valuesSyn[k]
					#	except:
					#		valuesSyn[k] = v
					#	else:
					#		valuesSyn[k].append(v)
					#else:
					distanceValuesSyn.append(k)
					ldValuesSyn.append(value)
					#print 'Syn: %s => %s' % (str(k), str(value))
					#f.write(str(k) + ',' + str(value) + '\n')
			
			for k,v in distanceAndLd[1].iteritems():
				for value in v:
					#print 'Distance Type: %s, LD Type: %s' % (k, v)
					#print value
					#if '-a' in sys.argv:
					#	try:
					#		valuesNonSyn[k]
					#	except:
					#		valuesNonSyn[k] = v
					#	else:
					#		for individualLdValue in v:
					#			valuesNonSyn[k].append(individualLdValue)
					#else:
					distanceValuesNonSyn.append(k)
					ldValuesNonSyn.append(value)
					#print 'nonSyn: %s => %s' % (str(k), str(value))
					#f2.write(str(k) + ',' + str(value) + '\n')
			sys.stdout.write("%s => Syn: %s NonSyn: %s \r" % (str(file), str(len(ldValuesSyn)), str(len(ldValuesNonSyn))))
			sys.stdout.flush()
			#print 'Syn: %s' % len(ldValuesSyn)
			#print 'Non Syn: %s' % len(ldValuesNonSyn)
			
			
			
				#print str(k) + ',' + str(v)
				
			
	#		if len(variation)==1: #If only one polymorphism
	#			print '%s --> %s polymorphism. Location is: %s' % (file[9:15], len(variation), str(variation))
	##		elif len(variation)==0: #If no polymorphisms. Comment out to not list all the uninteresting non-polymorphic files
	##			print '%s --> No Polymorphisms' % (file)
	#		elif len(variation) > 1: #If more than one polymorphism
	#			print '%s --> %s polymorphisms. Locations are: %s' % (file[9:15], len(variation), str(variation))
	#print 'Syn: %s , %s' % (str(distanceValuesSyn), str(ldValuesSyn))
	#print 'Non Syn: %s, %s' % (str(distanceValuesNonSyn), str(ldValuesNonSyn))
	#plt.plot(distanceValuesSyn, ldValuesSyn, 'ro')
	#plt.plot(distanceValuesNonSyn, ldValuesNonSyn, 'bo')
	#plt.ylabel('Linkage Disequilibrium')
	#plt.xlabel('Distance (nucleotides)')
	#plt.show()
	#while True:
	#	print 'Current limits are ' + str(plt.axis())
	#	xmin = raw_input('Minimum on X-axis--> ')
	#	xmax = raw_input('Maximum on X-axis--> ')
	#	ymin = raw_input('Minimum on Y-axis--> ')
	#	ymax = raw_input('Maximum on Y-axis--> ')
	#	print '%s, %s, %s, %s' % (xmin, xmax, ymin, ymax)
	#	yorn = raw_input('Happy with these values? (y or n)')
	#	if yorn[0] == 'y':
	#		break
	#plt.plot(distanceValuesSyn, ldValuesSyn, 'ro')
	#plt.plot(distanceValuesNonSyn, ldValuesNonSyn, 'bo')
	#plt.ylabel('Linkage Disequilibrium')
	#plt.xlabel('Distance (nucleotides)')
	#plt.show()
	
       	#synYPos = raw_input('Y Position of top text label --> ')
	#nSynYPos = raw_input('Y Position of bottom text label --> ')
	print ''
	#for value in ldValuesSyn:
	#	print type(value)
	#for value in ldValuesNonSyn:
	#	print type(value)
	#print valuesSyn
	#print valuesNonSyn
	#if '-a' in sys.argv:
		#drawAverageGraph(distanceValuesSyn, ldValuesSyn, distanceValuesNonSyn, ldValuesNonSyn, bacteriaName, ext)
		#print 'Saving data to %s...' % ('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/')
		#storeAverageData(distanceValuesSyn, ldValuesSyn, distanceValuesNonSyn, ldValuesNonSyn, bacteriaName)
		#print ' Done.'
	if distanceValuesSyn and ldValuesSyn or distanceValuesNonSyn and ldValuesNonSyn:
		print 'Drawing average graph...'
		
		drawAverageGraph(distanceValuesSyn, ldValuesSyn, distanceValuesNonSyn, ldValuesNonSyn, bacteriaName, ext)
		print 'Done.'
		print 'Drawing normal graph...'
		n = pynotify.Notification("%s - Drawing normal graph" % bacteriaName)
		n.show()
		drawGraph(distanceValuesSyn, ldValuesSyn, distanceValuesNonSyn, ldValuesNonSyn, bacteriaName, ext)
		print 'Done.'
		print 'Saving data to %s...' % ('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/')
		storeData(distanceValuesSyn, ldValuesSyn, distanceValuesNonSyn, ldValuesNonSyn, bacteriaName)
		print ' Done.'
	else:
		print 'No polymorphisms found, not writing a graph.'
	print ' Done.'
	#print 'Uploading to dropbox...',
	#pdfPath = '"' + dir[:-1] + '.' + 'pdf' + '"'
	#pngPath = '"' + dir[:-1] + '.' + 'png' + '"'
	#dropboxUpload(pdfPath, '"/Biomedicine/Yr 3/FYP/Results for Adam/"')
	#dropboxUpload(pngPath, '"/Biomedicine/Yr 3/FYP/Results for Adam/"')
	#print ' Done.'
	#print 'Saving to usb...',
	#os.sys('cp /home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '.pdf /media/robert/ROBUSB/')
	#print ' Done.'
	return ((distanceValuesSyn, ldValuesSyn), (distanceValuesNonSyn, ldValuesNonSyn))

def storeAverageData(listSynDistance, listSynLd, listNonSynDistance, listNonSynLd, bacteriaName):
	f = open('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/averageSynDistance.txt', 'w')
	f.write(str(listSynDistance))
	f.close()
	f = open('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/averageSynLd.txt', 'w')
	f.write(str(listSynLd))
	f.close()
	f = open('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/averageNonSynDistance.txt', 'w')
	f.write(str(listNonSynDistance))
	f.close()
	f = open('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/averageNonSynLd.txt', 'w')
	f.write(str(listNonSynLd))
	f.close()
	return 
	
def loadAverageData(bacteriaName):
	f = open('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/averageSynDistance.txt', 'r')
	averageSynDistance = eval(f.read())
	f.close()
	f = open('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/averageSynLd.txt', 'r')
	averageSynLd = eval(f.read())
	f.close()
	f = open('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/averageNonSynDistance.txt', 'r')
	averageNonSynDistance = eval(f.read())
	f.close()
	f = open('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/averageNonSynLd.txt', 'r')
	averageNonSynLd = eval(f.read())
	f.close()
	return (averageSynDistance, averageSynLd, averageNonSynDistance, averageNonSynLd)
	
def drawAverageGraph(distanceValuesSyn, ldValuesSyn, distanceValuesNonSyn, ldValuesNonSyn, bacteriaName, ext):
	'''Takes all the info needed to plot an average graph for the results'''
	listSynDistance = []
	listSynLd = []
	listNonSynDistance = []
	listNonSynLd = []
	print 'Calculating averages...'
	#n = pynotify.Notification("%s - Calculating averages" % bacteriaName)
	#n.show()
	#print valuesSyn
	#print valuesNonSyn
	argumentsList = sys.argv
	if '-a' in argumentsList:
		aIndex = argumentsList.index('-a')
		numIndex = aIndex + 1
		try:
			float(argumentsList[numIndex])
			typeNumIndex = True
		except:
			typeNumIndex = False
		if typeNumIndex == True:
			try:
				rangeNumber = int(argumentsList[numIndex])
			except:
				print '------ERROR------ You must use an integer for the range'
				quit()
		else:
			rangeNumber = 0
	else:
		rangeNumber = 0
	if len(distanceValuesSyn) > 0:
		maxDistanceValuesSyn = max(distanceValuesSyn)
		#ldValuesList = []
		for distanceValue in distanceValuesSyn:
			#sys.stdout.write("%s bacteria to do. Averages: %s => %s - %s / %s \r" % (str('1') ,str(bacteriaName), str('Syn'), str(distanceValue), str(maxDistanceValuesSyn)))
			ldValuesList = []
			#sys.stdout.flush()
			distanceValueAlreadyCalculated = False
			#print "Range for this one is from %s to %s with %s being the dv. Type %s" % (int(distanceValue)-11, int(distanceValue)+11, distanceValue, type(distanceValue))
			rangeLimit = rangeNumber + 1
			for x in range(int(distanceValue)-rangeLimit, int(distanceValue)+rangeLimit):
				if x in listSynDistance:
					distanceValueAlreadyCalculated = True
			if distanceValue not in listSynDistance:
				if distanceValueAlreadyCalculated == False:
					sys.stdout.write("%s bacteria to do. Averages: %s => %s - %s / %s \r" % (str('1'), str(bacteriaName), str('Syn'), str(distanceValue), str(maxDistanceValuesSyn)))
					sys.stdout.flush()
					#print '---Unqiue'
					#If its not already in the list its a novel distanceValue, this means we should do stuff
					#ldValuesList = []
					for testDistanceValue in distanceValuesSyn:
						if distanceValue-rangeNumber <= testDistanceValue <= distanceValue+rangeNumber: #They're the same value so we should add the corresponding LD value to the list
							#print '---Adding an LD value'
							ldIndex = distanceValuesSyn.index(distanceValue)
							ldValuesList.append(float(ldValuesSyn[ldIndex]))
					if len(ldValuesList) > 0:
						#print len(ldValuesList)
						#print ldValuesList
						ldValuesAverage = sum(ldValuesList) / float(len(ldValuesList))
						#print "%s from %s / %s" % (ldValuesAverage, sum(ldValuesList), float(len(ldValuesList)))
						listSynDistance.append(float(distanceValue))
						listSynLd.append(float(ldValuesAverage))
	if len(distanceValuesNonSyn) > 0:
		maxDistanceValuesNonSyn = max(distanceValuesNonSyn)
		#ldValuesList = []
		for distanceValue in distanceValuesNonSyn:
			#sys.stdout.write("Averages: %s => %s - %s / %s \r" % (str('1'), str('Non Syn'), str(distanceValue), str(maxDistanceValuesNonSyn)))
			ldValuesList = []
			#sys.stdout.flush()
			distanceValueAlreadyCalculated = False
			#print "Range for this one is from %s to %s with %s being the dv. Type %s" % (int(distanceValue)-11, int(distanceValue)+11, distanceValue, type(distanceValue))
			rangeLimit = rangeNumber + 1
			for x in range(int(distanceValue)-rangeLimit, int(distanceValue)+rangeLimit):
				if x in listNonSynDistance:
					distanceValueAlreadyCalculated = True
			if distanceValue not in listNonSynDistance:
				if distanceValueAlreadyCalculated == False:
					sys.stdout.write("%s bacteria to do. Averages: %s => %s - %s / %s \r" % (str('1'), str(bacteriaName), str('Non Syn'), str(distanceValue), str(maxDistanceValuesNonSyn)))	
					sys.stdout.flush()			
					#print '---Unqiue'
					#If its not already in the list its a novel distanceValue, this means we should do stuff
					#ldValuesList = []
					for testDistanceValue in distanceValuesNonSyn:
						if distanceValue-rangeNumber <= testDistanceValue <= distanceValue+rangeNumber: #They're the same value so add the corresponding LD value to the list
							#print '---Adding an LD value'
							ldIndex = distanceValuesNonSyn.index(testDistanceValue)
							ldValuesList.append(float(ldValuesNonSyn[ldIndex]))
					if len(ldValuesList) > 0:
						#print len(ldValuesList)
						ldValuesAverage = sum(ldValuesList) / float(len(ldValuesList))
						#print "%s from %s / %s" % (ldValuesAverage, sum(ldValuesList), float(len(ldValuesList)))
						listNonSynDistance.append(float(distanceValue))
						listNonSynLd.append(float(ldValuesAverage))
	#print listSynDistance
	#time.sleep(1)
	#print listSynLd
	#time.sleep(1)
	#print listNonSynDistance
	#time.sleep(1)
	#print listNonSynLd
	#time.sleep(1)
	'''for k,v in valuesSyn.items():
		print sum(v)
		print type(sum(v))
		print len(v)
		ldAverage = float(sum(v)) / float(len(v))
		listSynDistance.append(float(k))
		listSynLd.append(ldAverage)
	for k,v in valuesNonSyn.items():
		print sum(v)
		print type(sum(v))
		print len(v)
		ldAverage = float(sum(v)) / float(len(v))
		listNonSynDistance.append(float(k))
		listNonSynLd.append(ldAverage)
	#print listSynDistance
	#time.sleep(2)
	#print listNonSynDistance
	#time.sleep(2)
	#print listSynLd
	#time.sleep(2)
	#print listNonSynLd
	#time.sleep(2)
	try:
		listSynDistance[50]
	except:
		pass
	else:
		print type(listSynDistance[50])
	'''
	print ' Done.'
	print 'Drawing average graph...'
	#n = pynotify.Notification("%s - Drawing average graph" % bacteriaName)
	#Do both graphs
	#plt.plot(listSynDistance, listSynLd, 'ro')
	#plt.plot(listNonSynDistance, listNonSynLd, 'bo')
	#plt.show()
	#print len(listSynDistance)
	#print len(listSynLd)
	#print len(listNonSynDistance)
	#print len(listNonSynLd)
	plt.plot(listSynDistance, listSynLd, 'ro')
	plt.plot(listNonSynDistance, listNonSynLd, 'bo')
	plt.ylabel('Average Linkage Disequilibrium')
	plt.xlabel('Distance (nucleotides)')
	plotAxes = plt.axis()
	#plotAxes = (plotAxes[0], plotAxes[1], 0.0, 1.0)
	#plt.axis(plotAxes)
	xPos = 0.6*plotAxes[1]
	synYPos = 0.9*plotAxes[3]
	nSynYPos = 0.8*plotAxes[3]
	plt.text(xPos, synYPos, 'Synonymous', color='r')
	plt.text(xPos, nSynYPos, 'Non Synonymous', color='b')
	plt.title(bacteriaName)
	plt.suptitle('A graph of Average Linkage Disequilibrium Against Loci Distance')
	print ' Done.'
	print 'Writing graph file to %s%s%s...' % ('/home/robert/FYP/Sequence Data/goodBacteria/Average ', bacteriaName, '.' + ext)
	#plt.show()
	#plt.savefig('/home/robert/FYP/Sequence Data/goodBacteria/Average ' + bacteriaName + '.' + 'pdf')
	if rangeNumber > 0:
		plt.savefig('/home/robert/FYP/Sequence Data/goodBacteria/Average ' + str(rangeNumber) + 'bp ' + bacteriaName + '.' + 'png')
	elif rangeNumber == 0:
		plt.savefig('/home/robert/FYP/Sequence Data/goodBacteria/Average ' + bacteriaName + '.' + 'png')
		print 'Storing average data...'
		storeAverageData(listSynDistance, listSynLd, listNonSynDistance, listNonSynLd, bacteriaName)
		print ' Done.'
	#plt.show()
	plt.close()
	'''if listSynDistance and listSynLd and len(listSynDistance) == len(listSynLd):
		#Do the Synonymous graph
		plt.figure(2)
		#plt.subplot(211)
		print len(listSynDistance)
		print len(listSynLd)
		plt.plot(listSynDistance, listSynLd, 'ro')
		plt.ylabel('Average Linkage Disequilibrium')
		plt.xlabel('Distance (nucleotides)')
		plotAxes = plt.axis()
		#plotAxes = (plotAxes[0], plotAxes[1], 0.0, 1.0)
		#plt.axis(plotAxes)
		xPos = 0.6*plotAxes[1]
		synYPos = 0.9*plotAxes[3]
		nSynYPos = 0.8*plotAxes[3]
		plt.text(xPos, synYPos, 'Synonymous', color='r')
		plt.title(bacteriaName)
		plt.suptitle('A graph of Average Linkage Disequilibrium Against Loci Distance')
		
		plt.savefig('/home/robert/FYP/Sequence Data/goodBacteria/Average Synonymous' + bacteriaName + '.' + 'png')
		print ' Done.'
	#print 'Writing graph file to %s%s%s...' % ('/home/robert/FYP/Sequence Data/goodBacteria/Average ', bacteriaName, '.' + ext),
	#plt.show()
	#plt.close()
	if listNonSynDistance and listNonSynLd and len(listNonSynDistance) == len(listNonSynLd):
		#Do the Non Synonymous graph
		plt.figure(3)
		plt.plot(listNonSynDistance, listSynLd, 'bo')
		plt.ylabel('Average Linkage Disequilibrium')
		plt.xlabel('Distance (nucleotides)')
		plotAxes = plt.axis()
		#plotAxes = (plotAxes[0], plotAxes[1], 0.0, 1.0)
		#plt.axis(plotAxes)
		xPos = 0.6*plotAxes[1]
		synYPos = 0.9*plotAxes[3]
		nSynYPos = 0.8*plotAxes[3]
		plt.text(xPos, synYPos, 'Non Synonymous', color='b')
		plt.title(bacteriaName)
		plt.suptitle('A graph of Average Linkage Disequilibrium Against Loci Distance')
		
		plt.savefig('/home/robert/FYP/Sequence Data/goodBacteria/Average Non Synonymous' + bacteriaName + '.' + 'png')
		print ' Done.'''
	#print 'Writing graph file to %s%s%s...' % ('/home/robert/FYP/Sequence Data/goodBacteria/Average ', bacteriaName, '.' + ext),
	#plt.savefig('/home/robert/FYP/Sequence Data/goodBacteria/Average Non Synonymous' + bacteriaName + '.' + 'pdf')
	#plt.savefig('/home/robert/FYP/Sequence Data/goodBacteria/Average Non Synonymous' + bacteriaName + '.' + 'png')
	
	print '----------------Done.'



def drawGraph(distanceValuesSyn, ldValuesSyn, distanceValuesNonSyn, ldValuesNonSyn, bacteriaName, ext):
	'''Takes all the info required to draw a graph for the results, draws it and saves it to the specified path'''
	#print distanceValuesSyn
	#print distanceValuesNonSyn
	print 'Drawing graph...',
	plt.plot(distanceValuesSyn, ldValuesSyn, 'ro')
	plt.plot(distanceValuesNonSyn, ldValuesNonSyn, 'bo')
	plt.ylabel('Linkage Disequilibrium')
	plt.xlabel('Distance (nucleotides)')
	#chooseChange = raw_input('Change the scale? (y/n) --> ')
       	#if chooseChange == 'y':
	#	axesScale = raw_input('Whats the borders then? [xmin, xmax, ymin, ymax] --> ')
	#	plt.axis(axesScale)
	plotAxes = plt.axis()
	#plotAxes[2] = 0.0
	#plotAxes[3] = 1.0
	#print_r(plotAxes)
	plotAxes = (plotAxes[0], plotAxes[1], 0.0, 1.0)
	plt.axis(plotAxes)
	#print '--------------------------------------------- ' + str(plotAxes)
	#time.sleep(3)
	xPos = 0.6*plotAxes[1]
	synYPos = 0.9*plotAxes[3]
	nSynYPos = 0.8*plotAxes[3]
	plt.text(xPos, synYPos, 'Synonymous', color='r')
	plt.text(xPos, nSynYPos, 'Non Synonymous', color='b')
	plt.title(bacteriaName)
	plt.suptitle('A graph of Linkage Disequilibrium Against Loci Distance')
	print ''
	print 'Writing graph file to %s%s%s...' % ('/home/robert/FYP/Sequence Data/goodBacteria/', bacteriaName, '.' + ext),
	plt.savefig('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '.' + 'pdf')
	plt.savefig('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '.' + 'png')
	print ' Done.'
	
	return '/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '.' + ext

def storeData(distanceValuesSyn, ldValuesSyn, distanceValuesNonSyn, ldValuesNonSyn, bacteriaName):
	'''Takes all the info that needs to be stored plus the bacteria name and writes the 
	lists to individual files in the bacteriaName folder'''
	f = open('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/distanceValuesSyn.txt', 'w')
	f.write(str(distanceValuesSyn))
	f.close()
	f = open('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/ldValuesSyn.txt', 'w')
	f.write(str(ldValuesSyn))
	f.close()
	f = open('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/distanceValuesNonSyn.txt', 'w')
	f.write(str(distanceValuesNonSyn))
	f.close()
	f = open('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/ldValuesNonSyn.txt', 'w')
	f.write(str(ldValuesNonSyn))
	f.close()
	return 
	
		
def chooseSequenceData():
	'''Prompts the user to select a bacterium from a list of folder in a location on the
	hard drive, then returns a path to the folder containing the bacterium files'''
	if os.name == 'nt': #On Windows
		dirPrefix = 'F:\\USER FILES\\Dropbox\\Dropbox\\Biomedicine\\Yr 3\\FYP\\Sequence Data\\' #Choose this folder to get the sequence data from
	elif os.name == 'posix': #On Mac/Linux
		dirPrefix = '/Users/robert/Dropbox/Biomedicine/Yr 3/FYP/Sequence Data/' #Choose this folder to get all the sequence data from
	folderList = os.listdir(dirPrefix) #Gets all the file names in the prefix dir
	bacteriaList = [] #Initialises a list for putting the bacteria names in
	for folder in folderList: 
		if folder[0] != '.':#Folder is not a hidden system folder
			bacteriaList.append(folder) #Append to the list of bacteria
	for bacteria in bacteriaList:
		print '%s. %s' % (bacteriaList.index(bacteria), bacteria) #Print out a list of the bacterium for the user to choose from
	print 'Which bacterium would you like to load?'
	while True: #This section is to catch user input and ensure it is a number and in range, if not it asks again
		chosenBacteriaIndex = raw_input('---> ') #Get user input
		if chosenBacteriaIndex.isalpha() == False: #Only if its a number
			chosenBacteriaIndex = int(chosenBacteriaIndex) #str -> int
			if chosenBacteriaIndex <= len(bacteriaList): #If its in range
				break #Input is all good, break to next part of program
		print 'Try again...'#Only printed if input is not good
	chosenBacteria = bacteriaList[chosenBacteriaIndex] #Set chosenBacteria to the name of the bacteria chosen by the user
	chosenBacteriaDir = dirPrefix + chosenBacteria + '/'
	return chosenBacteriaDir

def dropboxUpload(pdfPath, uploadPath):
	'''Takes a path to a pdf file and a path to upload it to on Dropbox then runs dropbox_uploader.sh from the ~/ directory'''
	os.system('~/dropbox_uploader.sh upload ' + pdfPath + ' ' + uploadPath)
	
	
def main():
	#chosenBacteriaDir = chooseSequenceData()
	#print chosenBacteriaDir
	#if os.name == 'nt': #If on windows, open this tester file
		#fileObject = openFile('F:\USER FILES\Dropbox\Dropbox\Biomedicine\Yr 3\FYP\Bifidobacterium animalis lactis\ortholog_000397.nt_ali.Bifidobacteriumanimalissubsplactis.fasta')
	#elif os.name == 'posix': #If on mac, open this tester file
	#fileObject = openFile('/Users/robert/Dropbox/Biomedicine/Yr 3/FYP/Sequence Data/Chlamydia psittaci/ortholog_000004.nt_ali.fasta')
	#fileObject = openFile('/home/robert/FYPsequences/Propionibacterium acnes/ortholog_000553.nt_ali.fasta')
	#sequences = getSequences(fileObject)
	#print sequences
	#for sequence in sequences:
	#	print len(sequence)
	#allSeqCodons = getAllCodons(sequences)
	#allPolymorphisms = categorisePolymorphisms(allSeqCodons)
	#synPolymorphisms = allPolymorphisms[0]
	#nonSynPolymorphisms = allPolymorphisms[1]
	#calcLinkageDisequilibrium(allPolymorphisms, sequences, file)
	#variation = getVariation(sequences)
	#print variation
	#for variationIndex in variation:
	#	currentIndexPoly = []
	#	for sequence in sequences:
	#		if sequence[variationIndex] not in currentIndexPoly:
	#			currentIndexPoly.append(sequence[variationIndex])
#	#	print '%s: %s' % (variationIndex, str(currentIndexPoly))
		
	#allCodons = getAllCodons(sequences)
	#print allCodons
	#allAminoAcids = getAllAminoAcids(allCodons)
	#print allAminoAcids
	#nonSynPoly = getNonSynonymousPolymorphisms(allAminoAcids)
	#if nonSynPoly:
	#print nonSynPoly	
	#print '-------------------------'
	#synPoly = getSynonymousPolymorphisms(allCodons)
	#print synPoly
	#print getNonSynonymousPolymorphisms(allSeqCodons)
	#print fileObject
	#print getVariation(sequences)
	#print seqNuc
	#print seqHeader
	#storeData(['distanceValuesSyn'], ['ldValuesSyn'], ['distanceValuesNonSyn'], ['ldValuesNonSyn'], 'Chlamydia psittaci')
	#time.sleep(10)
	pynotify.init("icon-summary-body")
	if len(sys.argv) == 1:
		print '--------------------------------------------------'
		print 'This program takes specified fasta files, a directory containing fasta files or a directory containing several folders which all contain fasta files'
		print 'It then calculates the Linkage Disequilibrium for every polymorphic site found within the files and outputs a graph of Linkage Disequilibrium'
		print 'against loci distance.'
		print 'Usage:'
		print 'ldcalc.py [Path To Folder Containing Single Bacteria Files] => Calculate the LD for a single bacteria'
		print 'ldcalc.py -e [Extension] => Specifiy the extension of the graph file to create. Depending on the modules you have installed the possible extensions are pdf, png and jpg'
		print 'ldcalc.py -s => Shutdown the script after running (Requires sudo)'
		print 'ldcalc.py -r [Path to directory containing folders containing fasta files] => Recursive mode, iterates through bacteria creating a graph for each one'

	if '-s' in sys.argv:
		if os.geteuid() != 0:
			exit("You need to have root privileges to run this script and shutdown after.\nPlease try again, this time using 'sudo'.")
	if '-r' in sys.argv:
		if '-e' in sys.argv:
			if sys.argv[sys.argv.index('-e') + 1]:
				ext = sys.argv[sys.argv.index('-e') + 1]
				print 'Using the %s extension to save the files' % ext
				iterateFolders('/home/robert/FYP/Sequence Data/goodBacteria/', ext)
			else:
				print 'If you supply the -e flag you must specify the extension to write the graph files to after'
				quit()
		else:
			iterateFolders('/home/robert/FYP/Sequence Data/goodBacteria/', 'png')
	else:
		if '-' not in sys.argv[-1]:
			bacteriaName = sys.argv[-1]
		else:
			bacteriaName = 'Alteromonas macleodii (All Strains)'
		dataPath = '/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/'
		bacteriaFiles = os.listdir(dataPath)
		if 'distanceValuesNonSyn.txt' in bacteriaFiles and 'distanceValuesSyn.txt' in bacteriaFiles and 'ldValuesNonSyn.txt' in bacteriaFiles and 'ldValuesSyn.txt' in bacteriaFiles and '-forcecalc' not in sys.argv:
			print 'LD and distances already calculated, loading now...'
			f = open(dataPath + 'distanceValuesNonSyn.txt')
			distanceValuesNonSyn = eval(f.read())
			f.close()
			f = open(dataPath + 'distanceValuesSyn.txt')
			distanceValuesSyn = eval(f.read())
			f.close()
			f = open(dataPath + 'ldValuesNonSyn.txt')
			ldValuesNonSyn = eval(f.read())
			f.close()
			f = open(dataPath + 'ldValuesSyn.txt')
			ldValuesSyn = eval(f.read())
			f.close()
			print ' Done.'
			drawAverageGraph(distanceValuesSyn, ldValuesSyn, distanceValuesNonSyn, ldValuesNonSyn, bacteriaName, 'png')
		else:
			iterateFiles('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/', bacteriaName, 'png')
	if '-s' in sys.argv:
		os.system('sudo shutdown -h 60')
		os.system('sudo shutdown -h 60')
	#if sys.argv[1]:
	#	bacteriaName = sys.argv[1]
	#else:
	#	bacteriaName = raw_input('Bacteria Name: ')
	#iterateFiles('/home/robert/FYP/Sequence Data/goodBacteria/' + bacteriaName + '/', bacteriaName) #Calculates for all the files in the chosen folder
	#ld = newCalcLinkageDisequilibrium(variation, sequences)
	#print sequences0
	#distanceValuesSyn = ld[0][0]
	#ldValuesSyn = ld[0][1]
	#distanceValuesNonSyn = ld[1][0]
	#ldValuesNonSyn = ld[1][1]
	
	#plt.acorr( #Add some correlation somehow!)
	#plt.show()
	
	
if __name__ == '__main__':
	main()
	
