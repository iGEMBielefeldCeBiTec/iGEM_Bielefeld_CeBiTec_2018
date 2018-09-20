#!/usr/b-in/env python
# -*- coding: utf-8 -*-
#

''' siRNA Designer
    Copyright (C) 2018  Vanessa Krämer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
	
	
'''


import math
import sys
from operator import itemgetter
import re
from sets import Set

def read_genesequence(fasta_file):
	
	geneseqdic = {}
	
	head = ""
	geneseq = ""
	
	try:
	
		with open (fasta_file) as f:
		
			line  = f.readline()
			
			if line[0] != '>':

				sys.exit("The Input has the wrong format")
		
			while line:
		
				if line[0] == '>':
					if not geneseqdic:
						
						head = line.strip()
						
					else:
						
						geneseqdic.update({head:geneseq})
						
						head = line.strip()
			
			
				else:
					geneseq = geneseq + line.strip()
			

				line = f.readline()
		
		geneseqdic.update({head:geneseq})
		
	except (IOError, OSError) as e:
		
		if e.errno==EPERM or e.errno==EACCES:
			sys.exit("You dont have permission to open the file.")
		elif e.errno==ENOENT:
			sys.exit("The input file does not exist.")
		elif IOError:
			sys.exit("Unexpected error occured.")
		elif OSError:
			sys.exit("Unexpected error occured.")
		
	return geneseqdic


def revcomp( seq ):
	"""! @brief construct reverse complement of sequence """
	
	new_seq = []
	
	bases = { 'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n', '-':'-' }
	for nt in seq.lower():
		try:
			new_seq.append( bases[nt] )
		except:
			new_seq.append( 'n' )
	return ''.join( new_seq[::-1] ).upper()


def comp( seq ):

	new_seq = []
	
	bases = { 'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n', '-':'-' }
	for nt in seq.lower():
		try:
			new_seq.append( bases[nt] )
		except:
			new_seq.append( 'n' )
	return ''.join(new_seq).upper()
	

#TODO Ns beachten???
def calculate_GC_content(sequence):
	
	at_count = 0
	
	gc_count = 0
	
	for nuc in sequence:
		
		if nuc == 'A' or nuc == 'T':
			
			at_count +=1
		
		else:
			
			gc_count +=1
	
	gc_content = gc_count/len(sequence)
	
	return gc_content
	

#standard 50mM Na+
def calculate_Tm(gc_content, seq_len):
	
	tm = 79.8 + (18.5* math.log10(0.05)) + (58.4 * gc_content)+(11.8* math.pow((gc_content),2))-(820/seq_len)
	
	return tm


def number_of_AT_region(sequence):
	
	score = 0
	
	i=14
	
	while i < 19:
		
		if sequence[i] == 'A' or sequence[i] == 'T':
			score += 1 
		
		i+=1
	
	return score


def number_of_AT_region_Ui_Tei(sequence):
	
	score = 0
	
	i=12
	
	while i < 19:
		
		if sequence[i] == 'A' or sequence[i] == 'T':
			score += 1 
		
		i+=1
	
	return score


def find_GC_stretch(sequence):
	
	if "GCGCGCGCGC" in sequence:
		
		return 0
	
	else:
		
		return 1


def test_rational_siRNA(sequence):

	score = 0
	
	GC = calculate_GC_content(sequence)
	
	TM = calculate_Tm(GC, len(sequence))
	
	AT_score = number_of_AT_region(sequence)
	
	if GC >= 0.3 and GC < 0.52:
		score +=1
	#print score
	if AT_score >=3:
		score = score + AT_score
	#print score
	if TM < 20:
		score +=1
	#print score	
	if sequence[18] == 'A':
		score +=1
	#print score	
	if sequence[2] == 'A':
		score +=1
	#print score
	if sequence[9] == 'T':
		score +=1
	#print score
	if sequence[18] == 'A' or sequence[18] == 'T':
		score -=1
	#print score
	if sequence[12] == 'A' or sequence[12] == 'T' or sequence[12] == 'C':
		score -=1
	
	#print score
	return score
	

def test_Ui_Tei_siRNA(sequence):
	
	score = 0
	
	AT_score = number_of_AT_region_Ui_Tei(sequence)
	GC_score = find_GC_stretch(sequence)
	
	#print "AT-score " + str(AT_score)
	
	if AT_score >=5:
		score += 1
	
	#print score
		
	if sequence[18] == 'A' or sequence[18] == 'T':
		score +=1
		
	#print score
	
	if sequence[0] == 'G' or sequence[0] == 'C':
		score +=1
	
	#print score
	
	#print "GC-score" + str(GC_score)
	if GC_score == 1:
		score +=1
	
	#print score
	return score
	

def get_rational_siRNA(gene_sequence):
	
	all_possible_siRNAs = {}
	
	overhang_possible_siRNAs = {}
	
	k = 19
	
	#make k-mers and check siRNA
	
	#offset of 50 nucleotides at start and end
	i = 50
	j = len(gene_sequence)-50-k
	
	while i < j:
		
		pos = i+k
		
		k_mer = revcomp(gene_sequence[i:pos])
		
		
		if "AAAA" in k_mer or "TTTT" in k_mer or "GGGG" in k_mer or "CCCC" in k_mer:
			pass
		
		else:
			
			k_mer_score = test_rational_siRNA(k_mer)
		
			#print k_mer_score
		
			if k_mer_score > 6:
			
				all_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k))})
				
			
				if gene_sequence[pos+2] != 'C' and gene_sequence[pos+3] != 'T':
					
					overlap = "A"+ "G"  + gene_sequence[pos+1]+ gene_sequence[pos]
					
					print overlap
					
					overhang_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k), overlap)})
		
		i +=1
		
	return overhang_possible_siRNAs


def get_rational_siRNA_silincing(gene_sequence):
	
	all_possible_siRNAs = {}
	
	overhang_possible_siRNAs = {}
	
	k = 19
	
	#make k-mers and check siRNA
	
	#offset of 50 nucleotides at start and end
	i = 50
	j = len(gene_sequence)-50-k
	
	while i < j:
		
		pos = i+k
		
		k_mer = revcomp(gene_sequence[i:pos])
		
		
		
		
		if "AAAA" in k_mer or "TTTT" in k_mer or "GGGG" in k_mer or "CCCC" in k_mer:
			pass
		
		else:
			
			k_mer_score = test_rational_siRNA(k_mer)
		
			#print k_mer_score
		
			if k_mer_score > 6:
			
				all_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k))})
					
				overhang_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k))})
		
		i +=1
		
	return overhang_possible_siRNAs
	
	
def get_rational_siRNA_short_sequence(gene_sequence):
	
	all_possible_siRNAs = {}
	
	overhang_possible_siRNAs = {}
	
	k = 19
	
	i = 0
	j = len(gene_sequence)-4
	
	while i < j-19:
		
		pos = i+k
		
		k_mer = revcomp(gene_sequence[i:pos])
		
		if "AAAA" in k_mer or "TTTT" in k_mer or "GGGG" in k_mer or "CCCC" in k_mer:
			pass
		
		else:
			k_mer_score = test_rational_siRNA(k_mer)
		
			#print k_mer_score
		
			if k_mer_score > 6:
			
				all_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k))})
			
				if gene_sequence[pos+3] != 'C' and gene_sequence[pos+4] != 'T':
					
					overlap = gene_sequence[pos+1] + gene_sequence[pos+2]+ "G" + "A"
					
					overhang_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k), overlap)})
		
		i +=1
		
	
	
	return overhang_possible_siRNAs


def get_rational_siRNA_short_sequence_silincing(gene_sequence):
	
	all_possible_siRNAs = {}
	
	overhang_possible_siRNAs = {}
	
	k = 19
	
	i = 0
	j = len(gene_sequence)-4
	
	while i < j-19:
		
		pos = i+k
		
		k_mer = revcomp(gene_sequence[i:pos])
		
		if "AAAA" in k_mer or "TTTT" in k_mer or "GGGG" in k_mer or "CCCC" in k_mer:
			pass
		
		else:
			k_mer_score = test_rational_siRNA(k_mer)
		
			#print k_mer_score
		
			if k_mer_score > 6:
			
				all_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k))})
			
					
				overhang_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k))})
		
		i +=1
		
	
	
	return overhang_possible_siRNAs
	

def get_Ui_Tei_siRNA(gene_sequence):
	
	all_possible_siRNAs = {}
	
	overhang_possible_siRNAs = {}
	
	
	k = 19
	
	#make k-mers and check siRNA
	
	#offset of 50 nucleotides at start and end
	i = 50
	j = len(gene_sequence)-50-k
	
	while i < j:
		
		pos = i+k
		
		k_mer = revcomp(gene_sequence[i:pos])
		if "AAAA" in k_mer or "TTTT" in k_mer or "GGGG" in k_mer or "CCCC" in k_mer:
			pass
		else:
			k_mer_score = test_Ui_Tei_siRNA(k_mer)
		
			#print k_mer_score
		
			if k_mer_score == 4:
				
				all_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k))})
			
				if gene_sequence[pos+2] != 'C' and gene_sequence[pos+3] != 'T':
					
					overlap = "A"+ "G"  + gene_sequence[pos+1]+ gene_sequence[pos]
					
					overhang_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k), overlap)})
		
		i +=1

	return overhang_possible_siRNAs


def get_Ui_Tei_siRNA_silincing(gene_sequence):
	
	all_possible_siRNAs = {}
	
	overhang_possible_siRNAs = {}
	
	
	k = 19
	
	#make k-mers and check siRNA
	
	#offset of 50 nucleotides at start and end
	i = 50
	j = len(gene_sequence)-50-k
	
	while i < j:
		
		pos = i+k
		
		k_mer = revcomp(gene_sequence[i:pos])
		if "AAAA" in k_mer or "TTTT" in k_mer or "GGGG" in k_mer or "CCCC" in k_mer:
			pass
		else:
			k_mer_score = test_Ui_Tei_siRNA(k_mer)
		
			#print k_mer_score
		
			if k_mer_score == 4:
				
				all_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k))})
					
				overhang_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k))})
		
		i +=1

	return overhang_possible_siRNAs


def get_Ui_Tei_siRNA_short_sequence(gene_sequence):
	
	all_possible_siRNAs = {}
	
	overhang_possible_siRNAs = {}
	
	k = 19
	
	i = 0
	j = len(gene_sequence)-4
	
	while i < j-19:
		
		pos = i+k
		
		k_mer = revcomp(gene_sequence[i:pos])
		if "AAAA" in k_mer or "TTTT" in k_mer or "GGGG" in k_mer or "CCCC" in k_mer:
			pass
		else:
			
			k_mer_score = test_Ui_Tei_siRNA(k_mer)
		
			#print k_mer_score
		
			if k_mer_score == 4:
			
				all_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k))})
			
				if gene_sequence[pos+2] != 'C' and gene_sequence[pos+3] != 'T':
					
					overlap = "A" + "G" + revcomp(gene_sequence[pos+1] + gene_sequence[pos+2]+ "G" + "A")
					
					overhang_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k), overlap)})
		
		i +=1

	return overhang_possible_siRNAs


def get_Ui_Tei_siRNA_short_sequence_silincing(gene_sequence):
	
	all_possible_siRNAs = {}
	
	overhang_possible_siRNAs = {}
	
	k = 19
	
	i = 0
	j = len(gene_sequence)-4
	
	while i < j-19:
		
		pos = i+k
		
		k_mer = revcomp(gene_sequence[i:pos])
		if "AAAA" in k_mer or "TTTT" in k_mer or "GGGG" in k_mer or "CCCC" in k_mer:
			pass
		else:
			
			k_mer_score = test_Ui_Tei_siRNA(k_mer)
		
			#print k_mer_score
		
			if k_mer_score == 4:
			
				all_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k))})
			
				overhang_possible_siRNAs.update({k_mer:(k_mer_score,str(i),str(i+k))})
		
		i +=1

	return overhang_possible_siRNAs


def get_probability(siRNA):
	
	P_eff_90 = 0.1
	P_eff_80 = 0.2
	
	P_inf_90 = 1-P_eff_90
	P_inf_80 = 1-P_eff_80
	
	independent_effective = {"A" : [0.12,0.3,0.318,0.218,0.271,0.294,0.247,0.298,0.271,0.229,0.25,0.283,0.259,0.282,0.27,0.232,0.288,0.313,0.312],"G" : [0.557,0.27,0.229,0.291,0.262,0.239,0.304,0.253,0.224,0.257,0.279,0.247,0.255,0.276,0.242,0.287,0.24,0.239,0.208],"C" : [0.208,0.242,0.208,0.294,0.248,0.208,0.263,0.229,0.244,0.275,0.256,0.234,0.247,0.217,0.202,0.251,0.235,0.178,0.233],"T" : [0.115,0.187,0.245,0.197,0.218,0.259,0.186,0.22,0.261,0.239,0.216,0.235,0.239,0.224,0.286,0.23,0.236,0.27,0.248], "N" : [0.115, 0.187, 0.208, 0.197, 0.218, 0.208, 0.186, 0.22, 0.224, 0.229, 0.216, 0.234, 0.239, 0.217, 0.202, 0.23, 0.235, 0.178, 0.208]}
	
	dependent_nucleotides = {"AA" : [0.16,0.328,0.2,0.247,0.292,0.229,0.282,0.278,0.204,0.267,0.274,0.195,0.245,0.264,0.218,0.249,0.3,0.3], "AG" : [0.35,0.272,0.336,0.341,0.265,0.363,0.311,0.286,0.323,0.356,0.284,0.314,0.329,0.306,0.329,0.275,0.288,0.235], "AC" : [0.28,0.192,0.264,0.22,0.212,0.257,0.214,0.254,0.257,0.257,0.25,0.263,0.236,0.179,0.258,0.218,0.154,0.231], "AT" : [0.21,0.208,0.2,0.192,0.23,0.151,0.194,0.181,0.217,0.12,0.192,0.229,0.19,0.251,0.196,0.259,0.254,0.235 ], "AN" : [0.16, 0.192, 0.2, 0.192, 0.212, 0.151, 0.194, 0.181, 0.204, 0.12, 0.192, 0.195, 0.19, 0.179, 0.196, 0.218, 0.154, 0.231], "GA" : [0.33,0.364,0.241,0.314,0.394,0.342,0.316,0.313,0.337,0.271,0.358,0.364,0.373,0.33,0.267,0.285,0.37,0.387], "GG" : [0.284,0.169,0.293,0.194,0.206,0.241,0.241,0.218,0.193,0.271,0.263,0.204,0.274,0.226,0.231,0.238,0.2,0.166], "GC" : [0.213,0.227,0.251,0.26,0.22,0.246,0.213,0.223,0.203,0.285,0.19,0.204,0.193,0.165,0.249,0.247,0.185,0.196], "GT" : [0.172,0.24,0.215,0.231,0.179,0.171,0.229,0.246,0.267,0.173,0.19,0.228,0.16,0.278,0.151,0.23,0.245,0.236], "GN" : [0.172, 0.169, 0.215, 0.194, 0.179, 0.171, 0.213, 0.218, 0.193, 0.173, 0.19, 0.204, 0.16, 0.165, 0.151, 0.23, 0.185, 0.166], "CA" : [0.358,0.371,0.277,0.331,0.329,0.295,0.37,0.33,0.261,0.306,0.305,0.282,0.325,0.331,0.222,0.421,0.352,0.412], "CG" : [0.116,0.149,0.127,0.171,0.13,0.185,0.137,0.131,0.192,0.131,0.155,0.169,0.146,0.16,0.107,0.148,0.133,0.108], "CC" : [0.295,0.198,0.318,0.245,0.203,0.243,0.224,0.199,0.271,0.262,0.239,0.241,0.214,0.204,0.187,0.167,0.204,0.23], "CT" : [0.231,0.282,0.277,0.253,0.338,0.277,0.269,0.34,0.276,0.301,0.3,0.308,0.316,0.304,231,0.263,0.306,0.23], "CN" : [0.116, 0.149, 0.127, 0.171, 0.13, 0.185, 0.137, 0.131, 0.192, 0.131, 0.155, 0.169, 0.146, 0.16, 0.107, 0.148, 0.133, 0.108], "TA" : [0.198,0.167,0.172,0.146,0.137,0.144,0.187,0.153,0.134,0.146,0.172,0.204,0.181,0.144,0.151,0.188,0.228,0.183], "TG" : [0.396,0.353,0.368,0.409,0.368,0.389,0.361,0.246,0.304,0.382,0.294,0.321,0.357,0.262,0.396,0.307,0.325,0.272], "TC" : [0.25,0.218,0.353,0.268,0.192,0.301,0.284,0.301,0.359,0.216,0.267,0.281,0.226,0.273,0.236,0.313,0.173,0.263], "TT" : [0.156,0.263,0.108,0.177,0.302,0.167,0.168,0.301,0.203,0.256,0.267,0.194,0.236,0.321,0.276,0.193,0.274,0.277], "TN" : [0.156, 0.167, 0.108, 0.146, 0.137, 0.144, 0.168, 0.153, 0.134, 0.146, 0.172, 0.194, 0.181, 0.144, 0.151, 0.188, 0.173, 0.183]}

	independent_ineffective = {"A" : [0.312,0.247,0.211,0.247,0.254,0.231,0.273,0.26,0.235,0.295,0.251,0.243,0.226,0.235,0.203,0.261,0.262,0.182,0.084], "G" : [0.185,0.229,0.256,0.262,0.237,0.259,0.236,0.254,0.286,0.279,0.231,0.257,0.323,0.266,0.293,0.26,0.255,0.301,0.319], "C" : [0.215,0.253,0.296,0.285,0.266,0.321,0.283,0.253,0.298,0.242,0.269,0.256,0.247,0.244,0.289,0.269,0.262,0.319,0.426], "T" : [0.288,0.272,0.236,0.207,0.243,0.189,0.208,0.234,0.182,0.184,0.248,0.43,0.204,0.255,0.215,0.21,0.221,0.198,0.171], "N" : [0.185, 0.229, 0.211, 0.207, 0.237, 0.189, 0.208, 0.234, 0.182, 0.184, 0.231, 0.243, 0.204, 0.235, 0.203, 0.21, 0.221, 0.182, 0.084]}
		
	deductive_ineffective = {"A" : [0.293,0.233,0.227,0.261,0.243,0.235,0.251,0.234,0.243,0.257,0.25,0.239,0.247,0.239,0.243,0.256,0.237,0.229,0.229],"G" : [0.148,0.243,0.257,0.236,0.246,0.254,0.232,0.249,0.259,0.248,0.24,0.251,0.248,0.241,0.253,0.238,0.253,0.254,0.264],"C" : [0.264,0.253,0.264,0.235,0.251,0.264,0.246,0.257,0.252,0.242,0.248,0.255,0.251,0.261,0.266,0.25,0.255,0.274,0.256],"T" : [0.295,0.271,0.252,0.268,0.261,0.247,0.271,0.26,0.246,0.254,0.261,0.255,0.254,0.259,0.238,0.257,0.255,0.243,0.251], "N" : [0.148, 0.233, 0.227, 0.235, 0.243, 0.235, 0.232, 0.234, 0.243, 0.242, 0.24, 0.239, 0.247, 0.239, 0.238, 0.238, 0.237, 0.229, 0.229]}
	
	produkt_eff = calculate_produkt_dependent(siRNA, dependent_nucleotides, independent_effective)
	produkt_ineff = calculate_produkt_independent(siRNA, independent_ineffective)
	
	probability_90 = (P_eff_90 * produkt_eff)/((P_eff_90 * produkt_eff)+(P_inf_90*produkt_ineff))
	
	probability_80 = (P_eff_80 * produkt_eff)/((P_eff_80 * produkt_eff)+(P_inf_80*produkt_ineff))
	
	return (probability_80, probability_90)
	
	
def calculate_produkt_dependent(si_seq, eff_probs, independent_effective):
	
	produkt = 0 
	
	i = 0
	
	while i < len(si_seq):
		
		if i == 0:
			
			nuc = si_seq[i]
			produkt = produkt + independent_effective[nuc][i]
			
		else:
			
			nuc = si_seq[i]
			nuc_before = si_seq[i-1]
			
			dic_code = nuc + nuc_before
			
			produkt = produkt * eff_probs[dic_code][i-1]
			
		i += 1
	
	return produkt


def calculate_produkt_independent(si_seq, independent):

	produkt = 0
	
	i = 0
		
	while i < len(si_seq):
		
		if i == 0:
			
			nuc = si_seq[i]
			produkt = produkt + independent[nuc][i]
			
		else:
			
			nuc = si_seq[i]
			produkt = produkt * independent[nuc][i]
		
		i += 1
	
	return produkt		
	
	
def get_siRNAs_probabilities(siRNA_dic):

		all_siRNAs_probs = []
		
		for siRNA in siRNA_dic:
			
			prob80, prob90 = get_probability(siRNA)
			
			overlap = siRNA_dic[siRNA][3]
			
			reverse_overlap = overlap[::-1]
			
			final_siRNA1 = overlap + revcomp(siRNA)
			
			final_siRNA2 = revcomp(final_siRNA1)
			
			all_siRNAs_probs.append((round(prob80, 2), final_siRNA1.upper(), final_siRNA2.upper()))
		
		return all_siRNAs_probs
		

def get_siRNAs_probabilities_bielefeld(siRNA_dic):
	
		MicC_overlap =  "GAAA"
		Plasmid_overlap = "TAGC"

		all_siRNAs_probs = []
		
		for siRNA in siRNA_dic:
			
			prob80, prob90 = get_probability(siRNA)
			
			overlap = siRNA_dic[siRNA][3]
			
			final_siRNA1 = overlap + siRNA
			
			final_siRNA2 = revcomp(final_siRNA1)
			
			final_siRNA1 = Plasmid_overlap + final_siRNA1
			
			final_siRNA2 = MicC_overlap + final_siRNA2
			
			all_siRNAs_probs.append((round(prob80, 2), final_siRNA1.upper(), final_siRNA2.upper()))
		
		return all_siRNAs_probs


def get_siRNAs_probabilities_only_silincing(siRNA_dic):

		all_siRNAs_probs = []
		ompA_mRNA = "GCCAGGGGUGCUCGGCAUAAGCCGAAGAUAUCGGUAGAGUUAAUAUUGAGCAGAUCCCCCGGUGAAGGAUUUAACCGUGUUAUCUCGUUGGAGAUAUUCAUGGCGUAUUUUGGAUGA"
		ompA = "GCCAGGGGTGCTCGGCATAAGCCGAAGATATCGGTAGAGTTAATATTGAGCAGATCCCCCGGTGAAGGATTTAACCGTGTTATCTCGTTGGAGATATTCATGGCGTATTTTGGATGA"
		
		for siRNA in siRNA_dic:
			
			prob80, prob90 = get_probability(siRNA)
			
			final_siRNA1 = siRNA
			
			final_siRNA2 = revcomp(final_siRNA1)
			
			all_siRNAs_probs.append((round(prob80, 2), final_siRNA1.upper(), final_siRNA2.upper()))
		
		return all_siRNAs_probs


def get_siRNAs_probabilities_only_silincing_bilefeld(siRNA_dic):
	
		omp = "ATGA"
		hfq = "GAAA"
		vec = "CCGG"

		all_siRNAs_probs = []
		
		for siRNA in siRNA_dic:
			
			prob80, prob90 = get_probability(siRNA)
			
			final_siRNA1 = omp + siRNA
			
			final_siRNA2 = vec + revcomp(siRNA)
			
			all_siRNAs_probs.append((round(prob80, 2), final_siRNA1.upper(), final_siRNA2.upper()))
		
		return all_siRNAs_probs


__usage__ = """
	Usage:
	python siRNA_commandline.py
	
	required:
		-i <FULL PATH TO INPUT FILE>		gene sequnce
		-m <siRNA | RNAi>			desifn method
		-c					check siRNA
		
	if the check function is used a second input is required:
		-si <siRNA sequence | FULL PATH TO INPUT FILE>		siRNA sequence
	
	optionally:		  
		-v (Bielefeld Method)			if the Bielefeld vector is used
		-o <FULL PATH TO OUTPUT FILE>		output file
	
	All input files should be in FASTA format:
		> header
		gene/siRNA sequence
	
	
	 """
		
	


if __name__ == '__main__':
	
	allowed_chars = Set('ATGCNatgcn')
	
	input_file = ""
	
	output_file = ""
	
	method = ""
	
	vector= False
	
	arguments = sys.argv
	
	siRNA = ""
	
	check = False
	
	if len(arguments) > 1:
		
		i = 1
		
		while i < len(arguments):
			
			if arguments[i] == "-i":
				
				input_file = arguments[i+1]
				
				i=i+2
			
			elif arguments[i] == "-o":
				
				output_file = arguments[i+1]
				
				i=i+2
			
			elif arguments[i] == "-m":
				
				method = arguments[i+1]
				
				i=i+2
			
			elif arguments[i] == "-v":
				
				vector = True
				
				i=i+1
			
			elif arguments[i] == "-c":
				
				check = True
				
				i=i+1
			
			elif arguments[i] == "-si":
				
				siRNA = arguments[i+1]
				
				i=i+2
			
			elif arguments[i] == "help":
				
				sys.exit( __usage__ )				
			
	else:
		
		sys.exit(__usage__)	 
	
	geneseqs = read_genesequence(input_file)
	
	siRNAs = {}
	
	if not Set(siRNA).issubset(allowed_chars):
		
		siRNAs = read_genesequence(input_file)
		
		
	if check:
		
		if len(siRNAs.keys()) ==1:
			
			siRNAseq = siRNAs[siRNAs.keys()[0]]
			
			if Set(siRNAseq).issubset(allowed_chars):
				
				if len(geneseqs.keys()) ==1:
		
					geneseq = geneseqs[geneseqs.keys()[0]]
				
					if  Set(geneseq).issubset(allowed_chars):
				
						if vector:
						
							if method == "siRNA":
								
								siRNA = siRNAseq
								
								if len(siRNA) == 23:
								
									sicore = revcomp(siRNA[4:])
									matches = [m.start() for m in re.finditer('(?=' + sicore +')', geneseq.upper())]
									
									output = "The siRNA matches " + str(len(matches)) + " time/s within the entered gene sequence.\n \nIt matches at position: \n "
									
									for match in matches:
										
										output = output + str(match) + '\n'
									
									rational_score = test_rational_siRNA(siRNA[4:])
									
									vector_overlap = siRNA[:4]
									
									if vector_overlap == "ATGA":
									
										output = output + "\nThe entered siRNA has the correct vector overlap.\n \n"
									
									if rational_score > 6:
							
										output = output + "The entered siRNA fullfills the requirements of the Rational Design \n \n"
					
									ui_score = test_Ui_Tei_siRNA(siRNA[4:])
					
									if ui_score == 4:
						
										output = output + "The entered siRNA fullfills the requirements of the Ui Tei Design \n \n"
									
									if output_file != "":
									
										with open (output_file, 'w') as out:
											
											out.write(output)
									
									sys.exit(output)
									
								else:
									
									sys.exit( "The siRNA has not the correct size.")
								
								
							elif method == "RNAi":
								
								siRNA = siRNAseq
								
								if len(siRNA) == 27:
								
									sicore = revcomp(siRNA[8:])
									matches = [m.start() for m in re.finditer('(?=' + sicore +')', geneseq.upper())]
									
									output = "The siRNA matches " + str(len(matches)) + " time/s with the entered gene sequence.\n \nIt matches at position: \n"
									
									for match in matches:
										
										output = output + str(match) + '\n'
									
									vector_overlap = siRNA[:4]
									
									rnai_overlap = siRNA[4:8]
									
									if vector_overlap == "TAGC":
									
										output = output + "\nThe entered siRNA has the correct vector overlap.\n \n"
									
									if rnai_overlap[0] == 'A' and rnai_overlap[1] == 'G':
									
										output = output + "The entered siRNA has the right RNAi overlap.\n \n"
									
									rational_score = test_rational_siRNA(siRNA[8:])
					
									if rational_score > 6:
							
										output = output + "The entered siRNA fullfills the requirements of the Rational Design \n \n"
					
									ui_score = test_Ui_Tei_siRNA(siRNA[8:])
					
									if ui_score == 4:
						
										output = output + "The entered siRNA fullfills the requirements of the Ui Tei Design \n \n"
									
									if output_file != "":
									
										with open (output_file, 'w') as out:
											
											out.write(output)
									
									sys.exit(output)
									
								else:
									sys.exit( "The siRNA has not the correct size.")
							
							else:
									
								sys.exit("Invalid or missing -m argument.")
						
						else:
							
							if method == "siRNA":
								
								siRNA = siRNAseq
								
								if len(siRNA) == 19:
								
									sicore = revcomp(siRNA)
									matches = [m.start() for m in re.finditer('(?=' + sicore +')', geneseq.upper())]
									
									output = "The siRNA matches " + str(len(matches)) + " time/s with the entered gene sequence.\n \nIt matches at position: \n"
									
									for match in matches:
										
										output = output + str(match) + '\n'
									
									rational_score = test_rational_siRNA(siRNA)	
									
									if rational_score > 6:
							
										output = output + "\nThe entered siRNA fullfills the requirements of the Rational Design \n \n"
					
									ui_score = test_Ui_Tei_siRNA(siRNA)
					
									if ui_score == 4:
						
										output = output + "The entered siRNA fullfills the requirements of the Ui Tei Design \n \n"
									
									if output_file != "":
									
										with open (output_file, 'w') as out:
											
											out.write(output)
									
									sys.exit(output)
								
								else:
									sys.exit( "The siRNA has not the correct size.")
							
							elif method == "RNAi":
								
								siRNA = siRNAseq
								
								if len(siRNA) == 23:
								
									sicore = revcomp(siRNA[4:])
									matches = [m.start() for m in re.finditer('(?=' + sicore +')', geneseq.upper())]
									
									output = "The siRNA matches " + str(len(matches)) + " time/s with the entered gene sequence.\n \nIt matches at position: \n"
									
									for match in matches:
										
										output = output + str(match) + '\n'
									
									rational_score = test_rational_siRNA(siRNA[4:])
									
									rnai_overlap = siRNA[:4]
									
									if rnai_overlap[0] == 'A' and rnai_overlap[1] == 'G':
									
										output = output + "\nThe entered siRNA has the right RNAi overlap.\n \n"
									
									if rational_score > 6:
							
										output = output + "The entered siRNA fullfills the requirements of the Rational Design \n \n"
					
									ui_score = test_Ui_Tei_siRNA(siRNA[4:])
					
									if ui_score == 4:
						
										output = output + "The entered siRNA fullfills the requirements of the Ui Tei Design \n \n"
									
									if output_file != "":
									
										with open (output_file, 'w') as out:
											
											out.write(output)
									
									sys.exit(output)
								
								else:
									sys.exit( "The siRNA has not the correct size.")
							
							else:
									
								sys.exit("Invalid or missing -m argument.")
				
					else:
			
						sys.exit("The entered file contain invalid characters.")
						
				
				else:
					
					sys.exit("Multiple sequences as input is not yet supported.")
					
		
		elif Set(siRNA).issubset(allowed_chars):
			
			if len(geneseqs.keys()) ==1:
		
				geneseq = geneseqs[geneseqs.keys()[0]]
			
				if  Set(geneseq).issubset(allowed_chars):
		
					if vector:
						
						if method == "siRNA":
							
							
							if len(siRNA) == 23:
								sicore = revcomp(siRNA[4:])
								matches = [m.start() for m in re.finditer('(?=' + sicore +')', geneseq.upper())]
								
								output = "The siRNA matches " + str(len(matches)) + " time/s within the entered gene sequence.\n \nIt matches at position: \n"
								
								for match in matches:
									
									output = output + str(match) + '\n'
								
								rational_score = test_rational_siRNA(siRNA[4:])
								
								vector_overlap = siRNA[:4]
								
								if vector_overlap == "ATGA":
								
									output = output + "\nThe entered siRNA has the correct vector overlap.\n \n"
								
								if rational_score > 6:
						
									output = output + "The entered siRNA fullfills the requirements of the Rational Design \n \n"
				
								ui_score = test_Ui_Tei_siRNA(siRNA[4:])
				
								if ui_score == 4:
					
									output = output + "The entered siRNA fullfills the requirements of the Ui Tei Design \n \n"
								
								if output_file != "":
									
									with open (output_file, 'w') as out:
										
										out.write(output)
								
								sys.exit(output)
								
							else:
								
								sys.exit( "The siRNA has not the correct size.")
							
							
						elif method == "RNAi":
							
							
							if len(siRNA) == 27:
							
								sicore = revcomp(siRNA[8:])
								matches = [m.start() for m in re.finditer('(?=' + sicore +')', geneseq.upper())]
								
								output = "The siRNA matches " + str(len(matches)) + " time/s with the entered gene sequence.\n \nIt matches at position: \n"
								
								for match in matches:
									
									output = output + str(match) + '\n'
								
								vector_overlap = siRNA[:4]
								
								rnai_overlap = siRNA[4:8]
								
								if vector_overlap == "TAGC":
								
									output = output + "\nThe entered siRNA has the correct vector overlap.\n \n"
								
								if rnai_overlap[0] == 'A' and rnai_overlap[1] == 'G':
								
									output = output + "The entered siRNA has the right RNAi overlap.\n \n"
								
								rational_score = test_rational_siRNA(siRNA[8:])
				
								if rational_score > 6:
						
									output = output + "The entered siRNA fullfills the requirements of the Rational Design \n \n"
				
								ui_score = test_Ui_Tei_siRNA(siRNA[8:])
				
								if ui_score == 4:
					
									output = output + "The entered siRNA fullfills the requirements of the Ui Tei Design \n \n"
								
								if output_file != "":
									
									with open (output_file, 'w') as out:
										
										out.write(output)
								
								sys.exit(output)
								
							else:
								sys.exit( "The siRNA has not the correct size.")
						
						else:
								
							sys.exit("Invalid or missing -m argument.")
					
					else:
						
						if method == "siRNA":
							
							if len(siRNA) == 19:
							
								sicore = revcomp(siRNA)
								matches = [m.start() for m in re.finditer('(?=' + sicore +')', geneseq.upper())]
								
								output = "The siRNA matches " + str(len(matches)) + " time/s with the entered gene sequence.\n \nIt matches at position: \n"
								
								for match in matches:
									
									output = output + str(match) + '\n'
								
								rational_score = test_rational_siRNA(siRNA)	
								
								if rational_score > 6:
						
									output = output + "\nThe entered siRNA fullfills the requirements of the Rational Design \n \n"
				
								ui_score = test_Ui_Tei_siRNA(siRNA)
				
								if ui_score == 4:
					
									output = output + "The entered siRNA fullfills the requirements of the Ui Tei Design \n \n"
								
								if output_file != "":
									
									with open (output_file, 'w') as out:
										
										out.write(output)
								
								sys.exit(output)
							
							else:
								sys.exit( "The siRNA has not the correct size.")
						
						elif method == "RNAi":
							
							if len(siRNA) == 23:
							
								sicore = revcomp(siRNA[4:])
								matches = [m.start() for m in re.finditer('(?=' + sicore +')', geneseq.upper())]
								
								output = "The siRNA matches " + str(len(matches)) + " time/s with the entered gene sequence.\n \nIt matches at position: \n"
								
								for match in matches:
									
									output = output + str(match) + '\n'
								
								rational_score = test_rational_siRNA(siRNA[4:])
								
								rnai_overlap = siRNA[:4]
								
								if rnai_overlap[0] == 'A' and rnai_overlap[1] == 'G':
								
									output = output + "\nThe entered siRNA has the right RNAi overlap.\n \n"
								
								if rational_score > 6:
						
									output = output + "The entered siRNA fullfills the requirements of the Rational Design \n \n"
				
								ui_score = test_Ui_Tei_siRNA(siRNA[4:])
				
								if ui_score == 4:
					
									output = output + "The entered siRNA fullfills the requirements of the Ui Tei Design \n \n"
								
								if output_file != "":
									
									with open (output_file, 'w') as out:
										
										out.write(output)
								
								sys.exit(output)
							
							else:
								sys.exit( "The siRNA has not the correct size.")
						
						else:
								
							sys.exit("Invalid or missing -m argument.")
	
				else:
		
					sys.exit("The entered file contain invalid characters.")
						
				
			else:
					
				sys.exit("Multiple sequences as input is not yet supported.")
	
	
	elif len(geneseqs.keys()) ==1:
		
		geneseq = geneseqs[geneseqs.keys()[0]]
	
		if  Set(geneseq).issubset(allowed_chars):
		
			if vector:
				
				if method == "siRNA":
					
					if len(geneseq) >=250:
						rational_candidates = get_rational_siRNA(geneseq)
						uitei_candidates = get_Ui_Tei_siRNA(geneseq)
				
						all_siRNAS = rational_candidates.copy()
						all_siRNAS.update(uitei_candidates)
				
						final_list = get_siRNAs_probabilities_only_silincing_bilefeld(all_siRNAS)
						final_list.sort()
						
						output = "Sense siRNA \t Antisense siRNA \t Probability \n" 
						
						for el in reversed(final_list):
							
							output = output + el[1] + '\t' + el[2] + '\t' + str(el[0]) + '\n'
						
						
						if output_file != "":
							
							with open (output_file, 'w') as out:
								
								out.write(output)
						
						sys.exit("Results: \n" + output)
					
					
					elif len(geneseq) <250 and len(geneseq) > 30:
				
						rational_candidates = get_rational_siRNA_short_sequence(geneseq)
						uitei_candidates = get_Ui_Tei_siRNA_short_sequence(geneseq)
				
						all_siRNAS = rational_candidates.copy()
						all_siRNAS.update(uitei_candidates)
				
						final_list = get_siRNAs_probabilities_only_silincing_bilefeld(all_siRNAS)
						final_list.sort()
						
						output = "Sense siRNA \t Antisense siRNA \t Probability \n" 
						
						for el in reversed(final_list):
							
							output = output + el[1] + '\t' + el[2] + '\t' + str(el[0]) + '\n'
						
						if output_file != "":
							
							with open (output_file, 'w') as out:
								
								out.write(output)
						
						sys.exit("Results: \n" + output)
					
					else:
						
						sys.exit("The input genesequence is not long enough to find possible siRNAs.")
				
				elif method == "RNAi":
					
					if len(geneseq) >=250:
						rational_candidates = get_rational_siRNA(geneseq)
						uitei_candidates = get_Ui_Tei_siRNA(geneseq)
		
						all_siRNAS = rational_candidates.copy()
						all_siRNAS.update(uitei_candidates)
		
						final_list = get_siRNAs_probabilities_bielefeld(all_siRNAS)
						final_list.sort()
						
						output = "Sense siRNA \t Antisense siRNA \t Probability \n" 
						
						for el in reversed(final_list):
							
							output = output + el[1] + '\t' + el[2] + '\t' + str(el[0]) + '\n'
						
						if output_file != "":
							
							with open (output_file, 'w') as out:
								
								out.write(output)
						
						sys.exit("Results: \n" + output)
					
					elif len(geneseq) <250 and len(geneseq) > 30:
			
						rational_candidates = get_rational_siRNA_short_sequence(geneseq)
						uitei_candidates = get_Ui_Tei_siRNA_short_sequence(geneseq)
		
						all_siRNAS = rational_candidates.copy()
						all_siRNAS.update(uitei_candidates)
		
						final_list = get_siRNAs_probabilities_bielefeld(all_siRNAS)
						final_list.sort()
						
						output = "Sense siRNA \t Antisense siRNA \t Probability \n" 
						
						for el in reversed(final_list):
							
							output = output + el[1] + '\t' + el[2] + '\t' + str(el[0]) + '\n'
						
						if output_file != "":
							
							with open (output_file, 'w') as out:
								
								out.write(output)
						
						sys.exit("Results: \n" + output)
					
					else:
						
						sys.exit("The input genesequence is not long enough to find possible siRNAs.")
				
				else:
						
					sys.exit("Invalid or missing -m argument.")
			
			else:
				
				if method == "siRNA":
					
					if len(geneseq) >=250:
						
						rational_candidates = get_rational_siRNA_silincing(geneseq)
						uitei_candidates = get_Ui_Tei_siRNA_silincing(geneseq)
				
						all_siRNAS = rational_candidates.copy()
						all_siRNAS.update(uitei_candidates)
				
						final_list = get_siRNAs_probabilities_only_silincing(all_siRNAS)
						final_list.sort()
						
						output = "Sense siRNA \t Antisense siRNA \t Probability \n" 
						
						for el in reversed(final_list):
							
							output = output + el[1] + '\t' + el[2] + '\t' + str(el[0]) + '\n'
						
						if output_file != "":
							
							with open (output_file, 'w') as out:
								
								out.write(output)
						
						sys.exit("Results: \n" + output)
					
					elif len(geneseq) <250 and len(geneseq) > 30:
						
						rational_candidates = get_rational_siRNA_short_sequence_silincing(geneseq)
						uitei_candidates = get_Ui_Tei_siRNA_short_sequence_silincing(geneseq)
				
						all_siRNAS = rational_candidates.copy()
						all_siRNAS.update(uitei_candidates)
				
						final_list = get_siRNAs_probabilities_only_silincing(all_siRNAS)
						final_list.sort()
						
						output = "Sense siRNA \t Antisense siRNA \t Probability \n" 
						
						for el in reversed(final_list):
							
							output = output + el[1] + '\t' + el[2] + '\t' + str(el[0]) + '\n'
						
						if output_file != "":
							
							with open (output_file, 'w') as out:
								
								out.write(output)
						
						sys.exit("Results: \n" + output)
							
					
					else:
						
						sys.exit("The input genesequence is not long enough to find possible siRNAs.")
				
				
				elif method == "RNAi":
					
					if len(geneseq) >=250:
						rational_candidates = get_rational_siRNA(geneseq)
						uitei_candidates = get_Ui_Tei_siRNA(geneseq)
				
						all_siRNAS = rational_candidates.copy()
						all_siRNAS.update(uitei_candidates)
				
						final_list = get_siRNAs_probabilities(all_siRNAS)
						final_list.sort()
						
						output = "Sense siRNA \t Antisense siRNA \t Probability \n" 
						
						for el in reversed(final_list):
							
							output = output + el[1] + '\t' + el[2] + '\t' + str(el[0]) + '\n'
						
						if output_file != "":
							
							with open (output_file, 'w') as out:
								
								out.write(output)
						
						sys.exit("Results: \n" + output)
					
					elif len(geneseq) <250 and len(geneseq) > 30:
			
						rational_candidates = get_rational_siRNA_short_sequence(geneseq)
						uitei_candidates = get_Ui_Tei_siRNA_short_sequence(geneseq)
				
						all_siRNAS = rational_candidates.copy()
						all_siRNAS.update(uitei_candidates)
				
						final_list = get_siRNAs_probabilities(all_siRNAS)
						final_list.sort()
						
						output = "Sense siRNA \t Antisense siRNA \t Probability \n" 
						
						for el in reversed(final_list):
							
							output = output + el[1] + '\t' + el[2] + '\t' + str(el[0]) + '\n'
						
						if output_file != "":
							
							with open (output_file, 'w') as out:
								
								out.write(output)
						
						sys.exit("Results: \n" + output)
					
					else:
						
						sys.exit("The input genesequence is not long enough to find possible siRNAs.")
				
				else:
						
					sys.exit("Invalid or missing -m argument.")
		
		else:
			
			sys.exit("The entered file contain invalid characters.")
	
	else:
		
		sys.exit("Multiple sequences as input is not yet supported.")
