'''Dini Alice 830931'''

#! /usr/bin/python
import sys 
import numpy as np
import itertools

masa={'A':115, 'L':170, 'R':225.0, 'K':200, 'N':160.0, 'M':185.0, 'D':150.0, 'F':210.0, 'C':135.0, 'P':145.0, 'Q':180.0, 'S':115.0, 'E':190.0, 'T':140.0, 'G':75.0, 'W':255.0, 'H':195.0, 'Y':230.0, 'I':175.0, 'V':155.0}


def parse_dssp(dsspfile):
	dssp={}
	fdssp=open(dsspfile)
	c=0
	for line in fdssp:
		if line.find('  #  RESIDUE')==0:
			c=1
			continue
		if c==1: 
			num=line[5:10].strip() 
			res=line[13]
			ch = line[11]
			ss=line[16]
			asa=float(line[35:38])
			phi=float(line[103:109])
			psi=float(line[109:115])
			if ss==' ': ss='C' 
			if ch!=' ':
				dssp[ch] = dssp.get(ch,{})
				dssp[ch][num]=[res,ss,asa,phi,psi]
	return dssp

def rasa(dssp, chain):
	'''Computes the relative solvent accessibility'''
	ks = dssp.get(chain, dssp).keys()
	ks.sort()
	rasa_list=[]
	for k in ks:
		r = dssp[chain][k][2]/masa[dssp[chain][k][0]]
		rasa_list.append([dssp[chain][k][0], k, round(min(1,r),2)])
	return rasa_list

	
def get_asa(dssp, chain): 
	'''To extract individual surfaces of a chain from trimer/tetramer/dimer'''
	s = sum(dssp[chain][j][2] for j in dssp[chain])
	return s



if __name__== '__main__':
	if len(sys.argv)>=2: 
		dsspfile=sys.argv[1]                    #dssp file in input of the tetramer or trimer
		dssp_chain_file = sys.argv[2]           #dssp file of the chain or of the trimer, smaller than the previous
		chain = sys.argv[3]                     #selected chain to be compared
		dssp=parse_dssp(dsspfile)
		ch = parse_dssp(dssp_chain_file)
		
		'''With the function rasa we get the relative solvent accessibility, given a certain chain from a complex.
		To compute the rasa of a chain alone, just input as chain file (3rd argument in the terminal) the specific
		chain dssp and select the chain itself.'''
		
		rasa1 = rasa(dssp, chain) 
		rasa2 = rasa(ch, chain)
		
		'''To compare complex vs monomer rasa, select in the terminal complex.dssp chainX.dssp X 
		where X is the chain.'''
		
		c = 0
		for i in range(len(rasa1)):
			if rasa1[i] != rasa2[i] and abs(rasa1[i][2] - rasa2[i][2])>0.10:
				c+=1
				print rasa1[i][1], rasa1[i][0], rasa1[i][2], rasa2[i][2], abs(rasa1[i][2] - rasa2[i][2]), str(round(100*abs(rasa1[i][2] - rasa2[i][2]),2))+'%'
				
		print 'Total number of changing residues between chain',chain,'in complex vs monomeric:',c
		print 'Difference in surface accessible area:', abs(get_asa(dssp, chain) - get_asa(ch, chain)) 
		

