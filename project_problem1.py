'''Dini Alice 830931'''

#! /usr/bin/python
import sys 
import numpy as np
import itertools

def parse_pdb(pdbfile):
	'''Parses pdb file to extract information'''
	fpdb=open(pdbfile)
	dict_pdb ={}
	dict_hetatm={}
	for line in fpdb:
		if line[0:4]=='ATOM' or line[0:6]=='HETATM':
			x=float(line[30:38])
			y=float(line[38:46])
			z=float(line[46:54])
			at=line[12:16].strip() 
			ch = line[21]
			res=line[17:20]
			n=line[22:26].strip()
			coord=[x,y,z]
			if (line[0:4]=='ATOM'):
				dict_pdb[ch]=dict_pdb.get(ch,{})
				dict_pdb[ch][(res,n)] = dict_pdb[ch].get((res,n), {})
				dict_pdb[ch][(res,n)][at] = coord
			if (line[0:6]=='HETATM') and (res != 'HOH'):
				dict_hetatm[ch]=dict_hetatm.get(ch,{})
				dict_hetatm[ch][(res,n)] = dict_hetatm[ch].get((res,n), {})
				dict_hetatm[ch][(res,n)][at] = coord
	return dict_pdb, dict_hetatm

def dist(a1,a2):
	'''Returns the distance between two given atoms'''
	d=np.sqrt((a1[0]-a2[0])**2+(a1[1]-a2[1])**2+(a1[2]-a2[2])**2)
	return d

def heme_oxy_dist(dict_pdb, dict_hetatm, chain):
	'''Returns close-by residues to oxygen groups and heme groups'''
	distances={}
	for i in dict_pdb[chain]:
		for j in dict_hetatm[chain]:
			for k in dict_pdb[chain][i]:
				for w in dict_hetatm[chain][j]:
					d=dist(dict_pdb[chain][i][k], dict_hetatm[chain][j][w])
					if d <3.5:
                                                if distances.get((i,j), 0)==0:
							distances[(i,j)] = [[k, w]]
						else:
							distances[(i,j)].append([k, w])
							
	h=0
	o=0
	hk = sorted(distances.keys())
	for i in hk:
		if i[1][0] == 'HEM':
			h+=1
			print chain, i[1][0], i[0][0], i[0][1], ', '.join(('-'.join(k) for k in distances[i]))
	for j in hk:
		if j[1][0] == 'OXY':
			o+=1
			print chain, j[1][0], j[0][0], j[0][1], ', '.join(('-'.join(k) for k in distances[j]))
	return 'HEM:',h,'OXY:',o


def interacting_res(dict_pdb,chain1, chain2):
	'''Returns residues interacting between two chains given in input'''
	residues={}
	for i in dict_pdb[chain1]:
		for j in dict_pdb[chain2]:
			if i != j:
				for x in dict_pdb[chain1][i]:
					for y in dict_pdb[chain2][j]:
						d = dist(dict_pdb[chain1][i][x], dict_pdb[chain2][j][y])
						if d <3.5:
							if residues.get((i,j), 0)==0:
								residues[(i,j)] = [[x, y]]
							else:
								residues[(i,j)].append([x, y])
							
	rk=sorted(residues.keys())
	for i in rk:
		print chain1, i[0][1], i[0][0],  chain2, i[1][1], i[1][0], ', '.join(('-'.join(k) for k in residues[i]))
	return


if __name__== '__main__':
	if len(sys.argv)>=2:
		'''Requires as input a pdb file, preferentially of the entire complex, and obtains the interacting 
		residues in between two chains, given as well as input. Furthermore, it also prints the interactions
		that each chain has with heme and oxygen groups.'''
		
		pdbfile = sys.argv[1]
		chain1 = sys.argv[2]
		chain2 = sys.argv[3]
		dict_pdb, dict_hetatm=parse_pdb(pdbfile)
		print(heme_oxy_dist(dict_pdb, dict_hetatm, chain1))
		print(heme_oxy_dist(dict_pdb, dict_hetatm, chain2))
		interacting_res(dict_pdb,chain1,chain2)
		
	
		
		
		
		
		
		
