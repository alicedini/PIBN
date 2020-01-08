'''Dini Alice 830931'''

#! /usr/bin/python
import sys
import networkx as nx
import numpy as np

def search_mitab(mitabfile):
	list_int=[]
	f=open(mitabfile)
	for line in f:
		v=line.split('\t') 
		int1=v[0].replace('uniprotkb:','')
		int2=v[1].replace('uniprotkb:','')
		texp=v[6].replace('psi-mi:','')
		tax1=v[9].replace('taxid:','')
		tax2=v[10].replace('taxid:','')
		cexp=v[11].replace('psi-mi:','')
		pint=(int1, int2, texp, tax1, tax2, cexp)
		list_int.append(pint)
	return list_int

def filter_organism(list_int, organism):
	organism=organism.lower()
	filter_list=[pint for pint in list_int if pint[3].lower().find(organism)>-1 and pint[4].lower().find(organism)>-1]
	return filter_list

def clean(list_int):
	new = {}
	for pint in list_int:
		k = [pint[0], pint[1]]
		k.sort()
		new[tuple(k)] = new.get(tuple(k),0)+1
	return new
				
def stats(ints, p1, p2): 
	g=nx.Graph() 
	g.add_edges_from(ints)
	n=g.number_of_nodes()
	subnet1 = nx.single_source_shortest_path_length(g,source=p1,cutoff=2) 
	subnet2 = nx.single_source_shortest_path_length(g,source=p2,cutoff=2) 
	total_subnet = set(subnet1.keys()).union(set(subnet2.keys()))
	sg = g.subgraph(total_subnet)
	
	#'''Degree of alpha and beta'''
	print 'Degree of',u'\u03b1 subunit:',sg.degree(p1)
	print 'Degree of',u'\u03b2 subunit: ',sg.degree(p2)

	'''Betweenness centrality of alpha and beta in the reduced subnet'''
	betw = nx.betweenness_centrality(sg)	#Computed once for all
	print 'Betweenness centrality of',u'\u03b1 subunit:', betw[p1]
	print 'Betweenness centrality of',u'\u03b2 subunit:', betw[p2]
	
	#'''Transitivity of alpha and beta in the reduced subnet'''
	print 'Transitivity of ',u'\u03b1 subunit:', nx.clustering(sg, p1)
	print 'Transitivity of',u'\u03b2 subunit', nx.clustering(sg, p2)
	
	'''Direct interactions subnet'''
	'''We seek directly interacting proteins to find their betweenness centrality in the overall system
	We are also excluding the two imput and interacting proteins from the results'''
	
	set_prot2 = set([k for k in subnet2.keys() if subnet2[k]==1 and k!=p2])
	set_prot1 = set([k for k in subnet1.keys() if subnet1[k]==1 and k!=p1])
	total_network_dir=set_prot1.union(set_prot2)
	print 'The proteins directly interacting with both the subunits are:'
	for i in set_prot2.intersection(set_prot1):
		print i
	print 'Protein			Degree			Betweenness centrality					Transitivity'
	for x in total_network_dir:
		print x,'			',sg.degree(x),'			',betw[x],'					',nx.clustering(sg, x)
	return 
	
		
if __name__ == '__main__':
        '''Instructions: mitab file to be parsed (eg:intact.txt) protein 1 protein 2 organism'''
        
	mitabfile=sys.argv[1]
	p1=sys.argv[2]
	p2 = sys.argv[3]
	organism=sys.argv[4]
	list_int=search_mitab(mitabfile)
	filter_list=filter_organism(list_int, organism)
	cleaned_list = clean(filter_list)
	stats(cleaned_list, p1, p2)

