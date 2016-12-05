# a Python script to download all the PDB files corresponding
# to the human X-ray structures of a given gene

# Part of the Supplemental Information for:
# Drug Target Interaction Prediction using LINCS Data
# Nicolas A. Pabon, Yan Xia, Sam Estabrooks, Amanda K. Herbrand, Evelyn SuB, Ricardo Biondi, Jeffrey L.  Brodsky, Carlos J. Camacho,  and Ziv Bar-Joseph

# Code written by Nicolas Pabon
# 12/5/2016

import sys
import urllib
import urllib2
import os
import shutil
import xml.etree.ElementTree as ET
import subprocess as sp
import tempfile
import prody
import argparse
import glob
import generate_pdb_models

def get_pdb_ids(gene_name):
	'''
	Query the PDB for human x-ray structures of a given gene_name

	Arguments:
	gene_name -- the UNIPROT gene name (ex. PDPK1)

	Return:
	pdb_ids -- a list of the PDB ids for human x-ray structures of gene_name
	'''

	# fill in the RESTful search template
	# organism (human) and experimental method (X-ray) are hard-coded in template
	query_text = """
	<orgPdbCompositeQuery version="1.0">
	 <queryRefinement>
	  <queryRefinementLevel>0</queryRefinementLevel>
	  <orgPdbQuery>
	    <version>head</version>
	    <queryType>org.pdb.query.simple.UniprotGeneNameQuery</queryType>
	    <description>UniProt Gene Name:  {gene_name}</description>
	    <query>{gene_name}</query>
	  </orgPdbQuery>
	 </queryRefinement>
	 <queryRefinement>
	  <queryRefinementLevel>1</queryRefinementLevel>
	  <conjunctionType>and</conjunctionType>
	  <orgPdbQuery>
	    <version>head</version>
	    <queryType>org.pdb.query.simple.TreeEntityQuery</queryType>
	    <description>TaxonomyTree Search for Homo sapiens (human)</description>
	    <t>1</t>
	    <n>9606</n>
	    <nodeDesc>Homo sapiens (human)</nodeDesc>
	  </orgPdbQuery>
	 </queryRefinement>
	 <queryRefinement>
	  <queryRefinementLevel>2</queryRefinementLevel>
	  <conjunctionType>and</conjunctionType>
	  <orgPdbQuery>
	    <version>head</version>
	    <queryType>org.pdb.query.simple.ExpTypeQuery</queryType>
	    <description>Experimental Method is X-RAY</description>
	    <mvStructure.expMethod.value>X-RAY</mvStructure.expMethod.value>
	  </orgPdbQuery>
	 </queryRefinement>
	</orgPdbCompositeQuery>
	""".format(gene_name = gene_name)

	# query the pdb and compile results into a list of pdb IDs

	print "Querying the PDB for gene: " + gene_name

	url = 'http://www.rcsb.org/pdb/rest/search'
	request = urllib2.Request(url, data=query_text)
	results = urllib2.urlopen(request)
	pdb_ids = []
	for result in results:
		pdb_id = result.strip()
		pdb_ids.append(pdb_id)

	# sanitize the PDB IDs, sometimes they come back as something like 2FUC:1
	# we want to get rid of the ':1'
	pdb_ids = [ i.split(':')[0] for i in pdb_ids ]

	print 'Found ' + str(len(pdb_ids)) + ' PDB IDs for ' + gene_name
	return pdb_ids

def find_representative_gene_models(pdb_models):
	'''
	This function will identify a few representative models of the target
	gene from among all the models generated from the PDB structures associated
	with the gene.

	Arguments:
	pdb_models -- list of all the representative models for the gene (all PDBs)

	Returns:
	gene_models -- the list of representative models for the gene
	'''

	hi_res_cutoff = 2.0 # Angstroms

	hi_res_iso_models = []
	lo_res_iso_models = []

	# find the resolution for each model
	# sort by model type and resolution
	for pdb_model in pdb_models:
		pdb_model_resolution = get_pdb_resolution(pdb_model)
		if pdb_model_resolution <= hi_res_cutoff:
			hi_res_iso_models.append(pdb_model)
		else:
			lo_res_iso_models.append(pdb_model)

	# find representative iso models
	print 'Finding representative gene ISO models...'
	max_num_rep_iso_models = 3
	max_clust_rad = 10.0

	# start with a standard clustering radius of 4A RMSD
	print 'Starting with a default RMSD clustering radius of 4.0 Anstroms'
	clust_rad = 4.0
	rep_gene_iso_models = find_rep_gene_iso_models(hi_res_iso_models, lo_res_iso_models, clust_rad)

	# try to minimize the number of rep models by increasing the clustering radius
	while len(rep_gene_iso_models) > max_num_rep_iso_models and clust_rad < max_clust_rad:
		clust_rad += 1.0
		print 'Max number of rep gene iso models exceeded...'
		print 'Increasing RMSD clustering radius to', clust_rad, 'Anstroms'
		rep_gene_iso_models = find_rep_gene_iso_models(hi_res_iso_models, lo_res_iso_models, clust_rad)

	rep_gene_models = rep_gene_iso_models 


	# move representative gene models to their own directory
	new_rep_gene_models = []
	if len(rep_gene_models) > 0: 
		gene_dir = os.path.abspath(os.path.join(rep_gene_models[0], os.pardir+'/'+os.pardir+'/'+os.pardir+'/'+os.pardir))
		rep_model_dir = gene_dir + '/representative_gene_models'
		if os.path.exists(rep_model_dir):
			shutil.rmtree(rep_model_dir)
		os.mkdir(rep_model_dir)

		for rep in rep_gene_models:
			rep_pardir = os.path.abspath(os.path.join(rep, os.pardir))
			rep_name = os.path.basename(rep_pardir)
			new_path = rep_model_dir + '/' + rep_name
			shutil.copytree(rep_pardir, new_path)

			# record new location
			new_name = new_path + '/' + rep_name + '.pdb'
			new_rep_gene_models.append(new_name)


	return new_rep_gene_models

def get_pdb_resolution(pdb_model):
	'''
	This function takes in a representative pdb model and returns the 
	resolution of the crystal structure read from the pdbinfo XML file.
	Input model must be a representative 
	'''

	pdb_dir = os.path.abspath(os.path.join(pdb_model, os.pardir+'/'+os.pardir+'/'+os.pardir))
	pdbinfo_file = glob.glob(pdb_dir + '/*pdbinfo.xml')[0]
	
	# prase the xml description
	pdbinfo_xml = open(pdbinfo_file).read()
	pdbinfo_tree = ET.fromstring(pdbinfo_xml)

	pdb_entry = pdbinfo_tree.getchildren()[0]
	pdb_res = pdb_entry.get('resolution')
	if pdb_res != None:
		pdb_res = float(pdb_res)
	else:
		pdb_res = 10.0

	return pdb_res

def find_rep_gene_iso_models(hi_res_iso_models, lo_res_iso_models, rep_rmsd_cutoff):
	'''
	Function to pick representative iso gene models from a pool of representative
	pdb iso models. Uses a greedy algorithm to cover as much of the gene sequence
	as possible using first high resolution models and then filling any gaps 
	with low resolution models

	'''

	# only make a model a representative model if is at least 15 residues
	# long and if it includes at least 10
	# residues that have never been seen in previous models or if it has
	# a significantly different conformation than previous models
	rep_gene_iso_models = []
	min_length = 20
	num_new_residue_cutoff =10
	rep_overlap_cutoff = 10

	# tag each model with it's sequence coverage 
	hi_res_iso_models = [ [m, get_seq_range(m)] for m in hi_res_iso_models ]
	lo_res_iso_models = [ [m, get_seq_range(m)] for m in lo_res_iso_models ]

	# sort lists of models by length of sequence coverage 
	sorted_hi_res_iso_models = sorted(hi_res_iso_models, key=lambda m: -1*len(m[1]))
	sorted_lo_res_iso_models = sorted(lo_res_iso_models, key=lambda m: -1*len(m[1]))
	sorted_iso_models = sorted_hi_res_iso_models + sorted_lo_res_iso_models

	# use greedy algorithm to try to cover full gene sequence
	gene_coverage = []

	# start with large hi res models, end with small lo res models
	for model in sorted_iso_models:
		model_file = model[0]
		model_coverage = model[1]
		
		# discrard structures that have too few number of residues
		if len(model_coverage) >= min_length:
			intersection = list(set(model_coverage) & set(gene_coverage))
			num_new_residues = len(model_coverage) - len(intersection)

			# if rep model list is empty, make it a rep model
			if len(rep_gene_iso_models) == 0:
				rep_gene_iso_models.append(model_file)
				gene_coverage += model_coverage

			# otherwise, if this model has enough new residues, add it to the representatives list
			elif num_new_residues >= num_new_residue_cutoff:
				rep_gene_iso_models.append(model_file)
				gene_coverage += model_coverage
				gene_coverage = list(set(gene_coverage))

			# otherwise check if it has a unique conformation
			else:
				model_struct = prody.parsePDB(model_file)
				redundant = False
				for rep_gene_iso_model in rep_gene_iso_models:
					rep_struct = prody.parsePDB(rep_gene_iso_model) # get structure
					# calc RMSD between model and rep
					alignment = prody.matchAlign(model_struct,rep_struct,overlap=rep_overlap_cutoff)
					if alignment != None:
						rmsd = prody.calcRMSD(alignment[1], alignment[2])
						if rmsd <= rep_rmsd_cutoff:
							redundant = True # we already have a representative for this segment
							break
				# if the model does not match any of our representative models,
				# then it is unique - add it to the representative models list
				if not redundant:
					rep_gene_iso_models.append(model_file)
					gene_coverage += model_coverage
					gene_coverage = list(set(gene_coverage))

	return rep_gene_iso_models

def get_seq_range(model):
	'''
	Takes in a pdb model file and returns a list containing the
	sequence indices that the structure spans
	'''

	struct = prody.parsePDB(model)
	unique_indeces = []
 
	if struct is not None:
		for i in struct.getResnums():
			if i not in unique_indeces:
				unique_indeces.append(i)

	return unique_indeces

##### MAIN #####

# initialize argument parser
parser = argparse.ArgumentParser(description="Mine the PDB for structural models of a gene")
parser.add_argument("gene", type=str, help="UNIPROT primary gene name")
parser.add_argument("dest", type=str, help="Destination parent directory for the gene")
args = parser.parse_args()

gene_name = args.gene # gene name (must be uniprot primary name)
destination = args.dest # gene target parent directory

# get PDB IDs for all structures of that gene
pdb_ids = get_pdb_ids(gene_name)
if len(pdb_ids) == 0:
	print 'No structures found. Query Complete'
	sys.exit()

# extract representative models from each PDB structure
rep_pdb_models = []
for pdb_id in pdb_ids:
	rep_pdb_models += fetch_iso_pdb_models.fetch_iso_pdb_models(gene_name, destination, pdb_id)

# compare different representative PDB models and select representative models for the gene itself
rep_gene_models = find_representative_gene_models(rep_pdb_models)
print '\nFound', len(rep_gene_models), 'representative gene models:\n'
for i, model in enumerate(rep_gene_models):
	print str(i+1) + ')\t' + os.path.basename(model)
print '\n'

# if no rep iso models were generated, write down the gene in a log file
if len(rep_gene_models) == 0:
	logfile_name = 'genes_without_rep_models.txt'
	logfile = open(logfile_name, "a")
	logfile.write(gene_name + '\n')
	logfile.close()




