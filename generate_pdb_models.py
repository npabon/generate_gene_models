# a Python script to generate representative structures of a gene
# from a specific entry in the Protein Data Bank

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

def fetch_iso_pdb_models(gene_name, destination, pdb_id):
	'''
	Top level function. Generates representative iso 
	models of a particular gene from a PDB structure of that gene.
	'''

	struct_file, molinfo_file, pdbinfo_file  = get_pdb(gene_name, pdb_id, destination)
	rep_pdb_models = generate_pdb_models(gene_name, struct_file, molinfo_file)

	# sanitize the pdb models
	for model in rep_pdb_models:
		sanitize_pdb_model(model)
	
	print '\nFound', len(rep_pdb_models), 'representative pdb models:\n'
	for i, model in enumerate(rep_pdb_models):
		print str(i+1) + ')\t' + os.path.basename(model)

	return rep_pdb_models


def get_pdb(gene_name, pdb_id, dest):
	'''
	Download a pdb file and its description from the website

	Arguments:
	gene_name -- the name of the gene
	pdbid -- the pdb id
	dest -- where you want it saved (directory)

	Return:
	output_struct -- saved pdb structure file name
	output_desc -- saved pdb info file name
	'''

	model_name = gene_name + '_' + pdb_id
	model_dir = dest + '/' + model_name
	output_struct = model_dir + '/' + model_name + '.pdb'
	output_molinfo = model_dir + '/' + model_name + '_molinfo.xml'
	output_pdbinfo = model_dir + '/' + model_name + '_pdbinfo.xml'

	# get the structure file
	if not os.path.exists(output_struct):
		#down load the file
		struct_url = 'http://www.rcsb.org/pdb/files/{pdb}.pdb'.format(pdb=pdb_id)
		print '\ngetting:', struct_url
		pdb_file = urllib.urlretrieve(struct_url)

		#move it to where you want it
		if not os.path.exists(dest):
			os.mkdir(dest)
		if not os.path.exists(model_dir):
			os.mkdir(model_dir)
		shutil.move(pdb_file[0], output_struct)
	else:
		print '\n', output_struct, 'already exists, skipping'

	# get the molecule description and write it to a file
	if not os.path.exists(output_molinfo):
		# get the XML description of the molecule
		molinfo_url = 'http://www.rcsb.org/pdb/rest/describeMol?structureId={pdb}'.format(pdb=pdb_id)
		print 'getting:', molinfo_url
		molinfo_xml = urllib2.urlopen(molinfo_url).read()

		# write XML description to file
		molinfo_file = open(output_molinfo, 'w')
		molinfo_file.write(molinfo_xml)
		molinfo_file.close()
	else:
		print output_molinfo, 'already exists, skipping'

	# get the pdb description and write it to a file
	if not os.path.exists(output_pdbinfo):
		# get the XML description of the pdb
		pdbinfo_url = 'http://www.rcsb.org/pdb/rest/describePDB?structureId={pdb}'.format(pdb=pdb_id)
		print 'getting:', pdbinfo_url
		pdbinfo_xml = urllib2.urlopen(pdbinfo_url).read()

		# write XML description to file
		pdbinfo_file = open(output_pdbinfo, 'w')
		pdbinfo_file.write(pdbinfo_xml)
		pdbinfo_file.close()
	else:
		print output_pdbinfo, 'already exists, skipping'

	return output_struct, output_molinfo, output_pdbinfo

def generate_pdb_models(gene_name, struct_file, desc_file):
	'''
	Generate models to dock against from an original pdb structure
	of a gene target

	Arguments:
	gene_name -- uniprot primary name of the gene target
	struct_file -- the path to the original pdb file
	desc_file -- the path to the pdb description xml file

	Returns:
	models -- a list of the pdb models generated for the pdb
	'''

	pdb_models = []

	# get identify the chains corresponding to target gene
	gene_chains = identify_chains(gene_name, desc_file)	

	# save models for each relevant chain (just isolated chain)
	for chain_id in gene_chains:
		chain_models = save_models(struct_file, chain_id)
		pdb_models += chain_models

	# eliminate redundant models
	pruned_pdb_models = prune_pdb_models(pdb_models)

	return pruned_pdb_models

def identify_chains(gene_name, desc_file):
	'''
	Identifies which chains in a PDB file correspond to the target
	gene of interest

	Arguments:
	gene_name -- uniprot primary name of the gene target
	desc_file -- the path to the pdb description xml file

	Returns:
	gene_chains -- list of chain IDs corresponding to gene_name
	'''

	gene_chains = []

	# prase the xml description
	desc_xml = open(desc_file).read()
	desc_tree = ET.fromstring(desc_xml)

	# identify gene chains
	polymers = desc_tree.find('structureId').getchildren()
	for polymer in polymers:
		molecule = polymer.find('macroMolecule')
		# get accession numbers of polymer chains
		if molecule is not None:
			acc_id = molecule.find('accession').get('id')
			# lookup uniprot accession # to get gene name
			uniprot_url = 'http://www.uniprot.org/uniprot/{acc_id}.xml'.format(acc_id=acc_id)
			acc_id_xml = urllib2.urlopen(uniprot_url).read()
			root = ET.fromstring(acc_id_xml)
			acc_id_entry = root.find('{http://uniprot.org/uniprot}entry')
			acc_id_gene = acc_id_entry.find('{http://uniprot.org/uniprot}gene')
			if acc_id_gene is not None:
				acc_id_primary_gene_name = acc_id_gene.find("./*[@type='primary']")
				if acc_id_primary_gene_name is not None:
					acc_id_primary_gene_name = acc_id_primary_gene_name.text
					# if it's the gene we want, save the chain IDs
					if acc_id_primary_gene_name == gene_name:
						for chain in polymer.iter('chain'):
							gene_chains.append(chain.get('id'))

	return gene_chains

def save_models(struct_file, chain_id):
	'''
	Save target models from original pdb structure file using PyMol and the
	chain corresponding to our target gene of interest. Saves an isolated
	model of the chain

	Arguments:
	struct_file -- the path to the original pdb file
	chain_id -- pdb chain ID of the target gene chains 

	Returns:
	model_names -- list of the pdb models generated for the given chain
	'''

	model_names = []

	template_scpt ='''
from pymol import cmd
import sys
import os
import datetime

# load the raw pdb file
cmd.load("{struct_file}")

# define selections for isolated chain models
iso_sele = "{iso_sele}"

# save the isolated model
cmd.save("{iso_outfile}", iso_sele)
'''

	basename = os.path.splitext(os.path.basename(struct_file))[0] # ex. HDAC6_3PDH
	model_dir_template = os.path.splitext(struct_file)[0] # ex. HDAC6/HDAC6_3PDH/HDAC6_3PDH

	# Prepare directory and selection to save isolated chain model
	iso_sele = '(polymer or resn ZN+K+CA) and chain {chain_id}'.format(chain_id = chain_id)
	iso_model_dir = model_dir_template + '_' + chain_id + '_iso/' # ex. HDAC6/HDAC6_3PHD/HDAC6_3PHD_A_iso/
	if not os.path.exists(iso_model_dir):
		os.mkdir(iso_model_dir)
	iso_outfile = iso_model_dir + basename + '_' + chain_id + '_iso.pdb'
	model_names.append(iso_outfile)

	# modify the template script
	pymol_scpt = template_scpt.format(struct_file=struct_file, iso_sele=iso_sele, iso_outfile=iso_outfile)

	# write the script to a temporary file
	pyout = tempfile.NamedTemporaryFile(mode='w+b', suffix='.py', dir='/tmp/', delete=False)
	pyout.write(pymol_scpt)
	pyout.close()

	print 'Saving models: ', os.path.basename(iso_outfile)

	# call pymol and run the script
	save_models_cmd = 'pymol -qck -r ' + pyout.name
	if sp.call(save_models_cmd, shell=True) == 1:
		print '\nError When Running:' + gen_pse_cmd
		print 'Failed to generate models for:', struct_file, '\n'

	return model_names

def prune_pdb_models(pdb_models):
	'''
	This function takes a list of structural models corresponding to a single
	pdb ID (just isolated models). It prunes them to find representative
	models and eliminates redundant ones

	Arguments: 
	pdb_models -- full list of pdb models (iso)

	Returns:
	pruned_models -- list of pruned representative pdb models
	'''
	pruned_models = []

	# determine which files actually exist, delete parent dirs of those that don't
	iso_pdb_models = []

	for model in pdb_models:
		if not os.path.exists(model):
			print os.path.basename(model), 'does not exist! Deleting parent directory.'
			delete_model(model)
		else:
			iso_pdb_models.append(model)


	# find representative models
	rep_overlap_cutoff = 50 # percent seq overlap required (90% seq ID required)
	rep_rmsd_cutoff = 5 # models less than 4A apart are represented by a single model

	# find representative iso models
	print 'Finding representative PDB ISO models...'
	rep_iso_models = []

	for iso_model in iso_pdb_models:
		if len(rep_iso_models) == 0:
			rep_iso_models.append(iso_model)
		else:
			model = prody.parsePDB(iso_model) # get structure
			redundant = False
	
			for rep_iso_model in rep_iso_models:
				rep = prody.parsePDB(rep_iso_model) # get structure

				# calc RMSD between model and rep
				alignment = prody.matchAlign(model,rep,overlap=rep_overlap_cutoff)
				if alignment != None:
					rmsd = prody.calcRMSD(alignment[1], alignment[2])
					if rmsd <= rep_rmsd_cutoff:
						redundant = True # we already have a representative for this segment
						# take the larger structure as the representative
						if model.numResidues() > rep.numResidues():
							rep_iso_models.remove(rep_iso_model)
							rep_iso_models.append(iso_model)
						break

			# if the iso model does not match any of our representative models,
			# then add it to the representative models list
			if not redundant:
				rep_iso_models.append(iso_model)

	print 'Found', len(rep_iso_models), 'representative ISO models:', map(os.path.basename, rep_iso_models)


	# move representative models to their own directory
	if len(rep_iso_models) > 0: 
		pdb_dir = os.path.abspath(os.path.join(rep_iso_models[0], os.pardir+'/'+os.pardir))
		rep_model_dir = pdb_dir + '/representative_pdb_models/'
		if os.path.exists(rep_model_dir):
			shutil.rmtree(rep_model_dir)
		os.mkdir(rep_model_dir)

		for rep_iso_model in rep_iso_models:
			rep_iso_model_pardir = os.path.abspath(os.path.join(rep_iso_model, os.pardir))
			new_path = rep_model_dir + '/' + os.path.basename(rep_iso_model_pardir)
			shutil.copytree(rep_iso_model_pardir, new_path)
			# define new pathname to keep track of the models once we move them
			new_iso_model_path = rep_model_dir + os.path.basename(rep_iso_model_pardir) + '/' + os.path.basename(rep_iso_model)
			pruned_models.append(new_iso_model_path)


	# return all representative pdb models
	return pruned_models

def delete_model(model_file):
	'''
	This function deletes the parent directory containing a specific model

	Arguments:
	model_file -- path to the pdb model file to be deleted
	'''

	parent_dir = os.path.abspath(os.path.join(model_file, os.pardir))
	shutil.rmtree(parent_dir)

def sanitize_pdb_model(model):
	'''
	This function takes an input model corresponding to a PDB file
	and then manually removes problematic lines such as any ANISOU or 
	CONNECT lines. 

	Also makes sure each residue number is only 3 digits long to avoid
	smina parsing errors

	Arguments:
	model -- path to the pdb model file to be sanitized
	'''

	clean_lines = []

	model_file = open(model)
	for line in model_file:
		line = line.strip()
		# ignore all unwanted lines
		if not line.split()[0] == 'ANISOU':
			if not line.split()[0] == 'CONECT':
				if not line.split()[0] == 'MASTER':
					# reduce residue number to 3 digits
					newline = line[:22] + ' ' + line[23:]
					clean_lines.append(newline)

	# overwrite the old model with the sanitized model
	sanitized_model_file = open(model, 'w')
	for line in clean_lines:
		sanitized_model_file.write("%s\n" % line)
	sanitized_model_file.close()
				


def main():
	''' Parse command line args and run.'''
	
	parser = argparse.ArgumentParser(description="Generate docking models from a PDB structure of a given gene")
	parser.add_argument("gene", type=str, help="UNIPROT primary gene name")
	parser.add_argument("gene_dir", type=str, help="Parent directory for all the gene's models")
	parser.add_argument("pdb_id", type=str, help="PDB ID of the structure to generate models from")
	args = parser.parse_args()

	gene_name = args.gene # gene name (must be uniprot primary name)
	destination = args.gene_dir # gene target parent directory
	pdb_id = args.pdb_id

	fetch_iso_pdb_models(gene_name, destination, pdb_id)



if __name__ == '__main__':
   main()
