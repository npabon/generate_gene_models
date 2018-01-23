# Generate Gene Models

## SYNOPSIS

generate_gene_modles.py is a Python script for rapidly mining the Protein Data Bank and generating representative structural models of a given human gene. The script first downloads all PDB entries associated with the gene (using the accompanying script generate_pdb_models.py). Then the script compares all PDB structures and selects representative structures optimizing for sequence coverage, structural resolution, and structural diversity (i.e. multiple conformations).



## USAGE
```
usage: generate_gene_models.py [-h] gene dest

Mine the PDB for structural models of a gene

positional arguments:
  gene        UNIPROT primary gene name
  dest        Destination parent directory for the gene

optional arguments:
  -h, --help  show this help message and exit
```


## INPUT

```gene``` - The Uniprot Primary gene name for the human gene to generate representative structural models of

```dest``` - The top-level directory in which to store the structural models of gene (will be created if it doesn't already exist)



## OUTPUT

The output model files are stored in the 'dest' directory and are organized by PDB and then by chain. The naming system is <geneID>_<pdbID>_<chainID>. Assume you run the following command:

```python generate_gene_models.py geneX ~/geneX```

The directory ~/geneX would be organized as follows:

```
  ~/geneX/
    representative_gene_models/           * This directory contains the representative gene models
      geneX_pdbID_chainID_iso/
        geneX_pdbID_chainID_iso.pdb       * This is a representative gene model. It is chain <chainID> of the structure of geneX with PDB ID <pdbID>
      ...                                 * There may be multiple representative structural model of geneX
    geneX_pdbID/                          * The rest of the directories at this level contain structures specific to different PDB entries for geneX
      representative_pdb_models/          * This directory contains the representative structures (chains) extracted from PDB <pdbID> 
        geneX_pdbID_chainID_iso/
          geneX_pdbID_chainID_iso.pdb     * This is a representative pdb model. It is chain <chainID> of the structure of geneX with PDB ID <pdbID>
      geneX_pdbID.pdb                     * The original PDB file
      geneX_pdbID_chainID/
        geneX_pdbID_chainID.pdb           * This is chain <chainID> from PDB <pdbID>
      ...                                 * There may be multiple unique chains in <pdbID> that are all geneX
      geneX_pdbID_molinfo.xml             * metadata from the PDB
      geneX_pdbID_pdbinfo.xml             * metadata from the PDB
    ...                                   * There may be multiple pdbIDs corresponding to geneX
```






