#ABC ID PIPELINE

This pipeline is mean't to analyze the ABC transporter family across a large number (>150) sequenced arthropods. It uses the ABC_scan program to search the predicted protein sets of non*model species. More info on the pipelin can be found in the forthcoming publications (Denecke et. al 2021; BMC Genomics). 

## Required dependencies

* ABC_scan (https://github.com/shanedenecke/ABC_scan)
* R
	* tidyverse
	* data.table
	* ggtree
	* ggsci
	* agricolae
	* gplots
	* ape
	* treeio
* python3
	* Biopython
	* re 
	* os 
	* sys
	* pandas
	
* CAFEv5 (https://github.com/hahnlab/CAFE5/releases)
* Orthofinder (https://github.com/davidemms/OrthoFinder)
* HMMER package (http://hmmer.org/)
* MAFFT (https://mafft.cbrc.jp/alignment/software/)
* BUSCO (https://busco.ezlab.org/)
* Custom_Applications (https://github.com/shanedenecke/Custom_Applications)
* BLAST+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* MEME suite (https://meme-suite.org/meme/)
* The GENERAL_REFERENCE FOLDER (Various large documents such as proteomes available upon request)

## Instructions 



## Disclaimer

This pipeline has only been tested on Fedora31 and Debian 9 based linux systems from one user. It is not designed to run on other machines and has some paths specific to my setup.

If you are interested in adapting this pipeline, plese contact Shane Denecke shanedenecke@gmail.com.