#READ ME
HoP(Host of Phage) is a computational tool which integrates two modules respectively using the deep learning and the Markov chain model to identify the host of a given phage fragment from metagenome or metavirome data at the genus level.  

##How to use HoP

###Dependencies

The codes of HoP are implemented on Python 3.6. To use it, Python 3.6 together with the following packages are required.

- Python 3.6

 - [numpy 1.17.4](https://pypi.org/project/numpy/)
 - [pandas 0.25.3](https://pypi.org/project/pandas/)
 - [biopython 1.71](https://pypi.org/project/biopython/)
 - [numba 0.46.0](https://pypi.org/project/numba/)
 - [sklearn 0.0](https://scikit-learn.org/stable/index.html)
 - [pytorch 1.3.0](https://pytorch.org/)

We recommend to use [Anaconda](https://www.anaconda.com/) to install python and use `conda` to install all dependencies except pytorch. Just simply run:

	conda intstall numpy=1.17.4
	conda intstall pandas=0.25.3
	conda intstall biopython=1.71
	conda intstall numba=0.46.0
	conda intstall sklearn=0.0
	
The version of pytorch may affect the use of the model, and the version used by HoP is 1.3.0. You can just follow the If your machine has a GPU, please configure the corresponding CUDA, CUDNN. Then you can check your own cuda version by `nvidia-smi`, and [downloading corresponding pytorch and torchvision](https://download.pytorch.org/whl/cu100/torch_stable.html), and then install them by `pip`.

	pip install torch-1.3.0+cu100-cp36-cp36m-linux_x86_64.whl
	pip install torchvision-0.4.1+cu100-cp36-cp36m-linux_x86_64.whl


###Downloading
	
After you have configured the environment required to run HoP, you can download the relevant files in the following ways.

	git clone https://github.com/jie-tan/HoP.git

or

	wget https://github.com/jie-tan/HoP/archive/HoP-master.zip
	unzip HoP-master.zip
	cd HoP-master


###Data preparation

Before using HoP to predict host of phage fragments, you need to do some data preparation, 

- [PPR-Meta](http://cqb.pku.edu.cn/ZhuLab/PPR_Meta/)

	PPR-Meta is designed to identify metagenomic sequences as phages, chromosomes or plasmids. If you are getting metagenomic data rather than metavirome data, it is recommended to use PPR-Meta or other tools that can identify short phage sequence fragments to screen out the phage sequence. 

- [Prodigal](https://github.com/hyattpd/Prodigal)

	Prodigal is a fast and reliable protein-coding gene prediction tool. Although prodigal is designed for  prokaryotic genomes, we have found that it has the best performance on short phage sequence fragments during actual use, so we still recommend using it for gene prediction of phage fragments in metagenome or metavriome.

#
	prodigal -i my.metagenome.fna -o my.genes -a my.proteins.faa -p meta

- [taxize](https://www.rdocumentation.org/packages/taxize/versions/0.9.99)

	Taxize is an R package for taxonomic information from around the web. If you want to use HoP with your own candidate host range, a genera list are required to provide. Since you may lack the taxonomy information of prokaryotes and only have its GenBankID or RefSeqID, taxize is recommaded to use for batch processing.

#
	library(taxize)
	# Get NCBI taxonomy UID from GenBankID. 
	# GenBank accession number, gi number, RefSeq accession number can be the input of genbank2uid.
	uid <- genbank2uid(id='AJ748748') 
	# Retrieve the taxonomic hierarchy for a given taxon UID.
	temp <- classification(uid, db = 'ncbi')
	ind <- which(temp[[1]]$rank=='genus')        # genus as an example
	genus <- temp[[1]]$name[ind]
	 

###Usage

	python predict.py [-h] -q QUERY_PHAGE -c QUERY_PHAGE_CDS [-o OUTPUT_DIR]
                  [-w WEIGHT_HOP_S] [-g GENUS_RANGE] [--all]

#####Options

	-h, --help          show this help message and exit
	-q QUERY_PHAGE      A directory containing all single query phage fragments with .fasta/.fna suffix 
						OR a files with .fasta/.fna suffix which contains all query phage fragments
	-c QUERY_PHAGE_CDS  A directory containing all cds output files with .fasta/.fna suffix 
						of single query phage fragment predicted by Prodigal 
						OR a files with .fasta/.fna suffix 
						which contains all cds output of query phage fragments
	-o OUTPUT_DIR       Output directory. The default is the current path
	-w WEIGHT_HOP_S     Weight of HoP-S. Default = 0.5
	-g GENUS_RANGE      A file containing host genera names of interest
	--all               If set, scores of all genera are outputted

#####Example
	mkdir output_example
	python predict.py -q examples/phage_frag_10.fna -c examples/phage_cds_10.fna -o output_example --all

##Citation

##Contact
If you find any bugs or encounter any problems while using HoP, please feel free to contact <jie_tan@pku.edu.cn>.


