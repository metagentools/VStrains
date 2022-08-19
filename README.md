# vsAware: Graph Algorithms for De Novo Construction of Viral Strains

Manual
===========

Table of Contents
-----------------

1. [About vsAware](#sec1) </br>
2. [Installation](#sec2) </br>
   2.1. [Quick Install](#sec2.1) </br>
   2.1. [Manual Install](#sec2.2) </br>
3. [Running vsAware](#sec3) </br>
   3.1. [Quick Usage](#sec3.1) </br>
   3.2. [Support SPAdes](#sec3.2) </br>
   3.3. [Parameters](#sec3.4) </br>
4. [Citation](#sec4) </br>
5. [Feedback and bug reports](#sec5)</br>

<a name="sec1"></a>
# About vsAware

vsAware is a strain-aware assembly tools, which constructs full-length virus strain with aid from de Bruijn assembly graph and Strain-oblivious assembled contigs.

Please refer to our [paper](NULL) and [supplementary Material](NULL) for details methodology.

<a name="sec2"></a>
# Installation

vsAware requires a 64-bit Linux system or Mac OS and python (supported versions are python3: 3.2 and higher), 

<a name="sec2.1"></a>
## Quick Install (**recommended**)

Install [(mini)conda](https://conda.io/miniconda.html) as a light-weighted package management tool. Run following commands to initialize&setup the conda environment for vsAware

```bash
# add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment
conda create --name vsAware-env

# activate conda environment
conda activate vsAware-env

conda install -c bioconda -c conda-forge python=3 graph-tool minimap2 numpy pandas scikit-learn gfapy
```

<a name="sec2.2"></a>
## Manual Install

Manually install dependencies: 
- [minimap2](https://github.com/lh3/minimap2)  

And python modules:
- [graph-tool](https://graph-tool.skewed.de)
- [numpy](https://numpy.org)
- [scikit-learn](https://scikit-learn.org/stable/install.html)
- [pandas](https://pandas.pydata.org/docs/getting_started/install.html)

<a name="sec3"></a>
# Running vsAware

vsAware natively support assembly result from SPAdes (includes metaSPAdes, metaviralSPAdes).

<a name="sec3.1"></a>
## Quick Usage

```bash
usage: vsAware [-h] -a {spades} -g GFA_FILE [-p PATH_FILE] [-mc MIN_COV] [-ml MIN_LEN] [-r REF_FILE] [-o OUTPUT_DIR]

Construct full-length viral strains under deno vo approach from contigs and assembly graph, currently supports SPAdes

optional arguments:
  -h, --help            show this help message and exit
  -a {spades}, --assembler {spades}
                        name of the assembler used. [spades]
  -g GFA_FILE, --graph GFA_FILE
                        path to the assembly graph, (.gfa format)
  -p PATH_FILE, --path PATH_FILE
                        contig file from SPAdes (.paths format), only required for SPAdes. e.g., contigs.paths
  -mc MIN_COV, --minimum_coverage MIN_COV
                        minimum node coverage cutoff, [default: auto]
  -ml MIN_LEN, --minimum_contig_length MIN_LEN
                        minimum initial contig length for strains [default: 250]
  -r REF_FILE, --reference_fa REF_FILE
                        path to the reference strain, .fasta format, DEBUG_MODE only
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        path to the output directory [default: acc/]
```

vsAware takes as input assembly graph in Graphical Fragment Assembly (GFA) Format and reported contig information. Both input can be found in the output directory after running SPAdes assembler. Do not modify any contig/node name from the assembly result for consistency. Please refer to [SPAdes](https://github.com/ablab/spades) for further guideline. Example usage as below:

```bash
# SPAdes assembler example, pair-end reads
python spades.py -1 forward.fa -2 reverse.fa --careful -t 16 -o output_dir
```

<a name="sec3.2"></a>
## Support SPAdes

For SPAdes, we recommend to use `--careful` option for more accurate assembly graph and contigs result. Please use `assembly_graph_after_simplification.gfa` and `contigs.paths` as input, and set `-a` flag to `spades`. Example usage as below:

```bash
python vsAware.py -a spades -g assembly_graph_after_simplification.gfa -p contigs.paths -o output_dir
```

<a name="sec3.3"></a>
## Parameters

### Minimum Node Coverage

This sets the minimum node coverage for filtering the inaccurate nodes from initial assembly graph. By default, the node coverage is automatically set to 5% of median node coverage among the graph by tool, which demonstrates good result. Please use `-mc` flag to input the customized minimum node coverage if needed.

### Minimum Contig Length

Since SPAdes usually output all the nodes from assembly graph as contigs, short output contig may lead to less accuracy and confidence. By default, contig with length less than 250bp (as 250bp is single read length from pair-end read) is filtered out. Please use `-ml` flag to input the customized minimum node coverage if needed.

<a name="sec4"></a>
# Citation

TODO

<a name="sec5"></a>
# Feedback and bug reports

Thanks for using vsAware, please feel free to raise any concern via `Issues`.