# VStrains: De Novo Reconstruction of Viral Strains via Iterative Path Extraction From Assembly Graphs

![GitHub](https://img.shields.io/github/license/metagentools/VStrains)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Manual
===========

Table of Contents
-----------------

1. [About VStrains](#sec1) </br>
2. [Installation](#sec2) </br>
   2.1. [Quick Install](#sec2.1) </br>
   2.1. [Manual Install](#sec2.2) </br>
3. [Running VStrains](#sec3) </br>
   3.1. [Quick Usage](#sec3.1) </br>
   3.2. [Support SPAdes](#sec3.2) </br>
4. [Experiment](#sec4) </br>
5. [Citation](#sec5) </br>
6. [Feedback and bug reports](#sec6)</br>
<a name="sec1"></a>
# About VStrains

VStrains is a de novo approach for reconstructing strains from viral quasispecies.

<!-- Please refer to our [paper](NULL) and [supplementary Material](NULL) for details methodology. -->

<a name="sec2"></a>
# Installation

VStrains requires a 64-bit Linux system or Mac OS and python (supported versions are python3: 3.2 and higher), 

<a name="sec2.1"></a>
## Quick Install (**recommended**)

Install [(mini)conda](https://conda.io/miniconda.html) as a light-weighted package management tool. Run the following commands to initialize and setup the conda environment for VStrains

```bash
# add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment
conda create --name VStrains-env

# activate conda environment
conda activate VStrains-env

conda install -c bioconda -c conda-forge python=3 graph-tool minimap2 numpy gfapy matplotlib
```

<a name="sec2.2"></a>
## Manual Install

Manually install dependencies: 
- [minimap2](https://github.com/lh3/minimap2)  

And python modules:
- [graph-tool](https://graph-tool.skewed.de)
- [numpy](https://numpy.org)
- [gfapy](https://github.com/ggonnella/gfapy)
- [matplotlib](https://matplotlib.org)

<a name="sec3"></a>
# Running VStrains

VStrains supports assembly results from [SPAdes](https://github.com/ablab/spades) (includes metaSPAdes, metaviralSPAdes, etc).

<a name="sec3.1"></a>
## Quick Usage

```bash
usage: VStrains [-h] -a {spades} -g GFA_FILE [-p PATH_FILE] [-o OUTPUT_DIR] -fwd FWD -rve RVE

Construct full-length viral strains under de novo approach from contigs and assembly graph, currently supports
SPAdes

optional arguments:
  -h, --help            show this help message and exit
  -a {spades}, --assembler {spades}
                        name of the assembler used. [spades]
  -g GFA_FILE, --graph GFA_FILE
                        path to the assembly graph, (.gfa format)
  -p PATH_FILE, --path PATH_FILE
                        contig file from SPAdes (.paths format), only required for SPAdes. e.g., contigs.paths
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        path to the output directory [default: acc/]
  -fwd FWD, --fwd_file FWD
                        paired-end sequencing reads, forward strand (.fastq format)
  -rve RVE, --rve_file RVE
                        paired-end sequencing reads, reverse strand (.fastq format)
```

VStrains takes as input an assembly graph in Graphical Fragment Assembly (GFA) Format and associated contig information. Both inputs can be found in the output directory after running SPAdes assembler. Please also provide the raw reads in paired-end format (e.g., forward.fastq, reverse.fastq) together as inputs. Do not modify any contig/node name from the SPAdes assembly results for consistency. Please refer to [SPAdes](https://github.com/ablab/spades) for further guideline. Example usage as below:

```bash
# SPAdes assembler example, pair-end reads
python spades.py -1 forward.fastq -2 reverse.fastq --careful -t 16 -o output_dir
```

<a name="sec3.2"></a>
## Support SPAdes

For SPAdes, we recommend to use `--careful` option for more accurate assembly results. Please use `assembly_graph_after_simplification.gfa` and `contigs.paths` as input, and set `-a` flag to `spades`. Example usage as below:

```bash
python src/vstrains.py -a spades -g assembly_graph_after_simplification.gfa -p contigs.paths -o output_dir -fwd forward.fastq -rve reverse.fastq
```

<!-- <a name="sec3.3"></a> -->
<!-- ## Parameters -->

<!-- ### Minimum Node Coverage

This sets the minimum node coverage for filtering the inaccurate nodes from initial assembly graph. By default, the node coverage is automatically set based on coverage distribution, which demonstrates good result among all tested datasets. Please use `-mc` flag to input the customized minimum node coverage if needed.

### Minimum Contig Length

Since SPAdes normally output all the nodes from assembly graph as contigs, short or low coverage contig may lead to less accuracy and confidence. By default, single node contig with length less than 250bp or coverage less then `--mc` (defined above) is filtered out. Please use `-ml` flag to input the customized minimum contig length if needed. -->

<a name="sec4"></a>
# Experiment

VStrains is evaluated on both simulated and real datasets under default settings, and the source of the datasets can be found in the links listed below:
1. Simulated Dataset, can be found at [savage-benchmark](https://bitbucket.org/jbaaijens/savage-benchmarks/src/master/) (No preprocessing is required)
   - 6 Poliovirus (20,000x)
   - 10 HCV (20,000x)
   - 15 ZIKV (20,000x)
2. Real Dataset (please refer to [Supplementary Material](https://www.biorxiv.org/content/10.1101/2022.10.21.513181v3.supplementary-material) for preprocessing the real datasets)
   - 5 HIV labmix (20,000x) [SRR961514](https://www.ncbi.nlm.nih.gov/sra/?term=SRR961514), reference genome sequences are available at [link](https://github.com/cbg-ethz/5-virus-mix/blob/master/data/REF.fasta)
   - 2 SARS-COV-2 (4,000x) [SRR18009684](https://www.ncbi.nlm.nih.gov/sra/?term=SRR18009684), [SRR18009686](https://www.ncbi.nlm.nih.gov/sra/?term=SRR18009686), pre-processed reads and individually assemble ground-truth reference sequences can be found at [link](https://github.com/RunpengLuo/sarscov2-4000x)

<a name="sec5"></a>
# Citation
VStrains has been accepted at [RECOMB 2023](http://recomb2023.bilkent.edu.tr/program.html) and preprint is available at [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.10.21.513181v3).

If you use VStrains in your work, please cite the following publications.

Runpeng Luo and Yu Lin, VStrains: De Novo Reconstruction of Viral Strains via Iterative Path Extraction From Assembly Graphs

<a name="sec6"></a>
# Feedback and bug reports

Thanks for using VStrains. Please feel free to provide any feedback or raise any concern via `Issues`.
