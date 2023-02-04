<p align="center">
  <img src="VStrains_logo.png" width="500" title="VStrains logo" alt="VStrains logo">
</p>

# VStrains: De Novo Reconstruction of Viral Strains via Iterative Path Extraction From Assembly Graphs

![GitHub](https://img.shields.io/github/license/metagentools/VStrains)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Manual
===========

Table of Contents
-----------------

1. [About VStrains](#sec1) </br>
2. [Updates](#sec2) </br>
3. [Installation](#sec3) </br>
   3.1. [Quick Install](#sec3.1) </br>
   3.1. [Manual Install](#sec3.2) </br>
4. [Running VStrains](#sec4) </br>
   4.1. [Quick Usage](#sec4.1) </br>
   4.2. [Support SPAdes](#sec4.2) </br>
   4.3. [Output](#sec4.3) </br>
5. [Stand-alone binaries](#sec5) </br>
6. [Experiment](#sec6) </br>
7. [Citation](#sec7) </br>
8. [Feedback and bug reports](#sec8)</br>

<a name="sec1"></a>
# About VStrains

VStrains is a de novo approach for reconstructing strains from viral quasispecies.

<!-- Please refer to our [paper](NULL) and [supplementary Material](NULL) for details methodology. -->

<a name="sec2"></a>
# Updates

## VStrains 1.1.0 Release (03 Feb 2023)
* Replace the PE link inference module `VStrains_Alignment.py` with `VStrains_PE_Inference.py`
   
   `VStrains_PE_Inference.py` implements a hash table approach that produce efficient perfect match lookup, the new module leads to consistent evaluation results and substantially decrease the runtime and memory usage against previous alignment approach.

<!-- * support direct install for Conda -->

<a name="sec3"></a>
# Installation

VStrains requires a 64-bit Linux system or Mac OS and python (supported versions are python3: 3.2 and higher).

<a name="sec3.1"></a>
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

<a name="sec3.2"></a>
## Manual Install

Manually install dependencies: 
- [minimap2](https://github.com/lh3/minimap2)  

And python modules:
- [graph-tool](https://graph-tool.skewed.de)
- [numpy](https://numpy.org)
- [gfapy](https://github.com/ggonnella/gfapy)
- [matplotlib](https://matplotlib.org)

<a name="sec4"></a>
# Running VStrains

VStrains supports assembly results from [SPAdes](https://github.com/ablab/spades) (includes metaSPAdes, metaviralSPAdes, etc).

<a name="sec4.1"></a>
## Quick Usage

```
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

<a name="sec4.2"></a>
## Support SPAdes

For SPAdes, we recommend to use `--careful` option for more accurate assembly results. Please use `assembly_graph_after_simplification.gfa` and `contigs.paths` as input, and set `-a` flag to `spades`. Example usage as below:

```bash
python src/vstrains.py -a spades -g assembly_graph_after_simplification.gfa -p contigs.paths -o output_dir -fwd forward.fastq -rve reverse.fastq
```

<a name="sec4.3"></a>
## Output


VStrains stores all output files in `<output_dir>` , which is set by the user.

* `<output_dir>/aln/` directory contains paired-end (PE) linkage information, which is stored in `pe_info` and `st_info`.
* `<output_dir>/gfa/` directory contains iteratively simplified assembly graphs, where `graph_L0.gfa` contains the assembly graph produced by SPAdes after Strandedness Canonization, `split_graph_final.gfa` contains the assembly graph after Graph Disentanglement, and `graph_S_final.gfa` contains the assembly graph after Contig-based Path Extraction, the rests are intermediate results. All the assembly graphs are in [GFA 1.0 format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md).
* `<output_dir>/paf/` and `<output_dir>/tmp/` are temporary directories, feel free to ignore them.
* `<output_dir>/strain.fasta` contains resulting strains in `.fasta`, the headers for each strain has the form `NODE_<strain name>_<sequence length>_<coverage>` which is compatiable to SPAdes contigs format.
* `<output_dir>/strain.paths` contains paths in the assembly graph (input `GFA_FILE`) corresponding to `strain.fasta` using [Bandage](https://github.com/rrwick/Bandage) for further downstream analysis.
* `<output_dir>/vstrains.log` contains the VStrains log.
<!-- <a name="sec3.3"></a> -->
<!-- ## Parameters -->

<!-- ### Minimum Node Coverage

This sets the minimum node coverage for filtering the inaccurate nodes from initial assembly graph. By default, the node coverage is automatically set based on coverage distribution, which demonstrates good result among all tested datasets. Please use `-mc` flag to input the customized minimum node coverage if needed.

### Minimum Contig Length

Since SPAdes normally output all the nodes from assembly graph as contigs, short or low coverage contig may lead to less accuracy and confidence. By default, single node contig with length less than 250bp or coverage less then `--mc` (defined above) is filtered out. Please use `-ml` flag to input the customized minimum contig length if needed. -->

<a name="sec5"></a>
# Stand-alone binaries

`evals/quast_evaluation.py` is a wrapper script for strain-level experimental result analysis using [MetaQUAST](https://github.com/ablab/quast).

```
usage: quast_evaluation.py [-h] -quast QUAST [-cs FILES [FILES ...]] [-d IDIR] -ref REF_FILE -o OUTPUT_DIR

Use MetaQUAST to evaluate assembly result

options:
  -h, --help            show this help message and exit
  -quast QUAST, --path_to_quast QUAST
                        path to MetaQuast python script, version >= 5.2.0
  -cs FILES [FILES ...], --contig_files FILES [FILES ...]
                        contig files from different tools, separated by space
  -d IDIR, --contig_dir IDIR
                        contig files from different tools, stored in the directory, .fasta format
  -ref REF_FILE, --ref_file REF_FILE
                        ref file (single)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        output directory
```

<a name="sec6"></a>
# Experiment

VStrains is evaluated on both simulated and real datasets under default settings, and the source of the datasets can be found in the links listed below:
1. Simulated Dataset, can be found at [savage-benchmark](https://bitbucket.org/jbaaijens/savage-benchmarks/src/master/) (No preprocessing is required)
   - 6 Poliovirus (20,000x)
   - 10 HCV (20,000x)
   - 15 ZIKV (20,000x)
2. Real Dataset (please refer to [Supplementary Material](https://www.biorxiv.org/content/10.1101/2022.10.21.513181v3.supplementary-material) for preprocessing the real datasets)
   - 5 HIV labmix (20,000x) [SRR961514](https://www.ncbi.nlm.nih.gov/sra/?term=SRR961514), reference genome sequences are available at [5 HIV References](https://github.com/cbg-ethz/5-virus-mix/blob/master/data/REF.fasta)
   - 2 SARS-COV-2 (4,000x) [SRR18009684](https://www.ncbi.nlm.nih.gov/sra/?term=SRR18009684), [SRR18009686](https://www.ncbi.nlm.nih.gov/sra/?term=SRR18009686), pre-processed reads and individually assemble ground-truth reference sequences can be found at [2 SARS-COV-2 Dataset](https://github.com/RunpengLuo/sarscov2-4000x)

<a name="sec7"></a>
# Citation
VStrains has been accepted at [RECOMB 2023](http://recomb2023.bilkent.edu.tr/program.html) and preprint is available at [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.10.21.513181v3).

If you use VStrains in your work, please cite the following publications.

Runpeng Luo and Yu Lin, VStrains: De Novo Reconstruction of Viral Strains via Iterative Path Extraction From Assembly Graphs

<a name="sec8"></a>
# Feedback and bug reports

Thanks for using VStrains. If any bugs be experienced during execution, please re-run the program with additional `-d` flag and provide the `vstains.log` together with user cases via `Issues`
