# CircosRNAfold
An python script tool for visualization of Circos plots with high-throughput sequencing data.

The goal of this project is to develop a simple, circular visualization tool of binding sites and their raw data on folded 2D RNA structures. Generating the plots allows the user to get an overview of their experimental data without the need to read in Circos or similar.


The pyCircos package (https://github.com/ponnhide/pyCircos) is used to create the diagram and Viennarna (https://github.com/ViennaRNA/ViennaRNA) is used to predict the fold.


# Installation

You can install this script directly from github with:
```
git clone https://github.com/mloskant/CircosRNAfold.git
```

## Prerequisites
You can install the dependencies for using the script via the ```requirements.text``` file by running the following code:
```
pip install -r requirements.txt
```
In addtition, Python 3.10 or higher is required.


# Usage
The script can be started with various parameters. 
The required parameters are: 

- `-g` or `--genes` : The genes in BED format
- `-fi` or `--fasta` : The chromosomal sequences as a fasta file

Optional parameters can also be used:
- `-b` or `--bindings` : The binding sites of proteins in BED format (max. 4 files)
- `-bw` : BigWig file(s) with the raw data (max. 4 files)
- `-m` : The mature miRNAs in BED format
- `-color` : The colors for the iCLIP plots. If none or too few/too many colors are specified than there are specified bw files, default colors are used (max. 4 colors)
- `-o` or `--output` : Name of the folder in which all pngs are saved

## Running CircosRNAfold with sample data: 
You must first <em>sample_data/example_data_TAIR10_chr_all.fas.gz</em> unpack and then:

```
python3 circosrnafold.py --genes sample_data/example_data_genes.bed --bindings sample_data/example_data_binding_sites.bed --fasta sample_data/example_data_TAIR10_chr_all.fas -bw sample_data/example_data_bw_minus.bw -m sample_data/example_data_mature_rnas.bed --output circoplots
```

# Setting raw data color(s)
A maximum of four colors can be specified for the maximum of four bigWig files. There are the following four default colors: #800080 (purple), #006400 (green), #FF6347 (red), #4682B4 (blue). Depending on the number of bigWig files specified, the appropriate number of default colors are used in sequence as listed above. The default colors are used if no colors or too many/too few are specified by you. If this is the case, a message appears on the terminal. To specify the colors yourself, use the `-color` parameter followed by the colors in the hexadecimal system (without a preceding #).

# Gallery


# Input file format overview
This overview clearly sets out all the input file types and their requirements.


| Input | Description |
|--------------------|--------------------------------------------------------|
|Genome data|FASTA file conforming to the standard|
|Genes|BED7+ files conforming to the standard|
|iCLIP raw data| bigWig files conforming to the standard|
|Binding Site data|BED7+ files conforming to the standard|
|Mature RNAs|BED7+ files conforming to the standard|

# License
MIT License
