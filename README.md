# The FASTGenomics Data Package Format

## Introduction

Single-cell RNA-seq datasets typically consist of several data types that are 
all required for the understanding of an experiment. The FASTGenomics 
ecosystem for single-cell RNA-seq analyses provides functionality to make 
these analyses as easy and convenient as possible. To enable the data analysis 
with [FASTGenomics](https://fastgenomics.org), the dataset must be provided in 
a defined format which is detailed below. Briefly, a dataset consists of files 
containing expression data, metadata describing cells and genes as well as the 
experimental conditions. To reduce disk space, all these files are bundled 
into one ZIP file.


### Components of a FASTGenomics Data Package

The following table gives an overview of files that have to be included in a 
FASTGenomics dataset package:

* `manifest.yml`	Detailed description of the dataset package, including file definitions and dataset description.
* `expression_data.tsv`	Expression data for Entrez ID-coded genes in sparse FASTGenomics format used for data analysis.
* `cell_metadata.tsv`	Tab-separated file with metadata about cells used in the analysis.
* `gene_metadata.tsv`	File containing the gene IDs used for analysis.


### Data Package Description

The data package description is supplied via the manifest.yml file, which has 
the following structure:

* schema_version: 2.1
* data:
  * cell_metadata:
    * file: <cell metadata file name, e.g. cell_metadata.tsv>
    * organism: <NCBI taxonomy ID, e.g. 10090 for mouse>
    * \[optional\] batch_column: <column name in cell metadata file>
  * gene_metadata:
    * file: <file name containing the gene IDs used for analysis, e.g. gene_metadata.tsv>
  * expression_data: 
    * file: <file name for the single-cell RNA-seq dataset in sparse FASTGenomics format, e.g. expression_data.tsv>
  * \[optional\] supplemental: 
    * unconsidered_genes:
      * expression_data:
        * file: <file containting data of unconsidered genes, e.g. unconsidered_expression_data.tsv>
      * gene_metadata:
        * file: <file containing gene IDs excluded from analysis, e.g. unconsidered_gene_metadata.tsv>
* metadata:
  * title: <title of the dataset, composed of last name of first author et al., journal name, followed by publication year in brackets, e.g. Mass et al., Science (2016)>
  * technology: <single-cell RNA-seq technology used to generate the dataset, e.g. MARS-seq>
  * version: <version of the dataset, starting from 1>
  * contact: <contact person and/or institution, e.g. Comma Soft AG>
  * description: <description of the datset, including experimental setting, link to public repository and ID of dataset, as well as link to publication, where applicable>
  * short_description: <short description of dataset, e.g. the experimental setting>
  * preprocessing:
    * notes: <description how the dataset has been prepared, including cell and gene exclusion criteria, etc.>
    * tools:
      * <tools used for dataset preparation, e.g. FGpackageR>
      * <another tool used for dataset preparation>
    * \[optional\] image: <optional file name for image shown in the FASTGenomics data store, e.g. image.png>


### File Format Specifications

#### Sparse Expression File Format

Single-cell RNA-seq data is typically zero-inflated, i.e. for many genes no 
gene expression values are available. In a dense gene expression matrix, where
each column represents a cell and each row represents a gene and each 
row/column combination holds the expression value for a particular gene in a 
particular cell, expression values of unexpressed and uncaptured genes are 
represented as zeros. Depending on the technology used to generate the 
single-cell expression dataset, the proportion of zeros in a dense matrix may 
be higher than 90%. To save disk space, FASTGenomics therefore uses a sparse 
matrix data format to store expression data. A FASTGenomics sparse expression 
matrix file is a simple text file storing data in three tab-separated columns with a header line. This case-sensitive header 
line stores mandatory column names (`cellId`, `entrezId` and `expressionValue`) 
and data type information (`Integer`, `Number` and `String`) separated by an asterisk 
(`*`). The first column contains zero-based integer values identifying cells. 
The second column contains the identifier representing a gene (the Entrez ID 
for genes analyzed in [FASTGenomics](https://fastgenomics.org)). The third column holds the expression 
value of this gene in a particular cell (see example below).

```
cellId*Integer	entrezId*Integer	expressionValue*Number 
0	12544	4.0
0	67608	1.0
0	12390	1.0
1	12544	5.0
1	67608	1.2
1	12390	3.3
2	12544	4.5
2	67608	10.0
2	12390	1.2
```

Gene expression values of cells examined in [FASTGenomics](https://fastgenomics.org) are stored in 
a file defined in the manifest describing the data package, e.g. 
`expression_data.tsv`. Expression data excluded from analysis (see below) can 
be included in the dataset to increase transparency. 


#### Gene IDs

The [FASTGenomics](https://fastgenomics.org) pipeline works with Entrez IDs for genes. Many 
single-cell expression datasets stored in public databases (e.g. [Gene Expression Omnibus](http://www.ncbi-nlm.nih.gov/geo)) 
encode genes with other ID types (e.g. ENSEMBL IDs or gene symbols). In this 
case, IDs must be mapped to Entrez IDs. Genes for which no Entrez ID or 
several Entrez IDs are available should be excluded from data analysis.


#### Gene Metadata

Gene metadata is supplied in tab-separated text files with a header line. Gene 
metadata must be defined for all genes used in the analysis; the respective 
file must have the name `entrezId*Integer` in its first column, further columns, 
e.g. containing mappings to gene symbols, etc. are optional. The optional file 
containing metadata on excluded genes lists their IDs (in the first column 
with appropriate type definition, e.g. `gene*String`) along with optional 
further information (e.g. reason for exclusion).

#### Cell Metadata

Cell metadata is provided in a tab-separated text file that contains a header line. 

#### Image

If an image is being supplied in the manifest, use a common-known image format such as png or jpg (e.g. dropseq.png).

### Package Bundling

All files to be included in the FASTGenomics data package have to be packed in
one ZIP file.
