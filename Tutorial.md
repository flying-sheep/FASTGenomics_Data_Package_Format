# How to create a FASTGenomics Data Package

## Prerequisites
To run this tutorial, you need to have R installed, which can be obtained from 
the [R project website](https://www.r-project.org/). Furthermore, you need 
several add-on packages that help to conveniently create files in specific 
formats. If you have not yet installed the following packages, open R and
execute the following commands:

```
source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg,db")
install.packages("Matrix", dependencies = TRUE)
install.packages("yaml", dependencies = TRUE)
```

## Loading the data

First of all, we need to read the experiment data into R. Here, we're using the
example file available on GitHub, but you can also use your own data, but then
make sure that you set the working directory using `setwd`.

```
# load the data, returns a data.frame object
dense_matrix <- read.table("https://raw.githubusercontent.com/FASTGenomics/FASTGenomics_Data_Package_Format/master/example/dense_matrix_HGNC.txt",
                           header = TRUE,
		                   sep = "\t")
# convert data.frame to matrix
dense_matrix <- as.matrix(dense_matrix)

# look at the expression data
dense_matrix
```

Furthermore, there may be data available about the experiment, that we also
need to load into R.

```
# load the metadata table
dense_meta <- read.table("https://raw.githubusercontent.com/FASTGenomics/FASTGenomics_Data_Package_Format/master/example/dense_matrix_metadata.txt",
                         header = TRUE,
		                 sep = "\t")
						 
# look at the metadata
dense_meta
```

## Preparing cell metadata
The FASTGenomics platform handles cells not with their IDs, but with numbers
starting from 0, therefore we need to provided a mapping. As you see in our
example, the columns in the expression matrix are not in the same order as in
the metadata table, so we first adjust the ordering.

```
# sort the metadata table
index <- match(colnames(dense_matrix), dense_meta$sampleName)
dense_meta <- dense_meta[index, ]

# look at the metadata
dense_meta
```

In your own data, ensure that metadata is available for all cells. Now prepare
a FASTGenomics metadata table and adapt the expression data to it.

```
# make FASTGenomics metadata table
FG_cell_meta <- data.frame(cellId = 0:(ncol(dense_matrix)-1),
                           dense_meta,
						   stringsAsFactors = FALSE)

# assign FASTGenomics-compatible column names
colnames(dense_matrix) <- FG_cell_meta$cellId
```

## Initial quality control
You may want to apply some initial quality filters to your dataset and for
instance remove cells that do not express a certain number of genes and remove
genes that are not expressed in a minimum number of cells.

```
# count the genes per cell that have at least one read
genes_per_cell <- apply(dense_matrix, 2, function(x) sum(x>0))

# find cells that express at least one gene
cell_index <- genes_per_cell>0

# for documentation, save data from excluded cells
excluded_cell_dense_matrix <- dense_matrix[, !cell_index]
excluded_cell_dense_meta <- FG_cell_meta[!cell_index, ]

# remove the excluded cells from dataset
dense_matrix <- dense_matrix[, cell_index]
FG_cell_meta <- FG_cell_meta[cell_index, ]
colnames(FG_cell_meta) <- c("cellId*Integer", "sampleName*String", "cellType*String", "batch*String", "quality*String")

# count the number of cells, in which a gene is expressed (i.e. has at least one read)
cells_per_gene <- apply(dense_matrix, 1, function(x) sum(x>0))

# find genes that are expressed in at least one cell
gene_index <- cells_per_gene>0

# save data from excluded genes
excluded_gene_dense_matrix <- dense_matrix[!gene_index, ]

# remove excluded genes from dataset
dense_matrix <- dense_matrix[gene_index, ]
```

## ID mapping
The FASTGenomics platform works with Entrez IDs as gene identifiers. Our
example has however HGNC symbols, i.e. we need to platform an ID mapping to
translate the ID types from one into another. The following example should
make you aware of some pitfalls that may occur during that process.

```
# load the library providing the ID mapping
library(org.Hs.eg.db)

# translate the HGNC symbols to Entrez IDs
entrez_IDs <- mapIds(org.Hs.eg.db,
                     keys = rownames(dense_matrix),
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "list")
					 
# find HGNC symbols that were mapped to more than one Entrez ID
multiMappingIndex <- sapply(entrez_IDs, length)>1

# find HGNC symbols that were not mapped to any Entrez ID
unassignedIndex <- sapply(sapply(entrez_IDs, na.omit), length)==0

# combine names of multimapping Entrez IDs into one
entrez_IDs <- sapply(entrez_IDs , function(x) paste(x, collapse = " // "))

# check if some Entrez IDs are duplicated
duplicatedTargetIDs <- unique(entrez_IDs[duplicated(entrez_IDs[!unassignedIndex])])

# ... and find those that are concerned
duplicatedTargetIndex <- if (length(duplicatedTargetIDs)>0) {
  entrez_IDs %in% duplicatedTargetIDs
} else {
  rep(FALSE, length(entrez_IDs))
}

# create a log for the ID translation
mappingLog <- data.frame(entrezID = entrez_IDs,
						 originalID = rownames(dense_matrix),
                         log = "successfully mapped to unique Entrez ID",
						 stringsAsFactors= FALSE)
mappingLog$log[multiMappingIndex] <- "mapped to multiple Entrez IDs"
mappingLog$log[unassignedIndex] <- "no Entrez ID defined"
mappingLog$log[duplicatedTargetIndex] <- "multiple original IDs mapped to same Entrez ID"

# find genes to keep in dataset
entrezIndex <- !multiMappingIndex & !unassignedIndex & !duplicatedTargetIndex

# make Entrez ID-based dense expression matrix
dense_entrez <- dense_matrix[entrezIndex, ]
rownames(dense_entrez) <- entrez_IDs[entrezIndex]

# create gene metadata table
FG_gene_meta <- mappingLog[entrezIndex,  ]
colnames(FG_gene_meta) <- c("entrezID*Integer", "originalID*String", "log*String")
```

## Conversion of dense expression matrix to FASTGenomics format
Since single-cell RNA-seq data is sparse, a lot of memory (RAM during analyis
and disk space for the data files) can be saved when data is not handled in
dense, but sparse format. Here we show you how to do this.

```
# load the library to provide sparse matrix functionality
library(Matrix)

# convert dense to sparse matrix
sparse_entrez <- Matrix(dense_entrez,
                        sparse = TRUE)

# take a look
sparse_entrez
```

As you can see, all the zeros have disappeared from the matrix. Now let's
transform this to the FASTGenomics data format.

```
# make a temporary object
tmp <- summary(sparse_entrez)

# make FASTGenomics sparse expression matrix
FG_expression <- data.frame(cellId = colnames(sparse_entrez)[tmp$j],
                            entrezID = rownames(sparse_entrez)[tmp$i],
						    expressionValue = tmp$x)
colnames(FG_expression) <- c("cellId*Integer", "entrezID*Integer", "expressionValue*Number")
```

## Setting up the dataset manifest

```
# load yaml file parser
library(yaml)

# if you set up the manifest, fill in as much information as possible
FG_manifest <- list(schema_version = 2.1,
                    data = list(cell_metadata = list(file = "cell_metadata.tsv",
                                                     organism = 9606,
                                                     batch_column = "batch"),
                                gene_metadata = list(file = "gene_metadata.tsv"),
                                expression_data = list(file = "expression_data.tsv"),
                                supplemental = list(unconsidered_genes = list(expression_data = list(file = NULL),
                                                                              gene_metadata = list(file = NULL)))),
                    metadata = list(title = "FASTGenomics Data Package Example",
                                    technology = "in silico-generated data",
                                    version = 1,
                                    contact = "Comma Soft AG",
                                    description = "in silico-generated dataset for demonstrating how to create a FASTGenomics Data Package; available from https://github.com/FASTGenomics/FASTGenomics_Data_Package_Format",
                                    short_description = "In Silico Example Data Package",
                                    preprocessing = list(notes = "expression defined as at least one read present; cells with 0 expressed genes and genes expressed in 0 cells were removed; HGNC symbols were translated to Entrez IDs using org.Hs.eg.db version 3.4.2",
                                                      tools = "R with Matrix and yaml packages as well as the Bioconductor package org.Hs.eg.db."),
                                    image = "image.png"))
```

Finally, all the data needs to be saved according to the specifications made in
the FASTGenomics manifest.

```
write.table(FG_expression, 
            "expression_data.tsv",
			row.names = FALSE,
			col.names = TRUE,
			quote = FALSE,
			sep = "\t")

write.table(FG_cell_meta,
            "cell_metadata.tsv",
			row.names = FALSE,
			col.names = TRUE,
			quote = FALSE,
			sep = "\t")

write.table(FG_gene_meta,
            "gene_metadata.tsv",
			row.names = FALSE,
			col.names = TRUE,
			quote = FALSE,
			sep = "\t")

write(as.yaml(FG_manifest),
      "manifest.yml")
```

Now compress and assemble these files into one ZIP file and you are done.
