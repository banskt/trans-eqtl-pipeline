# Comparison of different methods in trans-eQTL

(Currently in development)

The following methods will be included:
* [x] MatrixEQTL
* [ ] GNetLMM
* [ ] CPMA (implemented as JPA within TEJAAS)
* [ ] TEJAAS

We also want to compare:
* [ ] effect of different pre-filtering methods
* [ ] kNN
* [ ] effect of sparsity in TEJAAS

And, finally we plot everything together:
* [ ] Plot

## Method
We use the gene expression of two different tissues within the same population.
We find trans-eQTLs using different methods,
and then compare the methods using tissue-consistent trans-eQTLs (which are found in both tissues).

## Input
The pipeline expects the following input files:
* Genotype (in gzipped dosage format)
* Expression (tab-separated text file, gene name in column1, expression for `N` patients in the next `N` columns, header line starting with `gene_id` in first column and sample ids in the next `N` columns)
* Sample (a dummy [sample file in Oxford format](http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html))

## How to run
1. Update the file paths in `main/PATH`.
2. Create a `CONFIG` file (see example in `configs/CONFIG`).
3. Run the different scripts from within `main` directory.
