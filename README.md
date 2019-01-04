# Comparison of different methods in trans-eQTL

(Currently in development)

The following methods will be included:
* [x] [MatrixEQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/)
* [ ] [GNetLMM](https://github.com/PMBio/GNetLMM) (see [Installation instructions](https://github.com/banskt/trans-eqtl-pipeline/wiki/Install-GNetLMM-in-GWDG-cluster))
* [x] [CPMA](https://github.com/cotsapaslab/CPMAtranseqtl) (The authors did not provide updated software, the pipeline uses JPA scores from TEJAAS)
* [x] TEJAAS

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
* Expression (tab-separated text file: genes in rows, samples in columns. Header row with sample-ids, first column with gene names. In the header row, the first column is named `gene_id`). 
* Sample (a dummy [sample file in Oxford format](http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html))
* GENCODE file
* gene position file (for MatrixEQTL)
* MAF file from 1000Genomes

## Required softwares
* Python >3.6 (numpy, mpmath)
* TEJAAS
* LDSTORE
* GNeTLMM
* R v3.4.1 (MatrixEQTL)

## How to run
1. Within `bsubfiles` folder, change the job submission criteria and module loadings as per your requirements (GWDG users, skip this step)
2. Modify `main/utils/submit_job` to your own job scheduling mechanism (`bsub` users, skip this step)
3. Update the path of external programs `main/EXTERNAL` 
4. Update the path of your datasets in `main/DATA`.
5. Create a `CONFIG` file (see example in `configs/CONFIG`).
6. Run the pipeline from within `main` directory.
```
cd main
./01_validation_pipeline.sh configs/CONFIG
```
