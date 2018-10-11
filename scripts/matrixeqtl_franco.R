
# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Script for benchmarking by F. Simonetti



library("optparse")
# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)
 
option_list = list(
    make_option(c("-g", "--genotype"), type="character", default=NULL, 
              help="path to genotype file name", metavar="character"),
    make_option(c("-c", "--chrom"), type="numeric", default=NULL, 
              help="Chromosome number", metavar="number"),
    make_option(c("-d", "--donors"), type="character", default=NULL, 
              help="path to donor ids file name", metavar="character"),
    make_option(c("-u", "--selectdonors"), type="character", default=NULL, 
              help="path to selected donor ids file", metavar="character"),
    make_option(c("-i", "--geneinfo"), type="character", default=NULL, 
              help="path to gene info file name", metavar="character"),
    make_option(c("-s", "--dataset"), type="character", default=NULL, 
              help="dataset to process [gtex, cardiogenics]", metavar="character"),
    make_option(c("-x", "--expression"), type="character", default=NULL, 
              help="path to gene expression file name", metavar="character"),
    make_option(c("-p", "--pvalcis"), type="numeric", default=1e-3, 
              help="pvalue threshold for cis eQTLs [default= %default]", metavar="number"),
    make_option(c("-t", "--pvaltrans"), type="numeric", default=1e-3, 
              help="pvalue threshold for cis eQTLs [default= %default]", metavar="number"),
    make_option(c("-o", "--outdir"), type="character", default=".", 
              help="path to output directory [default \"%default\"]", metavar="character"),
    make_option(c("-m", "--model"), type="character", default="modelLINEAR", 
              help="Model to use from modelANOVA, modelLINEAR, or modelLINEAR_CROSS [default \"%default\"]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# opt$dataset    ="gtex"
# opt$genotype   ="/cbscratch/franco/datasets/gtex/prefiltered/GTEx_450Indiv_filtered_chr22.gz"
# opt$donors     ="/cbscratch/franco/datasets/gtex/donor_ids.fam"
# opt$expression ="/cbscratch/franco/datasets/gtex/Whole_Blood_Analysis.v6p.normalized.expression.txt"
# Gene expression samples: 338
# Gene expression genes: 23152
# Genotype SNPs: 72607




# opt$dataset    ="cardiogenics"
# opt$genotype   ="/cbscratch/franco/datasets/cardiogenics/genotypes/prefiltered/CG_dosages_filtered_22.imputed.gz"
# opt$donors     ="/cbscratch/franco/datasets/cardiogenics/genotypes/CG.sample"
# opt$expression ="/cbscratch/franco/datasets/cardiogenics/cardio_mono_expr.txt"
# opt$outdir     ="/cbscratch/franco/tejaas_output/tests"
# opt$geneinfo   ="/cbscratch/franco/datasets/gtex/genepos.gencode.v19.txt"
# opt$chrom      = 22
# opt$selectdonors="/usr/users/fsimone/Tejaas_Pipeline/devtools/mono_usersamples.txt"
# Gene expression samples: 744
# Gene expression genes: 15304
# Genotype SNPs: 71200


if (is.null(opt$genotype) | is.null(opt$expression) | is.null(opt$dataset) | is.null(opt$geneinfo)) {
  print_help(opt_parser)
  stop("Mandatory arguments missing.n", call.=FALSE)
}

if (opt$dataset!="gtex" & opt$dataset!="cardiogenics") {
    print_help(opt_parser)
    stop(paste("Available datasets are \"gtex\", \"cardiogenics\". Input is: ", opt$dataset, sep=""))
}

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
if (opt$model == "modelLINEAR")
    {useModel = modelLINEAR;}
if (opt$model == "modelLINEAR")
    {useModel = modelLINEAR;}
if (opt$model == "modelLINEAR")
    {useModel = modelLINEAR;
} else { stop("Invalid model selection, must be one of modelANOVA, modelLINEAR, or modelLINEAR_CROSS") }


if(is.null(opt$geneinfo)) {
    stop("Gene info file is missing");    
} else {
    genepos = read.table(opt$geneinfo, header = TRUE, stringsAsFactors = FALSE, colClasses = c("character", "character", "numeric", "numeric"));        
}

# Output file name
output_file_name_cis = paste(opt$outdir, "/",opt$dataset, "_MatrixEQTL_chr", opt$chrom, ".cisout", sep="");
output_file_name_tra = paste(opt$outdir, "/",opt$dataset, "_MatrixEQTL_chr", opt$chrom, ".transout", sep="");
dir.create(file.path(opt$outdir), showWarnings = FALSE, recursive=TRUE)

## Location of the package with the data files.
# base.dir = find.package('MatrixEQTL');


# Genotype file name
SNP_file_name = opt$genotype;
donors_file_name = opt$donors;

## Load genotype data
# crop matrix for genotypes only
readGTEX<-function(SNP_file_name, donors_file_name, maf_filter=T, maf_thres=0.1) {
    message("Reading Genotype ...")
    message(SNP_file_name)
    snps_mat = read.csv(file=SNP_file_name, sep=" ", stringsAsFactors=F, header=F)
    row_names = snps_mat[,2]
    donor_ids = read.csv(file=donors_file_name, sep=" ", stringsAsFactors=F, header=F)[,1]

    # get SNPs positions for cis and trans analysis (before cropping the snp matrix)
    snpspos = snps_mat[,c(2,1,3)]
    snpspos[,2] = paste("chr", snpspos[,2], sep="")
    colnames(snpspos) = c("snpid","chr","pos")
    snpspos[,3] = as.numeric(snpspos[,3])

    mafs = snps_mat[,6]

    snps_mat = snps_mat[,7:ncol(snps_mat)]
    rownames(snps_mat) = row_names
    colnames(snps_mat) = donor_ids

    if (maf_filter) {
        mafs2remove = mafs < maf_thres | mafs > (1-maf_thres)
        
        mafs      = mafs[!mafs2remove]
        snps_mat  = snps_mat[!mafs2remove, ]
        snpspos   = snpspos[!mafs2remove,]
    }    

    snps_mat = as.matrix(snps_mat)
    return(list(snps_mat, snpspos))
}

readOxford<-function(SNP_file_name, donors_file_name, chrom, maf_filter=T, maf_thres=0.1) {
    message("Reading Genotype ...")
    message(SNP_file_name)
    snps_mat = read.csv(file=SNP_file_name, sep=" ", stringsAsFactors=F, header=F)
    row_names = snps_mat[,2]
    donor_ids = read.csv(file=donors_file_name, sep=" ", stringsAsFactors=F, header=F, comment.char="#")[,1]

    snpspos = snps_mat[,c(2,1,3)]
    snpspos[,2] = rep(paste("chr", chrom, sep=""), nrow(snpspos))
    colnames(snpspos) = c("snpid","chr","pos")
    snpspos[,3] = as.numeric(snpspos[,3])

    snps_mat_gt = snps_mat[,6:ncol(snps_mat)]
    AAindex = seq(1, ncol(snps_mat_gt),3)
    ABindex = seq(2, ncol(snps_mat_gt),3)
    BBindex = seq(3, ncol(snps_mat_gt),3)

    dosages_mat = snps_mat_gt[,ABindex] + 2*snps_mat_gt[,BBindex] 

    rownames(dosages_mat) = row_names
    colnames(dosages_mat) = donor_ids

    if (maf_filter) {
        mafs        = apply(dosages_mat,1,sum) / 2 / ncol(dosages_mat)
        mafs2remove = mafs < maf_thres | mafs > (1-maf_thres)
        
        mafs        = mafs[!mafs2remove]
        dosages_mat = dosages_mat[!mafs2remove, ]
        snpspos     = snpspos[!mafs2remove,]
    }

    dosages_mat = as.matrix(dosages_mat)
    return(list(dosages_mat, snpspos))
}

if (is.null(opt$chrom)) {
    warning("Chromosome number has not been set")
}



if (opt$dataset == "gtex"){
    res = readGTEX(SNP_file_name, donors_file_name)
}

if (opt$dataset == "cardiogenics") {

    if (is.null(opt$chrom)) {
        stop("Chromosome number required for reading Cardiogenics data")
    }

    #trim genepos ensembl ids
    ensembls = sapply(genepos$geneid, function(x){strsplit(x, '.', fixed=T)})
    ensembls = as.character(sapply(ensembls, "[[", 1))
    genepos$geneid = ensembls

    # read genotype
    # res = readOxford(SNP_file_name, donors_file_name, opt$chrom)    
    # after filtering the genotypes, it is all in GTEX format
    res = readGTEX(SNP_file_name, donors_file_name)
}

snps_mat = res[[1]] #genotype matrix
snpspos  = res[[2]] #SNP position info


snps = SlicedData$new();
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$CreateFromMatrix( snps_mat ) 


## Load gene expression data
expression_file_name = opt$expression;
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);


# match columns of samples from genotype and expression

if (! is.null(opt$selectdonors)) {
    message("Selecting donors from user list")
    selectdonors = read.csv(file=opt$selectdonors, sep=" ", stringsAsFactors=F, header=F, comment.char="#")
    col_index1 = match(selectdonors[[1]], colnames(gene))
    col_index2 = match(selectdonors[[1]], colnames(snps))
} else {
    common_donors = intersect(colnames(gene), colnames(snps))
    col_index1 = match(common_donors, colnames(gene))
    col_index2 = match(common_donors, colnames(snps))
}


if (all(colnames(snps)[col_index2] == colnames(gene)[col_index1])) {
    gene$ColumnSubsample(col_index1)
    snps$ColumnSubsample(col_index2)
}



# Covariates file name. Set to character() for no covariates
if (is.null(opt$covariates)) {
    cvrt = SlicedData$new();
} else {

    ## Load covariates
    cvrt = SlicedData$new();
    cvrt$fileDelimiter = "\t";      # the TAB character
    cvrt$fileOmitCharacters = "NA"; # denote missing values;
    cvrt$fileSkipRows = 1;          # one row of column labels
    cvrt$fileSkipColumns = 1;       # one column of row labels
    if(length(covariates_file_name)>0) {
        cvrt$LoadFile(opt$covariates);
    }
}



# Only associations significant at this level will be saved
pvOutputThreshold_cis = opt$pvalcis;
pvOutputThreshold_tra = opt$pvaltrans;

# Distance for local gene-SNP pairs
cisDist = 1e6;

# Error covariance matrix. Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

message(paste("Gene expression samples: ", ncol(gene), sep=""))
message(paste("Gene expression genes: ",   nrow(gene), sep=""))
message(paste("Genotype SNPs: ",   nrow(snps), sep=""))

## Run the analysis
me = Matrix_eQTL_main(
snps = snps,
gene = gene,
cvrt = cvrt, 
output_file_name = output_file_name_tra,
pvOutputThreshold = pvOutputThreshold_tra,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
useModel = useModel, 
errorCovariance = errorCovariance, 
verbose = TRUE,
cisDist = cisDist,
snpspos = snpspos, 
genepos = genepos,
pvalue.hist = TRUE,
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat(paste('Detected eQTLs:', me$all$neqtls, '\n', sep=""));


# show(me$all$eqtls)

## Plot the histogram of all p-values

# png(paste(opt$outdir, "/",opt$dataset, "_MatrixEQTL.pvalue.hist.png", sep=""))
# plot(me)
# dev.off()

# me$all
# $ntests
# [1] 3697304944
# $neqtls
# [1] 71141786
# $hist.bins
# $hist.counts
