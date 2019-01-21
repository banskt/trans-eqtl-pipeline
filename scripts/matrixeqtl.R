# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Script for benchmarking by F. Simonetti
# Minor modifications by Saikat Banerjee to include file in bash pipeline, remove input of chromosome

library("optparse")
library(MatrixEQTL)

option_list = list(
    make_option(c("-g", "--genotype"), type="character", default=NULL,
              help="path to genotype file name", metavar="character"),
    make_option(c("-d", "--donors"), type="character", default=NULL,
              help="path to donor ids file name", metavar="character"),
    make_option(c("-u", "--selectdonors"), type="character", default=NULL,
              help="path to selected donor ids file", metavar="character"),
    make_option(c("-i", "--geneinfo"), type="character", default=NULL,
              help="path to gene info file name", metavar="character"),
    make_option(c("-s", "--datatype"), type="character", default=NULL,
              help="type of dataset [gtex, cardiogenics]", metavar="character"),
    make_option(c("-e", "--expression"), type="character", default=NULL,
              help="path to gene expression file name", metavar="character"),
    make_option(c("-p", "--pvalcis"), type="numeric", default=1e-3,
              help="pvalue threshold for cis eQTLs [default= %default]", metavar="number"),
    make_option(c("-t", "--pvaltrans"), type="numeric", default=1e-3,
              help="pvalue threshold for cis eQTLs [default= %default]", metavar="number"),
    make_option(c("-o", "--outfilecis"), type="character", default=NULL,
              help="path to output file for cis-eQTLs", metavar="character"),
    make_option(c("-q", "--outfiletrans"), type="character", default=NULL,
              help="path to output file for trans-eQTLs", metavar="character"),
    make_option(c("-m", "--model"), type="character", default="modelLINEAR",
              help="Model to use from modelANOVA, modelLINEAR, or modelLINEAR_CROSS [default \"%default\"]", metavar="character"),
    make_option(c("-r", "--randomize"), action="store_true", default=FALSE,
              help="Randomize the gene expression"),
    make_option(c("-R", "--shufflewith"), type="character", default=NULL,
              help="file with shuffled donor ids for genotype shuffling", metavar="character")
);

opt_parser = OptionParser(option_list=option_list, add_help_option = TRUE);
opt = parse_args(opt_parser);

read_genotype <- function(SNP_file_name, donors_file_name, maf_filter=T, maf_thres=0.1) {
    message("Reading Genotype ...")
    message(SNP_file_name)
    snps_mat = read.csv(file=SNP_file_name, sep=" ", stringsAsFactors=F, header=F)
    row_names = snps_mat[,2]
    donor_ids = read.csv(file=donors_file_name, sep=" ", stringsAsFactors=F, header=F, skip=2)[,1]

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

if(is.null(opt$geneinfo)) {
    stop("Gene info file is missing");
} else {
    genepos = read.table(opt$geneinfo, header = TRUE, stringsAsFactors = FALSE, colClasses = c("character", "character", "numeric", "numeric"));
}

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Output file name
output_file_name_cis = opt$outfilecis;
output_file_name_tra = opt$outfiletrans;

SNP_file_name = opt$genotype;
donors_file_name = opt$donors;

res = read_genotype(SNP_file_name, donors_file_name)
snps_mat = res[[1]] #genotype matrix
snpspos  = res[[2]] #SNP position info

if (!is.null(opt$shufflewith)) {
    message("Shuffling genotype using supplied donor IDs");
    shuffled_colnames = scan(opt$shufflewith, what="", sep="\n")
    colnames(snps_mat) = shuffled_colnames;
}

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

if (opt$randomize) {
    message("Randomizing gene expression values by permuting samples")
    expr_mat = as.matrix(gene);
    shuffled = t(apply(expr_mat, 1, sample));
    colnames(shuffled) = colnames(gene);
    gene$CreateFromMatrix(shuffled);
}

# match columns of samples from genotype and expression
if (!is.null(opt$selectdonors)) {
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

if (opt$datatype == "cardiogenics") {
    #trim genepos ensembl ids
    ensembls = sapply(genepos$geneid, function(x){strsplit(x, '.', fixed=T)})
    ensembls = as.character(sapply(ensembls, "[[", 1))
    genepos$geneid = ensembls
}


# Only associations significant at this level will be saved
pvOutputThreshold_cis = opt$pvalcis;
pvOutputThreshold_tra = opt$pvaltrans;

# Distance for local gene-SNP pairs
cisDist = 1e6;

# Error covariance matrix. Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

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
