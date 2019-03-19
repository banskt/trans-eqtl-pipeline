# Reads an expression file in as a dataframe, transposes it, and
# saves it as txt table (RDS commented out).  If a covariate file is present (as 3rd
# input argument), then it will be used to produce new expression data
# by making a linear model expression ~ covariate, and then pulling the
# residuals as the new expression data.
#
# The input expression file is expected to be tab-delimted, with people
# as columns and genes as row.

argv <- commandArgs(trailingOnly = TRUE)
expressionfile <- argv[1]
outfile <- argv[2]
# Presence of covariate file suggests to correct for PEER factors, etc.
covariatefile <- ifelse(length(argv) == 3, argv[3], NA)

if (is.na(covariatefile)) {
  stop("missing required covariate file")
}

expression <- read.table(expressionfile, stringsAsFactors = FALSE, header = TRUE, row.names = 1, check.names=F)
message(paste(sep=" ", "Read",nrow(expression), "genes and", ncol(expression), "samples"))


if (!is.na(covariatefile)) {
  # Correct expression data for covariates.
  covariate <- read.table(covariatefile, stringsAsFactors = FALSE, header = TRUE, row.names = 1, check.names=F)
  message(paste(sep=" ", "Read", nrow(covariate), "covariates for", ncol(covariate), "samples"))

  if ( all(colnames(covariate) == colnames(expression))) {
    # Transpose expression.
    expression <- t(expression)
    covariate <- t(covariate)

    for (i in 1:length(colnames(expression))) {
      fit <- lm(expression[,i] ~ covariate)
      expression[,i] <- fit$residuals
    }
  } else {
    stop("Samples are not ordered accross covariates and expression")
  }
}

# saveRDS(expression, RDSout)
lm_expr = as.data.frame(t(expression))
lm_expr["gene_id"] <- rownames(lm_expr)                     # add gene_id as a column
lm_expr = lm_expr[,c(ncol(lm_expr), 1:(ncol(lm_expr)-1))]   # re-arrange columns
write.table(lm_expr, file=outfile, row.names = F, col.names = T, quote = F, sep="\t")