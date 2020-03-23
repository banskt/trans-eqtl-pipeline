import pandas as pd
import argparse
import gzip, os

def parse_args():

    parser = argparse.ArgumentParser(description='Convert rpkm for a given tissue into inverse quantile normalized gene expression')

    parser.add_argument('--input',
                        type=str,
                        dest='input',
                        metavar='FILE',
                        help='input covariate file (tissue-specific)')

    parser.add_argument('--age',
                        type=str,
                        dest='age_file',
                        metavar='FILE',
                        help='age covariate file for all samples')

    parser.add_argument('--output',
                        type=str,
                        dest="outfile",
                        metavar='STR',
                        help='output file')


    opts = parser.parse_args()
    return opts

opts = parse_args()

age = pd.read_table(opts.age_file, header=None, index_col=0)
age.columns = ["AGE", "PMI"]

covars = pd.read_table(opts.input, header=0, index_col=0)

# Fix last id trailing space
ids = [x.rstrip() for x in covars.columns]
covars.columns = ids
ids_w_age = [x for x in ids if x in age.index]

if len(ids_w_age) != len(ids):
    ids_not_age = [x for x in ids if x not in age.index]
    print("Following subjects don't have AGE:", ids_not_age)
    raise Exception("No AGE found for some SUBJECTS")

covariate_age_table = pd.concat([covars, age.loc[ids].T]) #.drop("PMI")

'''
We don't want this for now
AGE2 = covariate_age_table.loc["AGE"]**2
AGE3 = covariate_age_table.loc["AGE"]**3

ages23 = pd.concat([AGE2, AGE3], axis=1)
ages23.columns = ["AGE2", "AGE3"]

# covariates_ages123_table = pd.concat([covariate_age_table, ages23.T])
'''
covariates_ages123_table = covariate_age_table
covariates_ages123_table.index.name = "ID"
# covariates_ages123_table = covariates_ages123_table.rename_axis("ID", axis=1)
covariates_ages123_table.to_csv(opts.outfile, sep="\t", header=True, index=True)