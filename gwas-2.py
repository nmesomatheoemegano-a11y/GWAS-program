
import argparse
import pysam
import pandas as pd
import numpy as np
import statsmodels.api as sm
import gzip

# Your code goes here
def run_gwas(bcf_path, pheno_path, covar_path, output_path):
  print(f"Starting GWAS analysis....")
  print(f"Input BCF file: {bcf_path}")
  print(f"Input phenotype file: {pheno_path}")
  print(f"Input covariates file: {covar_path}")
  print(f"Output file: {output_path}")

#opening the genotype data from BCF using Pysam
  ordered_bcf_samples = None
  with pysam.VariantFile(bcf_path, "r") as ifile:
    ordered_bcf_samples = list(ifile.header.samples)
    n_samples = len(ordered_bcf_samples)

  #sanity checks
  assert n_samples == 2504, "Number of samples is wrong."

  #loading phenotype and covariates data from text file
  try:
     df_hdlphenotypes = pd.read_csv(pheno_path, sep = "\t")
     df_covariates = pd.read_csv(covar_path, sep = "\t")
  except FileNotFoundError as e:
      print(f"Error loading files: {e}")
      return

  #making sure the order of the samples in phenotype and covariate data frame correspond with those in BCF files.
  df_hdlphenotypes = df_hdlphenotypes.set_index('IID').reindex(ordered_bcf_samples).reset_index()
  df_covariates = df_covariates.set_index('IID').reindex(ordered_bcf_samples).reset_index()

  #sanity checks
  assert len(df_hdlphenotypes) == n_samples, "Number of samples in hdl phenotype data frame is wrong."
  assert len(df_covariates) == n_samples, "Number of samples in covariates data frame is wrong."
  assert all(x == y for x, y in zip(df_hdlphenotypes ['IID'].values, ordered_bcf_samples)), "Sample order in BCF file and Hdl phenotype data frame do not match."
  assert all(x == y for x, y in zip(df_covariates ['IID'].values, ordered_bcf_samples)), "Sample order in BCF file and covariates data frame do not match."

  #combining all variables and covariates into a single data frame
  df_data = pd.DataFrame({
      "HDL" : df_hdlphenotypes['HDL'],
      "genotype": np.full(n_samples, np.nan, dtype = np.float32),
      "sex": df_covariates['SEX'],
      "age": df_covariates['AGE'],
      "PC1": df_covariates['PC1'],
      "PC2": df_covariates['PC2'],
      "PC3": df_covariates['PC3'],
      "PC4": df_covariates['PC4'],
      "PC5": df_covariates['PC5'],
      "PC6": df_covariates['PC6'],
      "PC7": df_covariates['PC7'],
      "PC8": df_covariates['PC8'],
      "PC9": df_covariates['PC9'],
      "PC10": df_covariates['PC10']
      })

  #sanity check
  assert len(df_data) == len(df_hdlphenotypes), "Length of the final data frame does not match the original data"
  assert df_data.columns.tolist() == ['HDL', 'genotype', 'sex', 'age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'], "Wrong final data frame created"
  assert df_data["HDL"].dtype == "float64", "HDL column has wrong data type"
  assert df_data["genotype"].dtype == "float32", "genotype column has wrong data type."
  assert df_data["sex"].dtype == "int64", "sex column has wrong data type."
  assert df_data["age"].dtype == "int64", "age column has wrong data type."
  assert df_data["PC1"].dtype == "float64", "PC1 column has wrong data type."
  assert df_data["PC2"].dtype == "float64", "PC2 column has wrong data type."
  assert df_data["PC3"].dtype == "float64", "PC3 column has wrong data type."
  assert df_data["PC4"].dtype == "float64", "PC4 column has wrong data type."
  assert df_data["PC5"].dtype == "float64", "PC5 column has wrong data type."
  assert df_data["PC6"].dtype == "float64", "PC6 column has wrong data type."
  assert df_data["PC7"].dtype == "float64", "PC7 column has wrong data type."
  assert df_data["PC8"].dtype == "float64", "PC8 column has wrong data type."
  assert df_data["PC9"].dtype == "float64", "PC9 column has wrong data type."
  assert df_data["PC10"].dtype == "float64", "PC10 column has wrong data type."

  y = df_data["HDL"] #dependent variable
  X = df_data[["genotype", "sex", "age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"]] #independent variable
  X = sm.add_constant(X) #constant term for intercept

  #initializing Main loop to read variant one by one and compute the association regression
  with pysam.VariantFile(bcf_path, "r") as ifile, gzip.open(output_path, "wt") as ofile:

      #writing header of output file
      ofile.write("CHROM\tPOSITION\tREF_ALLELE\tALT_ALLELE\tN\tEFFECT\tPVALUE\n")

      for bcf_record in ifile:
          np_genotypes = np.full(n_samples, np.nan, dtype = np.float32)

          for i in range(0, n_samples):
              genotype = bcf_record.samples[i]['GT']
              if None in genotype:
                  continue
              np_genotypes[i] = sum(genotype)

          #skip missing values
          if np.isnan(np_genotypes).all():
              continue

        #Update the genotype column in existing X data frame
          X["genotype"] = np_genotypes

        #association testing using OLS
          try:
              regression_model = sm.OLS(y, X, missing='drop')
              results = regression_model.fit()

              effect_size = results.params["genotype"]
              p_value = results.pvalues["genotype"]

              output_line = [
                  str(bcf_record.chrom),
                  str(bcf_record.pos),
                  str(bcf_record.ref),
                  str(bcf_record.alts[0]),
                  str(n_samples),
                  str(effect_size),
                  str(p_value)
                  ]
              #printing the designed format for output file

              ofile.write("\t".join(output_line) + "\n")
          except Exception as e:
              continue

  print (f"GWAS analysis complete. Results saved to {output_path}")

def parse_args():
    parser = argparse.ArgumentParser(description = "GWAS.py is a tool designed to run genome wide association on input bcf files, with necessary input files.")
    parser.add_argument("-g", "--bcf_path", required = True, help="Path to input BCF file in Google drive")
    parser.add_argument("-p", "--pheno_path", required = True, help="Path to input phenotype file in Google drive")
    parser.add_argument("-c", "--covar_path", required = True, help="Path to input covariates file in Google drive")
    parser.add_argument("-o", "--output_path", required = True, help="Output file zipped")
    return parser.parse_args()

def main():
    args = parse_args()
    run_gwas(args.bcf_path, args.pheno_path, args.covar_path, args.output_path)

if __name__ == "__main__":
    main()

    exit(0)
