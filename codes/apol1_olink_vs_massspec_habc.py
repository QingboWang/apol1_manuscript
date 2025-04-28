import pandas as pd
import numpy as np
import glob
import time as tm
from matplotlib import pyplot as plt
from scipy import stats
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

tb = pd.read_csv("habc_peptide_inferred_genotypes.tsv.gz", sep="\t", compression="gzip", index_col=0) #We actually don't even need the inferred genotype here, but just since this df has everything.
tbe = tb[tb.population==0]
tbe = tbe[~(tbe['APOL1_tr|A5PL32|A5PL32_HUMAN'].isna()|tbe['APOL1_olink'].isna())]
plt.scatter(tbe['APOL1_tr|A5PL32|A5PL32_HUMAN'], tbe['APOL1_olink'], alpha=0.8)
pearsonr = stats.pearsonr(tbe['APOL1_tr|A5PL32|A5PL32_HUMAN'], tbe['APOL1_olink'])
plt.xlabel("APOL1 mass spec.")
plt.ylabel("APOL1 Olink")
plt.title("APOL1 expression correlation\nin HealthABC EUR (r={0:.3f})".format(pearsonr[0]))
plt.show()
tba = tb[tb.population==1]
tba = tba[~(tba['APOL1_tr|A5PL32|A5PL32_HUMAN'].isna()|tba['APOL1_olink'].isna())]
plt.scatter(tba['APOL1_tr|A5PL32|A5PL32_HUMAN'], tba['APOL1_olink'], alpha=0.8)
pearsonr = stats.pearsonr(tba['APOL1_tr|A5PL32|A5PL32_HUMAN'], tba['APOL1_olink'])
plt.xlabel("APOL1 mass spec.")
plt.ylabel("APOL1 Olink")
plt.title("APOL1 expression correlation\nin HealthABC AFR (r={0:.3f})".format(pearsonr[0]))
plt.show()