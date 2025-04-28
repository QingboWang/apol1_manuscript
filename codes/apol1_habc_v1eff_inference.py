import pandas as pd
import numpy as np
import glob
import time as tm
from matplotlib import pyplot as plt
from scipy import stats
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

#Inference of V1 genotype based on the peptide readout:
mat = pd.read_csv("/Users/qsw/Desktop/habc/perpeptide_apol1_expression_mat.tsv.gz", sep="\t", compression="gzip", index_col=0)
mat = mat.loc[:,mat.isna().mean()<0.2] #More than 20% NAs => low quality eptide, remove
mat = mat.astype(float)

#Utilizing the fact K.LEDNIRR.L and K.SELEDNIRR.L intersects with V1, the fraction can be used for "pseudo" genotyping:
annot = x.copy(deep=True)
annot[(y <= x+1) & (y >= x-0.5)] = "REF/ALT" #Manually set threshold
annot[(y > x+1)|(x==0)] = "REF/REF" #Manually set threshold
annot[(y < x-0.5)|(y==0)] = "ALT/ALT" #Manually set threshold
annot[(x==0)&(y==0)] = "REF/ALT"
plt.scatter(x[annot=="ALT/ALT"], y[annot=="ALT/ALT"], alpha=0.5, color="tab:red")
plt.scatter(x[annot=="REF/ALT"], y[annot=="REF/ALT"], alpha=0.5, color="tab:orange")
plt.scatter(x[annot=="REF/REF"], y[annot=="REF/REF"], alpha=0.5, color="tab:grey")
plt.plot(np.arange(-4,3), np.arange(-4,3)+1, color="black", linestyle="--", linewidth=0.5)
plt.plot(np.arange(-4,3), np.arange(-4,3)-0.5, color="black", linestyle="--", linewidth=0.5)
plt.show()

#Now we do the plot for AFR and EUR separately, just for sanity check:
meta = pd.read_csv("~/Desktop/habc/HABC_CALICO_BARCODE_MAP.tsv", sep="\t")
meta.set_index("HABCID", inplace=True)
meta.index = meta.index.astype(str)
pop = meta.RACE[~meta.index.duplicated(keep="first")]
pop = pop[nonna[nonna].index]
x1 = x[nonna][pop==0]
y1 = y[nonna][pop==0]
plt.figure(figsize=(4.5,4.5))
plt.title("Inferred rs2239785 genotype in HealthABC (EUR)", fontsize=12)
plt.scatter(x1[annot=="REF/REF"], y1[annot=="REF/REF"], alpha=0.5, color="tab:grey", label="G/G (REF/REF)")
plt.scatter(x1[annot=="REF/ALT"], y1[annot=="REF/ALT"], alpha=0.5, color="tab:orange", label="G/A (REF/ALT)")
plt.scatter(x1[annot=="ALT/ALT"], y1[annot=="ALT/ALT"], alpha=0.5, color="tab:red", label="A/A (ALT/ALT)")
plt.plot([-6,2.2], [-6,2.2], color="black", linestyle="--", linewidth=0.5, zorder=-2)
plt.axvline(x=0, color="black", linestyle="--", linewidth=0.5, zorder=-2)
plt.axhline(y=0, color="black", linestyle="--", linewidth=0.5, zorder=-2)
plt.xlabel("peptide: K.LEDNIRR.L")
plt.ylabel("peptide: K.SELEDNIRR.L")
plt.xlim([-6,2.2])
plt.ylim([-10,2.2])
plt.legend()
plt.tight_layout()
plt.savefig("habc_v1_inferred_scatter_EUR.png", dpi=400)
plt.savefig("habc_v1_inferred_scatter_EUR.pdf", dpi=400)
plt.close()

pop = meta.RACE[~meta.index.duplicated(keep="first")]
pop = pop[nonna[nonna].index]
x1 = x[nonna][pop==1]
y1 = y[nonna][pop==1]
plt.figure(figsize=(4.5,4.5))
plt.title("Inferred rs2239785 genotype in HealthABC (AFR)", fontsize=12)
plt.scatter(x1[annot=="REF/REF"], y1[annot=="REF/REF"], alpha=0.5, color="tab:grey", label="G/G (REF/REF)")
plt.scatter(x1[annot=="REF/ALT"], y1[annot=="REF/ALT"], alpha=0.5, color="tab:orange", label="G/A (REF/ALT)")
plt.scatter(x1[annot=="ALT/ALT"], y1[annot=="ALT/ALT"], alpha=0.5, color="tab:red", label="A/A (ALT/ALT)")
plt.plot([-6,2.2], [-6,2.2], color="black", linestyle="--", linewidth=0.5, zorder=-2)
plt.axvline(x=0, color="black", linestyle="--", linewidth=0.5, zorder=-2)
plt.axhline(y=0, color="black", linestyle="--", linewidth=0.5, zorder=-2)
plt.xlabel("peptide: K.LEDNIRR.L")
plt.ylabel("peptide: K.SELEDNIRR.L")
plt.xlim([-6,2.2])
plt.ylim([-10,2.2])
plt.legend()
plt.tight_layout()
plt.savefig("habc_v1_inferred_scatter_AFR.png", dpi=400)
plt.savefig("habc_v1_inferred_scatter_AFR.pdf", dpi=400)
plt.close()


#Now look at their association on APOL1 expression (the right peptide):
#Especially focusing on the difference between mass spec vs olink
tb = pd.read_csv("gt_v1_inferred.tsv", sep="\t", index_col=0) #Saved the above result and reading again here
prot = pd.read_csv("proteomics_cape_export.v1.tsv.gz", sep="\t", compression="gzip")
prot.set_index("Gene", inplace=True)
prot = prot.loc[prot.index.intersection(["APOL1", "HP", "HPR","KLKB1","F12", "APOA1"]), :] #The guys we care about
prot["prot_id"] = prot.index + "_" + prot.Protein
mat12 = prot.groupby(["prot_id", "HABCID"])["value"].mean().unstack(level=0) #Ave, accounting for NAs
tb = tb.join(mat12['APOL1_tr|A5PL32|A5PL32_HUMAN'], how="left")
#Olink:
ol0 = pd.read_csv("HABC_Olink_model_dataframe.tsv.gz", sep="\t")
ol0.set_index("HABCID", inplace=True)
tb = tb.join(ol0.groupby("HABCID").APOL1.mean(), how="left")
tb.rename(columns={"APOL1":"APOL1_olink"}, inplace=True)

#Want a few more covariates to plug in for association test:
meta = pd.read_csv("HABC_CALICO_BARCODE_MAP.tsv", sep="\t")
meta.set_index("HABCID", inplace=True)
meta.index = meta.index.astype(int)
tb = tb.join(meta.sort_values(by="timepoint").groupby("HABCID")[["GENDER", "CV1AGE", "BMI"]].head(1), how="left")

#Look at the linear regression results on each of the measurements:
import statsmodels.api as sm
def linear_regression(y: pd.Series, x: pd.Series, covariates: pd.DataFrame):
    model = sm.OLS(y, pd.concat([x, covariates], axis=1), missing="drop").fit()
    results_df = pd.DataFrame({'variable': ['x'] + list(covariates.columns),'beta': model.params,'se': model.bse,'p_value': model.pvalues})
    return results_df
tb["inferred_v1_asdosage"] = tb.inferred_v1.replace({"REF/REF":0, "REF/ALT":1, "ALT/ALT":2})
v1effs = []
#EUR, mass spec
tbe = tb[tb.population==0]
x = tbe.inferred_v1_asdosage
y = tbe["APOL1_tr|A5PL32|A5PL32_HUMAN"]
Z = tbe.loc[:,["GENDER","CV1AGE","BMI", "population"]] #covariate. population = intercept in reality in this case
rmna = (~np.isnan(x)) & (~np.isnan(y)) & (~np.isnan(Z).any(axis=1))
x = x[rmna]
y = y[rmna]
Z = Z[rmna]
v1eff = linear_regression(y, x, Z)
v1eff["population"] = "EUR"
v1eff["measurement"] = "mass"
v1effs.append(v1eff)
#EUR, olink:
x = tbe.inferred_v1_asdosage
y = tbe["APOL1_olink"]
Z = tbe.loc[:,["GENDER","CV1AGE","BMI", "population"]] #covariate. population = intercept in reality in this case
rmna = (~np.isnan(x)) & (~np.isnan(y)) & (~np.isnan(Z).any(axis=1))
x = x[rmna]
y = y[rmna]
Z = Z[rmna]
v1eff = linear_regression(y, x, Z)
v1eff["population"] = "EUR"
v1eff["measurement"] = "olink"
v1effs.append(v1eff)
#AFR:
tba = tb[tb.population==1]
x = tba.inferred_v1_asdosage
y = tba["APOL1_tr|A5PL32|A5PL32_HUMAN"]
Z = tba.loc[:,["GENDER","CV1AGE","BMI", "population"]] #covariate. population = intercept in reality in this case
rmna = (~np.isnan(x)) & (~np.isnan(y)) & (~np.isnan(Z).any(axis=1))
x = x[rmna]
y = y[rmna]
Z = Z[rmna]
v1eff = linear_regression(y, x, Z)
v1eff["population"] = "AFR"
v1eff["measurement"] = "mass"
v1effs.append(v1eff)
#EUR, olink:
x = tba.inferred_v1_asdosage
y = tba["APOL1_olink"]
Z = tba.loc[:,["GENDER","CV1AGE","BMI", "population"]] #covariate. population = intercept in reality in this case
rmna = (~np.isnan(x)) & (~np.isnan(y)) & (~np.isnan(Z).any(axis=1))
x = x[rmna]
y = y[rmna]
Z = Z[rmna]
v1eff = linear_regression(y, x, Z)
v1eff = linear_regression(y, x, Z)
v1eff["population"] = "AFR"
v1eff["measurement"] = "olink"
v1effs.append(v1eff)
v1effs = pd.concat(v1effs)
v1effs = v1effs.loc["inferred_v1_asdosage",:]
#And plot this:
plt.figure(figsize=(4.5,4.5))
for i in range(4):
    plt.text(x=0.02, y=3-i, s="p={0:.3f}".format(v1effs.p_value.values[i]), color="black") #Color based on pval = TODO
plt.errorbar(x=v1effs.beta.values[0],y=3, xerr=v1effs.se.values[0], fmt='o', color="darkgreen")
plt.errorbar(x=v1effs.beta.values[1],y=2, xerr=v1effs.se.values[1], fmt='o', color="tab:green")
plt.errorbar(x=v1effs.beta.values[2],y=1, xerr=v1effs.se.values[2], fmt='o', color="tab:blue")
plt.errorbar(x=v1effs.beta.values[3],y=0, xerr=v1effs.se.values[3], fmt='o', color="tab:cyan")
plt.axvline(x=0, color="black", linestyle="--", linewidth=0.5)
plt.title("Inferred effect of rs2239785 on APOL1 expression in HealthABC", fontsize=12)
plt.xlabel("Effect size of rs2239785 ALT allele (A) (Inferred)")
plt.yticks([0,1,2,3],["AFR, Olink\n(n={0})".format(sum(~tba["APOL1_olink"].isna())),
                      "AFR, mass. spec.\n(n={0})".format(sum(~tba["APOL1_tr|A5PL32|A5PL32_HUMAN"].isna())),
                      "EUR, Olink\n(n={0})".format(sum(~tbe["APOL1_olink"].isna())),
                      "EUR, mass. spec.\n(n={0})".format(sum(~tbe["APOL1_tr|A5PL32|A5PL32_HUMAN"].isna()))])
ax = plt.gca()
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlim([-0.4,0.1])
plt.tight_layout()
plt.show()
