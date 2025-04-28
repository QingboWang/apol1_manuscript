import pandas as pd
import numpy as np
import glob
import time as tm
from matplotlib import pyplot as plt
from scipy import stats
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
import statsmodels.api as sm
def linear_regression(y: pd.Series, x: pd.Series, covariates: pd.DataFrame):
    model = sm.OLS(y, pd.concat([x, covariates], axis=1), missing="drop").fit()
    results_df = pd.DataFrame({'variable': ['x'] + list(covariates.columns),'beta': model.params,'se': model.bse,'p_value': model.pvalues})
    return results_df
from scipy.stats import norm, rankdata
def inverse_normal_transform(series):
  ranks = rankdata(series.dropna(), method='average')
  normalized_ranks = (ranks - 0.5) / len(ranks)
  transformed_values = norm.ppf(normalized_ranks)
  result = pd.Series(index=series.index)
  result.loc[series.dropna().index] = transformed_values
  return result

### Imputation of G1/G2 from peptide level data:

#We essentially build a model to predict the peptide in EUR using all the other peptides
#And then use the deviation from the prediction as the score (either inv-normed or not. Start with not)
tb = pd.read_csv("habc_peptide_inferred_genotypes.tsv.gz", sep="\t", compression="gzip", index_col=0)
mat = pd.read_csv("perpeptide_apol1_expression_mat.tsv.gz", sep="\t", compression="gzip", index_col=0)
mat = mat.loc[:,mat.T.protein_id=="tr|A5PL32|A5PL32_HUMAN"]
mat = mat.iloc[2:,:]
mat = mat.loc[:,mat.isna().mean()<0.3]
mat = mat.astype(float)
mat.index = mat.index.astype(int)
#And join at this stage:
tb = tb.join(mat, how="left")
tba = tb[tb.population==1]
tbe = tb[tb.population==0]
collist = tba.columns[-16:]
collist_feat = collist.difference(["K.LNILNNNYK.I"])

#Train a linear model in EUR:
from sklearn.linear_model import Lasso
model = Lasso(alpha=0.01)
X = tbe[collist_feat].copy(deep=True)
X.fillna(X.mean(axis=0), inplace=True)
y = tbe["K.LNILNNNYK.I"].copy(deep=True)
nans = ( y.isna() ) #NA and outlier values - not use for training.
X = X[~nans]
y = y[~nans]
model.fit(X, y)
y_pred = model.predict(X) #Predicted y within EUR
# Make predictions on the new data = AFR => What is deviating is the existence of 
Xafr = tba[collist_feat].copy(deep=True)
Xafr.fillna(X.mean(axis=0), inplace=True)
afr_pred = model.predict(Xafr)
#Fill NAs back with NA for AFR:
afr_pred = pd.Series(afr_pred, index=Xafr.index)
afr_pred[tba["K.LNILNNNYK.I"].isna()] = np.nan
#Plot obs vs exp for sanity check:
plt.scatter(tba["K.LNILNNNYK.I"], afr_pred, color="tab:blue")
plt.axvline(0, color="black", linestyle="--", linewidth=0.5)
plt.axhline(0, color="black", linestyle="--", linewidth=0.5)
plt.plot([-1,1],[-1,1], color="black", linestyle="--", linewidth=0.5)
plt.xlabel("Observed K.LNILNNNYK.I")
plt.ylabel("Predicted K.LNILNNNYK.I")
plt.show()

#Have a (manually defined) decision boundary in a linear model:
#And plot that information too in the plot:
model = Lasso(alpha=0.01)
X = tbe[collist_feat].copy(deep=True)
X.fillna(X.mean(axis=0), inplace=True)
y = tbe["K.LNILNNNYK.I"].copy(deep=True)
nans = ( y.isna() ) #NA and outlier values - not use for training.
X = X[~nans]
y = y[~nans]
model.fit(X, y)
y_pred = model.predict(X) #For later plot
tba["g12_guess"] = 0
tba.loc[tba["K.LNILNNNYK.I"]<tba.predicted_peptide_read*5/3, "g12_guess"] = 1 #Manually defined decision boundary
tba.loc[tba["K.LNILNNNYK.I"]<tba.predicted_peptide_read*5/3-1, "g12_guess"] = 2 #Manually defined decision boundary
tba.loc[tba["K.LNILNNNYK.I"]==0, "g12_guess"] = 2 #Manually defined decision boundary
plt.figure(figsize=(5, 5))
plt.scatter(y_pred, y, color="tab:blue", alpha=0.2, label="EUR, 0")
plt.scatter(tba[tba.g12_guess==0].predicted_peptide_read, tba[tba.g12_guess==0]["K.LNILNNNYK.I"], color="tab:grey", alpha=0.5, label="AFR, 0")
plt.scatter(tba[tba.g12_guess==1].predicted_peptide_read, tba[tba.g12_guess==1]["K.LNILNNNYK.I"], color="tab:orange", alpha=0.5, label="AFR, 1")
plt.scatter(tba[tba.g12_guess==2].predicted_peptide_read, tba[tba.g12_guess==2]["K.LNILNNNYK.I"], color="tab:red", alpha=0.5, label="AFR, 2")
plt.axhline(0, color="black", linestyle="--", linewidth=0.5)
plt.axvline(0, color="black", linestyle="--", linewidth=0.5)
plt.plot([-10,10], [-10*5/3,10*5/3], color="black", linestyle="--", linewidth=0.5)
plt.plot([-10,10], [-10*5/3 -1,10*5/3 -1], color="black", linestyle="--", linewidth=0.5)
plt.ylabel("Observed K.LNILNNNYK.I")
plt.xlabel("Predicted K.LNILNNNYK.I\n(based on other peptides)")
plt.legend(title = "Population, G1/G2\ndosage estimation", fontsize=11)
plt.title("Inferred G1/G2 dosage in HealthABC")
plt.xlim([-1.1,0.3])
plt.ylim([-4,1])
plt.tight_layout()
plt.savefig("apol_habc_g12_estimation.png", dpi=400)
plt.savefig("apol_habc_g12_estimation.pdf", dpi=400)
plt.close()

#And the effect size:
g12effs = []
tba.loc[tba.g12_mass.isna(),"g12_guess"] = np.nan
tba.loc[tba.g12_mass==0,"g12_guess"] = np.nan #Also.
x, y, Z = tba.g12_guess, tba["APOL1_tr|A5PL32|A5PL32_HUMAN"], tba.loc[:,["GENDER","CV1AGE","BMI", "inferred_v1_asdosage", "population"]] #covariate. population = intercept in reality in this case
rmna = (~np.isnan(x)) & (~np.isnan(y)) & (~np.isnan(Z).any(axis=1))
x, y, Z = x[rmna], y[rmna], Z[rmna]
g12eff = linear_regression(y, x, Z)
g12eff["population"] = "AFR"
g12eff["measurement"] = "mass"
g12effs.append(g12eff)
#AFR, olink:
x, y, Z = tba.g12_guess, tba["APOL1_olink"], tba.loc[:,["GENDER","CV1AGE","BMI", "inferred_v1_asdosage", "population"]] #covariate. population = intercept in reality in this case
rmna = (~np.isnan(x)) & (~np.isnan(y)) & (~np.isnan(Z).any(axis=1))
x, y, Z = x[rmna], y[rmna], Z[rmna]
g12eff = linear_regression(y, x, Z)
g12eff["population"] = "AFR"
g12eff["measurement"] = "olink"
g12effs.append(g12eff)
g12effs = pd.concat(g12effs)
g12effs = g12effs.loc["g12_guess",:]

#eff size difference p-value:
effect_size_1 = g12effs.loc[g12effs.measurement=="mass", 'beta']
SE_1 = g12effs.loc[g12effs.measurement=="mass", 'se']
effect_size_2 = g12effs.loc[g12effs.measurement=="olink", 'beta']
SE_2 = g12effs.loc[g12effs.measurement=="olink", 'se']
z_statistic = (effect_size_1 - effect_size_2) / (SE_1**2 + SE_2**2)**0.5
p_value = 2 * (1 - stats.norm.cdf(abs(z_statistic)))  # Two-tailed test

plt.figure(figsize=(5, 2.5))
plt.errorbar(x=g12effs.beta.values[0],y=1, xerr=g12effs.se.values[0], fmt='o', color="black")
plt.errorbar(x=g12effs.beta.values[1],y=0, xerr=g12effs.se.values[1], fmt='o', color="black")
plt.axvline(x=0, color="black", linestyle="--", linewidth=0.5)
plt.title("Inferred effect of G1/G2\non APOL1 expression in HealthABC", fontsize=12)
plt.xlabel("Effect size of G1/G2 (Per dosage, inferred)")
plt.yticks([0,1],["Olink\n(AFR n={0})".format(sum(~(tba["APOL1_olink"].isna()|tba.g12_guess.isna()))),
                      "mass. spec.\n(AFR n={0})".format(sum(~tba["APOL1_tr|A5PL32|A5PL32_HUMAN"].isna()|tba.g12_guess.isna()))])
plt.text(.92, 0.4, "p={0:.4f}".format(p_value[0]), fontsize=12)
ax = plt.gca()
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlim([-0.04,1.1])
plt.ylim([-0.2,1.2])
plt.tight_layout()
plt.savefig("apol_habc_g1g2_effsize_comparison.png", dpi=400)
plt.savefig("apol_habc_g1g2_effsize_comparison.pdf", dpi=400)
plt.close()