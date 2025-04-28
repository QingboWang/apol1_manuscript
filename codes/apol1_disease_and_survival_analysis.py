#Survival analysis for APOL1_olink

#env: source ~/miniconda3/bin/activate
#conda activate ml_env
from lifelines import CoxPHFitter
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

###Health ABC:
df = pd.read_excel("HABC_DEEPOMICS_Y1_Y2_1111_FASTED_SAMPLES_with_BATCH.xlsx")
df["smoking_binary"] = df.SMK1=="Former"
df["is_male"] = df.GENDER==0
# Now add APOl1 expression:
#From Olink within habc
ol = pd.read_csv("HABC_Olink_model_dataframe.tsv.gz", sep="\t")
ol.set_index("HABCID", inplace=True)
#From mass spec:
prot = pd.read_csv("proteomics_cape_export.v1.tsv.gz", sep="\t", compression="gzip")
ap = prot[prot.Gene.fillna("").str.contains("APOL1")]
ap = ap[ap.Protein=="tr|A5PL32|A5PL32_HUMAN"]
#Take the timepoint mean:
ap = ap.groupby("HABCID")["value"].mean()

#Join:
df.set_index("HABCID", inplace=True)
df = df.join(ap, how="left")
df.rename(columns={"value": "APOL1_mass"}, inplace=True)
ol_apol1 = ol.groupby(ol.index).APOL1.mean()
df = df.join(ol_apol1, how="left")
df.rename(columns={"APOL1": "APOL1_olink"}, inplace=True)

#Add G1/G2/V1 (inferred):
gt = pd.read_csv("v1_g12_guess_from_mass.tsv", sep="\t", index_col=0)
df = df.join(gt, how="left")
#Impute those missing with mean:
df["v1_guess"].fillna(df["v1_guess"].mean(), inplace=True)
df["g12_guess"].fillna(df["g12_guess"].mean(), inplace=True)
#Mass spec:
cols = ["APOL1_mass", "CV1AGE", "RACE", "smoking_binary", "is_male", "BMI", "v1_guess", "g12_guess", "SurvTime", "event"]
dffilt = df[cols]
dffilt = dffilt[~dffilt.APOL1_mass.isna()]
for feat in cols[:-2]:
   dffilt[feat] = (dffilt[feat]-dffilt[feat].mean())/dffilt[feat].std() #Standardize
cph = CoxPHFitter()
cph.fit(dffilt, duration_col='SurvTime', event_col='event')
cph.print_summary()
#Plot
cph.plot()
plt.ylabel("Factors")
plt.title("HealthABC, all N, per SD")
plt.tight_layout()
plt.savefig("habc_survival_APOL1mass_withgenotype_perSD.pdf", dpi=400)
plt.show() 
plt.close()
#Output the summary
cph.summary.to_csv("habc_survival_APOL1mass_withgenotype_perSD.tsv", sep="\t")

###UKB death:
from datetime import datetime
recdate = pd.read_csv("field_53.csv.gz", sep=",", index_col=0)
age = pd.read_csv("field_21022.csv.gz", sep=",", index_col=0)
sex = pd.read_csv("field_31.csv.gz", sep=",", index_col=0)
bmi = pd.read_csv("field_21001.csv.gz", sep=",", index_col=0)
smoking = pd.read_csv("field_20116.csv.gz", sep=",", index_col=0)
death = pd.read_csv("field_40000.csv.gz", sep=",", index_col=0)
esrd = pd.read_csv("field_42026.csv.gz", sep=",", index_col=0)
#Engineer:
recdate = recdate[recdate.instance_id==0]["value"] #Only at the first arrival
recdate.name = "recdate"
death_date = death["value"] #Only at the first arrival (confirmed by death.instance_id.value_counts())
death_date.name = "death_date"
esrd_date = esrd["value"] #Only at the first arrival (confirmed by esrd.instance_id.value_counts())
esrd_date.name = "esrd_date"
age = age["value"]
age.name = "age"
sex = sex["value"] #Sex has only one instance: (sex.instance_id==0).value_counts()
sex.name = "is_male"
bmi = bmi[bmi.instance_id==0]["value"] #Only at the first arrival
bmi.name = "bmi"
smoking = smoking[smoking.instance_id==0]["value"] #Only at the first arrival
smoking.name = "smoking_ever"
mapper = {0: 0, 1: 1, 2: 1, -3:np.nan} #Encoding: never-ever style https://biobank.ctsu.ox.ac.uk/ukb/coding.cgi?id=90
smoking= smoking.map(mapper)

#Day until death:
end_date = datetime.strptime("2024-09-01", "%Y-%m-%d").date()
def date_difference_in_years(date1, date2):
  difference_in_days = (date2 - date1).days #Signed
  # Approximate the difference in years
  difference_in_years = round(difference_in_days / 365.2425, 1) 
  return difference_in_years
dates = pd.DataFrame(recdate).join(death_date, how="left").join(esrd_date, how="left")
dates["current_date"] = end_date
dates["recdate"] = dates.recdate.apply(lambda x: datetime.strptime(x, "%Y-%m-%d").date() if pd.notnull(x) else x) #Formatting
dates["esrd_date"] = dates.esrd_date.apply(lambda x: datetime.strptime(x, "%Y-%m-%d").date() if pd.notnull(x) else x) #Formatting
dates["death_date"] = dates.death_date.apply(lambda x: datetime.strptime(x, "%Y-%m-%d").date() if pd.notnull(x) else x) #Formatting
dates["years_until_esrd"] = dates.apply(lambda x: date_difference_in_years(x["recdate"], x["esrd_date"]) if pd.notnull(x["esrd_date"]) else np.nan, axis=1) #Negative => Exclude later
dates["years_until_death"] = dates.apply(lambda x: date_difference_in_years(x["recdate"], x["death_date"]) if pd.notnull(x["death_date"]) else np.nan, axis=1)
dates["years_until_now"] = dates.apply(lambda x: date_difference_in_years(x["recdate"], x["current_date"]), axis=1)
#If died, that is the end of cencus date:
dates["years_until_final_censor"] = dates.apply(lambda x: x["years_until_death"] if pd.notnull(x["years_until_death"]) else x["years_until_now"], axis=1)
#annotate event and annotate years until the event:
dates["event_indicator"] = dates["years_until_death"].apply(lambda x: 1 if pd.notnull(x) else 0)
dates["time_to_event"] = dates["years_until_death"].fillna(dates["years_until_final_censor"])

#Concatenate with other features:
df = dates[["event_indicator", "time_to_event"]].join(age, how="left").join(sex, how="left").join(bmi, how="left").join(smoking, how="left")

#PPP APOL1 data (already filtered previously):
p_eur = pd.read_csv("ppp_eur_expression_mat_naincluded.tsv.gz", sep="\t", index_col=0, usecols=["SampleID", "OID30708"])
p_afr = pd.read_csv("ppp_afr_expression_mat_naincluded.tsv.gz", sep="\t", index_col=0, usecols=["SampleID", "OID30708"])
apol1 = pd.concat([p_eur, p_afr])
#Population data (already filtered previously):
sample_eur = pd.read_csv("samples_to_use_eur_genotypefilt.txt", sep="\t", index_col=0)
sample_eur["is_afr"] = 0
sample_afr = pd.read_csv("samples_to_use_afr_genotypefilt.txt", sep="\t", index_col=0)
sample_afr["is_afr"] = 1
populations = pd.concat([sample_eur, sample_afr])
#Join these two:
df2 = pd.DataFrame(populations).join(apol1, how="inner")
df2.rename(columns={"OID30708":"apol1"}, inplace=True)
df2 = df2[~df2.apol1.isna()]
df = df.join(df2, how="inner")
df.event_indicator.value_counts() #Very few cases...
df.fillna(df.mean(axis=0), inplace=True) #Filling missing BMI etc with mean

#Now add exome info:
exome = pd.read_csv("ukb_exome_g1g2v1_dosages.tsv.gz", sep="\t", index_col=0)
df = df.join(exome, how="inner") #Require exome info.
del df['G1G'] #G1M is enough
df["G1M"].fillna(df["G1M"].mode()[0], inplace=True)
df["G2"].fillna(df["G2"].mode()[0], inplace=True)
df["V1"].fillna(df["V1"].mode()[0], inplace=True)

#In local, read, flip V1 and do it:
#Flip V1 since we define it as REF:
df["V1"] = 2 - df.V1
#Do the coxph:
from lifelines import CoxPHFitter
cph = CoxPHFitter()
cph.fit(df, duration_col='time_to_event', event_col='event_indicator')
cph.print_summary()
from matplotlib import pyplot as plt
cph = CoxPHFitter()
cph.fit(df, duration_col='time_to_event', event_col='event_indicator')
cph.plot()
plt.title("Death Harzard Ratio in UKB EUR/AFR")
plt.ylabel("Factor")
plt.tight_layout()
plt.savefig("/Users/qsw/Downloads/ukb_apol1_survival_withgenotype_v1ref.pdf", dpi=400)
plt.show()
plt.close()

###Linear regression and logistic regression in UKB CKD and eGFR:
import pandas as pd
import numpy as np
from datetime import datetime

#In local:
df = pd.read_csv("apol1_logr_input.tsv.gz", sep="\t", compression="gzip", index_col=0)
#Add exome:
exome = pd.read_csv("~/Desktop/ukb_data/ukb_exome_g1g2v1_dosages.tsv.gz", sep="\t", index_col=0)
df = df.join(exome, how="inner")
del df["G1G"]
#Fill with mode (as there are very few missing)
df["G1M"].fillna(df["G1M"].mode()[0], inplace=True)
df["G2"].fillna(df["G2"].mode()[0], inplace=True)
df["V1"].fillna(df["V1"].mode()[0], inplace=True)
#Flip V1:
df["V1"] = 2 - df["V1"]

#Logistic regression weight:
import statsmodels.api as sm
def logistic_regression_odds_ratio(df, feature_cols, outcome_col):
  #Remove those where outcome is NA:
  df = df[~df[outcome_col].isna()]
  # Add a constant term to the features for the intercept
  X = sm.add_constant(df[feature_cols])
  y = df[outcome_col]
  # Fit the logistic regression model
  logit_model = sm.Logit(y, X)
  result = logit_model.fit()
  # Extract odds ratios and p-values
  odds_ratios = pd.DataFrame({'OddsRatio': result.params,'Pval': result.pvalues, 'Lower CI': result.conf_int()[0], 'Upper CI': result.conf_int()[1]})
  # Remove the intercept row and return
  return odds_ratios.iloc[1:]
features = ["is_afr", "age", "is_male", "bmi", "smoking_ever", "G1M", "G2", "V1", "apol1"]
outcome = "esrd"
pvals = logistic_regression_odds_ratio(df, features, outcome)
print(pvals)
features = ["is_afr", "age", "is_male", "bmi", "smoking_ever", "G1M", "G2", "V1", "apol1"]
outcome = "ckd"
pvals = logistic_regression_odds_ratio(df, features, outcome)
print(pvals)

#Also egfr:
def linear_regression_p_values(df, feature_cols, outcome_col):
  #Remove those where outcome is NA:
  df = df[~df[outcome_col].isna()]  
  # Add a constant term to the features for the intercept
  X = sm.add_constant(df[feature_cols])
  y = df[outcome_col]
  # Fit the linear regression model
  model = sm.OLS(y, X)
  result = model.fit()
  # Extract coefficients and p-values
  results_df = pd.DataFrame({'beta': result.params,'P-value': result.pvalues, 'Lower CI': result.conf_int()[0], 'Upper CI': result.conf_int()[1]})
  return results_df.iloc[1:]
df = pd.read_csv("apol1_logr_input.tsv.gz", sep="\t", compression="gzip", index_col=0)
#Add exome:
exome = pd.read_csv("ukb_exome_g1g2v1_dosages.tsv.gz", sep="\t", index_col=0)
df = df.join(exome, how="inner")
del df["G1G"]
df["G1M"].fillna(df["G1M"].mean(), inplace=True)
df["G2"].fillna(df["G2"].mean(), inplace=True)
df["V1"].fillna(df["V1"].mean(), inplace=True)
#Flip V1:
df["V1"] = 2 - df["V1"]
#Join egfr
egfr = pd.read_csv("egfr_3formulas.txt", sep="\t", index_col=0)
df = df.join(egfr, how="left")
#Standardize features:
features = ["is_afr", "age", "is_male", "bmi", "smoking_ever", "G1M", "G2", "V1", "apol1"]
for feat in features:
   df[feat] = (df[feat] - df[feat].mean()) / df[feat].std()
outcome = "esrd"
pvals = logistic_regression_odds_ratio(df, features, outcome)
print(pvals)
features = ["is_afr", "age", "is_male", "bmi", "smoking_ever", "G1M", "G2", "V1", "apol1"]
outcome = "ckd"
pvals = logistic_regression_odds_ratio(df, features, outcome)
print(pvals)
pvals.to_csv("~/Downloads/ukb_apol1_ckd_logr_persd.tsv", sep="\t")
features = ["is_afr", "age", "is_male", "bmi", "smoking_ever", "G1M", "G2", "V1", "apol1"]
for outcome in ["egfr_cr2021", "egfr_crcys2021", "egfr_cys2012"]:
    pvals = linear_regression_p_values(df, features, outcome)
    print (outcome)
    print(pvals)
    pvals.to_csv("~/Downloads/ukb_apol1_"+outcome+"_linr_persd.tsv", sep="\t")
#Plot
pvals = pd.read_csv("~/Downloads/ukb_apol1_ckd_logr_persd.tsv", sep="\t", index_col=0)
pvals.sort_values(by="OddsRatio", ascending=True, inplace=True)
plt.figure(figsize=(8,3.5))
plt.title("Odds ratio for CKD in UKB EUR/AFR")
plt.errorbar(x = pvals.OddsRatio, xerr = [pvals.OddsRatio-pvals["Lower CI"], pvals["Upper CI"]-pvals.OddsRatio], 
             y = np.arange(pvals.shape[0]), fmt="o", color="black")
plt.xlabel("Per SD Log(Odds Ratio)")
plt.ylabel("Factor")
plt.yticks(np.arange(pvals.shape[0]), pvals.index)
plt.axvline(x=0, linestyle="--", color="tab:gray", zorder=-2)
plt.tight_layout()
plt.savefig("ukb_apol1_ckd_logr_withgenotype_v1ref_persd.pdf", dpi=400)
plt.show()
plt.close()
#Plot EGFR assoc:
for outcome in ["egfr_cr2021", "egfr_crcys2021", "egfr_cys2012"]:
    pvals = pd.read_csv("~/Downloads/ukb_apol1_"+outcome+"_linr_persd.tsv", sep="\t", index_col=0)
    pvals.sort_values(by="beta", ascending=True, inplace=True)
    plt.figure(figsize=(8,3.5))
    plt.title("Per SD Beta for {0} in UKB EUR/AFR".format(outcome))
    plt.errorbar(x = pvals.beta, xerr = [pvals.beta-pvals["Lower CI"], pvals["Upper CI"]-pvals.beta], 
                y = np.arange(pvals.shape[0]), fmt="o", color="black")
    plt.xlabel("Beta")
    plt.ylabel("Factor")
    plt.yticks(np.arange(pvals.shape[0]), pvals.index)
    plt.axvline(x=0, linestyle="--", color="tab:gray", zorder=-2)
    plt.tight_layout()
    plt.savefig("ukb_apol1_{0}_linr_withgenotype_v1ref_persd.pdf".format(outcome), dpi=400)
    plt.show()
    plt.close()

#per SD for the coxPH as well:
df = pd.read_csv("ukb_apol1_df_for_survival_withgenotype.tsv.gz", sep="\t", compression="gzip", index_col=0)
#Flip V1:
df["V1"] = 2 - df.V1
#SD:
for feat in df.columns.difference(["event_indicator", "time_to_event"]):
   df[feat] = (df[feat] - df[feat].mean())/df[feat].std()
#Do the coxph:
from lifelines import CoxPHFitter
cph = CoxPHFitter()
cph.fit(df, duration_col='time_to_event', event_col='event_indicator')
cph.print_summary()
from matplotlib import pyplot as plt
cph = CoxPHFitter()
cph.fit(df, duration_col='time_to_event', event_col='event_indicator')
cph.plot()
plt.title("Death Harzard Ratio in UKB EUR/AFR (per SD)")
plt.ylabel("Factor")
plt.tight_layout()
plt.savefig("ukb_apol1_survival_withgenotype_v1ref_persd.pdf", dpi=400)
plt.show()
plt.close()
cph.summary.to_csv("ukb_apol1_survival_withgenotype_v1ref_persd.tsv", sep="\t")




