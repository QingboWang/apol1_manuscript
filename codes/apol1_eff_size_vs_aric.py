#Comparing the effect sizes of APOL1 cis-variants in UKB and ARIC

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

#AFR:
ss = pd.read_parquet("156_APOL1_OID30708.cis_qtl_pairs.chr22.parquet")
ss["ref_is_minor"] = ss.af>0.5
aa = pd.read_csv("apol_finemapping_AA.csv", sep=",")
aa.index = "chr" + aa.CHR.astype(str) + "_" + aa["position (hg38)"].astype(str) + "_" + aa["Tested Allele"]
aa.columns = aa.columns + "_aric"
ss.index = ss.variant_id.str.split("_").str[0]+"_"+ss.variant_id.str.split("_").str[1] + "_" + ss.variant_id.str.split("_").str[3].str.split(";").str[0]
ss1 = ss.join(aa, how="left")
ss.index = ss.variant_id.str.split("_").str[0]+"_"+ss.variant_id.str.split("_").str[1] + "_" + ss.variant_id.str.split("_").str[2]
ss2 = ss.join(aa, how="left")
ss2["Beta_aric"] = -ss2.Beta_aric #For those where ref/alt are flipped, we flip the beta
ss1 = ss1[~ss1.Beta_aric.isna()]
ss2 = ss2[~ss2.Beta_aric.isna()]
df = pd.concat([ss1,ss2], axis=0)
g1 = df[df.variant_id.str.startswith("chr22_36265988_T_G")]
g2_proxy = df[df.variant_id.str.startswith("chr22_36269923_C_T")] #Proxy of G2, in good level of LD
v1 = df[df.variant_id.str.startswith("chr22_36265284_G_A")]#rs2239785
plt.errorbar(x=df.slope, y=df.Beta_aric,yerr=df.SE_aric,xerr=df.slope_se, fmt='o', color="tab:blue", label="Others")
plt.errorbar(x=g1.slope, y=g1.Beta_aric,yerr=g1.SE_aric,xerr=g1.slope_se, fmt='o', color="tab:orange", label="G1 (rs60910145)")
plt.errorbar(x=g2_proxy.slope, y=g2_proxy.Beta_aric,yerr=g2_proxy.SE_aric,xerr=g2_proxy.slope_se, fmt='o', color="tab:green", label="Proxy of G2\n(rs6000222)")
plt.errorbar(x=v1.slope, y=v1.Beta_aric,yerr=v1.SE_aric,xerr=v1.slope_se, fmt='o', color="tab:red", label="V1 (rs2239785)")
plt.axvline(x=0, linestyle="--", color="black")
plt.axhline(y=0, linestyle="--", color="black")
plt.legend(title="Variant:")
plt.title("Effect sizes of cis-variants on APOL1")
plt.xlabel("Beta in UKB AFR")
plt.ylabel("Beta in ARIC AA")
plt.show()

#EUR:
ss = pd.read_parquet("156_APOL1_OID30708.cis_qtl_pairs.chr22.parquet")
ss["ref_is_minor"] = ss.af>0.5
aa = pd.read_csv("apol_finemapping_EA.csv", sep=",")
aa.index = "chr" + aa.CHR.astype(str) + "_" + aa["position (hg38)"].astype(str) + "_" + aa["Tested Allele"]
aa.columns = aa.columns + "_aric"
ss.index = ss.variant_id.str.split("_").str[0]+"_"+ss.variant_id.str.split("_").str[1] + "_" + ss.variant_id.str.split("_").str[3].str.split(";").str[0]
ss1 = ss.join(aa, how="left")
ss.index = ss.variant_id.str.split("_").str[0]+"_"+ss.variant_id.str.split("_").str[1] + "_" + ss.variant_id.str.split("_").str[2]
ss2 = ss.join(aa, how="left")
ss2["Beta_aric"] = -ss2.Beta_aric  #For those where ref/alt are flipped, we flip the beta
ss1 = ss1[~ss1.Beta_aric.isna()]
ss2 = ss2[~ss2.Beta_aric.isna()]
df = pd.concat([ss1,ss2], axis=0)
v1 = df[df.variant_id.str.contains("rs2239785")]
v2 = df[df.variant_id.str.contains("rs6000220")]
plt.errorbar(x=df.slope, y=df.Beta_aric,yerr=df.SE_aric,xerr=df.slope_se, fmt='o', color="tab:blue", label="Others")
plt.errorbar(x=v1.slope, y=v1.Beta_aric,yerr=v1.SE_aric,xerr=v1.slope_se, fmt='o', color="tab:green", label='"v1" (rs2239785)')
plt.errorbar(x=v2.slope, y=v2.Beta_aric,yerr=v2.SE_aric,xerr=v2.slope_se, fmt='o', color="tab:orange", label='"v2" (rs6000220)')
plt.axvline(x=0, linestyle="--", color="black")
plt.axhline(y=0, linestyle="--", color="black")
plt.legend(title="Variant:")
plt.title("Effect sizes of cis-variants on APOL1")
plt.xlabel("Beta in UKB EUR")
plt.ylabel("Beta in ARIC EA")
plt.show()
