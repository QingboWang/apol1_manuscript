import pandas as pd
import numpy as np
import glob
import time as tm
from matplotlib import pyplot as plt
from scipy import stats
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"


#G1, G2 and V1 causal eff size comparison (posterior beta from FINEMAP):
#AFR:
fm = pd.read_csv("~/Desktop/ukb_data/0919dl_mainfig_source_apol1_afr/apol1_ukb_afr_geno_ex_hybrid_n947_cis.snp", sep=" ")
sus = pd.read_csv("~/Desktop/ukb_data/0919dl_mainfig_source_apol1_afr/ppp_apol1_afr_geno_ex_hybrid_susie_pip.tsv", sep="\t")
fm.set_index("rsid", inplace=True)
sus.set_index("rsid", inplace=True)
fm = fm.join(sus, how="left")
fm["min_pip"] = np.minimum(fm.prob, fm.pip)
fm.sort_values(by="min_pip", ascending=False).head(1).T #G2
fm.sort_values(by="min_pip", ascending=False).head(3).tail(1).T #G1
fm[fm.index=="22:36265284:G:A"].T #V1
#EUR (V1):
fme = pd.read_csv("~/Desktop/ukb_data/eur_ss_for_apol1_mainfig/156_APOL1_OID30708/156_APOL1_OID30708.snp", sep=" ")
fme.set_index("rsid", inplace=True)

#Combine:
eff = pd.concat([fm.loc[["22:36265995:AATAATT:A", "22:36265988:T:G", "22:36265284:G:A"],["beta", "se", "mean_incl", "sd_incl"]], fme[fme.position==36265284][["beta", "se", "mean_incl", "sd_incl"]]])
eff.index = ["G2", "G1", "V1", "V1_eur"]
eff = eff.loc[["G2", "G1", "V1"],:]
colors = ["tab:purple", "tab:red", "tab:green"]
plt.figure(figsize=(10,2.2))
for i in range(3):
    plt.errorbar([i], abs(eff.mean_incl.values[i]), yerr=eff.sd_incl.values[i], fmt="o", color=colors[i])
    plt.axhline(y=abs(eff.mean_incl.values[i]), color=colors[i], linestyle="--", linewidth=0.5)
ax = plt.gca()  # Get the current axes
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)  
plt.ylim([0,1.5])
plt.xlim([-0.5,2.5])
plt.xticks([0,1,2], ["G2", "G1M", "rs2239785"])
plt.ylabel("|Causal effect size|\n(Olink NPX unit)")
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/v1_beta_vsg1g2.png", dpi=500)
plt.savefig("/Users/qsw/Desktop/misc_plots/v1_beta_vsg1g2.pdf", dpi=500)
plt.show()
plt.close()