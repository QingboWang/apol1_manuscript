#Checking the cis-pQTL effect of KLKB1, F12 etc on each other and on APOL1, HPR, KNG1 etc:
import pandas as pd
import numpy as np
import glob
import time as tm
from matplotlib import pyplot as plt
from scipy import stats
from sklearn.linear_model import LinearRegression
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

#Extracted relevant variants from the full exome data:
dos = pd.read_csv("ppp_apol1_eur_causalvars_and_apol1.tsv.gz", sep="\t", compression="gzip", index_col=0)   
dos.rename(columns = {'5:177409531:A:G':"f12", '4:186236880:G:A':"klkb1", '3:186743735:C:T':"kng_lof",'22:36253920:C:T':"v2"}, inplace=True)
dos["kng_lof_any"] = (dos.kng_lof > 0).astype(int)
#HPR:
hpr = pd.read_csv("hpr_representative_variant_genotype.tsv.gz", sep="\t", compression="gzip", index_col=0).squeeze()
hpr.name = "hpr" #This is rs9930957
hpr = hpr.round(0)
dos = dos.join(hpr, how="left")

#Expression data (in EUR) for F12 and KLKB1 as well:
exp = pd.read_csv("ppp_eur_expression_mat_naincluded.tsv.gz", sep="\t", usecols=["SampleID", "OID30760", "OID30739"], compression="gzip", index_col=0)
exp.rename(columns={"OID30760":"KLKB1", "OID30739":"F12"}, inplace=True)#KLKB1, F12
#Violin plot:
dos = dos.join(exp, how="left")
fig, ax = plt.subplots(nrows=5, ncols=5, figsize=(12,12), sharex=True, sharey=True)
for i, prot in zip([1,2,4], ["KLKB1", "F12", "apol1"]):
    for j, var in enumerate(['kng_lof', 'klkb1', 'f12', "hpr", 'v1']):
        vio = []
        for cnt in range(3):
            y = dos[dos[var]==cnt][prot]
            y.dropna(inplace=True)
            vio.append(y.tolist())
        parts = ax[i][j].violinplot(vio, showmeans=False, showextrema=False, widths=0.7, points=1000)
        mu = dos.groupby(var)[prot].mean()
        sig = dos.groupby(var)[prot].sem()
        ax[i][j].errorbar([1,2,3], mu, yerr=sig, fmt="o", color="black")
        #Linear regression:
        model = LinearRegression()
        d = dos[[var, prot]]
        d.dropna(inplace=True)
        model.fit(pd.DataFrame(d[var]), d[prot])
        slope = model.coef_[0]
        intercept = model.intercept_
        pval = stats.linregress(d[var], d[prot])[3]
        ax[i][j].plot([1,3],np.array([0,2])*slope+intercept, linestyle="--", linewidth=1, color="black")
        ax[i][j].set_title("slope={:.3f}, p={:.2e}".format(slope, max(pval, 1e-320)), fontsize=12)
        #Change color if significant:
        col = "tab:grey"
        if pval < 5e-8:
            if slope>0: col = "tab:red"
            else: col = "tab:blue"
        for pc in parts['bodies']:
            pc.set_facecolor(col)  # Set the desired color
            pc.set_edgecolor('black')
            pc.set_alpha(0.5)
for i in [1,2,4]:
    ax[i][0].set_ylabel("{0} expression (NPX)".format(["KNG1", "KLKB1",  "F12", "HPR", "APOL1"][i]))
for i in [0,3]:
    ax[i][0].set_ylabel("{0} expression (NPX):\nNot included in Olink".format(["KNG1", "KLKB1",  "F12", "HPR", "APOL1"][i]))
for j in range(5):
    ax[4][j].set_xticks([1,2,3], ["REF/REF", "REF/ALT", "ALT/ALT"], fontsize=11)
    ax[4][j].set_xlabel("{0} cis-variant".format(["KNG1", "KLKB1",  "F12", "HPR", "APOL1"][j]))
ax[0][0].set_ylim([-1.5, 1.5])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/klkb1_etc_violins_lofver.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/klkb1_etc_violins_lofver.pdf", dpi=400)
plt.close()


#Just go with mean and sem:
fig, ax = plt.subplots(nrows=5, ncols=5, figsize=(14,12), sharex=True, sharey=False)
for i, prot in zip([1,2,4], ["KLKB1", "F12", "apol1"]):
    for j, var in enumerate(['kng1_causal', 'klkb1', 'f12', "hpr", 'v1']):
        #Linear regression:
        model = LinearRegression()
        d = dos[[var, prot]]
        d.dropna(inplace=True)
        model.fit(pd.DataFrame(d[var]), d[prot])
        slope = model.coef_[0]
        intercept = model.intercept_
        pval = stats.linregress(d[var], d[prot])[3]
        #Change color if significant:
        col = "tab:grey"
        if pval < 5e-8:
            if slope>0: col = "tab:red"
            else: col = "tab:blue"
        #Scatter plot:
        mu = dos.groupby(var)[prot].mean()
        sig = dos.groupby(var)[prot].sem()
        ax[i][j].errorbar([1,2,3], mu, yerr=sig, fmt="o", color=col)
        #Linear regression plot:
        ax[i][j].plot([1,3],np.array([0,2])*slope+intercept, linestyle="--", linewidth=1, color=col)
        #Title:
        ax[i][j].set_title("slope={:.3f}, p={:.2e}".format(slope, max(pval, 1e-320)), fontsize=12)
for i in [1,2,4]:
    ax[i][0].set_ylabel("{0} expression (NPX)".format(["KNG1", "KLKB1",  "F12", "HPR", "APOL1"][i]))
for i in [0,3]:
    ax[i][0].set_ylabel("{0} expression (NPX):\nNot included in Olink".format(["KNG1", "KLKB1",  "F12", "HPR", "APOL1"][i]))
for j in range(5):
    ax[4][j].set_xticks([1,2,3], ["REF/REF", "REF/ALT", "ALT/ALT"], fontsize=11)
    ax[4][j].set_xlabel("{0} cis-variant".format(["KNG1", "KLKB1",  "F12", "HPR", "APOL1"][j]))
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/klkb1_etc_dot.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/klkb1_etc_dot.pdf", dpi=400)
plt.close()
