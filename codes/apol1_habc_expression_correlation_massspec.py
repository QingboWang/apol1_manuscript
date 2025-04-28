import pandas as pd
import numpy as np
import glob
import time as tm
from matplotlib import pyplot as plt
from scipy import stats
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

#Load healthABC data:
prot = pd.read_csv("~/Desktop/habc/proteomics_cape_export.v1.tsv.gz", sep="\t", compression="gzip")
prot.set_index("Gene", inplace=True)
prot = prot.loc[prot.index.intersection(["APOL1", "HP", "HPR","KLKB1","F12", "KNG1", "APOA1"]), :] #The guys we care about
prot.index.unique() #They actually exist!

#Make it a matrix per timepoint, individual, gene and mean of timepoint:
prot["prot_id"] = prot.index + "_" + prot.Protein
mat1 = prot[prot.timepoint==1].groupby(["prot_id", "HABCID"])["value"].mean().unstack(level=0)
mat2 = prot[prot.timepoint==2].groupby(["prot_id", "HABCID"])["value"].mean().unstack(level=0)
mat12 = prot.groupby(["prot_id", "HABCID"])["value"].mean().unstack(level=0) #Ave, accounting for NAs
mat12min = prot.groupby(["prot_id", "HABCID"])["value"].min().unstack(level=0) #for later error bar
mat12max = prot.groupby(["prot_id", "HABCID"])["value"].max().unstack(level=0) #for later error bar

#Take the correlation of the key genes of interest and plot:
colors = ["tab:olive", "tab:green", "tab:cyan", "darkgreen", "tab:blue", "tab:purple", "tab:purple"]
y = mat12["APOL1_tr|A5PL32|A5PL32_HUMAN"]
ymin = mat12min["APOL1_tr|A5PL32|A5PL32_HUMAN"]
ymax = mat12max["APOL1_tr|A5PL32|A5PL32_HUMAN"]
xtextpos = [-2,-3,1,-4,-1.5, -2, -0.75]
for i, gene in enumerate(["F12_sp|P00748|FA12_HUMAN", "HPR_sp|P00739|HPTR_HUMAN", "HP_sp|P00738|HPT_HUMAN", "HP_tr|H0Y300|H0Y300_HUMAN", "KLKB1_tr|H0YAC1|H0YAC1_HUMAN", "KNG1_tr|D3DNU8|D3DNU8_HUMAN", "KNG1_sp|P01042|KNG1_HUMAN"]):
    x = mat12[gene]
    xmin = mat12min[gene]
    xmax = mat12max[gene]
    plt.figure(figsize=(4.5,4.5))
    plt.axvline(x=0, color="black", linestyle="--", linewidth=0.5, zorder=-2)
    plt.axhline(y=0, color="black", linestyle="--", linewidth=0.5, zorder=-2)
    plt.scatter(x, y, label=gene.split("_")[0], color=colors[i], alpha=0.5)
    plt.title("Protein expression in HealthABC (Mass Spec.)", fontsize=12)
    plt.ylabel("APOL1 (A5PL32)")
    plt.xlabel(gene.split("_")[0] + " (" + gene.split("|")[1] + ")")
    plt.text(x=xtextpos[i], y=1.05, s="r={0:.3f}".format(corr.loc["APOL1_tr|A5PL32|A5PL32_HUMAN", gene]), fontsize=15)
    plt.tight_layout()
    plt.savefig("/Users/qsw/Desktop/misc_plots/habc_protein_expression_{0}vs{1}.png".format("APOL1", gene), dpi=400)
    plt.savefig("/Users/qsw/Desktop/misc_plots/habc_protein_expression_{0}vs{1}.pdf".format("APOL1", gene), dpi=400)
    plt.close()

#Also p-value for correlation:
y = mat12["APOL1_tr|A5PL32|A5PL32_HUMAN"]
for i, gene in enumerate(["F12_sp|P00748|FA12_HUMAN", "HPR_sp|P00739|HPTR_HUMAN", "HP_sp|P00738|HPT_HUMAN", "HP_tr|H0Y300|H0Y300_HUMAN", "KLKB1_tr|H0YAC1|H0YAC1_HUMAN", "KNG1_tr|D3DNU8|D3DNU8_HUMAN", "KNG1_sp|P01042|KNG1_HUMAN"]):
    x = mat12[gene]
    isna = (x.isna() | y.isna())
    x1 = x[~isna]
    y1 = y[~isna]
    print (gene)
    print (stats.pearsonr(x1, y1))