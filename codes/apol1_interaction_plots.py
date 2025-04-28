import pandas as pd
import numpy as np
import glob
import time as tm
from matplotlib import pyplot as plt
from scipy import stats
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
import seaborn as sns

#The interaction p-values as a heatmap:
pvals = pd.read_csv("apol1_cis_trans_interaction_pvals_eur.tsv", sep="\t", index_col=0)
pvals.index = ["APOL1\n(rs6000220)", "APOL1\n(rs148296684)", "APOL1\n(rs2239785)", 
               "F12\n(rs1801020)", "KLKB1\n(rs3733402)", "KNG1\n(rs5030062)", "KNG1\n(rs76438938)", 
               "HPR\n(rs116891509)", "HPR\n(rs217181)", "HPR\n(rs9930957)"]
pvals.columns = pvals.index
plt.figure(figsize=(10,10))
ax = sns.heatmap(-np.log10(pvals), cmap="copper", annot=pvals, fmt=".1e", vmin=0, vmax=15, linewidths=1, cbar_kws={'label': '-log10(p) for interaction effect\non APOL1 expression (NPX)', 'shrink': 0.75})
ax.set_facecolor('gray')
plt.xlabel("Gene\n(variant)")
plt.ylabel("Gene\n(variant)")
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.savefig("2way_intr_heatmap_full.png", dpi=400)
plt.savefig("2way_intr_heatmap_full.pdf", dpi=400)
plt.close()

#Simplified main figure version:
touse = ["KNG1\n(rs5030062)", "KLKB1\n(rs3733402)", "F12\n(rs1801020)", "HPR\n(rs9930957)", "APOL1\n(rs2239785)"]
pvals = pvals.loc[touse, touse]
pvals.index = ["KNG1", "KLKB1", "F12", "HPR", "APOL1"]
pvals.columns = pvals.index
plt.figure(figsize=(5,5))
ax = sns.heatmap(-np.log10(pvals), cmap="copper", annot=pvals, fmt=".1e", vmin=0, vmax=15, linewidths=1, cbar_kws={'label': '-log10(p) for interaction effect\non APOL1 expression (NPX)', 'shrink': 0.75})
ax.set_facecolor('gray')
plt.xlabel("Gene")
plt.ylabel("Gene")
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.savefig("2way_intr_heatmap_main.png", dpi=400)
plt.savefig("2way_intr_heatmap_main.pdf", dpi=400)
plt.close()


#And all the interaction 3d plots:
dos = pd.read_csv("cis_trans_interaction_variants_eur_apol1.tsv.gz", sep="\t", index_col=0)
cols = ["kl", "f12", "kng_intron", "v1", "hpr3"]
colnames = ["KLKB1", "F12", "KNG1", "APOL1", "HPR"]
snpnames = ["rs3733402", "rs1801020", "rs5030062", "rs2239785", "rs9930957"]
refs = ["G", "A", "A", "G", "C"]
alts = ["A", "G", "C", "A", "T"]
X = [[0,1,2],[0,1,2],[0,1,2]]
Y = [[0,0,0],[1,1,1],[2,2,2]]
for cnt1 in range(5):
    for cnt2 in range(5):
        mu = dos.groupby([cols[cnt1], cols[cnt2]]).apol1.mean().unstack().sort_index().T.sort_index().T
        sig = dos.groupby([cols[cnt1], cols[cnt2]]).apol1.sem().unstack().sort_index().T.sort_index().T
        mu_expected = mu*0
        for i in range(3):
            mu_expected.iloc[i,0] = mu.iloc[i,0]
            mu_expected.iloc[0,i] = mu.iloc[0,i]
        for i in range(1,3):
            for j in range(1,3):
                mu_expected.iloc[i,j] = mu.iloc[i,0] + mu.iloc[0,j] - mu.iloc[0,0]
        Z = mu.values
        Zexp = mu_expected.values
        sig_values = sig.values
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(projection='3d')
        ax.view_init(elev=20, azim=-50)
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))  # Set pane color to white with 0.0 alpha
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.xaxis.grid(False)
        ax.grid(color='black', linestyle='-', linewidth=0.5)  # Adjust parameters as needed
        ax.plot_wireframe(X, Y, Zexp, color="tab:grey", linestyle="--", linewidth=0.5)
        ax.plot_wireframe(X, Y, Z, color="tab:orange", linestyle="--", linewidth=0.5)
        ax.plot_surface(X, Y, Zexp, alpha=0.1, color="tab:grey")
        ax.plot_surface(X, Y, Z, alpha=0.1, color="tab:orange")
        ax.set_xticks([0,1,2])
        ax.set_xticklabels(["{0}/{0}".format(refs[cnt2]),"{0}/{1}".format(refs[cnt2], alts[cnt2]),"{0}/{0}".format(alts[cnt2])])
        ax.set_yticks([0,1,2])
        ax.set_yticklabels(["{0}/{0}".format(refs[cnt1]),"{0}/{1}".format(refs[cnt1], alts[cnt1]),"{0}/{0}".format(alts[cnt1])])
        ax.set_xlabel("{0} cis-variant\n({1}, REF={2})".format(colnames[cnt2], snpnames[cnt2], refs[cnt2]))        
        ax.set_ylabel("{0} cis-variant\n({1}, REF={2})".format(colnames[cnt1], snpnames[cnt1], refs[cnt1]))
        ax.set_zlabel("APOL1 expression (Olink NPX)")
        for i in range(len(X)):
            for j in range(len(X[i])):
                ax.plot([X[i][j], X[i][j]], [Y[i][j], Y[i][j]], [Z[i][j] + sig_values[i][j], Z[i][j] - sig_values[i][j]], marker="_", color="tab:orange")
                ax.scatter(X[i][j], Y[i][j], Zexp[i][j], color="tab:grey", marker=".")
                ax.scatter(X[i][j], Y[i][j], Z[i][j], color="tab:orange")
        #Legend:
        i, j = 0,0
        ax.scatter(X[i][j], Y[i][j], Zexp[i][j], color="tab:grey", marker=".", label="Expected (additive)")
        ax.scatter(X[i][j], Y[i][j], Z[i][j], color="tab:orange", label="Observed")
        plt.gca().tick_params(axis='x', pad=0.3)
        plt.gca().tick_params(axis='y', pad=0.3)
        plt.legend()
        plt.savefig("/Users/qsw/Desktop/misc_plots/3d_plots_update1129/3d_plot_{0}_vs_{1}.png".format(colnames[cnt1], colnames[cnt2]), dpi=400)
        plt.savefig("/Users/qsw/Desktop/misc_plots/3d_plots_update1129/3d_plot_{0}_vs_{1}.pdf".format(colnames[cnt1], colnames[cnt2]), dpi=400)
        plt.close()