#Summary of the manhattan plots we generated

import pandas as pd
import numpy as np
import glob
import time as tm
from matplotlib import pyplot as plt
from scipy import stats
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"


#The AFR manhattan plot highlighting G1 and G2:
ss = pd.read_parquet("apol1_geno_ex_hybrid_nominal_n947.cis_qtl_pairs.chr22.parquet")
ssc1 = pd.read_parquet("apol1_geno_ex_hybrid_g2adj_n947.cis_qtl_pairs.chr22.parquet")
ssc2 = pd.read_parquet("apol1_geno_ex_hybrid_g1g2adj_n947.cis_qtl_pairs.chr22.parquet")
fm = pd.read_csv("apol1_ukb_afr_geno_ex_hybrid_n947_cis.snp", sep=" ")
sus = pd.read_csv("ppp_apol1_afr_geno_ex_hybrid_susie_pip.tsv", sep="\t")
fm.set_index("rsid", inplace=True)
sus.set_index("rsid", inplace=True)
fm = fm.join(sus, how="left")
fm["min_pip"] = np.minimum(fm.prob, fm.pip)
ss.set_index("variant_id", inplace=True)
fm = fm.join(ss.start_distance, how="left")
ss.reset_index(inplace=True)    
fm = fm[(-0.2*10**5<fm.start_distance)&(fm.start_distance<0.4*10**5)]
from matplotlib.cm import viridis
colors = ["#000000", "#555555", "#777777"]
ytextposs = [17, 9.5]
texts = ["Conditioning on G2", "Conditioning on G2 and G1"]
topylims = [151, 1.3, 0, 22, 12]#manual
sss = [ss[(-0.2*10**5<ss.start_distance)&(ss.start_distance<0.4*10**5)], ssc1[(-0.2*10**5<ssc1.start_distance)&(ssc1.start_distance<0.4*10**5)],
       ssc2[(-0.2*10**5<ssc2.start_distance)&(ssc2.start_distance<0.4*10**5)]]
fig, ax = plt.subplots(nrows=5, ncols=1, figsize=(10,5), sharex=False, sharey=False, gridspec_kw={'height_ratios': [1.5, 0.5, 0.2, 1, 1]})
i = 0 #First main plot
ax[i].scatter(sss[i].start_distance, -np.log10(sss[i].pval_nominal), color=colors[i])
ax[i].scatter(sss[i][sss[i].variant_id=="22:36265995:AATAATT:A"].start_distance, -np.log10(sss[i][sss[i].variant_id=="22:36265995:AATAATT:A"].pval_nominal), color="tab:purple", label="G2")
ax[i].scatter(sss[i][sss[i].variant_id=="22:36265988:T:G"].start_distance, -np.log10(sss[i][sss[i].variant_id=="22:36265988:T:G"].pval_nominal), color="tab:red", label="G1M")
ax[i].scatter(sss[i][sss[i].variant_id=="22:36265860:A:G"].start_distance, -np.log10(sss[i][sss[i].variant_id=="22:36265860:A:G"].pval_nominal), color="tab:orange", label="G1G")    
ax[i].axhline(y=-np.log10(5*10**-8), linestyle="--", linewidth=1, color="tab:blue")
i = 1 #PIP - we plot the SuSiE version for main fig.
ax[i].scatter(fm.start_distance, fm.min_pip, color="navy")
ax[i].scatter(fm[fm.index=="22:36265995:AATAATT:A"].start_distance, fm[fm.index=="22:36265995:AATAATT:A"].pip, color="tab:purple")
ax[i].scatter(fm[fm.index=="22:36265988:T:G"].start_distance, fm[fm.index=="22:36265988:T:G"].min_pip, color="tab:red")
ax[i].scatter(fm[fm.index=="22:36265860:A:G"].start_distance, fm[fm.index=="22:36265860:A:G"].min_pip, color="tab:orange")
for i in [3,4]:
        ax[i].scatter(sss[i-2].start_distance, -np.log10(sss[i-2].pval_nominal), color=colors[i-2])
        ax[i].scatter(sss[i-2][sss[i-2].variant_id=="22:36265995:AATAATT:A"].start_distance, 
                      -np.log10(sss[i-2][sss[i-2].variant_id=="22:36265995:AATAATT:A"].pval_nominal), color="tab:purple")
        ax[i].scatter(sss[i-2][sss[i-2].variant_id=="22:36265988:T:G"].start_distance, 
                      -np.log10(sss[i-2][sss[i-2].variant_id=="22:36265988:T:G"].pval_nominal), color="tab:red")
        ax[i].scatter(sss[i-2][sss[i-2].variant_id=="22:36265860:A:G"].start_distance, 
                      -np.log10(sss[i-2][sss[i-2].variant_id=="22:36265860:A:G"].pval_nominal), color="tab:orange")
        ax[i].axhline(y=-np.log10(5*10**-8), linestyle="--", linewidth=1, color="tab:blue")
        ax[i].text(x=-20000, y=ytextposs[i-3], s=texts[i-3], ha="left", va="center", color=colors[i-2])
ax[0].set_title("cis-pQTL effect on APOL1 (Olink) in UKB AFR")
ax[0].set_ylabel("-log10(p)")
ax[1].set_ylabel("PIP")
ax[3].set_ylabel("-log10(p)")
ax[4].set_ylabel("-log10(p)")
ax[4].set_xlabel("Distance to TSS of APOL1")
for cnt, a in enumerate(ax.flat):
    ax[cnt].set_ylim(top=topylims[cnt])
    ax[cnt].set_xlim([-0.21*10**5, 0.41*10**5])
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)    
    if cnt==2:
        a.spines['bottom'].set_visible(False)
        a.spines['left'].set_visible(False)
        a.set_yticks([])
    if cnt!=4:
         a.set_xticks([])
#fig.text(0, 0.55, ', -log10(p)', va='center', rotation='vertical', fontsize=12) #Shared y axis    
ax[0].legend()
#ax[0].set_ylim([-1,80])
plt.tight_layout()
fig.subplots_adjust(hspace=0.2)
plt.savefig("afr_manhattan_main.png", dpi=500)
plt.savefig("afr_manhattan_main.pdf", dpi=500)
plt.show()
plt.close()


#G1 vs G2 plots:
d_g12_geno = pd.read_csv("ukb_exome_g1g2v1_dosages.tsv.gz", sep="\t", index_col=0)
apol1_afr = pd.read_csv("ppp_afr_expression_mat_naincluded.tsv.gz", sep="\t", index_col=0)
apol1_afr = apol1_afr.OID30708.squeeze()
apol1_afr.dropna(inplace=True)
d_g12 = d_g12_geno.join(apol1_afr, how="inner")
d_g12.rename(columns={"OID30708":"apol1"}, inplace=True)
from sklearn.linear_model import LinearRegression
np.random.seed(3)
d_g12["rd"] = np.random.rand(d_g12.shape[0])
#colors = ["tab:grey", "tab:orange", "tab:red"]
from matplotlib.cm import viridis
colors = [viridis(0.3), viridis(0.6), viridis(0.9)]
labels = ["REF/REF","REF/ALT", "ALT/ALT"]
d_g12 = d_g12[(d_g12.G1M==d_g12.G1G)]
#d_g12 = d_g12[(d_g12.V1==0)]
#Plot:
plt.figure(figsize=(7.5,3))
for i in range(3):
    for j in range(3):
        d = d_g12[(d_g12.G1M==i)&(d_g12.G2==j)]
        if i==0:
            plt.scatter((d.rd-0.5)/5+i+(j-1)/4, d.apol1, color=colors[j], label=labels[j], alpha=0.75)
        else:
            plt.scatter((d.rd-0.5)/5+i+(j-1)/4, d.apol1, color=colors[j], alpha=0.75)
        #mu = d.apol1.mean()
        #plt.plot([i+(j-1)/4 - 1/5, i+(j-1)/4 + 1/5], [mu,mu], color=colors[j], linestyle="--", linewidth=2)        
#Regression line:
for i in range(2):
    model = LinearRegression()
    model.fit(pd.DataFrame(d_g12[(d_g12.G1M==i)].G2), d_g12[(d_g12.G1M==i)].apol1)
    slope = model.coef_[0]
    intercept = model.intercept_
    plt.plot([i-1/4,i+1/4],np.array([0,2])*slope+intercept, linestyle="--", linewidth=2, color="black")
    plt.text(i-1/4, 4.3, f"$\\beta$={slope:.3f}")
plt.legend(title="G2 genotype", bbox_to_anchor=(1.05,1))
plt.xlabel("G1 genotype")
plt.ylabel("APOL1 (Olink NPX)")
plt.xticks([0,1,2],labels)
plt.ylim([-0.9,4.9])
plt.tight_layout()
plt.savefig("g1g2_intr.png", dpi=500)
plt.savefig("g1g2_intr.pdf", dpi=500)
plt.show()
plt.close()

#The other side:
from matplotlib.cm import plasma
colors = [plasma(0.1), plasma(0.3), plasma(0.6)]
plt.figure(figsize=(7.5,3))
for i in range(3):
    for j in range(3):
        d = d_g12[(d_g12.G1M==j)&(d_g12.G2==i)]
        if i==0:
            plt.scatter((d.rd-0.5)/5+i+(j-1)/4, d.apol1, color=colors[j], label=labels[j], alpha=0.75)
        else:
            plt.scatter((d.rd-0.5)/5+i+(j-1)/4, d.apol1, color=colors[j], alpha=0.75)
#Regression line:
for i in range(2):
    model = LinearRegression()
    model.fit(pd.DataFrame(d_g12[(d_g12.G2==i)].G1M), d_g12[(d_g12.G2==i)].apol1)
    slope = model.coef_[0]
    intercept = model.intercept_
    plt.plot([i-1/4,i+1/4],np.array([0,2])*slope+intercept, linestyle="--", linewidth=2, color="black")            
    plt.text(i-1/4, 4.3, f"$\\beta$={slope:.3f}")
plt.legend(title="G1 genotype", bbox_to_anchor=(1.05,1))
plt.xlabel("G2 genotype")
plt.ylabel("APOL1 (Olink NPX)")
plt.xticks([0,1,2],labels)
plt.ylim([-0.9,4.9])
plt.tight_layout()
plt.savefig("g1g2_intr2.png", dpi=500)
plt.savefig("g1g2_intr2.pdf", dpi=500)
plt.show()
plt.close()

#Manhattan and fine-mapping results for EUR:
ss = pd.read_parquet("156_APOL1_OID30708.cis_qtl_pairs.chr22.parquet")
fm = pd.read_csv("156_APOL1_OID30708.snp", sep=" ")
sus = pd.read_csv("156_APOL1_OID30708_susie_pip.txt", sep="\t")
ss["abs_z"] = abs(ss.slope/ss.slope_se)
ss.set_index("variant_id", inplace=True)
fm.set_index("rsid", inplace=True)
sus.set_index("rsid", inplace=True)
fm = fm.join(sus, how="left")
fm["min_pip"] = np.minimum(fm.prob, fm.pip)
fm.index = fm.index.str.replace("APOL1_OID30708_","")
fm = fm.join(ss[["start_distance","abs_z"]], how="left")
ss.reset_index(inplace=True)    
ss2 = pd.read_parquet("apol1_eur_v1adj.cis_qtl_pairs.chr22.parquet")
ss2["abs_z"] = abs(ss2.slope/ss2.slope_se)
ss2.set_index("variant_id", inplace=True)
ss3 = pd.read_parquet("apol1_eur_v1v2adj.cis_qtl_pairs.chr22.parquet")
ss3["abs_z"] = abs(ss3.slope/ss3.slope_se)
ss3.set_index("variant_id", inplace=True)
ss4 = pd.read_parquet("apol1_eur_v1v2r1adj.cis_qtl_pairs.chr22.parquet")
ss4["abs_z"] = abs(ss4.slope/ss4.slope_se)
ss4.set_index("variant_id", inplace=True)

fm = fm[(-0.2*10**5<fm.start_distance)&(fm.start_distance<0.4*10**5)]
ss2 = ss2[(-0.2*10**5<ss2.start_distance)&(ss2.start_distance<0.4*10**5)]
ss3 = ss3[(-0.2*10**5<ss3.start_distance)&(ss3.start_distance<0.4*10**5)]
ss4 = ss4[(-0.2*10**5<ss4.start_distance)&(ss4.start_distance<0.4*10**5)]
topylims = [45, 1.3, 0, 28, 15]#manual
fig, ax = plt.subplots(nrows=5, ncols=1, figsize=(10,5), sharex=False, sharey=False, gridspec_kw={'height_ratios': [1.5, 0.5, 0.2, 1, 1]})
i = 0 #First main plot
ax[i].scatter(fm.start_distance, fm.abs_z, color="black")
ax[i].scatter(fm[fm.index.str.startswith("chr22_36265284_G_A")].start_distance, fm[fm.index.str.startswith("chr22_36265284_G_A")].abs_z, color="tab:green", label="rs2239785 (missense)")
ax[i].scatter(fm[fm.index.str.startswith("chr22_36253920_C_T")].start_distance, fm[fm.index.str.startswith("chr22_36253920_C_T")].abs_z, color="tab:cyan", label="rs6000220 (5'UTR)")
ax[i].axhline(y=5.327, linestyle="--", linewidth=1, color="tab:blue")
i = 1 #PIP - we plot the SuSiE version for main fig.
ax[i].scatter(fm.start_distance, fm.pip, color="navy")
ax[i].scatter(fm[fm.index.str.startswith("chr22_36265284_G_A")].start_distance, fm[fm.index.str.startswith("chr22_36265284_G_A")].pip, color="tab:green")
ax[i].scatter(fm[fm.index.str.startswith("chr22_36253920_C_T")].start_distance, fm[fm.index.str.startswith("chr22_36253920_C_T")].pip, color="tab:cyan")
i = 3 #conditional
ax[i].scatter(ss2.start_distance, ss2.abs_z, color="#444444")
ax[i].scatter(ss2[ss2.index.str.startswith("chr22_36265284_G_A")].start_distance, ss2[ss2.index.str.startswith("chr22_36265284_G_A")].abs_z, color="tab:green")
ax[i].scatter(ss2[ss2.index.str.startswith("chr22_36253920_C_T")].start_distance, ss2[ss2.index.str.startswith("chr22_36253920_C_T")].abs_z, color="tab:cyan")
ax[i].axhline(y=5.327, linestyle="--", linewidth=1, color="tab:blue")
i = 4 #conditional 2
ax[i].scatter(ss3.start_distance, ss3.abs_z, color="#777777")
ax[i].scatter(ss3[ss3.index.str.startswith("chr22_36265284_G_A")].start_distance, ss3[ss3.index.str.startswith("chr22_36265284_G_A")].abs_z, color="tab:green")
ax[i].scatter(ss3[ss3.index.str.startswith("chr22_36253920_C_T")].start_distance, ss3[ss3.index.str.startswith("chr22_36253920_C_T")].abs_z, color="tab:cyan")
ax[i].axhline(y=5.327, linestyle="--", linewidth=1, color="tab:blue")
#i=5 #conditional 3 - didn't really help. It is a rare variant.
ax[0].set_title("cis-pQTL effect on APOL1 (Olink) in UKB EUR", fontsize=12)
ax[0].set_ylabel("|z-score|")
ax[1].set_ylabel("PIP")
ax[3].set_ylabel("|z-score|")
ax[4].set_ylabel("|z-score|")
ax[4].set_xlabel("Distance to TSS of APOL1")
ax[3].text(x=-20000, y=25, s="Conditioning on rs2239785", ha="left", va="center", color="#444444")
ax[4].text(x=-20000, y=13, s="Conditioning on rs2239785 and rs6000220", ha="left", va="center", color="#777777")
for cnt, a in enumerate(ax.flat):
    ax[cnt].set_ylim(top=topylims[cnt])
    ax[cnt].set_xlim([-0.21*10**5, 0.41*10**5])
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)    
    if cnt==2:
        a.spines['bottom'].set_visible(False)
        a.spines['left'].set_visible(False)
        a.set_yticks([])
    if cnt!=4:
         a.set_xticks([])
#fig.text(0, 0.55, ', -log10(p)', va='center', rotation='vertical', fontsize=12) #Shared y axis    
ax[0].legend()
plt.tight_layout()
fig.subplots_adjust(hspace=0.2)
plt.savefig("eur_manhattan_main.png", dpi=500)
plt.savefig("eur_manhattan_main.pdf", dpi=500)
plt.show()
plt.close()



#Trans players:

#KLKB1:
ss = pd.read_csv("apol1_eur_trans_from_KLKB1.tsv.gz", sep="\t", index_col=0)
fm = pd.read_csv("KLKB1_trans_apol1.snp", sep=" ", index_col=0)
sus = pd.read_csv("KLKB1_trans_on_apol1_susie_pip.txt", sep="\t", index_col=0)
ss.set_index("variant_id", inplace=True)
fm.set_index("rsid", inplace=True)
sus.set_index("rsid", inplace=True)
fm.index = fm.index.str.split("OID30708_").str[1]
sus.index = sus.index.str.split("OID30708_").str[1]
ss = ss.join(fm.prob, how="left")
ss = ss.join(sus.pip, how="left")
ss.rename(columns={"prob":"pip_fm", "pip":"pip_sus"}, inplace=True)
ss["min_pip"] = np.minimum(ss.pip_fm, ss.pip_sus)
ss["position"] = ss.index.str.split("_").str[1].astype(int)
tss = 186208979
ss["dist_tss"] = ss.position - tss
ss = ss[(-0.5*10**5<ss.dist_tss)&(ss.dist_tss<1.5*10**5)]
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,3), sharex=True, sharey=False, gridspec_kw={'height_ratios': [1.5, 0.5]})
i = 0 #First main plot
ax[i].scatter(ss.dist_tss, -np.log10(ss.pval), color="black", edgecolors="#ffffff", linewidth=0.1)
ax[i].scatter(ss[ss.position==186236880].dist_tss, -np.log10(ss[ss.position==186236880].pval), color="tab:green", label="rs3733402\n(missense)")
ax[i].axhline(y=-np.log10(5*10**-8), linestyle="--", linewidth=1, color="tab:blue")
i = 1 #PIP
ax[i].scatter(ss.dist_tss, ss.min_pip, color="navy", edgecolors="#ffffff", linewidth=0.1)
ax[i].scatter(ss[ss.position==186236880].dist_tss, ss[ss.position==186236880].min_pip, color="tab:green")
ax[0].set_title("trans-pQTL effect of KLKB1 cis-variant on APOL1 (Olink) in UKB EUR", fontsize=12)
ax[0].set_ylabel("-log10(p)")
ax[1].set_ylabel("PIP")
ax[1].set_xlabel("Distance to TSS of KLKB1")
for cnt, a in enumerate(ax.flat):
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)    
ax[0].legend()
ax[0].set_ylim([-5,200])
ax[1].set_ylim([-0.01, 0.2])
ax[1].set_xlim([-0.51*10**5, 1.51*10**5])
plt.tight_layout()
fig.subplots_adjust(hspace=0.2)
plt.savefig("kng1_manhattan_main.png", dpi=500)
plt.savefig("kng1_manhattan_main.pdf", dpi=500)
plt.show()
plt.close()

#F12:
ss = pd.read_csv("apol1_eur_trans_from_F12.tsv.gz", sep="\t", index_col=0)
fm = pd.read_csv("F12_trans_apol1.snp", sep=" ", index_col=0)
sus = pd.read_csv("F12_trans_on_apol1_susie_pip.txt", sep="\t", index_col=0)
ss.set_index("variant_id", inplace=True)
fm.set_index("rsid", inplace=True)
sus.set_index("rsid", inplace=True)
fm.index = fm.index.str.split("OID30708_").str[1]
sus.index = sus.index.str.split("OID30708_").str[1]
ss = ss.join(fm.prob, how="left")
ss = ss.join(sus.pip, how="left")
ss.rename(columns={"prob":"pip_fm", "pip":"pip_sus"}, inplace=True)
ss["min_pip"] = np.minimum(ss.pip_fm, ss.pip_sus)
ss["position"] = ss.index.str.split("_").str[1].astype(int)
tss = 177416583
ss["dist_tss"] = ss.position - tss
ss = ss[(-10**5<ss.dist_tss)&(ss.dist_tss<10**5)]
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,3), sharex=True, sharey=False, gridspec_kw={'height_ratios': [1.5, 0.5]})
i = 0 #First main plot
ax[i].scatter(ss.dist_tss, -np.log10(ss.pval), color="black", edgecolors="#ffffff", linewidth=0.1)
ax[i].scatter(ss[ss.position==177409531].dist_tss, -np.log10(ss[ss.position==177409531].pval), color="tab:green", label="rs1801020\n(5'UTR)")
ax[i].axhline(y=-np.log10(5*10**-8), linestyle="--", linewidth=1, color="tab:blue")
i = 1 #PIP
ax[i].scatter(ss.dist_tss, ss.min_pip, color="navy", edgecolors="#ffffff", linewidth=0.1)
ax[i].scatter(ss[ss.position==177409531].dist_tss, ss[ss.position==177409531].min_pip, color="tab:green")
ax[0].set_title("trans-pQTL effect of F12 cis-variant on APOL1 (Olink) in UKB EUR", fontsize=12)
ax[0].set_ylabel("-log10(p)")
ax[1].set_ylabel("PIP")
ax[1].set_xlabel("Distance to TSS of F12")
for cnt, a in enumerate(ax.flat):
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)    
ax[0].legend()
ax[0].set_ylim([-3,150])
ax[1].set_ylim([-0.02, 1.0])
ax[1].set_xlim([-1.01*10**5, 1.01*10**5])
plt.tight_layout()
fig.subplots_adjust(hspace=0.2)
plt.savefig("f12_manhattan_main.png", dpi=500)
plt.savefig("f12_manhattan_main.pdf", dpi=500)
plt.show()
plt.close()


#KNG1:
ss = pd.read_csv("apol1_eur_trans_from_KNG1.tsv.gz", sep="\t", index_col=0)
fm = pd.read_csv("KNG1_trans_apol1.snp", sep=" ", index_col=0)
sus = pd.read_csv("KNG1_trans_on_apol1_susie_pip.txt", sep="\t", index_col=0)
ss.set_index("variant_id", inplace=True)
fm.set_index("rsid", inplace=True)
sus.set_index("rsid", inplace=True)
fm.index = fm.index.str.split("OID30708_").str[1]
sus.index = sus.index.str.split("OID30708_").str[1]
ss = ss.join(fm.prob, how="left")
ss = ss.join(sus.pip, how="left")
ss.rename(columns={"prob":"pip_fm", "pip":"pip_sus"}, inplace=True)
ss["min_pip"] = np.minimum(ss.pip_fm, ss.pip_sus)
ss["position"] = ss.index.str.split("_").str[1].astype(int)
tss = 186717348
ss["dist_tss"] = ss.position - tss
ss = ss[(-10**5<ss.dist_tss)&(ss.dist_tss<1.3*10**5)]
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,3), sharex=True, sharey=False, gridspec_kw={'height_ratios': [1.5, 0.5]})
i = 0 #First main plot
ax[i].scatter(ss.dist_tss, -np.log10(ss.pval), color="black", edgecolors="#ffffff", linewidth=0.1)
ax[i].scatter(ss[ss.position==186736391].dist_tss, -np.log10(ss[ss.position==186736391].pval), color="tab:green", label="rs5030062\n(intronic)")
ax[i].scatter(ss[ss.position==186743735].dist_tss, -np.log10(ss[ss.position==186743735].pval), color="tab:orange", label="rs76438938\n(3'UTR)")
ax[i].axhline(y=-np.log10(5*10**-8), linestyle="--", linewidth=1, color="tab:blue")
i = 1 #PIP
ax[i].scatter(ss.dist_tss, ss.min_pip, color="navy", edgecolors="#ffffff", linewidth=0.1)
ax[i].scatter(ss[ss.position==186736391].dist_tss, ss[ss.position==186736391].min_pip, color="tab:green")
ax[i].scatter(ss[ss.position==186743735].dist_tss, ss[ss.position==186743735].min_pip, color="tab:orange", label="rs76438938\n(3'UTR)")
ax[0].set_title("trans-pQTL effect of KNG1 cis-variant on APOL1 (Olink) in UKB EUR", fontsize=12)
ax[0].set_ylabel("-log10(p)")
ax[1].set_ylabel("PIP")
ax[1].set_xlabel("Distance to TSS of KNG1")
for cnt, a in enumerate(ax.flat):
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)    
ax[0].legend(loc="upper right")
ax[0].set_ylim([-1,32])
ax[1].set_ylim([-0.02, 1.0])
ax[1].set_xlim([-1.01*10**5, 1.31*10**5])
plt.tight_layout()
fig.subplots_adjust(hspace=0.2)
plt.savefig("kng1_manhattan_main.png", dpi=500)
plt.savefig("kng1_manhattan_main.pdf", dpi=500)
plt.show()
plt.close()

#HPR EUR:
ss = pd.read_csv("apol1_eur_trans_from_HPR.tsv.gz", sep="\t", index_col=0)
fm = pd.read_csv("HPR_trans_apol1.snp", sep=" ", index_col=0)
sus = pd.read_csv("HPR_trans_on_apol1_susie_pip.txt", sep="\t", index_col=0)
ss.set_index("variant_id", inplace=True)
fm.set_index("rsid", inplace=True)
sus.set_index("rsid", inplace=True)
fm.index = fm.index.str.split("OID30708_").str[1]
sus.index = sus.index.str.split("OID30708_").str[1]
ss = ss.join(fm.prob, how="left")
ss = ss.join(sus.pip, how="left")
ss.rename(columns={"prob":"pip_fm", "pip":"pip_sus"}, inplace=True)
ss["min_pip"] = np.minimum(ss.pip_fm, ss.pip_sus)
ss["position"] = ss.index.str.split("_").str[1].astype(int)
tss = 72063148
ss["dist_tss"] = ss.position - tss
ss["abs_z"] = abs(ss.b/ss.b_se)
#ss = ss[(-10**5<ss.dist_tss)&(ss.dist_tss<1.3*10**5)]
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(11,3), sharex=True, sharey=False, gridspec_kw={'height_ratios': [1.5, 0.5]})
i = 0 #First main plot
ax[i].scatter(ss.dist_tss, ss.abs_z, color="black", edgecolors="#ffffff", linewidth=0.1)
ax[i].scatter(ss[ss.position==72080103].dist_tss, ss[ss.position==72080103].abs_z, color="tab:green", label="rs217181\n(intronic)")
ax[i].scatter(ss[ss.position==72116024].dist_tss, ss[ss.position==72116024].abs_z, color="tab:orange", label="rs9930957\n(intronic)")
ax[i].scatter(ss[ss.position==72071661].dist_tss, ss[ss.position==72071661].abs_z, color="tab:cyan", label="rs116891509\n(intronic)")
ax[i].axhline(y=5.327, linestyle="--", linewidth=1, color="tab:blue")
i = 1 #PIP
ax[i].scatter(ss.dist_tss, ss.min_pip, color="navy", edgecolors="#ffffff", linewidth=0.1)
ax[i].scatter(ss[ss.position==72080103].dist_tss, ss[ss.position==72080103].min_pip, color="tab:green")
ax[i].scatter(ss[ss.position==72116024].dist_tss, ss[ss.position==72116024].min_pip, color="tab:orange")
ax[i].scatter(ss[ss.position==72071661].dist_tss, ss[ss.position==72071661].min_pip, color="tab:cyan")
ax[0].set_title("trans-pQTL effect of HPR cis-variant on APOL1 (Olink) in UKB EUR", fontsize=12)
ax[0].set_ylabel("|z-score|")
ax[1].set_ylabel("PIP")
ax[1].set_xlabel("Distance to TSS of HPR")
for cnt, a in enumerate(ax.flat):
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)    
ax[0].legend(bbox_to_anchor=(1.02,1))
ax[0].set_ylim([-1,70])
ax[1].set_ylim([-0.05, 1.15])
ax[1].set_xlim([-2.05*10**5, 2.2*10**5])
plt.tight_layout()
fig.subplots_adjust(hspace=0.2)
plt.savefig("hpr_eur_manhattan_main.png", dpi=500)
plt.savefig("hpr_eur_manhattan_main.pdf", dpi=500)
plt.show()
plt.close()

#pval version:
ss["mlog10p"] = -np.log10(ss.pval.replace(0,1e-320))
#ss = ss[(-10**5<ss.dist_tss)&(ss.dist_tss<1.3*10**5)]
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(11,3), sharex=True, sharey=False, gridspec_kw={'height_ratios': [1.5, 0.5]})
i = 0 #First main plot
ax[i].scatter(ss.dist_tss, ss.mlog10p, color="black", edgecolors="#ffffff", linewidth=0.1)
ax[i].scatter(ss[ss.position==72080103].dist_tss, ss[ss.position==72080103].mlog10p, color="tab:green", label="rs217181\n(intronic)")
ax[i].scatter(ss[ss.position==72116024].dist_tss, ss[ss.position==72116024].mlog10p, color="tab:orange", label="rs9930957\n(intronic)")
ax[i].scatter(ss[ss.position==72071661].dist_tss, ss[ss.position==72071661].mlog10p, color="tab:cyan", label="rs116891509\n(intronic)")
ax[i].axhline(y=-np.log10(5e-8), linestyle="--", linewidth=1, color="tab:blue")
i = 1 #PIP
ax[i].scatter(ss.dist_tss, ss.min_pip, color="navy", edgecolors="#ffffff", linewidth=0.1)
ax[i].scatter(ss[ss.position==72080103].dist_tss, ss[ss.position==72080103].min_pip, color="tab:green")
ax[i].scatter(ss[ss.position==72116024].dist_tss, ss[ss.position==72116024].min_pip, color="tab:orange")
ax[i].scatter(ss[ss.position==72071661].dist_tss, ss[ss.position==72071661].min_pip, color="tab:cyan")
ax[0].set_title("trans-pQTL effect of HPR cis-variant on APOL1 (Olink) in UKB EUR", fontsize=12)
ax[0].set_ylabel("-log10(p)")
ax[1].set_ylabel("PIP")
ax[1].set_xlabel("Distance to TSS of HPR")
for cnt, a in enumerate(ax.flat):
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)    
ax[0].legend(bbox_to_anchor=(1.02,1))
#ax[0].set_ylim([-1,70])
ax[1].set_ylim([-0.05, 1.15])
ax[1].set_xlim([-2.05*10**5, 2.2*10**5])
plt.tight_layout()
fig.subplots_adjust(hspace=0.2)
plt.savefig("hpr_eur_manhattan_pval.png", dpi=500)
plt.savefig("hpr_eur_manhattan_pval.pdf", dpi=500)
plt.show()
plt.close()


#HPR AFR:
ss = pd.read_csv("apol1_afr_trans_from_HPR.tsv.gz", sep="\t", index_col=0)
fm = pd.read_csv("HPR_trans_apol1_afr.snp", sep=" ", index_col=0)
sus = pd.read_csv("HPR_trans_on_apol1_afr_susie_pip.txt", sep="\t", index_col=0)
ss.set_index("variant_id", inplace=True)
fm.set_index("rsid", inplace=True)
sus.set_index("rsid", inplace=True)
ss.index = ss.index.str.replace("cnv","dummy_72063148")
fm.index = fm.index.str.replace("chr16_1_cnv_cnv","dummy_72063148")
sus.index = sus.index.str.replace("chr16_1_cnv_cnv","dummy_72063148")
fm.index = fm.index.str.split("OID30708_").str[1]
sus.index = sus.index.str.split("OID30708_").str[1]
ss = ss.join(fm.prob, how="left")
ss = ss.join(sus.pip, how="left")
ss.rename(columns={"prob":"pip_fm", "pip":"pip_sus"}, inplace=True)
ss["min_pip"] = np.minimum(ss.pip_fm, ss.pip_sus)
ss["position"] = ss.index.str.split("_").str[1].astype(int)
tss = 72063148
ss["dist_tss"] = ss.position - tss
#ss = ss[(-10**5<ss.dist_tss)&(ss.dist_tss<1.3*10**5)]
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(11,3), sharex=True, sharey=False, gridspec_kw={'height_ratios': [1.5, 0.5]})
i = 0 #First main plot
ax[i].scatter(ss.dist_tss, -np.log10(ss.pval), color="black", edgecolors="#ffffff", linewidth=0.1)
ax[i].scatter(ss[ss.dist_tss==0].dist_tss, -np.log10(ss[ss.dist_tss==0].pval), color="tab:red", marker="D", label="HPR CNV", s=30)
ax[i].scatter(ss[ss.position==72072066].dist_tss, -np.log10(ss[ss.position==72072066].pval), color="tab:green", label="rs217184\n(intronic)")
ax[i].axhline(y=-np.log10(5*10**-8), linestyle="--", linewidth=1, color="tab:blue")
i = 1 #PIP
ax[i].scatter(ss.dist_tss, ss.min_pip, color="navy", edgecolors="#ffffff", linewidth=0.1)
ax[i].scatter(ss[ss.dist_tss==0].dist_tss, ss[ss.dist_tss==0].min_pip, color="tab:red", marker="D", label="HPR CNV", s=30)
ax[i].scatter(ss[ss.position==72072066].dist_tss, ss[ss.position==72072066].min_pip, color="tab:green")
ax[0].set_title("trans-pQTL effect of HPR cis-variant on APOL1 (Olink) in UKB AFR", fontsize=12)
ax[0].set_ylabel("-log10(p)")
ax[1].set_ylabel("PIP")
ax[1].set_xlabel("Distance to TSS of HPR")
for cnt, a in enumerate(ax.flat):
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)    
ax[0].legend(bbox_to_anchor=(1.02,1))
ax[0].set_ylim([-0.2,10])
ax[1].set_ylim([-0.05, 1.15])
ax[1].set_xlim([-2.05*10**5, 2.2*10**5])
plt.tight_layout()
fig.subplots_adjust(hspace=0.2)
plt.savefig("hpr_afr_manhattan_main.png", dpi=500)
plt.savefig("hpr_afr_manhattan_main.pdf", dpi=500)
plt.show()
plt.close()