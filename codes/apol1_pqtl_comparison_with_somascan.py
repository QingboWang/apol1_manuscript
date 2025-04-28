import pandas as pd
import numpy as np
import glob
import time as tm
from matplotlib import pyplot as plt
from scipy import stats
plt.rcParams.update({'font.size': 16})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"



#Fig. 1a: miami plot EUR
#Fig. 1 sub: effect size scatter plot for each region
#Fig. 1b: miami plot AFR
#Fig. 1c: effect size scatter plot for each region

#To make the position to abs:
chrlen = {
1:	248956422,2:	242193529,3:	198295559,4:	190214555,5:	181538259,6:	170805979,7:	159345973,8:	145138636,
9:	138394717,10:	133797422,11:	135086622,12:	133275309,13:	114364328,14:	107043718,15:	101991189,16:	90338345,
17:	83257441,18:	80373285,19:	58617616,20:	64444167,21:	46709983,22:	50818468,23:	156040895,24:	57227415}
offsets = {}
for i in range(25):
    offsets["chr{0}".format(i)] = sum([chrlen[j] for j in range(1,i)])

#Our Olink EUR:
df = pd.read_csv("/Users/qsw/Desktop/ukb_data/apol1_qtl_sumstats_eur_p005.tsv.gz", sep="\t", compression="gzip", index_col=0) #Filtered to nominal pass
df["minus_log10_pval"] = -np.log10(df.pval)
df.loc[np.isinf(df.minus_log10_pval),"minus_log10_pval"] = 320 #Cap at 320
df["chr"] = df.index.str.split("_").str[0]
df["pos"] = df.index.str.split("_").str[1].astype(int)
df["pos_abs"] = df.pos + df.chr.map(offsets)
df["minus_log10_pval"] = -np.log10(df.pval)

#SomaScan EUR
soma = pd.read_csv("~/Desktop/resources/decode_pqtl/11510_31_APOL1_Apo_L1.txt.gz", sep="\t")
soma = soma[soma.Pval<0.05] #Filtering to nominal pass
soma.loc[np.isinf(soma.minus_log10_pval),"minus_log10_pval"] = 320 #Cap at 320
soma = soma[(soma.Chrom!="chrX")&(soma.Chrom!="chrY")]
soma["pos_abs"] = soma.Pos + soma.Chrom.str.replace("X","23").str.replace("Y","24").map(offsets)

#Make the position absolute
xpos = []
for i in range(1,23):
    xpos.append((offsets["chr{0}".format(i)] + chrlen[i]/2))
xmax = offsets["chr22"] + chrlen[22]
fig, ax = plt.subplots(2, 1, figsize=(16, 7), sharex=True, sharey=False)
c0 = "#222222"
for i in list(range(1,23)):
    x = df[df.chr=="chr{0}".format(i)].pos_abs
    y = df[df.chr=="chr{0}".format(i)].minus_log10_pval
    ax[0].scatter(x, y, color=c0)
    x = soma[soma.Chrom=="chr{0}".format(i)].pos_abs
    y = soma[soma.Chrom=="chr{0}".format(i)].minus_log10_pval
    ax[1].scatter(x, y, color=c0)
    if c0=="#222222":
        c0 = "tab:grey"
    else: 
        c0 = "#222222"
#KNG1 in a specific color:
cond = (df.chr=="chr3")&(186717348-10**6<df.pos)&(df.pos<186744410+10**6)
x = df[cond].pos_abs
y = df[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:purple", label="KNG1 (chr3)")
cond = (soma.Chrom=="chr3")&(186717348-10**6<soma.Pos)&(soma.Pos<186744410+10**6)
x = soma[cond].pos_abs
y = soma[cond].minus_log10_pval
ax[1].scatter(x, y, color="tab:purple")
#KLKB1 in another specific color:
cond = (df.chr=="chr4")&(186208979-10**6<df.pos)&(df.pos<186258471+10**6)
x = df[cond].pos_abs
y = df[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:blue", label="KLKB1 (chr4)")
cond = (soma.Chrom=="chr4")&(186208979-10**6<soma.Pos)&(soma.Pos<186258471+10**6)
x = soma[cond].pos_abs
y = soma[cond].minus_log10_pval
ax[1].scatter(x, y, color="tab:blue")
#F12 in another specific color:
cond = (df.chr=="chr5")&(177402141-10**6<df.pos)&(df.pos<177409564+10**6)
x = df[cond].pos_abs
y = df[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:olive", label="F12 (chr5)")
cond = (soma.Chrom=="chr5")&(177402141-10**6<soma.Pos)&(soma.Pos<177409564+10**6)
x = soma[cond].pos_abs
y = soma[cond].minus_log10_pval
ax[1].scatter(x, y, color="tab:olive")
#HPR in another specific color:
cond = (df.chr=="chr16")&(72063148-10**6<df.pos)&(df.pos<72077246+10**6)
x = df[cond].pos_abs
y = df[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:green", label="HPR (chr16)")
cond = (soma.Chrom=="chr16")&(72063148-10**6<soma.Pos)&(soma.Pos<72077246+10**6)
x = soma[cond].pos_abs
y = soma[cond].minus_log10_pval
ax[1].scatter(x, y, color="tab:green")
#APOL1 in another specific color:
cond = (df.chr=="chr22")&(36253071-10**6<df.pos)&(df.pos<36267530+10**6)
x = df[cond].pos_abs
y = df[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:orange", label="APOL1 (chr22)")
cond = (soma.Chrom=="chr22")&(36253071-10**6<soma.Pos)&(soma.Pos<36267530+10**6)
x = soma[cond].pos_abs
y = soma[cond].minus_log10_pval

ax[1].scatter(x, y, color="tab:orange")
ax[1].invert_yaxis()  # Invert the y-axis
ax[0].legend(title="Gene (chromosome):", loc="upper left", fontsize=16)
ax[0].set_ylabel("-log10(p)")
ax[1].set_ylabel("-log10(p)")
ax[1].set_xlabel("Position (Chr)")
ax[1].set_xticks(xpos)
ax[1].set_xticklabels(list(range(1,23)))
ax[1].set_xlim([0 - 3*10**7, xmax + 3*10**7])
ax[0].set_ylim(bottom=0)
ax[0].set_title("pQTL effect on APOL1 in two studies with different technologies", fontsize=16)
ax[0].text(offsets["chr7"], 300, "Study: UKB EUR\nTechnology: Olink", ha="left", va="top", fontsize=16)
ax[1].text(offsets["chr7"], 50, "Study: deCODE\nTechnology: SomaScan", ha="left", va="top", fontsize=16)
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_miami_eur.pdf", dpi=250)
plt.close()


#Eff. size comparison scatter plto:
#Note: beta is the beta of "effectAllele" in soma
df = []
for chr in range(1,23): #Out olink data
    fns = glob.glob("/Users/qsw/Downloads/apol1_trans_eur/apol1_trans_eur_ss_chr{0}*.tsv.gz".format(chr))
    print ("starting chr{0}, {1} files, {2}".format(chr, len(fns), tm.ctime()))
    for cnt, fn in enumerate(fns):
        dfsub = pd.read_csv(fn, sep="\t", index_col=0)
        df.append(dfsub)
        if cnt%1000==0:
            print ("done {0} files {1}".format(cnt, tm.ctime()))
    print ("done chr{0}, {1}".format(chr, tm.ctime()))                   
df = pd.concat(df)#Keeping only p<0.05 for simplicity
df = df[~df.index.str.split("_").str[1].str.contains(r'[^0-9]')] #Removing non-canonical chrs
#SomaScan EUR
soma = pd.read_csv("~/Desktop/resources/decode_pqtl/11510_31_APOL1_Apo_L1.txt.gz", sep="\t")
soma.loc[np.isinf(soma.minus_log10_pval),"minus_log10_pval"] = 320 #Cap at 320
soma = soma[(soma.Chrom!="chrX")&(soma.Chrom!="chrY")]

#Parse for joining later:
soma.index = soma.Chrom+"_"+soma.Pos.astype(str)+"_"+soma.otherAllele+"_"+soma.effectAllele
df.index = df.index.str.split(";").str[0]
df["chr"] = df.index.str.split("_").str[0]
df["pos"] = df.index.str.split("_").str[1].astype(int)
soma_flipped = soma.copy(deep=True)
soma_flipped.index = soma_flipped.index = soma_flipped.Chrom+"_"+soma_flipped.Pos.astype(str)+"_"+soma_flipped.effectAllele+"_"+soma_flipped.otherAllele
soma_flipped["Beta"] = soma_flipped["Beta"]*-1
soma_flipped = soma_flipped[soma_flipped.effectAllele!=soma_flipped.otherAllele] #To avoid duplication
soma_full = pd.concat([soma, soma_flipped])

#APOL1 cis region:
cond = (df.chr=="chr22")&(36253071-10**6<df.pos)&(df.pos<36267530+10**6)
apo_olink = df[cond]
cond = (soma_full.Chrom=="chr22")&(36253071-10**6<soma_full.Pos)&(soma_full.Pos<36267530+10**6)
apo_soma = soma_full[cond]
#Join
apo = apo_olink[["beta","se", 'pval']].join(apo_soma[["Beta","SE", 'Pval']], how="inner", rsuffix="_soma")
apo = apo[(apo.pval < 0.05)&(apo.Pval < 0.05)] #Otherwise too noisy
#Pearson R in this filtered data:
r = stats.pearsonr(apo.beta, apo.Beta)[0]
#Plot:
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=apo.beta, y=apo.Beta, xerr=apo.se, yerr=apo.SE, alpha=0.4, fmt='o', color="tab:orange")
plt.title("Effect of APOL1 cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB EUR (Olink)")
plt.ylabel("Beta in deCODE (SomaScan)")
plt.text(x=0.1, y=0.3, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.4,0.4])
plt.ylim([-0.4,0.4])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_cis_effsize_eur.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_cis_effsize_eur.pdf", dpi=400)
plt.close()

#KNG1:
cond = (df.chr=="chr3")&(186717348-10**6<df.pos)&(df.pos<186744410+10**6)
kng1_olink = df[cond]
cond = (soma_full.Chrom=="chr3")&(186717348-10**6<soma_full.Pos)&(soma_full.Pos<186744410+10**6)
kng1_soma = soma_full[cond]
kng1 = kng1_olink[["beta","se", 'pval']].join(kng1_soma[["Beta","SE", 'Pval']], how="inner", rsuffix="_soma")
kng1 = kng1[(kng1.pval < 0.05)|(kng1.Pval < 0.05)] #Otherwise too noisy
#Pearson R in this filtered data:
r = stats.pearsonr(kng1.beta, kng1.Beta)[0]
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=kng1.beta, y=kng1.Beta, xerr=kng1.se, yerr=kng1.SE, alpha=0.4, fmt='o', color="tab:purple")
plt.title("Effect of KNG1 cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB EUR (Olink)")
plt.ylabel("Beta in deCODE (SomaScan)")
plt.text(x=-0.1, y=0.09, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.11,0.11])
plt.ylim([-0.11,0.11])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_kng1_effsize_eur.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_kng1_effsize_eur.pdf", dpi=400)
plt.close()

#HPR cis region:
cond = (df.chr=="chr16")&(72063148-10**6<df.pos)&(df.pos<72077246+10**6)
hpr_olink = df[cond]
cond = (soma_full.Chrom=="chr16")&(72063148-10**6<soma_full.Pos)&(soma_full.Pos<72077246+10**6)
hpr_soma = soma_full[cond]
hpr = hpr_olink[["beta","se", 'pval']].join(hpr_soma[["Beta","SE", 'Pval']], how="inner", rsuffix="_soma")
hpr = hpr[(hpr.pval < 0.05)&(hpr.Pval < 0.05)] #Otherwise too noisy
#Pearson R in this filtered data:
r = stats.pearsonr(hpr.beta, hpr.Beta)[0]
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=hpr.beta, y=hpr.Beta, xerr=hpr.se, yerr=hpr.SE, alpha=0.4, fmt='o', color="tab:green")
plt.title("Effect of HPR cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB EUR (Olink)")
plt.ylabel("Beta in deCODE (SomaScan)")
plt.text(x=-0.5, y=0.45, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.6,0.6])
plt.ylim([-0.6,0.6])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_hpr_effsize_eur.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_hpr_effsize_eur.pdf", dpi=400)
plt.close()

#F12 cis region:
cond = (df.chr=="chr5")&(177402141-10**6<df.pos)&(df.pos<177409564+10**6)
f12_olink = df[cond]
cond = (soma_full.Chrom=="chr5")&(177402141-10**6<soma_full.Pos)&(soma_full.Pos<177409564+10**6)
f12_soma = soma_full[cond]
f12 = f12_olink[["beta","se","pval"]].join(f12_soma[["Beta","SE","Pval"]], how="inner")
f12 = f12[(f12.pval < 0.05)|(f12.Pval < 0.05)] #Otherwise too noisy, but "and" is too few
#Pearson R in this filtered data:
r = stats.pearsonr(f12.beta, f12.Beta)[0]
#Remove a few outlier for visual purpose:
f12 = f12[(abs(f12.beta)<0.2)&(abs(f12.Beta)<0.2)]
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=f12.beta, y=f12.Beta, xerr=f12.se, yerr=f12.SE, alpha=0.4, fmt='o', color="tab:olive")
plt.title("Effect of F12 cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB EUR (Olink)")
plt.ylabel("Beta in deCODE (SomaScan)")
plt.text(x=-0.18, y=0.15, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.2,0.2])
plt.ylim([-0.2,0.2])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_f12_effsize_eur.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_f12_effsize_eur.pdf", dpi=400)
plt.close()

#KLKB1 cis region:
cond = (df.chr=="chr4")&(186208979-10**6<df.pos)&(df.pos<186258471+10**6)
klkb1_olink = df[cond]
cond = (soma_full.Chrom=="chr4")&(186208979-10**6<soma_full.Pos)&(soma_full.Pos<186258471+10**6)
klkb1_soma = soma_full[cond]
klkb1 = klkb1_olink[["beta","se", "pval"]].join(klkb1_soma[["Beta","SE", "Pval"]], how="inner")
klkb1 = klkb1[(klkb1.pval < 0.05)|(klkb1.Pval < 0.05)] #Otherwise too noisy, but "and" is too few
r = stats.pearsonr(f12.beta, f12.Beta)[0]
#Remove a few outlier for visual purpose:
klkb1 = klkb1[(abs(klkb1.beta)<0.2)&(abs(klkb1.Beta)<0.2)]
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=klkb1.beta, y=klkb1.Beta, xerr=klkb1.se, yerr=klkb1.SE, alpha=0.4, fmt='o', color="tab:blue")
plt.title("Effect of KLKB1 cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB EUR (Olink)")
plt.ylabel("Beta in DECODE (SomaScan)")
plt.text(x=-0.18, y=0.15, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.2,0.2])
plt.ylim([-0.2,0.2])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_klkb1_effsize_eur.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_klkb1_effsize_eur.pdf", dpi=400)
plt.close()


#Doing the same for AFR:
df = []
for chr in range(1,23):
    dfsub = pd.read_csv("~/Desktop/ukb_data/apol1_trans_afr/apol1_trans_afr_ss_chr{0}.tsv.gz".format(chr), sep="\t", index_col=0)
    df.append(dfsub[dfsub.pval<0.05])
df = pd.concat(df)#Keeping only p<0.05 for simplicity
df = df[~df.index.str.split("_").str[1].str.contains(r'[^0-9]')] #Removing non-canonical chrs
#soma = pd.read_csv("~/Desktop/resources/somascan_afr_2022/GCST90233273_buildGRCh38.tsv.gz", sep="\t") #corresponds to 31
#soma = pd.read_csv("~/Desktop/resources/somascan_afr_2022/GCST90239378_buildGRCh38.tsv.gz", sep="\t") #corresponds to the noise guy
soma = pd.read_csv("~/Desktop/resources/somascan_afr_2022/GCST90233274_buildGRCh38.tsv.gz", sep="\t") #corresponds to 51
#Filter to pval<0.05:
soma = soma[soma.p_value<0.05]

#Parse:
df["minus_log10_pval"] = -np.log10(df.pval)
df.loc[np.isinf(df.minus_log10_pval),"minus_log10_pval"] = 320 #Cap at 320
df["chr"] = df.index.str.split("_").str[0]
df["pos"] = df.index.str.split("_").str[1].astype(int)
df["pos_abs"] = df.pos + df.chr.map(offsets)
soma["minus_log10_pval"] = -np.log10(soma.p_value)
soma.loc[np.isinf(soma.minus_log10_pval),"minus_log10_pval"] = 320 #Cap at 320
soma = soma[(soma.chromosome!="X")&(soma.chromosome!="Y")]
soma["pos_abs"] = soma.base_pair_location + ("chr"+soma.chromosome.astype(str)).map(offsets)

#Plot:

#First, manhattan:
#Do some pruning first:
dfsub = df.copy(deep=True)
dfsub.sort_values(by="pos_abs", inplace=True)
dfsub["to_skip"] = False
shape_before = dfsub.shape[0]
for iter in range(10):
    for i in range(dfsub.shape[0]-1):
        if (dfsub.pval.iloc[i]>0.05*10**-3) & (dfsub.pos.iloc[i+1]-dfsub.pos.iloc[i]<10**4) & (dfsub.chr.iloc[i+1]==dfsub.chr.iloc[i]):
            dfsub.to_skip.iloc[i] = True #Two-minute ish in the first iter.
    dfsub = dfsub[~dfsub.to_skip]
    if dfsub.shape[0]==shape_before:
        break
    shape_before = dfsub.shape[0]
somasub = soma.copy(deep=True) #For now
somasub.sort_values(by="pos_abs", inplace=True)
somasub["to_skip"] = False
shape_before = somasub.shape[0]
for iter in range(10):
    for i in range(somasub.shape[0]-1):
        if (somasub.p_value.iloc[i]>0.05*10**-3) & (somasub.base_pair_location.iloc[i+1]-somasub.base_pair_location.iloc[i]<10**4) & (somasub.chromosome.iloc[i+1]==somasub.chromosome.iloc[i]):
            somasub.to_skip.iloc[i] = True #Two-minute ish
    somasub = somasub[~somasub.to_skip]        
    if somasub.shape[0]==shape_before:
        break
    shape_before = somasub.shape[0]
plt.rcParams.update({'font.size': 16})
xpos = []
for i in range(1,23):
    xpos.append((offsets["chr{0}".format(i)] + chrlen[i]/2))
xmax = offsets["chr22"] + chrlen[22]
fig, ax = plt.subplots(2, 1, figsize=(16, 7), sharex=True, sharey=False, gridspec_kw={'height_ratios': [2.2, 1]})
c0 = "#222222"
for i in list(range(1,23)):
    x = dfsub[dfsub.chr=="chr{0}".format(i)].pos_abs
    y = dfsub[dfsub.chr=="chr{0}".format(i)].minus_log10_pval
    ax[0].scatter(x, y, color=c0)
    x = somasub[somasub.chromosome==i].pos_abs
    y = somasub[somasub.chromosome==i].minus_log10_pval
    ax[1].scatter(x, y, color=c0)
    if c0=="#222222":
        c0 = "tab:grey"
    else: 
        c0 = "#222222"
#KNG1:
cond = (dfsub.chr=="chr3")&(186717348-10**6<dfsub.pos)&(dfsub.pos<186744410+10**6)
x = dfsub[cond].pos_abs
y = dfsub[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:purple", label="KNG1 (chr3)")
cond = (somasubsub.chromosome==3)&(186717348-10**6<somasubsub.base_pair_location)&(somasubsub.base_pair_location<186744410+10**6)
x = somasubsub[cond].pos_abs
y = somasubsub[cond].minus_log10_pval
ax[1].scatter(x, y, color="tab:purple")
#KLKB1
cond = (dfsub.chr=="chr4")&(186208979-10**6<dfsub.pos)&(dfsub.pos<186258471+10**6)
x = dfsub[cond].pos_abs
y = dfsub[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:blue", label="KLKB1 (chr4)")
cond = (somasubsub.chromosome==4)&(186208979-10**6<somasubsub.base_pair_location)&(somasubsub.base_pair_location<186258471+10**6)
x = somasubsub[cond].pos_abs
y = somasubsub[cond].minus_log10_pval
ax[1].scatter(x, y, color="tab:blue")
#F12
cond = (dfsub.chr=="chr5")&(177402141-10**6<dfsub.pos)&(dfsub.pos<177409564+10**6)
x = dfsub[cond].pos_abs
y = dfsub[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:olive", label="F12 (chr5)")
cond = (somasubsub.chromosome==5)&(177402141-10**6<somasubsub.base_pair_location)&(somasubsub.base_pair_location<177409564+10**6)
x = somasubsub[cond].pos_abs
y = somasubsub[cond].minus_log10_pval
ax[1].scatter(x, y, color="tab:olive")
#HPR
cond = (dfsub.chr=="chr16")&(72063148-10**6<dfsub.pos)&(dfsub.pos<72077246+10**6)
x = dfsub[cond].pos_abs
y = dfsub[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:green", label="HPR (chr16)")
cond = (somasubsub.chromosome==16)&(72063148-10**6<somasubsub.base_pair_location)&(somasubsub.base_pair_location<72077246+10**6)
x = somasubsub[cond].pos_abs
y = somasubsub[cond].minus_log10_pval
ax[1].scatter(x, y, color="tab:green")
#APOL1:
cond = (dfsub.chr=="chr22")&(36253071-10**6<dfsub.pos)&(dfsub.pos<36267530+10**6)
x = dfsub[cond].pos_abs
y = dfsub[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:orange", label="APOL1 (chr22)")
cond = (somasubsub.chromosome==22)&(36253071-10**6<somasubsub.base_pair_location)&(somasubsub.base_pair_location<36267530+10**6)
x = somasubsub[cond].pos_abs
y = somasubsub[cond].minus_log10_pval
ax[1].scatter(x, y, color="tab:orange")
ax[1].invert_yaxis()  # Invert the y-axis
ax[0].legend(title="Gene (chromosomeosome):", loc="upper left", fontsize=14)
ax[0].set_ylabel("-log10(p)")
ax[1].set_ylabel("-log10(p)")
ax[1].set_xlabel("Position (Chr)")
ax[1].set_xticks(xpos)
ax[1].set_xticklabels(list(range(1,23)))
ax[1].set_xlim([0 - 3*10**7, xmax + 3*10**7])
ax[0].set_ylim(bottom=0)
ax[1].set_ylim(top=0)
ax[0].set_title("pQTL effect on APOL1 in two studies with different technologies", fontsize=16)
ax[0].text(offsets["chr7"], 80, "Study: UKB AFR\nTechnology: Olink", ha="left", va="top", fontsize=18)
ax[1].text(offsets["chr7"], 20, "Study: AASK AFR\nTechnology: Somascan\n(Aptamer ID: 11510_51)", ha="left", va="top", fontsize=16)
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_miami_afr_ID_11510_51.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_miami_afr_ID_11510_51.pdf", dpi=400)
plt.close()


#Effect size (not filtering to p005 for this a priori):
plt.rcParams.update({'font.size': 12})
df = []
for chr in range(1,23):
    dfsub = pd.read_csv("~/Desktop/ukb_data/apol1_trans_afr/apol1_trans_afr_ss_chr{0}.tsv.gz".format(chr), sep="\t", index_col=0)
    df.append(dfsub[dfsub.pval<0.05])
df = pd.concat(df)#Keeping only p<0.05 for simplicity
df = df[~df.index.str.split("_").str[1].str.contains(r'[^0-9]')] #Removing non-canonical chrs
#soma = pd.read_csv("~/Desktop/resources/somascan_afr_2022/GCST90233273_buildGRCh38.tsv.gz", sep="\t") #corresponds to 31
#soma = pd.read_csv("~/Desktop/resources/somascan_afr_2022/GCST90239378_buildGRCh38.tsv.gz", sep="\t") #corresponds to the noise guy
soma = pd.read_csv("~/Desktop/resources/somascan_afr_2022/GCST90233274_buildGRCh38.tsv.gz", sep="\t") #corresponds to 51

#Parse:
df["minus_log10_pval"] = -np.log10(df.pval)
df.loc[np.isinf(df.minus_log10_pval),"minus_log10_pval"] = 320 #Cap at 320
df["chr"] = df.index.str.split("_").str[0]
df["pos"] = df.index.str.split("_").str[1].astype(int)
df["pos_abs"] = df.pos + df.chr.map(offsets)
soma["minus_log10_pval"] = -np.log10(soma.p_value)
soma.loc[np.isinf(soma.minus_log10_pval),"minus_log10_pval"] = 320 #Cap at 320
soma = soma[(soma.chromosome!="X")&(soma.chromosome!="Y")]
soma["pos_abs"] = soma.base_pair_location + ("chr"+soma.chromosome.astype(str)).map(offsets)

#Prep. for join:
soma.index = "chr"+soma.chromosome.astype(str)+"_"+soma.base_pair_location.astype(str)+"_"+soma.other_allele+"_"+soma.effect_allele
df.index = df.index.str.split(";").str[0]
df["chr"] = df.index.str.split("_").str[0]
df["pos"] = df.index.str.split("_").str[1].astype(int)
soma_flipped = soma.copy(deep=True)
soma_flipped.index = "chr"+soma_flipped.chromosome.astype(str)+"_"+soma_flipped.base_pair_location.astype(str)+"_"+soma_flipped.effect_allele+"_"+soma_flipped.other_allele
soma_flipped["beta"] = soma_flipped["beta"]*-1
soma_flipped = soma_flipped[soma_flipped.effect_allele!=soma_flipped.other_allele] #To avoid duplication
soma_full = pd.concat([soma, soma_flipped])

#Plot:

#KNG1 cis region:
cond = (df.chr=="chr3")&(186717348-10**6<df.pos)&(df.pos<186744410+10**6)
apo_olink = df[cond]
cond = (soma_full.chromosome==3)&(186717348-10**6<soma_full.base_pair_location)&(soma_full.base_pair_location<186744410+10**6)
apo_soma = soma_full[cond]
#Join
apo = apo_olink[["beta","se", 'pval']].join(apo_soma[["beta", "standard_error", "p_value"]], how="inner", rsuffix="_soma")
apo = apo[(apo.pval < 0.05)|(apo.p_value < 0.05)] #Otherwise too noisy
#Pearson R in this filtered data:
r = stats.pearsonr(apo.beta, apo.beta_soma)[0]
#Plot:
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=apo.beta, y=apo.beta_soma, xerr=apo.se, yerr=apo.standard_error, alpha=0.4, fmt='o', color="tab:purple")
plt.title("Effect of KNG1 cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB AFR (Olink)")
plt.ylabel("Beta in AASK (SomaScan)")
plt.text(x=-0.4, y=0.42, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.5,.5])
plt.ylim([-0.5,.5])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol1_kng1_effsize_afr_ID_11510_51.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol1_kng1_effsize_afr_ID_11510_51.pdf", dpi=400)
plt.close()

#APOL1 cis region:
cond = (df.chr=="chr22")&(36253071-10**6<df.pos)&(df.pos<36267530+10**6)
apo_olink = df[cond]
cond = (soma_full.chromosome==22)&(36253071-10**6<soma_full.base_pair_location)&(soma_full.base_pair_location<36267530+10**6)
apo_soma = soma_full[cond]
#Join
apo = apo_olink[["beta","se", 'pval']].join(apo_soma[["beta", "standard_error", "p_value"]], how="inner", rsuffix="_soma")
apo = apo[(apo.pval < 0.05)&(apo.p_value < 0.05)] #Otherwise too noisy
#Pearson R in this filtered data:
r = stats.pearsonr(apo.beta, apo.beta_soma)[0]
#Plot:
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=apo.beta, y=apo.beta_soma, xerr=apo.se, yerr=apo.standard_error, alpha=0.4, fmt='o', color="tab:orange")
plt.title("Effect of APOL1 cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB AFR (Olink)")
plt.ylabel("Beta in AASK (SomaScan)")
plt.text(x=0.1, y=1.1, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.8,1.3])
plt.ylim([-0.8,1.3])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_cis_effsize_afr_ID_11510_51.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_cis_effsize_afr_ID_11510_51.pdf", dpi=400)
plt.close()

#HPR cis region:
cond = (df.chr=="chr16")&(72063148-10**6<df.pos)&(df.pos<72077246+10**6)
hpr_olink = df[cond]
cond = (soma_full.chromosome==16)&(72063148-10**6<soma_full.base_pair_location)&(soma_full.base_pair_location<72077246+10**6)
hpr_soma = soma_full[cond]
hpr = hpr_olink[["beta","se", 'pval']].join(hpr_soma[["beta", "standard_error", "p_value"]], how="inner", rsuffix="_soma")
hpr = hpr[(hpr.pval < 0.05)&(hpr.p_value < 0.05)] #Otherwise too noisy
#Pearson R in this filtered data:
r = stats.pearsonr(hpr.beta, hpr.beta_soma)[0]
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=hpr.beta, y=hpr.beta_soma, xerr=hpr.se, yerr=hpr.standard_error, alpha=0.4, fmt='o', color="tab:green")
plt.title("Effect of HPR cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB AFR (Olink)")
plt.ylabel("Beta in AASK (SomaScan)")
plt.text(x=-0.5, y=0.6, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.6,0.6])
plt.ylim([-0.7,0.7])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_hpr_effsize_afr_ID_11510_51.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_hpr_effsize_afr_ID_11510_51.pdf", dpi=400)
plt.close()

#F12 cis region:
cond = (df.chr=="chr5")&(177402141-10**6<df.pos)&(df.pos<177409564+10**6)
f12_olink = df[cond]
cond = (soma_full.chromosome==5)&(177402141-10**6<soma_full.base_pair_location)&(soma_full.base_pair_location<177409564+10**6)
f12_soma = soma_full[cond]
f12 = f12_olink[["beta","se","pval"]].join(f12_soma[["beta", "standard_error", "p_value"]], how="inner", rsuffix="_soma")
#Pearson R in this filtered data:
r = stats.pearsonr(f12.beta, f12.beta_soma)[0]
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=f12.beta, y=f12.beta_soma, xerr=f12.se, yerr=f12.standard_error, alpha=0.4, fmt='o', color="tab:olive")
plt.title("Effect of F12 cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB AFR (Olink)")
plt.ylabel("Beta in AASK (SomaScan)")
plt.text(x=-0.48, y=0.4, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.5,0.5])
plt.ylim([-0.5,0.5])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_f12_effsize_afr_ID_11510_51.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_f12_effsize_afr_ID_11510_51.pdf", dpi=400)
plt.close()


#KLKB1 cis region:
cond = (df.chr=="chr4")&(186208979-10**6<df.pos)&(df.pos<186258471+10**6)
klkb1_olink = df[cond]
cond = (soma_full.chromosome==4)&(186208979-10**6<soma_full.base_pair_location)&(soma_full.base_pair_location<186258471+10**6)
klkb1_soma = soma_full[cond]
klkb1 = klkb1_olink[["beta","se", "pval"]].join(klkb1_soma[["beta", "standard_error", "p_value"]], how="inner", rsuffix="_soma")
r = stats.pearsonr(klkb1.beta, klkb1.beta_soma)[0]
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=klkb1.beta, y=klkb1.beta_soma, xerr=klkb1.se, yerr=klkb1.standard_error, alpha=0.4, fmt='o', color="tab:blue")
plt.title("Effect of KLKB1 cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB AFR (Olink)")
plt.ylabel("Beta in AASK (SomaScan)")
plt.text(x=-0.53, y=0.45, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.55,0.55])
plt.ylim([-0.55,0.55])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_klkb1_effsize_afr_ID_11510_51.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_klkb1_effsize_afr_ID_11510_51.pdf", dpi=400)
plt.close()


#The same thing with the other soma data in AFR (..31):
df = []
for chr in range(1,23):
    dfsub = pd.read_csv("~/Desktop/ukb_data/apol1_trans_afr/apol1_trans_afr_ss_chr{0}.tsv.gz".format(chr), sep="\t", index_col=0)
    df.append(dfsub[dfsub.pval<0.05])
df = pd.concat(df)#Keeping only p<0.05 for simplicity
df = df[~df.index.str.split("_").str[1].str.contains(r'[^0-9]')] #Removing non-canonical chrs
soma = pd.read_csv("~/Desktop/resources/somascan_afr_2022/GCST90233273_buildGRCh38.tsv.gz", sep="\t") #corresponds to 31
#soma = pd.read_csv("~/Desktop/resources/somascan_afr_2022/GCST90239378_buildGRCh38.tsv.gz", sep="\t") #corresponds to the noise guy
#soma = pd.read_csv("~/Desktop/resources/somascan_afr_2022/GCST90233274_buildGRCh38.tsv.gz", sep="\t") #corresponds to 51
#Filter to pval<0.05:
soma = soma[soma.p_value<0.05]

#Parse:
df["minus_log10_pval"] = -np.log10(df.pval)
df.loc[np.isinf(df.minus_log10_pval),"minus_log10_pval"] = 320 #Cap at 320
df["chr"] = df.index.str.split("_").str[0]
df["pos"] = df.index.str.split("_").str[1].astype(int)
df["pos_abs"] = df.pos + df.chr.map(offsets)
soma["minus_log10_pval"] = -np.log10(soma.p_value)
soma.loc[np.isinf(soma.minus_log10_pval),"minus_log10_pval"] = 320 #Cap at 320
soma = soma[(soma.chromosome!="X")&(soma.chromosome!="Y")]
soma["pos_abs"] = soma.base_pair_location + ("chr"+soma.chromosome.astype(str)).map(offsets)

#Make it sparse:
dfsub = df.copy(deep=True)
dfsub.sort_values(by="pos_abs", inplace=True)
dfsub["to_skip"] = False
shape_before = dfsub.shape[0]
for iter in range(10):
    for i in range(dfsub.shape[0]-1):
        if (dfsub.pval.iloc[i]>0.05*10**-3) & (dfsub.pos.iloc[i+1]-dfsub.pos.iloc[i]<10**4) & (dfsub.chr.iloc[i+1]==dfsub.chr.iloc[i]):
            dfsub.to_skip.iloc[i] = True #Two-minute ish in the first iter.
    dfsub = dfsub[~dfsub.to_skip]
    if dfsub.shape[0]==shape_before:
        break
    shape_before = dfsub.shape[0]
somasub = soma.copy(deep=True) #For now
somasub.sort_values(by="pos_abs", inplace=True)
somasub["to_skip"] = False
shape_before = somasub.shape[0]
for iter in range(10):
    for i in range(somasub.shape[0]-1):
        if (somasub.p_value.iloc[i]>0.05*10**-1) & (somasub.base_pair_location.iloc[i+1]-somasub.base_pair_location.iloc[i]<10**4) & (somasub.chromosome.iloc[i+1]==somasub.chromosome.iloc[i]):
            somasub.to_skip.iloc[i] = True #Two-minute ish
    somasub = somasub[~somasub.to_skip]        
    if somasub.shape[0]==shape_before:
        break
    shape_before = somasub.shape[0]
#Plot:
plt.rcParams.update({'font.size': 16})
xpos = []
for i in range(1,23):
    xpos.append((offsets["chr{0}".format(i)] + chrlen[i]/2))
xmax = offsets["chr22"] + chrlen[22]
fig, ax = plt.subplots(2, 1, figsize=(16, 7), sharex=True, sharey=False, gridspec_kw={'height_ratios': [2.2, 1]})
c0 = "#222222"
for i in list(range(1,23)):
    x = dfsub[dfsub.chr=="chr{0}".format(i)].pos_abs
    y = dfsub[dfsub.chr=="chr{0}".format(i)].minus_log10_pval
    ax[0].scatter(x, y, color=c0)
    x = somasub[somasub.chromosome==i].pos_abs
    y = somasub[somasub.chromosome==i].minus_log10_pval
    ax[1].scatter(x, y, color=c0)
    if c0=="#222222":
        c0 = "tab:grey"
    else: 
        c0 = "#222222"
#KNG1:
cond = (dfsub.chr=="chr3")&(186717348-10**6<dfsub.pos)&(dfsub.pos<186744410+10**6)
x = dfsub[cond].pos_abs
y = dfsub[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:purple", label="KNG1 (chr3)")
cond = (somasub.chromosome==3)&(186717348-10**6<somasub.base_pair_location)&(somasub.base_pair_location<186744410+10**6)
x = somasub[cond].pos_abs
y = somasub[cond].minus_log10_pval
ax[1].scatter(x, y, color="tab:purple")
#KLKB1
cond = (dfsub.chr=="chr4")&(186208979-10**6<dfsub.pos)&(dfsub.pos<186258471+10**6)
x = dfsub[cond].pos_abs
y = dfsub[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:blue", label="KLKB1 (chr4)")
cond = (somasub.chromosome==4)&(186208979-10**6<somasub.base_pair_location)&(somasub.base_pair_location<186258471+10**6)
x = somasub[cond].pos_abs
y = somasub[cond].minus_log10_pval
ax[1].scatter(x, y, color="tab:blue")
#F12
cond = (dfsub.chr=="chr5")&(177402141-10**6<dfsub.pos)&(dfsub.pos<177409564+10**6)
x = dfsub[cond].pos_abs
y = dfsub[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:olive", label="F12 (chr5)")
cond = (somasub.chromosome==5)&(177402141-10**6<somasub.base_pair_location)&(somasub.base_pair_location<177409564+10**6)
x = somasub[cond].pos_abs
y = somasub[cond].minus_log10_pval
ax[1].scatter(x, y, color="tab:olive")
#HPR
cond = (dfsub.chr=="chr16")&(72063148-10**6<dfsub.pos)&(dfsub.pos<72077246+10**6)
x = dfsub[cond].pos_abs
y = dfsub[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:green", label="HPR (chr16)")
cond = (somasub.chromosome==16)&(72063148-10**6<somasub.base_pair_location)&(somasub.base_pair_location<72077246+10**6)
x = somasub[cond].pos_abs
y = somasub[cond].minus_log10_pval
ax[1].scatter(x, y, color="tab:green")
#APOL1:
cond = (dfsub.chr=="chr22")&(36253071-10**6<dfsub.pos)&(dfsub.pos<36267530+10**6)
x = dfsub[cond].pos_abs
y = dfsub[cond].minus_log10_pval
ax[0].scatter(x, y, color="tab:orange", label="APOL1 (chr22)")
cond = (somasub.chromosome==22)&(36253071-10**6<somasub.base_pair_location)&(somasub.base_pair_location<36267530+10**6)
x = somasub[cond].pos_abs
y = somasub[cond].minus_log10_pval
ax[1].scatter(x, y, color="tab:orange")
ax[1].invert_yaxis()  # Invert the y-axis
ax[0].legend(title="Gene (chromosomeosome):", loc="upper left", fontsize=14)
ax[0].set_ylabel("-log10(p)")
ax[1].set_ylabel("-log10(p)")
ax[1].set_xlabel("Position (Chr)")
ax[1].set_xticks(xpos)
ax[1].set_xticklabels(list(range(1,23)))
ax[1].set_xlim([0 - 3*10**7, xmax + 3*10**7])
ax[0].set_ylim(bottom=0)
#ax[0].set_ylim([0, 40]) #test
ax[1].set_ylim(top=0)
ax[0].set_title("pQTL effect on APOL1 in two studies with different technologies", fontsize=16)
ax[0].text(offsets["chr7"], 80, "Study: UKB AFR\nTechnology: Olink", ha="left", va="top", fontsize=18)
ax[1].text(offsets["chr7"], 8, "Study: AASK AFR\nTechnology: Somascan (Aptamer ID 11510_31)", ha="left", va="top", fontsize=16)
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_miami_afr_ID_11510_31.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_miami_afr_ID_11510_31.pdf", dpi=400)
plt.close()


#Effect size (not filtering to p005 for this a priori):
df = []
for chr in range(1,23):
    dfsub = pd.read_csv("~/Desktop/ukb_data/apol1_trans_afr/apol1_trans_afr_ss_chr{0}.tsv.gz".format(chr), sep="\t", index_col=0)
    df.append(dfsub[dfsub.pval<0.05])
df = pd.concat(df)#Keeping only p<0.05 for simplicity
df = df[~df.index.str.split("_").str[1].str.contains(r'[^0-9]')] #Removing non-canonical chrs
soma = pd.read_csv("~/Desktop/resources/somascan_afr_2022/GCST90233273_buildGRCh38.tsv.gz", sep="\t") #corresponds to 31
#soma = pd.read_csv("~/Desktop/resources/somascan_afr_2022/GCST90239378_buildGRCh38.tsv.gz", sep="\t") #corresponds to the noise guy
#soma = pd.read_csv("~/Desktop/resources/somascan_afr_2022/GCST90233274_buildGRCh38.tsv.gz", sep="\t") #corresponds to 51
#Filter to pval<0.05:
#soma = soma[soma.p_value<0.05]

#Parse:
df["minus_log10_pval"] = -np.log10(df.pval)
df.loc[np.isinf(df.minus_log10_pval),"minus_log10_pval"] = 320 #Cap at 320
df["chr"] = df.index.str.split("_").str[0]
df["pos"] = df.index.str.split("_").str[1].astype(int)
df["pos_abs"] = df.pos + df.chr.map(offsets)
soma["minus_log10_pval"] = -np.log10(soma.p_value)
soma.loc[np.isinf(soma.minus_log10_pval),"minus_log10_pval"] = 320 #Cap at 320
soma = soma[(soma.chromosome!="X")&(soma.chromosome!="Y")]
soma["pos_abs"] = soma.base_pair_location + ("chr"+soma.chromosome.astype(str)).map(offsets)

#Prep. for join:
soma.index = "chr"+soma.chromosome.astype(str)+"_"+soma.base_pair_location.astype(str)+"_"+soma.other_allele+"_"+soma.effect_allele
df.index = df.index.str.split(";").str[0]
df["chr"] = df.index.str.split("_").str[0]
df["pos"] = df.index.str.split("_").str[1].astype(int)
soma_flipped = soma.copy(deep=True)
soma_flipped.index = "chr"+soma_flipped.chromosome.astype(str)+"_"+soma_flipped.base_pair_location.astype(str)+"_"+soma_flipped.effect_allele+"_"+soma_flipped.other_allele
soma_flipped["beta"] = soma_flipped["beta"]*-1
soma_flipped = soma_flipped[soma_flipped.effect_allele!=soma_flipped.other_allele] #To avoid duplication
soma_full = pd.concat([soma, soma_flipped])

#Plot:
plt.rcParams.update({'font.size': 12})
#KNG1 cis region:
cond = (df.chr=="chr3")&(186717348-10**6<df.pos)&(df.pos<186744410+10**6)
apo_olink = df[cond]
cond = (soma_full.chromosome==3)&(186717348-10**6<soma_full.base_pair_location)&(soma_full.base_pair_location<186744410+10**6)
apo_soma = soma_full[cond]
#Join
apo = apo_olink[["beta","se", 'pval']].join(apo_soma[["beta", "standard_error", "p_value"]], how="inner", rsuffix="_soma")
apo = apo[(apo.pval < 0.05)|(apo.p_value < 0.05)] #Otherwise too noisy
#Pearson R in this filtered data:
r = stats.pearsonr(apo.beta, apo.beta_soma)[0]
#Plot:
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=apo.beta, y=apo.beta_soma, xerr=apo.se, yerr=apo.standard_error, alpha=0.4, fmt='o', color="tab:purple")
plt.title("Effect of KNG1 cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB AFR (Olink)")
plt.ylabel("Beta in AASK\n(SomaScan aptamer 2)")
plt.text(x=-0.35, y=0.35, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.4,.4])
plt.ylim([-0.4,.4])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol1_kng1_effsize_afr_ID_11510_31.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol1_kng1_effsize_afr_ID_11510_31.pdf", dpi=400)
plt.close()


#APOL1 cis region:
cond = (df.chr=="chr22")&(36253071-10**6<df.pos)&(df.pos<36267530+10**6)
apo_olink = df[cond]
cond = (soma_full.chromosome==22)&(36253071-10**6<soma_full.base_pair_location)&(soma_full.base_pair_location<36267530+10**6)
apo_soma = soma_full[cond]
#Join
apo = apo_olink[["beta","se", 'pval']].join(apo_soma[["beta", "standard_error", "p_value"]], how="inner", rsuffix="_soma")
apo = apo[(apo.pval < 0.05)&(apo.p_value < 0.05)] #Otherwise too noisy
#Pearson R in this filtered data:
r = stats.pearsonr(apo.beta, apo.beta_soma)[0]
#Plot:
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=apo.beta, y=apo.beta_soma, xerr=apo.se, yerr=apo.standard_error, alpha=0.4, fmt='o', color="tab:orange")
plt.title("Effect of APOL1 cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB AFR (Olink)")
plt.ylabel("Beta in AASK\n(SomaScan aptamer 2)")
plt.text(x=0.1, y=0.7, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.8,1.3])
plt.ylim([-0.8,0.8])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_cis_effsize_afr_ID_11510_31.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_cis_effsize_afr_ID_11510_31.pdf", dpi=400)
plt.close()


#HPR cis region:
cond = (df.chr=="chr16")&(72063148-10**6<df.pos)&(df.pos<72077246+10**6)
hpr_olink = df[cond]
cond = (soma_full.chromosome==16)&(72063148-10**6<soma_full.base_pair_location)&(soma_full.base_pair_location<72077246+10**6)
hpr_soma = soma_full[cond]
hpr = hpr_olink[["beta","se", 'pval']].join(hpr_soma[["beta", "standard_error", "p_value"]], how="inner", rsuffix="_soma")
hpr = hpr[(hpr.pval < 0.05)&(hpr.p_value < 0.05)] #Otherwise too noisy
#Pearson R in this filtered data:
r = stats.pearsonr(hpr.beta, hpr.beta_soma)[0]
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=hpr.beta, y=hpr.beta_soma, xerr=hpr.se, yerr=hpr.standard_error, alpha=0.4, fmt='o', color="tab:green")
plt.title("Effect of HPR cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB AFR (Olink)")
plt.ylabel("Beta in AASK\n(SomaScan aptamer 2)")
plt.text(x=-0.5, y=0.6, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.6,0.6])
plt.ylim([-0.75,0.75])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_hpr_effsize_afr_ID_11510_31.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_hpr_effsize_afr_ID_11510_31.pdf", dpi=400)
plt.close()


#F12 cis region:
cond = (df.chr=="chr5")&(177402141-10**6<df.pos)&(df.pos<177409564+10**6)
f12_olink = df[cond]
cond = (soma_full.chromosome==5)&(177402141-10**6<soma_full.base_pair_location)&(soma_full.base_pair_location<177409564+10**6)
f12_soma = soma_full[cond]
f12 = f12_olink[["beta","se","pval"]].join(f12_soma[["beta", "standard_error", "p_value"]], how="inner", rsuffix="_soma")
#Pearson R in this filtered data:
r = stats.pearsonr(f12.beta, f12.beta_soma)[0]
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=f12.beta, y=f12.beta_soma, xerr=f12.se, yerr=f12.standard_error, alpha=0.4, fmt='o', color="tab:olive")
plt.title("Effect of F12 cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB AFR (Olink)")
plt.ylabel("Beta in AASK\n(SomaScan aptamer2)")
plt.text(x=-0.48, y=0.4, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.5,0.5])
plt.ylim([-0.5,0.5])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_f12_effsize_afr_ID_11510_31.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_f12_effsize_afr_ID_11510_31.pdf", dpi=400)
plt.close()


#KLKB1 cis region:
cond = (df.chr=="chr4")&(186208979-10**6<df.pos)&(df.pos<186258471+10**6)
klkb1_olink = df[cond]
cond = (soma_full.chromosome==4)&(186208979-10**6<soma_full.base_pair_location)&(soma_full.base_pair_location<186258471+10**6)
klkb1_soma = soma_full[cond]
klkb1 = klkb1_olink[["beta","se", "pval"]].join(klkb1_soma[["beta", "standard_error", "p_value"]], how="inner", rsuffix="_soma")
r = stats.pearsonr(klkb1.beta, klkb1.beta_soma)[0]
plt.figure(figsize=(4.5,4.5))
plt.axvline(x=0, color="black", linestyle="--", linewidth=.5)
plt.axhline(y=0, color="black", linestyle="--", linewidth=.5)
plt.errorbar(x=klkb1.beta, y=klkb1.beta_soma, xerr=klkb1.se, yerr=klkb1.standard_error, alpha=0.4, fmt='o', color="tab:blue")
plt.title("Effect of KLKB1 cis-variants on APOL1 expression", fontsize=12)
plt.xlabel("Beta in UKB AFR (Olink)")
plt.ylabel("Beta in AASK (SomaScan aptamer 2)")
plt.text(x=-0.53, y=0.45, s="r = {0:.3f}".format(r), fontsize=15)
ax = plt.gca()
ax.set_aspect('auto')
plt.xlim([-0.55,0.55])
plt.ylim([-0.55,0.55])
plt.tight_layout()
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_klkb1_effsize_afr_ID_11510_31.png", dpi=400)
plt.savefig("/Users/qsw/Desktop/misc_plots/apol_klkb1_effsize_afr_ID_11510_31.pdf", dpi=400)
plt.close()