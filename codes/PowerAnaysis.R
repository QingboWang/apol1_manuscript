# Power analysis given odds ratio of ESRD based on different risk alleles and plot

# Functions

# Two proportion Z-test: calculate the test statistic Z, 
# part of Z: a,b,c
calc_Zstat = function(n1, n2, p1, p2) {
a = 1/n1 + 1/n2
b = (n1*p1 + n2*p2)/(n1+n2)
c = 1 - b
Z = (p1 - p2) / sqrt(a*b*c) 
return(Z) }


# Calculate p1 given odds ratio, number of individuals, number of cases, risk allele frequency
cal_p1=function(N, n_case, f, R){ 
  # p1. is the risk genotype frequency in the total population, assume reccessive
  p1.=f^2
  # p.1 is the disease prevalence in the total population
  p.1=n_case/N
  # p11 is the proportion of people who have risk genotype and sick
  # p1 is the genotype frequency in the case
  # R is the odds ratio
  S = sqrt((1+(p1.+p.1)*(R-1))^2 + 4*R*(1-R)*p1.*p.1)
  p11 = (1+(p1.+p.1)*(R-1) - S) / (2*(R-1))
  p1 = p11 / p.1 
  return(p1)
}

#cal_p1(N=N, n_case=n_case, f=af_a, R=OR)

# Caclulate p2 given OR and p1
cal_p2<-function(OR, p1){
  p2=p1/(OR-OR*p1+p1)
  return(p2)
}

# Calculate power: given critical value and Z stat, get probability of P[X>x].
power.OR =function(N, n_case, OR, f){
  # get critical value given probability, one-side
  alpha = 0.05
  critical_value = qnorm(p=1-alpha, lower.tail=TRUE) #~1.65
  # calculate p1
  p1 = cal_p1(N=N, n_case=n_case, f=f, R=OR)
  # calculate p2
  p2 = cal_p2(OR=OR, p1=p1)
  # calculate z stat analytically
  Z_analy = calc_Zstat(n1=n_case, 
                     n2=N-n_case, 
                     p1=p1, 
                     p2=p2)
  power = pnorm(q=critical_value,  
                sd =1, mean = Z_analy, 
                lower.tail=FALSE)
  return(power)}


# plot the power-OR with different risk allele frequency
plot_power.OR=function(N, n_case, frq_G1, frq_G2, frq_G1G2, frq_V1){
# create the table with RAF, OR, calculate the power given these two
RAF = c(frq_G1, frq_G2, frq_G1G2, frq_V1)
OR <- seq(from = 1.01, to = 5, length = 50)
OR_power <- expand.grid(OR = OR, RAF = RAF)
OR_power$power <- NA

for(i in 1:nrow(OR_power)){
   OR_power$power[i] <- power.OR(N=N,
                                 n_case=n_case, 
                                 OR=OR_power$OR[i],
                                 f=OR_power$RAF[i])
}  

# plot the OR vs power
library(ggplot2)
# Basic line plot with points
p=ggplot(data=OR_power, aes(x=OR, y=power, group=RAF, color=as.factor(RAF))) +
 scale_color_manual(labels = c(paste0("G2, RAF=", frq_G2),  paste0("G1M, RAF=", frq_G1), paste0("G1/G2, RAF=", frq_G1G2), paste0("rs2239785, RAF=", frq_V1)), 
                     values = c("#9467bd",  "#d62728", "#0D98BA", "#2ca02c")) +
  geom_line(linewidth=0.7)+
  geom_point(size=0.9)+
  theme_bw() +
  xlab("Odds ratio") + 
  ylab("Power") + 
  ylim(c(0, 1)) +
  scale_x_continuous(breaks=seq(1,5,1)) +
  guides(color=guide_legend(title="Risk Allele")) # change legend title
return(p)
}



# --------------UKBB AFR-----------------

# ukbb data Risk allele frequency
frq_G1_ukbb=0.29
frq_G2_ukbb=0.15
frq_V1_ukbb=0.70
frq_G1G2_ukbb=0.44
N_ukbb=6407
n_case_ukbb=68

p1=plot_power.OR(N=N_ukbb, 
	n_case=n_case_ukbb, 
	frq_G1=frq_G1_ukbb, 
	frq_G2=frq_G2_ukbb, 
	frq_G1G2=frq_G1G2_ukbb, 
	frq_V1=frq_V1_ukbb)
p_ukbb_AFR = p1 + ggtitle("Power analysis for ESRD in UKBB AFR") 
#ggsave(plot=p_ukbb_AFR, filename="../Figure/PowerAnalysis_ESRD_UKBB_AFR.png", width=6, height=4, dpi=300, units = "in")
ggsave(plot=p_ukbb_AFR, filename="../Figure/PowerAnalysis_ESRD_UKBB_AFR.pdf", width=6, height=4, dpi=500, units = "in")

# --------------PMBB AFR-----------------

frq_G1_pmbb=0.22
frq_G2_pmbb=0.14
frq_V1_pmbb=0.63
frq_G1G2_pmbb=0.36
N_pmbb=10908
n_case_pmbb=911

p2=plot_power.OR(N=N_pmbb, 
	n_case=n_case_pmbb, 
	frq_G1=frq_G1_pmbb, 
	frq_G2=frq_G2_ukbb, 
	frq_G1G2=frq_G1G2_pmbb, 
	frq_V1=frq_V1_pmbb)
p_pmbb_AFR = p2 + ggtitle("Power analysis for ESRD in PMBB AFR") 
#ggsave(plot=p_pmbb_AFR, filename="../Figure/PowerAnalysis_ESRD_PMBB_AFR.png", width=6, height=4, dpi=300, units = "in")
ggsave(plot=p_pmbb_AFR, filename="../Figure/PowerAnalysis_ESRD_PMBB_AFR.pdf", width=6, height=4, dpi=500, units = "in")


# --------------UKBB EUR V1-----------------

frq_V1=0.19 #ALLELE FREQ OF V1 in EUR
N=410372 # EUR in UKBB
n_case=1336 #ESRD in EUR in UKBB

# create the table with RAF, OR, calculate the power given these two
RAF = c(frq_V1)
OR <- seq(from = 1.01, to = 2, length = 50)
OR_power <- expand.grid(OR = OR, RAF = RAF)
OR_power$power <- NA

for(i in 1:nrow(OR_power)){
   OR_power$power[i] <- power.OR(N=N,
                                 n_case=n_case, 
                                 OR=OR_power$OR[i],
                                 f=OR_power$RAF[i])
}  

# plot the OR vs power
library(ggplot2)
# Basic line plot with points
p_ukb_esrd_eur=ggplot(data=OR_power, aes(x=OR, y=power, group=RAF, color=as.factor(RAF))) +
 scale_color_manual(labels = c("rs2239785, RAF=0.19"), 
                     values = c("#2ca02c")) +
  geom_line()+
  geom_point(size=1)+
  theme_bw() +
  xlab("Odds ratio") + 
  ylab("Power") + 
  ggtitle("Power analysis for ESRD in UKBB EUR") +
  ylim(c(0, 1)) +
  scale_x_continuous(breaks = c(seq(1,2,0.2))) +
  guides(color=guide_legend(title="Risk Allele")) # change legend title
#ggsave(plot=p_ukb_esrd_eur, filename="../Figure/PowerAnalysis_ESRD_UKBB_EUR_V1.png", width=6, height=4, dpi=300, units = "in")
ggsave(plot=p_ukb_esrd_eur, filename="../Figure/PowerAnalysis_ESRD_UKBB_EUR_V1.pdf", width=6, height=4, dpi=500, units = "in")


# --------------AoU AFR-----------------
# for AoU data, using the freq from https://genetics.opentargets.org/
# given a convservative estimate for disease prevalence
frq_G1_aou=0.22
frq_G2_aou=0.14
frq_V1_aou=0.74
frq_G1G2_aou=0.36
N_aou=50064
#p_case_lower=4400/254700 #prevalence in the total population (who have EHR), as a lower bound for conservative estimation, assuming prevalence is the same across populations
#n_case_aou=round(N*p_case_lower, digits = 0) #ESRD CASES
n_case_aou=911

p3=plot_power.OR(N=N_aou, 
	n_case=n_case_aou, 
	frq_G1=frq_G1_aou, 
	frq_G2=frq_G2_aou, 
	frq_G1G2=frq_G1G2_aou, 
	frq_V1=frq_V1_aou)
p_aou_AFR = p3 + ggtitle("Power analysis for ESRD in AoU AFR") 
#ggsave(plot=p_aou_AFR, filename="../Figure/PowerAnalysis_ESRD_AoU_AFR.png", width=6, height=4, dpi=300, units = "in")
ggsave(plot=p_aou_AFR, filename="../Figure/PowerAnalysis_ESRD_AoU_AFR.pdf", width=6, height=4, dpi=500, units = "in")



