# Script to model and analyze consumer-resource size and size-ratio data across communities
# Note that predictions for feasible size ratios will change a bit every time this script is run because if the random seed is not set for the beta distribution sampling. 

require(grDevices) # for colours
library(ggplot2)
library(reshape2)

# options(error = browser) # for interactive debugging

rm(list=ls()) # clear objects
graphics.off() #close all open figures and graphics objects
IDir <- "../Data/"
ODir <- "../Tables+Figs/"
DataName <- "AllWebCRs"
# DataName <- "AllWebCRsReduc"

############################################
################ Functions #################
############################################

scal_pars <- function(data_m_C,data_m_R,z_0,r_0,beta_r,B_0_2D,B_0_3D,x_0_2D, x_0_3D, beta_x, h_02D, h_03D,beta_a, ForStr,a_02D,a_03D,p_v,p_d,gamma_2D,gamma_3D) { # returns scaling parameters given scalar-, vector- or matrix-valued consumer and resource mass data (inputs "data_m_C", "data_m_R", both of same length)
	
	z <- z_0 * data_m_C^(beta_r-1) # mortality rate
	r <- r_0 * data_m_R^(beta_r-1) # production rate
	LossRate_2D <- B_0_2D * data_m_C^beta_r # Biomass loss rate
	LossRate_3D <- B_0_3D * data_m_C^beta_r
	x_2D <- x_0_2D * (data_m_R ^ -beta_x) #~ Numerical abundance (m^-2) according to Peters; Also see Leaper and Raffaeli (1999) Ecol Lett and Cyr et al (1997). Oikos, 79(2), 333-346
	x_3D <- x_0_3D * (data_m_R ^ -beta_x) # 

	K_2D <- x_2D * data_m_R # Biomass carrying capacity; i.e., K = x_0 * modl_m_R^(1-beta_x)
	K_3D <- x_3D * data_m_R

	data_k_RC <- data_m_R/data_m_C

	h_2D <- h_02D * data_m_R * data_m_C^-beta_a # per-capita handling time
	h_3D <- h_03D * data_m_R * data_m_C^-beta_a #

	if (ForStr == 'A') { #~ Active
		a_2D <- (a_02D * data_m_C ^ (p_v + (2 * p_d * (2 - 1)))) * (1 + data_k_RC ^ (2 * p_v)) ^ 0.5 * data_k_RC ^ ((2 - 1) * p_d) / (1 + data_k_RC ^ gamma_2D) 
		a_3D <- (a_03D * data_m_C ^ (p_v + (2 * p_d * (3 - 1)))) * (1 + data_k_RC ^ (2 * p_v)) ^ 0.5 * data_k_RC ^ ((3 - 1) * p_d) / (1 + data_k_RC^gamma_3D)
	} else if (ForStr == 'S') {#~ Sit wait, mass-specific
		a_2D <- (a_02D * data_m_C ^ (p_v + (2 * p_d * (2 - 1)))) * data_k_RC ^ (p_v + (p_d * (2 - 1))) / (1 + data_k_RC ^ gamma_2D)
		a_3D <- (a_03D * data_m_C ^ (p_v + (2 * p_d * (3 - 1)))) * data_k_RC ^ (p_v + (p_d * (3 - 1))) / (1 + data_k_RC ^ gamma_3D)
	} else if (ForStr == 'G') {#~ Grazer, mass-specific
		
		a_2D <- (a_02D * data_m_C ^ (p_v + (2 * p_d * (2 - 1)))) * data_k_RC ^ (p_d * (2 - 1)) / (1 + data_k_RC ^ gamma_2D)   

		a_3D <- a_03D * (data_m_C ^ (p_v + (2 * p_d * (3 - 1)))) * (data_k_RC ^ (p_d * (3 - 1))) / (1 + data_k_RC ^ gamma_3D)
		}

	par_list <- list("z" = z, "r" = r, "LossRate_2D" = LossRate_2D, "LossRate_3D" = LossRate_3D,"x_2D" = x_2D,"x_3D" = x_3D,"K_2D" = K_2D,"K_3D" = K_3D,"h_2D" = h_2D,"h_3D" = h_3D,"a_2D" = a_2D,"a_3D" = a_3D)

	return(par_list) 

}

get_sizehist <- function(data, brks){ #returns normalized histogram for size and size ratio
	tmpHist <- hist(data,breaks = brks,plot = FALSE)
	tmpHist$counts <- tmpHist$counts/max(tmpHist$counts) #normalize [0,1]
	
	return(tmpHist)
}

############################################
############################################

## Assign model parameter values

Bet_a <- 1 #Beta distribution parameter 1
Bet_b <- 1 #Beta distribution parameter 2 (setting both = 1 gives uniform, which is what we use)

r_0 <- (exp(23.82)*exp(-0.65/((8.617*10^-5)*(273.15+20))))/(1000 * (24*3600)) #constant for scaling of intrinsic growth rate, obtained from MTE papers assuming t_ref of 20

# B_0_2D <- 3/(7*10^06) # Intercept for scaling of metabolic rate (originally expressed in J/s, coverted to kg/s) from Peters, also see Carbone et al 2014 Ecol Lett and and Rizzuto et al 2017

B_0_2D <- 9*1000/(24*3600*7*10^6) # Nagy 2005 for field, corrected for temperature and converted to mass units

B_0_3D <- B_0_2D

z_0 <- B_0_2D

h_02D <- 10000 #constant for handling time
h_03D <- 1000 #constant for handling time
eff <- 0.3 #Consumer's biomass conversion efficiency
beta_a <- 0.9 #Body mass-metabolism scaling exponent for active/maximal state
beta_r <- 0.75 #Body mass-metabolism scaling exponent for resting state
beta_x <- 0.75 #~exponent for abundance scaling

p_v <- 0.26 #velocity scaling exponent
p_d <- 0.20 #detection scaling exponent

x_0_2D <- 1 #~ 2D Resource carrying capacity constant; 1 is approx intercept at 1 kg for biomass abundance scaling from Peters book, terrestrial inverts 
x_0_3D <- 300 # 3D Resource carrying capacity constant; 300 is approx intercept at 1 kg for biomass abundance scaling from Peters book, aquat inverts

a_02D <- 10^-3.5 #2D constant for scaling of search rate
a_03D <- 10^-1.8 #3D constant for scaling of search rate
gamma_2D <- 2
gamma_3D <- 2

ForStr <- 'G'

Resol <- 300 # for plotting
RdmPts <- 10000 # for simulations

StabCrit <- 1

#Create string containing parameter names and values to include in output file name:
ParaName <- paste(ForStr,'_','z0_',as.character(z_0),'_','r0_',as.character(r_0),'_','h0_',as.character(h_02D),'and',as.character(h_03D), '_', 'gamma_',as.character(gamma_2D),'and', as.character(gamma_3D), '_', 'x0_', as.character(x_0_2D),'and',as.character(x_0_3D), sep = '')

Data <- as.matrix(read.table(paste(IDir,DataName,".csv", sep = ""), header = TRUE, fill = TRUE,sep = ","))

ComNames <- unique(Data[,"Community"]) 
m_C_all <- as.matrix(as.numeric(Data[,"Cmass"]))#Extract consumer sizes
m_R_all <- as.matrix(as.numeric(Data[,"Rmass"]))#Extract resource sizes

# Set histogram limits based on real data

kRCLims <- c(floor(min(log10(m_R_all/m_C_all))),ceiling(max(log10(m_R_all/m_C_all))))
kBins <- seq(kRCLims[1],kRCLims[2],len=80) #pre-designate bins for overall as well as each community

mLims_all <- c(floor(min(log10(c(m_C_all,m_R_all)))),ceiling(max(log10(c(m_C_all,m_R_all)))))
mBins <- seq(mLims_all[1],mLims_all[2],len=30) #pre-designate bins for overall as well as each community

ReSummNames <- c("Community", "Taxa", "2D Taxa", "3D Taxa","Links","2D Links","3D Links", #Headers for for tabulated results
	"D Overlap (Consumers)", "D Overlap (Resources)",
	"Mean 2D log10(kRC)","95pctCI Mean 2D log10(kRC)", 
	"Mean 3D log10(kRC)","95pctCI Mean 3D log10(kRC)", 
	"Displacement Mean log10(kRC)", "P-val Rmz Displacement Mean log10(kRC)","Predict Displacement Mean log10(kRC) - stab", 
	"Predict Displacement Mean log10(kRC) - coex", 
	"Predict Displacement Mean log10(kRC) - ener", 
	"Median 2D log10(kRC)", "Median 3D log10(kRC)", 
	"Displacement Median log10(kRC)", "P-val Displacement Median log10(kRC)",
	"Displacement estimate Median log10(kRC)", "95pctCI Displacement estimate Median log10(kRC)",
	"P-val Rmz Displacement Median log10(kRC)",
	"Predict Displacement Median log10(kRC) - stab", 
	"Predict Displacement Median log10(kRC) - coex", 
	"Predict Displacement Median log10(kRC) - ener", 
	"Mean 2D log10(Body mass)", "95pctCI Mean 2D log10(Body mass)", 
	"Mean 3D log10(Body mass)", "95pctCI Mean 3D log10(Body mass)", 
	"Median 2D log10(Body mass)","Median 3D log10(Body mass)",
	"Displacement Median log10(Body mass)", "P-val Displacement Median log10(Body mass)", 
	"Displacement estimate Median log10(Body mass)", "95pctCI Displacement estimate Median log10(Body mass)")
				
ReSumm <- matrix(data = "",length(ComNames)+1,length(ReSummNames),dimnames = list(NULL,ReSummNames)) #Preallocate empty array for tabulated results

pdf(paste(ODir,DataName,"_Fig1.pdf",sep=""), 11.7, 8.3)

par(mfcol=c(length(ComNames)+1,2),lwd=1.5,mar=c(0,2,1,1),bty="l",oma=c(0,0,0,0)) #initialize multi-paneled plot

set.seed(123)

for(i in 1:(length(ComNames)+1)){ #loop to run analysis and plotting for each community
	if (i == 1){tmpData <- Data; ComName <- "All"} #all communities together
	else {tmpData <- Data[which(Data[,"Community"]==ComNames[i-1]),] ;ComName <- ComNames[i-1]} #Otherwise extract data for current community

	m_C <- as.matrix(as.numeric(tmpData[,"Cmass"]))
	m_R <- as.matrix(as.numeric(tmpData[,"Rmass"]))
	m_C2D <- m_C[which(tmpData[,"D"]=="2")] #Extract 2D size data 
	m_R2D <- m_R[which(tmpData[,"D"]=="2")] #Extract 2D size data 
	m_C3D <- m_C[which(tmpData[,"D"]=="3")] #Extract 3D size data 
	m_R3D <- m_R[which(tmpData[,"D"]=="3")] #Extract 3D size data 
	kRCs <- as.matrix(m_R/m_C)#Extract size-ratios
	kRCs2D <- m_R2D/m_C2D #Extract 2D size-ratios 
	kRCs3D <- m_R3D/m_C3D #Extract 3D size-ratios
	Con2D <- unique(tmpData[which(tmpData[,"D"]=="2"),"Consumer"]) #get 2D consumer species
	Res2D <- unique((tmpData[which(tmpData[,"D"]=="2"),"Resource"])) #get 2D resource species
	Con3D <- unique((tmpData[which(tmpData[,"D"]=="3"),"Consumer"])) #get 3D consumer species
	Res3D <- unique((tmpData[which(tmpData[,"D"]=="3"),"Resource"])) #get 3D resource species
	Perm_kRC <- matrix(data = "",RdmPts,2,dimnames = list(NULL,c("Mean","Median"))) #Preallocate empty array for randomized deviations
	for (j in 1:RdmPts){ #generate randomized lists to calculate significance of observed 2D-3D differences
		Perm2Dk_RC <- sample(c(m_R,m_C),length(m_R2D),replace = TRUE)/sample(c(m_R,m_C),length(m_C2D),replace = TRUE) #
		Perm3Dk_RC <- sample(c(m_R,m_C),length(m_R3D),replace = TRUE)/sample(c(m_R,m_C),length(m_C3D),replace = TRUE) #
		Perm_kRC[j,"Mean"] <- mean(log10(Perm3Dk_RC)) - mean(log10(Perm2Dk_RC))
		Perm_kRC[j,"Median"] <- median(log10(Perm3Dk_RC)) - median(log10(Perm2Dk_RC))
	}

	#Size limits for current community
	mLims <- c(floor(min(log10(c(m_C,m_R)))),ceiling(max(log10(c(m_C,m_R)))))#calculate min and max log10 sizes
	m_RLims <- c(floor(min(log10(m_R))),ceiling(max(log10(m_R))))#calculate resource log10 min-max	sizes 
	m_CLims <- c(floor(min(log10(m_C))),ceiling(max(log10(m_C))))#calculate consumer log10 min-max	sizes

	# Store results for this part of analysis

	ReSumm[i,"Community"] <- ComName
	ReSumm[i,"Taxa"] <- length(unique(c(tmpData[,"Consumer"],tmpData[,"Resource"])))
	ReSumm[i,"2D Taxa"] <- length(unique(c(tmpData[which(tmpData[,"D"]=="2"),"Consumer"],tmpData[which(tmpData[,"D"]=="2"),"Resource"])))
	ReSumm[i,"3D Taxa"] <- length(unique(c(tmpData[which(tmpData[,"D"]=="3"),"Consumer"],tmpData[which(tmpData[,"D"]=="3"),"Resource"])))
	ReSumm[i,"Links"] <- length(kRCs)
	ReSumm[i,"2D Links"] <- length(kRCs2D)
	ReSumm[i,"3D Links"] <- length(kRCs3D)
	ReSumm[i,"D Overlap (Consumers)"] <- length(intersect(Con2D,Con3D))/length(union(Con2D,Con3D))#Jaccard overlap
	ReSumm[i,"D Overlap (Resources)"] <- length(intersect(Res2D,Res3D))/length(union(Res2D,Res3D)) #Jaccard overlap
	ReSumm[i,"Mean 2D log10(kRC)"] <- mean(log10(kRCs2D))
	ReSumm[i,"95pctCI Mean 2D log10(kRC)"] <- qnorm(1-0.05/2)*sd(log10(kRCs2D))/sqrt(length(kRCs2D)) #95pctCI
	ReSumm[i,"Mean 3D log10(kRC)"] <- mean(log10(kRCs3D))
	ReSumm[i,"95pctCI Mean 3D log10(kRC)"] <- qnorm(1-0.05/2)*sd(log10(kRCs3D))/sqrt(length(kRCs3D)) #95pctCI
	ReSumm[i,"Displacement Mean log10(kRC)"] <- mean(log10(kRCs3D))-mean(log10(kRCs2D)) #displacement	
	ReSumm[i,"P-val Rmz Displacement Mean log10(kRC)"] <- length(which(as.numeric(Perm_kRC[,"Mean"]) < as.numeric(ReSumm[i,"Displacement Mean log10(kRC)"])))/RdmPts 
	ReSumm[i,"Median 2D log10(kRC)"] <- median(log10(kRCs2D))
	ReSumm[i,"Median 3D log10(kRC)"] <- median(log10(kRCs3D))
	if (length(kRCs2D)>0 && length(kRCs3D)>0){ #if there were both 2D and 3D interactions
		ReSumm[i,"Displacement Median log10(kRC)"] <- median(log10(kRCs3D))-median(log10(kRCs2D)) #displacement
		sigtest <- wilcox.test(log10(kRCs3D), log10(kRCs2D), conf.int = TRUE) #two-sample wilcox test
		ReSumm[i,"P-val Displacement Median log10(kRC)"] <- sigtest$p.value
		ReSumm[i,"Displacement estimate Median log10(kRC)"] <- sigtest$estimate
		ReSumm[i,"95pctCI Displacement estimate Median log10(kRC)"] <- sigtest$conf.int[2]-sigtest$estimate #95pct CI of location displacement (pseudomedian)
		ReSumm[i,"P-val Rmz Displacement Median log10(kRC)"] <- length(which(as.numeric(Perm_kRC[,"Median"]) < as.numeric(ReSumm[i,"Displacement Median log10(kRC)"])))/RdmPts 
	}
	else {ReSumm[i,c("Displacement Median log10(kRC)","P-val Displacement Median log10(kRC)", "Displacement estimate Median log10(kRC)",
	"95pctCI Displacement estimate Median log10(kRC)","95pctCI Displacement estimate Median log10(kRC)")] <- "NA"}

	############ Perform body size analysis and store results  ############
	
	ReSumm[i,"Mean 2D log10(Body mass)"] <- mean(log10(unique(c(m_C2D,m_R2D))))
	ReSumm[i,"95pctCI Mean 2D log10(Body mass)"] <- qnorm(1-0.05/2)*sd(log10(unique(c(m_C2D,m_R2D))))/sqrt(length(unique(c(m_C2D,m_R2D)))) #95pctCI
	ReSumm[i,"Mean 3D log10(Body mass)"] <- mean(log10(unique(c(m_C3D,m_R3D))))
	ReSumm[i,"95pctCI Mean 3D log10(Body mass)"] <- qnorm(1-0.05/2)*sd(log10(unique(c(m_C3D,m_R3D))))/sqrt(length(unique(c(m_C3D,m_R3D)))) #95pctCI
	ReSumm[i,"Median 2D log10(Body mass)"] <- median(log10(unique(c(m_C2D,m_R2D))))
	ReSumm[i,"Median 3D log10(Body mass)"] <- median(log10(unique(c(m_C3D,m_R3D))))
	if (length(kRCs2D)>0 && length(kRCs3D)>0){
		ReSumm[i,"Displacement Median log10(Body mass)"] <- median(log10(unique(c(m_C3D,m_R3D))))-median(log10(unique(c(m_C2D,m_R2D)))) #displacement

		tmp3Dtst <- unique(c(m_C3D,m_R3D))
		tmp2Dtst <- unique(c(m_C2D,m_R2D))
		intersectTst <- intersect(tmp3Dtst,tmp2Dtst)
		mark <- rep(0,length(tmp3Dtst)) #mark masses overlapping between 2D and 3D (cannot test for differences otherwise)
		for (k in 1:length(tmp3Dtst)){
			for (l in 1:length(intersectTst)){
				if(tmp3Dtst[k] == intersectTst[l]) {
					mark[k] = 1
				}
			} 
		}
		tmp3Dtst <- tmp3Dtst[which(mark == 0)] #remove overlapping masses

		mark <- rep(0,length(tmp2Dtst))
		for (k in 1:length(tmp2Dtst)){
			for (l in 1:length(intersectTst)){
				if(tmp2Dtst[k] == intersectTst[l]) {
					mark[k] = 1
				}
			} 
		}
		tmp2Dtst <- tmp2Dtst[which(mark == 0)]#remove overlapping masses
   
    if (length(tmp2Dtst)>0 && length(tmp3Dtst)>0){ #only if there are still separate sets of masses after removing overlap
      sigtest <- wilcox.test(log10(tmp3Dtst), log10(tmp2Dtst), conf.int = TRUE) #two-sample wilcox test
      }
    else {
    ReSumm[i,c("Displacement Median log10(Body mass)","P-val Displacement Median log10(Body mass)", "Displacement estimate Median log10(Body mass)","95pctCI Displacement estimate Median log10(Body mass)")] <- "NA"
    }
 
		ReSumm[i,"P-val Displacement Median log10(Body mass)"] <- sigtest$p.value
		ReSumm[i,"Displacement estimate Median log10(Body mass)"] <- sigtest$estimate
		ReSumm[i,"95pctCI Displacement estimate Median log10(Body mass)"] <- sigtest$conf.int[2]-sigtest$estimate #95pct CI of location displacement (pseudomedian)
	} 
	else {ReSumm[i,c("Displacement Median log10(Body mass)", "P-val Displacement Median log10(Body mass)", "Displacement estimate Median log10(Body mass)","95pctCI Displacement estimate Median log10(Body mass)")] <- "NA"}

	############################################
	######### Random interaction model #########
	############################################

	modl_m_C <- 10^as.matrix((m_CLims[1] + ((m_CLims[2]-m_CLims[1])*rbeta(RdmPts,Bet_a,Bet_b)))) #generate random consumer vector
	modl_m_R <- 10^as.matrix((m_RLims[1] + ((m_RLims[2]-m_RLims[1])*rbeta(RdmPts,Bet_a,Bet_b)))) #generate random resource vector

	par_list <- scal_pars(modl_m_C,modl_m_R,z_0,r_0,beta_r,B_0_2D,B_0_3D,x_0_2D, x_0_3D, beta_x, h_02D, h_03D,beta_a, ForStr,a_02D,a_03D,p_v,p_d,gamma_2D,gamma_3D)
	
	z <- par_list$z; r <- par_list$r; LossRate_2D <- par_list$LossRate_2D; LossRate_3D <- par_list$LossRate_3D; x_2D <- par_list$x_2D; x_3D <- par_list$x_3D; K_2D <- par_list$K_2D; K_3D <- par_list$K_3D; h_2D <- par_list$h_2D; h_3D <- par_list$h_3D; a_2D <- par_list$a_2D; a_3D <- par_list$a_3D

  #%%%%%%%%%%%%%%%%%%%%% Energetic bounds %%%%%%%%%%%%%%%%%%%%% 

	modl_k_RC <- modl_m_R/modl_m_C #Size-ratio

  	CRate_2D <- modl_m_R * a_2D * x_2D / (1 + h_2D * a_2D * x_2D) #~Search-constrained biomass consumption rate 
	CRate_3D <- modl_m_R * a_3D * x_3D / (1 + h_3D * a_3D * x_3D) 

	NetRate_2D <- eff * CRate_2D - LossRate_2D
	NetRate_3D <- eff * CRate_3D - LossRate_3D

# %%%%%%%%%%%%%%%%%%% Coexistence Bounds %%%%%%%%%%%%%%%%%%%% 

	a_2D_ <- a_2D/modl_m_C # get mass-specific search rate (denoted by a' in manuscript)
	a_3D_ <- a_3D/modl_m_C 
	Coex2D <- K_2D - (z / (a_2D_* (eff - z * h_2D))) # coexistence condition
	Coex3D <- K_3D - (z / (a_3D_* (eff - z * h_3D)))
	
# %%%%%%%%%%%%%%%%%%%%% Stability bounds %%%%%%%%%%%%%%%%%%%%% 

	Stab2D <- matrix(0,RdmPts,1); Stab3D <- matrix(0,RdmPts,1)
	J11_2D <- (r*z*(K_2D*(h_2D^2)*z*a_2D_ + eff + h_2D*(z - K_2D*a_2D_*eff)))/(eff*K_2D*a_2D_*(h_2D*z - eff)) #1st element of jacobian matrix 
	J12_2D <- -z/eff # 2nd element
	J21_2D <- r*(eff - z/(K_2D*a_2D_)- h_2D*z ) # 3rd element
	J22_2D <- matrix(0,RdmPts,1) # 4th element
	J11_3D <- (r*z*(K_3D*(h_3D^2)*z*a_3D_ + eff + h_3D*(z - K_3D * a_3D_ * eff))) / (K_3D * a_3D_ * (h_3D*z - eff) * eff)
	J12_3D <- -z/eff
	J21_3D <- r*(-h_3D*z - z/(K_3D*a_3D_)+ eff)
	J22_3D <- matrix(0,RdmPts,1)
	for (k in 1:RdmPts){
		tmpEig2D <- eigen(rbind(c(J11_2D[k],J12_2D[k]),c(J21_2D[k],J22_2D[k])),symmetric = FALSE,only.values = TRUE)
		tmpEig3D <- eigen(rbind(c(J11_3D[k],J12_3D[k]),c(J21_3D[k],J22_3D[k])),symmetric = FALSE,only.values = TRUE)
		if (max(Re(tmpEig2D$values)) < 0 && all(Im(tmpEig2D$values)== 0)) {Stab2D[k] = 4} #Stable point
		else if (max(Re(tmpEig2D$values)) < 0 && any(Im(tmpEig2D$values)>0)) {Stab2D[k] = 3} #Stable point with damped oscillations
		else if (all(Re(tmpEig2D$values) == 0) && any(Im(tmpEig2D$values)>0)) {Stab2D[k] = 2} # neutrally stable with oscillations
		else if (any(Re(tmpEig2D$values)>0)	&& any(Im(tmpEig2D$values)>0)) {Stab2D[k] = 1} # unstable with oscillations
		else if (any(Re(tmpEig2D$values)>0)	&& all(Im(tmpEig2D$values)==0)) {Stab2D[k] = 0} # unstable
		else if (all(Re(tmpEig2D$values)==0) && all(Im(tmpEig2D$values)==0)) {Stab2D[k] = 0} #?
		
		if (max(Re(tmpEig3D$values)) < 0 && all(Im(tmpEig3D$values)== 0)) {Stab3D[k] = 4} #Stable point
		else if (max(Re(tmpEig3D$values)) < 0 && any(Im(tmpEig3D$values)> 0)) {Stab3D[k] = 3} #Stable point with oscillations
		else if (all(Re(tmpEig3D$values) == 0) && any(Im(tmpEig3D$values)>0)) {Stab3D[k] = 2} #neutrally stable with oscillations
		else if (any(Re(tmpEig3D$values)> 0)	&& any(Im(tmpEig3D$values)>0)) {Stab3D[k] = 1} #unstable with oscillations
		else if (any(Re(tmpEig3D$values)>0)	&& all(Im(tmpEig3D$values)==0)) {Stab3D[k] = 0} #unstable
		else if (all(Re(tmpEig3D$values)==0)	&& all(Im(tmpEig3D$values)==0)) {Stab3D[k] = 0} #?
	}

	modl_k_RC_2DFilt0 <- modl_k_RC[NetRate_2D > 0]
	modl_k_RC_3DFilt0 <- modl_k_RC[NetRate_3D > 0]

  modl_k_RC_2DFilt <- modl_k_RC[Coex2D > 0]
  modl_k_RC_3DFilt <- modl_k_RC[Coex3D > 0]

	modl_k_RC_2DFilt1 <- modl_k_RC[Stab2D > StabCrit] 
	modl_k_RC_3DFilt1 <- modl_k_RC[Stab3D > StabCrit]

	## Store predictions

	ReSumm[i,"Predict Displacement Mean log10(kRC) - stab"] <- mean(log10(modl_k_RC_3DFilt1)) - mean(log10(modl_k_RC_2DFilt1))

	ReSumm[i,"Predict Displacement Mean log10(kRC) - coex"] <- mean(log10(modl_k_RC_3DFilt)) - mean(log10(modl_k_RC_2DFilt))
	
	ReSumm[i,"Predict Displacement Mean log10(kRC) - ener"] <- mean(log10(modl_k_RC_3DFilt0)) - mean(log10(modl_k_RC_2DFilt0))
 
	ReSumm[i,"Predict Displacement Median log10(kRC) - stab"] <- median(log10(modl_k_RC_3DFilt1)) - median(log10(modl_k_RC_2DFilt1))

	ReSumm[i,"Predict Displacement Median log10(kRC) - coex"] <- median(log10(modl_k_RC_3DFilt)) - median(log10(modl_k_RC_2DFilt))	
	
	ReSumm[i,"Predict Displacement Median log10(kRC) - ener"] <- median(log10(modl_k_RC_3DFilt0)) - median(log10(modl_k_RC_2DFilt0))
	
	###################################################################
	####################### Plot distributions ########################
	###################################################################
		
	if (length(kRCs2D)>0 && length(kRCs3D)>0) { # only if there were 2D as well as 3D interactions
		kRCHist <- get_sizehist(log10(c(kRCs2D,kRCs3D)), kBins)
		kRCHist2D <- get_sizehist(log10(kRCs2D), kBins)
		kRCHist3D <- get_sizehist(log10(kRCs3D), kBins)

		mHist <- get_sizehist(log10(c(m_C2D,m_C3D)), kBins)		
		mHist2D <- get_sizehist(log10(m_C2D), kBins)
		mHist3D <- get_sizehist(log10(m_C3D), kBins)		
	}
	else if (length(kRCs3D)==0) {
		kRCHist <- get_sizehist(log10(kRCs2D), kBins)
		kRCHist2D <- get_sizehist(log10(kRCs2D), kBins)
		kRCHist3D <- get_sizehist(0, kBins)
		
		mHist <- get_sizehist(log10(m_C2D), kBins)
		mHist2D <- get_sizehist(log10(m_C2D), kBins)
		mHist3D <- get_sizehist(0, kBins)
	}
	else if (length(kRCs2D)==0) {
		kRCHist <- get_sizehist(log10(kRCs3D), kBins)
		kRCHist2D <- get_sizehist(0, kBins)
		kRCHist3D <- get_sizehist(log10(kRCs3D), kBins)

		mHist <- get_sizehist(log10(m_C3D), kBins)
		mHist2D <- get_sizehist(0, kBins)
		mHist3D <- get_sizehist(log10(m_C3D), kBins)
	}

	# plot 
	par(mfg = c(i,1)); plot(kRCHist2D,xlim=kRCLims,ylim=c(0,1),col = rgb(1, 0, 0, 0.5),tck = 0.05, main = NULL)
	par(mfg = c(i,1)); plot(kRCHist3D,xlim=kRCLims,ylim = c(0,1),col = rgb(0, 0, 1, 0.5), tck = 0.05, main = ComName)
		
	par(mfg = c(i,1)); arrows(as.numeric(ReSumm[i,"Median 2D log10(kRC)"]), c(0,1)[2], as.numeric(ReSumm[i,"Median 2D log10(kRC)"]), 0, length = 0.2, angle = 30, code = 2, col = "red")
	par(mfg = c(i,1)); arrows(as.numeric(ReSumm[i,"Median 3D log10(kRC)"]), c(0,1)[2], as.numeric(ReSumm[i,"Median 3D log10(kRC)"]), 0, length = 0.2, angle = 30, code = 2, col = "blue")
	par(mfg = c(i,1)); text(as.numeric(ReSumm[i,"Median 3D log10(kRC)"]), c(0,1)[2],format(as.numeric(ReSumm[i,"Displacement Median log10(kRC)"]),digits = 3))

	par(mfg = c(i,2)); plot(mHist2D,xlim=mLims,ylim=c(0,1),col = rgb(1, 0, 0, 0.5),tck = 0.05, main = NULL)
	par(mfg = c(i,2)); plot(mHist3D,xlim=mLims,ylim = c(0,1),col = rgb(0, 0, 1, 0.5),tck = 0.05,main = ComName)

}

dev.off()

## Write all results to file
write.csv(ReSumm,paste(ODir,DataName,"_Table2_",ParaName,'.csv',sep = ""),row.names = FALSE)

#############################################################
##################### Additional Plots ######################
#############################################################

###### Raw size-ratio data histograms #####

DataToPlot <-  data.frame("Web" = Data[,"Community"], "SizeRatios" = as.numeric(Data[,"Log10CRBR"]))

p2 <- ggplot(DataToPlot, aes(x = SizeRatios)) + 
		geom_histogram(binwidth = .5) +
		xlab("Log10(Body mass ratio)") 
#~		facet_wrap( ~ Web, scales = "free",ncol = 3)
	
pdf(paste(ODir,DataName,"_AllSizeRatios.pdf",sep = ""))
print(p2)
dev.off()

p3 <- ggplot(DataToPlot, aes(x = SizeRatios)) + 
#~		geom_histogram(breaks = seq(-14,10,len=50)) +
		geom_histogram(binwidth = .5) +
		xlab("Log10(Body mass ratio)") +
		facet_wrap( ~ Web, scales = "free",ncol = 3)
			
pdf(paste(ODir,DataName,"_CommunitySizeRatios.pdf",sep = ""))
print(p3)
dev.off()

##### m_C vs m_R data with overall predicted bounds #####

## Prepare data

DataToPlot <- as.data.frame(Data[,c("Cmass","Rmass","D")],stringsAsFactors = FALSE) #The data to plot
colnames(DataToPlot) <- c('m_C', 'm_R','D')
#~DataToPlot$PointColor <- rep('',nrow(DataToPlot)) #Add color column
DataToPlot$m_C <- log10(as.numeric(DataToPlot$m_C))
DataToPlot$m_R <- log10(as.numeric(DataToPlot$m_R))
DataToPlot$D[DataToPlot$D == 2] <- '2D' #factors have to have exactly the same name as ModelToPlot
#~DataToPlot$PointColor[DataToPlot$D == '2D'] <- 'red'
DataToPlot$D[DataToPlot$D == 3] <- '3D' 
#~DataToPlot$PointColor[DataToPlot$D == '3D'] <- 'blue'

#####
m_C_vec <- as.matrix(10^seq(mLims_all[1],mLims_all[2],len = Resol))

CMassMat <- m_C_vec[,rep(1,Resol)] #Initialize mass matrices
RMassMat <- t(CMassMat)
kMat <- RMassMat/CMassMat

par_list <- scal_pars(CMassMat,RMassMat,z_0,r_0,beta_r,B_0_2D,B_0_3D,x_0_2D, x_0_3D, beta_x, h_02D, h_03D,beta_a, ForStr,a_02D,a_03D,p_v,p_d,gamma_2D,gamma_3D)

z <- par_list$z; r <- par_list$r; LossRate_2D <- par_list$LossRate_2D; LossRate_3D <- par_list$LossRate_3D; x_2D <- par_list$x_2D; x_3D <- par_list$x_3D; K_2D <- par_list$K_2D; K_3D <- par_list$K_3D; h_2D <- par_list$h_2D; h_3D <- par_list$h_3D; a_2D <- par_list$a_2D; a_3D <- par_list$a_3D

######### Energetic bounds ##########

CRateMat_2D <- RMassMat * a_2D * x_2D / (1 + h_2D * a_2D * x_2D) #~Search-constrained consumption rate 
CRateMat_3D <- RMassMat * a_3D * x_3D / (1 + h_3D * a_3D * x_3D) #~Search-constrained consumption rate 

NetRateMat_2D <- CRateMat_2D - LossRate_2D
NetRateMat_3D <- CRateMat_3D - LossRate_3D

NetRateVec_2D <- melt(CRateMat_2D - LossRate_2D)
NetRateVec_3D <- melt(CRateMat_3D - LossRate_3D)
colnames(NetRateVec_2D) <- c('m_C', 'm_R','Value')
colnames(NetRateVec_3D) <- c('m_C', 'm_R','Value')

NetRateVec_2D$Value[NetRateVec_2D$Value <=0] <- 0 #threshold consumption rates
NetRateVec_3D$Value[NetRateVec_3D$Value <=0] <- 0
NetRateVec_2D$Value[NetRateVec_2D$Value > 0] <- 1
NetRateVec_3D$Value[NetRateVec_3D$Value > 0] <- 1

NetRateVec_2D$m_C <- melt(log10(CMassMat))$value
NetRateVec_2D$m_R <- melt(log10(RMassMat))$value
NetRateVec_2D$D <- '2D'
NetRateModToPlot_2D <- NetRateVec_2D[NetRateVec_2D$Value>0,] 

NetRateVec_3D$m_C <- melt(log10(CMassMat))$value
NetRateVec_3D$m_R <- melt(log10(RMassMat))$value
NetRateVec_3D$D <- '3D'
NetRateModToPlot_3D <- NetRateVec_3D[NetRateVec_3D$Value>0,]

NetRateModToPlot <- rbind(NetRateModToPlot_2D,NetRateModToPlot_3D)

################### Coexistence Bounds ################### 

a_2D_ <- a_2D/CMassMat # get mass-specific search rate
a_3D_ <- a_3D/CMassMat 
CoexMat_2D <- K_2D - (z / (eff * a_2D_)) #coexistence condition
CoexMat_3D <- K_3D - (z / (eff * a_3D_))

CoexVec_2D <- melt(CoexMat_2D)
CoexVec_3D <- melt(CoexMat_3D)
colnames(CoexVec_2D) <- c('m_C', 'm_R','Value')
colnames(CoexVec_3D) <- c('m_C', 'm_R','Value')

CoexVec_2D$Value[CoexVec_2D$Value <=0] <- 0 #threshold
CoexVec_3D$Value[CoexVec_3D$Value <=0] <- 0
CoexVec_2D$Value[CoexVec_2D$Value > 0] <- 1
CoexVec_3D$Value[CoexVec_3D$Value > 0] <- 1

CoexVec_2D$m_C <- melt(log10(CMassMat))$value
CoexVec_2D$m_R <- melt(log10(RMassMat))$value
CoexVec_2D$D <- '2D'
CoexModToPlot_2D <- CoexVec_2D[CoexVec_2D$Value>0,]

CoexVec_3D$m_C <- melt(log10(CMassMat))$value
CoexVec_3D$m_R <- melt(log10(RMassMat))$value
CoexVec_3D$D <- '3D'
CoexModToPlot_3D <- CoexVec_3D[CoexVec_3D$Value>0,]

CoexModToPlot <- rbind(CoexModToPlot_2D,CoexModToPlot_3D)

#%%%%%%%%%%%%%%%%%%%%% Stability bounds %%%%%%%%%%%%%%%%%%%%% 

StabVec_2D <- as.data.frame(matrix(0,Resol^2,3)); StabVec_3D <- StabVec_2D
colnames(StabVec_2D) <- c('m_C', 'm_R','Value'); 
colnames(StabVec_3D) <- c('m_C', 'm_R','Value')
StabVec_2D$D <- '2D'; StabVec_3D$D <- '3D'
StabVec_2D$m_C <- melt(log10(CMassMat))$value
StabVec_2D$m_R <- melt(log10(RMassMat))$value
StabVec_3D$m_C <- melt(log10(CMassMat))$value
StabVec_3D$m_R <- melt(log10(RMassMat))$value

J11_2D <- melt((r*z*(K_2D*(h_2D^2)*z*a_2D_ + eff + h_2D*(z - K_2D*a_2D_*eff)))/(eff*K_2D*a_2D_*(h_2D*z - eff)))$value #1st elements of jacobian matrix 
J12_2D <- melt(-z/eff)$value # 2nd elements
J21_2D <- melt(r*(eff - z/(K_2D*a_2D_)- h_2D*z ))$value # 3rd elements
J22_2D <- melt(matrix(0,Resol,1))$value # 4th elements
JacMat_2D <- cbind(J11_2D,J12_2D,J21_2D,J22_2D) # combine into Resol x 4 matrix

J11_3D <- melt((r*z*(K_3D*(h_3D^2)*z*a_3D_ + eff + h_3D*(z - K_3D*a_3D_*eff)))/(eff*K_3D*a_3D_*(h_3D*z - eff)))$value #1st elements of jacobian matrix 
J12_3D <- melt(-z/eff)$value # 2nd elements
J21_3D <- melt(r*(eff - z/(K_3D*a_3D_)- h_3D*z ))$value # 3rd elements
J22_3D <- melt(matrix(0,Resol,1))$value # 4th elements
JacMat_3D <- cbind(J11_3D,J12_3D,J21_3D,J22_3D) # combine into Resol x 4 matrix

for (k in 1:Resol^2){
	tmpEig2D <- eigen(rbind(c(JacMat_2D[k,1:2]),c(JacMat_2D[k,3:4])), symmetric = FALSE,only.values = TRUE)
	tmpEig3D <- eigen(rbind(c(JacMat_3D[k,1:2]),c(JacMat_3D[k,3:4])), symmetric = FALSE,only.values = TRUE)
	if (max(Re(tmpEig2D$values)) < 0 && all(Im(tmpEig2D$values)== 0)) {StabVec_2D$Value[k] = 4} #Stable point
	else if (max(Re(tmpEig2D$values)) < 0 && any(Im(tmpEig2D$values)>0)) {StabVec_2D$Value[k] = 3} #Stable point with damped oscillations
	else if (all(Re(tmpEig2D$values) == 0) && any(Im(tmpEig2D$values)>0)) {StabVec_2D$Value[k] = 2} # neutrally stable with oscillations
	else if (any(Re(tmpEig2D$values)>0)	&& any(Im(tmpEig2D$values)>0)) {StabVec_2D$Value[k] = 1} # unstable with oscillations
	else if (any(Re(tmpEig2D$values)>0)	&& all(Im(tmpEig2D$values)==0)) {StabVec_2D$Value[k] = 0} # unstable
	else if (all(Re(tmpEig2D$values)==0) && all(Im(tmpEig2D$values)==0)) {StabMat_2D[k] = 0} #?
	
	if (max(Re(tmpEig3D$values)) < 0 && all(Im(tmpEig3D$values)== 0)) {StabVec_3D$Value[k] = 4} #Stable point
	else if (max(Re(tmpEig3D$values)) < 0 && any(Im(tmpEig3D$values)> 0)) {StabVec_3D$Value[k] = 3} #Stable point with oscillations
	else if (all(Re(tmpEig3D$values) == 0) && any(Im(tmpEig3D$values)>0)) {StabVec_3D$Value[k] = 2} #neutrally stable with oscillations
	else if (any(Re(tmpEig3D$values)> 0) && any(Im(tmpEig3D$values)>0)) {StabVec_3D$Value[k] = 1} #unstable with oscillations
	else if (any(Re(tmpEig3D$values)>0)	&& all(Im(tmpEig3D$values)==0)) {StabVec_3D$Value[k] = 0} #unstable
	else if (all(Re(tmpEig3D$values)==0) && all(Im(tmpEig3D$values)==0)) {StabVec_3D$Value[k] = 0} #?
	# print(k)
}

# StabModToPlot <- rbind(StabVec_2D,StabVec_3D)

StabModToPlot_2D <- StabVec_2D[StabVec_2D$Value > 0,]
StabModToPlot_3D <- StabVec_3D[StabVec_3D$Value > 0,]
StabModToPlot <- rbind(StabModToPlot_2D,StabModToPlot_3D)

##### Plotting #####   

ax_lims <- data.frame(x1=mLims_all[1],x2=mLims_all[2],y1=mLims_all[1],y2=mLims_all[2])
p <- ggplot() + theme_bw() +
  	geom_tile(data=NetRateModToPlot,aes(x = as.numeric(m_C),y = as.numeric(m_R), fill=Value)) + #this will yield large files
	scale_fill_gradient(low="white", high="green") +
	geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), linetype = "solid",colour = "black", size = .15, data = ax_lims) +		
	geom_point(data = DataToPlot, aes(x = m_C,y = m_R,colour = D), size = I(2), alpha = .60) +
	guides(fill=FALSE, colour = "none") + # turn off all marker and fill legends
	facet_wrap( ~ D, ncol = 2) + 
	theme(aspect.ratio = 1) +
	# xlim(mLims_all[1], mLims_all[2]) + ylim(mLims_all[1], mLims_all[2])+
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
	xlab(expression(log[10](italic(m[C])))) + 
	ylab(expression(log[10](italic(m[R]))))

ggsave(filename = paste(ODir, "mCmR_energetic_bounds_",ParaName,'.pdf', sep = ""), plot = p, height = 11.7, width = 8.3)

p1 <- ggplot() + theme_bw() +
  	geom_tile(data=CoexModToPlot,aes(x = as.numeric(m_C),y = as.numeric(m_R), fill=Value)) + #this will yield large files
	scale_fill_gradient(low="white", high="green") + 
	geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), linetype = "solid",colour = "black", size = .15, data = ax_lims) +	
	geom_point(data = DataToPlot, aes(x = m_C,y = m_R,colour = D), size = I(2), alpha = .60)+
	guides(fill=FALSE, colour = "none") + # turn off all marker and fill legends
	facet_wrap( ~ D, ncol = 2) + 
	theme(aspect.ratio = 1) +
	# xlim(mLims_all[1], mLims_all[2]) + ylim(mLims_all[1], mLims_all[2])+
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
	xlab(expression(log[10](italic(m[C])))) + 
	ylab(expression(log[10](italic(m[R])))) 

ggsave(filename = paste(ODir, "mCmR_coex_bounds_",ParaName,'.pdf', sep = ""), plot = p1, height = 11.7, width = 8.3)

p <- ggplot() + theme_bw() +
  	geom_tile(data=StabModToPlot,aes(x = as.numeric(m_C),y = as.numeric(m_R), fill=Value)) + #this will yield large files
	scale_fill_gradient(low="yellow", high="orange") + 
	geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), linetype = "solid",colour = "black", size = .15, data = ax_lims) +	
	geom_point(data = DataToPlot, aes(x = m_C,y = m_R,colour = D), size = I(1), alpha = .75) +
	guides(fill=FALSE, colour = "none") + # turn off all marker and fill legends
	facet_wrap( ~ D, ncol = 2) + 
	theme(aspect.ratio = 1) +
	# xlim(mLims_all[1], mLims_all[2]) + ylim(mLims_all[1], mLims_all[2])+
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
	xlab(expression(log[10](italic(m[C])))) + 
	ylab(expression(log[10](italic(m[R]))))

ggsave(filename = paste(ODir, "mCmR_stab_bounds_",ParaName,'.pdf', sep = ""), plot = p, height = 11.7, width = 8.3)
