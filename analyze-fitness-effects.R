#24 August 2017
#TNS - tyler.n.starr@gmail.com
#script to read in data from bulk growth assays and isogenic growths to derive and analyze selection coefficients of mutations between historical states in Hsp90

setwd("/path/to/source/file/directory") #also, make an empty folder in this location called "plots" to deposit all figures

library(mixtools)
library(seqinr)
library(ggplot2)
library(Biostrings)

################################################################################
#collating data for Sc-reversions experiment
#read in bulk growth data for individual reversion mutations to Sc: gives log(freq(mut)/freq(wt)) versus time-in-#-generations
Sc.rep1 <- read.csv(file="data_in/Sc-reversions_rep1.csv",header=T, stringsAsFactors=F)
Sc.rep2 <- read.csv(file="data_in/Sc-reversions_rep2.csv",header=T, stringsAsFactors=F)

#read in table of # bc sequence counts at each of the same timepoints for WT barcodes in the Sc reversions bulk competition
Sc.WT.rep1 <- read.csv(file="data_in/Sc-reversions_WT_rep1.csv",header=F,stringsAsFactors=F)
Sc.WT.rep2 <- read.csv(file="data_in/Sc-reversions_WT_rep2.csv",header=F,stringsAsFactors=F)
#reversion mutants were represented by a median of 763 reads in the first timepoint, 75% represented by >572 reads --> 
# #find pooled groups of wt barcodes to create independent samples that each start with ~500 reads; this will give a conservative SEM for wildtype, as most variants have higher sequencing depth than this

#rep1 groups: [1:26],[27:55],[56:87],[88:112]
Sc.WT.rep1.pooled <- data.frame(position=c("WT1","WT2","WT3","WT4"),codon=c("WT1","WT2","WT3","WT4"),aa=c("WT1","WT2","WT3","WT4"),gen0=rep(NA,4),gen2.8=rep(NA,4),gen4.2=rep(NA,4),gen6=rep(NA,4),gen9.3=rep(NA,4),gen13.6=rep(NA,4),gen15.7=rep(NA,4),gen19.3=rep(NA,4))
Sc.WT.rep1.pooled[1,4:11] <- colSums(Sc.WT.rep1[1:26,2:9])
Sc.WT.rep1.pooled[2,4:11] <- colSums(Sc.WT.rep1[27:55,2:9])
Sc.WT.rep1.pooled[3,4:11] <- colSums(Sc.WT.rep1[56:87,2:9])
Sc.WT.rep1.pooled[4,4:11] <- colSums(Sc.WT.rep1[88:112,2:9])

#rep2 groups: [1:18],[19:36],[37:58],[59:80],[81:93],[94:108]
Sc.WT.rep2.pooled <- data.frame(position=c("WT1","WT2","WT3","WT4","WT5","WT6"),codon=c("WT1","WT2","WT3","WT4","WT5","WT6"),aa=c("WT1","WT2","WT3","WT4","WT5","WT6"),gen0=rep(NA,6),gen1.2=rep(NA,6),gen2.8=rep(NA,6),gen4.2=rep(NA,6),gen6=rep(NA,6),gen9.3=rep(NA,6),gen11.4=rep(NA,6),gen13.6=rep(NA,6),gen15.7=rep(NA,6),gen19.3=rep(NA,6))
Sc.WT.rep2.pooled[1,4:13] <- colSums(Sc.WT.rep2[1:18,2:11])
Sc.WT.rep2.pooled[2,4:13] <- colSums(Sc.WT.rep2[19:36,2:11])
Sc.WT.rep2.pooled[3,4:13] <- colSums(Sc.WT.rep2[37:58,2:11])
Sc.WT.rep2.pooled[4,4:13] <- colSums(Sc.WT.rep2[59:80,2:11])
Sc.WT.rep2.pooled[5,4:13] <- colSums(Sc.WT.rep2[81:93,2:11])
Sc.WT.rep2.pooled[6,4:13] <- colSums(Sc.WT.rep2[94:108,2:11])

#convert read counts to log(reads_i/reads_WT), where reads_WT is the total number of WT reads at time x
Sc.WT.rep1.pooled$gen0 <- log(Sc.WT.rep1.pooled$gen0/2345)
Sc.WT.rep1.pooled$gen2.8 <- log(Sc.WT.rep1.pooled$gen2.8/3098)
Sc.WT.rep1.pooled$gen4.2 <- log(Sc.WT.rep1.pooled$gen4.2/4489)
Sc.WT.rep1.pooled$gen6 <- log(Sc.WT.rep1.pooled$gen6/4340)
Sc.WT.rep1.pooled$gen9.3 <- log(Sc.WT.rep1.pooled$gen9.3/4042)
Sc.WT.rep1.pooled$gen13.6 <- log(Sc.WT.rep1.pooled$gen13.6/4069)
Sc.WT.rep1.pooled$gen15.7 <- log(Sc.WT.rep1.pooled$gen15.7/905)
Sc.WT.rep1.pooled$gen19.3 <- log(Sc.WT.rep1.pooled$gen19.3/3570)

Sc.WT.rep2.pooled$gen0 <- log(Sc.WT.rep2.pooled$gen0/3069)
Sc.WT.rep2.pooled$gen1.2 <- log(Sc.WT.rep2.pooled$gen1.2/6284)
Sc.WT.rep2.pooled$gen2.8 <- log(Sc.WT.rep2.pooled$gen2.8/3901)
Sc.WT.rep2.pooled$gen4.2 <- log(Sc.WT.rep2.pooled$gen4.2/6847)
Sc.WT.rep2.pooled$gen6 <- log(Sc.WT.rep2.pooled$gen6/4671)
Sc.WT.rep2.pooled$gen9.3 <- log(Sc.WT.rep2.pooled$gen9.3/4844)
Sc.WT.rep2.pooled$gen11.4 <- log(Sc.WT.rep2.pooled$gen11.4/5860)
Sc.WT.rep2.pooled$gen13.6 <- log(Sc.WT.rep2.pooled$gen13.6/7218)
Sc.WT.rep2.pooled$gen15.7 <- log(Sc.WT.rep2.pooled$gen15.7/6793)
Sc.WT.rep2.pooled$gen19.3 <- log(Sc.WT.rep2.pooled$gen19.3/5854)

#append WT reads to Sc.rep2 data
Sc.rep1 <- rbind(Sc.rep1,Sc.WT.rep1.pooled)
Sc.rep2 <- rbind(Sc.rep2,Sc.WT.rep2.pooled)

#for each variant, calculate slope of log(f(mut)/f(wt)) ~ time-in-generations: this gives raw selection coefficient from competitive growth assay (equivalent to sT in Chevin paper, per-generation selection coefficient)
timepoints1 <- c(0,2.8,4.2,6,9.3,13.6,15.7,19.3)
timepoints2 <- c(0,1.2,2.8,4.2,6,9.3,11.4,13.6,15.7,19.3)
for(i in 1:nrow(Sc.rep1)){
  fit <- lm(as.numeric(Sc.rep1[i,4:11]) ~ timepoints1)
  Sc.rep1$slope[i] <- summary(fit)$coefficients[2,1]
  Sc.rep1$slope.se[i] <- summary(fit)$coefficients[2,2]
  Sc.rep1$mean.res[i] <- mean(abs(summary(fit)$residuals))
}
for(i in 1:nrow(Sc.rep2)){
  fit <- lm(as.numeric(Sc.rep2[i,4:13]) ~ timepoints2)
  Sc.rep2$slope[i] <- summary(fit)$coefficients[2,1]
  Sc.rep2$slope.se[i] <- summary(fit)$coefficients[2,2]
  Sc.rep2$mean.res[i] <- mean(abs(summary(fit)$residuals))
}

#linear transform values such that relative fitness of null allele is 0
#convert raw s (slope) to fitness given w = e^s
Sc.rep1$raw.fitness <- exp(Sc.rep1$slope)
Sc.rep2$raw.fitness <- exp(Sc.rep2$slope)
#stretch relative fitness space to be from 0 (null) to 1 (wt)
Sc.rep1$rel.fitness <- (Sc.rep1$raw.fitness - Sc.rep1$raw.fitness[2])/(Sc.rep1$raw.fitness[1]-Sc.rep1$raw.fitness[2])
Sc.rep2$rel.fitness <- (Sc.rep2$raw.fitness - Sc.rep2$raw.fitness[2])/(Sc.rep2$raw.fitness[1]-Sc.rep2$raw.fitness[2])
#convert back to selection coefficient s = ln(w)
Sc.rep1$s <- log(Sc.rep1$rel.fitness)
Sc.rep2$s <- log(Sc.rep2$rel.fitness)

#save wt values
Sc.WT <- rbind(Sc.rep1[Sc.rep1$position %in% c("WT1","WT2","WT3","WT4"),c("position","rel.fitness","s")],Sc.rep2[Sc.rep2$position %in% c("WT1","WT2","WT3","WT4","WT5","WT6"),c("position","rel.fitness","s")])

#reformat s values into data table with other stats on each individual ancestral reversion; calculate mean.s's per replicate and across all 2 or 4 observations
Sc.data <- read.csv(file="./data_in/Sc-reversions_summary.csv",header=T, stringsAsFactors=F)
for(i in 1:nrow(Sc.data)){
  subset.r1 <- Sc.rep1[Sc.rep1$position==Sc.data[i,"Sc.position"] & Sc.rep1$aa==as.character(Sc.data[i,"anc.AA"]),]
  subset.r2 <- Sc.rep2[Sc.rep2$position==Sc.data[i,"Sc.position"] & Sc.rep2$aa==as.character(Sc.data[i,"anc.AA"]),]
  if(nrow(subset.r1)==0){
    Sc.data[i,"rep1.s1"] <- NA; Sc.data[i,"rep1.s2"] <- NA
  }else if(nrow(subset.r1)==1){
    Sc.data[i,"rep1.s1"] <- as.numeric(subset.r1[1,"s"]); Sc.data[i,"rep1.s2"] <- NA
  }else if(nrow(subset.r1)==2){
    Sc.data[i,"rep1.s1"] <- as.numeric(subset.r1[1,"s"]); Sc.data[i,"rep1.s2"] <- as.numeric(subset.r1[2,"s"])
  }
  if(nrow(subset.r2)==0){
    Sc.data[i,"rep2.s1"] <- NA; Sc.data[i,"rep2.s2"] <- NA
  }else if(nrow(subset.r2)==1){
    Sc.data[i,"rep2.s1"] <- as.numeric(subset.r2[1,"s"]); Sc.data[i,"rep2.s2"] <- NA
  }else if(nrow(subset.r2)==2){
    Sc.data[i,"rep2.s1"] <- as.numeric(subset.r2[1,"s"]); Sc.data[i,"rep2.s2"] <- as.numeric(subset.r2[2,"s"])
  }
  Sc.data[i,"rep1.mean.s"] <- mean(c(Sc.data[i,"rep1.s1"],Sc.data[i,"rep1.s2"]),na.rm=T)
  Sc.data[i,"rep2.mean.s"] <- mean(c(Sc.data[i,"rep2.s1"],Sc.data[i,"rep2.s2"]),na.rm=T)
  Sc.data[i,"mean.s"] <- mean(c(Sc.data[i,"rep1.s1"],Sc.data[i,"rep1.s2"],Sc.data[i,"rep2.s1"],Sc.data[i,"rep2.s2"]),na.rm=T)
  Sc.data[i,"SE.mean.s"] <- sd(c(Sc.data[i,"rep1.s1"],Sc.data[i,"rep1.s2"],Sc.data[i,"rep2.s1"],Sc.data[i,"rep2.s2"]),na.rm=T)/sqrt(sum(!is.na(c(Sc.data[i,"rep1.s1"],Sc.data[i,"rep1.s2"],Sc.data[i,"rep2.s1"],Sc.data[i,"rep2.s2"]))))
}

#indicator if measurement was made in parallel bulk assays (versus isogenic stuff added in below)
Sc.data$bulk <- FALSE
Sc.data[!is.na(Sc.data$mean.s),"bulk"] <- TRUE

#one mutation (T5S) was missed in library prep, growth rate was determined in isogenic culture in replicate
#for isogenic growths, read in ln(OD) versus time (in hrs) measures
Sc.iso.rep1 <- read.csv(file="data_in/Sc-reversions_isogenic-growths_rep1.csv",header=T, stringsAsFactors=F)
Sc.iso.rep2 <- read.csv(file="data_in/Sc-reversions_isogenic-growths_rep2.csv",header=T, stringsAsFactors=F)
Sc.iso.rep3 <- read.csv(file="data_in/Sc-reversions_isogenic-growths_rep3.csv",header=T, stringsAsFactors=F)
timepoints1 <- c(4,5.5,8.5,9.5,13.5)
timepoints2 <- c(4.5,6.5,8.5,12,14.5)
timepoints3 <- c(3,6,9,12,15.5)

for(i in 1:nrow(Sc.iso.rep1)){
  fit <- lm(as.numeric(Sc.iso.rep1[i,2:6]) ~ timepoints1)
  Sc.iso.rep1$slope[i] <- summary(fit)$coefficients[2,1] #slope of ln(OD) versus time is malthusian parameter
  Sc.iso.rep1$slope.se[i] <- summary(fit)$coefficients[2,2]
}
#calculate dimensionless per-generation selection coefficient, assuming that seleciton is density-independent, as (slope_x - slope_wt)/slope_wt*ln2 (last part equivalent to multiplying by generation time of WT)
Sc.iso.rep1$s <- (Sc.iso.rep1$slope-Sc.iso.rep1$slope[1])/(Sc.iso.rep1$slope[1])*log(2)
#linear transform values such that relative fitness of null allele is 0
#convert raw s to fitness given w=e^s
Sc.iso.rep1$raw.fitness <- exp(Sc.iso.rep1$s)
#stretch relative fitness space to be from 0 (null) to 1 (wt)
Sc.iso.rep1$rel.fitness <- (Sc.iso.rep1$raw.fitness - Sc.iso.rep1$raw.fitness[2])/(Sc.iso.rep1$raw.fitness[1] - Sc.iso.rep1$raw.fitness[2])
#convert back to selection coefficient s=ln(w)
Sc.iso.rep1$s <- log(Sc.iso.rep1$rel.fitness)

for(i in 1:nrow(Sc.iso.rep2)){
  fit <- lm(as.numeric(Sc.iso.rep2[i,2:6]) ~ timepoints2)
  Sc.iso.rep2$slope[i] <- summary(fit)$coefficients[2,1] #slope of ln(OD) versus time is malthusian parameter
  Sc.iso.rep2$slope.se[i] <- summary(fit)$coefficients[2,2]
}
#calculate dimensionless per-generation selection coefficient, assuming that seleciton is density-independent, as (slope_x - slope_wt)/slope_wt*ln2 (last part equivalent to multiplying by generation time of WT)
Sc.iso.rep2$s <- (Sc.iso.rep2$slope-Sc.iso.rep2$slope[1])/(Sc.iso.rep2$slope[1])*log(2)
#linear transform values such that relative fitness of null allele is 0
#convert raw s to fitness given w=e^s
Sc.iso.rep2$raw.fitness <- exp(Sc.iso.rep2$s)
#stretch relative fitness space to be from 0 (null) to 1 (wt)
Sc.iso.rep2$rel.fitness <- (Sc.iso.rep2$raw.fitness - Sc.iso.rep2$raw.fitness[2])/(Sc.iso.rep2$raw.fitness[1] - Sc.iso.rep2$raw.fitness[2])
#convert back to selection coefficient s=ln(w)
Sc.iso.rep2$s <- log(Sc.iso.rep2$rel.fitness)

for(i in 1:nrow(Sc.iso.rep3)){
  fit <- lm(as.numeric(Sc.iso.rep3[i,2:6]) ~ timepoints3)
  Sc.iso.rep3$slope[i] <- summary(fit)$coefficients[2,1] #slope of ln(OD) versus time is malthusian parameter
  Sc.iso.rep3$slope.se[i] <- summary(fit)$coefficients[2,2]
}
#calculate dimensionless per-generation selection coefficient, assuming that seleciton is density-independent, as (slope_x - slope_wt)/slope_wt*ln2 (last part equivalent to multiplying by generation time of WT)
Sc.iso.rep3$s <- (Sc.iso.rep3$slope-Sc.iso.rep3$slope[1])/(Sc.iso.rep3$slope[1])*log(2)
#linear transform values such taht relative fitness of null allele is 0
#convert raw s to fitness given w=e^s
Sc.iso.rep3$raw.fitness <- exp(Sc.iso.rep3$s)
#stretch relative fitness space to be from 0 (null) to 1 (wt)
Sc.iso.rep3$rel.fitness <- (Sc.iso.rep3$raw.fitness - Sc.iso.rep3$raw.fitness[2])/(Sc.iso.rep3$raw.fitness[1] - Sc.iso.rep3$raw.fitness[2])
#convert back to selection coefficient s=ln(w)
Sc.iso.rep3$s <- log(Sc.iso.rep3$rel.fitness)

for(mut in Sc.data[is.na(Sc.data$mean.s),"shorthand"]){
  Sc.data[Sc.data$shorthand==mut,"mean.s"] <- mean(Sc.iso.rep1[Sc.iso.rep1$ID==mut,"s"],Sc.iso.rep1[Sc.iso.rep1$ID==mut,"s"])
}

################################################################################################################################

#collating data for ancAmo-fwds experiment
#read in bulk growth data for individual fwd sub mutations to ancAmo: gives log(freq(mut)/freq(wt)) versus time in #-generations, where wt is ancAmo+L378i
ancAmo.rep1 <- read.csv(file="data_in/ancAmo-fwd-subs_rep1.csv",header=T, stringsAsFactors=F)
ancAmo.rep2 <- read.csv(file="data_in/ancAmo-fwd-subs_rep2.csv",header=T, stringsAsFactors=F)

#read in table of # bc sequence counts versus time (in hrs) for WT barcodes in the Sc reversions bulk competition, rep2
ancAmo.WT.rep1 <- read.csv(file="data_in/ancAmo-fwd-subs_WT_rep1.csv",header=F,stringsAsFactors=F)
ancAmo.WT.rep2 <- read.csv(file="data_in/ancAmo-fwd-subs_WT_rep2.csv",header=F,stringsAsFactors=F)
#forward mutants were represented by a median of 29,000 reads in first timepoint, 75% represented by >17500 reads
#find pooled groups of wt barcodes to create independent samples that each start with ~17000 reads
#rep1 groups: [c(8,9)](17343),[c(2,7)](16861),[c(3,5,6)](17119),[c(1,4)](18165)
ancAmo.WT.rep1.pooled <- data.frame(position=c("WT1","WT2","WT3","WT4"),codon=c("WT1","WT2","WT3","WT4"),aa=c("WT1","WT2","WT3","WT4"),gen0=rep(NA,4),gen1.25=rep(NA,4),gen2.5=rep(NA,4),gen3.75=rep(NA,4),gen6.67=rep(NA,4),gen9.17=rep(NA,4),gen10.83=rep(NA,4),gen12.92=rep(NA,4))
ancAmo.WT.rep1.pooled[1,4:11] <- colSums(ancAmo.WT.rep1[c(8,9),2:9])
ancAmo.WT.rep1.pooled[2,4:11] <- colSums(ancAmo.WT.rep1[c(2,7),2:9])
ancAmo.WT.rep1.pooled[3,4:11] <- colSums(ancAmo.WT.rep1[c(3,5,6),2:9])
ancAmo.WT.rep1.pooled[4,4:11] <- colSums(ancAmo.WT.rep1[c(1,4),2:9])

#rep2 groups: [c(1,3,7)](14380),[c(6,8)](14560),[c(2,4,5,9)](14444)
ancAmo.WT.rep2.pooled <- data.frame(position=c("WT1","WT2","WT3"),codon=c("WT1","WT2","WT3"),aa=c("WT1","WT2","WT3"),gen0=rep(NA,3),gen1.25=rep(NA,3),gen2.5=rep(NA,3),gen3.75=rep(NA,3),gen6.67=rep(NA,3),gen9.17=rep(NA,3),gen10.83=rep(NA,3),gen12.92=rep(NA,3))
ancAmo.WT.rep2.pooled[1,4:11] <- colSums(ancAmo.WT.rep2[c(1,3,7),2:9])
ancAmo.WT.rep2.pooled[2,4:11] <- colSums(ancAmo.WT.rep2[c(6,8),2:9])
ancAmo.WT.rep2.pooled[3,4:11] <- colSums(ancAmo.WT.rep2[c(2,4,5,9),2:9])

#convert read counts to log(reads_i/reads_WT), where reads_WT is the total number of WT reads at time gen.x
ancAmo.WT.rep1.pooled$gen0 <- log(ancAmo.WT.rep1.pooled$gen0/69488)
ancAmo.WT.rep1.pooled$gen1.25 <- log(ancAmo.WT.rep1.pooled$gen1.25/107075)
ancAmo.WT.rep1.pooled$gen2.5 <- log(ancAmo.WT.rep1.pooled$gen2.5/102544)
ancAmo.WT.rep1.pooled$gen3.75 <- log(ancAmo.WT.rep1.pooled$gen3.75/85927)
ancAmo.WT.rep1.pooled$gen6.67 <- log(ancAmo.WT.rep1.pooled$gen6.67/112269)
ancAmo.WT.rep1.pooled$gen9.17 <- log(ancAmo.WT.rep1.pooled$gen9.17/108183)
ancAmo.WT.rep1.pooled$gen10.83 <- log(ancAmo.WT.rep1.pooled$gen10.83/91876)
ancAmo.WT.rep1.pooled$gen12.92 <- log(ancAmo.WT.rep1.pooled$gen12.92/116837)

ancAmo.WT.rep2.pooled$gen0 <- log(ancAmo.WT.rep2.pooled$gen0/43384)
ancAmo.WT.rep2.pooled$gen1.25 <- log(ancAmo.WT.rep2.pooled$gen1.25/48057)
ancAmo.WT.rep2.pooled$gen2.5 <- log(ancAmo.WT.rep2.pooled$gen2.5/47685)
ancAmo.WT.rep2.pooled$gen3.75 <- log(ancAmo.WT.rep2.pooled$gen3.75/44944)
ancAmo.WT.rep2.pooled$gen6.67 <- log(ancAmo.WT.rep2.pooled$gen6.67/52136)
ancAmo.WT.rep2.pooled$gen9.17 <- log(ancAmo.WT.rep2.pooled$gen9.17/51360)
ancAmo.WT.rep2.pooled$gen10.83 <- log(ancAmo.WT.rep2.pooled$gen10.83/55415)
ancAmo.WT.rep2.pooled$gen12.92 <- log(ancAmo.WT.rep2.pooled$gen12.92/52959)


#append WT reads to ancAmo.repx data
ancAmo.rep1 <- rbind(ancAmo.rep1,ancAmo.WT.rep1.pooled)
ancAmo.rep2 <- rbind(ancAmo.rep2,ancAmo.WT.rep2.pooled)

#for each variant, calculate slope of log(f(mut)/f(wt)) ~ time-in-generations: this gives raw selection coefficient from competitive growth assay
timepoints1 <- c(0,1.25,2.5,3.75,6.67,9.17,10.83,12.92)
timepoints2 <- c(0,1.25,2.5,3.75,6.67,9.17,10.83,12.92)
for(i in 1:nrow(ancAmo.rep1)){
  fit <- lm(as.numeric(ancAmo.rep1[i,4:11]) ~ timepoints1)
  ancAmo.rep1$slope[i] <- summary(fit)$coefficients[2,1]
  ancAmo.rep1$slope.se[i] <- summary(fit)$coefficients[2,2]
  ancAmo.rep1$mean.res[i] <- mean(abs(summary(fit)$residuals))
}
for(i in 1:nrow(ancAmo.rep2)){
  fit <- lm(as.numeric(ancAmo.rep2[i,4:11]) ~ timepoints2)
  ancAmo.rep2$slope[i] <- summary(fit)$coefficients[2,1]
  ancAmo.rep2$slope.se[i] <- summary(fit)$coefficients[2,2]
  ancAmo.rep2$mean.res[i] <- mean(abs(summary(fit)$residuals))
}

#linear transform values such that relative fitness of null allele is 0
#convert raw s (slope) to fitness given w = e^s
ancAmo.rep1$raw.fitness <- exp(ancAmo.rep1$slope)
ancAmo.rep2$raw.fitness <- exp(ancAmo.rep2$slope)
#stretch relative fitness space to be from 0 (null) to 1 (wt ancAmo)
ancAmo.rep1$rel.fitness <- (ancAmo.rep1$raw.fitness - ancAmo.rep1$raw.fitness[2])/(ancAmo.rep1$raw.fitness[1]-ancAmo.rep1$raw.fitness[2])
ancAmo.rep2$rel.fitness <- (ancAmo.rep2$raw.fitness - ancAmo.rep2$raw.fitness[2])/(ancAmo.rep2$raw.fitness[1]-ancAmo.rep2$raw.fitness[2])
#convert back to selection coefficient s = ln(w)
ancAmo.rep1$s <- log(ancAmo.rep1$rel.fitness)
ancAmo.rep2$s <- log(ancAmo.rep2$rel.fitness)

#save wt values
ancAmo.WT <- rbind(ancAmo.rep1[ancAmo.rep1$position %in% c("WT1","WT2","WT3","WT4"),c("position","rel.fitness","s")],ancAmo.rep2[ancAmo.rep2$position %in% c("WT1","WT2","WT3"),c("position","rel.fitness","s")])

#fitness of ancAmoNTD and ancAmo+L378i (ancAmo-wt) relative to w.Sc=1, and Sc fitness relative to w.ancAmo-wt=1
ancAmo.w1 <- (ancAmo.rep1$rel.fitness[5]-ancAmo.rep1$rel.fitness[2])/(ancAmo.rep1$rel.fitness[3]-ancAmo.rep1$rel.fitness[2])
ancAmo.w2 <- (ancAmo.rep2$rel.fitness[5]-ancAmo.rep2$rel.fitness[2])/(ancAmo.rep2$rel.fitness[3]-ancAmo.rep2$rel.fitness[2])
ancAmo.mean.w <- mean(c(ancAmo.w1,ancAmo.w2))
ancAmo.se.w <- sd(c(ancAmo.w1,ancAmo.w2))/sqrt(2)
rm(ancAmo.w1);rm(ancAmo.w2)

ancAmo.i378.w1 <- (ancAmo.rep1$rel.fitness[1]-ancAmo.rep1$rel.fitness[2])/(ancAmo.rep1$rel.fitness[3]-ancAmo.rep1$rel.fitness[2])
ancAmo.i378.w2 <- (ancAmo.rep2$rel.fitness[1]-ancAmo.rep2$rel.fitness[2])/(ancAmo.rep2$rel.fitness[3]-ancAmo.rep2$rel.fitness[2])
ancAmo.i378.mean.w <- mean(c(ancAmo.i378.w1,ancAmo.i378.w2))
ancAmo.i378.se.w <- sd(c(ancAmo.i378.w1,ancAmo.i378.w2))/sqrt(2)
rm(ancAmo.i378.w1);rm(ancAmo.i378.w2)

Sc.mean.w <- mean(c(ancAmo.rep1[3,"rel.fitness"],ancAmo.rep2[3,"rel.fitness"]))
Sc.se.w <- sd(c(ancAmo.rep1[3,"rel.fitness"],ancAmo.rep2[3,"rel.fitness"]))/sqrt(2)

#reformat s values into data table with other stats on each individual ancestral reversion; calculate mean.s's per replicate and across all 2 or 4 observations (for ancAmo, all 2, once per rep)
ancAmo.data <- read.csv(file="./data_in/ancAmo-fwd-subs_summary.csv",header=T, stringsAsFactors=F)
for(i in 1:nrow(ancAmo.data)){
  subset.r1 <- ancAmo.rep1[ancAmo.rep1$position==ancAmo.data[i,"Sc.position"] & ancAmo.rep1$aa==as.character(ancAmo.data[i,"der.AA"]),]
  subset.r2 <- ancAmo.rep2[ancAmo.rep2$position==ancAmo.data[i,"Sc.position"] & ancAmo.rep2$aa==as.character(ancAmo.data[i,"der.AA"]),]
  if(nrow(subset.r1)==0){
    ancAmo.data[i,"rep1.s"] <- NA
  }else if(nrow(subset.r1)==1){
    ancAmo.data[i,"rep1.s"] <- as.numeric(subset.r1[1,"s"])
  }
  if(nrow(subset.r2)==0){
    ancAmo.data[i,"rep1.s"] <- NA
  }else if(nrow(subset.r2)==1){
    ancAmo.data[i,"rep2.s"] <- as.numeric(subset.r2[1,"s"])
  }
  ancAmo.data[i,"mean.s"] <- mean(c(ancAmo.data[i,"rep1.s"],ancAmo.data[i,"rep2.s"]),na.rm=T)
  ancAmo.data[i,"SE.mean.s"] <- sd(c(ancAmo.data[i,"rep1.s"],ancAmo.data[i,"rep2.s"]),na.rm=T)/sqrt(sum(!is.na(c(ancAmo.data[i,"rep1.s"],ancAmo.data[i,"rep2.s"]))))
}

#indicator if measurement was made in parallel bulk assays (versus isogenic stuff added in below)
ancAmo.data$bulk <- FALSE
ancAmo.data[!is.na(ancAmo.data$mean.s),"bulk"] <- TRUE

#several mutations were missed in library prep, growth rate was determined in isogenic culture
#for isogenic growths, read in ln(OD) versus time (in hrs) measures
ancAmo.iso.rep1 <- read.csv(file="data_in/ancAmo-fwd-subs_isogenic-growths_rep1.csv",header=T, stringsAsFactors=F)

timepoints1 <- c(21,23,25,30)

for(i in 1:nrow(ancAmo.iso.rep1)){
  fit <- lm(as.numeric(ancAmo.iso.rep1[i,2:5]) ~ timepoints1)
  ancAmo.iso.rep1$slope[i] <- summary(fit)$coefficients[2,1] #slope of ln(OD) versus time is Malthusian parameter
  ancAmo.iso.rep1$slope.se[i] <- summary(fit)$coefficients[2,2]
}
#calculate dimensionless per-generation selection coefficient, assuming that seleciton is density-independent, as (slope_x - slope_wt)/slope_wt*ln2 (last part equivalent to multiplying by generation time of WT)
ancAmo.iso.rep1$s <- (ancAmo.iso.rep1$slope-ancAmo.iso.rep1$slope[1])/(ancAmo.iso.rep1$slope[1])*log(2)
#linear transform values such that relative fitness of null allele is 0
#convert raw s to fitness given w=e^s
ancAmo.iso.rep1$raw.fitness <- exp(ancAmo.iso.rep1$s)
#stretch relative fitness space to be from 0 (null) to 1 (wt)
ancAmo.iso.rep1$rel.fitness <- (ancAmo.iso.rep1$raw.fitness - ancAmo.iso.rep1$raw.fitness[2])/(ancAmo.iso.rep1$raw.fitness[1] - ancAmo.iso.rep1$raw.fitness[2])
#convert back to selection coefficient s=ln(w)
ancAmo.iso.rep1$s <- log(ancAmo.iso.rep1$rel.fitness)

for(mut in ancAmo.data[is.na(ancAmo.data$mean.s),"shorthand"]){
  ancAmo.data[ancAmo.data$shorthand==mut,"mean.s"] <- ancAmo.iso.rep1[ancAmo.iso.rep1$ID==mut,"s"]
}

################################################################################################################################
#plot: mean.s between replicates for each bulk experiment
pdf(file="plots/replicate-correlations.pdf",height=5,width=9,useDingbats=F)
par(mfrow=c(1,2))
#ancAmo forwards:
plot(ancAmo.data$rep1.s[ancAmo.data$bulk==T],ancAmo.data$rep2.s[ancAmo.data$bulk==T],ylim=c(-0.6,0.03),xlim=c(-0.6,0.03),main="AncAmo forward mutations",xlab="mean s, replicate 1",ylab="mean s, replicate 2",pch=20); abline(lm(ancAmo.data$rep2.s[ancAmo.data$bulk==T]~ancAmo.data$rep1.s[ancAmo.data$bulk==T])); summary(lm(ancAmo.data$rep2.s[ancAmo.data$bulk==T]~ancAmo.data$rep1.s[ancAmo.data$bulk==T]))
text(-0.6,0.0,expression(paste(R^2,"=0.943")),adj=c(0,0))
text(-0.6,-0.04,"slope=1.06",adj=c(0,0))
#Sc reversions
plot(Sc.data$rep1.mean.s[Sc.data$bulk==T],Sc.data$rep2.mean.s[Sc.data$bulk==T],ylim=c(-0.6,0.03),xlim=c(-0.6,0.03),main="Sc reverse mutations",xlab="mean s, replicate 1",ylab="mean s, replicate 2",pch=20); abline(lm(Sc.data$rep2.mean.s[Sc.data$bulk==T]~Sc.data$rep1.mean.s[Sc.data$bulk==T])); summary(lm(Sc.data$rep2.mean.s[Sc.data$bulk==T]~Sc.data$rep1.mean.s[Sc.data$bulk==T]))
text(-0.6,0.0,expression(paste(R^2,"=0.984")),adj=c(0,0))
text(-0.6,-0.04,"slope=0.91",adj=c(0,0))
dev.off()

#zoom in on nearly neutral region for each plot, for inset
pdf(file="plots/replicate-correlations_zoom.pdf",height=5,width=9,useDingbats=F)
par(mfrow=c(1,2))
#ancAmo forwards:
plot(ancAmo.data$rep1.s[ancAmo.data$bulk==T],ancAmo.data$rep2.s[ancAmo.data$bulk==T],ylim=c(-0.08,0.03),xlim=c(-0.08,0.03),main="AncAmo forward mutations",xlab="mean s, replicate 1",ylab="mean s, replicate 2",pch=20); abline(lm(ancAmo.data$rep2.s[ancAmo.data$bulk==T]~ancAmo.data$rep1.s[ancAmo.data$bulk==T]))
abline(v=0,lty=2)
abline(h=0,lty=2)
#Sc reversions
plot(Sc.data$rep1.mean.s[Sc.data$bulk==T],Sc.data$rep2.mean.s[Sc.data$bulk==T],ylim=c(-0.08,0.03),xlim=c(-0.08,0.03),main="Sc reverse mutations",xlab="mean s, replicate 1",ylab="mean s, replicate 2",pch=20); abline(lm(Sc.data$rep2.mean.s[Sc.data$bulk==T]~Sc.data$rep1.mean.s[Sc.data$bulk==T]))
abline(v=0,lty=2)
abline(h=0,lty=2)
dev.off()

################################################################################################################################
#average reversion, forward mutation is deleterious
#bulk determination only
#ancAmo, bulk competition values
wilcox.test(ancAmo.data$mean.s[ancAmo.data$bulk==TRUE],conf.int=T)
#ancAmo, bulk plus filled in missing values from isogenic growths
wilcox.test(ancAmo.data$mean.s, conf.int=T)
#ancAmo, excluding strongly deleterious outliers
wilcox.test(ancAmo.data$mean.s[ancAmo.data$bulk==TRUE & ancAmo.data$mean.s>-0.1],conf.int=T)
#ancAmo, excluding derived states that are not unambiguously reconstructed
wilcox.test(ancAmo.data[ancAmo.data$max.PP ==1 & ancAmo.data$bulk==TRUE,"mean.s"])
#ancAmo, excluding states that are not unambiguously reconstructed ancestrally in ancAmo
wilcox.test(ancAmo.data[ancAmo.data$ancAmo.PP==1 & ancAmo.data$bulk==TRUE,"mean.s"])
#ancAmo, excluding states whose ancestral or derved PP is not 1.0
wilcox.test(ancAmo.data[ancAmo.data$max.PP ==1 & ancAmo.data$ancAmo.PP==1,"mean.s"],conf.int=T)

#Sc, bulk competition values
wilcox.test(Sc.data$mean.s[Sc.data$bulk==TRUE],conf.int=T)#$p.value
#Sc, bulk plus filled in missing values from isogenic growths
wilcox.test(Sc.data$mean.s,conf.int=T)#$p.value
#Sc, excluding strongly deleterious outliers
wilcox.test(Sc.data$mean.s[Sc.data$bulk==TRUE & Sc.data$mean.s > -0.1],conf.int=T)#$p.value
#Sc, exclude ancestral states that don't achieve posterior probability of 1.0 in any ancestor of the trajectory
wilcox.test(Sc.data[Sc.data$max.PP ==1 ,"mean.s"],conf.int=T)

#relating the average s of reversions as measured in this experiment to an estimate of the s as measured at wt level Hsp90 expression
#from Jiang et al, have some relevant information:
#the strain in which experiments were performed has Hsp90 expression levels 0.008 (0.8%) that of "Wt" Scer (approximated by the GPD promoter); Julia says might be more like 5% with updated data, but this is the minimum so the most conservative extreme estimate
#function between relative growth rate (wt ~ 1) and amount of Hsp90 activity (expression) is y = (x)/(x+0.014), from Jiang et al. PLOS Genetics 2013
s_exp <- c(-0.01129891, -0.00964616 , -0.008057479) #95% CI and median s of reversions, excluding strongly deleterious muts
sGR_exp <- s_exp * log(2) #convert per generation selection coefficient to relative growth rate difference metric used in microbial growth curve comparisons
x_ref <- 0.008 #relative expression of Hsp90 in our experimental reference strain is 0.008 that of GPD (~wt)
y_ref <- (x_ref)/(x_ref+0.014) #relative growth rate of the strain experiments were conducted in
y_exp <- y_ref + sGR_exp #express growth rate of mean reversion relative to GPD strain instead of experimental strain
x_exp <- (0.014*y_exp)/(1-y_exp) #relative Hsp90 activity of the mean reversion
delta_x_exp <- x_exp - 0.008 #difference in Hsp90 functional activity of the average reversion --> this is scale-independent metric of effect of the mutations!
#what is the relative growth rate difference associated with this same difference in Hsp90 function at the GPD (~wt) expression level?
x_gpd <- 1 #defined as 1
y_gpd <- (x_gpd)/(x_gpd + 0.014) #corresponding relative growth rate (should be 1, but not quite because of model fit)
y_exp_gpd <- (1-delta_x_exp)/((1-delta_x_exp)+0.014) #what the relative growth rate would be for an identical functional defect of reversion but at GPD expression levels
sGR_gpd <- y_exp_gpd - y_gpd #what the relative growth rate difference would be versus GPD
s_gpd <- sGR_gpd/log(2) #what the per-generation selection coefficient would be versus GPD-level expressing wildtype strain, given the observed functional defect of reversions

################################################################################################################################
#calculate the fraction of reversions, forward mutations that are deleterious via comparison to the wt sampling distribution
#Sc data
vals <- Sc.data$mean.s[Sc.data$bulk==TRUE & Sc.data$mean.s > -0.04] #mixtools has troubles with outliers -- infer mixture from the bulk part of the distribution, any below 0.04 are unambiguously deleterious anyway
hist(vals,breaks=20,col="gray50",xlab="Selection coefficient (relative to ScHsp90)",main="",freq=T,xlim=c(-0.03, 0.01));abline(v=0,lty=2)
hist(Sc.WT$s,add=T,col="#00FFFF",freq=T,breaks=8)
# #for summary figure, after inferring model can run following commented lines to illustrate mixture over distributions
# curve((mix3$lambda[1]*dnorm(x,mix3$mu[1],mix3$sigma[1]))*0.25,add=T,lty=2,col="darkred",lwd=2)
# curve((mix3$lambda[2]*dnorm(x,mix3$mu[2],mix3$sigma[2]))*0.25,add=T,lty=2,col="darkgreen",lwd=2)
# curve((mix3$lambda[3]*dnorm(x,mix3$mu[3],mix3$sigma[3]))*0.25,add=T,lty=2,col="darkblue",lwd=2)
# curve(dmix(x)*0.25,add=T,lwd=2)

#use mixtools to fit gaussian mixture models, one component of which is defined by the wt sampling distributions, others vary. Figure out best number of components via AIC
n <- 2; mix2 <- normalmixEM(vals, k=n, lambda=rep(1/n, n), mu=c(mean(Sc.WT$s),-0.01), sigma=rep(sd(Sc.WT$s),n), mean.constr=c(mean(Sc.WT$s), rep(NA,n-1)), sd.constr=c(sd(Sc.WT$s),rep(NA,n-1)), arbmean=T, arbvar=T)
n <- 3; mix3 <- normalmixEM(vals, k=n, lambda=rep(1/n, n), mu=c(mean(Sc.WT$s),-0.01,-0.005), sigma=rep(sd(Sc.WT$s),n), mean.constr=c(mean(Sc.WT$s), rep(NA,n-1)), sd.constr=c(sd(Sc.WT$s),rep(NA,n-1)), arbmean=T, arbvar=T)
n <- 4; mix4 <- normalmixEM(vals, k=n, lambda=rep(1/n, n), mu=c(mean(Sc.WT$s),-0.01,-0.005,0.005), sigma=rep(sd(Sc.WT$s),n), mean.constr=c(mean(Sc.WT$s), rep(NA,n-1)), sd.constr=c(sd(Sc.WT$s),rep(NA,n-1)), arbmean=T, arbvar=T)
n <- 5; mix5 <- normalmixEM(vals, k=n, lambda=rep(1/n, n), mu=c(mean(Sc.WT$s),-0.01,-0.005,-0.003,0.005), sigma=rep(sd(Sc.WT$s),n), mean.constr=c(mean(Sc.WT$s), rep(NA,n-1)), sd.constr=c(sd(Sc.WT$s),rep(NA,n-1)), arbmean=T, arbvar=T)
summary.mixEM(mix2) #loglik 334.3584, k = 3, AIC = -662.7168
summary.mixEM(mix3) #loglik 340.9462, k = 6, AIC = -669.8924
summary.mixEM(mix4) #loglik 341.8108, k = 9, AIC = -665.6216
summary.mixEM(mix5) #loglik 343.1684, k = 12, AIC = -662.3368

#3-component mixture favored by AIC
#define a density function for this 3-component mixture model
dmix <- function(x){
  return(mix3$lambda[1]*dnorm(x, mix3$mu[1], mix3$sigma[1])+mix3$lambda[2]*dnorm(x, mix3$mu[2], mix3$sigma[2])+mix3$lambda[3]*dnorm(x, mix3$mu[3], mix3$sigma[3]))
}

#qqplot to judge model fit
set.seed(10101)
qqplot(c(rnorm(round(100000*mix3$lambda[1]),mix3$mu[1],mix3$sigma[1]),rnorm(round(100000*mix3$lambda[2]),mix3$mu[2],mix3$sigma[2]),rnorm(round(100000*mix3$lambda[3]),mix3$mu[3],mix3$sigma[3])),vals,xlab="Theoretical quantiles (mixture)", ylab="Sample quantiles",xlim=c(-0.045, 0.025),ylim=c(-0.045,0.025))
abline(a=0,b=1)

#define the posterior probability of being deleterious, neutral, beneficial for each mutation -- comes from the relative pdf of deleterious versus wt (versus beneficial, though there is no beneficial component in this mixture model) at the given s for each mutation
Sc.data$PP.neut <- NA
Sc.data$PP.ben <- NA
Sc.data$PP.del <- NA
Sc.data[Sc.data$mean.s < -0.04,"PP.neut"] <- 0
Sc.data[Sc.data$mean.s > -0.04 & Sc.data$bulk==T,"PP.neut"] <- mix3$posterior[,1]
#no inferred beneficial component in the mixture model, so PP.ben is 0 and PP.del = 1-PP.neut for all variants
Sc.data$PP.del <- 1-Sc.data$PP.neut
Sc.data[Sc.data$bulk==T,"PP.ben"] <- 0

sum(Sc.data$PP.del,na.rm=T)/sum(!is.na(Sc.data$PP.del)) # 0.9310345
sum(Sc.data$PP.neut,na.rm=T)/sum(!is.na(Sc.data$PP.neut)) # 0.06896547
sum(Sc.data$PP.ben,na.rm=T)/sum(!is.na(Sc.data$PP.ben)) # 0

#plot of distribution of selection coefficients of reversions with an indication of the expected fraction in each bin that are neutral
pdf(file="plots/Sc-rev_mean-s-histogram_bulk_indicate-neutral-bg.pdf",width=5.5,height=5,useDingbats=F)
hist(Sc.data$mean.s[Sc.data$bulk==TRUE],breaks=seq(-0.6,0.01,0.01),col="gray50",xlab="selection coefficient of reverse mutation",main="",xlim=c(-0.6,0.05),ylim=c(0,55));abline(v=0,lty=2)
breaks <- seq(-0.6,0.01,0.01)
heights <- vector(length=length(breaks),mode="numeric") 
for(i in 2:length(breaks)){
  heights[i] <- sum(Sc.data[Sc.data$bulk==TRUE & Sc.data$mean.s >= breaks[i-1] & Sc.data$mean.s < breaks[i],"PP.neut"])
}
heights[heights==0] <- NA
points(breaks-0.005, heights, pch=19,col="white",cex=0.5) #this gives proportional stack for histogram representation, have to modify in illustrator to make stacked histogram breaks at points
dev.off()

#zoom in on ~nearly neutral region
pdf(file="plots/Sc-rev_mean-s-histogram_bulk_nearly-neutral_indicate-neutral-bg.pdf",width=5.5,height=5,useDingbats=F)
hist(Sc.data$mean.s[Sc.data$bulk==TRUE & Sc.data$mean.s > -0.1],breaks=seq(-0.05,0.01,0.0025),col="gray50",xlab="selection coefficient of reverse mutation",main="",xlim=c(-0.05,0.01),ylim=c(0,25),right=T);abline(v=0,lty=2)
breaks <- seq(-0.05, 0.01, 0.0025)
heights <- vector(length=length(breaks),mode="numeric") 
for(i in 2:length(breaks)){
  heights[i] <- sum(Sc.data[Sc.data$bulk==TRUE & Sc.data$mean.s >= breaks[i-1] & Sc.data$mean.s < breaks[i],"PP.neut"])
}
heights[heights==0] <- NA
points(breaks-0.00125, heights, pch=19,col="white",cex=0.5) #modify in Illustrator to make stacked histogram bars
dev.off()

# #bootstrap estimate of the proportion deleterious, neutral, beneficial
# B <- 10000 #number bootstrap replicates
# p.del.boot <- vector(mode="numeric",length=B) #vectors to store proportions in each category in each bootstrap
# p.ben.boot <- vector(mode="numeric",length=B)
# p.neut.boot <- vector(mode="numeric",length=B)
# set.seed(1990)
# for(i in 1:B){
#   #sample mean.s values with replacement
#   boot <- sample(Sc.data[Sc.data$bulk==T,"mean.s"], size=length(Sc.data[Sc.data$bulk==T,"mean.s"]),replace=T)
#   #for inferring mixture model, take out outliers
#   boot.bulk <- boot[boot>-0.04]
#   mix.boot <- normalmixEM(boot.bulk, k=3, lambda=rep(1/3, 3), mu=c(mean(Sc.WT$s),-0.01,-0.005), sigma=rep(sd(Sc.WT$s),3), mean.constr=c(mean(Sc.WT$s), rep(NA,2)), sd.constr=c(sd(Sc.WT$s),rep(NA,2)), arbmean=T, arbvar=T)
#   p.del <- sum(mix.boot$posterior[,mix.boot$mu < mean(Sc.WT$s)])/sum(mix.boot$posterior) #proportion of main distribution mutations that are deleterious
#   p.ben <- sum(mix.boot$posterior[,mix.boot$mu > -mean(Sc.WT$s)])/sum(mix.boot$posterior) #proportion of main distribution mutations that are beneficial
#   p.neut <- sum(mix.boot$posterior[,mix.boot$mu >= mean(Sc.WT$s) & mix.boot$mu <= -mean(Sc.WT$s)])/sum(mix.boot$posterior) #proportion of main distribution mutations that are neutral
#   #adjust proportion deleterious for the unambiguously deleterious "outliers" that were not fit in the mixture model
#   p.del.boot[i] <- (p.del*length(boot.bulk) + length(boot)-length(boot.bulk))/length(boot)
#   p.ben.boot[i] <- p.ben*length(boot.bulk)/length(boot)
#   p.neut.boot[i] <- p.neut*length(boot.bulk)/length(boot)
# }
# hist(p.del.boot,breaks=40)
# abline(v=0.9310345)
# quantile(p.del.boot,c(0.025,0.975)) #0.830, 1.000

#ancAmo data
vals <- ancAmo.data$mean.s[ancAmo.data$bulk==TRUE & ancAmo.data$mean.s > -0.04]
hist(vals,breaks=15,col="gray50",xlab="selection coefficient of forward mutation",main="",freq=T);abline(v=0,lty=2)
hist(ancAmo.WT$s,add=T,col="#00FFFF",freq=T,breaks=5)
# #for summary figure, after inferring model can run following commented lines to illustrate mixture over distributions
# curve((mix5$lambda[1]*dnorm(x,mix5$mu[1],mix5$sigma[1]))*0.405,add=T,lty=2,col="darkred",lwd=2)
# curve((mix5$lambda[2]*dnorm(x,mix5$mu[2],mix5$sigma[2]))*0.405,add=T,lty=2,col="darkgreen",lwd=2)
# curve((mix5$lambda[3]*dnorm(x,mix5$mu[3],mix5$sigma[3]))*0.405,add=T,lty=2,col="darkblue",lwd=2)
# curve((mix5$lambda[4]*dnorm(x,mix5$mu[4],mix5$sigma[4]))*0.405,add=T,lty=2,col="orange",lwd=2)
# curve((mix5$lambda[5]*dnorm(x,mix5$mu[5],mix5$sigma[5]))*0.405,add=T,lty=2,col="yellow",lwd=2)
# curve(dmix5(x)*0.405,add=T,lwd=2)

#use mixtools to fit gaussian mixture models, one component of which is defined by the wt sampling distributions, others vary. Figure out best number of components via AIC
n <- 2; mix2 <- normalmixEM(vals, k=n, lambda=rep(1/n, n), mu=c(mean(ancAmo.WT$s),-0.01), sigma=rep(sd(ancAmo.WT$s),n), mean.constr=c(mean(ancAmo.WT$s), rep(NA,n-1)), sd.constr=c(sd(ancAmo.WT$s),rep(NA,n-1)), arbmean=T, arbvar=T)
n <- 3; mix3 <- normalmixEM(vals, k=n, lambda=rep(1/n, n), mu=c(mean(ancAmo.WT$s),-0.02,-0.005), sigma=rep(sd(ancAmo.WT$s),n), mean.constr=c(mean(ancAmo.WT$s), rep(NA,n-1)), sd.constr=c(sd(ancAmo.WT$s),rep(NA,n-1)), arbmean=T, arbvar=T)
n <- 4; mix4 <- normalmixEM(vals, k=n, lambda=rep(1/n, n), mu=c(mean(ancAmo.WT$s),-0.01,-0.005,0.02), sigma=rep(sd(ancAmo.WT$s),n), mean.constr=c(mean(ancAmo.WT$s), rep(NA,n-1)), sd.constr=c(sd(ancAmo.WT$s),rep(NA,n-1)), arbmean=T, arbvar=T)
n <- 5; mix5 <- normalmixEM(vals, k=n, lambda=rep(1/n, n), maxit=2000, mu=c(mean(ancAmo.WT$s),-0.02,-0.001,0.001,0.02), sigma=rep(sd(ancAmo.WT$s),n), mean.constr=c(mean(ancAmo.WT$s), rep(NA,n-1)), sd.constr=c(sd(ancAmo.WT$s),rep(NA,n-1)), arbmean=T, arbvar=T)
n <- 6; mix6 <- normalmixEM(vals, k=n, lambda=rep(1/n, n), maxit=2000, mu=c(mean(ancAmo.WT$s),-0.02,-0.001,-0.001,0.001,0.02), sigma=rep(sd(ancAmo.WT$s),n), mean.constr=c(mean(ancAmo.WT$s), rep(NA,n-1)), sd.constr=c(sd(ancAmo.WT$s),rep(NA,n-1)), arbmean=T, arbvar=T)

summary.mixEM(mix2) #loglik 243.8518 , k = 3, AIC = -481.7036
summary.mixEM(mix3) #loglik 245.4855 , k = 6, AIC = -478.971
summary.mixEM(mix4) #loglik 246.215, k = 9, AIC = -474.43
summary.mixEM(mix5) #loglik 252.74, k = 12, AIC = -481.48
summary.mixEM(mix6) #loglik 252.74 , k = 15, AIC = -475.48

#2-component and 5-component similarly favored by AIC -- tentatively would favor 5-component, it's more conservative (for # deleterious) and just "looks" like a better fit...
#define density functions for these 2- and 5-component mixture model
dmix2 <- function(x){
  return(mix2$lambda[1]*dnorm(x, mix2$mu[1], mix2$sigma[1])+mix2$lambda[2]*dnorm(x, mix2$mu[2], mix2$sigma[2]))
}
dmix5 <- function(x){
  return(mix5$lambda[1]*dnorm(x, mix5$mu[1], mix5$sigma[1])+mix5$lambda[2]*dnorm(x, mix5$mu[2], mix5$sigma[2])+mix5$lambda[3]*dnorm(x, mix5$mu[3], mix5$sigma[3])+mix5$lambda[4]*dnorm(x, mix5$mu[4], mix5$sigma[4])+mix5$lambda[5]*dnorm(x, mix5$mu[5], mix5$sigma[5]))
}

#qqplots to judge model fits
set.seed(101010)
qqplot(c(rnorm(round(100000*mix2$lambda[1]),mix2$mu[1],mix2$sigma[1]),rnorm(round(100000*mix2$lambda[2]),mix2$mu[2],mix2$sigma[2])),vals,xlab="Theoretical quantiles (mixture)", ylab="Sample quantiles",xlim=c(-0.065,0.06),ylim=c(-0.065,0.06))
abline(a=0,b=1)

qqplot(c(rnorm(round(100000*mix5$lambda[1]),mix5$mu[1],mix5$sigma[1]),rnorm(round(100000*mix5$lambda[2]),mix5$mu[2],mix5$sigma[2]),rnorm(round(100000*mix5$lambda[3]),mix5$mu[3],mix5$sigma[3]),rnorm(round(100000*mix5$lambda[4]),mix5$mu[4],mix5$sigma[4]),rnorm(round(100000*mix5$lambda[5]),mix5$mu[5],mix5$sigma[5])),vals,xlab="Theoretical quantiles (mixture)", ylab="Sample quantiles",xlim=c(-0.045,0.03),ylim=c(-0.045,0.03))
abline(a=0,b=1)

#define the posterior probability of being deleterious, neutral, beneficial for each mutation -- comes from the relative pdf of deleterious versus wt versus beneficial components at the given s for each mutation
ancAmo.data$PP.neut <- NA
ancAmo.data$PP.ben <- NA
ancAmo.data$PP.del <- NA
ancAmo.data[ancAmo.data$mean.s < -0.04 & ancAmo.data$bulk == T,"PP.neut"] <- 0
ancAmo.data[ancAmo.data$mean.s < -0.04 & ancAmo.data$bulk == T,"PP.ben"] <- 0
ancAmo.data[ancAmo.data$mean.s < -0.04 & ancAmo.data$bulk == T,"PP.del"] <- 1
ancAmo.data[ancAmo.data$mean.s > -0.04 & ancAmo.data$bulk == T,"PP.neut"] <- mix5$posterior[,1]
ancAmo.data[ancAmo.data$mean.s > -0.04 & ancAmo.data$bulk == T,"PP.ben"] <- mix5$posterior[,4]+mix5$posterior[,5]
ancAmo.data[ancAmo.data$mean.s > -0.04 & ancAmo.data$bulk == T,"PP.del"] <- mix5$posterior[,2]+mix5$posterior[,3]

sum(ancAmo.data$PP.del,na.rm=T)/sum(!is.na(ancAmo.data$PP.del)) # 0.5280663
sum(ancAmo.data$PP.ben,na.rm=T)/sum(!is.na(ancAmo.data$PP.ben)) # 0.151148
sum(ancAmo.data$PP.neut,na.rm=T)/sum(!is.na(ancAmo.data$PP.neut)) # 0.3207857

#plot of distribution of selection coefficients of reversions with an indication of the expected fraction in each bin that are neutral, beneficial
pdf(file="plots/ancAmo-rev_mean-s-histogram_bulk_indicate-neutral-bg.pdf",width=5.5,height=5,useDingbats=F)
hist(ancAmo.data$mean.s[ancAmo.data$bulk==TRUE],breaks=seq(-0.3, 0.05, 0.01),col="gray50",xlab="selection coefficient of forward sub in ancAmo",main="",xlim=c(-0.3,0.05),ylim=c(0,35));abline(v=0,lty=2)
breaks <- seq(-0.3, 0.05, 0.01)
heights.neut <- vector(length=length(breaks), mode="numeric")
heights.ben <- vector(length=length(breaks), mode="numeric")
for(i in 2:length(breaks)){
  heights.neut[i] <- sum(ancAmo.data[ancAmo.data$bulk==TRUE & ancAmo.data$mean.s >= breaks[i-1] & ancAmo.data$mean.s < breaks[i],"PP.neut"])
  heights.ben[i] <- sum(ancAmo.data[ancAmo.data$bulk==TRUE & ancAmo.data$mean.s >= breaks[i-1] & ancAmo.data$mean.s < breaks[i],"PP.ben"])
}
heights.ben <- heights.neut + heights.ben
heights.neut[heights.neut==0] <- NA
heights.ben[heights.ben==0] <- NA
points(breaks-0.004, heights.neut, pch=19,col="white",cex=0.5)
points(breaks-0.006, heights.ben, pch=19,col="cyan",cex=0.5) #modify in Illustrator to make stacked histogram bars
dev.off()

#zoom in on ~nearly neutral region
pdf(file="plots/ancAmo-rev_mean-s-histogram_bulk_nearly-neutral_indicate-neutral-bg.pdf",width=5.5,height=5,useDingbats=F)
hist(ancAmo.data$mean.s[ancAmo.data$bulk==TRUE & ancAmo.data$mean.s > -0.1],breaks=seq(-0.075,0.02,0.005),col="gray50",xlab="selection coefficient of forward mutation",main="",xlim=c(-0.075,0.02),ylim=c(0,20),right=T);abline(v=0,lty=2)
breaks <- seq(-0.075,0.02,0.005)
heights.neut <- vector(length=length(breaks), mode="numeric")
heights.ben <- vector(length=length(breaks), mode="numeric")
for(i in 2:length(breaks)){
  heights.neut[i] <- sum(ancAmo.data[ancAmo.data$bulk==TRUE & ancAmo.data$mean.s >= breaks[i-1] & ancAmo.data$mean.s < breaks[i],"PP.neut"])
  heights.ben[i] <- sum(ancAmo.data[ancAmo.data$bulk==TRUE & ancAmo.data$mean.s >= breaks[i-1] & ancAmo.data$mean.s < breaks[i],"PP.ben"])
}
heights.ben <- heights.neut + heights.ben
heights.neut[heights.neut==0] <- NA
heights.ben[heights.ben==0] <- NA
points(breaks-0.003, heights.neut, pch=19,col="white",cex=0.5)
points(breaks-0.002, heights.ben, pch=19,col="cyan",cex=0.5) #modify in Illustrator to make stacked colors in histogram bins
dev.off()

# #bootstrap estimate of the proportion deleterious, neutral, beneficial
# B <- 10000
# p.del.boot <- vector(mode="numeric",length=B)
# p.ben.boot <- vector(mode="numeric",length=B)
# p.neut.boot <- vector(mode="numeric",length=B)
# set.seed(10990)
# for(i in 1:B){
#   #sample mean.s values with replacement
#   boot <- sample(ancAmo.data[ancAmo.data$bulk==T,"mean.s"], size=length(ancAmo.data[ancAmo.data$bulk==T,"mean.s"]),replace=T)
#   #for inferring mixture model, take out outliers
#   boot.bulk <- boot[boot>-0.04]
#   mix.boot <- normalmixEM(boot.bulk, k=5, lambda=rep(1/5, 5), maxit=2000, mu=c(mean(ancAmo.WT$s),-0.02,-0.001,0.001,0.02), sigma=rep(sd(ancAmo.WT$s),5), mean.constr=c(mean(ancAmo.WT$s), rep(NA,4)), sd.constr=c(sd(ancAmo.WT$s),rep(NA,4)), arbmean=T, arbvar=T)
#   p.del <- sum(mix.boot$posterior[,mix.boot$mu < mean(ancAmo.WT$s)])/sum(mix.boot$posterior) #proportion of main distribution mutations that are deleterious
#   p.ben <- sum(mix.boot$posterior[,mix.boot$mu > -mean(ancAmo.WT$s)])/sum(mix.boot$posterior) #proportion of main distribution mutations that are beneficial
#   p.neut <- sum(mix.boot$posterior[,mix.boot$mu >= mean(ancAmo.WT$s) & mix.boot$mu <= -mean(ancAmo.WT$s)])/sum(mix.boot$posterior) #proportion of main distribution mutations that are neutral
#   #adjust proportion deleterious for the unambiguously deleterious "outliers" that were not fit in the mixture model
#   p.del.boot[i] <- (p.del*length(boot.bulk) + length(boot)-length(boot.bulk))/length(boot)
#   p.ben.boot[i] <- p.ben*length(boot.bulk)/length(boot)
#   p.neut.boot[i] <- p.neut*length(boot.bulk)/length(boot)
# }
# hist(p.del.boot)
# abline(v=0.5280663)
# quantile(p.del.boot,c(0.025,0.975)) #0.274, 0.955
# 
# hist(p.ben.boot)
# abline(v=0.151148)
# quantile(p.ben.boot,c(0.025,0.975)) #0.034, 0.573
# 
# hist(p.neut.boot)
# abline(v=0.151148)
# quantile(p.neut.boot,c(0.025,0.975)) #0.000, 0.593


################################################################################################################################
#one metric of extent of epistasis is to compare the expected fitness of a genotypes with multiple mutations to the expected fitness if mutations were additive
alignment <- read.fasta("data_in/hsp90a_ancestors_Saccer-numbering.fas", seqtype="AA", as.string=FALSE)

#expected s/w of ancAmo in absence of epistasis given s of individual states in ancAmo as determined in Sc background
s.ancAmo <- 0
s.se.ancAmo <- 0
for(i in 1:nrow(Sc.data)){
  if(alignment[[36]][Sc.data[i,"Sc.position"]]==Sc.data[i,"anc.AA"]){
    if(Sc.data[i,"bulk"] == TRUE){
      s.ancAmo <- s.ancAmo+Sc.data[i,"mean.s"]
      s.se.ancAmo <- sqrt(s.se.ancAmo^2+Sc.data[i,"SE.mean.s"]^2)
    }
  }
}
exp(s.ancAmo)
exp(s.ancAmo-1.96*s.se.ancAmo);exp(s.ancAmo+1.96*s.se.ancAmo)
#fitness of ancAmo in absence of epistasis would be 0.233 (95% CI 0.213, 0.256)

s.ancAsco <- 0
s.se.ancAsco <- 0
for(i in 1:nrow(Sc.data)){
  if(alignment[[25]][Sc.data[i,"Sc.position"]]==Sc.data[i,"anc.AA"]){
    if(Sc.data[i,"bulk"] == TRUE){
      s.ancAsco <- s.ancAsco+Sc.data[i,"mean.s"]
      s.se.ancAsco <- sqrt(s.se.ancAsco^2+Sc.data[i,"SE.mean.s"]^2)
    }
  }
}
exp(s.ancAsco)
exp(s.ancAsco-1.96*s.se.ancAsco);exp(s.ancAsco+1.96*s.se.ancAsco)
#fitness of ancAsco in absence of epistasis would be 0.647 (95% CI 0.612, 0.686)

#actual estimate of ancAmo fitness as experimentally measured?
ancAmo.mean.w; ancAmo.mean.w - 1.96*ancAmo.se.w; ancAmo.mean.w + 1.96*ancAmo.se.w #0.429 (95% CI 0.424, 0.435)
#and with additional L378i reversion to the ancAmo state
ancAmo.i378.mean.w;ancAmo.i378.mean.w - 1.96*ancAmo.i378.se.w;ancAmo.i378.mean.w + 1.96*ancAmo.i378.se.w #0.9646 (95% CI 0.9629, 0.9663)

#observed fitness of ancAsco as calculated in isogenic cultures
ancAsco.iso <- read.csv(file="data_in/ancAsco-complementation.csv",header=T,stringsAsFactors=F)
timepoints <- c(16,18,20,21,23.5,28)
for(i in 1:nrow(ancAsco.iso)){
  fit <- lm(as.numeric(ancAsco.iso[i,2:7]) ~ timepoints)
  ancAsco.iso$slope[i] <- summary(fit)$coefficients[2,1] #slope of ln(OD) versus time is malthusian parameter
  ancAsco.iso$slope.se[i] <- summary(fit)$coefficients[2,2]
}
#calculate dimensionless per-generation selection coefficient, assuming that seleciton is density-independent, as (slope_x - slope_wt)/slope_wt*ln2 (last part equivalent to multiplying by generation time of WT)
ancAsco.iso$s <- (ancAsco.iso$slope-ancAsco.iso$slope[1])/(ancAsco.iso$slope[1])*log(2)
#linear transform values such that relative fitness of null allele is 0
#convert raw s to fitness given w=e^s
ancAsco.iso$raw.fitness <- exp(ancAsco.iso$s)
#stretch relative fitness space to be from 0 (null) to 1 (wt)
ancAsco.iso$rel.fitness <- (ancAsco.iso$raw.fitness - ancAsco.iso$raw.fitness[2])/(ancAsco.iso$raw.fitness[1] - ancAsco.iso$raw.fitness[2])
#convert back to selection coefficient s=ln(w)
ancAsco.iso$s <- log(ancAsco.iso$rel.fitness)

#make figure: expected fitness of ancestors versus observed
ancestors <- read.csv(file="data_in/ancestors_info.csv",header=T,stringsAsFactors = F)
#predicted fitness of ancestors given effects of mutations in Sc background?
for(i in 1:nrow(ancestors)){
  s.anc <- 0
  s.se.anc <- 0
  s.anc.2 <- 0 #if excluding large effect (s < -0.1)
  s.se.anc.2 <- 0
  for(j in 1:nrow(Sc.data)){
    if(alignment[[i]][Sc.data[j,"Sc.position"]]==Sc.data[j,"anc.AA"]){
      if(Sc.data[j,"bulk"]==TRUE){
        s.anc <- s.anc + Sc.data[j,"mean.s"]
        s.se.anc <- sqrt(s.se.anc^2 + Sc.data[j,"SE.mean.s"]^2)
      }
      if(Sc.data[j,"bulk"]==TRUE & Sc.data[j,"mean.s"] > -0.1){
        s.anc.2 <- s.anc.2 + Sc.data[j,"mean.s"]
        s.se.anc.2 <- sqrt(s.se.anc.2^2 + Sc.data[j,"SE.mean.s"]^2)
      }
    }
  }
  ancestors$s.predicted.from.Sc[i] <- s.anc
  ancestors$s.SE.predicted.from.Sc[i] <- s.se.anc
  ancestors$s.predicted.from.Sc.small.effect[i] <- s.anc.2
  ancestors$s.SE.predicted.from.Sc.small.effect[i] <- s.se.anc.2
}

pdf(file="plots/predicted_v_observed_anc-fitness_Sc-values_log-fitness.pdf",5,5,useDingbats = F)
plot(ancestors$distance.from.hsp82,ancestors$s.predicted.from.Sc,pch="",xlab="evolutionary distance from ScHsp82 (expected # subs/site)",ylab="selection coefficient (log-fitness)")
arrows(ancestors$distance.from.hsp82,ancestors$s.predicted.from.Sc.small.effect - 1.96*ancestors$s.SE.predicted.from.Sc.small.effect, ancestors$distance.from.hsp82,ancestors$s.predicted.from.Sc.small.effect + 1.96*ancestors$s.SE.predicted.from.Sc.small.effect,length=0.025,angle=90,code=3,col="gray60")
points(ancestors$distance.from.hsp82,ancestors$s.predicted.from.Sc.small.effect,pch=20,col="gray60")
arrows(ancestors$distance.from.hsp82,ancestors$s.predicted.from.Sc - 1.96*ancestors$s.SE.predicted.from.Sc, ancestors$distance.from.hsp82,ancestors$s.predicted.from.Sc + 1.96*ancestors$s.SE.predicted.from.Sc,length=0.025,angle=90,code=3,col="gray30")
points(ancestors$distance.from.hsp82,ancestors$s.predicted.from.Sc,pch=20)
arrows(ancestors$distance.from.hsp82[36],log(ancAmo.i378.mean.w-1.96*ancAmo.i378.se.w),ancestors$distance.from.hsp82[36],log(ancAmo.i378.mean.w+1.96*ancAmo.i378.se.w),length=0.025,angle=90,code=3,col="blue")
points(ancestors$distance.from.hsp82[36],log(ancAmo.i378.mean.w),pch=20,col="blue",cex=1.4)
points(ancestors$distance.from.hsp82[25],log(0.9905),pch=20,col="green",cex=1.4) #fitness of ancAsco measured in isogenic growth relative to Sc
dev.off()

#predicted fitness of ancestors given effects of mutations in ancAmo background?
for(i in 1:nrow(ancestors)){
  s.anc <- 0
  s.se.anc <- 0
  s.anc.2 <- 0 #if excluding large effect (s < -0.1)
  s.se.anc.2 <- 0
  for(j in 1:nrow(ancAmo.data)){
    if(alignment[[i]][ancAmo.data[j,"Sc.position"]]==ancAmo.data[j,"der.AA"]){
      if(ancAmo.data[j,"bulk"]==TRUE){
        s.anc <- s.anc + ancAmo.data[j,"mean.s"]
        s.se.anc <- sqrt(s.se.anc^2 + ancAmo.data[j,"SE.mean.s"]^2)
      }
      if(ancAmo.data[j,"bulk"]==TRUE & ancAmo.data[j,"mean.s"] > -0.1){
        s.anc.2 <- s.anc.2 + ancAmo.data[j,"mean.s"]
        s.se.anc.2 <- sqrt(s.se.anc.2^2 + ancAmo.data[j,"SE.mean.s"]^2)
      }
    }
  }
  ancestors$s.predicted.from.ancAmo[i] <- s.anc
  ancestors$s.SE.predicted.from.ancAmo[i] <- s.se.anc
  ancestors$s.predicted.from.ancAmo.small.effect[i] <- s.anc.2
  ancestors$s.SE.predicted.from.ancAmo.small.effect[i] <- s.se.anc.2
}

pdf(file="plots/predicted_v_observed_anc-fitness_ancAmo-values_log-fitness.pdf",5,5,useDingbats = F)
plot(ancestors$distance.from.ancAmo,ancestors$s.predicted.from.ancAmo,pch="",xlab="evolutionary distance from ancAmoHsp90 (expected # subs/site)",ylab="selection coefficient (log-fitness)",ylim=c(-1,0.1))
arrows(ancestors$distance.from.ancAmo,ancestors$s.predicted.from.ancAmo.small.effect - 1.96*ancestors$s.SE.predicted.from.ancAmo.small.effect, ancestors$distance.from.ancAmo,ancestors$s.predicted.from.ancAmo.small.effect + 1.96*ancestors$s.SE.predicted.from.ancAmo.small.effect,length=0.025,angle=90,code=3,col="gray60")
points(ancestors$distance.from.ancAmo,ancestors$s.predicted.from.ancAmo.small.effect,pch=20,col="gray60")
arrows(ancestors$distance.from.ancAmo,ancestors$s.predicted.from.ancAmo - 1.96*ancestors$s.SE.predicted.from.ancAmo, ancestors$distance.from.ancAmo,ancestors$s.predicted.from.ancAmo + 1.96*ancestors$s.SE.predicted.from.ancAmo,length=0.025,angle=90,code=3,col="gray30")
points(ancestors$distance.from.ancAmo,ancestors$s.predicted.from.ancAmo,pch=20)
arrows(ancestors$distance.from.ancAmo[1],log(Sc.mean.w-1.96*Sc.se.w),ancestors$distance.from.ancAmo[1],log(Sc.mean.w+1.96*Sc.se.w),length=0.025,angle=90,code=3,col="red")
points(ancestors$distance.from.ancAmo[1],log(Sc.mean.w),pch=20,col="red",cex=1.4)
points(ancestors$distance.from.ancAmo[25],log(1.0269),pch=20,col="green",cex=1.4) #fitness of ancAsco measured in isogenic growth relative to Sc=1.036698
dev.off()

################################################################################################################################
#compare effect of equivalent mutation between backgrounds
#add column to ancAmo.data that reflects the selection coefficient of the equivalent mutation in the reversion data (i.e. change the directionality of selection coefficient if needed)
for (i in 1:nrow(ancAmo.data)){
  if(nrow(ancAmo.data[ancAmo.data$Sc.position==ancAmo.data$Sc.position[i],])==1){ #if there's only one alt state... it's either reverse changes or the same change
    if(ancAmo.data$shorthand[i] %in% Sc.data$shorthand){ #if it's the same change
      ancAmo.data$mean.s.rev[i] <- Sc.data[Sc.data$shorthand==ancAmo.data$shorthand[i],"mean.s"] #the mean.s is already polarized in the same direction
      ancAmo.data$SE.mean.s.rev[i] <- Sc.data[Sc.data$shorthand==ancAmo.data$shorthand[i],"SE.mean.s"]
      ancAmo.data$direction[i] <- 0 #factor for difference in direction between fwd, rev sub: 0 means same sub (parallel); 1 means rev sub (reversion); 2 means sub to same state from different states (convergent)
      ancAmo.data$Sc.PP.neut[i] <- Sc.data[Sc.data$shorthand==ancAmo.data$shorthand[i],"PP.neut"]
      ancAmo.data$Sc.PP.del[i] <- Sc.data[Sc.data$shorthand==ancAmo.data$shorthand[i],"PP.del"]
      ancAmo.data$Sc.PP.ben[i] <- Sc.data[Sc.data$shorthand==ancAmo.data$shorthand[i],"PP.ben"]
    } else { #otherwise
      ancAmo.data$mean.s.rev[i] <- -Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i],"mean.s"] #the mean.s is in the opposite direction, take the inverse
      ancAmo.data$SE.mean.s.rev[i] <- Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i],"SE.mean.s"]
      ancAmo.data$direction[i] <- 1
      ancAmo.data$Sc.PP.neut[i] <- Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$ancAmo.AA[i],"PP.neut"]
      ancAmo.data$Sc.PP.ben[i] <- Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$ancAmo.AA[i],"PP.del"]
      ancAmo.data$Sc.PP.del[i] <- Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$ancAmo.AA[i],"PP.ben"]
    }
  } else { #if there's more than one alt state
    if(ancAmo.data$shorthand[i] %in% Sc.data$shorthand){ #if Sc and ancAmo have the same state, they'll already all be polarized in the right direction
      ancAmo.data$mean.s.rev[i] <- Sc.data[Sc.data$shorthand==ancAmo.data$shorthand[i],"mean.s"]
      ancAmo.data$SE.mean.s.rev[i] <- Sc.data[Sc.data$shorthand==ancAmo.data$shorthand[i],"SE.mean.s"]
      ancAmo.data$direction[i] <- 0
      ancAmo.data$Sc.PP.neut[i] <- Sc.data[Sc.data$shorthand==ancAmo.data$shorthand[i],"PP.neut"]
      ancAmo.data$Sc.PP.del[i] <- Sc.data[Sc.data$shorthand==ancAmo.data$shorthand[i],"PP.del"]
      ancAmo.data$Sc.PP.ben[i] <- Sc.data[Sc.data$shorthand==ancAmo.data$shorthand[i],"PP.ben"]
    } else if(ancAmo.data$der.AA[i] == ancAmo.data$Sc.AA[i]){ #if Sc and ancAmo have different states, one of the options will be the reverse mutation ...
      ancAmo.data$mean.s.rev[i] <- -Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$ancAmo.AA[i],"mean.s"]
      ancAmo.data$SE.mean.s.rev[i] <- Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$ancAmo.AA[i],"SE.mean.s"]
      ancAmo.data$direction[i] <- 1
      ancAmo.data$Sc.PP.neut[i] <- Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$ancAmo.AA[i],"PP.neut"]
      ancAmo.data$Sc.PP.ben[i] <- Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$ancAmo.AA[i],"PP.del"]
      ancAmo.data$Sc.PP.del[i] <- Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$ancAmo.AA[i],"PP.ben"]
    } else { #and the others will use a shared intermediate state
      ancAmo.data$mean.s.rev[i] <- -Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$ancAmo.AA[i],"mean.s"] + Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$der.AA[i],"mean.s"]
      ancAmo.data$SE.mean.s.rev[i] <- sqrt((Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$ancAmo.AA[i],"SE.mean.s"])^2 + (Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$der.AA[i],"SE.mean.s"])^2)
      ancAmo.data$direction[i] <- 2
      ancAmo.data$Sc.PP.neut[i] <- Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$der.AA[i],"PP.neut"] #note, this is NOT the PP that this mutation in this column is neut (or del, ben below) as it is for the class 0 and class 1 subs, but rather the PP that the intermediate state j was neutral if reversed in Sc (for calculation below, see ***)
      ancAmo.data$Sc.PP.del[i] <- Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$der.AA[i],"PP.del"]
      ancAmo.data$Sc.PP.ben[i] <- Sc.data[Sc.data$Sc.position==ancAmo.data$Sc.position[i] & Sc.data$Sc.AA==ancAmo.data$Sc.AA[i] & Sc.data$anc.AA == ancAmo.data$der.AA[i],"PP.ben"]
    }
  }
}
ancAmo.data$direction <- as.factor(ancAmo.data$direction)

#take out isolated measurements in ancAmo background which have no estimates of error
ancAmo.data.bulk <- ancAmo.data[ancAmo.data$bulk==TRUE,]

#fraction of all subs that are neither contingent NOR irreversible?
#need to do separately for each class
#first, class 0 (i->j->i); PP.del is prob contingent, Sc.PP.del is prob entrenched
x <- sum(ancAmo.data.bulk[ancAmo.data.bulk$direction==0,"Sc.PP.neut"]*ancAmo.data.bulk[ancAmo.data.bulk$direction==0,"PP.neut"] + ancAmo.data.bulk[ancAmo.data.bulk$direction==0,"Sc.PP.neut"]*ancAmo.data.bulk[ancAmo.data.bulk$direction==0,"PP.ben"] + ancAmo.data.bulk[ancAmo.data.bulk$direction==0,"Sc.PP.ben"]*ancAmo.data.bulk[ancAmo.data.bulk$direction==0,"PP.neut"] + ancAmo.data.bulk[ancAmo.data.bulk$direction==0,"Sc.PP.ben"]*ancAmo.data.bulk[ancAmo.data.bulk$direction==0,"PP.ben"],na.rm=T)
#next, class 1 (i->j); PP.del is prob contingent, Sc.PP.ben is prob entrenched
x <- x + sum(ancAmo.data.bulk[ancAmo.data.bulk$direction==1,"Sc.PP.neut"]*ancAmo.data.bulk[ancAmo.data.bulk$direction==1,"PP.neut"] + ancAmo.data.bulk[ancAmo.data.bulk$direction==1,"Sc.PP.neut"]*ancAmo.data.bulk[ancAmo.data.bulk$direction==1,"PP.ben"] + ancAmo.data.bulk[ancAmo.data.bulk$direction==1,"Sc.PP.del"]*ancAmo.data.bulk[ancAmo.data.bulk$direction==1,"PP.neut"] + ancAmo.data.bulk[ancAmo.data.bulk$direction==1,"Sc.PP.del"]*ancAmo.data.bulk[ancAmo.data.bulk$direction==1,"PP.ben"],na.rm=T)
#*** finally, class 2 (i->j->k); PP.del is prob that k is contingent, Sc.PP.del is prob k is irreversible (remember this Sc.PP.del is not for the substituion in that row, but rather for that intermediate state j, coming from different starting points in ancAmo versus Sc)
x <- x + sum(ancAmo.data.bulk[ancAmo.data.bulk$direction==2,"Sc.PP.neut"]*ancAmo.data.bulk[ancAmo.data.bulk$direction==2,"PP.neut"] + ancAmo.data.bulk[ancAmo.data.bulk$direction==2,"Sc.PP.neut"]*ancAmo.data.bulk[ancAmo.data.bulk$direction==2,"PP.ben"] + ancAmo.data.bulk[ancAmo.data.bulk$direction==2,"Sc.PP.ben"]*ancAmo.data.bulk[ancAmo.data.bulk$direction==2,"PP.neut"] + ancAmo.data.bulk[ancAmo.data.bulk$direction==2,"Sc.PP.ben"]*ancAmo.data.bulk[ancAmo.data.bulk$direction==2,"PP.ben"],na.rm=T)
x/(nrow(ancAmo.data.bulk)-1) #5.01%

#for each mutation, compute joint PP of various combinations of being contingent and irreversible
for(i in 1:nrow(ancAmo.data.bulk)){
  if(ancAmo.data.bulk$direction[i]==0){ #first, class 0 (i->j->i); PP.del is prob contingent, Sc.PP.del is prob entrenched
    ancAmo.data.bulk$PP.del.del[i] <- ancAmo.data.bulk[i,"PP.del"]*ancAmo.data.bulk[i,"Sc.PP.del"] #PP that ancestral state preferred in ancAmo, derived state preferred in Sc
    ancAmo.data.bulk$PP.del.neut[i] <- ancAmo.data.bulk[i,"PP.del"]*ancAmo.data.bulk[i,"Sc.PP.neut"] #PP that ancestral state preferred in ancAmo, neither state preferred in Sc
    ancAmo.data.bulk$PP.del.ben[i] <- ancAmo.data.bulk[i,"PP.del"]*ancAmo.data.bulk[i,"Sc.PP.ben"] #PP that ancestral state preferred in ancAmo, ancestral state preferred in Sc
    ancAmo.data.bulk$PP.neut.del[i] <- ancAmo.data.bulk[i,"PP.neut"]*ancAmo.data.bulk[i,"Sc.PP.del"] #PP that neither state preferred in ancAmo, derived state preferred in Sc
    ancAmo.data.bulk$PP.neut.neut[i] <- ancAmo.data.bulk[i,"PP.neut"]*ancAmo.data.bulk[i,"Sc.PP.neut"] #PP that neither state preferred in ancAmo, neither state preferred in Sc
    ancAmo.data.bulk$PP.neut.ben[i] <- ancAmo.data.bulk[i,"PP.neut"]*ancAmo.data.bulk[i,"Sc.PP.ben"] #PP that neither state preferred in ancAmo, ancestral state preferred in Sc
    ancAmo.data.bulk$PP.ben.del[i] <- ancAmo.data.bulk[i,"PP.ben"]*ancAmo.data.bulk[i,"Sc.PP.del"] #PP that derived state preferred in ancAmo, derived state preferred in Sc
    ancAmo.data.bulk$PP.ben.neut[i] <- ancAmo.data.bulk[i,"PP.ben"]*ancAmo.data.bulk[i,"Sc.PP.neut"] #PP that derived state preferred in ancAmo, neither state preferred in Sc
    ancAmo.data.bulk$PP.ben.ben[i] <- ancAmo.data.bulk[i,"PP.ben"]*ancAmo.data.bulk[i,"Sc.PP.ben"] #PP that derived state preferred in ancAmo, ancestral state preferred in Sc
  }else if(ancAmo.data.bulk$direction[i]==1){ #next, class 1 (i->j); PP.del is prob contingent, Sc.PP.ben is prob entrenched
    ancAmo.data.bulk$PP.del.del[i] <- ancAmo.data.bulk[i,"PP.del"]*ancAmo.data.bulk[i,"Sc.PP.ben"] #PP that ancestral state preferred in ancAmo, derived state preferred in Sc
    ancAmo.data.bulk$PP.del.neut[i] <- ancAmo.data.bulk[i,"PP.del"]*ancAmo.data.bulk[i,"Sc.PP.neut"] #PP that ancestral state preferred in ancAmo, neither state preferred in Sc
    ancAmo.data.bulk$PP.del.ben[i] <- ancAmo.data.bulk[i,"PP.del"]*ancAmo.data.bulk[i,"Sc.PP.del"] #PP that ancestral state preferred in ancAmo, ancestral state preferred in Sc
    ancAmo.data.bulk$PP.neut.del[i] <- ancAmo.data.bulk[i,"PP.neut"]*ancAmo.data.bulk[i,"Sc.PP.ben"] #PP that neither state preferred in ancAmo, derived state preferred in Sc
    ancAmo.data.bulk$PP.neut.neut[i] <- ancAmo.data.bulk[i,"PP.neut"]*ancAmo.data.bulk[i,"Sc.PP.neut"] #PP that neither state preferred in ancAmo, neither state preferred in Sc
    ancAmo.data.bulk$PP.neut.ben[i] <- ancAmo.data.bulk[i,"PP.neut"]*ancAmo.data.bulk[i,"Sc.PP.del"] #PP that neither state preferred in ancAmo, ancestral state preferred in Sc
    ancAmo.data.bulk$PP.ben.del[i] <- ancAmo.data.bulk[i,"PP.ben"]*ancAmo.data.bulk[i,"Sc.PP.ben"] #PP that derived state preferred in ancAmo, derived state preferred in Sc
    ancAmo.data.bulk$PP.ben.neut[i] <- ancAmo.data.bulk[i,"PP.ben"]*ancAmo.data.bulk[i,"Sc.PP.neut"] #PP that derived state preferred in ancAmo, neither state preferred in Sc
    ancAmo.data.bulk$PP.ben.ben[i] <- ancAmo.data.bulk[i,"PP.ben"]*ancAmo.data.bulk[i,"Sc.PP.del"] #PP that derived state preferred in ancAmo, ancestral state preferred in Sc
  }else if(ancAmo.data.bulk$direction[i]==2){ #*** finally, class 2 (i->j->k); PP.del is prob that k is contingent, Sc.PP.del is prob k is irreversible (remember this Sc.PP.del is not for the substituion in that row, but rather for that intermediate state j, coming from different starting points in ancAmo versus Sc)
    ancAmo.data.bulk$PP.del.del[i] <- ancAmo.data.bulk[i,"PP.del"]*ancAmo.data.bulk[i,"Sc.PP.del"] #PP that ancestral state preferred in ancAmo, derived state preferred in Sc
    ancAmo.data.bulk$PP.del.neut[i] <- ancAmo.data.bulk[i,"PP.del"]*ancAmo.data.bulk[i,"Sc.PP.neut"] #PP that ancestral state preferred in ancAmo, neither state preferred in Sc
    ancAmo.data.bulk$PP.del.ben[i] <- ancAmo.data.bulk[i,"PP.del"]*ancAmo.data.bulk[i,"Sc.PP.ben"] #PP that ancestral state preferred in ancAmo, ancestral state preferred in Sc
    ancAmo.data.bulk$PP.neut.del[i] <- ancAmo.data.bulk[i,"PP.neut"]*ancAmo.data.bulk[i,"Sc.PP.del"] #PP that neither state preferred in ancAmo, derived state preferred in Sc
    ancAmo.data.bulk$PP.neut.neut[i] <- ancAmo.data.bulk[i,"PP.neut"]*ancAmo.data.bulk[i,"Sc.PP.neut"] #PP that neither state preferred in ancAmo, neither state preferred in Sc
    ancAmo.data.bulk$PP.neut.ben[i] <- ancAmo.data.bulk[i,"PP.neut"]*ancAmo.data.bulk[i,"Sc.PP.ben"] #PP that neither state preferred in ancAmo, ancestral state preferred in Sc
    ancAmo.data.bulk$PP.ben.del[i] <- ancAmo.data.bulk[i,"PP.ben"]*ancAmo.data.bulk[i,"Sc.PP.del"] #PP that derived state preferred in ancAmo, derived state preferred in Sc
    ancAmo.data.bulk$PP.ben.neut[i] <- ancAmo.data.bulk[i,"PP.ben"]*ancAmo.data.bulk[i,"Sc.PP.neut"] #PP that derived state preferred in ancAmo, neither state preferred in Sc
    ancAmo.data.bulk$PP.ben.ben[i] <- ancAmo.data.bulk[i,"PP.ben"]*ancAmo.data.bulk[i,"Sc.PP.ben"] #PP that derived state preferred in ancAmo, ancestral state preferred in Sc
  }
}
sum(ancAmo.data.bulk$PP.del.ben,ancAmo.data.bulk$PP.del.neut,ancAmo.data.bulk$PP.del.del,ancAmo.data.bulk$PP.neut.ben,ancAmo.data.bulk$PP.neut.neut,ancAmo.data.bulk$PP.neut.del,ancAmo.data.bulk$PP.ben.ben,ancAmo.data.bulk$PP.ben.neut,ancAmo.data.bulk$PP.ben.del,na.rm=T)
sum(ancAmo.data.bulk$PP.del.del,na.rm=T)/86 # 0.502 proportion ancestral state preferred in ancAmo, derived state preferred in Sc (contingent and entrenched)
sum(ancAmo.data.bulk$PP.del.neut,na.rm=T)/86 # 0.025 proportion ancestral state preferred in ancAmo, neither state preferred in Sc (contingent)
sum(ancAmo.data.bulk$PP.del.ben,na.rm=T)/86 # 0 proportion ancestral state preferred in both ancAmo and Sc
sum(ancAmo.data.bulk$PP.neut.del,na.rm=T)/86 # 0.301 proportion neither state preferred in ancAmo, derived state preferred in Sc (irreversible, entrenched)
sum(ancAmo.data.bulk$PP.neut.neut,na.rm=T)/86 # 0.019 proportion neither state preferred in both ancAmo and Sc
sum(ancAmo.data.bulk$PP.neut.ben,na.rm=T)/86 # 0 proportion neither state preferred in ancAmo, ancestral state preferred in Sc
sum(ancAmo.data.bulk$PP.ben.del,na.rm=T)/86 # 0.122 proportion derived state preferred in both ancAmo and Sc (irreversible, though not entrenched)
sum(ancAmo.data.bulk$PP.ben.neut,na.rm=T)/86 # 0.031 proportion derived state preferred in ancAmo, neither state preferred in Sc
sum(ancAmo.data.bulk$PP.ben.ben,na.rm=T)/86 # 0 proportion derived state preferred in ancAmo, ancestral state preferred in Sc

sum(ancAmo.data.bulk$PP.del.ben+ancAmo.data.bulk$PP.neut.neut+ancAmo.data.bulk$PP.ben.del,na.rm=T)/86 #fraction nonepistatic: 0.141 (though some of these labeled "nonepistatic" can still be contingent or irreversible)
1-sum(ancAmo.data.bulk$PP.del.ben+ancAmo.data.bulk$PP.neut.neut+ancAmo.data.bulk$PP.ben.del,na.rm=T)/86 #fraction epistatic: 0.859


pdf(file="plots/joint-probabilities_del-neut-ben_all.pdf",width=5,height=5,useDingbats=F)
barplot(rev(c(sum(ancAmo.data.bulk$PP.del.del,na.rm=T)/86, # 0.502 proportion ancestral state preferred in ancAmo, derived state preferred in Sc (contingent and irreversible)
          sum(ancAmo.data.bulk$PP.del.neut,na.rm=T)/86, # 0.025 proportion ancestral state preferred in ancAmo, neither state preferred in Sc (contingent)
          sum(ancAmo.data.bulk$PP.del.ben,na.rm=T)/86, # 0 proportion ancestral state preferred in both ancAmo and Sc (contingent)
          sum(ancAmo.data.bulk$PP.neut.del,na.rm=T)/86, # 0.301 proportion neither state preferred in ancAmo, derived state preferred in Sc (irreversible)
          sum(ancAmo.data.bulk$PP.neut.neut,na.rm=T)/86, # 0.019 proportion neither state preferred in both ancAmo and Sc
          sum(ancAmo.data.bulk$PP.neut.ben,na.rm=T)/86, # 0 proportion neither state preferred in ancAmo, ancestral state preferred in Sc
          sum(ancAmo.data.bulk$PP.ben.del,na.rm=T)/86, # 0.122 proportion derived state preferred in both ancAmo and Sc (irreversible)
          sum(ancAmo.data.bulk$PP.ben.neut,na.rm=T)/86, # 0.031 proportion derived state preferred in ancAmo, neither state preferred in Sc
          sum(ancAmo.data.bulk$PP.ben.ben,na.rm=T)/86)),xlab="proportion substitutions",horiz=T) # 0 proportion derived state preferred in ancAmo, ancestral state preferred in Sc
dev.off()

#correlation in mutational effects between backgrounds, for i->j substitutions
cor.test(ancAmo.data.bulk[ancAmo.data.bulk$mean.s > -0.1 & ancAmo.data.bulk$mean.s < 0.1 & ancAmo.data.bulk$mean.s.rev > -0.1 & ancAmo.data.bulk$mean.s.rev < 0.1 & ancAmo.data.bulk$direction==1,"mean.s.rev"], ancAmo.data.bulk[ancAmo.data.bulk$mean.s > -0.1 & ancAmo.data.bulk$mean.s < 0.1 & ancAmo.data.bulk$mean.s.rev > -0.1 & ancAmo.data.bulk$mean.s.rev < 0.1 & ancAmo.data.bulk$direction==1,"mean.s"])
#r = -0.36, P = 0.009

pdf(file="plots/Sc-s_vs_AncAmo-s_only-i-j-points.pdf",width=5.5,height=5,useDingbats=F)
p3 <- ggplot(data=ancAmo.data.bulk[ancAmo.data.bulk$mean.s > -0.1 & ancAmo.data.bulk$mean.s < 0.1 & ancAmo.data.bulk$mean.s.rev > -0.1 & ancAmo.data.bulk$mean.s.rev < 0.1 & ancAmo.data.bulk$direction==1,],aes(x=mean.s,y=mean.s.rev))
p3+geom_errorbar(aes(x=ancAmo.data.bulk[ancAmo.data.bulk$mean.s > -0.1 & ancAmo.data.bulk$mean.s < 0.1 & ancAmo.data.bulk$mean.s.rev > -0.1 & ancAmo.data.bulk$mean.s.rev < 0.1 & ancAmo.data.bulk$direction==1,"mean.s"],ymin=ancAmo.data.bulk[ancAmo.data.bulk$mean.s > -0.1 & ancAmo.data.bulk$mean.s < 0.1 & ancAmo.data.bulk$mean.s.rev > -0.1 & ancAmo.data.bulk$mean.s.rev < 0.1 & ancAmo.data.bulk$direction==1,"mean.s.rev"]-ancAmo.data.bulk[ancAmo.data.bulk$mean.s > -0.1 & ancAmo.data.bulk$mean.s < 0.1 & ancAmo.data.bulk$mean.s.rev > -0.1 & ancAmo.data.bulk$mean.s.rev < 0.1 & ancAmo.data.bulk$direction==1,"SE.mean.s.rev"],ymax=ancAmo.data.bulk[ancAmo.data.bulk$mean.s > -0.1 & ancAmo.data.bulk$mean.s < 0.1 & ancAmo.data.bulk$mean.s.rev > -0.1 & ancAmo.data.bulk$mean.s.rev < 0.1 & ancAmo.data.bulk$direction==1,"mean.s.rev"]+ancAmo.data.bulk[ancAmo.data.bulk$mean.s > -0.1 & ancAmo.data.bulk$mean.s < 0.1 & ancAmo.data.bulk$mean.s.rev > -0.1 & ancAmo.data.bulk$mean.s.rev < 0.1 & ancAmo.data.bulk$direction==1,"SE.mean.s.rev"]),width=0.0005,colour="gray70")+
  geom_errorbarh(aes(y=ancAmo.data.bulk[ancAmo.data.bulk$mean.s > -0.1 & ancAmo.data.bulk$mean.s < 0.1 & ancAmo.data.bulk$mean.s.rev > -0.1 & ancAmo.data.bulk$mean.s.rev < 0.1 & ancAmo.data.bulk$direction==1,"mean.s.rev"],xmin=ancAmo.data.bulk[ancAmo.data.bulk$mean.s > -0.1 & ancAmo.data.bulk$mean.s < 0.1 & ancAmo.data.bulk$mean.s.rev > -0.1 & ancAmo.data.bulk$mean.s.rev < 0.1 & ancAmo.data.bulk$direction==1,"mean.s"]-ancAmo.data.bulk[ancAmo.data.bulk$mean.s > -0.1 & ancAmo.data.bulk$mean.s < 0.1 & ancAmo.data.bulk$mean.s.rev > -0.1 & ancAmo.data.bulk$mean.s.rev < 0.1 & ancAmo.data.bulk$direction==1,"SE.mean.s"],xmax=ancAmo.data.bulk[ancAmo.data.bulk$mean.s > -0.1 & ancAmo.data.bulk$mean.s < 0.1 & ancAmo.data.bulk$mean.s.rev > -0.1 & ancAmo.data.bulk$mean.s.rev < 0.1 & ancAmo.data.bulk$direction==1,"mean.s"]+ancAmo.data.bulk[ancAmo.data.bulk$mean.s > -0.1 & ancAmo.data.bulk$mean.s < 0.1 & ancAmo.data.bulk$mean.s.rev > -0.1 & ancAmo.data.bulk$mean.s.rev < 0.1 & ancAmo.data.bulk$direction==1,"SE.mean.s"]),height=0.0005,colour="gray70")+
  geom_point(size=3.5)+
  scale_colour_manual(values="black")+
  theme_classic()+
  theme(text=element_text(size=16))+
  xlab(expression(italic(s)[ij]~","~AncAmoHsp90))+ylab(expression(italic(s)[ij]~","~ScHsp90))
dev.off()

########################################################################################################################
#want to look for features that might be associated with the degree to which a mutation is deleterious in each of the datasets
#some properties already present in data table, add remaining ones I want to examine

#overall evolutionary rate of the site. Do epistatically constrained sites evolve more slowly?
#read in relative evolutionary rate at each site, from PAML (in alignment numbering)
evol.rate <- read.table(file="data_in/paml-rates.txt",header=T,sep="\t")
for(i in 1:nrow(Sc.data)){
  Sc.data$evol.rate[i] <- evol.rate[evol.rate$site==Sc.data$alignment.position[i],"rate"]
}
for(i in 1:nrow(ancAmo.data)){
  ancAmo.data$evol.rate[i] <- evol.rate[evol.rate$site==ancAmo.data$alignment.position[i],"rate"]
}

#RSA is known to be a strong determinant of evolutionary rate. Do deleterious states have lower RSA?
#read in table of per-site RSA from the 2CG9 crystal structure, calculated using the script getRelativeSolventAccessibility-NWV.py by D. Allan Drummond and Nicholas VanKuren
RSA <- read.table(file="data_in/RSA.txt",header=T,sep="\t")
for(i in 1:nrow(Sc.data)){
  Sc.data$RSA[i] <- mean(RSA[RSA$site==Sc.data$Sc.position[i],"rsa"])
}
for(i in 1:nrow(ancAmo.data)){
  ancAmo.data$RSA[i] <- mean(RSA[RSA$site==ancAmo.data$Sc.position[i],"rsa"])
}

#distance from the ATP gamma-phosphate was shown by Mishra et al. Cell Reports to constrain mutational tolerance in this Hsp90 NTD
#read in table of distance from each variable site to each ATP gamma phosphate, calculated using pairwisedistances.py script by Pietro Gatti-Lafranconi
dist.to.gP <- read.table(file="data_in/distance-to-ATPgP.txt",header=T,sep="\t")
for(i in 1:nrow(Sc.data)){
  Sc.data$dist.to.gP[i] <- min(dist.to.gP[dist.to.gP$site==Sc.data$Sc.position[i],"distance"])
  if(Sc.data$dist.to.gP[i]==Inf){Sc.data$dist.to.gP[i] <- NA}
}
for(i in 1:nrow(ancAmo.data)){
  ancAmo.data$dist.to.gP[i] <- min(dist.to.gP[dist.to.gP$site==ancAmo.data$Sc.position[i],"distance"])
  if(ancAmo.data$dist.to.gP[i]==Inf){ancAmo.data$dist.to.gP[i] <- NA}
}

#characteristics of actual amino acid states involved in deleterious v neutral mutations?
#overall difference in BLOSUM62 score of aa sub?
data(BLOSUM62)
for(i in 1:nrow(Sc.data)){
  Sc.data$BLOSUM62[i] <- BLOSUM62[Sc.data$Sc.AA[i],Sc.data$anc.AA[i]]
}
for(i in 1:nrow(ancAmo.data)){
  ancAmo.data$BLOSUM62[i] <- BLOSUM62[ancAmo.data$ancAmo.AA[i],ancAmo.data$der.AA[i]]
}

#hydrophobicity: kyte and doolittle scale
data(aaindex)
which(sapply(aaindex, function(x) length(grep("Kyte", x$A)) != 0))
aaindex[[151]]$I
for(i in 1:nrow(Sc.data)){
  Sc.data$diff.hydrophobicity[i] <- aaindex[[151]]$I[aaa(Sc.data[i,"Sc.AA"])]-aaindex[[151]]$I[aaa(Sc.data[i,"anc.AA"])]
}
for(i in 1:nrow(ancAmo.data)){
  ancAmo.data$diff.hydrophobicity[i] <- aaindex[[151]]$I[aaa(ancAmo.data[i,"ancAmo.AA"])]-aaindex[[151]]$I[aaa(ancAmo.data[i,"der.AA"])]
}

#abs value of diff in hydrophobicity?
Sc.data$abs.diff.hydrophobicity <- abs(Sc.data$diff.hydrophobicity)
ancAmo.data$abs.diff.hydrophobicity <- abs(ancAmo.data$diff.hydrophobicity)

#grantham volume
which(sapply(aaindex, function(x) length(grep("Grantham", x$A)) != 0))
aaindex[[112]]$I
for(i in 1:nrow(Sc.data)){
  Sc.data$diff.volume[i] <- aaindex[[112]]$I[aaa(Sc.data[i,"Sc.AA"])]-aaindex[[112]]$I[aaa(Sc.data[i,"anc.AA"])]
}
for(i in 1:nrow(ancAmo.data)){
  ancAmo.data$diff.volume[i] <- aaindex[[112]]$I[aaa(ancAmo.data[i,"ancAmo.AA"])]-aaindex[[112]]$I[aaa(ancAmo.data[i,"der.AA"])]
}

#abs value of diff in volume?
Sc.data$abs.diff.volume <- abs(Sc.data$diff.volume)
ancAmo.data$abs.diff.volume <- abs(ancAmo.data$diff.volume)

#mutational tolerance, from NTD-wide EMPIRIC experiment of Mishra, Flynn et al. Cell Reports 2015
#read in csv giving relative growth rate differences for all mutations to the NTD measured in bulk growth assays
EMPIRIC <- read.csv(file="data_in/Mishra_NTD-EMPIRIC.csv",header=T,stringsAsFactors = F)

for(i in 1:nrow(Sc.data)){
  Sc.data$avg.mut.effect[i] <- mean(EMPIRIC[EMPIRIC$Sc.position==Sc.data$Sc.position[i],"s"],na.rm=T)
}
for(i in 1:nrow(ancAmo.data)){
  ancAmo.data$avg.mut.effect[i] <- mean(EMPIRIC[EMPIRIC$Sc.position==ancAmo.data$Sc.position[i],"s"],na.rm=T)
}

#do I see a rank correlation (Spearman) between properties and mean.s measurements?
#Sc, all muts
var <- "avg.mut.effect"; cor.test(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"])
var <- "evol.rate"; cor.test(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"])
var <- "RSA"; cor.test(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"])
var <- "dist.to.gP"; cor.test(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"])
var <- "BLOSUM62"; cor.test(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"])
var <- "diff.hydrophobicity"; cor.test(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"])
var <- "abs.diff.hydrophobicity"; cor.test(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"])
var <- "diff.volume"; cor.test(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"])
var <- "abs.diff.volume"; cor.test(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T, var],Sc.data[Sc.data$bulk==T, "mean.s"])

#Sc, exclude strongly deleterious outliers
var <- "avg.mut.effect"; cor.test(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"])
var <- "evol.rate"; cor.test(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"])
var <- "RSA"; cor.test(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"])
var <- "dist.to.gP"; cor.test(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"])
var <- "BLOSUM62"; cor.test(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"])
var <- "diff.hydrophobicity"; cor.test(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"])
var <- "abs.diff.hydrophobicity"; cor.test(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"])
var <- "diff.volume"; cor.test(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"])
var <- "abs.diff.volume"; cor.test(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"], method="spearman");plot(Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, var],Sc.data[Sc.data$bulk==T & Sc.data$mean.s > -0.1, "mean.s"])

#ancAmo, all muts
var <- "avg.mut.effect"; cor.test(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"])
var <- "evol.rate"; cor.test(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"])
var <- "RSA"; cor.test(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"])
var <- "dist.to.gP"; cor.test(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"])
var <- "BLOSUM62"; cor.test(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"])
var <- "diff.hydrophobicity"; cor.test(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"])
var <- "abs.diff.hydrophobicity"; cor.test(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"])
var <- "diff.volume"; cor.test(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"])
var <- "abs.diff.volume"; cor.test(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T, var],ancAmo.data[ancAmo.data$bulk==T, "mean.s"])

#ancAmo, exclude strongly deleterious outliers
var <- "avg.mut.effect"; cor.test(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"])
var <- "evol.rate"; cor.test(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"])
var <- "RSA"; cor.test(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"])
var <- "dist.to.gP"; cor.test(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"])
var <- "BLOSUM62"; cor.test(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"])
var <- "diff.hydrophobicity"; cor.test(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"])
var <- "abs.diff.hydrophobicity"; cor.test(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"])
var <- "diff.volume"; cor.test(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"])
var <- "abs.diff.volume"; cor.test(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"], method="spearman");plot(ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, var],ancAmo.data[ancAmo.data$bulk==T & ancAmo.data$mean.s > -0.1, "mean.s"])

########################################################################################################################
#compare full EMPIRIC data (effects of all mutations on fitness as measured by Mishra, Flynn et al. 2015) to result here
#evidence for same epistatic effect?
hist(EMPIRIC[EMPIRIC$shorthand %in% Sc.data$shorthand,"s"],col="gray70",main="",xlab="Selection coefficient from Mishra et al.",ylab="Number of ancestral reversions")
median(EMPIRIC[EMPIRIC$shorthand %in% Sc.data$shorthand,"s"])
wilcox.test(EMPIRIC[EMPIRIC$shorthand %in% Sc.data$shorthand,"s"],conf.int=T,paired=F)
#same average 1% fitness defect, but not a significant difference (exptl noise)

########################################################################################################################
#analyze epistatic compensations from isogenic growths
#open up and get rel.fitness for remaining isogenic growths, same as before
Sc.iso.rep4 <- read.csv(file="data_in/Sc-reversions_isogenic-growths_rep4.csv",header=T, stringsAsFactors=F)
Sc.iso.rep5 <- read.csv(file="data_in/Sc-reversions_isogenic-growths_rep5.csv",header=T, stringsAsFactors=F)
timepoints4 <- c(17,22,26,30,38.5)
timepoints5 <- c(3,6,11,18,28)

for(i in 1:nrow(Sc.iso.rep4)){
  fit <- lm(as.numeric(Sc.iso.rep4[i,2:6]) ~ timepoints4)
  Sc.iso.rep4$slope[i] <- summary(fit)$coefficients[2,1] #slope of ln(OD) versus time is malthusian parameter
  Sc.iso.rep4$slope.se[i] <- summary(fit)$coefficients[2,2]
}
#calculate dimensionless per-generation selection coefficient, assuming that seleciton is density-independent, as (slope_x - slope_wt)/slope_wt*ln2 (last part equivalent to multiplying by generation time of WT)
Sc.iso.rep4$s <- (Sc.iso.rep4$slope-Sc.iso.rep4$slope[1])/(Sc.iso.rep4$slope[1])*log(2)
#linear transform values such that relative fitness of null allele is 0
#convert raw s to fitness given w=e^s
Sc.iso.rep4$raw.fitness <- exp(Sc.iso.rep4$s)
#stretch relative fitness space to be from 0 (null) to 1 (wt)
Sc.iso.rep4$rel.fitness <- (Sc.iso.rep4$raw.fitness - Sc.iso.rep4$raw.fitness[2])/(Sc.iso.rep4$raw.fitness[1] - Sc.iso.rep4$raw.fitness[2])
#convert back to selection coefficient s=ln(w)
Sc.iso.rep4$s <- log(Sc.iso.rep4$rel.fitness)

for(i in 1:nrow(Sc.iso.rep5)){
  fit <- lm(as.numeric(Sc.iso.rep5[i,2:6]) ~ timepoints5)
  Sc.iso.rep5$slope[i] <- summary(fit)$coefficients[2,1] #slope of ln(OD) versus time is malthusian parameter
  Sc.iso.rep5$slope.se[i] <- summary(fit)$coefficients[2,2]
}
#calculate dimensionless per-generation selection coefficient, assuming that seleciton is density-independent, as (slope_x - slope_wt)/slope_wt*ln2 (last part equivalent to multiplying by generation time of WT)
Sc.iso.rep5$s <- (Sc.iso.rep5$slope-Sc.iso.rep5$slope[1])/(Sc.iso.rep5$slope[1])*log(2)
#linear transform values such that relative fitness of null allele is 0
#convert raw s to fitness given w=e^s
Sc.iso.rep5$raw.fitness <- exp(Sc.iso.rep5$s)
#stretch relative fitness space to be from 0 (null) to 1 (wt)
Sc.iso.rep5$rel.fitness <- (Sc.iso.rep5$raw.fitness - Sc.iso.rep5$raw.fitness[2])/(Sc.iso.rep5$raw.fitness[1] - Sc.iso.rep5$raw.fitness[2])
#convert back to selection coefficient s=ln(w)
Sc.iso.rep5$s <- log(Sc.iso.rep5$rel.fitness)

#V23F/L378I interaction
V23F <- c(Sc.iso.rep4[Sc.iso.rep4$ID=="V23F","rel.fitness"],Sc.iso.rep5[Sc.iso.rep5$ID=="V23F","rel.fitness"])
mean(V23F);sd(V23F)/sqrt(2)

L378I <- c(Sc.iso.rep4[Sc.iso.rep4$ID=="L378I","rel.fitness"],Sc.iso.rep5[Sc.iso.rep5$ID=="L378I","rel.fitness"])
mean(L378I);sd(L378I)/sqrt(2)

V23F.L378I <- Sc.iso.rep4[Sc.iso.rep4$ID=="V23F/L378I","rel.fitness"] #only have one observation

#E7A/T13N/N151A interaction
#single revertants
E7A <- c(Sc.iso.rep1[Sc.iso.rep1$ID=="E7A","rel.fitness"],Sc.iso.rep2[Sc.iso.rep2$ID=="E7A","rel.fitness"],Sc.iso.rep3[Sc.iso.rep3$ID=="E7A","rel.fitness"])
mean(E7A);sd(E7A)/sqrt(3)
T13N <- c(Sc.iso.rep1[Sc.iso.rep1$ID=="T13N","rel.fitness"],Sc.iso.rep2[Sc.iso.rep2$ID=="T13N","rel.fitness"],Sc.iso.rep3[Sc.iso.rep3$ID=="T13N","rel.fitness"])
mean(T13N);sd(T13N)/sqrt(3)
N151A <- c(Sc.iso.rep1[Sc.iso.rep1$ID=="N151A","rel.fitness"],Sc.iso.rep2[Sc.iso.rep2$ID=="N151A","rel.fitness"],Sc.iso.rep3[Sc.iso.rep3$ID=="N151A","rel.fitness"])
mean(N151A);sd(N151A)/sqrt(3)
#double revertants
E7A.T13N <- c(Sc.iso.rep1[Sc.iso.rep1$ID=="E7A/T13N","rel.fitness"],Sc.iso.rep2[Sc.iso.rep2$ID=="E7A/T13N","rel.fitness"],Sc.iso.rep3[Sc.iso.rep3$ID=="E7A/T13N","rel.fitness"])
mean(E7A.T13N);sd(E7A.T13N)/sqrt(3)
E7A.N151A <- c(Sc.iso.rep1[Sc.iso.rep1$ID=="E7A/N151A","rel.fitness"],Sc.iso.rep2[Sc.iso.rep2$ID=="E7A/N151A","rel.fitness"],Sc.iso.rep3[Sc.iso.rep3$ID=="E7A/N151A","rel.fitness"])
mean(E7A.N151A);sd(E7A.N151A)/sqrt(3)
T13N.N151A <- c(Sc.iso.rep3[Sc.iso.rep3$ID=="T13N/N151A","rel.fitness"],Sc.iso.rep5[Sc.iso.rep5$ID=="T13N/N151A","rel.fitness"])
mean(T13N.N151A);sd(T13N.N151A)/sqrt(2)
#triple revertant
E7A.T13N.N151A <- c(Sc.iso.rep1[Sc.iso.rep1$ID=="E7A/T13N/N151A","rel.fitness"],Sc.iso.rep2[Sc.iso.rep2$ID=="E7A/T13N/N151A","rel.fitness"],Sc.iso.rep3[Sc.iso.rep3$ID=="E7A/T13N/N151A","rel.fitness"])
mean(E7A.T13N.N151A);sd(E7A.T13N.N151A)/sqrt(3)

t.test(E7A,E7A.T13N)
t.test(E7A,E7A.N151A)
t.test(E7A,E7A.T13N.N151A)
t.test(E7A.T13N,E7A.T13N.N151A)
t.test(T13N.N151A,E7A.T13N.N151A)
t.test(E7A.N151A,E7A.T13N.N151A)

#test whether states that interact with E7A interact with V23F, and vice-versa
V23F.T13N <- Sc.iso.rep5[Sc.iso.rep5$ID=="V23F/T13N","rel.fitness"]
V23F.N151A <- Sc.iso.rep5[Sc.iso.rep5$ID=="V23F/N151A","rel.fitness"]
V23F.T13N.N151A <- Sc.iso.rep5[Sc.iso.rep5$ID=="V23F/T13N/N151A","rel.fitness"]
E7A.L378I <- Sc.iso.rep5[Sc.iso.rep5$ID=="E7A/L378I","rel.fitness"]

#plot: for V23F/L378 interaction, cycle of mutants rel fitness to Sc wildtype
pdf(file="plots/V23F.L378I.pdf",2.5,3,useDingbats = F)
plot(1,log(1),xlim=c(0.5,3.5),ylim=c(log(0.67),log(1.03)),pch=20,ylab="log-fitness relative to ScHsp82")
arrows(2,mean(log(V23F)) - sd(log(V23F))/sqrt(length(V23F)), 2, mean(log(V23F)) + sd(log(V23F))/sqrt(length(V23F)),length=0.025,angle=90,code=3,col="black")
points(2,mean(log(V23F)), pch=20,col="black")
arrows(2,mean(log(L378I)) - sd(log(L378I))/sqrt(length(L378I)), 2,mean(log(L378I)) + sd(log(L378I))/sqrt(length(L378I)),length=0.025,angle=90,code=3,col="black")
points(2,mean(log(L378I)),pch=20,col="black")
points(3,mean(log(V23F.L378I)),pch=20,col="black")
dev.off()

#plot: for E7A/T13N/N151A interaction, cycle of mutants rel fitness to Sc wildtype
pdf(file="plots/E7A.T13N.N151A.pdf",3,3,useDingbats = F)
plot(1,log(1),xlim=c(0.5,4.5),ylim=c(log(0.7),log(1.03)),pch=20,ylab="log-fitness relative to ScHsp82")
arrows(2,mean(log(E7A)) - sd(log(E7A))/sqrt(length(E7A)), 2, mean(log(E7A)) + sd(log(E7A))/sqrt(length(E7A)),length=0.025,angle=90,code=3,col="black")
points(2,mean(log(E7A)), pch=20,col="black")
arrows(2,mean(log(T13N)) - sd(log(T13N))/sqrt(length(T13N)), 2,mean(log(T13N)) + sd(log(T13N))/sqrt(length(T13N)),length=0.025,angle=90,code=3,col="black")
points(2,mean(log(T13N)),pch=20,col="black")
arrows(2,mean(log(N151A)) - sd(log(N151A))/sqrt(length(N151A)), 2,mean(log(N151A)) + sd(log(N151A))/sqrt(length(N151A)),length=0.025,angle=90,code=3,col="black")
points(2,mean(log(N151A)),pch=20,col="black")
arrows(3,mean(log(E7A.T13N)) - sd(log(E7A.T13N))/sqrt(length(E7A.T13N)), 3, mean(log(E7A.T13N)) + sd(log(E7A.T13N))/sqrt(length(E7A.T13N)),length=0.025,angle=90,code=3,col="black")
points(3,mean(log(E7A.T13N)), pch=20,col="black")
arrows(3,mean(log(E7A.N151A)) - sd(log(E7A.N151A))/sqrt(length(E7A.N151A)), 3, mean(log(E7A.N151A)) + sd(log(E7A.N151A))/sqrt(length(E7A.N151A)),length=0.025,angle=90,code=3,col="black")
points(3,mean(log(E7A.N151A)), pch=20,col="black")
arrows(3,mean(log(T13N.N151A)) - sd(log(T13N.N151A))/sqrt(length(T13N.N151A)), 3, mean(log(T13N.N151A)) + sd(log(T13N.N151A))/sqrt(length(T13N.N151A)),length=0.025,angle=90,code=3,col="black")
points(3,mean(log(T13N.N151A)), pch=20,col="black")
arrows(4,mean(log(E7A.T13N.N151A)) - sd(log(E7A.T13N.N151A))/sqrt(length(E7A.T13N.N151A)), 4, mean(log(E7A.T13N.N151A)) + sd(log(E7A.T13N.N151A))/sqrt(length(E7A.T13N.N151A)),length=0.025,angle=90,code=3,col="black")
points(4,mean(log(E7A.T13N.N151A)), pch=20,col="black")
dev.off()

#plot for E7A/L378 (lack of) interaction
pdf(file="plots/E7A.L378I.pdf",2.5,3,useDingbats = F)
plot(1,log(1),xlim=c(0.5,3.5),ylim=c(log(0.7),log(1.05)),pch=20,ylab="log-fitness relative to ScHsp82")
arrows(2,mean(log(E7A)) - sd(log(E7A))/sqrt(length(E7A)), 2, mean(log(E7A)) + sd(log(E7A))/sqrt(length(E7A)),length=0.025,angle=90,code=3,col="black")
points(2,mean(log(E7A)), pch=20,col="black")
arrows(2,mean(log(L378I)) - sd(log(L378I))/sqrt(length(L378I)), 2,mean(log(L378I)) + sd(log(L378I))/sqrt(length(L378I)),length=0.025,angle=90,code=3,col="black")
points(2,mean(log(L378I)),pch=20,col="black")
points(3,mean(log(E7A.L378I)),pch=20,col="black")
dev.off()

#plot for V23F/T13N/N151A (lack of) interaction
pdf(file="plots/V23F.T13N.N151A.pdf",3,3,useDingbats = F)
plot(1,log(1),xlim=c(0.5,4.5),ylim=c(log(0.65),log(1.03)),pch=20,ylab="log-fitness relative to ScHsp82")
arrows(2,mean(log(V23F)) - sd(log(V23F))/sqrt(length(V23F)), 2, mean(log(V23F)) + sd(log(V23F))/sqrt(length(V23F)),length=0.025,angle=90,code=3,col="black")
points(2,mean(log(V23F)), pch=20,col="black")
arrows(2,mean(log(T13N)) - sd(log(T13N))/sqrt(length(T13N)), 2,mean(log(T13N)) + sd(log(T13N))/sqrt(length(T13N)),length=0.025,angle=90,code=3,col="black")
points(2,mean(log(T13N)),pch=20,col="black")
arrows(2,mean(log(N151A)) - sd(log(N151A))/sqrt(length(N151A)), 2,mean(log(N151A)) + sd(log(N151A))/sqrt(length(N151A)),length=0.025,angle=90,code=3,col="black")
points(2,mean(log(N151A)),pch=20,col="black")
arrows(3,mean(log(V23F.T13N)) - sd(log(V23F.T13N))/sqrt(length(V23F.T13N)), 3, mean(log(V23F.T13N)) + sd(log(V23F.T13N))/sqrt(length(V23F.T13N)),length=0.025,angle=90,code=3,col="black")
points(3,mean(log(V23F.T13N)), pch=20,col="black")
arrows(3,mean(log(V23F.N151A)) - sd(log(V23F.N151A))/sqrt(length(V23F.N151A)), 3, mean(log(V23F.N151A)) + sd(log(V23F.N151A))/sqrt(length(V23F.N151A)),length=0.025,angle=90,code=3,col="black")
points(3,mean(log(V23F.N151A)), pch=20,col="black")
arrows(3,mean(log(T13N.N151A)) - sd(log(T13N.N151A))/sqrt(length(T13N.N151A)), 3, mean(log(T13N.N151A)) + sd(log(T13N.N151A))/sqrt(length(T13N.N151A)),length=0.025,angle=90,code=3,col="black")
points(3,mean(log(T13N.N151A)), pch=20,col="black")
arrows(4,mean(log(V23F.T13N.N151A)) - sd(log(V23F.T13N.N151A))/sqrt(length(V23F.T13N.N151A)), 4, mean(log(V23F.T13N.N151A)) + sd(log(V23F.T13N.N151A))/sqrt(length(V23F.T13N.N151A)),length=0.025,angle=90,code=3,col="black")
points(4,mean(log(V23F.T13N.N151A)), pch=20,col="black")
dev.off()

#correlation in fitness as determined in monoculture growths or bulk competitions for the six genotypes where we have measurements both ways
iso <- c(log(mean(E7A)),log(mean(L378I)),log(mean(N151A)),log(mean(T13N)),log(mean(V23F)),ancAsco.iso[4,"s"],ancAsco.iso[5,"s"])
bulk <- c(Sc.data[Sc.data$shorthand=="E7A","mean.s"],NA,Sc.data[Sc.data$shorthand=="N151A","mean.s"],Sc.data[Sc.data$shorthand=="T13N","mean.s"],Sc.data[Sc.data$shorthand=="V23F","mean.s"],log(ancAmo.mean.w),log(ancAmo.i378.mean.w))

plot(bulk,iso,pch=19,xlim=c(-1,0.05),ylim=c(-1,0.05),xlab="Selection coefficient, bulk competition", ylab="Selection coefficient, monoculture growth");abline(lm(iso~bulk),lty=2);summary(lm(iso~bulk))

#FIN