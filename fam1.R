# Processing the raw data of the MDC 7984 Family; Adjusting the data
# last modified 03.08.2021
# created 31.01.2020
# Florian Krause & Dr. Danny Arends

### Set the working directory ###
setwd("")


# Read the merged input data
genopheno <- read.table("fam1input.txt", sep = " ", header = TRUE, na.strings=c(""), colClasses = c("character"))

####Remove rows with NAs in the "Genotype" column

new_genopheno <- genopheno[- grep("NA", genopheno$Genotype),]
genopheno <- new_genopheno

### Read in the pheno MRI inputfile as phenomri ###

phenomri<- read.csv("fam1mri.txt", sep = "\t", header = TRUE, na.strings=c(""), colClasses = c("character"))

### Set the rownames of the animal ID's: paste "MDC-7984-" ###
rownames(phenomri) <- paste0('MDC-', phenomri$ID.Nr)

### Change rownames ###
rownames(phenomri) <- gsub("MDC-", "", rownames(phenomri))

genopheno <- cbind(genopheno, phenomri[rownames(genopheno),])

#genopheno

### remove NAs in MRI phenotypes ###

toremain <- which(is.na(genopheno[,"fat.lean70"])) 
genopheno <- genopheno[-toremain,]  # ToRem <- which(is.na(x[,"gt"]))

# Adjusting

genopheno[,"fat.lean70"] <- gsub("," , ".", genopheno[,"fat.lean70"])
genopheno[,"fat.lean70"] <- as.numeric(genopheno[,"fat.lean70"])
#genopheno[,"X70"] <- as.numeric(genopheno$X70)
anova(lm(genopheno[, "fat.lean70"]~genopheno[, "sex"])) #  Significant

genopheno$WG <- as.numeric(genopheno$WG)
anova(lm(genopheno[, "fat.lean70"]~genopheno[, "WG"]))

genopheno$WG <- as.numeric(genopheno$GenGood)
anova(lm(genopheno[, "fat.lean70"]~genopheno[, "GenGood"])) #  Significant

# Contains the residuals to add back in the mean of the timepoint

toremove <- which(is.na(genopheno[, "fat.lean70"]))
if(length(toremove)>0) genopheno <- genopheno[-toremove,]
fit <- lm(genopheno[, "fat.lean70"] ~ genopheno[, "sex"] + genopheno[, "GenGood"])
d70 <- resid(fit) + mean(genopheno[, "fat.lean70"])

nwo <- ordered(genopheno$Genotype, levels = c("BFMI/BFMI", "BFMI/MDC", "MDC/MDC", "MDC/B6N", "BFMI/B6N", "B6N/B6N"))

#################################################t-test#########################################
x <- d70
y <- nwo
x_name <- "FAT.LEAN.DAY.70"
y_name <- "ORDERED.GENOTYPES"

df <- data.frame(x,y)
names(df) <- c(x_name,y_name)

#BFMI/MDC vs MDC/B6N
BFMIMDC <- subset(df, ORDERED.GENOTYPES == "BFMI/MDC", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
MDCB6N <- subset(df, ORDERED.GENOTYPES == "MDC/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIMDCvsMDCB6N <- rbind(BFMIMDC, MDCB6N)
t.test(BFMIMDCvsMDCB6N$FAT.LEAN.DAY.70~BFMIMDCvsMDCB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#BFMI/MDC vs MDC/MDC
BFMIMDC <- subset(df, ORDERED.GENOTYPES == "BFMI/MDC", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
MDCMDC <- subset(df, ORDERED.GENOTYPES == "MDC/MDC", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIMDCvsMDCMDC <- rbind(BFMIMDC, MDCMDC)
t.test(BFMIMDCvsMDCMDC$FAT.LEAN.DAY.70~BFMIMDCvsMDCMDC$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#BFMI/MDC vs. BFMI/B6N
BFMIMDC <- subset(df, ORDERED.GENOTYPES == "BFMI/MDC", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIB6N <- subset(df, ORDERED.GENOTYPES == "BFMI/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIMDCvsBFMIB6N <- rbind(BFMIMDC, BFMIB6N)
t.test(BFMIMDCvsBFMIB6N$FAT.LEAN.DAY.70~BFMIMDCvsBFMIB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#BFMI/BFMI vs BFMI/B6N
BFMIBFMI <- subset(df, ORDERED.GENOTYPES == "BFMI/BFMI", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIB6N <- subset(df, ORDERED.GENOTYPES == "BFMI/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIBFMIvsBFMIB6N <- rbind(BFMIBFMI, BFMIB6N)
t.test(BFMIBFMIvsBFMIB6N$FAT.LEAN.DAY.70~BFMIBFMIvsBFMIB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="greater")

#BFMI/BFMI vs. MDC/MDC
BFMIBFMI <- subset(df, ORDERED.GENOTYPES == "BFMI/BFMI", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
MDCMDC <- subset(df, ORDERED.GENOTYPES == "MDC/MDC", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIBFMIvsMDCMDC <- rbind(BFMIBFMI, MDCMDC)
t.test(BFMIBFMIvsMDCMDC$FAT.LEAN.DAY.70~BFMIBFMIvsMDCMDC$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#MDC/MDC vs. B6N/B6N
MDCMDC <- subset(df, ORDERED.GENOTYPES == "MDC/MDC", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
B6NB6N <- subset(df, ORDERED.GENOTYPES == "B6N/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
MDCMDCvsB6NB6N <- rbind(MDCMDC, B6NB6N)
t.test(MDCMDCvsB6NB6N$FAT.LEAN.DAY.70~MDCMDCvsB6NB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#BFMI/BFMI vs B6N/B6N
BFMIBFMI <- subset(df, ORDERED.GENOTYPES == "BFMI/BFMI", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
B6NB6N <- subset(df, ORDERED.GENOTYPES == "B6N/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIBFMIvsB6NB6N <- rbind(BFMIBFMI, B6NB6N)
t.test(BFMIBFMIvsB6NB6N$FAT.LEAN.DAY.70~BFMIBFMIvsB6NB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="greater")

#MDC/B6N vs. BFMI/B6N
MDCB6N <- subset(df, ORDERED.GENOTYPES == "MDC/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIB6N <- subset(df, ORDERED.GENOTYPES == "BFMI/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
MDCB6NvsBFMIB6N <- rbind(MDCB6N, BFMIB6N)
t.test(MDCB6NvsBFMIB6N$FAT.LEAN.DAY.70~MDCB6NvsBFMIB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#MDC/MDC vs. MDC/B6N

MDCMDC <- subset(df, ORDERED.GENOTYPES == "MDC/MDC", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
MDCB6N <- subset(df, ORDERED.GENOTYPES == "MDC/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
MDCMDCvsMDCB6N <- rbind(MDCMDC, BFMIB6N)
t.test(MDCMDCvsMDCB6N$FAT.LEAN.DAY.70~MDCMDCvsMDCB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#BFMI/BFMI vs BFMI/MDC
BFMIBFMI <- subset(df, ORDERED.GENOTYPES == "BFMI/BFMI", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIMDC <- subset(df, ORDERED.GENOTYPES == "BFMI/MDC", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIBFMIvsBFMIMDC <- rbind(BFMIBFMI, BFMIMDC)
t.test(BFMIBFMIvsBFMIMDC$FAT.LEAN.DAY.70~BFMIBFMIvsBFMIMDC$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#BFMI/BFMI vs MDC/B6N
BFMIBFMI <- subset(df, ORDERED.GENOTYPES == "BFMI/BFMI", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
MDCB6N <- subset(df, ORDERED.GENOTYPES == "MDC/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIBFMIvsMDCB6N <- rbind(BFMIBFMI, MDCB6N)
t.test(BFMIBFMIvsMDCB6N$FAT.LEAN.DAY.70~BFMIBFMIvsMDCB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#BFMI/MDC vs B6N/B6N
BFMIMDC <- subset(df, ORDERED.GENOTYPES == "BFMI/MDC", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
B6NB6N <- subset(df, ORDERED.GENOTYPES == "B6N/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIMDCvsB6NB6N <- rbind(BFMIMDC, B6NB6N)
t.test(BFMIMDCvsB6NB6N$FAT.LEAN.DAY.70~BFMIMDCvsB6NB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#BFMI/B6N vs MDC/B6N
BFMIB6N <- subset(df, ORDERED.GENOTYPES == "BFMI/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
MDCB6N <- subset(df, ORDERED.GENOTYPES == "MDC/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIB6NvsMDCB6N <- rbind(BFMIB6N, MDCB6N)
t.test(BFMIB6NvsMDCB6N$FAT.LEAN.DAY.70~BFMIB6NvsMDCB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#BFMI/B6N vs B6N/B6N
BFMIB6N <- subset(df, ORDERED.GENOTYPES == "BFMI/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
B6NB6N <- subset(df, ORDERED.GENOTYPES == "B6N/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIB6NvsB6NB6N <- rbind(BFMIB6N, B6NB6N)
t.test(BFMIB6NvsB6NB6N$FAT.LEAN.DAY.70~BFMIB6NvsB6NB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#MDC/B6N vs B6N/B6N
MDCB6N <- subset(df, ORDERED.GENOTYPES == "MDC/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
B6NB6N <- subset(df, ORDERED.GENOTYPES == "B6N/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
MDCB6NvsB6NB6N <- rbind(MDCB6N, B6NB6N)
t.test(MDCB6NvsB6NB6N$FAT.LEAN.DAY.70~MDCB6NvsB6NB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#MDC/B6N vs BFMI/B6N
MDCB6N <- subset(df, ORDERED.GENOTYPES == "MDC/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIB6N <- subset(df, ORDERED.GENOTYPES == "BFMI/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
MDCB6NvsBFMIB6N <- rbind(MDCB6N, BFMIB6N)
t.test(MDCB6NvsBFMIB6N$FAT.LEAN.DAY.70~MDCB6NvsBFMIB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#BFMI/B6N vs MDC/MDC
BFMIB6N <- subset(df, ORDERED.GENOTYPES == "BFMI/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
MDCMDC <- subset(df, ORDERED.GENOTYPES == "MDC/MDC", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIB6NvsMDCMDC <- rbind(BFMIB6N, MDCMDC)
t.test(BFMIB6NvsMDCMDC$FAT.LEAN.DAY.70~BFMIB6NvsMDCMDC$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

mean(BFMIBFMI$FAT.LEAN.DAY.70)
mean(BFMIMDC$FAT.LEAN.DAY.70)
mean(MDCMDC$FAT.LEAN.DAY.70)
mean(BFMIB6N$FAT.LEAN.DAY.70)
mean(MDCB6N$FAT.LEAN.DAY.70)
mean(B6NB6N$FAT.LEAN.DAY.70)

sd(BFMIBFMI$FAT.LEAN.DAY.70)
sd(BFMIMDC$FAT.LEAN.DAY.70)
sd(MDCMDC$FAT.LEAN.DAY.70)
sd(BFMIB6N$FAT.LEAN.DAY.70)
sd(MDCB6N$FAT.LEAN.DAY.70)
sd(B6NB6N$FAT.LEAN.DAY.70)


#################################################t-test################################################

###Renaming MDC to final Name ###

lkup <- NA
lkup[df[, "ORDERED.GENOTYPES"]=="MDC/B6N"] <- "I8∆1/B6N"
lkup[df[, "ORDERED.GENOTYPES"]=="MDC/MDC"] <- "I8∆1/I8∆1"
lkup[df[, "ORDERED.GENOTYPES"]=="BFMI/MDC"] <- "BFMI/I8∆1"
lkup[df[, "ORDERED.GENOTYPES"]=="BFMI/BFMI"] <- "BFMI/BFMI"
lkup[df[, "ORDERED.GENOTYPES"]=="BFMI/B6N"] <- "BFMI/B6N"
lkup[df[, "ORDERED.GENOTYPES"]=="B6N/B6N"] <- "B6N/B6N"

### Add the lkup as "Genotype" vector to geno

#as.data.frame(lkup)
df = cbind(df, ORDERED.GENOTYPES = lkup)

df <- df[, -2]

nwo <- ordered(df$ORDERED.GENOTYPES, levels = c("BFMI/BFMI", "BFMI/I8∆1", "I8∆1/I8∆1", "I8∆1/B6N", "BFMI/B6N", "B6N/B6N"))

myplot <- boxplot(df$FAT.LEAN.DAY.70~nwo,



main = "I8∆1 Family",
xlab = "Genotypes",
ylab = "Adjusted fat-to-lean ratio at week 10", 
yaxt="none",
main="Turn off y-axis",
ylim = c(0, 0.6)


)

myplot

myplot$n
text(1:length(myplot$n),rep(0,length(myplot$n)),myplot$n)

axis(2, las=2)

write.table(genopheno, file = "final_table_7984.txt")

