# Processing the raw data of the MDC 8015 Family; Adjusting the data
# last modified 03.08.2021
# created 31.01.2020
# Florian Krause & Dr. Danny Arends

### Set the working directory ###
setwd("")

# Read the merged input data
genopheno <- read.table("fam2input_groups.txt", sep = " ", header = TRUE, na.strings=c(""), colClasses = c("character"))

####Remove rows with NAs in the "Genotype" column

new_genopheno <- genopheno[- grep("NA", genopheno$Genotype),]
genopheno <- new_genopheno

### Read in the pheno MRI inputfile as phenomri ###

phenomri<- read.csv("fam2mri.txt", sep = "\t", header = TRUE, na.strings=c(""), colClasses = c("character"))

### Set the rownames of the animal ID's: paste "MDC-8015-" ###
rownames(phenomri) <- paste0('MDC-', phenomri$ID.Nr)

### Change rownames ###
rownames(phenomri) <- gsub("MDC-", "", rownames(phenomri))

#genopheno <- cbind(genopheno, phenomri[rownames(phenomri),])
genopheno <- cbind(genopheno, phenomri[rownames(genopheno),])

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

nwo <- ordered(genopheno$Genotype, levels = c("BFMI/BFMI", "y/MDC", "x/B6N"))

#################################################t-test#########################################
x <- d70
y <- nwo
x_name <- "FAT.LEAN.DAY.70"
y_name <- "ORDERED.GENOTYPES"

df <- data.frame(x,y)
names(df) <- c(x_name,y_name)

#y/MDC vs x/B6N
yMDC <- subset(df, ORDERED.GENOTYPES == "y/MDC", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
xB6N <- subset(df, ORDERED.GENOTYPES == "x/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
yMDCvsxB6N <- rbind(yMDC, xB6N)
t.test(yMDCvsxB6N$FAT.LEAN.DAY.70~yMDCvsxB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")

#BFMI/BFMI vs x/B6N
BFMIBFMI <- subset(df, ORDERED.GENOTYPES == "BFMI/BFMI", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
xB6N <- subset(df, ORDERED.GENOTYPES == "x/B6N", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIBFMIvsxB6N <- rbind(BFMIBFMI, xB6N)
t.test(BFMIBFMIvsxB6N$FAT.LEAN.DAY.70~BFMIBFMIvsxB6N$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="greater")

#BFMI/BFMI vs. y/MDC
BFMIBFMI <- subset(df, ORDERED.GENOTYPES == "BFMI/BFMI", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
yMDC <- subset(df, ORDERED.GENOTYPES == "y/MDC", select = c("ORDERED.GENOTYPES", "FAT.LEAN.DAY.70"))
BFMIBFMIvsyMDC <- rbind(BFMIBFMI, yMDC)
t.test(BFMIBFMIvsyMDC$FAT.LEAN.DAY.70~BFMIBFMIvsyMDC$ORDERED.GENOTYPES, var.equal = FALSE, alternative ="two.sided")


mean(yMDC$FAT.LEAN.DAY.70)
mean(xB6N$FAT.LEAN.DAY.70)

sd(yMDC$FAT.LEAN.DAY.70)
sd(xB6N$FAT.LEAN.DAY.70)

#################################################t-test################################################

###Renaming MDC to final Name ###

lkup <- NA
lkup[df[, "ORDERED.GENOTYPES"]=="x/B6N"] <- "x/B6N"
lkup[df[, "ORDERED.GENOTYPES"]=="y/MDC"] <- "y/I8∆2"
lkup[df[, "ORDERED.GENOTYPES"]=="y/MDC"] <- "y/I8∆2"
lkup[df[, "ORDERED.GENOTYPES"]=="BFMI/BFMI"] <- "BFMI/BFMI"
lkup[df[, "ORDERED.GENOTYPES"]=="x/B6N"] <- "x/B6N"
lkup[df[, "ORDERED.GENOTYPES"]=="x/B6N"] <- "x/B6N"

### Add the lkup as "Genotype" vector to geno

#as.data.frame(lkup)
df = cbind(df, ORDERED.GENOTYPES = lkup)

df <- df[, -2]

nwo <- ordered(df$ORDERED.GENOTYPES, levels = c("BFMI/BFMI", "y/I8∆2", "x/B6N"))

myplot <- boxplot(df$FAT.LEAN.DAY.70~nwo,

main = "I8∆2 Family",
xlab = "Genotype groups",
ylab = "Adjusted fat-to-lean ratio at week 10", 
yaxt="none",
main="Turn off y-axis",
ylim = c(0, 0.6)


)

myplot

myplot$n
text(1:length(myplot$n),rep(0,length(myplot$n)),myplot$n)

axis(2, las=2)
