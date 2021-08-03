# t-Test for both Families
# last modified 03.08.2021
# created 27.07.2021
# Florian Krause

### Set the working directory ###
setwd("")


samples <- read.table("both_families_ttest.csv", sep = ",", header = TRUE, na.strings=c(""), colClasses = c("character"))

samples$fatlean <- as.numeric(samples$fatlean)
fatlean <- samples$fatlean
animals <- samples$animals


t.test(samples$fatlean~samples$animals, var.equal = TRUE, alternative = "two.sided")

