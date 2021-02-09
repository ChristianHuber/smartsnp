# Simulate data

# Data parameters
cellN <- 1e06
sampleN <- 100
snpN <- cellN/sampleN

# Random sample of genotypes /0|1|2/ and missing values /9/
sampleSNP <- sample(x = c(0,1,2,9), size = cellN, replace = TRUE, prob=c(0.33,0.33,0.33,0.01))

# Build dataset
dataSNP <- matrix(sampleSNP, nrow = snpN, ncol = sampleN)
dim(dataSNP)

# Save data to working directory
write.table(dataSNP, file = "data-raw/dataSNP", col.names = FALSE, row.names = FALSE)

# Make data available to package
library(usethis)
usethis::use_data(dataSNP)
