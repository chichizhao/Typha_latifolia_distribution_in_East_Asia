#!/Rscript
# RDA analysis
# set the working directory
# reference: https://popgen.nescent.org/2018-03-27_RDA_GEA.html
setwd("/home/chichi/data/china/china3/pca/")
# Load the data
# the samples names
samples <- read.table("samples.txt", header = FALSE)
library(data.table)
genotype_data <- fread("genetic_matrix.txt", header = FALSE)
# generate the genetic_data
# remove columns 1 and 2 of the genotype_data
genotype_data <- genotype_data[, -c(1, 2)]

genotype_data[genotype_data == "0/0"] <- 0
genotype_data[genotype_data == "0|0"] <- 0
genotype_data[genotype_data == "./."] <- 0
genotype_data[genotype_data == ".|."] <- 0
genotype_data[genotype_data == "0/1"] <- 1
genotype_data[genotype_data == "0|1"] <- 1
genotype_data[genotype_data == "1/1"] <- 2
genotype_data[genotype_data == "1|1"] <- 2

genetic_matrix_t <- t(genotype_data)
# count the rows of the genetic_matrix_t
rows <- nrow(genetic_matrix_t)
# give the genetic_matrix_t the samples names index
rownames(genetic_matrix_t) <- samples$V1[1:rows]

write.csv(genetic_matrix_t, "genetic_matrix.csv", row.names = TRUE)
# read the genetic_matrix_t
genetic_matrix_t <- read.csv("genetic_matrix.csv", header = TRUE, row.names = 1)


library(vegan)
library(psych)

#bio_factors <- read.csv("bio_factors.csv")
bio_factors <- read.table("/home/chichi/data/china/china3/future/lfmm/env_data4.csv", header=T, sep=",", row.names=1)
# confirm the row names of the bio_factors are the same as the genetic_matrix_t
identical(rownames(bio_factors), rownames(genetic_matrix_t))
# read the samples classification
samples_classification <- read.table("pop3.txt", header = FALSE)
# rename the column names
colnames(samples_classification) <- c("site", "pop")
# add the samples classification to the bio_factors, matching index and samples names
sorted_samples_classification <- samples_classification[match(rownames(bio_factors), samples_classification$site), ]
# add the samples classification to the bio_factors
bio_factors <- cbind(bio_factors, sorted_samples_classification)
# rename the column name
# colnames(bio_factors)[ncol(sorted_samples_classification$V2)] <- "pop"

# convert the pop column to factor
bio_factors$pop <- as.factor(bio_factors$pop)
# set the pop column into levels
bio_factors$pop <- factor(bio_factors$pop, levels = c("WA", "WB","WC","BA","EC","EB","EA","MA"))
# reorder the bio_factors to match the genetic_matrix_t
bio_factors <- bio_factors[rownames(genetic_matrix_t), ]
identical(rownames(bio_factors), rownames(genetic_matrix_t))
pairs.panels(bio_factors, scale=T)
pdf("bio_factors.pdf", width=12, height=12)
pairs.panels(bio_factors, scale=T)
dev.off()
# here the pattern are strong similar are those 
# bio19, bio17, bio14
# bio16, bio18, bio13, bio12
# bio11, bio9, bio6
# bio8, bio10, bio5
# bio7, bio4
bio_factors_selected <- subset(bio_factors, select=-c(bio_17, bio_14, bio_16,bio_18, bio_13, bio_9, bio_6,bio_8, bio_5, bio_7,bio_1))

pdf("bio_factors_selected.pdf", width=12, height=12)
pairs.panels(bio_factors_selected, scale=T)
dev.off()

bio_factors2 <- read.table("/home/chichi/Typha_latifolia_in_China/script/PCA_East_Asia_climate_TOP_3_8pops.csv", header=T, sep=",")
pairs.panels(bio_factors2, scale=T)
pdf("bio_factors2.pdf", width=12, height=12)
pairs.panels(bio_factors2, scale=T)

dev.off()
bio_factors2_selected <- subset(bio_factors2, select=-c(bio17,bio14,bio16,bio13, bio12,bio11,bio9,bio7,bio5,bio8,bio1))
pdf("bio_factors2_selected.pdf", width=12, height=12)
pairs.panels(bio_factors2_selected, scale=T)
dev.off()
# select the same samples the row 3 to 108 of the genetic_matrix.csv
genetic_matrix <- genetic_matrix_t
# Perform RDA
str(genetic_matrix)
genetic_matrix <- apply(genetic_matrix, 2, as.numeric)
identical(rownames(bio_factors), rownames(genetic_matrix_t))

# remove the columns with na NaN null
genetic_matrix <- na.omit(genetic_matrix)

str(genetic_matrix)
summary(genetic_matrix)
non_numeric_indices <- which(is.na(apply(genetic_matrix, 2, as.numeric)), arr.ind = TRUE)
non_numeric_values <- genetic_matrix[non_numeric_indices]
print(non_numeric_values)
# Option 1: Remove rows with non-numeric values
genetic_matrix <- genetic_matrix[complete.cases(apply(genetic_matrix, 2, as.numeric)), ]

# Option 2: Replace non-numeric values with a specific value (e.g., 0)
genetic_matrix[is.na(apply(genetic_matrix, 2, as.numeric))] <- 0

bio_factors_selected <- na.omit(bio_factors_selected)
str(bio_factors_selected)
rows <- nrow(bio_factors_selected)
rda_result <- rda(genetic_matrix ~ ., data = bio_factors_selected, scale = TRUE)

RsquareAdj(rda_result)

# Summarize and interpret the results
summary(rda_result)
# save the results
save(rda_result, file = "rda_result.RData")

screeplot(rda_result)

pdf("rda_screeplot.pdf", width=6, height=6)
# only plot the first 5, mark each rda axis use bars # set the ylim to 500000
screeplot(rda_result, npcs=5, type="bar", main="RDA Results", ylab="Variance")
dev.off()




plot(rda_result, scaling = 3)
pdf("rda_plot.pdf", width=12, height=12)
plot(rda_result, scaling = 3)
dev.off()
png("rda_plot.png", width=5, height=5, units="in", res=300)
plot(rda_result, scaling = 3,
     xlim = c(-22, 19), ylim = c(-16, 9))
dev.off()
png("rda_plot2.png", width=5, height=5, units="in", res=300)
plot(rda_result,choices = c(1, 3), scaling=3)
dev.off() 


load.rda <- scores(rda_result, choices=c(1:3), display="species")
png("rda_plot3.png", width=5, height=5, units="in", res=300)
hist(load.rda[,1], main="Loadings on RDA1")
dev.off()

# read rda_result
load("rda_result.RData")
library(parallel)
num_cores <- detectCores()
# check the significanct of the RDA
options(mc.cores = num_cores/2)
signif.full <- anova.cca(rda_result, parallel = options("mc.cores"))

# save the results
save(signif.full, file = "signif_full.RData")
signif.axis <- anova.cca(rda_result, by = "axis", parallel = options("mc.cores"))

pop <- bio_factors_selected$pop
# change the bio name, remove the _ in the name column
library(dplyr)
bio_factors_selected <- bio_factors_selected %>% rename_with(~ gsub("_", "", .x))
# also remove the _ in rda_result,change the bp column


#bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c","#fb9a99","#fdbf6f")
# colors = { 'WA':'#F3740B', 'WB':'#EFA38A', 'WC':'#D2352C','BA':'#808080','EC':'#ADD9ED', 'EB':'#30BBCE','EA':'#4A7DB4','MA':'#68AC56'}
bg <- c("#F3740B", "#EFA38A", "#D2352C", "#808080", "#ADD9ED", "#30BBCE", "#4A7DB4", "#68AC56")
png("rda_plot4.png", width=5, height=5, units="in", res=300)
#pdf("rda_plot4.pdf", width=5, height=5)
par(mar = c(4, 4, 1, 1))
plot(rda_result, scaling = 3, type="n",xlim = c(-24, 23), ylim = c(-16, 8), xlab="RDA1 ***", ylab="RDA2 ***")

points(rda_result, display="species", pch=20, cex=1, col="gray32", scaling=3)           # the SNPs
points(rda_result, display="sites", pch=21, cex=1, col="gray32", scaling=3, bg=bg[pop]) # the wolves
# Extract the labels for the "bp" display
bp_labels <- labels(rda_result, display = "bp", scaling = 3)

# Replace underscores with spaces in the labels
bp_labels <- gsub("_", "", bp_labels)

# Plot the modified labels
text(rda_result, scaling = 3, display = "bp", col = "#0868ac", cex = 1, labels = bp_labels)
#text(rda_result, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topright", legend=levels(pop), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

dev.off()

# rda 1 and 3 
png("rda_plot5.png", width=5, height=5, units="in", res=300)

par(mar = c(4, 4, 1, 1))
plot(rda_result,choices = c(1, 3),xlim = c(-26, 20), ylim = c(-15, 7),scaling=3, xlab="RDA1 ***", ylab="RDA3 .")

points(rda_result,choices = c(1, 3), display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(rda_result,choices = c(1, 3), display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[pop]) # the wolves

text(rda_result,choices = c(1, 3), scaling = 3, display = "bp", col = "#0868ac", cex = 1)
legend("bottomleft", legend=levels(pop), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

dev.off()