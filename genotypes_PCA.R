library(VariantAnnotation)
library(ggfortify)
library(vcfR)
library(usethis)
library(ggfortify)
edit_r_environ()

vcf_file <- "ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
param <- ScanVcfParam(info = c("AF"), geno=NA)

data <- read.vcfR(vcf_file)
head(getFIX(data))

gt <- extract.gt(data, element = "GT")

filtered_genotypes_df <- as.data.frame(gt)

View(filtered_genotypes_df)

for (row in seq(1, nrow(filtered_genotypes_df), by = 1000)) {
  for (col in 1:ncol(filtered_genotypes_df)) {
    if (filtered_genotypes_df[row,col]=="0|0") {
      filtered_genotypes_df[row,col] = 0
    } else if (filtered_genotypes_df[row,col] == "0|1" | filtered_genotypes_df[row,col] == "1|0") {
      filtered_genotypes_df[row,col] = 1
    } else if (filtered_genotypes_df[row,col]=="1|1"){
      filtered_genotypes_df[row,col] = 2
    }
  }
} 

subset_df <- filtered_genotypes_df[seq(1, nrow(filtered_genotypes_df), by = 1000),]
View(alt_alleles_df)

alt_alleles_df <- matrix(as.numeric(unlist(subset_df)), 
                         ncol = ncol(subset_df))

rownames(alt_alleles_df) <- rownames(subset_df)
colnames(alt_alleles_df) <- colnames(subset_df)

pc <- prcomp(na.omit(alt_alleles_df))
pc

autoplot(pc)

popu_data <- read.csv("igsr_populations.tsv", header = TRUE, sep = "\t")

fl <- "phase1_integrated_calls.20101123.ALL.panel"

sample_vec <- c()
con <- file(fl, open = "r")
sampleline <- readLines(con, n = ncol(filtered_genotypes_df))

for (line in 1:ncol(filtered_genotypes_df)) {
  sample <- strsplit(sampleline[line], "\t")[[1]][2]
  sample_vec <- c(sample_vec, sample)
}

alt_alleles_df_t <- as.data.frame(t(alt_alleles_df))
View(alt_alleles_df_t)
alt_alleles_df_t["popu"] <- sample_vec
alt_alleles_df_t["popu"]

pc <- prcomp(na.omit(alt_alleles_df_t[,1:(ncol(alt_alleles_df_t)-2)]))
pc

library(ggplot2)
library(dplyr)

autoplot(pc, data = alt_alleles_df_t, color="popu") +
  labs(title="PCA Plot of Populations")


data.frame(pc$x, popu=alt_alleles_df_t$popu) %>%
ggplot(aes(x=PC1,y=PC2,fill=popu)) +
  geom_point(aes(col=popu)) + #Size and alpha just for fun
  scale_color_manual(values=c("dodgerblue2", 
                              "green4", # red
                              "#E31A1C",
                              "#6A3D9A", # purple
                              "darkorange", # orange
                              "cyan3", 
                              "deeppink4",
                              "skyblue2", 
                              "#FB9A99", # lt pink
                              "darkblue",
                              "black", # lt purple
                              "brown", # lt orange
                              "azure4", 
                              "darkolivegreen",
                              "maroon")) +
  labs(title="PCA Plot of Populations") +
  theme_classic()

ggsave("PCAPopu.png", width = 7, height = 7, units = "in")

sample_vec2 <- c()
con <- file(fl, open = "r")
sampleline <- readLines(con, n = ncol(filtered_genotypes_df))

for (line in 1:ncol(filtered_genotypes_df)) {
  sample <- strsplit(sampleline[line], "\t")[[1]][3]
  sample_vec2 <- c(sample_vec2, sample)
}

alt_alleles_df_t["superpopu"] <- sample_vec2

data.frame(pc$x, superpopu=alt_alleles_df_t$superpopu) %>%
  ggplot(aes(x=PC1,y=PC2,fill=superpopu)) +
  geom_point(aes(col=superpopu)) + #Size and alpha just for fun
  scale_color_manual(values=c("dodgerblue2", 
                              "#E31A1C", # red
                              "green4",
                              "#6A3D9A")) +
  labs(title="PCA Plot of Superpopulations") +
  theme_classic()

ggsave("PCASuper.png", width=7, height=7, units = "in")
