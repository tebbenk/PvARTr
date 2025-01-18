library(edgeR)
library(ggplot2)
library(pcaMethods)
library(mixOmics)
library(tidyverse)
library(dplyr)
library(ggrepel)

##Pv
##Read in counts
pv_counts <- read.delim("pv_counts_9724_unstranded.txt", comment.char="#")
pv_counts_countsonly <- pv_counts[,c(7:118)]
rownames(pv_counts_countsonly) <- pv_counts$Geneid

for ( col in 1:ncol(pv_counts_countsonly)){
  colnames(pv_counts_countsonly)[col] <-  sub("_pv_subset_rmdup_sorted.bam*", "", colnames(pv_counts_countsonly)[col])
}

exclude <- c("PQRC21_046_D0_2HR", "PQRC21_021_R3_D0", "PQRC21_021_R3_D0_2hr", "PQRC21_009_D0", "PQRC21_008_D0_8HR_B")

pv_counts_countsonly <- pv_counts_countsonly[,!(colnames(pv_counts_countsonly) %in% exclude)]

gff <- read.delim("/Volumes/projects-t3/SerreDLab-3/kieran.tebben/mali/genomes/PlasmoDB-54_PvivaxP01.gff", comment.char="#", header = FALSE)
rRNA <- subset(gff, V3 == "rRNA")
rRNA <- separate(data = rRNA, col = V9, into = c("V9", "V10", "V11"), sep = ";")
rRNA <- separate(data = rRNA, col = V10, into = c("ID", "Gene"), sep = "=")
pv_counts_countsonly <- pv_counts_countsonly[!(rownames(pv_counts_countsonly) %in% rRNA$Gene),]

#Normalize count data and remove lowly expressed genes
pv_counts_countsonly[] <- sapply(pv_counts_countsonly, as.numeric)
dglist_pv <- DGEList(counts=pv_counts_countsonly, genes=pv_counts_countsonly[c(1)])

keep <- rowSums(cpm(dglist_pv)>10) >= (ncol(dglist_pv$counts)/2)
table(keep)
dglist_pv_filter <- dglist_pv[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_pv_filter))
dim(dglist_pv_filter)

dglist_pv_filter <- calcNormFactors(dglist_pv_filter, method="TMM")

#CIBERSORT Input File
#norm_counts <- as.data.frame(cpm(dglist_pv_filter$counts))
#norm_counts <- na.omit(norm_counts)
#norm_counts$Gene <- rownames(norm_counts)
#write.table(norm_counts, file = "normalized_pv_counts_ARTr_92324.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Read in meta data
pheno <- read.delim("AR_phenotypes.txt") 
rownames(pheno) <- pheno$Sample
pheno$alt_phenotype <- ifelse(pheno$Phenotype == "Average", "Average", "Resistant")

#Pv PCA - all samples
count <- as.data.frame(cpm(dglist_pv_filter))
count <- count[,colnames(count) %in% pheno$Sample]
count <- as.data.frame(t(count))
pca <- prcomp(count)
pcs <- as.data.frame(pca$x)
loadings <- as.data.frame(pca$rotation)

#Plot PCA with ggplot
pcs <- merge(pheno, pcs, by = 0)
ggplot(pcs, aes(x = PC1, y = PC2, color = alt_phenotype, label = ID_Deidentified)) + 
  geom_point(size = 4) +
  theme_linedraw() + 
  labs(color = "Phenotype") + 
  theme(text = element_text(size = 20))+ 
  geom_text(hjust=1.5, vjust=1.5) + 
  scale_color_manual(values = c("black", "grey", "purple", "gold"))


#Baseline PCA
pheno_sub <- subset(pheno, hour == 0)
counts_baseline <- pv_counts_countsonly[,colnames(pv_counts_countsonly) %in% pheno_sub$Sample]
pheno_sub <- pheno_sub[pheno_sub$Sample %in% colnames(counts_baseline),]
pheno_sub <- pheno_sub[order(match(pheno_sub$Sample, colnames(counts_baseline))),]
counts_baseline[] <- sapply(counts_baseline, as.numeric)

dglist_baseline <- DGEList(counts=counts_baseline, genes=rownames(counts_baseline))

keep <- rowSums(cpm(dglist_baseline)>10) >= (ncol(dglist_baseline$counts)/2)
table(keep)
dglist_baseline_filter <- dglist_baseline[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_baseline_filter))
dim(dglist_baseline_filter)

dglist_baseline_filter <- calcNormFactors(dglist_baseline_filter, method="TMM")

count <- as.data.frame(cpm(dglist_baseline_filter))
count <- count[,colnames(count) %in% pheno_sub$Sample]
count <- as.data.frame(t(count))
pca <- prcomp(count)
pcs <- as.data.frame(pca$x)
loadings <- as.data.frame(pca$rotation)

pcs <- merge(pheno_sub, pcs, by = 0)
ggplot(pcs, aes(x = PC1, y = PC2, color = alt_phenotype, label = Jean)) + 
  geom_point(size = 4) +
  theme_linedraw() + 
  labs(color = "Phenotype") + 
  theme(text = element_text(size = 20))+ 
  geom_text(hjust=1.5, vjust=1.5) + 
  scale_color_manual(values = c("black", "grey"), label = c("Fast-clearing", "Slow-clearing"))

#Correlation of slope half life and PCs
summary(lm(PC1~ART_halflife, data = pcs))
ggplot(pcs, aes(x = PC1, y = ART_halflife)) + 
  geom_point(size = 4) + 
  theme_linedraw() + 
  theme(text = element_text(size = 20)) + 
  ylab("Slope Half-life of ART") 

summary(lm(PC2~ART_halflife, data = pcs))
ggplot(pcs, aes(x = PC2, y = ART_halflife)) + 
  geom_point(size = 4) + 
  theme_linedraw() + 
  theme(text = element_text(size = 20)) + 
  ylab("Slope Half-life of ART")

#Baseline differences in gene expression
pheno_sub <- subset(pheno, hour == 0)

counts_baseline <- pv_counts_countsonly[,colnames(pv_counts_countsonly) %in% pheno_sub$Sample]
pheno_sub <- pheno_sub[pheno_sub$Sample %in% colnames(counts_baseline),]
pheno_sub <- pheno_sub[order(match(pheno_sub$Sample, colnames(counts_baseline))),]
counts_baseline[] <- sapply(counts_baseline, as.numeric)

dglist_baseline <- DGEList(counts=counts_baseline, genes=rownames(counts_baseline))

keep <- rowSums(cpm(dglist_baseline)>10) >= (ncol(dglist_baseline$counts)/2)
table(keep)
dglist_baseline_filter <- dglist_baseline[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_baseline_filter))
dim(dglist_baseline_filter)

dglist_baseline_filter <- calcNormFactors(dglist_baseline_filter, method="TMM")

id <- as.factor(pheno_sub$ID)
phenotype <- as.factor(pheno_sub$alt_phenotype)

celltypes <- read.delim("pv_ART_cibersort_92724.txt")
celltypes <- celltypes[celltypes$Mixture %in% colnames(dglist_baseline_filter),]
celltypes <- celltypes[ order(match(celltypes$Mixture, colnames(dglist_baseline_filter))), ]

ring <- as.numeric(celltypes$ring)
troph <- as.numeric(celltypes$troph)
schiz <- as.numeric(celltypes$schiz)
M <- as.numeric(celltypes$M)

design_baseline <- model.matrix(~phenotype + ring + troph + schiz + M)
rownames(design_baseline) <- colnames(dglist_baseline_filter)
design_baseline

y_baseline <- estimateDisp(dglist_baseline_filter, design_baseline)
fit_baseline <- glmFit(y_baseline, design_baseline)
lrt_baseline <- glmLRT(fit_baseline, coef = 2)

topTags(lrt_baseline)
summary(decideTests(lrt_baseline))

results_baseline <- as.data.frame(lrt_baseline$table)
results_baseline$FDR <- p.adjust(results_baseline$PValue, method="fdr")

results_baseline$deexpressed <- ifelse(results_baseline$FDR <= 0.1 & results_baseline$logFC >= 0, "UP", 
                                       ifelse(results_baseline$FDR <= 0.1 & results_baseline$logFC < 0, "DOWN", "NO"))
results_baseline$gene <- rownames(results_baseline)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(results_baseline) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.1)) +
  scale_color_manual(values=mycolors) + 
  theme_classic() + theme(text = element_text(size = 30), legend.position = "none")

#write.table(results_baseline, file = "baseline_DEG_92724.txt", sep = "\t", quote = FALSE)
#sig_baseline <- subset(results_baseline, FDR <= 0.1)
celltypes$Sample <- celltypes$Mixture
df <- merge(pheno_sub, celltypes, by = "Sample")

ggplot(df, aes(x = alt_phenotype, y = schiz, fill = alt_phenotype)) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab("Proportion of Schizonts") + 
  xlab("Phenotype") + 
  theme(legend.position = "none", text = element_text(size = 25)) + 
  ylim(0,1) + 
  stat_compare_means(method = "t.test", size = 8) + 
  scale_fill_manual(values = c("darkgrey", "white")) + 
  scale_x_discrete(labels = c("Fast-clearing", "Slow-clearing"))

##Baseline differences in the proportion of rings
celltypes$Sample <- celltypes$Mixture
celltypes_vars <- merge(celltypes, pheno_sub, by = "Sample")

ggplot(celltypes_vars, aes(x = alt_phenotype, y = ring, fill = alt_phenotype)) + 
  geom_boxplot() +
  theme_classic() + 
  ylab("% Rings") +
  theme(legend.position = "none", text = element_text(size = 25)) + 
  stat_compare_means(method = "t.test", size = 8, label.x = 1, label.y = 1) + 
  scale_fill_manual(values = c("darkgrey", "white")) + 
  scale_x_discrete(labels = c("Fast-clearing", "Slow-clearing")) + 
  xlab("Phenotype")


## Baseline differences in hemoglobin gene expression
counts <- as.data.frame(cpm(dglist_baseline_filter))
counts_t <- as.data.frame(t(counts))
counts_t_merge <- merge(counts_t, pheno_sub, by = 0)

falcilysin <- as.data.frame(counts_t_merge[,(colnames(counts_t_merge) == "PVP01_1111200") | (colnames(counts_t_merge) == "alt_phenotype")])
PIBP <- as.data.frame(counts_t_merge[,(colnames(counts_t_merge) == "PVP01_0316400") | (colnames(counts_t_merge) == "alt_phenotype")])
VP3 <- as.data.frame(counts_t_merge[,(colnames(counts_t_merge) == "PVP01_0916100") | (colnames(counts_t_merge) == "alt_phenotype")])
VPS45 <- as.data.frame(counts_t_merge[,(colnames(counts_t_merge) == "PVP01_0408600") | (colnames(counts_t_merge) == "alt_phenotype")])

ggplot(falcilysin, aes(x = alt_phenotype, y = PVP01_1111200, fill = alt_phenotype)) + 
  geom_boxplot() +
  theme_classic() + 
  ylab("Falcilysin Expression") +
  theme(legend.position = "none", text = element_text(size = 25)) + 
  stat_compare_means(method = "t.test", size = 8) + 
  scale_fill_manual(values = c("darkgrey", "white")) + 
  scale_x_discrete(labels = c("Fast-clearing", "Slow-clearing")) + 
  xlab("Phenotype")

ggplot(PIBP, aes(x = alt_phenotype, y = PVP01_0316400, fill = alt_phenotype)) + 
  geom_boxplot() +
  theme_classic() + 
  ylab("PIBP Expression") +
  theme(legend.position = "none", text = element_text(size = 25)) + 
  stat_compare_means(method = "t.test", size = 8) + 
  scale_fill_manual(values = c("darkgrey", "white")) + 
  scale_x_discrete(labels = c("Fast-clearing", "Slow-clearing")) + 
  xlab("Phenotype")

ggplot(VP3, aes(x = alt_phenotype, y = PVP01_0916100, fill = alt_phenotype)) + 
  geom_boxplot() +
  theme_classic() + 
  ylab("VP3 Expression") +
  theme(legend.position = "none", text = element_text(size = 25)) + 
  stat_compare_means(method = "t.test", size = 8) + 
  scale_fill_manual(values = c("darkgrey", "white")) + 
  scale_x_discrete(labels = c("Fast-clearing", "Slow-clearing")) + 
  xlab("Phenotype")

ggplot(VPS45, aes(x = alt_phenotype, y = PVP01_0408600, fill = alt_phenotype)) + 
  geom_boxplot() +
  theme_classic() + 
  ylab("VPS45 Expression") +
  theme(legend.position = "none", text = element_text(size = 25)) + 
  stat_compare_means(method = "t.test", size = 8) + 
  scale_fill_manual(values = c("darkgrey", "white")) + 
  scale_x_discrete(labels = c("Fast-clearing", "Slow-clearing")) + 
  xlab("Phenotype")

######## WT SAMPLES #########
#Differential expression of paired samples at time 0 vs. time 1 hour
pheno_WT <- subset(pheno, Phenotype == "Average")

paired_subset_0v1 <- pheno_WT %>%
  group_by(ID_Deidentified) %>%
  filter(hour == 0 | hour == 1) 
pairs <- paired_subset_0v1 %>% group_by(ID) %>% count(ID) %>% filter(n == 2)

paired_subset_0v1 <- paired_subset_0v1[paired_subset_0v1$ID %in% pairs$ID,]

pv_counts <- pv_counts_countsonly[,colnames(pv_counts_countsonly) %in% paired_subset_0v1$Sample]

dglist_pv_0v1 <- DGEList(counts=pv_counts, genes=rownames(pv_counts))

keep <- rowSums(cpm(dglist_pv_0v1)>10) >= (ncol(dglist_pv_0v1$counts)/2)
table(keep)
dglist_pv_filter_0v1 <- dglist_pv_0v1[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_pv_filter_0v1))
dim(dglist_pv_filter_0v1)

dglist_pv_filter_0v1 <- calcNormFactors(dglist_pv_filter_0v1, method="TMM")

paired_subset_0v1 <- paired_subset_0v1[order(match(paired_subset_0v1$Sample, colnames(dglist_pv_filter_0v1))),]

celltypes <- read.delim("pv_ART_cibersort_92724.txt")
celltypes <- celltypes[celltypes$Mixture %in% colnames(dglist_pv_filter_0v1),]
celltypes <- celltypes[ order(match(celltypes$Mixture, colnames(dglist_pv_filter_0v1))), ]

ring <- as.numeric(celltypes$ring)
troph <- as.numeric(celltypes$troph)
schiz <- as.numeric(celltypes$schiz)
M <- as.numeric(celltypes$M)

id <- as.factor(paired_subset_0v1$ID)
hour <- as.factor(paired_subset_0v1$hour)

design_0v1 <- model.matrix(~id + hour)
rownames(design_0v1) <- colnames(dglist_pv_filter_0v1)
design_0v1

y0v1 <- estimateDisp(dglist_pv_filter_0v1, design_0v1)
fit0v1 <- glmFit(y0v1, design_0v1)
lrt0v1 <- glmLRT(fit0v1, coef = 15)

topTags(lrt0v1)
summary(decideTests(lrt0v1))

results0v1 <- as.data.frame(lrt0v1$table)
results0v1$FDR <- p.adjust(results0v1$PValue, method="fdr")

results0v1$deexpressed <- ifelse(results0v1$FDR <= 0.1 & results0v1$logFC >= 0, "UP", 
                                 ifelse(results0v1$FDR <= 0.1 & results0v1$logFC < 0, "DOWN", "NO"))
results0v1$gene <- rownames(results0v1)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

results0v1 <- subset(results0v1, gene != "PVP01_0202700" & gene != "PVP01_0202800")
#write.table(results0v1, file = "DEG_WT_0v1_92724.txt", sep = "\t", quote = FALSE)

ggplot(results0v1) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.1)) +
  scale_color_manual(values=mycolors) + 
  theme_classic() + theme(text = element_text(size = 30), legend.position = "none") + 
  xlim(-2, 2)

sig0v1_WT <- subset(results0v1, FDR <= 0.1)

celltypes$Sample <- celltypes$Mixture

df2 <- merge(celltypes, paired_subset_0v1, by = "Sample")
ggplot(df2, aes(x = as.factor(hour), y = schiz, fill = as.factor(hour))) + 
  geom_boxplot() + 
  theme_classic() + 
  xlab("Hour") + ylab("Proportion of Schizonts") +
  theme(text = element_text(size = 30), legend.position = "none") +
  stat_compare_means(method = "t.test", paired = TRUE, label.y = 1, size = 8) +
  scale_fill_manual(values = c("darkgrey", "white"))

##Differential expression of paired samples at time 1 hour vs time 2 hours
paired_subset_1v2 <- pheno_WT %>%
  group_by(ID_Deidentified) %>%
  filter(hour == 1 | hour == 2) 
pairs <- paired_subset_1v2 %>% group_by(ID) %>% count(ID) %>% filter(n == 2)

paired_subset_1v2 <- paired_subset_1v2[paired_subset_1v2$ID %in% pairs$ID,]

pv_counts <- pv_counts_countsonly[,colnames(pv_counts_countsonly) %in% paired_subset_1v2$Sample]
paired_subset_1v2 <- paired_subset_1v2[paired_subset_1v2$Sample %in% colnames(pv_counts),]

dglist_pv_1v2 <- DGEList(counts=pv_counts, genes=rownames(pv_counts))

keep <- rowSums(cpm(dglist_pv_1v2)>10) >= (ncol(dglist_pv_1v2$counts)/2)
table(keep)
dglist_pv_filter_1v2 <- dglist_pv_1v2[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_pv_filter_1v2))
dim(dglist_pv_filter_1v2)

dglist_pv_filter_1v2 <- calcNormFactors(dglist_pv_filter_1v2, method="TMM")

paired_subset_1v2 <- paired_subset_1v2[order(match(paired_subset_1v2$Sample, colnames(dglist_pv_filter_1v2))),]

id <- as.factor(paired_subset_1v2$ID)
hour <- as.factor(paired_subset_1v2$hour)

design_1v2 <- model.matrix(~id + hour)
rownames(design_1v2) <- colnames(dglist_pv_filter_1v2)
design_1v2

y1v2 <- estimateDisp(dglist_pv_filter_1v2, design_1v2)
fit1v2 <- glmFit(y1v2, design_1v2)
lrt1v2 <- glmLRT(fit1v2, coef = 15)

topTags(lrt1v2)
summary(decideTests(lrt1v2))

results1v2 <- as.data.frame(lrt1v2$table)
results1v2$FDR <- p.adjust(results1v2$PValue, method="fdr")

results1v2$deexpressed <- ifelse(results1v2$FDR <= 0.1 & results1v2$logFC >= 0, "UP", 
                                 ifelse(results1v2$FDR <= 0.1 & results1v2$logFC < 0, "DOWN", "NO"))
results1v2$gene <- rownames(results1v2)
results1v2 <- subset(results1v2, gene != "PVP01_0202700" & gene != "PVP01_0202800")
#write.table(results1v2, file = "DEG_WT_1v2_92724.txt", sep = "\t", quote = FALSE, row.names = FALSE)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(results1v2) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.1)) +
  scale_color_manual(values=mycolors) + 
  theme_classic() + theme(text = element_text(size = 30), legend.position = "none") + xlim(-2,2)

sig1v2_WT <- subset(results1v2, FDR <= 0.1)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

celltypes <- read.delim("pv_ART_cibersort_92724.txt")
celltypes <- celltypes[celltypes$Mixture %in% colnames(dglist_pv_filter_1v2),]
celltypes <- celltypes[ order(match(celltypes$Mixture, colnames(dglist_pv_filter_1v2))), ]
celltypes$Sample <- celltypes$Mixture

df2 <- merge(celltypes, paired_subset_1v2, by = "Sample")
ggplot(df2, aes(x = as.factor(hour), y = ring, fill = as.factor(hour))) + 
  geom_boxplot() + 
  theme_classic() + 
  xlab("Hour") + ylab("Proportion of Schizonts") +
  theme(text = element_text(size = 30), legend.position = "none") +
  stat_compare_means(method = "t.test", paired = TRUE, label.y = 1, size = 8) + 
  scale_fill_manual(values = c("darkgrey", "white"))


#Names on plot
mygenes_down <- c("RPN5", "RPN12", "RPN13", "RPN3", "RPN8", "RPN9", "RPS10", "RPS12", "RPS19", "RPS21", "RPL34", "RPL44", "GluRS", "SEC12", "SEC13")
mygenes_up <- c("eIK1", "PK4")

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

genenames <- read.delim("pv_codinggenes.txt")
results1v2_genenames <- merge(results1v2, genenames, by = "gene", all = TRUE)

results1v2_genenames <- subset(results1v2_genenames, deexpressed != "NA")

ggplot(results1v2_genenames) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.1)) +
  geom_label_repel(data = results1v2_genenames %>% 
                     filter(results1v2_genenames$name %in% mygenes_down), 
                   aes(label = name, 
                       x = logFC, 
                       y = -log10(PValue)), 
                   size = 5, max.overlaps = 12,
                   color = "darkblue") +
  geom_label_repel(data = results1v2_genenames %>% 
                     filter(results1v2_genenames$name %in% mygenes_up), 
                   aes(label = name, 
                       x = logFC, 
                       y = -log10(PValue)), 
                   size = 5, max.overlaps = 12,
                   color = "darkred") +
  scale_color_manual(values=mycolors) + 
  theme_classic() + theme(text = element_text(size = 30), legend.position = "none") + xlab("LogFC") + ylab("-log10(PValue)") + xlim(-2,2)



##Differential expression of paired samples at time 4 hours vs time 2 hours
paired_subset_2v4 <- pheno_WT %>%
  group_by(ID_Deidentified) %>%
  filter(hour == 2 | hour == 4) 
pairs <- paired_subset_2v4 %>% group_by(ID) %>% count(ID) %>% filter(n == 2)

paired_subset_2v4 <- paired_subset_2v4[paired_subset_2v4$ID %in% pairs$ID,]
paired_subset_2v4 <- subset(paired_subset_2v4, ID_Deidentified != "M")

pv_counts <- pv_counts_countsonly[,colnames(pv_counts_countsonly) %in% paired_subset_2v4$Sample]
paired_subset_2v4 <- paired_subset_2v4[paired_subset_2v4$Sample %in% colnames(pv_counts),]

dglist_pv_2v4 <- DGEList(counts=pv_counts, genes=rownames(pv_counts))

keep <- rowSums(cpm(dglist_pv_2v4)>10) >= (ncol(dglist_pv_2v4$counts)/2)
table(keep)
dglist_pv_filter_2v4 <- dglist_pv_2v4[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_pv_filter_2v4))
dim(dglist_pv_filter_2v4)

dglist_pv_filter_2v4 <- calcNormFactors(dglist_pv_filter_2v4, method="TMM")

paired_subset_2v4 <- paired_subset_2v4[order(match(paired_subset_2v4$Sample, colnames(dglist_pv_filter_2v4))),]

celltypes <- read.delim("pv_ART_cibersort_92724.txt")
celltypes <- celltypes[celltypes$Mixture %in% colnames(dglist_pv_filter_2v4),]
celltypes <- celltypes[ order(match(celltypes$Mixture, colnames(dglist_pv_filter_2v4))), ]
celltypes$Sample <- celltypes$Mixture

ring <- as.numeric(celltypes$ring)
troph <- as.numeric(celltypes$troph)
schiz <- as.numeric(celltypes$schiz)
M <- as.numeric(celltypes$M)

id <- as.factor(paired_subset_2v4$ID)
hour <- as.factor(paired_subset_2v4$hour)

design_2v4 <- model.matrix(~id + hour)
rownames(design_2v4) <- colnames(dglist_pv_filter_2v4)
design_2v4

y2v4 <- estimateDisp(dglist_pv_filter_2v4, design_2v4)
fit2v4 <- glmFit(y2v4, design_2v4)
lrt2v4 <- glmLRT(fit2v4, coef = 12)

topTags(lrt2v4)
summary(decideTests(lrt2v4))

results2v4 <- as.data.frame(lrt2v4$table)
results2v4$FDR <- p.adjust(results2v4$PValue, method="fdr")

results2v4$deexpressed <- ifelse(results2v4$FDR <= 0.1 & results2v4$logFC >= 0, "UP", 
                                 ifelse(results2v4$FDR <= 0.1 & results2v4$logFC < 0, "DOWN", "NO"))
results2v4$gene <- rownames(results2v4)
#write.table(results2v4, file = "DEG_WT_2v4_92724.txt", sep = "\t", quote = FALSE, row.names = FALSE)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(results2v4) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.1)) +
  scale_color_manual(values=mycolors) + 
  theme_classic() + theme(text = element_text(size = 30), legend.position = "none")

sig2v4_WT <- subset(results2v4, FDR <= 0.1)

df2 <- merge(celltypes, paired_subset_2v4, by = "Sample")
ggplot(df2, aes(x = as.factor(hour), y = schiz, fill = as.factor(hour))) + 
  geom_boxplot() + 
  theme_classic() + 
  xlab("Hour") + ylab("Proportion of Schizonts") +
  theme(text = element_text(size = 30), legend.position = "none") +
  stat_compare_means(method = "t.test", paired = TRUE, label.y = 1, size = 8) + 
  scale_fill_manual(values = c("darkgrey", "white"))


######## RESISTANT AND SLOW PARASITES DE ######
#Differential expression of paired samples at time 0 vs. time 1 hour
pheno_resistant <- subset(pheno, Phenotype == "Resistant" | Phenotype == "Slow")

paired_subset_0v1 <- pheno_resistant %>%
  group_by(ID_Deidentified) %>%
  filter(hour == 0 | hour == 1) 
pairs <- paired_subset_0v1 %>% group_by(ID) %>% count(ID) %>% filter(n == 2)

paired_subset_0v1 <- paired_subset_0v1[paired_subset_0v1$ID %in% pairs$ID,]
paired_subset_0v1 <- subset(paired_subset_0v1, ID_Deidentified != "V")

pv_counts <- pv_counts_countsonly[,colnames(pv_counts_countsonly) %in% paired_subset_0v1$Sample]
paired_subset_0v1 <- paired_subset_0v1[paired_subset_0v1$Sample %in% colnames(pv_counts),]

dglist_pv_0v1 <- DGEList(counts=pv_counts, genes=rownames(pv_counts))

keep <- rowSums(cpm(dglist_pv_0v1)>10) >= (ncol(dglist_pv_0v1$counts)/2)
table(keep)
dglist_pv_filter_0v1 <- dglist_pv_0v1[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_pv_filter_0v1))
dim(dglist_pv_filter_0v1)

dglist_pv_filter_0v1 <- calcNormFactors(dglist_pv_filter_0v1, method="TMM")

paired_subset_0v1 <- paired_subset_0v1[order(match(paired_subset_0v1$Sample, colnames(dglist_pv_filter_0v1))),]

celltypes <- read.delim("pv_ART_cibersort_92724.txt")
celltypes <- celltypes[celltypes$Mixture %in% colnames(dglist_pv_filter_0v1),]
celltypes <- celltypes[ order(match(celltypes$Mixture, colnames(dglist_pv_filter_0v1))), ]

ring <- as.numeric(celltypes$ring)
troph <- as.numeric(celltypes$troph)
schiz <- as.numeric(celltypes$schiz)
M <- as.numeric(celltypes$M)
F <- as.numeric(celltypes$F)
id <- as.factor(paired_subset_0v1$ID)
hour <- as.factor(paired_subset_0v1$hour)

design_0v1 <- model.matrix(~id + hour)
rownames(design_0v1) <- colnames(dglist_pv_filter_0v1)
design_0v1

y0v1 <- estimateDisp(dglist_pv_filter_0v1, design_0v1)
fit0v1 <- glmFit(y0v1, design_0v1)
lrt0v1 <- glmLRT(fit0v1, coef = 13)

topTags(lrt0v1)
summary(decideTests(lrt0v1))

results0v1 <- as.data.frame(lrt0v1$table)
results0v1$FDR <- p.adjust(results0v1$PValue, method="fdr")

results0v1$deexpressed <- ifelse(results0v1$FDR <= 0.1 & results0v1$logFC >= 0, "UP", 
                                 ifelse(results0v1$FDR <= 0.1 & results0v1$logFC < 0, "DOWN", "NO"))
results0v1$gene <- rownames(results0v1)
#write.table(results0v1, file = "RES_DEG_0v1_92824.txt", sep = "\t", quote = FALSE, row.names = FALSE)
sig0v1_res <- subset(results0v1, FDR <0.1)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(results0v1) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.1)) +
  scale_color_manual(values=mycolors) + 
  theme_classic() + theme(text = element_text(size = 30), legend.position = "none")

celltypes$Sample <- celltypes$Mixture
df1 <- merge(celltypes, paired_subset_0v1, by = "Sample")
ggplot(df1, aes(x = as.factor(hour), y = troph, fill = as.factor(hour))) + 
  geom_boxplot() + 
  theme_classic() + 
  xlab("Hour") + ylab("Proportion of Trophozoites") +
  theme(text = element_text(size = 30), legend.position = "none") +
  stat_compare_means(method = "t.test", paired = TRUE, label.y = 1, size = 8) + 
  scale_fill_manual(values = c("darkgrey", "white"))

##Differential expression of paired samples at time 1 hour vs time 2 hours
paired_subset_1v2 <- pheno_resistant %>%
  group_by(ID_Deidentified) %>%
  filter(hour == 1 | hour == 2) 
pairs <- paired_subset_1v2 %>% group_by(ID) %>% count(ID) %>% filter(n == 2)

paired_subset_1v2 <- paired_subset_1v2[paired_subset_1v2$ID %in% pairs$ID,]
paired_subset_1v2 <- subset(paired_subset_1v2, ID_Deidentified != "V")

pv_counts <- pv_counts_countsonly[,colnames(pv_counts_countsonly) %in% paired_subset_1v2$Sample]
paired_subset_1v2 <- paired_subset_1v2[paired_subset_1v2$Sample %in% colnames(pv_counts),]

dglist_pv_1v2 <- DGEList(counts=pv_counts, genes=rownames(pv_counts))

keep <- rowSums(cpm(dglist_pv_1v2)>10) >= (ncol(dglist_pv_1v2$counts)/2)
table(keep)
dglist_pv_filter_1v2 <- dglist_pv_1v2[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_pv_filter_1v2))
dim(dglist_pv_filter_1v2)

dglist_pv_filter_1v2 <- calcNormFactors(dglist_pv_filter_1v2, method="TMM")

paired_subset_1v2 <- paired_subset_1v2[order(match(paired_subset_1v2$Sample, colnames(dglist_pv_filter_1v2))),]

celltypes <- read.delim("pv_ART_cibersort_92724.txt")
celltypes <- celltypes[celltypes$Mixture %in% colnames(dglist_pv_filter_1v2),]
celltypes <- celltypes[ order(match(celltypes$Mixture, colnames(dglist_pv_filter_1v2))), ]

id <- as.factor(paired_subset_1v2$ID)
hour <- as.factor(paired_subset_1v2$hour)

design_1v2 <- model.matrix(~id + hour)
rownames(design_1v2) <- colnames(dglist_pv_filter_1v2)
design_1v2

y1v2 <- estimateDisp(dglist_pv_filter_1v2, design_1v2)
fit1v2 <- glmFit(y1v2, design_1v2)
lrt1v2 <- glmLRT(fit1v2, coef = 11)

topTags(lrt1v2)
summary(decideTests(lrt1v2))

results1v2 <- as.data.frame(lrt1v2$table)
results1v2$FDR <- p.adjust(results1v2$PValue, method="fdr")

results1v2$deexpressed <- ifelse(results1v2$FDR <= 0.1 & results1v2$logFC >= 0, "UP", 
                                 ifelse(results1v2$FDR <= 0.1 & results1v2$logFC < 0, "DOWN", "NO"))
results1v2$gene <- rownames(results1v2)
#write.table(results1v2, file = "DEG_res_1v2_92824.txt", sep = "\t", quote = FALSE, row.names = FALSE)
sig1v2_res <- subset(results1v2, FDR < 0.1)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

#Plot with labels
names(mycolors) <- c("DOWN", "UP", "NO")
results1v2$kelchlabel <- ifelse(results1v2$gene == "PVP01_0607800", "K10", "NA")


results1v2$hemoglobinlabel <-     ifelse(results1v2$gene == "PVP01_1111200", "falcilysin",
                                         ifelse(results1v2$gene == "PVP01_0316400", "PIBP",
                                                ifelse(results1v2$gene == "PVP01_0916100", "VP3",
                                                       ifelse(results1v2$gene == "PVP01_0408600", "VPS45", "NA"))))


ggplot(results1v2) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.5)) + 
  scale_colour_manual(values = mycolors) + 
  geom_label_repel(data = subset(results1v2, kelchlabel != "NA"), 
                   aes(logFC,-log10(PValue),label=kelchlabel),
                   size = 8,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   color = "darkred") + 
  geom_label_repel(data = subset(results1v2, hemoglobinlabel != "NA"), 
                   aes(logFC,-log10(PValue),label=hemoglobinlabel),
                   size = 8,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   color = "darkblue") + 
  theme_classic() + 
  theme(text = element_text(size = 30), legend.position = "none")


celltypes$Sample <- celltypes$Mixture
df1 <- merge(celltypes, paired_subset_1v2, by = "Sample")
ggplot(df1, aes(x = as.factor(hour), y = schiz, fill = as.factor(hour))) + 
  geom_boxplot() + 
  theme_classic() + 
  xlab("Hour") + ylab("Proportion of Schizonts") +
  theme(text = element_text(size = 30), legend.position = "none") +
  stat_compare_means(method = "t.test", paired = TRUE, label.y = 1, size = 8) + 
  scale_fill_manual(values = c("darkgrey", "white"))

hour1 <- subset(df1, hour == 1)
hour2 <- subset(df1, hour ==2)

##Differential expression of paired samples at time 4 hours vs time 2 hours
paired_subset_2v4 <- pheno_resistant %>%
  group_by(ID_Deidentified) %>%
  filter(hour == 2 | hour == 4) 
pairs <- paired_subset_2v4 %>% group_by(ID) %>% count(ID) %>% filter(n == 2)

paired_subset_2v4 <- paired_subset_2v4[paired_subset_2v4$ID %in% pairs$ID,]
#paired_subset_2v4 <- subset(paired_subset_2v4, ID_Deidentified != "M")

pv_counts <- pv_counts_countsonly[,colnames(pv_counts_countsonly) %in% paired_subset_2v4$Sample]
paired_subset_2v4 <- paired_subset_2v4[paired_subset_2v4$Sample %in% colnames(pv_counts),]

dglist_pv_2v4 <- DGEList(counts=pv_counts, genes=rownames(pv_counts))

keep <- rowSums(cpm(dglist_pv_2v4)>10) >= (ncol(dglist_pv_2v4$counts)/2)
table(keep)
dglist_pv_filter_2v4 <- dglist_pv_2v4[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_pv_filter_2v4))
dim(dglist_pv_filter_2v4)

dglist_pv_filter_2v4 <- calcNormFactors(dglist_pv_filter_2v4, method="TMM")

paired_subset_2v4 <- paired_subset_2v4[order(match(paired_subset_2v4$Sample, colnames(dglist_pv_filter_2v4))),]

celltypes <- read.delim("pv_ART_cibersort_92724.txt")
celltypes <- celltypes[celltypes$Mixture %in% colnames(dglist_pv_filter_2v4),]
celltypes <- celltypes[ order(match(celltypes$Mixture, colnames(dglist_pv_filter_2v4))), ]

id <- as.factor(paired_subset_2v4$ID)
hour <- as.factor(paired_subset_2v4$hour)

design_2v4 <- model.matrix(~id + hour)
rownames(design_2v4) <- colnames(dglist_pv_filter_2v4)
design_2v4

y2v4 <- estimateDisp(dglist_pv_filter_2v4, design_2v4)
fit2v4 <- glmFit(y2v4, design_2v4)
lrt2v4 <- glmLRT(fit2v4, coef = 11)

topTags(lrt2v4)
summary(decideTests(lrt2v4))

results2v4 <- as.data.frame(lrt2v4$table)
results2v4$FDR <- p.adjust(results2v4$PValue, method="fdr")

results2v4$deexpressed <- ifelse(results2v4$FDR <= 0.1 & results2v4$logFC >= 0, "UP", 
                                 ifelse(results2v4$FDR <= 0.1 & results2v4$logFC < 0, "DOWN", "NO"))
results2v4$gene <- rownames(results2v4)
#write.table(results2v4, file = "DEG_res_2v4_92824.txt", sep = "\t", quote = FALSE, row.names = FALSE)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

results2v4$kelchlabel <- ifelse(results2v4$gene == "PVP01_1304500", "kelch domain-containing protein", 
                                ifelse(results2v4$gene == "PVP01_0926500", "kelch domain-containing protein",
                                       ifelse(results2v4$gene == "PVP01_0926400", "kelch domain-containing protein",
                                              ifelse(results2v4$gene == "PVP01_0607800", "K10", "NA"))))


results2v4$hemoglobinlabel <-     ifelse(results2v4$gene == "PVP01_1111200", "falcilysin",
                                         ifelse(results2v4$gene == "PVP01_0316400", "PIBP",
                                                ifelse(results2v4$gene == "PVP01_0916100", "VP3", 
                                                       ifelse(results2v4$gene == "PVP01_1018600", "PIP3", "NA"))))


ggplot(results2v4) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.5)) + 
  scale_colour_manual(values = mycolors) + 
  geom_label_repel(data = subset(results2v4, kelchlabel != "NA"), 
                   aes(logFC,-log10(PValue),label=kelchlabel),
                   size = 5,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   color = "darkred",
                   max.overlaps = 4) + 
  geom_label_repel(data = subset(results2v4, hemoglobinlabel != "NA"), 
                   aes(logFC,-log10(PValue),label=hemoglobinlabel),
                   size = 4,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   color = "darkblue") + 
  theme_classic() + 
  theme(text = element_text(size = 30), legend.position = "none")

sig2v4_res <- subset(results2v4, FDR <= 0.1)

celltypes$Sample <- celltypes$Mixture
df2 <- merge(celltypes, paired_subset_2v4, by = "Sample")
ggplot(df2, aes(x = as.factor(hour), y = schiz, fill = as.factor(hour))) + 
  geom_boxplot() + 
  theme_classic() + 
  xlab("Hour") + ylab("Proportion of Schizonts") +
  theme(text = element_text(size = 30), legend.position = "none") +
  stat_compare_means(method = "t.test", paired = TRUE, label.y = 1, size = 8) + 
  scale_fill_manual(values = c("white", "lightgrey"))

hour2 <- subset(df2, hour ==2 )
hour4 <- subset(df2, hour == 4)

#Calculate overlap of genes between WT and resistan
overlap0v1 <- merge(sig0v1_WT, sig0v1_res, by = "gene")
overlap1v2 <- merge(sig1v2_WT, sig1v2_res, by = "gene")
overlap2v4 <- merge(sig2v4_WT, sig2v4_res, by = "gene")

#Schematic figures
#WT
wt <- read.delim("WT.txt")

ggplot(wt, aes(x = Position, y = DEG)) + 
  geom_bar(position="dodge", stat="identity", fill = "lightgrey") + theme_classic() + 
  scale_x_continuous(limits = c(0,4.5), breaks = c(0,1,2,4),labels = c("Baseline", "1 Hour", "2 Hours", "4 Hours")) +
  theme(text = element_text(size = 20)) + 
  ylab("Number of Differentially Expressed Genes") + 
  xlab(NULL) + 
  annotate("text", x = 0.5, y = 95, label = "178", size = 7) +
  annotate("text", x = 1.5, y = 500, label = "783", size = 7) + 
  annotate("text", x = 3, y = 15, label = "0", size = 7) + 
  annotate("text", x = 0.5, y = 195, label = "N = 14 Pairs", size = 6) +
  annotate("text", x = 1.5, y = 800, label = "N = 14 Pairs", size = 6) + 
  annotate("text", x = 3, y = 50, label = "N = 11 Pairs", size = 6)

res <- read.delim("res_DEG_table.txt")
res$Group <- factor(res$Group, levels = c("unique", "shared"))

ggplot(res, aes(x = Position, y = DEG, fill = Group)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
  scale_x_continuous(limits = c(0,4.5), breaks = c(0,1,2,4),labels = c("Baseline", "1 Hour", "2 Hours", "4 Hours")) +
  theme(text = element_text(size = 20)) + 
  ylab("Number of Differentially Expressed Genes") + 
  xlab(NULL) + 
  annotate("text", x = 0.5, y = 0, label = "1", size = 7) +
  annotate("text", x = 1.5, y = 1100, label = "1,227", size = 7) + 
  annotate("text", x = 1.5, y = 250, label = "476", size = 7) + 
  annotate("text", x = 3, y = 1500, label = "2,445", size = 7) +
  annotate("text", x = 0.5, y = 150, label = "N = 12 Pairs", size = 6) + 
  annotate("text", x = 1.5, y = 1775, label = "N = 10 Pairs", size = 6) + 
  annotate("text", x = 3, y = 2500, label = "N = 10 Pairs", size = 6) + 
  scale_fill_manual(name = "Origin of DEG",labels = c("Unique to resistant", "Shared with susceptible"), values = c("red", "lightgrey")) + 
  theme(legend.position = "bottom")


########## HUMAN ############
hu_counts <- read.delim("human_counts_9724_unstranded.txt", comment.char="#")
hu_counts <- na.omit(hu_counts)
hu_counts_countsonly <- hu_counts[,c(7:118)]
rownames(hu_counts_countsonly) <- hu_counts$Geneid

for ( col in 1:ncol(hu_counts_countsonly)){
  colnames(hu_counts_countsonly)[col] <-  sub("_subset_human_rmdup.bam*", "", colnames(hu_counts_countsonly)[col])
}

exclude <- c("PQRC21_046_D0_2HR", "PQRC21_021_R3_D0", "PQRC21_021_R3_D0_2hr", "PQRC21_009_D0", "PQRC21_008_D0_8HR_B")

hu_counts_countsonly <- hu_counts_countsonly[,!(colnames(hu_counts_countsonly) %in% exclude)]

#Normalize count data and remove lowly expressed genes
hu_counts_countsonly[] <- sapply(hu_counts_countsonly, as.numeric)
hu_counts_countsonly <- na.omit(hu_counts_countsonly)
dglist_hu <- DGEList(counts=hu_counts_countsonly, genes=rownames(hu_counts_countsonly))

keep <- rowSums(cpm(dglist_hu)>10) >= (ncol(dglist_hu$counts)/2)
table(keep)
dglist_hu_filter <- dglist_hu[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_hu_filter))
dim(dglist_hu_filter)

dglist_hu_filter <- calcNormFactors(dglist_hu_filter, method="TMM")

#CIBERSORT Input File
# norm_counts <- as.data.frame(cpm(dglist_hu_filter$counts))
# norm_counts <- na.omit(norm_counts)
# norm_counts$Gene <- rownames(norm_counts)
# write.table(norm_counts, file = "normalized_human_counts_ARTr_92324.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Read in meta data
pheno <- read.delim("AR_phenotypes.txt")
rownames(pheno) <- pheno$Sample
pheno$alt_phenotype <- ifelse(pheno$Phenotype == "Average", "Average", "Resistant")

#Human PCA - all samples
count <- as.data.frame(cpm(dglist_hu_filter))
count <- as.data.frame(t(count))
pca <- prcomp(count)
pcs <- as.data.frame(pca$x)

#Plot PCA with ggplot
pcs <- merge(pheno, pcs, by = 0)
ggplot(pcs, aes(x = PC1, y = PC2, color = Phenotype)) + 
  geom_point() +
  theme_minimal() + 
  ggtitle("Human")

#Human PCA - baseline
pheno_sub <- subset(pheno, hour == 0)
counts_baseline <- hu_counts_countsonly[,colnames(hu_counts_countsonly) %in% pheno_sub$Sample]
pheno_sub <- pheno_sub[pheno_sub$Sample %in% colnames(counts_baseline),]
pheno_sub <- pheno_sub[order(match(pheno_sub$Sample, colnames(counts_baseline))),]
counts_baseline[] <- sapply(counts_baseline, as.numeric)

dglist_baseline <- DGEList(counts=counts_baseline, genes=rownames(counts_baseline))

keep <- rowSums(cpm(dglist_baseline)>10) >= (ncol(dglist_baseline$counts)/2)
table(keep)
dglist_baseline_filter <- dglist_baseline[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_baseline_filter))
dim(dglist_baseline_filter)

dglist_baseline_filter <- calcNormFactors(dglist_baseline_filter, method="TMM")

count <- as.data.frame(cpm(dglist_baseline_filter))
count <- count[,colnames(count) %in% pheno_sub$Sample]
count <- as.data.frame(t(count))
pca <- prcomp(count)
pcs <- as.data.frame(pca$x)
loadings <- as.data.frame(pca$rotation)

pcs <- merge(pheno_sub, pcs, by = 0)
ggplot(pcs, aes(x = PC2, y = PC3, color = alt_phenotype, label = ID_Deidentified)) + 
  geom_point(size = 4) +
  theme_linedraw() + 
  labs(color = "Phenotype") + 
  theme(text = element_text(size = 20))+ 
  geom_text(hjust=1.5, vjust=1.5) + 
  scale_color_manual(values = c("black", "grey", "purple", "gold"))

#Baseline differences in gene expression - human
pheno_sub <- subset(pheno, hour == 0)

counts_baseline <- hu_counts_countsonly[,colnames(hu_counts_countsonly) %in% pheno_sub$Sample]
pheno_sub <- pheno_sub[pheno_sub$Sample %in% colnames(counts_baseline),]
pheno_sub <- pheno_sub[order(match(pheno_sub$Sample, colnames(counts_baseline))),]
counts_baseline[] <- sapply(counts_baseline, as.numeric)

dglist_baseline <- DGEList(counts=counts_baseline, genes=rownames(counts_baseline))

keep <- rowSums(cpm(dglist_baseline)>10) >= (ncol(dglist_baseline$counts)/2)
table(keep)
dglist_baseline_filter <- dglist_baseline[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_baseline_filter))
dim(dglist_baseline_filter)

dglist_baseline_filter <- calcNormFactors(dglist_baseline_filter, method="TMM")

id <- as.factor(pheno_sub$ID)
phenotype <- as.factor(pheno_sub$alt_phenotype)
parasitemia <- as.numeric(pheno_sub$Parasitemia)

celltypes <- read.delim("human_cibersort_101824.txt")
celltypes <- celltypes[celltypes$Mixture %in% colnames(dglist_baseline_filter),]
celltypes <- celltypes[ order(match(celltypes$Mixture, colnames(dglist_baseline_filter))), ]

celltypes$B <- celltypes$B.cells.memory + celltypes$B.cells.naive
celltypes$T <- celltypes$T.cells.CD4.memory.activated + celltypes$T.cells.CD4.memory.resting + celltypes$T.cells.CD4.naive + celltypes$T.cells.CD8 + celltypes$T.cells.follicular.helper + celltypes$T.cells.gamma.delta + celltypes$T.cells.regulatory..Tregs.
celltypes$APC <- celltypes$Macrophages.M0 + celltypes$Macrophages.M1 + celltypes$Macrophages.M2 + celltypes$Dendritic.cells.activated + celltypes$Dendritic.cells.resting + celltypes$Monocytes
celltypes$mast <- celltypes$Mast.cells.activated + celltypes$Mast.cells.resting
celltypes$NK <- celltypes$NK.cells.activated + celltypes$NK.cells.resting
celltypes <- celltypes[,c(1,4,22,23,27:31)]

B <- as.numeric(celltypes$B)
Plasma <- as.numeric(celltypes$Plasma.cells)
T_cell <- as.numeric(celltypes$T)
APC <- as.numeric(celltypes$APC)
mast <- as.numeric(celltypes$mast)
NK <- as.numeric(celltypes$NK)
Neut <- as.numeric(celltypes$Neutrophils)

design_baseline <- model.matrix(~phenotype + B + Plasma + T_cell + APC + mast + NK + Neut + parasitemia)
rownames(design_baseline) <- colnames(dglist_baseline_filter)
design_baseline

y_baseline <- estimateDisp(dglist_baseline_filter, design_baseline)
fit_baseline <- glmFit(y_baseline, design_baseline)
lrt_baseline <- glmLRT(fit_baseline, coef = 2)

topTags(lrt_baseline)
summary(decideTests(lrt_baseline))

results_baseline <- as.data.frame(lrt_baseline$table)
results_baseline$FDR <- p.adjust(results_baseline$PValue, method="fdr")

results_baseline$deexpressed <- ifelse(results_baseline$FDR <= 0.1 & results_baseline$logFC >= 0, "UP", 
                                       ifelse(results_baseline$FDR <= 0.1 & results_baseline$logFC < 0, "DOWN", "NO"))
results_baseline$gene <- rownames(results_baseline)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(results_baseline) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.1)) +
  scale_color_manual(values=mycolors) + 
  theme_classic() + theme(text = element_text(size = 30), legend.position = "none")

write.table(results_baseline, file = "baseline_DEG_human_121124.txt", sep = "\t", quote = FALSE)
sig_baseline <- subset(results_baseline, FDR <= 0.1)

celltypes$Sample <- celltypes$Mixture
df <- merge(pheno_sub, celltypes, by = "Sample")

ggplot(df, aes(x = alt_phenotype, y = NK, fill = alt_phenotype)) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab("Proportion of NK Cells") + 
  xlab("Phenotype") + 
  theme(legend.position = "none", text = element_text(size = 25)) + 
  ylim(0,1) + 
  stat_compare_means(method = "t.test", size = 8) + 
  scale_fill_manual(values = c("darkgrey", "white")) + 
  scale_x_discrete(labels = c("Fast-clearing", "Slow-clearing"))

#Human DEG - all samples, paired time points
paired_subset_0v1 <- pheno %>%
  group_by(ID_Deidentified) %>%
  filter(hour == 0 | hour == 1) 
pairs <- paired_subset_0v1 %>% group_by(ID) %>% count(ID) %>% filter(n == 2)

paired_subset_0v1 <- paired_subset_0v1[paired_subset_0v1$ID %in% pairs$ID,]
paired_subset_0v1 <- subset(paired_subset_0v1, ID_Deidentified != "V")

hu_counts <- hu_counts_countsonly[,colnames(hu_counts_countsonly) %in% paired_subset_0v1$Sample]

dglist_hu_0v1 <- DGEList(counts=hu_counts, genes=rownames(hu_counts))

keep <- rowSums(cpm(dglist_hu_0v1)>10) >= (ncol(dglist_hu_0v1$counts)/2)
table(keep)
dglist_hu_filter_0v1 <- dglist_hu_0v1[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_hu_filter_0v1))
dim(dglist_hu_filter_0v1)

dglist_hu_filter_0v1 <- calcNormFactors(dglist_hu_filter_0v1, method="TMM")

paired_subset_0v1 <- paired_subset_0v1[order(match(paired_subset_0v1$Sample, colnames(dglist_hu_filter_0v1))),]

id <- as.factor(paired_subset_0v1$ID)
hour <- as.factor(paired_subset_0v1$hour)

design_0v1 <- model.matrix(~hour + id)
rownames(design_0v1) <- colnames(dglist_hu_filter_0v1)
design_0v1

y0v1 <- estimateDisp(dglist_hu_filter_0v1, design_0v1)
fit0v1 <- glmFit(y0v1, design_0v1)
lrt0v1 <- glmLRT(fit0v1, coef = 2)

topTags(lrt0v1)
summary(decideTests(lrt0v1))

results0v1 <- as.data.frame(lrt0v1$table)
results0v1$FDR <- p.adjust(results0v1$PValue, method="fdr")

results0v1$deexpressed <- ifelse(results0v1$FDR <= 0.1 & results0v1$logFC >= 0, "UP", 
                                 ifelse(results0v1$FDR <= 0.1 & results0v1$logFC < 0, "DOWN", "NO"))
results0v1$gene <- rownames(results0v1)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(results0v1) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.1)) +
  scale_color_manual(values=mycolors) + 
  theme_classic() + theme(text = element_text(size = 30), legend.position = "none") + 
  xlim(-2, 2)

sig0v1 <- subset(results0v1, FDR <= 0.1)
#write.table(results0v1, file = "Human_DEG_0v1_101824.txt")

##Human - PCA
##Read in counts
hu_counts <- read.delim("human_counts_9724_unstranded.txt", comment.char="#")
hu_counts_countsonly <- hu_counts[,c(7:118)]
rownames(hu_counts_countsonly) <- hu_counts$Geneid

exclude <- c("PQRC21_046_D0_2HR", "PQRC21_021_R3_D0", "PQRC21_021_R3_D0_2hr", "PQRC21_009_D0", "PQRC21_008_D0_8HR_B")

hu_counts_countsonly <- hu_counts_countsonly[,!(colnames(hu_counts_countsonly) %in% exclude)]

#Normalize count data and remove lowly expressed genes
hu_counts_countsonly[] <- sapply(hu_counts_countsonly, as.numeric)
hu_counts_countsonly <- na.omit(hu_counts_countsonly)
dglist_hu <- DGEList(counts=hu_counts_countsonly, genes=rownames(hu_counts_countsonly))

keep <- rowSums(cpm(dglist_hu)>10) >= (ncol(dglist_hu$counts)/2)
table(keep)
dglist_hu_filter <- dglist_hu[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_hu_filter))
dim(dglist_hu_filter)

dglist_hu_filter <- calcNormFactors(dglist_hu_filter, method="TMM")

#Read in meta data
pheno <- read.delim("AR_phenotypes.txt")
rownames(pheno) <- pheno$Sample
pheno$alt_phenotype <- ifelse(pheno$Phenotype == "Average", "Average", "Resistant")

#Human PCA 
count <- as.data.frame(cpm(dglist_hu_filter))
count <- as.data.frame(t(count))
pca <- prcomp(count)
pcs <- as.data.frame(pca$x)

#Plot PCA with ggplot
pcs <- merge(pheno, pcs, by = 0)
ggplot(pcs, aes(x = PC1, y = PC2, color = Phenotype)) + 
  geom_point() +
  theme_minimal() + 
  ggtitle("Human")

## Susceptible
# WT 0v1 hour
#Differential expression of paired samples at time 0 vs. time 1 hour
pheno_WT <- subset(pheno, Phenotype == "Average")

paired_subset_0v1 <- pheno_WT %>%
  group_by(ID_Deidentified) %>%
  filter(hour == 0 | hour == 1) 
pairs <- paired_subset_0v1 %>% group_by(ID) %>% count(ID) %>% filter(n == 2)

paired_subset_0v1 <- paired_subset_0v1[paired_subset_0v1$ID %in% pairs$ID,]

hu_counts <- hu_counts_countsonly[,colnames(hu_counts_countsonly) %in% paired_subset_0v1$Sample]

dglist_hu_0v1 <- DGEList(counts=hu_counts, genes=rownames(hu_counts))

keep <- rowSums(cpm(dglist_hu_0v1)>10) >= (ncol(dglist_hu_0v1$counts)/2)
table(keep)
dglist_hu_filter_0v1 <- dglist_hu_0v1[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_hu_filter_0v1))
dim(dglist_hu_filter_0v1)

dglist_hu_filter_0v1 <- calcNormFactors(dglist_hu_filter_0v1, method="TMM")

paired_subset_0v1 <- paired_subset_0v1[order(match(paired_subset_0v1$Sample, colnames(dglist_hu_filter_0v1))),]

celltypes <- read.delim("human_cibersort_101824.txt")
celltypes <- celltypes[celltypes$Mixture %in% colnames(dglist_hu_filter_0v1),]
celltypes <- celltypes[ order(match(celltypes$Mixture, colnames(dglist_hu_filter_0v1))), ]

id <- as.factor(paired_subset_0v1$ID)
hour <- as.factor(paired_subset_0v1$hour)

design_0v1 <- model.matrix(~id + hour)
rownames(design_0v1) <- colnames(dglist_hu_filter_0v1)
design_0v1

y0v1 <- estimateDisp(dglist_hu_filter_0v1, design_0v1)
fit0v1 <- glmFit(y0v1, design_0v1)
lrt0v1 <- glmLRT(fit0v1, coef = 15)

topTags(lrt0v1)
summary(decideTests(lrt0v1))

results0v1 <- as.data.frame(lrt0v1$table)
results0v1$FDR <- p.adjust(results0v1$PValue, method="fdr")

results0v1$deexpressed <- ifelse(results0v1$FDR <= 0.1 & results0v1$logFC >= 0, "UP", 
                                 ifelse(results0v1$FDR <= 0.1 & results0v1$logFC < 0, "DOWN", "NO"))
results0v1$gene <- rownames(results0v1)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

results0v1 <- subset(results0v1, gene != "PVP01_0202700" & gene != "PVP01_0202800")

ggplot(results0v1) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.1)) +
  scale_color_manual(values=mycolors) + 
  theme_classic() + theme(text = element_text(size = 30), legend.position = "none") + 
  xlim(-2, 2)


sig0v1 <- subset(results0v1, FDR <= 0.1)

celltypes$B <- celltypes$B.cells.memory + celltypes$B.cells.naive
celltypes$T <- celltypes$T.cells.CD4.memory.activated + celltypes$T.cells.CD4.memory.resting + celltypes$T.cells.CD4.naive + celltypes$T.cells.CD8 + celltypes$T.cells.follicular.helper + celltypes$T.cells.gamma.delta + celltypes$T.cells.regulatory..Tregs.
celltypes$APC <- celltypes$Macrophages.M0 + celltypes$Macrophages.M1 + celltypes$Macrophages.M2 + celltypes$Dendritic.cells.activated + celltypes$Dendritic.cells.resting + celltypes$Monocytes
celltypes$mast <- celltypes$Mast.cells.activated + celltypes$Mast.cells.resting
celltypes$NK <- celltypes$NK.cells.activated + celltypes$NK.cells.resting
celltypes$Sample <- celltypes$Mixture

df2 <- merge(celltypes, paired_subset_0v1, by = "Sample")
ggplot(df2, aes(x = as.factor(hour), y = Eosinophils, fill = as.factor(hour))) + 
  geom_boxplot() + 
  theme_classic() + 
  xlab("Hour") + ylab("Proportion of Neutrophils") +
  theme(text = element_text(size = 20), legend.position = "none") +
  stat_compare_means(method = "t.test", paired = TRUE, label.y = 1, size = 8)

##Differential expression of paired samples at time 1 hour vs time 2 hours
paired_subset_1v2 <- pheno_WT %>%
  group_by(ID_Deidentified) %>%
  filter(hour == 1 | hour == 2) 
pairs <- paired_subset_1v2 %>% group_by(ID) %>% count(ID) %>% filter(n == 2)

paired_subset_1v2 <- paired_subset_1v2[paired_subset_1v2$ID %in% pairs$ID,]

hu_counts <- hu_counts_countsonly[,colnames(hu_counts_countsonly) %in% paired_subset_1v2$Sample]
paired_subset_1v2 <- paired_subset_1v2[paired_subset_1v2$Sample %in% colnames(hu_counts),]

dglist_hu_1v2 <- DGEList(counts=hu_counts, genes=rownames(hu_counts))

keep <- rowSums(cpm(dglist_hu_1v2)>10) >= (ncol(dglist_hu_1v2$counts)/2)
table(keep)
dglist_hu_filter_1v2 <- dglist_hu_1v2[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_hu_filter_1v2))
dim(dglist_hu_filter_1v2)

dglist_hu_filter_1v2 <- calcNormFactors(dglist_hu_filter_1v2, method="TMM")

paired_subset_1v2 <- paired_subset_1v2[order(match(paired_subset_1v2$Sample, colnames(dglist_hu_filter_1v2))),]

id <- as.factor(paired_subset_1v2$ID)
hour <- as.factor(paired_subset_1v2$hour)

design_1v2 <- model.matrix(~id + hour)
rownames(design_1v2) <- colnames(dglist_hu_filter_1v2)
design_1v2

y1v2 <- estimateDisp(dglist_hu_filter_1v2, design_1v2)
fit1v2 <- glmFit(y1v2, design_1v2)
lrt1v2 <- glmLRT(fit1v2, coef = 15)

topTags(lrt1v2)
summary(decideTests(lrt1v2))

results1v2 <- as.data.frame(lrt1v2$table)
results1v2$FDR <- p.adjust(results1v2$PValue, method="fdr")

results1v2$deexpressed <- ifelse(results1v2$FDR <= 0.1 & results1v2$logFC >= 0, "UP", 
                                 ifelse(results1v2$FDR <= 0.1 & results1v2$logFC < 0, "DOWN", "NO"))
results1v2$gene <- rownames(results1v2)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(results1v2) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.1)) +
  scale_color_manual(values=mycolors) + 
  theme_classic() + theme(text = element_text(size = 30), legend.position = "none")

sig1v2 <- subset(results1v2, FDR <= 0.1)

##Differential expression of paired samples at time 4 hours vs time 2 hours
paired_subset_2v4 <- pheno_WT %>%
  group_by(ID_Deidentified) %>%
  filter(hour == 2 | hour == 4) 
pairs <- paired_subset_2v4 %>% group_by(ID) %>% count(ID) %>% filter(n == 2)

paired_subset_2v4 <- paired_subset_2v4[paired_subset_2v4$ID %in% pairs$ID,]
paired_subset_2v4 <- subset(paired_subset_2v4, ID_Deidentified != "M")

hu_counts <- hu_counts_countsonly[,colnames(hu_counts_countsonly) %in% paired_subset_2v4$Sample]
paired_subset_2v4 <- paired_subset_2v4[paired_subset_2v4$Sample %in% colnames(hu_counts),]


dglist_hu_2v4 <- DGEList(counts=hu_counts, genes=rownames(hu_counts))

keep <- rowSums(cpm(dglist_hu_2v4)>10) >= (ncol(dglist_hu_2v4$counts)/2)
table(keep)
dglist_hu_filter_2v4 <- dglist_hu_2v4[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_hu_filter_2v4))
dim(dglist_hu_filter_2v4)

dglist_hu_filter_2v4 <- calcNormFactors(dglist_hu_filter_2v4, method="TMM")

paired_subset_2v4 <- paired_subset_2v4[order(match(paired_subset_2v4$Sample, colnames(dglist_hu_filter_2v4))),]

id <- as.factor(paired_subset_2v4$ID)
hour <- as.factor(paired_subset_2v4$hour)

design_2v4 <- model.matrix(~id + hour)
rownames(design_2v4) <- colnames(dglist_hu_filter_2v4)
design_2v4

y2v4 <- estimateDisp(dglist_hu_filter_2v4, design_2v4)
fit2v4 <- glmFit(y2v4, design_2v4)
lrt2v4 <- glmLRT(fit2v4, coef = 12)

topTags(lrt2v4)
summary(decideTests(lrt2v4))

results2v4 <- as.data.frame(lrt2v4$table)
results2v4$FDR <- p.adjust(results2v4$PValue, method="fdr")

results2v4$deexpressed <- ifelse(results2v4$FDR <= 0.1 & results2v4$logFC >= 0, "UP", 
                                 ifelse(results2v4$FDR <= 0.1 & results2v4$logFC < 0, "DOWN", "NO"))
results2v4$gene <- rownames(results2v4)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(results2v4) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.1)) +
  scale_color_manual(values=mycolors) + 
  theme_classic() + theme(text = element_text(size = 30), legend.position = "none")

sig2v4 <- subset(results2v4, FDR <= 0.1)

######## Re-do analyses with resistant and slow parasites
#Differential expression of paired samples at time 0 vs. time 1 hour
pheno_resistant <- subset(pheno, Phenotype == "Resistant" | Phenotype == "Slow")

paired_subset_0v1 <- pheno_resistant %>%
  group_by(ID_Deidentified) %>%
  filter(hour == 0 | hour == 1) 
pairs <- paired_subset_0v1 %>% group_by(ID) %>% count(ID) %>% filter(n == 2)

paired_subset_0v1 <- paired_subset_0v1[paired_subset_0v1$ID %in% pairs$ID,]
paired_subset_0v1 <- subset(paired_subset_0v1, ID_Deidentified != "V")

paired_subset_0v1 <- paired_subset_0v1[paired_subset_0v1$ID %in% pairs$ID,]

hu_counts <- hu_counts_countsonly[,colnames(hu_counts_countsonly) %in% paired_subset_0v1$Sample]

dglist_hu_0v1 <- DGEList(counts=hu_counts, genes=rownames(hu_counts))

keep <- rowSums(cpm(dglist_hu_0v1)>10) >= (ncol(dglist_hu_0v1$counts)/2)
table(keep)
dglist_hu_filter_0v1 <- dglist_hu_0v1[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_hu_filter_0v1))
dim(dglist_hu_filter_0v1)

dglist_hu_filter_0v1 <- calcNormFactors(dglist_hu_filter_0v1, method="TMM")

paired_subset_0v1 <- paired_subset_0v1[order(match(paired_subset_0v1$Sample, colnames(dglist_hu_filter_0v1))),]

id <- as.factor(paired_subset_0v1$ID)
hour <- as.factor(paired_subset_0v1$hour)

design_0v1 <- model.matrix(~id + hour)
rownames(design_0v1) <- colnames(dglist_hu_filter_0v1)
design_0v1

y0v1 <- estimateDisp(dglist_hu_filter_0v1, design_0v1)
fit0v1 <- glmFit(y0v1, design_0v1)
lrt0v1 <- glmLRT(fit0v1, coef = 13)

topTags(lrt0v1)
summary(decideTests(lrt0v1))

results0v1 <- as.data.frame(lrt0v1$table)
results0v1$FDR <- p.adjust(results0v1$PValue, method="fdr")

results0v1$deexpressed <- ifelse(results0v1$FDR <= 0.1 & results0v1$logFC >= 0, "UP", 
                                 ifelse(results0v1$FDR <= 0.1 & results0v1$logFC < 0, "DOWN", "NO"))
results0v1$gene <- rownames(results0v1)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

results0v1 <- subset(results0v1, gene != "PVP01_0202700" & gene != "PVP01_0202800")

ggplot(results0v1) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.1)) +
  scale_color_manual(values=mycolors) + 
  theme_classic() + theme(text = element_text(size = 30), legend.position = "none") + 
  xlim(-1, 1)


sig0v1 <- subset(results0v1, FDR <= 0.1)

##Differential expression of paired samples at time 1 hour vs time 2 hours
paired_subset_1v2 <- pheno_resistant %>%
  group_by(ID_Deidentified) %>%
  filter(hour == 1 | hour == 2) 
pairs <- paired_subset_1v2 %>% group_by(ID) %>% count(ID) %>% filter(n == 2)

paired_subset_1v2 <- paired_subset_1v2[paired_subset_1v2$ID %in% pairs$ID,]
paired_subset_1v2 <- subset(paired_subset_1v2, ID_Deidentified != "V")

hu_counts <- hu_counts_countsonly[,colnames(hu_counts_countsonly) %in% paired_subset_1v2$Sample]
paired_subset_1v2 <- paired_subset_1v2[paired_subset_1v2$Sample %in% colnames(hu_counts),]

dglist_hu_1v2 <- DGEList(counts=hu_counts, genes=rownames(hu_counts))

keep <- rowSums(cpm(dglist_hu_1v2)>10) >= (ncol(dglist_hu_1v2$counts)/2)
table(keep)
dglist_hu_filter_1v2 <- dglist_hu_1v2[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_hu_filter_1v2))
dim(dglist_hu_filter_1v2)

dglist_hu_filter_1v2 <- calcNormFactors(dglist_hu_filter_1v2, method="TMM")

paired_subset_1v2 <- paired_subset_1v2[order(match(paired_subset_1v2$Sample, colnames(dglist_hu_filter_1v2))),]

id <- as.factor(paired_subset_1v2$ID)
hour <- as.factor(paired_subset_1v2$hour)

design_1v2 <- model.matrix(~id + hour)
rownames(design_1v2) <- colnames(dglist_hu_filter_1v2)
design_1v2

y1v2 <- estimateDisp(dglist_hu_filter_1v2, design_1v2)
fit1v2 <- glmFit(y1v2, design_1v2)
lrt1v2 <- glmLRT(fit1v2, coef = 11)

topTags(lrt1v2)
summary(decideTests(lrt1v2))

results1v2 <- as.data.frame(lrt1v2$table)
results1v2$FDR <- p.adjust(results1v2$PValue, method="fdr")

results1v2$deexpressed <- ifelse(results1v2$FDR <= 0.1 & results1v2$logFC >= 0, "UP", 
                                 ifelse(results1v2$FDR <= 0.1 & results1v2$logFC < 0, "DOWN", "NO"))
results1v2$gene <- rownames(results1v2)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(results1v2) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.1)) +
  scale_color_manual(values=mycolors) + 
  theme_classic() + theme(text = element_text(size = 30), legend.position = "none")

sig1v2 <- subset(results1v2, FDR <= 0.1)

##Differential expression of paired samples at time 4 hours vs time 2 hours
paired_subset_2v4 <- pheno_resistant %>%
  group_by(ID_Deidentified) %>%
  filter(hour == 2 | hour == 4) 
pairs <- paired_subset_2v4 %>% group_by(ID) %>% count(ID) %>% filter(n == 2)

paired_subset_2v4 <- paired_subset_2v4[paired_subset_2v4$ID %in% pairs$ID,]

hu_counts <- hu_counts_countsonly[,colnames(hu_counts_countsonly) %in% paired_subset_2v4$Sample]
paired_subset_2v4 <- paired_subset_2v4[paired_subset_2v4$Sample %in% colnames(hu_counts),]


dglist_hu_2v4 <- DGEList(counts=hu_counts, genes=rownames(hu_counts))

keep <- rowSums(cpm(dglist_hu_2v4)>10) >= (ncol(dglist_hu_2v4$counts)/2)
table(keep)
dglist_hu_filter_2v4 <- dglist_hu_2v4[keep, , keep.lib.sizes=FALSE]
summary(cpm(dglist_hu_filter_2v4))
dim(dglist_hu_filter_2v4)

dglist_hu_filter_2v4 <- calcNormFactors(dglist_hu_filter_2v4, method="TMM")

paired_subset_2v4 <- paired_subset_2v4[order(match(paired_subset_2v4$Sample, colnames(dglist_hu_filter_2v4))),]

id <- as.factor(paired_subset_2v4$ID)
hour <- as.factor(paired_subset_2v4$hour)

design_2v4 <- model.matrix(~id + hour)
rownames(design_2v4) <- colnames(dglist_hu_filter_2v4)
design_2v4

y2v4 <- estimateDisp(dglist_hu_filter_2v4, design_2v4)
fit2v4 <- glmFit(y2v4, design_2v4)
lrt2v4 <- glmLRT(fit2v4, coef = 11)

topTags(lrt2v4)
summary(decideTests(lrt2v4))

results2v4 <- as.data.frame(lrt2v4$table)
results2v4$FDR <- p.adjust(results2v4$PValue, method="fdr")

results2v4$deexpressed <- ifelse(results2v4$FDR <= 0.1 & results2v4$logFC >= 0, "UP", 
                                 ifelse(results2v4$FDR <= 0.1 & results2v4$logFC < 0, "DOWN", "NO"))
results2v4$gene <- rownames(results2v4)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(results2v4) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.1)) +
  scale_color_manual(values=mycolors) + 
  theme_classic() + theme(text = element_text(size = 30), legend.position = "none")

sig2v4 <- subset(results2v4, FDR <= 0.1)
