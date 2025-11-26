library(stringr)
library(car)

SAMPLES_TO_EXCLUDE <- c("X90_M3")
N <- 14

bf <- 0.05 / N

stages <- c('F1', 'F2', 'F3', 'F4', 'F5')
#stages <- c('M1', 'M2', 'M3', 'M4')
expression_file <- "DataChapterTwo/DESeq/NormalisedData/normalised_counts_gonad_length_scaled.tsv"

df <- read.csv(expression_file, sep = '\t', header = FALSE)
df <- t(df)
colnames(df) <- df[1,]
df <- df[2:nrow(df),]
row.names(df) <- NULL
df <- as.data.frame(df)
df <- subset(df, !df$gene_name %in% SAMPLES_TO_EXCLUDE)
samples <- df$gene_name
df <- as.data.frame(lapply(df, as.numeric))

df$stage <- substring(str_extract(samples, pattern = '_[A-Z][1-9]'), 2, 3)

df <- subset(df, df$stage %in% stages)
df$stage <- as.factor(df$stage)

leveneTest(log(star + 1) ~ stage, data = df) # No
model <- aov(log(star + 1) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(cyp11a1.2) ~ stage, data = df) # Yes
model <- aov(log(cyp11a1.2) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(cyp17a1) ~ stage, data = df) # Yes
model <- aov(log(cyp17a1) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(hsd3b1) ~ stage, data = df) # Yes
model <- aov(log(hsd3b1) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(hsd17b1) ~ stage, data = df) # Yes
model <- aov(log(hsd17b1) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)
# =====================================================
leveneTest(log(cyp19a1a) ~ stage, data = df) # No
model <- aov(log(cyp19a1a) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(foxl2a) ~ stage, data = df) # Yes
model <- aov(log(foxl2a) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(cyp11b1) ~ stage, data = df) # No
model <- aov(log(cyp11b1) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(hsd11b2 + 1) ~ stage, data = df) # No
model <- aov(log(hsd11b2 + 1) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)
# =====================================================
leveneTest(log(fshr) ~ stage, data = df) # 
model <- aov(log(fshr) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(lhcgr_ChrAur1) ~ stage, data = df) # 
model <- aov(log(lhcgr_ChrAur1) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(ar) ~ stage, data = df) # 
model <- aov(log(ar) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(armt1_fusion_2) ~ stage, data = df) # 
model <- aov(log(armt1_fusion_2) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(esr2a ~ stage, data = df) # 
model <- aov(esr2a ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

# leveneTest(esr2b ~ stage, data = df) # 
# model <- aov(esr2b ~ stage, data = df)
# shapiro.test(resid(model))
# summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
# plot(model, which = 1)
# plot(model, which = 2)


stages <- c('M3', 'M4', 'M2')
#stages <- c('M1', 'M2', 'M3', 'M4')
expression_file <- "C:/Users/cfnsjm/Local Doccuments/PhD/Programming/R/DESeq_final/NormalisedData/normalised_counts_gonad_length_scaled.tsv"

df <- read.csv(expression_file, sep = '\t', header = FALSE)
df <- t(df)
colnames(df) <- df[1,]
df <- df[2:nrow(df),]
row.names(df) <- NULL
df <- as.data.frame(df)
df <- subset(df, !df$gene_name %in% SAMPLES_TO_EXCLUDE)
samples <- df$gene_name
df <- as.data.frame(lapply(df, as.numeric))

df$stage <- substring(str_extract(samples, pattern = '_[A-Z][1-9]'), 2, 3)

df <- subset(df, df$stage %in% stages)
df$stage <- as.factor(df$stage)

leveneTest(log(star) ~ stage, data = df) # Yes
model <- aov(log(star) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(cyp11a1.2) ~ stage, data = df) # Yes
model <- aov(log(cyp11a1.2) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(cyp17a1) ~ stage, data = df) # Yes
model <- aov(log(cyp17a1) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(hsd3b1) ~ stage, data = df) # Yes
model <- aov(log(hsd3b1) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(hsd17b1) ~ stage, data = df) # Yes
model <- aov(log(hsd17b1) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)
# ======================================================
leveneTest(cyp19a1a ~ stage, data = df) # No
model <- aov(cyp19a1a ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(foxl2a + 1) ~ stage, data = df) # No
model <- aov(log(foxl2a + 1) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(cyp11b1) ~ stage, data = df) # Yes
model <- aov(log(cyp11b1) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(hsd11b2 ~ stage, data = df) # Yes
model <- aov(hsd11b2 ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)
# =====================================================
leveneTest(log(fshr) ~ stage, data = df) # 
model <- aov(log(fshr) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(lhcgr_ChrAur1) ~ stage, data = df) # 
model <- aov(log(lhcgr_ChrAur1) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(ar) ~ stage, data = df) # 
model <- aov(log(ar) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(armt1_fusion_2) ~ stage, data = df) # 
model <- aov(log(armt1_fusion_2) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(esr2a ~ stage, data = df) # 
model <- aov(esr2a ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)

leveneTest(log(esr2b) ~ stage, data = df) # 
model <- aov(log(esr2b) ~ stage, data = df)
shapiro.test(resid(model))
summary(model); summary(model)[[1]][["Pr(>F)"]][1] < bf
plot(model, which = 1)
plot(model, which = 2)
