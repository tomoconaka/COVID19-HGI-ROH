setwd("/home/tomoko/scratch/ROH/")
data <- readRDS("pheno_roh.rds")

data %>% filter(B2 != -1) %>% group_by(B2) %>% 
  summarise(median=median(Froh),
            lower_IQR = quantile(Froh)[2],
            higher_IQR = quantile(Froh)[4])

data %>% group_by(study) %>% 
  summarise(median=median(Froh),
            lower_IQR = quantile(Froh)[2],
            higher_IQR = quantile(Froh)[4])

df_roh_pheno <- data %>% select(A2, B2, Froh, Fhat1, Fhat3, age_at_diagnosis, sex, nation, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
df_roh_pheno <- df_roh_pheno %>% mutate_at(.vars = c("A2", "B2"), .funs = funs(ifelse(.==-1,NA,.))) %>% drop_na(c(Froh, Fhat1, Fhat3, age_at_diagnosis, sex, nation, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))


#################################################
#INBREEDING DEPRESSION ANALYSIS ALL TOGETHER -------
df_roh_pheno <- df_roh_pheno
mod_1<-glm(B2~Froh+Fhat1+age_at_diagnosis+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "ALL:Froh+Fhat1"
res$outcome <- "B2"
res$Ncase <- dim(df_roh_pheno %>% filter(B2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(B2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = T, row.names = F, sep="\t")

mod_1<-glm(A2~Froh+Fhat1+age_at_diagnosis+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "ALL:Froh+Fhat1"
res$outcome <- "A2"
res$Ncase <- dim(df_roh_pheno %>% filter(A2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(A2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

mod_1<-glm(B2~Froh+Fhat3+age_at_diagnosis+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "ALL:Froh+Fhat3"
res$outcome <- "B2"
res$Ncase <- dim(df_roh_pheno %>% filter(B2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(B2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

mod_1<-glm(A2~Froh+Fhat3+age_at_diagnosis+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "ALL:Froh+Fhat3"
res$outcome <- "A2"
res$Ncase <- dim(df_roh_pheno %>% filter(A2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(A2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

#INBREEDING DEPRESSION ANALYSIS MEN -------
df_roh_pheno <- data %>% select(A2, B2, Froh, Fhat1, Fhat3, age_at_diagnosis, sex, nation, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
df_roh_pheno <- df_roh_pheno %>% mutate_at(.vars = c("A2", "B2"), .funs = funs(ifelse(.==-1,NA,.))) %>% drop_na(c(Froh, Fhat1, Fhat3, age_at_diagnosis, sex, nation, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))
df_roh_pheno <- df_roh_pheno %>% filter(sex == 0)

mod_1<-glm(B2~Froh+Fhat1+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "MEN:Froh+Fhat1"
res$outcome <- "B2"
res$Ncase <- dim(df_roh_pheno %>% filter(B2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(B2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]
write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")


mod_1<-glm(A2~Froh+Fhat1+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "MEN:Froh+Fhat1"
res$outcome <- "A2"
res$Ncase <- dim(df_roh_pheno %>% filter(A2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(A2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

mod_1<-glm(B2~Froh+Fhat3+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "MEN:Froh+Fhat3"
res$outcome <- "B2"
res$Ncase <- dim(df_roh_pheno %>% filter(B2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(B2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

mod_1<-glm(A2~Froh+Fhat3+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "MEN:Froh+Fhat3"
res$outcome <- "A2"
res$Ncase <- dim(df_roh_pheno %>% filter(A2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(A2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")


#INBREEDING DEPRESSION ANALYSIS WOMEN -------
df_roh_pheno <- data %>% select(A2, B2, Froh, Fhat1, Fhat3, age_at_diagnosis, sex, nation, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
df_roh_pheno <- df_roh_pheno %>% mutate_at(.vars = c("A2", "B2"), .funs = funs(ifelse(.==-1,NA,.))) %>% drop_na(c(Froh, Fhat1, Fhat3, age_at_diagnosis, sex, nation, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))
df_roh_pheno <- df_roh_pheno %>% filter(sex == 1)

mod_1<-glm(B2~Froh+Fhat1+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "WOMEN:Froh+Fhat1"
res$outcome <- "B2"
res$Ncase <- dim(df_roh_pheno %>% filter(B2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(B2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")


mod_1<-glm(A2~Froh+Fhat1+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "WOMEN:Froh+Fhat1"
res$outcome <- "A2"
res$Ncase <- dim(df_roh_pheno %>% filter(A2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(A2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

mod_1<-glm(B2~Froh+Fhat3+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "WOMEN:Froh+Fhat3"
res$outcome <- "B2"
res$Ncase <- dim(df_roh_pheno %>% filter(B2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(B2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

mod_1<-glm(A2~Froh+Fhat3+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "WOMEN:Froh+Fhat3"
res$outcome <- "A2"
res$Ncase <- dim(df_roh_pheno %>% filter(A2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(A2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

#INBREEDING DEPRESSION ANALYSIS YOUNG MEN -------
df_roh_pheno <- data %>% select(A2, B2, Froh, Fhat1, Fhat3, age_at_diagnosis, sex, nation, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
df_roh_pheno <- df_roh_pheno %>% mutate_at(.vars = c("A2", "B2"), .funs = funs(ifelse(.==-1,NA,.))) %>% drop_na(c(Froh, Fhat1, Fhat3, age_at_diagnosis, sex, nation, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))
df_roh_pheno <- df_roh_pheno %>% filter(sex == 0 & age_at_diagnosis <= 60) 

mod_1<-glm(B2~Froh+Fhat1+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "MEN≤60:Froh+Fhat1"
res$outcome <- "B2"
res$Ncase <- dim(df_roh_pheno %>% filter(B2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(B2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")


mod_1<-glm(A2~Froh+Fhat1+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "MEN≤60:Froh+Fhat1"
res$outcome <- "A2"
res$Ncase <- dim(df_roh_pheno %>% filter(A2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(A2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

mod_1<-glm(B2~Froh+Fhat3+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "MEN≤60:Froh+Fhat3"
res$outcome <- "B2"
res$Ncase <- dim(df_roh_pheno %>% filter(B2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(B2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

mod_1<-glm(A2~Froh+Fhat3+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "MEN≤60:Froh+Fhat3"
res$outcome <- "A2"
res$Ncase <- dim(df_roh_pheno %>% filter(A2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(A2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

#INBREEDING DEPRESSION ANALYSIS OLD MEN -------
df_roh_pheno <- data %>% select(A2, B2, Froh, Fhat1, Fhat3, age_at_diagnosis, sex, nation, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
df_roh_pheno <- df_roh_pheno %>% mutate_at(.vars = c("A2", "B2"), .funs = funs(ifelse(.==-1,NA,.))) %>% drop_na(c(Froh, Fhat1, Fhat3, age_at_diagnosis, sex, nation, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))
df_roh_pheno <- df_roh_pheno %>% filter(sex == 0 & age_at_diagnosis > 60) 

mod_1<-glm(B2~Froh+Fhat1+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "MEN>60:Froh+Fhat1"
res$outcome <- "B2"
res$Ncase <- dim(df_roh_pheno %>% filter(B2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(B2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

mod_1<-glm(A2~Froh+Fhat1+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "MEN>60:Froh+Fhat1"
res$outcome <- "A2"
res$Ncase <- dim(df_roh_pheno %>% filter(A2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(A2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

mod_1<-glm(B2~Froh+Fhat3+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "MEN>60:Froh+Fhat3"
res$outcome <- "B2"
res$Ncase <- dim(df_roh_pheno %>% filter(B2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(B2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

mod_1<-glm(A2~Froh+Fhat3+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "MEN>60:Froh+Fhat3"
res$outcome <- "A2"
res$Ncase <- dim(df_roh_pheno %>% filter(A2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(A2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

#INBREEDING DEPRESSION ANALYSIS YOUNG WOMEN -------
df_roh_pheno <- data %>% select(A2, B2, Froh, Fhat1, Fhat3, age_at_diagnosis, sex, nation, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
df_roh_pheno <- df_roh_pheno %>% mutate_at(.vars = c("A2", "B2"), .funs = funs(ifelse(.==-1,NA,.))) %>% drop_na(c(Froh, Fhat1, Fhat3, age_at_diagnosis, sex, nation, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))
df_roh_pheno <- df_roh_pheno %>% filter(sex == 1 & age_at_diagnosis <= 60) 

mod_1<-glm(B2~Froh+Fhat1+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "WOMEN≤60:Froh+Fhat1"
res$outcome <- "B2"
res$Ncase <- dim(df_roh_pheno %>% filter(B2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(B2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")


mod_1<-glm(A2~Froh+Fhat1+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "WOMEN≤60:Froh+Fhat1"
res$outcome <- "A2"
res$Ncase <- dim(df_roh_pheno %>% filter(A2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(A2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T,sep="\t")

mod_1<-glm(B2~Froh+Fhat3+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "WOMEN≤60:Froh+Fhat3"
res$outcome <- "B2"
res$Ncase <- dim(df_roh_pheno %>% filter(B2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(B2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T,sep="\t")

mod_1<-glm(A2~Froh+Fhat3+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "WOMEN≤60:Froh+Fhat3"
res$outcome <- "A2"
res$Ncase <- dim(df_roh_pheno %>% filter(A2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(A2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

#INBREEDING DEPRESSION ANALYSIS OLD WOMEN -------
df_roh_pheno <- data %>% select(A2, B2, Froh, Fhat1, Fhat3, age_at_diagnosis, sex, nation, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
df_roh_pheno <- df_roh_pheno %>% mutate_at(.vars = c("A2", "B2"), .funs = funs(ifelse(.==-1,NA,.))) %>% drop_na(c(Froh, Fhat1, Fhat3, age_at_diagnosis, sex, nation, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))
df_roh_pheno <- df_roh_pheno %>% filter(sex == 1 & age_at_diagnosis > 60) 

mod_1<-glm(B2~Froh+Fhat1+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "WOMEN>60:Froh+Fhat1"
res$outcome <- "B2"
res$Ncase <- dim(df_roh_pheno %>% filter(B2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(B2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

mod_1<-glm(A2~Froh+Fhat1+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "WOMEN>60:Froh+Fhat1"
res$outcome <- "A2"
res$Ncase <- dim(df_roh_pheno %>% filter(A2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(A2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

mod_1<-glm(B2~Froh+Fhat3+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "WOMEN>60:Froh+Fhat3"
res$outcome <- "B2"
res$Ncase <- dim(df_roh_pheno %>% filter(B2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(B2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")

mod_1<-glm(A2~Froh+Fhat3+age_at_diagnosis+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+nation, dat = df_roh_pheno, family = binomial);summary(mod_1)
mod_coef<-coef(summary(mod_1)); mod_coef<-mod_coef[-1,]
res <- data.frame(matrix(0, dim(mod_coef)[1], 9))
colnames(res) <- c("model","outcome","Ncase", "Ncontrol","covariates", "Estimate", "ci_min","ci_max","P_val")
res$model <- "WOMEN>60:Froh+Fhat3"
res$outcome <- "A2"
res$Ncase <- dim(df_roh_pheno %>% filter(A2 == 1))[1]
res$Ncontrol <- dim(df_roh_pheno %>% filter(A2 == 0))[1]
res$covariates <- rownames(mod_coef)
res$Estimate <- mod_coef[,1]
res$P_val <- mod_coef[,4]
res$ci_min <- res$Estimate + qnorm(0.025)*mod_coef[,2]
res$ci_max <- res$Estimate + qnorm(0.975)*mod_coef[,2]

write.table(res, file="roh_analyses.tsv", quote=F, col.names = F, row.names = F, append=T, sep="\t")
