setwd("/home/tomoko/scratch/ROH/")

phenoloh <- readRDS("pheno_roh.rds")
library(tidyr)
library(dplyr)
phenoloh <- phenoloh %>% mutate(FID = 0)

files_list <- c("EUR.fimm.hom", "EUR.brasil.hom", "EUR.canada.hom", 
                "EUR.kiel.hom", "EUR.spain_alarcon.hom", "EUR.spain_planas.hom", "EUR.germany_ludwig.hom"
)

dat <- data.frame()                          
for (i in 1:length(files_list)) {
  dat <- rbind(dat, read.table((files_list[i]),header=TRUE))
}

library(data.table)
dat <- as.data.table(dat)
dat$IID<-as.factor(dat$IID)
setkey(dat,"IID")

dat <- dat %>% rename(Start = POS1, End = POS2) %>% 
  filter(KB > 500)
getloh <- dat %>% rename(Sample = IID, Chromosome = CHR) %>% 
  filter(KB > 500) %>%
  select(Sample, Chromosome, Start, End)

sum(table(getloh$Sample)>0)
max(table(getloh$Sample))
hist(table(getloh$Sample), xlab="Number of ROH per individual", br=50, main="")

library(GenomicRanges)
library(arm)
library(parallel)

phenoname <- "death"

gwasloh <- data.frame()
for(chr in 1:22){
  ss <- which(getloh[,2] == chr)
  grchr <- GRanges(paste("chr", chr, sep=""), IRanges(start=getloh$Start[ss], end=getloh$End[ss]),
                   sub = getloh[ss]$Sample)
  mn <- min(start(grchr))
  mx <- max(start(grchr))
  blocks <- round(seq(mn,mx, by=500000))
  for(bb in c(1:(length(blocks)-1))){
    block <- GRanges(paste("chr", chr, sep=""), IRanges(start=blocks[bb],
                                                        end=blocks[bb+1]))
    selover <- data.frame(findOverlaps(block,grchr))[,2]
    selsubsloh <- unique(as.character(grchr$sub[selover]))
    dat <- phenoloh
    dat$loh <- rep(0,nrow(dat))
    if(length(selsubsloh) > 0){
      dat <- dat %>% mutate(loh = ifelse(IID %in% selsubsloh, 1, 0))
      dat$pheno <- dat[,phenoname]
      dat <- dat %>% mutate(pheno = ifelse(pheno==-1, NA, pheno))
      assocblock <- summary(bayesglm(pheno ~ loh+age_at_diagnosis*sex+
                                        PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                      family="binomial", data=dat))$coeff["loh",c(1,4)]
      tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                        OR=exp(assocblock[1]), P=assocblock[2], freq=sum(dat$loh==1))
      assocchr <- bind_rows(assocchr, tmp)
    }
    if(length(selsubsloh) == 0){
      tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                        OR=NA, P=NA, freq=0)
      assocchr <- bind_rows(assocchr, tmp)
    }
  }
  gwasloh <- bind_rows(gwasloh, assocchr)
 }

gwasloh <- unique(gwasloh)

rownames(gwasloh) <- NULL
fr005 <- sum(table(phenoloh[,phenoname]))*0.05 
fr005
## [1] 166.35
tb <- gwasloh[gwasloh$freq > fr005, ]
rownames(tb) <- NULL
tb
gwasloh$phenotype <- "death"

saveRDS(gwasloh, file="gwasloh_death.rds")

phenoname <- "icu_admit"

assocchr <- data.frame()
gwasloh <- data.frame()
for(chr in 1:22){
  ss <- which(getloh[,2] == chr)
  grchr <- GRanges(paste("chr", chr, sep=""), IRanges(start=getloh$Start[ss], end=getloh$End[ss]),
                   sub = getloh[ss]$Sample)
  mn <- min(start(grchr))
  mx <- max(start(grchr))
  blocks <- round(seq(mn,mx, by=500000))
  for(bb in c(1:(length(blocks)-1))){
    block <- GRanges(paste("chr", chr, sep=""), IRanges(start=blocks[bb],
                                                        end=blocks[bb+1]))
    selover <- data.frame(findOverlaps(block,grchr))[,2]
    selsubsloh <- unique(as.character(grchr$sub[selover]))
    dat <- phenoloh
    dat$loh <- rep(0,nrow(dat))
    if(length(selsubsloh) > 0){
      dat <- dat %>% mutate(loh = ifelse(IID %in% selsubsloh, 1, 0))
      dat$pheno <- dat[,phenoname]
      dat <- dat %>% mutate(pheno = ifelse(pheno==-1, NA, pheno))
      assocblock <- summary(bayesglm(pheno ~ loh+age_at_diagnosis*sex+
                                       PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                     family="binomial", data=dat))$coeff["loh",c(1,4)]
      tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                        OR=exp(assocblock[1]), P=assocblock[2], freq=sum(dat$loh==1))
      assocchr <- bind_rows(assocchr, tmp)
    }
    if(length(selsubsloh) == 0){
      tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                        OR=NA, P=NA, freq=0)
      assocchr <- bind_rows(assocchr, tmp)
    }
  }
  gwasloh <- bind_rows(gwasloh, assocchr)
}
gwasloh <- unique(gwasloh)
rownames(gwasloh) <- NULL
fr005 <- sum(table(phenoloh[,phenoname]))*0.05 
fr005
## [1] 166.35
tb <- gwasloh[gwasloh$freq > fr005, ]
rownames(tb) <- NULL
tb
gwasloh$phenotype <- "icu admission"

saveRDS(gwasloh, file="gwasloh_icu.rds")

phenoname <- "A2"

assocchr <- data.frame()
gwasloh <- data.frame()
for(chr in 1:22){
  ss <- which(getloh[,2] == chr)
  grchr <- GRanges(paste("chr", chr, sep=""), IRanges(start=getloh$Start[ss], end=getloh$End[ss]),
                   sub = getloh[ss]$Sample)
  mn <- min(start(grchr))
  mx <- max(start(grchr))
  blocks <- round(seq(mn,mx, by=500000))
  for(bb in c(1:(length(blocks)-1))){
    block <- GRanges(paste("chr", chr, sep=""), IRanges(start=blocks[bb],
                                                        end=blocks[bb+1]))
    selover <- data.frame(findOverlaps(block,grchr))[,2]
    selsubsloh <- unique(as.character(grchr$sub[selover]))
    dat <- phenoloh
    dat$loh <- rep(0,nrow(dat))
    if(length(selsubsloh) > 0){
      dat <- dat %>% mutate(loh = ifelse(IID %in% selsubsloh, 1, 0))
      dat$pheno <- dat[,phenoname]
      dat <- dat %>% mutate(pheno = ifelse(pheno==-1, NA, pheno))
      assocblock <- summary(bayesglm(pheno ~ loh+age_at_diagnosis*sex+
                                       PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                     family="binomial", data=dat))$coeff["loh",c(1,4)]
      tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                        OR=exp(assocblock[1]), P=assocblock[2], freq=sum(dat$loh==1))
      assocchr <- bind_rows(assocchr, tmp)
    }
    if(length(selsubsloh) == 0){
      tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                        OR=NA, P=NA, freq=0)
      assocchr <- bind_rows(assocchr, tmp)
    }
  }
  gwasloh <- bind_rows(gwasloh, assocchr)
}
gwasloh <- unique(gwasloh)
rownames(gwasloh) <- NULL
fr005 <- sum(table(phenoloh[,phenoname]))*0.05 
fr005
## [1] 166.35
tb <- gwasloh[gwasloh$freq > fr005, ]
rownames(tb) <- NULL
tb
gwasloh$phenotype <- "A2"

saveRDS(gwasloh, file="gwasloh_a2.rds")


phenoname <- "B2"

assocchr <- data.frame()
gwasloh <- data.frame()
for(chr in 1:22){
  ss <- which(getloh[,2] == chr)
  grchr <- GRanges(paste("chr", chr, sep=""), IRanges(start=getloh$Start[ss], end=getloh$End[ss]),
                   sub = getloh[ss]$Sample)
  mn <- min(start(grchr))
  mx <- max(start(grchr))
  blocks <- round(seq(mn,mx, by=500000))
  for(bb in c(1:(length(blocks)-1))){
    block <- GRanges(paste("chr", chr, sep=""), IRanges(start=blocks[bb],
                                                        end=blocks[bb+1]))
    selover <- data.frame(findOverlaps(block,grchr))[,2]
    selsubsloh <- unique(as.character(grchr$sub[selover]))
    dat <- phenoloh
    dat$loh <- rep(0,nrow(dat))
    if(length(selsubsloh) > 0){
      dat <- dat %>% mutate(loh = ifelse(IID %in% selsubsloh, 1, 0))
      dat$pheno <- dat[,phenoname]
      dat <- dat %>% mutate(pheno = ifelse(pheno==-1, NA, pheno))
      assocblock <- summary(bayesglm(pheno ~ loh+age_at_diagnosis*sex+
                                       PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                     family="binomial", data=dat))$coeff["loh",c(1,4)]
      tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                        OR=exp(assocblock[1]), P=assocblock[2], freq=sum(dat$loh==1))
      assocchr <- bind_rows(assocchr, tmp)
    }
    if(length(selsubsloh) == 0){
      tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                        OR=NA, P=NA, freq=0)
      assocchr <- bind_rows(assocchr, tmp)
    }
  }
  gwasloh <- bind_rows(gwasloh, assocchr)
}
gwasloh <- unique(gwasloh)
rownames(gwasloh) <- NULL
fr005 <- sum(table(phenoloh[,phenoname]))*0.05 
fr005
## [1] 166.35
tb <- gwasloh[gwasloh$freq > fr005, ]
rownames(tb) <- NULL
tb
gwasloh$phenotype <- "hospitalization"

saveRDS(gwasloh, file="gwasloh_hospitalization.rds")


###replication
phenoname <- "death"

assocchr <- data.frame()
gwasloh <- data.frame()
chr <- 3
ss <- which(getloh[,2] == chr)
grchr <- GRanges(paste("chr", chr, sep=""), IRanges(start=getloh$Start[ss], end=getloh$End[ss]),
                   sub = getloh[ss]$Sample)
blocks <- c(37530691, 38030691)

for(bb in c(1:(length(blocks)-1))){
    block <- GRanges(paste("chr", chr, sep=""), IRanges(start=blocks[bb],
                                                        end=blocks[bb+1]))
    selover <- data.frame(findOverlaps(block,grchr))[,2]
    selsubsloh <- unique(as.character(grchr$sub[selover]))
    dat <- phenoloh
    dat$loh <- rep(0,nrow(dat))
    if(length(selsubsloh) > 0){
      dat <- dat %>% mutate(loh = ifelse(IID %in% selsubsloh, 1, 0))
      dat$pheno <- dat[,phenoname]
      dat <- dat %>% mutate(pheno = ifelse(pheno==-1, NA, pheno))
      assocblock <- summary(bayesglm(pheno ~ loh+age_at_diagnosis*sex+
                                       PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                     family="binomial", data=dat))$coeff["loh",c(1,4)]
      tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                        OR=exp(assocblock[1]), P=assocblock[2], freq=sum(dat$loh==1))
      assocchr <- bind_rows(assocchr, tmp)
    }
    if(length(selsubsloh) == 0){
      tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                        OR=NA, P=NA, freq=0)
      assocchr <- bind_rows(assocchr, tmp)
    }
  }
gwasloh <- bind_rows(gwasloh, assocchr)
gwasloh$phenotype[1] <- "death"

assocchr <- data.frame()

chr <- 15
ss <- which(getloh[,2] == chr)
grchr <- GRanges(paste("chr", chr, sep=""), IRanges(start=getloh$Start[ss], end=getloh$End[ss]),
                 sub = getloh[ss]$Sample)
blocks <- c(48444756, 48944756)

for(bb in c(1:(length(blocks)-1))){
  block <- GRanges(paste("chr", chr, sep=""), IRanges(start=blocks[bb],
                                                      end=blocks[bb+1]))
  selover <- data.frame(findOverlaps(block,grchr))[,2]
  selsubsloh <- unique(as.character(grchr$sub[selover]))
  dat <- phenoloh
  dat$loh <- rep(0,nrow(dat))
  if(length(selsubsloh) > 0){
    dat <- dat %>% mutate(loh = ifelse(IID %in% selsubsloh, 1, 0))
    dat$pheno <- dat[,phenoname]
    dat <- dat %>% mutate(pheno = ifelse(pheno==-1, NA, pheno))
    assocblock <- summary(bayesglm(pheno ~ loh+age_at_diagnosis*sex+
                                     PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                   family="binomial", data=dat))$coeff["loh",c(1,4)]
    tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                      OR=exp(assocblock[1]), P=assocblock[2], freq=sum(dat$loh==1))
    assocchr <- bind_rows(assocchr, tmp)
  }
  if(length(selsubsloh) == 0){
    tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                      OR=NA, P=NA, freq=0)
    assocchr <- bind_rows(assocchr, tmp)
  }
}
gwasloh <- bind_rows(gwasloh, assocchr)
gwasloh$phenotype[2] <- "death"

assocchr <- data.frame()

chr <- 2
ss <- which(getloh[,2] == chr)
grchr <- GRanges(paste("chr", chr, sep=""), IRanges(start=getloh$Start[ss], end=getloh$End[ss]),
                 sub = getloh[ss]$Sample)
blocks <- c(188011944, 188511944)

for(bb in c(1:(length(blocks)-1))){
  block <- GRanges(paste("chr", chr, sep=""), IRanges(start=blocks[bb],
                                                      end=blocks[bb+1]))
  selover <- data.frame(findOverlaps(block,grchr))[,2]
  selsubsloh <- unique(as.character(grchr$sub[selover]))
  dat <- phenoloh
  dat$loh <- rep(0,nrow(dat))
  if(length(selsubsloh) > 0){
    dat <- dat %>% mutate(loh = ifelse(IID %in% selsubsloh, 1, 0))
    dat$pheno <- dat[,phenoname]
    dat <- dat %>% mutate(pheno = ifelse(pheno==-1, NA, pheno))
    assocblock <- summary(bayesglm(pheno ~ loh+age_at_diagnosis*sex+
                                     PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                   family="binomial", data=dat))$coeff["loh",c(1,4)]
    tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                      OR=exp(assocblock[1]), P=assocblock[2], freq=sum(dat$loh==1))
    assocchr <- bind_rows(assocchr, tmp)
  }
  if(length(selsubsloh) == 0){
    tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                      OR=NA, P=NA, freq=0)
    assocchr <- bind_rows(assocchr, tmp)
  }
}
gwasloh <- bind_rows(gwasloh, assocchr)
gwasloh$phenotype[3] <- "death"


phenoname <- "icu_admit"

assocchr <- data.frame()
chr <- 3
ss <- which(getloh[,2] == chr)
grchr <- GRanges(paste("chr", chr, sep=""), IRanges(start=getloh$Start[ss], end=getloh$End[ss]),
                 sub = getloh[ss]$Sample)
blocks <- c(37530691, 38030691)

for(bb in c(1:(length(blocks)-1))){
  block <- GRanges(paste("chr", chr, sep=""), IRanges(start=blocks[bb],
                                                      end=blocks[bb+1]))
  selover <- data.frame(findOverlaps(block,grchr))[,2]
  selsubsloh <- unique(as.character(grchr$sub[selover]))
  dat <- phenoloh
  dat$loh <- rep(0,nrow(dat))
  if(length(selsubsloh) > 0){
    dat <- dat %>% mutate(loh = ifelse(IID %in% selsubsloh, 1, 0))
    dat$pheno <- dat[,phenoname]
    dat <- dat %>% mutate(pheno = ifelse(pheno==-1, NA, pheno))
    assocblock <- summary(bayesglm(pheno ~ loh+age_at_diagnosis*sex+
                                     PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                   family="binomial", data=dat))$coeff["loh",c(1,4)]
    tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                      OR=exp(assocblock[1]), P=assocblock[2], freq=sum(dat$loh==1))
    assocchr <- bind_rows(assocchr, tmp)
  }
  if(length(selsubsloh) == 0){
    tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                      OR=NA, P=NA, freq=0)
    assocchr <- bind_rows(assocchr, tmp)
  }
}
gwasloh <- bind_rows(gwasloh, assocchr)
gwasloh$phenotype[4] <- "icu admit"

assocchr <- data.frame()

chr <- 15
ss <- which(getloh[,2] == chr)
grchr <- GRanges(paste("chr", chr, sep=""), IRanges(start=getloh$Start[ss], end=getloh$End[ss]),
                 sub = getloh[ss]$Sample)
blocks <- c(48444756, 48944756)

for(bb in c(1:(length(blocks)-1))){
  block <- GRanges(paste("chr", chr, sep=""), IRanges(start=blocks[bb],
                                                      end=blocks[bb+1]))
  selover <- data.frame(findOverlaps(block,grchr))[,2]
  selsubsloh <- unique(as.character(grchr$sub[selover]))
  dat <- phenoloh
  dat$loh <- rep(0,nrow(dat))
  if(length(selsubsloh) > 0){
    dat <- dat %>% mutate(loh = ifelse(IID %in% selsubsloh, 1, 0))
    dat$pheno <- dat[,phenoname]
    dat <- dat %>% mutate(pheno = ifelse(pheno==-1, NA, pheno))
    assocblock <- summary(bayesglm(pheno ~ loh+age_at_diagnosis*sex+
                                     PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                   family="binomial", data=dat))$coeff["loh",c(1,4)]
    tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                      OR=exp(assocblock[1]), P=assocblock[2], freq=sum(dat$loh==1))
    assocchr <- bind_rows(assocchr, tmp)
  }
  if(length(selsubsloh) == 0){
    tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                      OR=NA, P=NA, freq=0)
    assocchr <- bind_rows(assocchr, tmp)
  }
}
gwasloh <- bind_rows(gwasloh, assocchr)
gwasloh$phenotype[5] <- "icu admit"

assocchr <- data.frame()

chr <- 2
ss <- which(getloh[,2] == chr)
grchr <- GRanges(paste("chr", chr, sep=""), IRanges(start=getloh$Start[ss], end=getloh$End[ss]),
                 sub = getloh[ss]$Sample)
blocks <- c(188011944, 188511944)

for(bb in c(1:(length(blocks)-1))){
  block <- GRanges(paste("chr", chr, sep=""), IRanges(start=blocks[bb],
                                                      end=blocks[bb+1]))
  selover <- data.frame(findOverlaps(block,grchr))[,2]
  selsubsloh <- unique(as.character(grchr$sub[selover]))
  dat <- phenoloh
  dat$loh <- rep(0,nrow(dat))
  if(length(selsubsloh) > 0){
    dat <- dat %>% mutate(loh = ifelse(IID %in% selsubsloh, 1, 0))
    dat$pheno <- dat[,phenoname]
    dat <- dat %>% mutate(pheno = ifelse(pheno==-1, NA, pheno))
    assocblock <- summary(bayesglm(pheno ~ loh+age_at_diagnosis*sex+
                                     PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                   family="binomial", data=dat))$coeff["loh",c(1,4)]
    tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                      OR=exp(assocblock[1]), P=assocblock[2], freq=sum(dat$loh==1))
    assocchr <- bind_rows(assocchr, tmp)
  }
  if(length(selsubsloh) == 0){
    tmp <- data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                      OR=NA, P=NA, freq=0)
    assocchr <- bind_rows(assocchr, tmp)
  }
}
gwasloh <- bind_rows(gwasloh, assocchr)
gwasloh$phenotype[6] <- "icu admit"

rownames(gwasloh) <- NULL

saveRDS(gwasloh, file="gwasloh_replication.rds")

