setwd("/home/tomoko/scratch/ROH/")

results_kb <- function(class)  {
  this_iids_roh <- dat[class,]
  my_list<-c("mean"=mean(this_iids_roh$KB[this_iids_roh$KB>=1500]),
             "sum"=sum(this_iids_roh$KB[this_iids_roh$KB>=1500]),
             "length"=length(this_iids_roh$KB[this_iids_roh$KB>=1500]),
             "Froh"=(sum(this_iids_roh$KB[this_iids_roh$KB>=1500]))/2881033)
  return(my_list) }

files_list <- c("EUR.fimm.hom.indiv", "EUR.brasil.hom.indiv", "EUR.canada.hom.indiv", 
                "EUR.kiel.hom.indiv", "EUR.spain_alarcon.hom.indiv", "EUR.spain_planas.hom.indiv"
                )
dat <- data.frame()                          
for (i in 1:length(files_list)) {
  dat <- rbind(dat, read.table((files_list[i]),header=TRUE))
}

library(data.table)
dat <- as.data.table(dat)
dat$IID<-as.factor(dat$IID)
setkey(dat,"IID")
results <- c()
nLevels <- length(levels(dat$IID))
start <- proc.time()
pb <- txtProgressBar(min = 0, max = nLevels, style = 3)
for (i in 1:nLevels){
  this_iid <- levels(dat$IID)[i]
  results <- rbind(results,results_kb(this_iid))
  setTxtProgressBar(pb,i)
}
close(pb)
proc.time()-start
results<-data.frame(levels(dat$IID),results)
results$IID<-results$levels.dat.IID.
results[results==0] <- NA

files_list <- c("EUR.fimm.ibc", "EUR.brasil.ibc", "EUR.canada.ibc", 
                "EUR.kiel.ibc", "EUR.spain_alarcon.ibc", "EUR.spain_planas.ibc")
                
dat <- data.frame()                          
for (i in 1:length(files_list)) {
  dat <- rbind(dat, read.table((files_list[i]),header=TRUE))
}

library(data.table)
dat <- as.data.table(dat)
dat$IID<-as.factor(dat$IID)
dat <- dat %>% select(IID, Fhat1, Fhat3)
 
results <- results %>% merge(dat, by="IID")

write.csv(results, file="ROH_results.csv")

final <- readRDS("~/01.data.QC/curated_clinical/complication_dat_final.rds")
final <- final %>% filter(study != "UKB")
final <- final %>% filter(pop == "EUR")

#FIMM
fimm <- fread("~/covid19-hgi-clinical-values/hgi_fimm/pheno_geno_IDmap")
fimm <- fimm %>% select(anonymized_patient_id, genotypeID)

#Hostage
kiel <- fread("~/covid19-hgi-clinical-values/hgi_kiel/pheno_geno_IDmap")
kiel <- kiel %>% select(anonymized_patient_id, genotypeID)

#Canada
map <- fread("/home/tomoko/11.pca/canada.sampleid.studyid.map", header=F) %>% filter(grepl("JGH",V2)) %>% 
  mutate(anonymized_patient_id = paste0("CA_",as.numeric(str_split(V2, "_",simplify=TRUE)[,2])))
canada <- map %>% rename(genotypeID = V1) %>% select(anonymized_patient_id, genotypeID)
  
#spain_alarcon
map <- fread("/home/tomoko/covid19-hgi-clinical-values/EGA/results/SPGRX.fam") %>% 
  mutate(anonymized_patient_id = paste0("SA_",V2))
spain_alarcon <- map %>% rename(genotypeID = V2) %>% select(anonymized_patient_id, genotypeID)

#spain_planas
map <- fread("/home/tomoko/covid19-hgi-clinical-values/EGA/results/INMUNGEN_CoV2.fam") %>% 
  mutate(anonymized_patient_id = paste0("SP_",str_split(V2, ".CEL",simplify=TRUE)[,1]))
spain_planas <- map %>% rename(genotypeID = V2) %>% select(anonymized_patient_id, genotypeID)

map <- bind_rows(fimm, kiel, canada, spain_planas, spain_alarcon)

results <- results %>% merge(map, by.x="IID", by.y="genotypeID")

data <- results %>% merge(final, by="anonymized_patient_id")
data <- data %>% rename(A2 = a1, B2 = hospitalization)

data %>% saveRDS("pheno_roh.rds")




