final <- readRDS("~/01.data.QC/curated_clinical/complication_dat_final.rds")
final <- final %>% filter(study != "UKB")
final <- final %>% filter(pop == "EUR")

#FIMM
map <- fread("~/covid19-hgi-clinical-values/hgi_fimm/pheno_geno_IDmap")
fimm <- map %>% filter(anonymized_patient_id %in% final$anonymized_patient_id)
write.table(cbind(0, fimm$genotypeID), file="/home/tomoko/scratch/ROH/EUR.fimm.sample",
            quote=F, col.names = F, row.names = F)

#Hostage
map <- fread("~/covid19-hgi-clinical-values/hgi_kiel/pheno_geno_IDmap")
kiel <- map %>% filter(anonymized_patient_id %in% final$anonymized_patient_id)
write.table(cbind(0, kiel$genotypeID), file="/home/tomoko/scratch/ROH/EUR.kiel.sample",
            quote=F, col.names = F, row.names = F)

#Brasil
map <- fread("~/covid19-hgi-clinical-values/hgi_fimm/pheno_geno_IDmap")
brasil <- map %>% filter(anonymized_patient_id %in% final$anonymized_patient_id) %>% filter(study == "brazil")
write.table(cbind(0, brasil$genotypeID), file="/home/tomoko/scratch/ROH/EUR.brasil.sample",
            quote=F, col.names = F, row.names = F)

#Canada
map <- fread("/home/tomoko/11.pca/canada.sampleid.studyid.map", header=F) %>% filter(grepl("JGH",V2)) %>% 
  mutate(anonymized_patient_id = paste0("CA_",as.numeric(str_split(V2, "_",simplify=TRUE)[,2])))
canada <- map %>% filter(anonymized_patient_id %in% final$anonymized_patient_id)
write.table(cbind(0, canada$V1), file="/home/tomoko/scratch/ROH/EUR.canada.sample",
            quote=F, col.names = F, row.names = F)

#spain_alarcon
map <- fread("/home/tomoko/covid19-hgi-clinical-values/EGA/results/SPGRX.fam") %>% 
  mutate(anonymized_patient_id = paste0("SA_",V2))
spain_alarcon <- map %>% filter(anonymized_patient_id %in% final$anonymized_patient_id)
write.table(cbind(0, spain_alarcon$V2), file="/home/tomoko/scratch/ROH/EUR.spain_alarcon.sample",
            quote=F, col.names = F, row.names = F)

#spain_planas
map <- fread("/home/tomoko/covid19-hgi-clinical-values/EGA/results/INMUNGEN_CoV2.fam") %>% 
  mutate(anonymized_patient_id = paste0("SP_",str_split(V2, ".CEL",simplify=TRUE)[,1]))
spain_planas <- map %>% filter(anonymized_patient_id %in% final$anonymized_patient_id)
write.table(cbind(0, spain_planas$V2), file="/home/tomoko/scratch/ROH/EUR.spain_planas.sample",
            quote=F, col.names = F, row.names = F)




