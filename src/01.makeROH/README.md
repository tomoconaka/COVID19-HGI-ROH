
```
/home/aganna/plink2  \
--bfile ~/covid19-hgi-clinical-values/hgi_${cohort}/geno/${cohort} \
--make-king-table \
--king-table-filter 0.25 \
--out /home/tomoko/scratch/ROH/${cohort}

/home/aganna/plink \
--bfile ~/covid19-hgi-clinical-values/hgi_${cohort}/geno/${cohort} \
--keep /home/tomoko/scratch/ROH/EUR.${cohort}.sample \
--remove <(tail -n+2 /home/tomoko/scratch/ROH/${cohort}.kin0 | awk '{print 0,$4}') \
--maf 0.01 \
--homozyg-window-snp 30 \
--homozyg-window-het 1 \
--homozyg-snp 30 \
--homozyg-kb 300 \
--homozyg-density 30 \
--out /home/tomoko/scratch/ROH/EUR.${cohort}


/home/aganna/plink \
--bfile ~/covid19-hgi-clinical-values/hgi_${cohort}/geno/${cohort} \
--keep /home/tomoko/scratch/ROH/EUR.${cohort}.sample \
--remove <(tail -n+2 /home/tomoko/scratch/ROH/${cohort}.kin0 | awk '{print 0,$4}') \
--maf 0.01 \
--ibc \
--out /home/tomoko/scratch/ROH/EUR.${cohort}
```
