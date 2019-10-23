
# CODE for training
We first extracted the post clumping SNPs  (using default PRSice2 clumping options) from our training data set and rerun PRSice2 with the no-clump flag selected.

```bash
Rscript PRSice.R --cov-file ALL_COVS_ALS.txt --out Results -t Summary_Statistics_GWAS_2016 --beta --snp SNP --A1 A1 --A2 A2 --stat b --se se --pvalue p --ld workingData_BinaryPLINK --print-snp --score std --perm 1000 --prsice PRSice_linux -n 24 --binary-target T --maf 1E-2 --quantile 4 --prevalence 0.00005 --no-full --bar-levels 5E-2 --fastscore --gtf Homo_sapiens.GRCh37.75.gtf --msigdb c5.bp.v6.2.symbols.gmt.txt --multi-plot 10 --no-clump
```

# CODE for Replication

```bash
> R
library("data.table")
VanRheenenTemp <- fread("Summary_Statistics_GWAS_2016", header = T)
temp <- VanRheenenTemp[order(VanRheenenTemp$SNP,VanRheenenTemp$p),]
temp$dupe <- duplicated(temp$SNP)
VanRheenen <- subset(temp, dupe == F)
listOfPathways <- read.table("significant.txt", header = T)
names(listOfPathways) <- c("id")
for(i in 1:length(listOfPathways$id))
{
	instrumentId <- as.character(listOfPathways$id[i])
	thisDataTemp <- fread(file = instrumentId, header = F)
	thisData <- subset(thisDataTemp, V5 =="Y")
	merged <- merge(thisData, VanRheenen, by.x = "V1", by.y="SNP")
	merged$allele <- toupper(merged$A1)
	toPlink <- merged[,c("V1","allele","b")]
	write.table(toPlink, file = paste(instrumentId,".significant.txt",sep = ""), quote = F, sep = " ", col.names = F, row.names = F)
}	
```


```bash
> PLINK
cat significant.txt | while read LINE 
do
echo "THIS IS LINE" $LINE
plink --bfile REPLICATION --score $LINE.bpe.txt --out $LINE.significanttoreplication
done	
```	


```bash
> R
library("data.table")
listOfProfiles <- read.table("significant.txt", header = T)
names(listOfProfiles) <- c("id")
covs1 <- fread("/data/LNG/ALS_PATHWAYS/PHENO_COVARIATES_ALS_LNGonly_unrelated.txt", header = T)
covs2 <- fread("/data/LNG/ALS_PATHWAYS/ALS_Hits_toCOVs.raw", header = T)
covs <- merge (covs1, covs2, by ="FID")
covs$index <- paste(covs$FID,covs$FID, sep = "_")
covs$CASE <- covs$PHENO - 1
outPut <- matrix(ncol = 4, nrow = length(listOfProfiles$id), NA)
colnames(outPut) <- c("pathway","b","se","p")
for(i in 1:length(listOfProfiles$id))
{
	profileName <- 	as.character(listOfProfiles$id[i])
	profile <- fread(file = paste(profileName, ".significanttoreplication.profile", sep = ""), header = T)
	profile$index <- paste(profile$FID, profile$IID, sep = "_")
	data <- merge(covs, profile, by = "index")
	meanControls <- mean(data$SCORE[data$CASE == 0])
	sdControls <- sd(data$SCORE[data$CASE == 0])
	data$zSCORE <- (data$SCORE - meanControls)/sdControls	
	grsTest <- glm(CASE ~ zSCORE + age_at_onset + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, family="binomial", data = data)
	beta <- summary(grsTest)$coefficients["zSCORE","Estimate"]
	se <- summary(grsTest)$coefficients["zSCORE","Std. Error"]
	p <- summary(grsTest)$coefficients["zSCORE","Pr(>|z|)"]
	outPut[i,1] <- profileName
	outPut[i,2] <- beta
	outPut[i,3] <- se
	outPut[i,4] <- p
}
write.table(outPut, "REPLICATION_plink.tab", quote = F, sep = "\t", row.names = F)
```	

