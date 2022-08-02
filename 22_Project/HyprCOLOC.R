#install.packages("devtools", repos='http://cran.us.r-project.org')
library(devtools)
library(hyprcoloc, lib = Sys.getenv("R_LIBS_USER"))

args = commandArgs(trailingOnly=TRUE)
print(args[1])
dataRaw <- read.csv(file = args[1])
row.names(dataRaw) <- dataRaw[, c("SNP")]

betasT <- as.matrix(dataRaw[, c("BETA_GWAS", "BETA", "BETA.1", "BETA.2")])
seT <- as.matrix(dataRaw[, c("SE_GWAS", "SE", "SE.1", "SE.2")])

colnames(betasT) <- c("T1", "T2", "T3", "T4")  
colnames(seT) <- c("T1", "T2", "T3", "T4")  


# If test data are not from ind. studies
#trait.cor <- hyprcoloc::test.corr
#ld.matrix <- hyprcoloc::test.ld

traits <- paste0("T", 1:dim(betasT)[2])
rsid <- rownames(betasT)
res <- hyprcoloc(betasT, seT, trait.names=traits, snp.id=rsid)
res

# Need to label traits as either binary or continuous
#binary.traits = c(1,1,1,rep(0,dim(betasT)[2]-3))
binary.traits = c(0,0)
res <- hyprcoloc(betasT, seT, trait.names=traits, snp.id=rsid, binary.outcomes = binary.traits)

# Compute a credible set of snps for each cluster of colocalized traits
res <- hyprcoloc(betasT, seT, trait.names=traits, snp.id=rsid, snpscores = TRUE)

res[[1]]

# first few snp scores for cluster 1
# head(res[[2]][[1]])

# This function identifies a credible set of snps which explains 95% of the posterior
# probability of colocalization for each cluster of colocalized traits
# cred.sets(res, value = 0.95)

# Show that these 2 SNPs are in strong LD
# These 2 SNPs explain nearly 90% of the posterior probability of colocalization and
# therefore likely* that the traits do colocalize, however it is impossible
# to prioritise one of these variants over the other.

#ld.matrix["rs12117612","rs7532349"]


outFileName <- paste(substr(args[1], start = 0, stop=nchar(args[1])-4), "HyprCOLOC.csv", sep="_")

write.csv(res[[1]], outFileName, row.names = FALSE)


