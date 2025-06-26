
# generates SNP summary statistics (100 draws) for a given block

# NB: requires snpStats package for loading PLINK genotype data, see: https://bioconductor.org/packages/release/bioc/html/snpStats.html


# arguments:
# - [1] block ID
# - [2] rG value (scaled from 0 to 100)
# - [3+] named parameter flags (one of each must be set, unless it has a default)
#   - model type: fixed, random -> sets the type of generative model
#   - causal SNP proportion: causal10, causal100 -> set the percentage of SNPs in a block that are causal (10% or 100%) 
#   - shared causal SNPs: unshared, shared (default = shared) -> should causal SNPs be the same for both phenotypes (only if rG = 0)?
#   - h2 scaling: hscale1, hscale10, hscale25 (default = hscale1) -> power scaling, corresponds to signal being concentrated in 100%, 10% or 4% of genomic blocks

# conditions used for the manuscript were all combinations of: 
# - rG of 0, 25, 50, 75, 100, with 10% causal SNPs, and causal SNPs unshared for rG = 0
# - fixed and random model types
# - h2 scaling of 1, 10, and 25


# input files:
# - UK Biobank genotype data PLINK format files: stored per block, subset to the required number of individuals, containing only the SNPs to be used for that block
# - generate_functions.r script file: contains helper functions to generate the output

# output files:
# - .realized file, containing the realized values of the local heritabilities and rG (only relevant for the random effects model)
# - .info file, contains the SNP ID, A1/A2, and sample size columns needed for analysis later (stored separately for efficiency)
# - .stats1 and .stats2 files, contains the generated summary statistics for all 100 draws (as columns), for each phenotype



flags = commandArgs(T)

check.arg = function(input, ..., default=NULL) {
	allowed = c(..., default)
	values = input[input %in% allowed]
	if (length(values) == 0 && !is.null(default)) values = default
	if (length(values) != 1) stop(paste0(ifelse(length(values) > 1, "too many", "no"), " argument matches"))
	return(values)
}


block = as.numeric(flags[1])
rg = as.numeric(flags[2]) 

effect.type = check.arg(flags, "fixed", "random")
prop.causal = ifelse(check.arg(flags, "causal10", "causal100") == "causal100", 1, 0.1)
shared.causal = check.arg(flags, "unshared", default="shared") == "shared"
h2.scale = as.numeric(gsub("hscale", "", check.arg(flags, "hscale10", "hscale25", default="hscale1"))) 

if (rg != 0 && !shared.causal) stop("shared.causal must be TRUE if rG is non-zero")

condition.label = paste0(effect.type, if (h2.scale > 1) paste0("-hscale", h2.scale), "-rg", rg, "-prop", (100*prop.causal))
if (!shared.causal) condition.label = paste0(condition.label, "-unshared")


no.draws = 100
h2.total = c(0.2, 0.8)	# total across all causal blocks
total.snps = 1015592  	# across all chromosomes


root = "[ROOT DIRECTORY]"

blocks.dir = paste0(root, "/blocks")
data.dir = paste0(blocks.dir, "/ukb_main")  # contains UK Biobank genotype data, per blocks  
generate.dir = paste0(root, "/generate")  
output.dir = paste0(generate.dir, "/output/", condition.label, "/block", block)
dir.create(output.dir, recursive=T, showWarnings=F)  


func.file = paste0(generate.dir, "/generate_functions.r")
ukb.prefix = paste0(data.dir, "/ukb_main_block", block)
output.prefix = paste0(output.dir, "/", condition.label, "_block", block)


source(func.file)


# load and standardize genotype data (mean-imputes missing values)
ukb = load.genotypes(ukb.prefix); ukb.info = ukb$info; ukb = ukb$geno

no.indiv = nrow(ukb)
no.snps = nrow(ukb.info)
no.causal = round(prop.causal * no.snps)

h2.block = h2.total * no.snps / total.snps * h2.scale
residual.var = 1 - h2.block



realized = data.frame(h2.1=rep(NA, no.draws), h2.2=NA, rg=NA)
Y = list(matrix(NA, no.indiv, no.draws), matrix(NA, no.indiv, no.draws))
for (i in 1:no.draws) {
	# select the causal SNPs
	causal.snps = select.causal(no.snps, no.causal, shared=shared.causal)
	
	# generate the genetic components G
	if (effect.type == "fixed") components = generate.ffx(ukb, causal.snps, h2.block, rg/100)
	else components = generate.rfx(ukb, causal.snps, h2.block, rg/100)
	
	# compute realized values (will be identical to model parameters under fixed effects model)
	realized[i,1:2] = sapply(components, var)
	realized$rg[i] = cor(components[[1]], components[[2]])

	# add residual to create phenotype (continuous, with population variance of 1)
	for (p in 1:2) Y[[p]][,i] = components[[p]] + rnorm(no.indiv, sd=sqrt(residual.var[p]))
}


write.output = function(dat, suffix) {write.table(dat, file=paste0(output.prefix, ".", suffix), row.names=F, quote=F, sep="\t")}


ukb.info$N = no.indiv
write.output(ukb.info, "info")

write.output(realized, "realized")

sumstats = lapply(1:2, function(p) {compute.stats(ukb, Y[[p]])})	
for (p in 1:2) write.output(sumstats[[p]]$statistic, paste0("stats", p))





