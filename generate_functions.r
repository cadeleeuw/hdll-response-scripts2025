

standardize = function(x) {x = scale(x); x[is.na(x)] = 0; return(scale(x))}

# returns standardized genotypes, with mean-imputed missing values
load.genotypes = function(prefix) {
	dat = snpStats::read.plink(prefix)
	geno = standardize(as(dat$genotypes, "numeric"))
	info = dat$map[,c("snp.name", "allele.1", "allele.2")]; names(info) = c("SNP", "A1", "A2")
	return(list(geno=geno, info=info))
}


# return selected as boolean vectors (list of two, identical if shared=T)
select.causal = function(total, no.causal, shared=T) {
	draw = function(total, no.causal) {1:total %in% sample(total, no.causal)}
	out = list(draw(total, no.causal))
	out[[2]] = if (!shared) draw(total, no.causal) else out[[1]]
	return(out)
}


# generate SNP effects under random effects model, and return corresponding genetic components G
generate.rfx = function(x, causal.snps, h2, rg=0) {
	no.causal = sapply(causal.snps, sum)
	beta = lapply(causal.snps, function(cs) {rep(0, length(cs))})
	
	for (i in 1:2) beta[[i]][causal.snps[[i]]] = rnorm(no.causal[i])  # draw indepently from standard normal
	if (rg != 0) { # create requisite correlation of effects if rG is non-zero (selection of causal SNPs will be the same for both phenotypes in this case, as guaranteed by the main script) 
		if (abs(rg) < 1) beta[[2]][causal.snps[[2]]] = rg * beta[[1]][causal.snps[[1]]] + sqrt(1 - rg^2) * beta[[2]][causal.snps[[2]]]
		else beta[[2]][causal.snps[[2]]] = sign(rg) * beta[[1]][causal.snps[[1]]]  # special case if perfectly correlated 
	}

	# rescale variance to (h2 for the block) / (number of causal SNPs in block)
	components = lapply(1:2, function(i) {x %*% beta[[i]] * sqrt(h2[i] / no.causal[i])})
	return(components)
}

# generate genetic components G under fixed effects model, with h2's and rG exactly equal to specified values
# NB: since x is standardized, means of G are inherently zero and (co)variances reduce to scaled squared/crossproducts 
generate.ffx = function(x, causal.snps, h2, rg=0) {
	# helper function to create genetic component with randomly drawn SNP effects, rescaled to a variance of 1
	get.component = function(causal) {
		beta = rep(0, length(causal)); beta[causal] = rnorm(sum(causal))
		G = x %*% beta
		return(G / sd(G))
	}
	
	N = nrow(x)
	G1 = get.component(causal.snps[[1]])
	
	if (rg == 0) { # causal SNPs may differ for the two phenotypes in this case, so need special procedure to maintain selection of causal SNPs for G2
		# creates two candidate G2 draws, then takes their weighted sum to exactly cancel out any random covariance they individually had with G1 
		G2 = cbind(get.component(causal.snps[[2]]), get.component(causal.snps[[2]])) 
		g.prod = t(G2) %*% G1; wt = -g.prod[2] / g.prod[1]
		G2 = wt * G2[,1] + G2[,2]; G2 = G2 / sd(G2)  # such that cov(G1, G2) = 0
	} else { # in this case, causal SNPs are same for both pheno, so can combine components without changing pattern of causal SNPs
		if (abs(rg) < 1) {
			G2 = lm(get.component(causal.snps[[2]]) ~ G1)$residual; G2 = G2 / sd(G2) # create component fully linearly independent of G1, with variance of 1
			G2 = rg*G1 + sqrt(1 - rg^2) * G2 # induce required correlation
		} else G2 = sign(rg)*G1 # special case if perfectly correlated 
	}
	
	components = list(G1, G2)
	for (p in 1:2) components[[p]] = components[[p]] * sqrt(h2[p])  # rescale to required h2 values
	return(components)
}


	
# computes summary statistics for genotype matrix x and outcome vector/matrix y
# genotype matrix is assumed to be standardized, as is guaranteed by the main script
compute.stats = function(x, y) {
	if (is.null(dim(y))) y = matrix(y, ncol=1)
	stats = stats.core(x, y)
	stats$statistic = stats$beta.hat / sqrt(stats$sampling.var)
	for (var in names(stats)) colnames(stats[[var]]) = paste0("draw", 1:ncol(y))
	return(stats)
}


# core function for computing summary statistics
# x is assumed to be standardized, therefore xtx = N - 1 for every SNP
stats.core = function(x, y) {
	N = nrow(x); no.snps = ncol(x); no.pheno = ncol(y)
	y = sweep(y, 2, apply(y, 2, mean), FUN="-")       					# center y, so intercept can be left implicit in subsequent computation
	var.y = outer(rep(1, no.snps), apply(y, 2, var), FUN="*")   # no.snps x no.pheno matrix with var(Y) values 
	
	beta.hat = crossprod(x, y) / (N - 1)              	# means are all zero, so can omit explicit intercept
	residual.var = var.y - beta.hat^2 * N / (N - 2)   	# adjust for two degrees of freedom, to account for implicit intercept
	sampling.var = residual.var / (N - 1)               
	return(list(beta.hat = beta.hat, sampling.var = sampling.var))
}










