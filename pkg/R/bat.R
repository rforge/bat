#####BAT - Biodiversity Assessment Tools R package
#####required packages
library("vegan")

#####xTree function adapted from http://owenpetchey.staff.shef.ac.uk/Code/Code/calculatingfd_assets/Xtree.r
#####by Jens Schumacher (described in Petchey & Gaston 2002, 2006)
xTree <- function(h) {
	nSpecies <- nrow(as.data.frame(h['order']))
	H1 <- matrix(0, nSpecies, 2 * nSpecies - 2) 
	l <- vector("numeric", 2 * nSpecies - 2)
	for(i in 1:(nSpecies - 1)) {
		if(h$merge[i, 1] < 0) {
			l[2 * i - 1] <- h$height[order(h$height)[i]] 
			H1[ - h$merge[i, 1], 2 * i - 1] <- 1
		} else {
			l[2 * i - 1] <- h$height[order(h$height)[i]] - h$height[order(h$height)[h$merge[i, 1]]]
			H1[, 2 * i - 1] <- H1[, 2 * h$merge[i, 1] - 1] + H1[ , 2 * h$merge[i, 1]]
		} 
		if(h$merge[i, 2] < 0) {
			l[2 * i] <- h$height[order(h$height)[i]] 
			H1[ - h$merge[i, 2], 2 * i] <- 1
		} else {
			l[2 * i] <- h$height[order(h$height)[i]] - h$height[order(h$height)[h$merge[i, 2]]]
			H1[, 2 * i] <- H1[, 2 * h$merge[i, 2] - 1] + H1[, 2 *h$merge[i, 2]]
		}
	} 
	rownames(H1) <- h$labels
	list(l, H1)
}

#####auxiliary functions
prep <- function(comm, tree, abund = TRUE){
	len <- tree[[1]] 							## length of each branch
	A <- tree[[2]]									## matrix species X branches
	minBranch <- min(len[colSums(A)==1]) 	## minimum branch length of terminal branches
	BA <- comm%*%A 												## matrix samples X branches
	if (!abund)	BA = ifelse(BA >= 1, 1, 0)
	return (list(lenBranch = len, sampleBranch = BA, minBranch = minBranch))
}

rarefaction <- function(comm){
	n <- sum(comm)
	for (s in 1:nrow(comm))
		n <- min(n, sum(comm[s,]))
	return(n)
}

#####observed diversity
sobs <- function(comm, tree){
	if (missing(tree)){
		return(length(colSums(comm)[colSums(comm) > 0]))
	} else {
		data <- prep(comm, tree)
		value <- ifelse (colSums(data$sampleBranch) > 0, 1, 0) # vector of observed branches
		return (sum(value*data$lenBranch))
	}
}

#####diversity of rare species for abundance - singletons, doubletons, tripletons, etc
srare <- function(comm, tree, n = 1){
	if(missing(tree)){
		return(length(colSums(comm)[colSums(comm) == n]))
	} else {
		data <- prep(comm, tree)
		value <- ifelse (colSums(data$sampleBranch) == n, 1, 0) # vector of branches with given abundance
		return (sum(value*data$lenBranch))
	}
}

#####diversity of rare species for incidence - uniques, duplicates, triplicates, etc
qrare <- function(comm, tree, n = 1){
	if(missing(tree)){
		comm <- ifelse(comm > 0, 1, 0)
		return(length(colSums(comm)[colSums(comm) == n]))
	} else {
		data <- prep(comm, tree, FALSE)
		value <- ifelse (colSums(data$sampleBranch) == n, 1, 0) # vector of branches with given incidence
		return (sum(value*data$lenBranch))
	}
}

#####minimum terminal branch length, =1 in case of TD
minBranch <- function(comm, tree){
	if (missing(tree)){
		return(1)
	} else {
		data <- prep(comm, tree)
		return(data$minBranch)
	}
}

#####non-parametric estimators
chao <- function(obs, s1, s2, mb){
	return(obs + (s1*(s1-mb))/(2*(s2+mb)))
}

jack1ab <- function(obs, s1){
	return(obs + s1)
}

jack1in <- function(obs, q1, q){
	return(obs + q1 * ((q-1)/q))
}

jack2ab <- function(obs, s1, s2){
	return(obs + 2*s1 - s2)
}

jack2in <- function(obs, q1, q2, q){
	if (q > 1)	return(obs + (q1*(2*q-3)/q - q2*(q-2)^2/(q*(q-1))))
	else return(obs + 2*q1 - q2)
}

pcorr <- function(obs, s1){
	return(1+(s1/obs)^2)
}

#####observed beta (a = shared species/edges, b/c = species/edges exclusive to either site, comm is a 2sites x species matrix)
betaObs <- function(comm, tree, abund = FALSE, func = "jaccard"){
	if (!abund) {														##if incidence data
		obs1 <- sobs(comm[1,,drop=F], tree)
		obs2 <- sobs(comm[2,,drop=F], tree)
		obsBoth <- sobs(comm, tree)
		if(tolower(substr(func, 1, 1)) != "s")
			denominator <- obsBoth
		else
			denominator <- obs1 + obs2
		b <- obsBoth - obs2
		c <- obsBoth - obs1
	} else {																				##if abundance data
		a <- 0
		for (i in 1:ncol(comm))
			a <- a + min(comm[1,i], comm[2,i])
		b <- sum(comm[1,]) - a
		c <- sum(comm[2,]) - a
		denominator <- a + b + c
		if(tolower(substr(func, 1, 1)) == "s")
			denominator <- denominator + a
	}
	return(list(Btotal = (b+c)/denominator, Brepl = 2*min(b,c)/denominator, Brich = abs(b-c)/denominator))
}

#' Alpha diversity (TD, PD or FD).
#' @description Observed alpha diversity with possible rarefaction, multiple sites simultaneously.
#' @param comm A sites x species matrix, with either abundance or incidence data.
#' @param tree An hclust or phylo object (used only for PD or FD).
#' @param raref An integer specifying the number of individuals for rarefaction (individual based).
#' If raref < 1 no rarefaction is made.
#' If raref = 1 rarefaction is made by the minimum abundance among all sites.
#' If raref > 1 rarefaction is made by the abundance indicated.
#' @param runs Number of resampling runs for rarefaction.
#' @details TD is equivalent to species richness. Calculations of PD and FD are based on Faith (1992) and Petchey & Gaston (2002, 2006), which measure PD and FD of a community as the total branch length of a tree linking all species represented in such community.
#' PD and FD are calculated based on an ultrametric tree (hclust or phylo object).
#' The number and order of species in comm must be the same as in tree.
#' @return A matrix of sites x diversity values (either "Obs" OR "Avg, Min, LowerCL, UpperCL and Max").
#' @references Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @references Petchey, O.L. & Gaston, K.J. (2002) Functional diversity (FD), species richness and community composition. Ecology Letters, 5, 402-411.
#' @references Petchey, O.L. & Gaston, K.J. (2006) Functional diversity: back to basics and looking forward. Ecology Letters, 9, 741-758.
#' @examples comm <- matrix(c(0,0,1,1,0,0,2,1,0,0), nrow = 2, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' alpha(comm)
#' alpha(comm, raref = 0)
#' alpha(comm, tree, 2, 100)
alpha <- function(comm, tree, raref = 0, runs = 1000){
	if (!missing(tree))
		tree <- xTree(as.hclust(tree))
	nComm <- nrow(comm)
	if(raref < 1){						# no rarefaction if 0 or negative
		results <- matrix(0, nComm, 1)
		for (s in 1:nComm){
			results[s,1] <- sobs(comm[s,, drop=F], tree)
		}
		colnames(results) <- "Obs"
		return (results)
	}
	if (raref == 1)
		raref <- rarefaction(comm)				# rarefy by minimum n among all communities
	results <- matrix(0, nComm, 5)
	for (s in 1:nComm){
		res <- c()
		for (r in 1:runs){
			res <- c(res,sobs(rrarefy(comm[s,], raref), tree))
		}
		results[s,] <- c(mean(res), min(res), quantile(res, 0.025), quantile(res, 0.975), max(res))
	}
	colnames(results) <- c("Avg", "Min", "LowerCL", "UpperCL", "Max")
	return (results)
}

#' Alpha diversity accumulation curves (observed and estimated).
#' @description Estimation of alpha diversity of a single site with accumulation of samples.
#' @param comm A samples x species matrix, with either abundance or incidence data.
#' @param tree An hclust or phylo object (used only for PD or FD).
#' @param func The class of estimators to be used:
#' If func ~ "curve", TD, PD or FD are based on extrapolating the accumulation curve of observed diversity.
#' If func ~ "nonparametric", TD, PD or FD are based on non-parametric estimators.
#' If func ~ "completeness", PD or FD estimates are based on the completeness of TD (requires a tree to be used).
#' @param runs Number of random permutations to be made to the sample order.
#' @details Observed diversity often is an underestimation of true diversity. Several approaches have been devised to estimate species richness (TD) from incomplete sampling.
#' These include: (1) fitting asymptotic functions to randomised accumulation curves (Soberon & Llorente 1993; Flather 1996)
#' (2) the use of non-parametric estimators based on the incidence or abundance of rare species (Heltshe & Forrester 1983; Chao 1984, 1987; Colwell & Coddington 1994).
#' A correction to non-parametric estimators has also been recently proposed, based on the proportion of singleton or unique species
#' (species represented by a single individual or in a single sample respectively; Lopez et al. 2012).
#' Cardoso et al. (2014) have proposed a way of adapting these approaches to estimate PD and FD, also adding a third possible approach for
#' these dimensions of diversity: (3) correct PD and FD values based on the completeness of TD, where completeness equals the proportion of estimated true diversity that was observed.
#' Calculations of PD and FD are based on Faith (1992) and Petchey & Gaston (2002, 2006), which measure PD and FD of a community as the total branch length of a tree linking all species represented in such community.
#' PD and FD are calculated based on an ultrametric tree (hclust or phylo object).
#' The number and order of species in comm must be the same as in tree.
#' @return A matrix of samples x diversity values (samples, individuals, observed and estimated diversity).
#' @references Cardoso, P., Rigal, F., Borges, P.A.V. & Carvalho, J.C. (2014) A new frontier in biodiversity inventory: a proposal for estimators of phylogenetic and functional diversity. Methods in Ecology and Evolution, in press.
#' @references Chao, A. (1984) Nonparametric estimation of the number of classes in a population. Scandinavian Journal of Statistics, 11, 265-270.
#' @references Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. Biometrics 43, 783-791.
#' @references Colwell, R.K. & Coddington, J.A. (1994). Estimating terrestrial biodiversity through extrapolation. Phil. Trans. Roy. Soc. London B 345, 101–118.
#' @references Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @references Flather, C. (1996) Fitting species-accumulation functions and assessing regional land use impacts on avian diversity. Journal of Biogeography, 23, 155-168.
#' @references Heltshe, J. & Forrester, N.E. (1983) Estimating species richness using the jackknife procedure. Biometrics, 39, 1-11.
#' @references Lopez, L.C.S., Fracasso, M.P.A., Mesquita, D.O., Palma, A.R.T. & Riul, P. (2012) The relationship between percentage of singletons and sampling effort: a new approach to reduce the bias of richness estimates. Ecological Indicators, 14, 164-169.
#' @references Petchey, O.L. & Gaston, K.J. (2002) Functional diversity (FD), species richness and community composition. Ecology Letters, 5, 402-411.
#' @references Petchey, O.L. & Gaston, K.J. (2006) Functional diversity: back to basics and looking forward. Ecology Letters, 9, 741-758.
#' @references Soberon, M.J. & Llorente, J. (1993) The use of species accumulation functions for the prediction of species richness. Conservation Biology, 7, 480-488.
#' @examples comm <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' alpha.accum(comm)
#' alpha.accum(comm, func = "nonparametric")
#' alpha.accum(comm, tree, "compl")
#' alpha.accum(comm, tree, "curve", 100)
alpha.accum <- function(comm, tree, func = "nonparametric", runs = 1000){
	if (!missing(tree))
		tree <- xTree(as.hclust(tree))
	
	#####function options:
	#####nonparametric (TD/PD/FD with non-parametric estimators)
	#####completeness (PD/FD with TD completeness correction)
	#####curve (TD/PD/FD with curve fitting)
	func <- match.arg(func, c("nonparametric", "completeness", "curve"))
	
	#####nonparametric (TD/PD/FD with non-parametric estimators)
	switch(func, nonparametric = {
		resultsArray <- array(0, dim = c(nrow(comm), 19, runs))
		for (r in 1:runs){
			comm <- comm[sample(nrow(comm)),, drop=FALSE]			#shuffle rows (samples)
			data <- matrix(0,1,ncol(comm))
			runData <- matrix(0,0,19)
			for (q in 1:nrow(comm)){
				data <- rbind(data, comm[q,])
				n <- sum(rowSums(data))
				obs <- sobs(data, tree)
				s1 <- srare(data, tree, 1)
				s2 <- srare(data, tree, 2)
				q1 <- qrare(data, tree, 1)
				q2 <- qrare(data, tree, 2)
				mb <- minBranch(data, tree)
				j1ab <- jack1ab(obs, s1)
				j1abP <- j1ab * pcorr(obs, s1)
				j1in <- jack1in(obs, q1, q)
				j1inP <- j1in * pcorr(obs, q1)
				j2ab <- jack2ab(obs, s1, s2)
				j2abP <- j2ab * pcorr(obs, s1)
				j2in <- jack2in(obs, q1, q2, q)
				j2inP <- j2in * pcorr(obs, q1)
				c1 <- chao(obs, s1, s2, mb)
				c1P <- c1 * pcorr(obs, s1)
				c2 <- chao(obs, q1, q2, mb)
				c2P <- c2 * pcorr(obs, q1)
				runData <- rbind(runData, c(q, n, obs, s1, s2, q1, q2, j1ab, j1abP, j1in, j1inP, j2ab, j2abP, j2in, j2inP, c1, c1P, c2, c2P))
			}
			resultsArray[,,r] <- runData
		}
		
		#####calculate averages or medians of all runs
		results <- matrix(0,nrow(comm),19)
		v <- array(0, dim = c(runs))
		for (i in 1:nrow(comm)){
			for (j in 1:19){
				for (k in 1:runs){
					v[k] <- resultsArray[i,j,k]
				}
				if (j < 16 || missing(tree))
					results[i,j] <- mean(v)
				else
					results[i,j] <- median(v)
			}
		}
		
		#####completeness (PD/FD with TD completeness correction)
	}, completeness = {
		if (missing(tree))
			stop("Completeness option not available without a tree...")
		results <- alpha.accum(comm, , "nonparametric", runs)
		obs <- matrix(0,nrow(comm),1)
		for (r in 1:runs){
			comm <- comm[sample(nrow(comm)),, drop=FALSE]			#shuffle rows (samples)
			for (s in 1:nrow(comm)){
				obs[s,1] <- obs[s,1] + sobs(comm[1:s,], tree)
			}
		}
		obs <- obs / runs
		for (i in 8:19)
			results[,i] <- obs * (results[,i] / results[,3])
		results[,3] <- obs
		
		#####curve (TD/PD/FD with curve fitting)
	}, curve = {
		results <- matrix(0,nrow(comm),6)
		for (r in 1:runs){
			comm <- comm[sample(nrow(comm)),, drop=FALSE]		#shuffle rows (samples)
			runData <- matrix(0,0,6)
			for (s in 1:nrow(comm)){
				n <- sum(rowSums(comm[1:s,,drop=FALSE]))
				obs <- sobs(comm[1:s,,drop=FALSE], tree)
				runData <- rbind(runData, c(s,n,obs,obs,obs,obs))
			}
			results <- results + runData
		}
		results <- results / runs
		for (s in 3:nrow(results)){				##fit curves only with 3 or more samples
			x <- results[1:s,1]
			y <- results[1:s,3]
			#####Clench
			clench <- try(nls(y ~ (a*x)/(1+b*x), start = list(a = 100, b = 1), control = nls.control(maxiter = 10000)), silent = TRUE); # does not stop in the case of error
			if(class(clench) != "try-error"){
				a <- coef(clench)[1]
				b <- coef(clench)[2]
				results[s,4] <- a/b
			}
			#####Exponential
			exponential <- try(nls(y ~ (a/b)*(1-exp(-b*x)), start = list(a = 100, b = 1), control = nls.control(maxiter = 1000)), silent = TRUE); # does not stop in the case of error
			if(class(exponential) != "try-error"){
				a <- coef(exponential)[1]
				b <- coef(exponential)[2]
				results[s,5] <- a/b
			}
			#####Weibull
			weibull <- try(nls(y ~ a*(1-exp(-(b*(x-c))^d)), start = list(a = 1, b = 1, c = 1, d = 1), control = nls.control(maxiter = 10000)), silent = TRUE); # does not stop in the case of error
			if(class(weibull) != "try-error"){
				a <- coef(weibull)[1]
				results[s,6] <- a
			}
		}
		colnames(results) <- c("Samples", "Ind", "Obs", "Clench", "Exponential", "Weibull")
		return (results)
	})
	colnames(results) <- c("Samples", "Ind", "Obs", "S1", "S2", "Q1", "Q2", "Jack1ab", "Jack1abP", "Jack1in", "Jack1inP", "Jack2ab", "Jack2abP", "Jack2in", "Jack2inP", "Chao1", "Chao1P", "Chao2", "Chao2P")
	return(results)
}

#' Alpha diversity estimates.
#' @description Estimation of alpha diversity of multiple sites simultaneously.
#' @param comm A sites x species matrix, with either abundances or number of incidences.
#' @param tree An hclust or phylo object (used only for PD or FD).
#' @param func The class of estimators to be used:
#' If func ~ "nonparametric", TD, PD or FD are based on non-parametric estimators.
#' If func ~ "completeness", PD or FD estimates are based on the completeness of TD (requires a tree to be used).
#' @details Observed diversity often is an underestimation of true diversity.
#' Non-parametric estimators based on the incidence or abundance of rare species have been proposed to overcome the problem of undersampling (Heltshe & Forrester 1983; Chao 1984, 1987; Colwell & Coddington 1994).
#' A correction to non-parametric estimators has also been recently proposed, based on the proportion of singleton or unique species
#' (species represented by a single individual or in a single sample respectively; Lopez et al. 2012).
#' Cardoso et al. (2014) have proposed a way of adapting non-parametric species richness estimators to PD and FD. They have also proposed correcting PD and FD values based on the completeness of TD, where completeness equals the proportion of estimated true diversity that was observed.
#' Calculations of PD and FD are based on Faith (1992) and Petchey & Gaston (2002, 2006), which measure PD and FD of a community as the total branch length of a tree linking all species represented in such community.
#' PD and FD are calculated based on an ultrametric tree (hclust or phylo object).
#' The number and order of species in comm must be the same as in tree.
#' @return A matrix of sites x diversity values (individuals, observed and estimated diversity).
#' @references Cardoso, P., Rigal, F., Borges, P.A.V. & Carvalho, J.C. (2014) A new frontier in biodiversity inventory: a proposal for estimators of phylogenetic and functional diversity. Methods in Ecology and Evolution, in press.
#' @references Chao, A. (1984) Nonparametric estimation of the number of classes in a population. Scandinavian Journal of Statistics, 11, 265-270.
#' @references Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. Biometrics 43, 783-791.
#' @references Colwell, R.K. & Coddington, J.A. (1994). Estimating terrestrial biodiversity through extrapolation. Phil. Trans. Roy. Soc. London B 345, 101–118.
#' @references Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @references Heltshe, J. & Forrester, N.E. (1983) Estimating species richness using the jackknife procedure. Biometrics, 39, 1-11.
#' @references Lopez, L.C.S., Fracasso, M.P.A., Mesquita, D.O., Palma, A.R.T. & Riul, P. (2012) The relationship between percentage of singletons and sampling effort: a new approach to reduce the bias of richness estimates. Ecological Indicators, 14, 164-169.
#' @references Petchey, O.L. & Gaston, K.J. (2002) Functional diversity (FD), species richness and community composition. Ecology Letters, 5, 402-411.
#' @references Petchey, O.L. & Gaston, K.J. (2006) Functional diversity: back to basics and looking forward. Ecology Letters, 9, 741-758.
#' @examples comm <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' alpha.estimate(comm)
#' alpha.estimate(comm, tree)
#' alpha.estimate(comm, tree, func = "completeness")
alpha.estimate <- function(comm, tree, func = "nonparametric"){
	if (max(comm) == 1)
		stop("No estimates are possible without abundance or incidence frequency data")
	if (!missing(tree))
		tree <- xTree(as.hclust(tree))
	
	#####function options:
	#####nonparametric (TD/PD/FD with non-parametric estimators)
	#####completeness (PD/FD with TD completeness correction)
	func <- match.arg(func, c("nonparametric", "completeness"))
	
	#####nonparametric (TD/PD/FD with non-parametric estimators)
	switch(func, nonparametric = {
		results <- matrix(0,0,10)
		for (s in 1:nrow(comm)){
			data <- comm[s,,drop = F]
			obs <- sobs(data, tree)
			n <- sum(data)
			s1 <- srare(data, tree, 1)
			s2 <- srare(data, tree, 2)
			mb <- minBranch(data, tree)
			j1ab <- jack1ab(obs, s1)
			j1abP <- j1ab * pcorr(obs, s1)
			j2ab <- jack2ab(obs, s1, s2)
			j2abP <- j2ab * pcorr(obs, s1)
			c1 <- chao(obs, s1, s2, mb)
			c1P <- c1 * pcorr(obs, s1)
			results <- rbind(results, c(n, obs, s1, s2, j1ab, j1abP, j2ab, j2abP, c1, c1P))
		}
		
		#####completeness (PD/FD with TD completeness correction)
	}, completeness = {
		if (missing(tree))
			stop("Completeness option not available without a tree...")
		results <- alpha.estimate(comm, , "nonparametric")
		obs <- matrix(0,nrow(comm),1)
		for (s in 1:nrow(comm))
			obs[s,1] <- obs[s,1] + sobs(comm[s,], tree)
		for (i in 5:10)
			results[,i] <- obs[,1] * (results[,i] / results[,2])
		results[,2] <- obs[,1]
	})
	colnames(results) <- c("Ind", "Obs", "S1", "S2", "Jack1ab", "Jack1abP", "Jack2ab", "Jack2abP", "Chao1", "Chao1P")
	return(results)
}

#' Beta diversity (TD, PD or FD).
#' @description Beta diversity with possible rarefaction, multiple sites simultaneously.
#' @param comm A sites x species matrix, with either abundance or incidence data.
#' @param tree An hclust or phylo object (used only for PD or FD).
#' @param abund A boolean (T/F) indicating whether abundance data should be used or converted to incidence before analysis.
#' @param func Indicates whether the Jaccard or Soerensen family of beta diversity measures should be used.
#' @param raref An integer specifying the number of individuals for rarefaction (individual based).
#' If raref < 1 no rarefaction is made.
#' If raref = 1 rarefaction is made by the minimum abundance among all sites.
#' If raref > 1 rarefaction is made by the abundance indicated.
#' @param runs Number of resampling runs for rarefaction.
#' @details The beta diversity measures used here follow the partitioning framework independently developed by Podani & Schmera (2011) and Carvalho et al. (2012)
#' and later expanded to PD and FD by Cardoso et al. (2014), where Btotal = Brepl + Brich. Btotal = total beta diversity, reflecting both species replacement and loss/gain;
#' Brepl = beta diversity explained by replacement of species alone; Brich = beta diversity explained by species loss/gain (richness differences) alone.
#' PD and FD are calculated based on an ultrametric tree (hclust or phylo object).
#' The number and order of species in comm must be the same as in tree.
#' @return Three distance matrices between sites, one per each of the three beta diversity measures (either "Obs" OR "Avg, Min, LowerCL, UpperCL and Max").
#' @references Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. Journal of Biogeography, 41, 749-761.
#' @references Carvalho, J.C., Cardoso, P. & Gomes, P. (2012) Determining the relative roles of species replacement and species richness differences in generating beta-diversity patterns. Global Ecology and Biogeography, 21, 760-771.
#' @references Podani, J. & Schmera, D. (2011) A new conceptual and methodological framework for exploring and explaining pattern in presence–absence data. Oikos, 120, 1625-1638.
#' @examples comm <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' beta(comm)
#' beta(comm, func = "Soerensen")
#' beta(comm, tree)
#' beta(comm, raref = 1)
#' beta(comm, tree, T, "s", raref = 2)
beta <- function(comm, tree, abund = FALSE, func = "jaccard", raref = 0, runs = 1000){
	if (!missing(tree))
		tree <- xTree(as.hclust(tree))
	
	nComm <- nrow(comm)
	
	if(raref < 1){						# no rarefaction if 0 or negative
		results <- array(0, dim=c(nComm, nComm, 3))
		for (i in 1:nComm){
			for (j in 1:nComm){
				commBoth <- as.matrix(rbind(comm[i,], comm[j,]))
				betaValues <- betaObs(commBoth, tree, abund, func)
				for(k in 1:3)
					results[i,j,k] <- betaValues[[k]]
			}
		}
		results <- list(Btotal = as.dist(results[,,1]),Brepl = as.dist(results[,,2]),Brich = as.dist(results[,,3]))
		return (results)
	}
	if (raref == 1)
		raref <- rarefaction(comm)				# rarefy by minimum n among all communities
	results <- array(0, dim=c(nComm, nComm, 3, 5))
	
	for (i in 1:nComm){
		for (j in 1:nComm){
			run <- matrix(0, runs, 3)
			for (r in 1:runs){
				commBoth <- as.matrix(rbind(rrarefy(comm[i,], raref), rrarefy(comm[j,], raref)))
				betaValues <- betaObs(commBoth, tree, abund, func)
				run[r,1] <- betaValues$Btotal
				run[r,2] <- betaValues$Brepl
				run[r,3] <- betaValues$Brich
			}
			for (b in 1:3){
				results[i,j,b,1] <- mean(run[,b])
				results[i,j,b,2] <- min(run[,b])
				results[i,j,b,3] <- quantile(run[,b], 0.025)
				results[i,j,b,4] <- quantile(run[,b], 0.975)
				results[i,j,b,5] <- max(run[,b])
			}
		}
	}
	results.total <- list(Btotal = as.dist(results[,,1,1]), Btotal.min = as.dist(results[,,1,2]), Btotal.lowCL = as.dist(results[,,1,3]), Btotal.upCL = as.dist(results[,,1,4]), Btotal.max = as.dist(results[,,1,5]))
	results.repl <- list(Brepl = as.dist(results[,,2,1]), Brepl.min = as.dist(results[,,2,2]), Brepl.lowCL = as.dist(results[,,2,3]), Brepl.upCL = as.dist(results[,,2,4]), Brepl.max = as.dist(results[,,2,5]))
	results.rich <- list(Brich = as.dist(results[,,3,1]), Brich.min = as.dist(results[,,3,2]), Brich.lowCL = as.dist(results[,,3,3]), Brich.upCL = as.dist(results[,,3,4]), Brich.max = as.dist(results[,,3,5]))
	results <- c(results.total, results.repl, results.rich)
	return (results)
}

#' Beta diversity accumulation curves.
#' @description Beta diversity between two sites with accumulation of samples.
#' @param comm1 A samples x species matrix for the first site, with either abundance or incidence data.
#' @param comm2 A samples x species matrix for the second site, with either abundance or incidence data.
#' @param tree An hclust or phylo object (used only for PD or FD).
#' @param abund A boolean (T/F) indicating whether abundance data should be used or converted to incidence before analysis.
#' @param func Indicates whether the Jaccard or Soerensen family of beta diversity measures should be used.
#' @param runs Number of random permutations to be made to the sample order.
#' @details As widely recognized for species richness, beta diversity is also biased when communities are undersampled.
#' Beta diversity accumulation curves have been proposed by Cardoso et al. (2009) to test if beta diversity has approached an asymptote when comparing two undersampled sites.
#' The beta diversity measures used here follow the partitioning framework independently developed by Podani & Schmera (2011) and Carvalho et al. (2012)
#' and later expanded to PD and FD by Cardoso et al. (2014), where Btotal = Brepl + Brich. Btotal = total beta diversity, reflecting both species replacement and loss/gain;
#' Brepl = beta diversity explained by replacement of species alone; Brich = beta diversity explained by species loss/gain (richness differences) alone.
#' PD and FD are calculated based on an ultrametric tree (hclust or phylo object).
#' The number and order of species in comm1 and comm2 must be the same as in tree. Also, the number of samples should be similar in both sites.
#' @return Three matrices of samples x diversity values, one per each of the three beta diversity measures (samples, individuals and observed diversity).
#' @references Cardoso, P., Borges, P.A.V. & Veech, J.A. (2009) Testing the performance of beta diversity measures based on incidence data: the robustness to undersampling. Diversity and Distributions, 15, 1081-1090.
#' @references Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. Journal of Biogeography, 41, 749-761.
#' @references Carvalho, J.C., Cardoso, P. & Gomes, P. (2012) Determining the relative roles of species replacement and species richness differences in generating beta-diversity patterns. Global Ecology and Biogeography, 21, 760-771.
#' @references Podani, J. & Schmera, D. (2011) A new conceptual and methodological framework for exploring and explaining pattern in presence–absence data. Oikos, 120, 1625-1638.
#' @examples comm1 <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
#' comm2 <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' beta.accum(comm1, comm2)
#' beta.accum(comm1, comm2, func = "Soerensen")
#' beta.accum(comm2, comm3, tree)
#' beta.accum(comm2, comm3, abund = TRUE)
#' beta.accum(comm2, comm3, tree, F)
beta.accum <- function(comm1, comm2, tree, abund = FALSE, func = "jaccard", runs = 1000){
	if(nrow(comm1) < 2 || nrow(comm1) != nrow(comm2))
		stop("Both communities should have multiple and the same number of samples")
	if (!missing(tree))
		tree <- xTree(as.hclust(tree))
	
	nSamples <- nrow(comm1)
	results <- matrix(0,nSamples, 4)
	colnames(results) <- c("Samples", "Btotal", "Brepl", "Brich")
	
	for (r in 1:runs){
		comm1 <- comm1[sample(nSamples),, drop=FALSE]			#shuffle samples of first community
		comm2 <- comm2[sample(nSamples),, drop=FALSE]			#shuffle samples of second community
		for (q in 1:nSamples){
			commBoth <- as.matrix(rbind(colSums(comm1[1:q,,drop=F]),colSums(comm2[1:q,,drop=F])))
			results[q,1] <- results[q,1] + q
			betaValues <- betaObs(commBoth, tree, abund, func)
			results[q,2] <- results[q,2] + betaValues$Btotal
			results[q,3] <- results[q,3] + betaValues$Brepl
			results[q,4] <- results[q,4] + betaValues$Brich
		}
	}
	results <- results/runs
	return (results)
}

#' Beta diversity among multiple sites.
#' @description Beta diversity with possible rarefaction - multiple sites measure calculated as the average or variance of all pairwise values.
#' @param comm A sites x species matrix, with either abundance or incidence data.
#' @param tree An hclust or phylo object (used only for PD or FD).
#' @param abund A boolean (T/F) indicating whether abundance data should be used or converted to incidence before analysis.
#' @param func Indicates whether the Jaccard or Soerensen family of beta diversity measures should be used.
#' @param raref An integer specifying the number of individuals for rarefaction (individual based).
#' If raref < 1 no rarefaction is made.
#' If raref = 1 rarefaction is made by the minimum abundance among all sites.
#' If raref > 1 rarefaction is made by the abundance indicated.
#' @param runs Number of resampling runs for rarefaction.
#' @details Beta diversity of multiple sites simultaneously is calculated as either the average or the variance among all pairwise comparisons (Legendre, 2014).
#' The beta diversity measures used here follow the partitioning framework independently developed by Podani & Schmera (2011) and Carvalho et al. (2012)
#' and later expanded to PD and FD by Cardoso et al. (2014), where Btotal = Brepl + Brich. Btotal = total beta diversity, reflecting both species replacement and loss/gain;
#' Brepl = beta diversity explained by replacement of species alone; Brich = beta diversity explained by species loss/gain (richness differences) alone.
#' PD and FD are calculated based on an ultrametric tree (hclust or phylo object).
#' The number and order of species in comm must be the same as in tree.
#' @return A matrix of beta measures x diversity values (average and variance).
#' @references Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. Journal of Biogeography, 41, 749-761.
#' @references Carvalho, J.C., Cardoso, P. & Gomes, P. (2012) Determining the relative roles of species replacement and species richness differences in generating beta-diversity patterns. Global Ecology and Biogeography, 21, 760-771.
#' @references Legendre, P. (2014) Interpreting the replacement and richness difference components of beta diversity. Global Ecology and Biogeography, in press.
#' @references Podani, J. & Schmera, D. (2011) A new conceptual and methodological framework for exploring and explaining pattern in presence–absence data. Oikos, 120, 1625-1638.
#' @examples comm <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' beta.multi(comm)
#' beta.multi(comm, func = "Soerensen")
#' beta.multi(comm, tree)
#' beta.multi(comm, raref = 1)
#' beta.multi(comm, tree, T, "s", raref = 2)
beta.multi <- function(comm, tree, abund = FALSE, func = "jaccard", raref = 0, runs = 1000){
	pairwise <- beta(comm, tree, abund, func, raref, runs)
	Btotal.avg <- mean(pairwise$Btotal)
	Brepl.avg <- mean(pairwise$Brepl)
	Brich.avg <- mean(pairwise$Brich)
	Btotal.var <- sum(pairwise$Btotal)/(ncol(comm)*(ncol(comm)-1))
	Brepl.var <- sum(pairwise$Brepl)/(ncol(comm)*(ncol(comm)-1))
	Brich.var <- sum(pairwise$Brich)/(ncol(comm)*(ncol(comm)-1))
	results <- matrix(c(Btotal.avg, Brepl.avg, Brich.avg, Btotal.var, Brepl.var, Brich.var), nrow = 3, ncol = 2)
	colnames(results) <- c("Average", "Variance")
	rownames(results) <- c("Btotal", "Brepl", "Brich")
	return(results)
}

#' Scaled mean square error of accumulation curves.
#' @description Accuracy (scaled mean square error) of accumulation curves compared with a known true diversity value (target).
#' @param accum A matrix resulting from the alpha.accum or beta.accum functions (samples x diversity values).
#' @param target The true known diversity value, with which the curve will be compared. If = -1, it is the diversity observed with all samples.
#' @details Among multiple measures of accuracy (Walther & Moore 2005) the SMSE presents several advantages, as it is (Cardoso et al. 2014):
#' (i) scaled to true diversity, so that similar absolute differences are weighted according to how much they represent of the real value;
#' (ii) scaled to the number of samples, so that values are independent of sample size;
#' (iii) squared, so that small, mostly meaningless fluctuations around the true value are down-weighted; and
#' (iv) independent of positive or negative deviation from the real value, as such differentiation is usually not necessary.
#' @return Accuracy values for all observed and estimated curves.
#' @references Cardoso, P., Rigal, F., Borges, P.A.V. & Carvalho, J.C. (2014) A new frontier in biodiversity inventory: a proposal for estimators of phylogenetic and functional diversity. Methods in Ecology and Evolution, in press.
#' @references Walther, B.A. & Moore, J.L. (2005) The concepts of bias, precision and accuracy, and their use in testing the performance of species richness estimators, with a literature reviewof estimator performance. Ecography, 28, 815-829.
#' @examples comm1 <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
#' comm2 <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' acc.alpha = alpha.accum(comm1)
#' accuracy(acc.alpha)
#' accuracy(acc.alpha, 10)
#' acc.beta = beta.accum(comm1, comm2, tree)
#' accuracy(acc.beta)
#' accuracy(acc.beta, c(1,1,0))
accuracy <- function(accum, target = -1){
	if(ncol(accum) > 10){																#if alpha
		if (target == -1)
			target <- accum[nrow(accum), 3]
		smse <- rep(0, 13)
		for (i in 1:nrow(accum)){
			smse[1] <- smse[1] + (accum[i,3] - target)^2
			for (j in 2:13)
				smse[j] <- smse[j] + (accum[i,j+6] - target)^2
		}
		smse <- smse / (target^2*nrow(accum))
		smse <- list(Obs=smse[1], Jack1ab=smse[2], Jack1abP=smse[3], Jack1in=smse[4], Jack1inP=smse[5], Jack2ab=smse[6], Jack2abP=smse[7], Jack2in=smse[8], Jack2inP=smse[9], Chao1=smse[10], Chao1P=smse[11], Chao2=smse[12], Chao2P=smse[13])
	} else {																						#if beta
		if (target[1] == -1)
			target <- accum[nrow(accum), 2:4]
		smse <- rep(0, 3)
		for (i in 1:nrow(accum)){
			for (j in 1:3)
				smse[j] <- smse[j] + (accum[i,j+1] - target[j])^2
		}
		smse <- smse / nrow(accum)
		smse <- list(Btotal=smse[1], Brepl=smse[2], Brich=smse[3])
	}
	return(smse)
}

#' Slope of accumulation curves.
#' @description This is similar to the first derivative of the curves at each of its points.
#' @param accum A matrix resulting from the alpha.accum or beta.accum functions (samples x diversity values).
#' @details The slope of an accumulation curve, of either observed or estimated diversity, allows verifying if the asymptote has been reached (Cardoso et al. 2011).
#' This is an indication of either the completeness of the inventory (low final slopes of the observed curve indicate high completeness) or reliability of the estimators (stability of the slope around a value of 0 along the curve indicates reliability).
#' @return A matrix of samples x slope values.
#' @references Cardoso, P., Pekar, S., Jocque, R. & Coddington, J.A. (2011) Global patterns of guild composition and functional diversity of spiders. PLoS One, 6, e21710.
#' @examples comm1 <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
#' comm2 <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' acc.alpha = alpha.accum(comm1)
#' slope(acc.alpha)
#' acc.beta = beta.accum(comm1, comm2, tree)
#' slope(acc.beta)
slope <- function(accum){
	if(ncol(accum) > 10){																#if alpha
		sl <- accum[,-2]
		accum <- rbind(rep(0,ncol(accum)), accum)
		for (i in 1:nrow(sl)){
			sl[i,1] <- i
			for (j in 2:ncol(sl)){
				sl[i,j] <- (accum[i+1,j+1]-accum[i,j+1])/(accum[i+1,2]-accum[i,2])
			}
		}
	} else {																						#if beta..
		sl <- accum
		sl[1,] <- 0
		sl[1,1] <- 1
		for (i in 2:nrow(sl)){
			for (j in 2:ncol(sl)){
				sl[i,j] <- (accum[i,j]-accum[i-1,j])
			}
		}
	}
	return(sl)
}


#' Sample data of spiders in Arrabida (Portugal)
#'
#' A dataset containing the abundance of 338 spider species in each of 320 samples. Details are described in:
#' Cardoso, P., Gaspar, C., Pereira, L.C., Silva, I., Henriques, S.S., Silva, R.R. & Sousa, P. (2008) Assessing spider species richness and composition in Mediterranean cork oak forests. Acta Oecologica, 33: 114-127.
#' 
#' @docType data
#' @keywords datasets
#' @name arrabida
#' @usage data(arrabida)
#' @format A data frame with 320 rows and 338 variables
NULL

#' Sample data of spiders in Geres (Portugal)
#'
#' A dataset containing the abundance of 338 spider species in each of 320 samples. Details are described in:
#' Cardoso, P., Scharff, N., Gaspar, C., Henriques, S.S., Carvalho, R., Castro, P.H., Schmidt, J.B., Silva, I., Szuts, T., Castro, A. & Crespo, L.C. (2008) Rapid biodiversity assessment of spiders (Araneae) using semi-quantitative sampling: a case study in a Mediterranean forest. Insect Conservation and Diversity, 1: 71-84.
#' 
#' @docType data
#' @keywords datasets
#' @name geres
#' @usage data(geres)
#' @format A data frame with 320 rows and 338 variables
NULL

#' Sample data of spiders in Guadiana (Portugal)
#'
#' A dataset containing the abundance of 338 spider species in each of 320 samples. Details are described in:
#' Cardoso, P., Henriques, S.S., Gaspar, C., Crespo, L.C., Carvalho, R., Schmidt, J.B., Sousa, P. & Szuts, T. (2009) Species richness and composition assessment of spiders in a Mediterranean scrubland. Journal of Insect Conservation, 13: 45-55.
#' 
#' @docType data
#' @keywords datasets
#' @name guadiana
#' @usage data(guadiana)
#' @format A data frame with 192 rows and 338 variables
NULL

#' Functional tree for 338 species of spiders
#'
#' A dataset representing the functional tree for 338 species of spiders captured in Portugal.
#' For each species were recorded: average size, type of web, type of hunting, stenophagy, vertical stratification in vegetation and circadial activity. Details are described in:
#' Cardoso, P., Pekar, S., Jocque, R. & Coddington, J.A. (2011) Global patterns of guild composition and functional diversity of spiders. PLoS One, 6: e21710.
#' 
#' @docType data
#' @keywords datasets
#' @name functree
#' @usage data(functree)
#' @format An hclust object with 338 species
NULL

#' Phylogenetic tree for 338 species of spiders
#'
#' A dataset representing an approximation to the phylogenetic tree for 338 species of spiders captured in Portugal.
#' The tree is based on the linnean hierarchy, with different suborders separated by 1 unit, families by 0.75, genera by 0.5 and species by 0.25.
#' 
#' @docType data
#' @keywords datasets
#' @name phylotree
#' @usage data(phylotree)
#' @format An hclust object with 338 species
NULL
