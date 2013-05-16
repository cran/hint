###########################################################################
#
# Functions for calculating hypergeometric intersection distributions.
#
# Author: Alex T. Kalinka (alex.t.kalinka@gmail.com)
#
###########################################################################


#### Main user-callable functions that do the dispatching to user-invisible worker functions below. ####
#
# Split into two sets: distribution, quantile, etc functions, and a more formal test function that will
#  take more flexible user input for testing (objects, which will be placed into categories internally).
#
# 1. Density, distribution (probability mass function), quantile, and random generation functions:


.hint.check.params <- function(n, a, b, q)
	{
	if(!a>=0 || !b>=0 || !n>=0 || !q>=0 || !a<=n || !b<=(n+q) || !q<=n || !q>=0){
		stop("the following constraints must be met:\nn > 0\n0 <= a <= n\n0 <= b <= (n + q)\n0 <= q <= n\n", call. = FALSE)
		}
	vmin <- max(a+b-n-min(floor(b/2),q),0)
	vmax <- min(a,b)
	vrange <- vmin:vmax
	return(vrange)
	}



dhint <- function(n, a, b, q = 0, range = NULL, log = FALSE, verbose = TRUE)
	{
	# range is a vector giving intersection sizes for which the user wishes to retrieve probabilities.
	vrange <- .hint.check.params(n, a, b, q)
	if(is.null(range)){
		range <- vrange
		}
	if(q==0){
		rn <- intersect(range, vrange)
		if(length(rn)==0){
			dist <- data.frame(v=range, p=rep(0,length(range)))
			}
		dist <- .hint.symm.sing(n, a, b, range = rn)
	}else{
		rn <- intersect(range, vrange)
		if(length(rn)==0){
			dist <- data.frame(v=range, p=rep(0,length(range)))
			}
		dist <- .hint.asymm.dup(n, a, b, q, range = rn, verbose = verbose)
		}
	if(log){
		dist[,2] <- log(dist[,2])
		}
	return(dist)
	}


phint <- function(n, a, b, q = 0, vals, upper.tail = TRUE, log.p = FALSE)
	{
	# vals are the values of v for which we want cumulative probabilities.
	dist <- dhint(n, a, b, q)
	if(log.p){
		dist[,2] <- log(dist[,2])
		}
	vrange <- .hint.check.params(n, a, b, q)
	rn <- intersect(vals, vrange)
	if(length(rn)==0){
		pp <- NULL
		for(i in 1:length(vals)){
			if(vals[i] < min(vrange)){
				if(upper.tail){
					pp[i] <- 1
				}else{
					pp[i] <- 0
					}
			}else{
				if(upper.tail){
					pp[i] <- 0
				}else{
					pp[i] <- 1
					}
				}
			}
		pval <- data.frame(v=vals, cum.p=rep(pp, length(vals)))
	}else{
		pv <- NULL
		for(i in 1:length(rn)){
			if(upper.tail){
				inds <- which(dist[,1]==rn[i]):nrow(dist)
			}else{
				inds <- 1:which(dist[,1]==rn[i])
				}
			pv[i] <- sum(dist[inds, 2])
			}
		pval <- data.frame(v=rn, cum.p=pv)
		}
	if(log.p){
		pval[,2] <- log(pval[,2])
		}
	return(pval)
	}


qhint <- function(p, n, a, b, q = 0, upper.tail = TRUE, log.p = FALSE)
	{
	# p is a probability.
	if(!p>=0 || !p<=1){
		stop("p must be between 0 and 1\n", call. = FALSE)
		}
	vrange <- .hint.check.params(n, a, b, q)
	vals <- vrange
	dist <- dhint(n, a, b, q, range=vals)
	pvals <- phint(n, a, b, q, upper.tail=upper.tail, vals=vals)
	inds <- which(pvals[,2]<=p)
	pv <- sum(dist[inds, 2])
	qq <- pvals[max(inds),1]
	if(log.p){
		pv <- log(pv)
		}
	ret <- data.frame(v=qq, cum.p=pv)
	return(ret)
	}


rhint <- function(num = 5, n, a, b, q = 0)
	{
	vrange <- .hint.check.params(n, a, b, q)
	probs <- dhint(n, a, b, q, range=vrange)
	samp <- sample(vrange, num, prob = probs[,2], replace = TRUE)
	return(samp)
	}


############################################################
#
# 2. Functions to perform tests of enrichment or depletion.
#    Includes distance distribution.
#
##

print.hint.test <- function(x, ...)
	# S3 method for generic function "print".
	# x is a "hint.test".
	{
	md <- match("d",names(x$parameters))
	if(is.na(md)){
		rv <- "v"
	}else{
		rv <- "d"
		}
	if(x$alternative == "greater"){
		p.test <- paste("P(X >= ",rv,")",sep="")
	}else if(x$alternative == "less"){
		p.test <- paste("P(X <= ",rv,")",sep="")
	}else{
		p.test <- "P(two-sided)"
		}
	cat("\nHypergeometric intersection test\n\nParameters:\n")
	print(x$parameters)
	cat("\n",p.test," = ",x$p.value,"\n\n")
	}


hint.test <- function(cats, draw1, draw2, alternative = "greater")
	{
	# cats is a data frame or character matrix of category names and numbers for each urn:
	#    cats urn.1 urn.2
	#	a     1     1
	#	b     1     1
	#	c     1     1
	#	d     1     1
	#	e     1     1
	# draw1 and draw2 are samples of the categories from urn 1 and 2 respectively:
	# [1] "a" "b" "c" "d" "e"
	#
	# alternative can be: greater, less, or two.sided.
	# two-sided defined by doubling: 2*min(p.upper, p.lower)
	#
	n <- length(unique(cats[,1]))
	a <- length(draw1)
	b <- length(draw2)
	q <- length(which(cats[,3]==2))
	v <- length(intersect(draw1,draw2))
	if(alternative == "greater"){
		ut <- TRUE
	}else if(alternative == "less"){
		ut <- FALSE
	}else if(alternative == "two.sided"){
		p1 <- phint(n, a, b, q, vals = v, upper.tail = TRUE)
		p2 <- phint(n, a, b, q, vals = v, upper.tail = FALSE)
	}else{
		stop("alternative hypothesis must be one of: greater, less, or two.sided\n", call. = FALSE)
		}
	if(alternative == "two.sided"){
		pval <- 2*min(p1, p2)
	}else{
		pval <- phint(n, a, b, q, vals = v, upper.tail = ut)[,2]
		}
	ret <- list()
	params <- c(n,a,b,q,v)
	names(params) <- c("n","a","b","q","v")
	ret$parameters <- params
	ret$p.value <- pval
	ret$alternative <- alternative
	class(ret) <- "hint.test"
	return(ret)
	}


hint.dist.test <- function(d, n1, a1, b1, n2, a2, b2, q1 = 0, q2 = 0, alternative = "greater")
	{
	# d is the distance we wish to test.
	# alternative can be: greater, less, or two.sided.
	# two-sided defined by doubling: 2*min(p.upper, p.lower)
	#
	vrange1 <- .hint.check.params(n1, a1, b1, q1)
	vrange2 <- .hint.check.params(n2, a2, b2, q2)
	dist <- .hint.dist.distr(n1, a1, b1, n2, a2, b2, q1, q2)
	if(alternative == "greater"){
		inds <- which(dist[,1] >= d)
		if(length(inds)==0){
			pval <- 0
		}else{
			pval <- sum(dist[inds, 2])
			}
	}else if(alternative == "less"){
		inds <- which(dist[,1] <= d)
		if(length(inds)==0){
			pval <- 0
		}else{
			pval <- sum(dist[inds, 2])
			}
	}else if(alternative == "two.sided"){
		inds1 <- which(dist[,1] >= d)
		inds2 <- which(dist[,1] <= d)
		if(length(inds1)==0){
			p1 <- 0
		}else{
			p1 <- sum(dist[inds1, 2])
			}
		if(length(inds2)==0){
			p2 <- 0
		}else{
			p2 <- sum(dist[inds2, 2])
			}
		pval <- 2*min(p1, p2)
	}else{
		stop("alternative hypothesis must be one of: greater, less, or two.sided\n", call. = FALSE)
		}
	ret <- list()
	params <- c(d,n1,a1,b1,q1,n2,a2,b2,q2)
	names(params) <- c("d","n1","a1","b1","q1","n2","a2","b2","q2")
	ret$parameters <- params
	ret$p.value <- pval
	ret$alternative <- alternative
	class(ret) <- "hint.test"
	return(ret)
	}



#### Worker functions: ####


# Symmetrical, singletons case.
.hint.symm.sing <- function(n, a, b, range)
	{
	dist <- NULL
	for(i in range){
		dist <- append(dist, exp(-1*lgamma((i+1)) - lgamma((b-i+1)) - lgamma((a-i+1)) - lgamma((n-b-a+i+1)) - lgamma((n+1)) + lgamma((b+1)) + lgamma((n-b+1)) + lgamma((a+1)) + lgamma((n-a+1))) )
		}
	ret <- data.frame(v=range, p=dist)
	return(ret)
	}


# Asymmetrical, duplicates case.
.hint.asymm.dup <- function(n, a, b, q, range, verbose = TRUE)
	{
	len <- length(range)
	dist <- rep(0, len)
	out <- .C("inters_distr_dup_opt",as.integer(n),as.integer(a),as.integer(b),as.integer(q),dist=as.double(dist),as.integer(range),as.integer(len),as.logical(verbose))
	if(verbose){
		cat("\n")
		}
	dist <- out$dist
	ret <- data.frame(v=range, p=dist)
	return(ret)
	}


# Distribution of distances between independent intersection distributions.
.hint.dist.distr <- function(n1, a1, b1, n2, a2, b2, q1 = 0, q2 = 0)
	{
	dist <- NULL
	m1 <- min(a1,b1)
	m2 <- min(a2,b2)
	m3 <- max(a1+b1-n1-min(floor(b1/2),q1),0)
	m4 <- max(a2+b2-n2-min(floor(b2/2),q2),0)
	if(m1 > m2){
		mxdist <- m1 - m4
		mndist <- max(m3 - m2,0)
	}else{
		mxdist <- m2 - m3
		mndist <- max(m4 - m1,0)
		}
	for(d in mndist:mxdist){
		cat("   Completed...",round((d/mxdist)*100,3),"%\r",collapse="")
		# Enumerate intersection pairs that produce distance d.
		rm <- m2 - d
		sm <- m1 - d
		A <- min(m1,rm)+1
		B <- min(m2,sm)+1
		if(A<0 && B>0){
			Q <- B
		}else if(A>0 && B<0 || d == 0){
			Q <- A
		}else{
			Q <- A + B
			}
		if(Q <= 0){dist[(d+1)] <- 0; next}
		pairs <- matrix(0,Q,2)
		lim <- A
		one <- 0; if(d==0){two <- 0}else{two <- d}
		if(lim>0){
			for(i in 1:lim){
				pairs[i,] <- c(one,two)
				one <- one + 1
				two <- two + 1
				}
			}
		one <- d; two <- 0
		if(d!=0){
			lim2 <- B
			if(lim<Q && lim2>0){
				for(i in (lim+1):Q){
					pairs[i,] <- c(one,two)
					one <- one + 1
					two <- two + 1
					}
				}
			}
		# Calculate probability for this d.
		sd <- 0 

		for(i in 1:Q){
			# First distribution.
			v <- pairs[i,1]
			if(v < m3){
				p1 <- 0
			}else{
				p1 <- dhint(n1, a1, b1, q1, range = v, verbose = FALSE)[,2]
				}
			# Second distribution.
			v <- pairs[i,2]
			if(v < m4){
				p2 <- 0
			}else{
				p2 <- dhint(n2, a2, b2, q2, range = v, verbose = FALSE)[,2]
				}
			sd <- sum(sd, p1*p2)
			}
		dist <- append(dist, sd)
		}
	cat("\n")
	ret <- data.frame(d=mndist:mxdist, p=dist)
	return(ret)
	}


### Simulation ###

# Two urns.
# Variable numbers in each category:
sim.hypint <- function(n, a, b, sims = 10000, Na = rep(1,n), Nb = rep(1,n))
	{
	# Na and Nb are vectors with the numbers in each category (if they're variable).
	sdist <- NULL; na <- nb <- NULL
	for(i in 1:n){
		na <- append(na, rep(i,Na[i]))
		nb <- append(nb, rep(i,Nb[i]))
		}

	for(i in 1:sims){
		cat("\r   Completed...",round((i/sims)*100,3),"%",collapse="")
		s1 <- sample(na, a, replace = FALSE)
		s2 <- sample(nb, b, replace = FALSE)
		sdist[i] <- length(intersect(s1,s2))
		}
	cat("\n")
	return(sdist)
	}





