#
#	Watson-Williams test for homogeneity
#
#	Allows to compare mean angles in two or more samples.
#	Equivalent, for angles, of an ANOVA/Kruskal-Wallis test.
#
#	Based on
#		Circular statistics in biology, Batschelet, E (1981)
#		§6.2, p99
#
# (c) Copyright 2010 Jean-Olivier Irisson. GNU General Public License
#
#------------------------------------------------------------

# Generic function
watson.williams.test <- function(x, ...) {
	UseMethod("watson.williams.test", x)
}

# Default method, for an angle vector and a grouping vector
# similar to equal.kappa.test
watson.williams.test.default <- function(x, group, data.name=NULL) {	
	# generate a data name of not already done
	if (is.null(data.name)) {
		data.name = paste(deparse(substitute(x)), "by", deparse(substitute(group)))
	}
	
	# check arguments
	ok <- complete.cases(x, group)
	x <- x[ok]
	group <- group[ok]
	if (length(x)==0 | length(table(group)) < 2) {
		warning("No observations or no groups (at least after removing missing values)")
		return(NULL)
	}
		
	# compute concentration parameters and check assumptions
	kt <-  equal.kappa.test(x, group)
	if (kt$p.value < 0.05) {
		warning("Concentration parameters: ", format(kt$kappa, digits=5) ," not equal between groups. The test is not applicable")
	} else if ( any(kt$kappa < 2) ) {
		warning("Concentration parameter: ", format(kt$kappa.all, digits=5)," < 2. The test is not applicable")
	}
	# TODO : also check that distributions conform to Von Mises?

	# number of groups
	group <- as.factor(group)
	k <- nlevels(group)
	# total sample size
	n <- length(x)
	# sample size per group
	group <- sort(group)
	y <- group[-1] != group[-n]
	i <- c(which(y), n)
	ns <- diff(c(0, i))
	# correction factor
	g <- 1 + 3 / (8 * kt$kappa.all)
	# sum of resultant vectors lengths
	sRi <- sum(kt$rho * ns)
	# total resultant vector length
	R <- kt$rho.all * n
	
	statistic <- g * ((n - k) * (sRi - R)) / ((k - 1) * (n - sRi))
	
	p.value <- df(statistic, k-1, n-k)

	# compute estimates of means
	means <- tapply(x, group, mean.circular)

	# return result
	result <- list(
		method = "Watson-Williams test for homogeneity of means",
		data.name = data.name,
		parameter = c(df1=k-1, df2=n-k),
		statistic = c(F=statistic),
		p.value = p.value,
		estimate = means
	)
	class(result) = "htest"
	return(result)
}

# Method for a list
watson.williams.test.list <- function(x) {
	# get data name
	data.name <- deparse(substitute(x))
	# convert into x and group
	k <- length(x)
	ns <- apply(x, length)
	x <- do.call(c, x)
	# NB: unlist() removes the circular attributes here
	group <- rep(1:k, times=ns)
	watson.williams.test.default(x, group, data.name)
}

# Method for a formula
watson.williams.test.formula <- function(formula, data) {
	# convert into x and group
	d <- model.frame(as.formula(formula), data)
	# get data name
	data.name <- paste(names(d), collapse=" by ")
	watson.williams.test.default(d[,1], d[,2], data.name)
}

#------------------------------------------------------------
#	Test data
#------------------------------------------------------------
# 
# xn = x = circular( c(rep(c(-20, -10, 0), c(1,7,2)), rep(c(-10, 0, 10, 20), c(3,3,3,1))), units="degrees", template="geographics")
# group = rep(c(1,2), each=10)
# watson.williams.test(xn, group)
# 
# xl = split(xn, group)
# watson.williams.test(xl)
# 
# xd = data.frame(group=group, angles=x)
# watson.williams.test(angles ~ group, xd)
# 
