#
#	Watson-Williams test for homogeneity of means
#
#	Allows to compare mean angles in two or more samples.
#	Equivalent, for angles, of an ANOVA/Kruskal-Wallis test.
#
#	Based on
#		Circular statistics in biology, Batschelet, E (1981)
#		§6.2, p. 99
#		Biostatistical analysis, Zar, J H (1999)
#		§27.4, p. 634
#		Directional statistics, Mardia, K.V. and Jupp, P.E. (2000)
#		p. 135
#
# (c) Copyright 2010-211 Jean-Olivier Irisson
#     GNU General Public License, v.3
#
#------------------------------------------------------------

# Generic function
watson.williams.test <- function(x, ...) {
	UseMethod("watson.williams.test", x)
}

# Default method, for an angle vector and a grouping vector
watson.williams.test.default <- function(x, group, ...) {

	# get data name
	data.name <- paste(deparse(substitute(x)), "by", deparse(substitute(group)))

	# check arguments
	ok <- complete.cases(x, group)
	x <- x[ok]
	group <- group[ok]
	if (length(x)==0 | length(table(group)) < 2) {
		stop("No observations or no groups (at least after removing missing values)")
	}

	# convert everything to the radians/trigonometric case
	if (is.circular(x)) {
		dc <- circularp(x)
		x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
		attr(x, "class") <- attr(x, "circularp") <- NULL
	} else {
		dc <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
	}

	# compute concentration parameters and check assumptions
	kt <- EqualKappaTestRad(x, group)
	# equality of concentration parameters
	if (kt$p.value < 0.05) {
		warning("Concentration parameters: ", format(kt$kappa, digits=5) ," not equal between groups. The test might not be applicable")
	}
	# sufficiently large concentration
	# Batschelet provides kappa > 2 in the two sample case but no indication in the multisample one
	# Zar's conditions are sample size dependent
	# Mardia & Jupp cite Stephens 1972 to justify that kappa >= 1 is OK
	if ( kt$kappa.all < 1 ) {
		warning("Global concentration parameter: ", format(kt$kappa.all, digits=3)," < 1. The test is probably not applicable")
	}
	# TODO : also check that distributions conform to Von Mises?

	result <- WatsonWilliamsTestRad(x, group, kt)
	result$data.name <- data.name

	# convert means back in their original units
	result$estimate <-  conversion.circular(circular(result$estimate), units=dc$units, type=dc$type, template=dc$template, modulo=dc$modulo, zero=dc$zero, dc$rotation)

	return(result)
}

# Method for a list
watson.williams.test.list <- function(x, ...) {
	# fecth or fill list names
	k <- length(x)
	if (is.null(names(x))) {
		names(x) <- 1:k
	}
	# get data name
	data.name <- paste(names(x), collapse=" and ")

	# convert into x and group
	ns <- lapply(x, length)
	group <- rep(names(x), times=ns)
	x <- do.call("c", x)
	# NB: unlist() removes the circular attributes here

	# call default method
	result <- watson.williams.test.default(x, group)
	result$data.name <- data.name

	return(result)
}

# Method for a formula
watson.williams.test.formula <- function(formula, data, ...) {
	# convert into x and group
	d <- model.frame(as.formula(formula), data)

	# get data name
	data.name <- paste(names(d), collapse=" by ")

	# call default method
	result <- watson.williams.test.default(d[,1], d[,2])
	result$data.name <- data.name

	return(result)
}


# Computation in the usual trigonometric space
WatsonWilliamsTestRad <- function(x, group, kt) {

	# number of groups
	group <- as.factor(group)
	k <- nlevels(group)
	# total sample size
	n <- length(x)
	# sample size per group
	ns <- as.numeric(table(group))
	# correction factor
	g <- 1 + 3 / (8 * kt$kappa.all)
	# sum of resultant vectors lengths
	sRi <- sum(kt$rho * ns)
	# total resultant vector length
	R <- kt$rho.all * n

	statistic <- g * ((n - k) * (sRi - R)) / ((k - 1) * (n - sRi))

	p.value <- pf(statistic, k-1, n-k, lower.tail=FALSE)

	# compute estimates of means
	means <- tapply(x, group, MeanCircularRad)
	names(means) <- paste("mean of", names(means))

	# return result
	result <- list(
		method = "Watson-Williams test for homogeneity of means",
		parameter = c(df1=k-1, df2=n-k),
		statistic = c(F=statistic),
		p.value = p.value,
		estimate = means
	)
	class(result) <- "htest"

	return(result)
}