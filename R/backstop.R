
snp.rhs.tests154 = function (formula, family="binomial", link, weights, subset,
          data=parent.frame(), snp.data,
          tests=NULL, robust=FALSE,
          control=glm.test.control(maxit=20, epsilon=1.e-4, R2Max=0.98),
          allow.missing=0.01) {
  
  call <- match.call()

  # Family and link function
  
  fam <- pmatch(tolower(family),
                c("binomial", "poisson", "gaussian", "gamma"))
  if (is.na(fam))
    stop("Unrecognized family argument")
  if (fam==0)
    stop("Ambiguous family argument")

  if (missing(link))
    lnk <- fam # Canonical link is default
  else {
    lnk <- pmatch(tolower(lnk), c("logit", "log", "identity", "inverse"))
    if (is.na(fam))
      stop("Unrecognized link argument")
    if (lnk==0)
      stop("Ambiguous link argument")
  }
  
  if (missing(data))
    data <- environment(formula)
  temp <- c("", "formula", "data", "weights", "subset")
  m <- call[match(temp, names(call), nomatch=0)]
  special <- c("strata", "cluster")
    
  Terms <- if (missing(data))
    terms(formula, special)
  else terms(formula, special, data=data)

  # Separate intercept term from rest of model formula

  INTERCEPT <- attr(Terms, "intercept")
  if (!INTERCEPT) # put one back temporarily
    attr(Terms, "intercept") <- 1

  # Model frame for null model

  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())

  # Find common rows and mapping to line up  snp.data

  if (is(snp.data, "snp.matrix")) 
    snames <- rownames(snp.data)
  else
    stop("snp.data must be stored as a snp.data.frame")
  
  mnames <- rownames(m)
  map <- match(mnames, snames, nomatch=0)
  use <- map>0
  N <- sum(use)
  if (N==0)
    stop("no observations")
  

  # Response vector

  if(is.factor(m[[1]])){
    # This is just to avoid the warning that conversion
    # to numeric(=double) from factor is ignored.
    Y <- model.response(m, "any")
  } else {
    Y <- model.response(m, "numeric")
  }
  if (fam==1) { # Binomial family
    Y <- as.factor(Y)
    if (nlevels(Y)==2)
      Y <- as.numeric(Y)-1
    else
      stop("Response is not on two levels")
  }

  # Model matrix -- minus any cluster or strata variables

  strats <- attr(Terms, "specials")$strata
  clust <- attr(Terms, "specials")$cluster
  dropx <- NULL
  if (length(clust)) {
    # Some check on valid clusters?
    warning("cluster(...) not yet implemented")
    if (missing(robust))
      robust <- TRUE
    tempc <- untangle.specials(Terms, "cluster", 1)
    clust <- strata(m[, tempc$vars], shortlabel=TRUE)
    dropx <- tempc$terms
  }
  else
    clust <- NULL
  if (length(strats)) {
    temps <- untangle.specials(Terms, "strata", 1)
    strats <- strata(m[, temps$vars], shortlabel=TRUE)
    dropx <- c(dropx, temps$terms)
  }
  else
    strats <- NULL

  if (length(dropx)) 
    newTerms <- Terms[-dropx]
  else
    newTerms <- Terms
  X <- model.matrix(newTerms, m)

  # Remove unit column vector
  
  ones <- match("(Intercept)", colnames(X), nomatch=0)
  if (ones)
    X <- X[, -ones, drop=FALSE]
  if (ncol(X) == 0)
    X <- NULL
  
  # Prior weights

  if (!is.null(model.weights(m))) {
    cat("One\n")
    weights <- model.weights(m)
    if (any(weights<0))
      stop("negative prior weights not allowed")
  }
  else
    weights <- NULL

  # Offsets

  offset <- model.offset(m)
  if(!is.null(offset)) {
    offset <- offset
    stop("Can't deal with offsets at the moment")
  }

  # Sort and reshape arrays if necessary

  if (N != length(use)) {
    Y <- Y[use]
    if (!is.null(X))
      X <- X[use,,drop=FALSE]
    if (!is.null(weights))
      weights <- weights[use]
    if (!is.null(strats))
      strats <- strats[use]
    if (!is.null(clust))
      clust <- clust[use]
    snp.data <- snp.data[map,,drop=FALSE]
  }
  else if (any(map != 1:N))
    snp.data <- snp.data[map,,drop=FALSE]
    
    

  # Coerce tests argument to correct form #
  if(is.null(tests))
    tests <- as.integer(1:ncol(snp.data))

  # Private utility function to 
  # match character array  to list of names
  .col.numbers <- function(inc, names) {
    res <- match(inc, names, nomatch=NA)
    if (any(is.na(res))) {
      warning("Unmatched SNP names in tests argument")
      res <- tests[!is.na(res)]
    }
    res
  }
  
  if(is.character(tests)) {
    tests <- .col.numbers(tests, colnames(snp.data))
  } else {
    if(is.list(tests)) {
      first <- tests[[1]]
      if (is.character(first)) {
        tests <- lapply(tests, .col.numbers, colnames(snp.data))
      } else {
        if(!is.numeric(first))
          stop("illegal first element of tests argument")
      }
      tests<-lapply(tests,as.integer)
    } else {
      if(!(is.integer(tests) || is.null(tests)))
        stop("illegal tests argument")
    }
  }
  
  .Call("snp_rhs_score154", Y, fam, lnk, X, strats, snp.data, weights,
        tests, robust, clust, control, allow.missing, PACKAGE="GGBase")
}

