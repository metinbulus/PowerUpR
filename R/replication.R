# mdh: tau^2 / V
# tau^2: variance of effect sizes
# V: average estimate variance
power.rep <- function(k = 2L, mdh = 1/4, mdh.null = 0, alpha = .05) {

  # validity checks
  if(!is.numeric(alpha) | alpha <= 0 | alpha >= 1)
    stop("'alpha' argument must be a numeric value between 0 and 1", call. = FALSE)
  if(mdh < mdh.null)
    stop("'mdh.max' must be greater than or equal to 'mdh.null'", call. = FALSE)
  if(!is.numeric(mdh) |  !is.numeric(mdh.null) | mdh < 0 | mdh.null < 0)
    stop("MDH arguments must be numeric values greater than zero", call. = FALSE)
  if(!is.integer(k) | k <= 0)
    stop("'k' must be a positive integer, e.g. 'k = 2L'", call. = FALSE)

  # calculate power
  c.value <- qchisq(p = 1- alpha, df = k-1, ncp = (k-1) * mdh.null, lower.tail = TRUE)
  power <- 1 - pchisq(c.value, df = k-1, ncp = (k-1) * mdh, lower.tail = TRUE)
  power.out <-  list(fun = "power.rep",
                     parms =list(k=k, mdh=mdh, mdh.null=mdh.null, alpha=alpha),
                     df=k-1, power=power)

  class(power.out) <- c("rep", "power")
  cat(ifelse(mdh.null == 0, "\nExact", "\nApproximate"), " replication",
      "\n Number of studies: ", k,
      "\n MDH for null: ", round(mdh.null, 3),
      "\n MDH for alternative: ", round(mdh, 3), "\n\nStatistical power:",
      "\n-------------------------\n", round(power, 3),
      "\n-------------------------\nDegrees of freedom: ", k-1,
      "\nType I error rate: ", alpha,
      "\nType II error rate: ", round(1 - power, 3), "\n",
      sep = "")
  return(invisible(power.out))
}
# power.rep()

# mdh: tau^2 / V
# tau^2: variance of effect sizes
# V: average estimate variance
mdh.rep <- function(k = 2L, mdh.max = 15, alpha = .05, power = 0.80, mdh.null = 0, step = .001, plot = FALSE) {

  # validity checks
  if(step < 0 | step > .10 | !is.numeric(step))
    stop("Incorrect 'step' argument (0 < 'step' < 0.10)", call. = FALSE)
  if(!is.numeric(power) | !is.numeric(alpha) | power <= 0 | power >= 1 | alpha <= 0 | alpha >= 1)
    stop("'power' or 'alpha' arguments must be numeric values between 0 and 1", call. = FALSE)
  if(mdh.max < mdh.null)
    stop("'mdh.max' must be greater than or equal to'mdh.null'", call. = FALSE)
  if(!is.numeric(mdh.max) |  !is.numeric(mdh.null) | mdh.max < 0 | mdh.null < 0)
    stop("MDH arguments must be numeric values greater than zero", call. = FALSE)
  if(!is.integer(k) | k <= 0)
    stop("'k' must be a positive integer, e.g. 'k = 2L'", call. = FALSE)

  # MDH brute force search
  mdh.seq <- seq(mdh.null, mdh.max, step)
  c.value <- qchisq(p = 1- alpha, df = k-1, ncp = (k-1) * mdh.null, lower.tail = TRUE)
  pwr.seq <- 1 - pchisq(c.value, df = k-1, ncp = (k-1) * mdh.seq, lower.tail = TRUE)
  mdh <- mdh.seq[which(abs(power-pwr.seq) == min(abs(power-pwr.seq)))]

  if(mdh == mdh.max)
    stop("Increase 'mdh.max'", call. = FALSE)

  mdh.out <- list(fun = "mdh.rep",
                  parms =list(k=k, mdh.max=mdh.max, mdh.null=mdh.null, power=power, alpha=alpha),
                  df=k-1, mdh=mdh)

  if(plot) {
    plot(mdh.seq, pwr.seq, col = adjustcolor(4, alpha.f = 0.5), type = "l", lty = 1,
         xlab = "Minimum Detectable Heterogeneity", ylab = "Statistical Power")
    points(mdh, power, pch = 21, bg = adjustcolor(2, alpha.f = 0.5), cex = 1.5)
    abline(v = mdh, h = power, lty = 5, col = adjustcolor(2, alpha.f = 0.5))
  }

  class(mdh.out) <- c("rep", "mdh")
  cat(ifelse(mdh.null == 0, "\nExact", "\nApproximate"), " replication",
      "\n Number of studies: ", k,
      "\n MDH for null: ", round(mdh.null, 3),
      "\n\nMinimum detectable heterogeneity:",
      "\n--------------------------------\n", round(mdh, 3),
      "\n--------------------------------\nDegrees of freedom: ", k-1,
      "\nType I error rate: ", alpha,
      "\nType II error rate: ", round(1 - power, 3), "\n",
      sep = "")
  return(invisible(mdh.out))
}
# mdh.rep()


q.test <- function(es, sse, mdh.null = 0, alpha = .05) {

  # validity checks
  if(!is.numeric(alpha) | alpha <= 0 | alpha >= 1)
    stop("'power' or 'alpha' arguments must be numeric values between 0 and 1", call. = FALSE)
  if(!is.numeric(mdh.null) | mdh.null < 0)
    stop("'mdh.null' argument must be a numeric greater than zero", call. = FALSE)
  if(!is.vector(es) | !is.vector(sse) | !is.numeric(es) | !is.numeric(sse) | length(es) != length(sse))
    stop("'es' and 'sse' arguments must be numeric vectors of the same length", call. = FALSE)

  k <- length(es)
  es_star <- sum(es/sse^2) / sum(1/sse^2)
  #es_star <- weighted.mean(x=es, w=1/sse^2) # weighted mean
  Q <- sum((es - es_star)^2 / sse^2)
  c.value <- qchisq(p = 1 - alpha, df = k-1, ncp = (k-1) * mdh.null)
  p.value <- 1 - pchisq(Q, df = k-1, ncp = (k-1) * mdh.null)

  rep.out <- list(fun = "Q.test",
                  parms = list(mdh.null = mdh.null, es = es, sse = sse),
                  Q = Q, df = k-1, p.value = p.value)

  class(rep.out) <- c("rep", "q.test")
  cat(ifelse(mdh.null == 0, "\nExact", "\nApproximate"), " replication test",
      "\n Number of studies: ", k,
      "\n MDH for null: ", round(mdh.null, 3), "\n\n",
      sep = "")
  print(round(unlist(rep.out[c("Q", "df", "p.value")]),3))

  return(invisible(rep.out))
}
