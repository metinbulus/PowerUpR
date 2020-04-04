
mdes.ira1r1 <- function(power=.80, alpha=.05, two.tailed=TRUE,
                             p=.50, g1=0, r21=0, n){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- n-g1-2
  SSE <- sqrt((1-r21)/(p*(1-p)*n))

  mdes <- .mdes.fun(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.mdes(effect = "main", power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)
  mdes.out <- list(fun = "mdes.ira1r1",
                   parms = list(power=power, alpha=alpha, two.tailed=two.tailed,
                                p=p, r21=r21, g1=g1,
                                n=n),
                   df = df,
                   ncp = mdes[1]/SSE,
                   mdes = mdes)
  class(mdes.out) <- c("main", "mdes")
  return(invisible(mdes.out))
}
# example
# mdes.ira1r1(n=200)
mdes.ira <- mdes.ira1r1

power.ira1r1 <- function(es=.25, alpha=.05, two.tailed=TRUE,
                              p=.50, g1=0, r21=0, n){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- n-g1-2
  SSE <- sqrt((1-r21)/(p*(1-p)*n))

  power <- .power.fun(es = es, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.power(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, es = es)
  power.out <-  list(fun = "power.ira1r1",
                     parms = list(es=es, alpha=alpha, two.tailed=two.tailed,
                                 p=p, r21=r21, g1=g1,
                                 n=n),
                     df = df,
                     ncp = es/SSE,
                     power = power)
  class(power.out) <- c("main", "power")
  return(invisible(power.out))
}
# example
# power.ira1r1(n=200)
power.ira <- power.ira1r1

mrss.ira1r1 <- function(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
                        n0=10, tol=.10,
                        p=.50, g1=0, r21=0){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  i <- 0
  conv <- FALSE
  while(i<=100 & conv==FALSE){
    df <- n0-g1-2
    if(df<= 0 | is.infinite(df)){break}
    T1 <- ifelse(two.tailed==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    n1 <- (M/es)^2 * ((1-r21)/(p*(1-p)))
    if(abs(n1-n0)<tol){conv <- TRUE}
    n0 <- (n1+n0)/2
    i <- i+1
  }
  n <- round(ifelse(df>0,round(n0),NA))

  n.out <-  list(fun = "mrss.ira1r1",
                 parms = list(es=es, power=power, alpha=alpha, two.tailed=two.tailed,
                              n0=n0, tol=tol,
                              p=p, r21=r21, g1=g1),
                 df=df,
                 ncp = M,
                 n = n)
  class(n.out) <- c("main", "mrss")
  cat("n =", n, "\n")
  return(invisible(n.out))
}
# mrss.ira1r1()
mrss.ira <- mrss.ira1r1

mdes.bira2f1 <- function(power=.80, alpha=.05, two.tailed=TRUE,
                         p=.50, g1=0, r21=0, n, J){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- J * (n - 2) - g1
  SSE <- sqrt((1-r21)/(p*(1-p)*J*n))

  mdes <- .mdes.fun(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.mdes(effect = "main", power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)
  mdes.out <- list(fun = "mdes.bira2f1",
                   parms = list(power=power, alpha=alpha, two.tailed=two.tailed,
                                p=p, r21=r21, g1=g1,
                                n=n, J=J),
                   df = df,
                   ncp = mdes[1]/SSE,
                   mdes = mdes)
  class(mdes.out) <- c("main", "mdes")
  return(invisible(mdes.out))
}
# example
# mdes.bira2f1(n=55, J=3)

power.bira2f1 <- function(es=.25, alpha=.05, two.tailed=TRUE,
                          p=.50, g1=0, r21=0, n, J){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- J * (n - 2) - g1
  SSE <- sqrt((1-r21)/(p*(1-p)*J*n))

  power <- .power.fun(es = es, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.power(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, es = es)
  power.out <-  list(fun = "power.bira2f1",
                     parms = list(es=es, alpha=alpha, two.tailed=two.tailed,
                                  p=p, r21=r21, g1=g1,
                                  n=n, J=J),
                     df = df,
                     ncp = es/SSE,
                     power = power)
  class(power.out) <- c("main", "power")
  return(invisible(power.out))
}
# example
# power.bira2f1(n=55, J=3)

mrss.bira2f1 <- function(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
                         J, n0=10, tol=.10,
                         p=.50, g1=0, r21=0){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  i <- 0
  conv <- FALSE
  while(i<=100 & conv==FALSE){
    df <- J*(n0-2)-g1
    if(df<= 0 | is.infinite(df)){break}
    T1 <- ifelse(two.tailed==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    n1 <- (M/es)^2 * ((1-r21)/(p*(1-p)*J))
    if(abs(n1-n0)<tol){conv <- TRUE}
    n0 <- (n1+n0)/2
    i <- i+1
  }
  n <- ifelse(df>0,round(n0),NA)

  mrss.out <-  list(fun = "mrss.bira2f1",
                    parms = list(es=es, power=power, alpha=alpha, two.tailed=two.tailed,
                                 J=J, n0=n0, tol=tol,
                                 p=p, r21=r21, g1=g1),
                    df = df,
                    ncp = M,
                    n = n)
  class(mrss.out) <- c("main", "mrss")
  cat("n =", n, "(per block)\n")
  return(invisible(mrss.out))
}
# example
# mrss.bira2f1(J=5)

mdes.bira2c1 <- function(power=.80, alpha=.05, two.tailed=TRUE,
                         p=.50, g1=0, r21=0, n, J){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- J * (n - 1) - g1 - 1
  SSE <- sqrt((1-r21)/(p*(1-p)*J*n))

  mdes <- .mdes.fun(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.mdes(effect = "main", power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)
  mdes.out <- list(fun = "mdes.bira2c1",
                   parms = list(power=power, alpha=alpha, two.tailed=two.tailed,
                                p=p, r21=r21, g1=g1,
                                n=n, J=J),
                   df = df,
                   ncp = mdes[1]/SSE,
                   mdes = mdes)
  class(mdes.out) <- c("main", "mdes")
  return(invisible(mdes.out))
}
# example
# mdes.bira2c1(n=55, J=14)

power.bira2c1 <- function(es=.25, alpha=.05, two.tailed=TRUE,
                          p=.50, g1=0, r21=0, n, J){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- J * (n - 1) - g1 - 1
  SSE <- sqrt((1-r21)/(p*(1-p)*J*n))

  power <- .power.fun(es = es, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.power(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, es = es)
  power.out <-  list(fun = "power.bira2c1",
                     parms = list(es=es, alpha=alpha, two.tailed=two.tailed,
                                  p=p, r21=r21, g1=g1,
                                  n=n, J=J),
                     df = df,
                     ncp = es/SSE,
                     power = power)
  class(power.out) <- c("main", "power")
  return(invisible(power.out))
}
# example
# power.bira2c1(n=55, J=14)

mrss.bira2c1 <- function(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
                         J, n0=10, tol=.10,
                         p=.50, g1=0, r21=0){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  i <- 0
  conv <- FALSE
  while(i<=100 & conv==FALSE){
    df <- J*(n0-1)-g1-1
    if(df<= 0 | is.infinite(df)){break}
    T1 <- ifelse(two.tailed==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    n1 <- (M/es)^2 * ((1-r21)/(p*(1-p)*J))
    if(abs(n1-n0)<tol){conv <- TRUE}
    n0 <- (n1+n0)/2
    i <- i+1
  }
  n <- ifelse(df>0,round(n0),NA)

  mrss.out <-  list(fun = "mrss.bira2c1",
                    parms = list(es=es, power=power, alpha=alpha, two.tailed=two.tailed,
                                 J=J, n0=n0, tol=tol,
                                 p=p, r21=r21, g1=g1),
                    df = df,
                    ncp = M,
                    n = n)
  class(mrss.out) <- c("main", "mrss")
  cat("n =", n, "(per block)\n")
  return(invisible(mrss.out))
}
# example
# mrss.bira2c1(J=5)
