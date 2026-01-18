#' @import pracma
NULL
#' R function that calculates the probability density of minimum of NIND random variables (PDF)
#'
#' @param y_min the value of y_min
#' @param distribution the distribution of marker value follows Normal or Gamma
#' @param arg1 if distribution is normal input mean parameters of all subclasses in a vector, if gamma input shape parameters
#' @param arg2 if distribution is normal input variance parameter, if gamma input rate parameters
#' @return The probability density of minimum of NIND random variables
#' @importFrom pracma perms
#' @export
f_order_min <- function(y_min,distribution,arg1,arg2){
  k=length(arg1)
  perm <- perms(c(1:k))

  if(distribution=="Normal"){
    mu=arg1
    sd=arg2

  if(k>1){f_min <- function(y_min){
    1/factorial(k-1)*eval(parse(text=paste(sapply(1:dim(perm)[1], function(j) paste(paste(sapply(perm[j,1],function(a) substitute(dnorm(y_min,mu[i],sd[i]),list(i=a)))),
                                                                                    paste(sapply(perm[j,2:k],function(a) substitute(pnorm(y_min,mu[i],sd[i],lower.tail = F),list(i=a))),collapse="*"),sep='*')),collapse='+')))

  }}
  if(k==1){f_min <- function(y_min) dnorm(y_min,mu,sd)}
  }

  if(distribution=="Gamma"){
    shape=arg1
    rate=arg2

    if(k>1){f_min <- function(y_min){
      1/factorial(k-1)*eval(parse(text=paste(sapply(1:dim(perm)[1], function(j) paste(paste(sapply(perm[j,1],function(a) substitute(dgamma(y_min,shape[i],rate[i]),list(i=a)))),
                                                                                      paste(sapply(perm[j,2:k],function(a) substitute(pgamma(y_min,shape[i],rate[i],lower.tail = F),list(i=a))),collapse="*"),sep='*')),collapse='+')))

    }}
    if(k==1){f_min <- function(y_min) dgamma(y_min,shape,rate)}
    }
  return(f_min(y_min))
}

#' R function that calculates the probability density of maximum of NIND random variables (PDF)
#'
#' @param y_max the value of y_max
#' @param distribution the distribution of marker value follows Normal or Gamma
#' @param arg1 if distribution is normal input mean parameters of all subclasses in a vector, if gamma input shape parameters
#' @param arg2 if distribution is normal input variance parameter, if gamma input rate parameters
#' @return The probability density of maximum of random variables
#' @export
f_order_max <- function(y_max,distribution,arg1,arg2){
  k=length(arg1)
  perm <- perms(c(1:k))

  if(distribution=="Normal"){
    mu=arg1
    sd=arg2

  if(k>1){
    f_max <- function(y_max){
      1/factorial(k-1)*eval(parse(text=paste(sapply(1:dim(perm)[1], function(j) paste(paste(sapply(perm[j,1],function(a) substitute(dnorm(y_max,mu[i],sd[i]),list(i=a)))),
                                                                                      paste(sapply(perm[j,2:k],function(a) substitute(pnorm(y_max,mu[i],sd[i]),list(i=a))),collapse="*"),sep='*')),collapse='+')))}
  }
  if(k==1){
    f_max <- function(y_max) dnorm(y_max,mu,sd)
  }
}

  if(distribution=="Gamma"){
    shape=arg1
    rate=arg2

    if(k>1){
      f_max <- function(y_max){
        1/factorial(k-1)*eval(parse(text=paste(sapply(1:dim(perm)[1], function(j) paste(paste(sapply(perm[j,1],function(a) substitute(dgamma(y_max,shape[i],rate[i]),list(i=a)))),
                                                                                        paste(sapply(perm[j,2:k],function(a) substitute(pgamma(y_max,shape[i],rate[i]),list(i=a))),collapse="*"),sep='*')),collapse='+')))}
    }
    if(k==1){
      f_max <- function(y_max) dgamma(y_max,shape,rate)
    }
  }
  return(f_max(y_max))
}

#' R function that calculates the partial of joint probability of min and max over max of NIND random variables
#'
#' @param y_max the value of y_max
#' @param y_min the value of y_min
#' @param distribution the distribution of marker value follows Normal or Gamma
#' @param arg1 if distribution is normal input mean parameters of all subclasses in a vector, if gamma input shape parameters
#' @param arg2 if distribution is normal input variance parameter, if gamma input rate parameters
#' @return The partial of joint probablity of min and max over max
#' @export

cdf_min_max_partial <- function(y_min,y_max,distribution,arg1,arg2){
    k <- length(arg1)
    perm <- perms(c(1:k))

    if(distribution=="Normal"){
      mu=arg1
      sd=arg2

    (y_min<=y_max)*(1/factorial(k-1)*


      eval(parse(text=paste(sapply(1:dim(perm)[1], function(j) paste('(',paste(paste(sapply(perm[j,1],function(a) substitute(dnorm(y_max,mu[i],sd[i]),list(i=a)))),
                                                                               paste(sapply(perm[j,2:k],function(a) substitute(pnorm(y_max,mu[i],sd[i]),list(i=a))),collapse="*"),sep='*'),'-',
                                                                     paste(paste(sapply(perm[j,1],function(a) substitute(dnorm(y_max,mu[i],sd[i]),list(i=a)))),
                                                                           paste(sapply(perm[j,2:k],function(a) substitute((pnorm(y_max,mu[i],sd[i])-pnorm(y_min,mu[i],sd[i])),list(i=a))),collapse="*"),sep='*'),')',sep='')),collapse='+')

    )))
  }

   if(distribution=="Gamma"){
     shape=arg1
     rate=arg2

     (y_min<=y_max)*(1/factorial(k-1)*


       eval(parse(text=paste(sapply(1:dim(perm)[1], function(j) paste('(',paste(paste(sapply(perm[j,1],function(a) substitute(dgamma(y_max,shape[i],rate[i]),list(i=a)))),
                                                                                paste(sapply(perm[j,2:k],function(a) substitute(pgamma(y_max,shape[i],rate[i]),list(i=a))),collapse="*"),sep='*'),'-',
                                                                      paste(paste(sapply(perm[j,1],function(a) substitute(dgamma(y_max,shape[i],rate[i]),list(i=a)))),
                                                                            paste(sapply(perm[j,2:k],function(a) substitute((pgamma(y_max,shape[i],rate[i])-pgamma(y_min,shape[i],rate[i])),list(i=a))),collapse="*"),sep='*'),')',sep='')),collapse='+')

       )))
   }
}

#' R function that calculates the conditional probability of minimum greater than y_min given maximum equals to y_max of random variables (upper tail probability of minimum given maximum)
#'
#' @param y_min the value of y_min
#' @param y_max the value of y_max
#' @param distribution the distribution of marker value follows Normal or Gamma
#' @param arg1 if distribution is normal input mean parameters of all subclasses in a vector, if gamma input shape parameters
#' @param arg2 if distribution is normal input variance parameter, if gamma input rate parameters
#' @return The conditional probability of minimum given maximum of random variables
#' @export
cdf_min_given_max_partial_upper <- function(y_min,y_max,distribution,arg1,arg2){
  k <- length(arg1)
  perm <- perms(c(1:k))

  if(k==1){func <- function(y_min,y_max) 1}

  if(k>1){if(f_order_max(y_max,distribution,arg1,arg2)!=0) {func<- function(y_min,y_max) (1-cdf_min_max_partial(y_min,y_max,distribution,arg1,arg2)/f_order_max(y_max,distribution,arg1,arg2))}
    else {func<- function(y_min,y_max) 0}}

  return(func(y_min,y_max))

}

#' R function that calculates the probability of r-th order statistics of normal random variables (CDF of r-th order statistics)
#'
#' @param x the value of x
#' @param r r-th order statistics
#' @param distribution the distribution of marker value follows Normal or Gamma
#' @param arg1 if distribution is normal input mean parameters of all subclasses in a vector, if gamma input shape parameters
#' @param arg2 if distribution is normal input variance parameter, if gamma input rate parameters
#' @return The probability of r-th order statistics of random variables smaller or equal to x
#' @export
cdf_order_r <- function(x,r,distribution,arg1,arg2){
  k <- length(arg1)
  perm <- perms(c(1:k))

  if(distribution=="Normal"){
    mu=arg1
    sd=arg2

  index_k <- 1/factorial(k)
  if(r<k){
    F_r <- function(x){

      eval(parse(text=paste(paste(sapply(c(r:(k-1)), function(w) {

        index <- 1/(factorial(w)*factorial(k-w));

        paste(index,'*','(',
              paste(sapply(1:dim(perm)[1], function(j) paste(paste(sapply(perm[j,1:w],function(a) substitute(pnorm(x,mu[i],sd[i]),list(i=a))),collapse="*"),
                                                             paste(sapply(perm[j,(w+1):k],function(a) substitute(pnorm(x,mu[i],sd[i],lower.tail = F),list(i=a))),collapse="*"),sep='*')),collapse='+'),")",sep='')})

        ,collapse='+'),
        '+',paste(paste(sapply(perm[1,],function(a) substitute(pnorm(x,mu[i],sd[i]),list(i=a))),collapse="*"),sep=''))))
    }


  }
  if(r==k){
    F_r <- function(x){
      index_k <- 1/factorial(k)
      eval(parse(text=paste(paste(sapply(perm[1,],function(a) substitute(pnorm(x,mu[i],sd[i]),list(i=a))),collapse="*"),sep='')))
    }
  }
}

  if(distribution=="Gamma"){
    shape=arg1
    rate=arg2
    index_k <- 1/factorial(k)
    if(r<k){
      F_r <- function(x){

        eval(parse(text=paste(paste(sapply(c(r:(k-1)), function(w) {

          index <- 1/(factorial(w)*factorial(k-w));

          paste(index,'*','(',
                paste(sapply(1:dim(perm)[1], function(j) paste(paste(sapply(perm[j,1:w],function(a) substitute(pgamma(x,shape[i],rate[i]),list(i=a))),collapse="*"),
                                                               paste(sapply(perm[j,(w+1):k],function(a) substitute(pgamma(x,shape[i],rate[i],lower.tail = F),list(i=a))),collapse="*"),sep='*')),collapse='+'),")",sep='')})

          ,collapse='+'),
          '+',paste(paste(sapply(perm[1,],function(a) substitute(pgamma(x,shape[i],rate[i]),list(i=a))),collapse="*"),sep=''))))
      }


    }
    if(r==k){
      F_r <- function(x){
        index_k <- 1/factorial(k)
        eval(parse(text=paste(paste(sapply(perm[1,],function(a) substitute(pgamma(x,shape[i],rate[i]),list(i=a))),collapse="*"),sep='')))
      }
    }
  }
  return(F_r(x))
}
