#' @import pracma
NULL
#' R function that calculates the probability density of minimum of NIND normal random variables (PDF)
#'
#' @param y_min the value of y_min
#' @param mu the vector of mean parameters of normal random variables
#' @param sd the vector of variance parameters of normal random variables
#' @return The probability density of minimum of normal random variables
#' @importFrom pracma perms
#' @export
f_order_min_normal <- function(y_min,mu,sd){
  k=length(mu)
  perm <- perms(c(1:k))
  if(k>1){f_min <- function(y_min){
    1/factorial(k-1)*eval(parse(text=paste(sapply(1:dim(perm)[1], function(j) paste(paste(sapply(perm[j,1],function(a) substitute(dnorm(y_min,mu[i],sd[i]),list(i=a)))),
                                                                                    paste(sapply(perm[j,2:k],function(a) substitute(pnorm(y_min,mu[i],sd[i],lower.tail = F),list(i=a))),collapse="*"),sep='*')),collapse='+')))

  }}
  if(k==1){f_min <- function(y_min) dnorm(y_min,mu,sd)}

  return(f_min(y_min))
}

#' R function that calculates the probability density of maximum of NIND normal random variables (PDF)
#'
#' @param y_max the value of y_max
#' @param mu the vector of mean parameters of normal random variables
#' @param sd the vector of variance parameters of normal random variables
#' @return The probability density of maximum of normal random variables
#' @export
f_order_max_normal <- function(y_max,mu,sd){
  k=length(mu)
  perm <- perms(c(1:k))
  if(k>1){
    f_max <- function(y_max){
      1/factorial(k-1)*eval(parse(text=paste(sapply(1:dim(perm)[1], function(j) paste(paste(sapply(perm[j,1],function(a) substitute(dnorm(y_max,mu[i],sd[i]),list(i=a)))),
                                                                                      paste(sapply(perm[j,2:k],function(a) substitute(pnorm(y_max,mu[i],sd[i]),list(i=a))),collapse="*"),sep='*')),collapse='+')))}
  }
  if(k==1){
    f_max <- function(y_max) dnorm(y_max,mu,sd)
  }
  return(f_max(y_max))
}

#' R function that calculates the partial of joint probability of min and max over max of NIND normal random variables
#'
#' @param y_max the value of y_max
#' @param y_min the value of y_min
#' @param mu the vector of mean parameters of normal random variables
#' @param sd the vector of variance parameters of normal random variables
#' @return The partial of joint probablity of min and max over max
#' @export

F_min_max_partial_normal <- function(y_min,y_max,mu,sd){
    k <- length(mu)
    perm <- perms(c(1:k))

    (y_min<=y_max)*(1/factorial(k-1)*


      eval(parse(text=paste(sapply(1:dim(perm)[1], function(j) paste('(',paste(paste(sapply(perm[j,1],function(a) substitute(dnorm(y_max,mu[i],sd[i]),list(i=a)))),
                                                                               paste(sapply(perm[j,2:k],function(a) substitute(pnorm(y_max,mu[i],sd[i]),list(i=a))),collapse="*"),sep='*'),'-',
                                                                     paste(paste(sapply(perm[j,1],function(a) substitute(dnorm(y_max,mu[i],sd[i]),list(i=a)))),
                                                                           paste(sapply(perm[j,2:k],function(a) substitute((pnorm(y_max,mu[i],sd[i])-pnorm(y_min,mu[i],sd[i])),list(i=a))),collapse="*"),sep='*'),')',sep='')),collapse='+')

    )))

}

#' R function that calculates the conditional probability of minimum greater than y_min given maximum equals to y_max of normal random variables (upper tail probability of minimum given maximum)
#'
#' @param y_min the value of y_min
#' @param y_max the value of y_max
#' @param mu the vector of mean parameters of normal random variables
#' @param sd the vector of variance parameters of normal random variables
#' @return The conditional probability of minimum given maximum of normal random variables
#' @export
F_min_given_max_partial_normal_upper <- function(y_min,y_max,mu,sd){
  k <- length(mu)
  perm <- perms(c(1:k))

  if(k==1){func <- function(y_min,y_max) 1}

  if(k>1){if(f_order_max_normal(y_max,mu,sd)!=0) {func<- function(y_min,y_max) (1-F_min_max_partial_normal(y_min,y_max,mu,sd)/f_order_max_normal(y_max,mu=mu,sd=sd))}
    else {func<- function(y_min,y_max) 0}}


  return(func(y_min,y_max))

}

#' R function that calculates the probability of r-th order statistics of normal random variables (CDF of r-th order statistics)
#'
#' @param x the value of x
#' @param r r-th order statistics
#' @param mu the vector of mean parameters of normal random variables
#' @param sd the vector of variance parameters of normal random variables
#' @return The probability of r-th order statistics of normal random variables smaller or equal to x
#' @export
F_order_r_normal <- function(x,mu,sd,r){
  k <- length(mu)
  perm <- perms(c(1:k))
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
  return(F_r(x))
}

###gamma

#' R function that calculates the probability density of minimum of gamma random variables (PDF)
#'
#' @param y_min the value of y_min
#' @param shape the vector of shape parameters of gamma random variables
#' @param rate the vector of rate parameters of gamma random variables
#' @return The probability density of minimum of gamma random variables
#' @export

f_order_min_gamma <- function(y_min,shape,rate){
  k=length(shape)
  perm <- perms(c(1:k))
  if(k>1){f_min <- function(y_min){
    1/factorial(k-1)*eval(parse(text=paste(sapply(1:dim(perm)[1], function(j) paste(paste(sapply(perm[j,1],function(a) substitute(dgamma(y_min,shape[i],rate[i]),list(i=a)))),
                                                                                    paste(sapply(perm[j,2:k],function(a) substitute(pgamma(y_min,shape[i],rate[i],lower.tail = F),list(i=a))),collapse="*"),sep='*')),collapse='+')))

  }}
  if(k==1){f_min <- function(y_min) dgamma(y_min,shape,rate)}

  return(f_min(y_min))
}

#' R function that calculates the probability density of maximum of gamma random variables (PDF)
#'
#' @param y_max the value of y_max
#' @param shape the vector of shape parameters of gamma random variables
#' @param rate the vector of rate parameters of gamma random variables
#' @return The probability density of maximum of gamma random variables
#' @export

f_order_max_gamma <- function(y_max,shape,rate){
  k=length(shape)
  perm <- perms(c(1:k))
  if(k>1){
    f_max <- function(y_max){
      1/factorial(k-1)*eval(parse(text=paste(sapply(1:dim(perm)[1], function(j) paste(paste(sapply(perm[j,1],function(a) substitute(dgamma(y_max,shape[i],rate[i]),list(i=a)))),
                                                                                      paste(sapply(perm[j,2:k],function(a) substitute(pgamma(y_max,shape[i],rate[i]),list(i=a))),collapse="*"),sep='*')),collapse='+')))}
  }
  if(k==1){
    f_max <- function(y_max) dgamma(y_max,shape,rate)
  }
  return(f_max(y_max))
}

#' R function that calculates the probability of r-th order statistics of gamma random variables (CDF of r-th order statistics)
#'
#' @param x the value of x
#' @param r r-th order statistics
#' @param shape the vector of shape parameters of gamma random variables
#' @param rate the vector of rate parameters of gamma random variables
#' @return The probability of r-th order statistics of gamma random variables smaller or equal to x
#' @export

F_order_r_gamma <- function(x,shape,rate,r){
  k <- length(shape)
  perm <- perms(c(1:k))
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
  return(F_r(x))
}




#' R function that calculates the partial of joint probability of min and max over max of NIND gamma random variables
#'
#' @param y_max the value of y_max
#' @param y_min the value of y_min
#' @param shape the vector of shape parameters of gamma random variables
#' @param rate the vector of rate parameters of gamma random variables
#' @return The partial of joint probablity of min and max over max
#' @export

F_min_max_partial_gamma <- function(y_min,y_max,shape,rate){
  k <- length(shape)
  perm <- perms(c(1:k))

  (y_min<=y_max)*1/factorial(k-1)*


    eval(parse(text=paste(sapply(1:dim(perm)[1], function(j) paste('(',paste(paste(sapply(perm[j,1],function(a) substitute(dgamma(y_max,shape[i],rate[i]),list(i=a)))),
                                                                             paste(sapply(perm[j,2:k],function(a) substitute(pgamma(y_max,shape[i],rate[i]),list(i=a))),collapse="*"),sep='*'),'-',
                                                                   paste(paste(sapply(perm[j,1],function(a) substitute(dgamma(y_max,shape[i],rate[i]),list(i=a)))),
                                                                         paste(sapply(perm[j,2:k],function(a) substitute((pgamma(y_max,shape[i],rate[i])-pgamma(y_min,shape[i],rate[i])),list(i=a))),collapse="*"),sep='*'),')',sep='')),collapse='+')

    ))

}

#' R function that calculates the conditional probability of minimum greater than y_min given maximum equals to y_max of gamma random variables (upper tail of conditional probability of minimum given maximum)
#'
#' @param y_min the value of y_min
#' @param y_max the value of y_max
#' @param shape the vector of shape parameters of gamma random variables
#' @param rate the vector of rate parameters of gamma random variables
#' @return The conditional probability of minimum given maximum of gamma random variables
#' @export

F_min_given_max_partial_gamma_upper <- function(y_min,y_max,shape,rate){
  k <- length(shape)
  perm <- perms(c(1:k))
  if(k==1){func <- function(y_min,y_max) 1}

  if(k>1){if(f_order_max_gamma(y_max,shape,rate)!=0){func<- function(y_min,y_max) (1-F_min_max_partial_gamma(y_min,y_max,shape,rate)/f_order_max_gamma(y_max,shape=shape,rate=rate))}
    if(f_order_max_gamma(y_max,shape,rate)==0){func<- function(y_min,y_max) 0}}

  return(func(y_min,y_max))
}
