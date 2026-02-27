#' @import pracma
#' @import cubature
NULL
#' R function that calculates the true values of AUCo when distribution is known
#'
#' @param k1 number of subclasses in main class-1
#' @param k2 number of subclasses in main class-2
#' @param distribution the distribution of marker value follows Normal or Gamma
#' @param arg1 if distribution is normal input mean parameters of all subclasses in a vector, if gamma input arg1 parameters
#' @param arg2 if distribution is gamma input variance parameter, if gamma input arg2 parameters
#' @return The true value of AUCo under given distribution and parameters
#' @importFrom cubature hcubature
#' @export
auco_func <- function(k1,k2,distribution,arg1,arg2){

if(distribution=="Normal"){
  integrand <- function(s) {
    # Part 1: Genearg2 the first term (pnorm products for k2)
    part1 <- paste0(
      "(",
      paste(
        sapply(1:k2, function(a) {
          deparse(substitute(pnorm(s, arg1[k1 + i], arg2[k1 + i], lower.tail = F), list(i = a)))
        }),
        collapse = " * "
      ),
      ")"
    )

    # Part 2: Genearg2 the second term (sum of products for k1)

if(k1==1){
    part2 <- paste(dnorm(s,arg1[k1],arg2[k1]))
}
if(k1>1){
    part2 <- paste(
      sapply(1:k1, function(b) {
        # Products excluding the current index
        exclude_b <- paste(
          sapply(setdiff(1:k1, b), function(c) {
            deparse(substitute(pnorm(s, arg1[i], arg2[i]), list(i = c)))
          }),
          collapse = " * "
        )
        # Combine with dnorm for the current index
        paste0(exclude_b, " * ", deparse(substitute(dnorm(s, arg1[i], arg2[i]), list(i = b))))
      }),
      collapse = " + "
    )
}
    # Combine the two parts
    final_expr <- paste0(part1, " * (", part2, ")")

    # Return the parsed and evaluated function
    eval(parse(text = final_expr))
  }

  integral <- try(integrate(Vectorize(integrand), lower = -Inf, upper = Inf,subdivisions=100), silent = T)
}

##Gamma

if(distribution=="Gamma"){
    integrand <- function(s) {

      part1 <- paste0(
        "(",
        paste(
          sapply(1:k2, function(a) {
            deparse(substitute(pgamma(s, arg1[k1 + i], arg2[k1 + i], lower.tail = F), list(i = a)))
          }),
          collapse = " * "
        ),
        ")"
      )


if(k1==1){
  part2 <- paste(dgamma(s,arg1[k1],arg2[k2]))
}
if(k1>1){
      part2 <- paste(
        sapply(1:k1, function(b) {
          # Products excluding the current index
          exclude_b <- paste(
            sapply(setdiff(1:k1, b), function(c) {
              deparse(substitute(pgamma(s, arg1[i], arg2[i]), list(i = c)))
            }),
            collapse = " * "
          )
          # Combine with dnorm for the current index
          paste0(exclude_b, " * ", deparse(substitute(dgamma(s, arg1[i], arg2[i]), list(i = b))))
        }),
        collapse = " + "
      )
}
      # Combine the two parts
      final_expr <- paste0(part1, " * (", part2, ")")

      # Return the parsed and evaluated function
      eval(parse(text = final_expr))
    }

    integral <- try(integrate(Vectorize(integrand), lower = 0, upper = Inf,subdivisions=100), silent = T)
  }



  if(inherits(integral, "try-error")){
  AUCo = "NA"}else{
    AUCo= integral$value
  }
  return(AUCo)
}

#' R function that calculates the true values of  VUSC when distribution is known
#'
#' @param k1 number of subclasses in main class-1
#' @param k2 number of subclasses in main class-2
#' @param k3 number of subclasses in main class-3
#' @param distribution the distribution of marker value follows Normal or Gamma
#' @param arg1 if distribution is normal input mean parameters of all subclasses in a vector, if gamma input arg1 parameters
#' @param arg2 if distribution is gamma input variance parameter, if gamma input arg2 parameters
#' @return The true value of VUSc under given distribution and parameters
#' @export
cvus_func <- function(k1,k2,k3,distribution,arg1,arg2){
  ##Distribution: Normal
  if(distribution=="Normal"){
    if(k1==1){
      cut1 <-function(p1){
        qnorm(p1,arg1[1],arg2[1])
      }}
    if(k1>1){
      solve1 <- function(p1,c1){p1-eval(parse(text=paste(sapply(1:k1,function(a) substitute(pnorm(c1,arg1[i],arg2[i]),list(i=a))),collapse="*")))}
      cut1 <-function(p1){
        uniroot(solve1,interval=c(-50,500),p1=p1)$root
      }}
    Cut1 <- Vectorize(cut1)
    if(k3==1){
      cut2<-function(p3){
        qnorm(1-p3,arg1[k1+k2+1],arg2[k1+k2+1])
      }}
    if(k3>1){
      solve2 <- function(p3,c2){1-p3-eval(parse(text=paste(sapply(1:k3,function(a) substitute(pnorm(c2,arg1[k1+k2+i],arg2[k1+k2+i],lower.tail = F),list(i=a))),collapse="*")))}
      cut2 <-function(p3){
        uniroot(solve2,interval=c(-50,500),p3=p3)$root
      }}
    Cut2 <- Vectorize(cut2)

    if(k2==1){
      integrand.3q1 <- function(p1,p3){(pnorm(Cut2(p3),arg1[k1+1],arg2[k1+1])-pnorm(Cut1(p1),arg1[k1+1],arg2[k1+1]))*(Cut2(p3)>Cut1(p1))}
    }
    if(k2>1){
      integrand.3q1 <- function(p1,p3){eval(parse(text=paste(sapply(1:k2,function(a) substitute((pnorm(Cut2(p3),arg1[k1+i],arg2[k1+i])-pnorm(Cut1(p1),arg1[k1+i],arg2[k1+i])),list(i=a))),collapse="*")))*(Cut2(p3)>Cut1(p1))}
    }

  }
  ##Distribution: Gamma
  if(distribution=="Gamma"){
    if(k1==1){
      cut1 <-function(p1){
        qgamma(p1,arg1[1],arg2[1])
      }}
    if(k1>1){
      solve1 <- function(p1,c1){p1-eval(parse(text=paste(sapply(1:k1,function(a) substitute(pgamma(c1,arg1[i],arg2[i]),list(i=a))),collapse="*")))}
      cut1 <-function(p1){
        uniroot(solve1,interval=c(0.001,500),p1=p1)$root
      }}
    Cut1 <- Vectorize(cut1)
    if(k3==1){
      cut2<-function(p3){
        qgamma(1-p3,arg1[k1+k2+1],arg2[k1+k2+1])
      }}
    if(k3>1){
      solve2 <- function(p3,c2){1-p3-eval(parse(text=paste(sapply(1:k3,function(a) substitute(pgamma(c2,arg1[k1+k2+i],arg2[k1+k2+i],lower.tail = F),list(i=a))),collapse="*")))}
      cut2 <-function(p3){
        uniroot(solve2,interval=c(0.001,500),p3=p3)$root
      }}
    Cut2 <- Vectorize(cut2)

    if(k2==1){
      integrand.3q1 <- function(p1,p3){(pgamma(Cut2(p3),arg1[k1+1],arg2[k1+1])-pgamma(Cut1(p1),arg1[k1+1],arg2[k1+1]))*(Cut2(p3)>Cut1(p1))}
    }
    if(k2>1){
      integrand.3q1 <- function(p1,p3){eval(parse(text=paste(sapply(1:k2,function(a) substitute((pgamma(Cut2(p3),arg1[k1+i],arg2[k1+i])-pgamma(Cut1(p1),arg1[k1+i],arg2[k1+i])),list(i=a))),collapse="*")))*(Cut2(p3)>Cut1(p1))}
    }

  }
  integrand.3q <- function(p) as.matrix(integrand.3q1(p[1,],p[2,]))

  integral.3 <- try(hcubature(integrand.3q, lowerLimit = c(0,0),
                              upperLimit = c(1,1), vectorInterface = T, tol=1e-6), silent = T)

  if(inherits(integral.3, "try-error")){
    CVUS.cubintegrate = NA}else{
      CVUS.cubintegrate= integral.3$integral
    }
  return(CVUS.cubintegrate)
}


#' R function that calculates the true values of HUMcm when distribution is known
#'
#' @param distribution the distribution of marker value follows Normal or Gamma
#' @param arg1 if distribution is normal input mean parameters of all subclasses in a vector, if gamma input arg1 parameters
#' @param arg2 if distribution is gamma input variance parameter, if gamma input arg2 parameters
#' @param num_sub the vector of number of subclasses in each main class
#' @return The true value of HUMcm under given distribution and parameters
#' @export
humc_fourclass <- function(distribution,arg1,arg2,num_sub){


k1 = num_sub[1]
k2 = num_sub[2]
k3 = num_sub[3]
k4 = num_sub[4]

  if(distribution=="Normal"){

    integrand.full <- function(x1,x2,x3){
      (cdf_min_given_max_partial_upper(y_min=x1,y_max=x2,distribution="Normal",arg1[(k1+1):(k1+k2)],arg2[(k1+1):(k1+k2)]))*
        f_order_max(y_max=x1,distribution="Normal",arg1[1:k1],arg2[1:k1])*(cdf_min_given_max_partial_upper(y_min=x2,y_max=x3,distribution="Normal",arg1[(k1+k2+1):(k1+k2+k3)],arg2[(k1+k2+1):(k1+k2+k3)]))*
        f_order_max(y_max=x2,distribution="Normal",arg1[(k1+1):(k1+k2)],arg2[(k1+1):(k1+k2)])*(1-cdf_order_r(x3,distribution="Normal",arg1[(k1+k2+k3+1):(k1+k2+k3+k4)],arg2[(k1+k2+k3+1):(k1+k2+k3+k4)],r=1))*f_order_max(x3,distribution="Normal",arg1[(k1+k2+1):(k1+k2+k3)],arg2[(k1+k2+1):(k1+k2+k3)])
    }

    integral <- integrate(Vectorize(function(x3) {
      integrate(Vectorize(function(x2) {
        integrate(function(x1) {
          integrand.full(x1,x2,x3)
        },-Inf,x2,rel.tol = 1e-6)$value
      }),-Inf,x3,rel.tol = 1e-6)$value
    }),-Inf,Inf,rel.tol = 1e-6)$value

  }

  if(distribution=="Gamma"){

    integrand.full <- function(x1,x2,x3){
      (cdf_min_given_max_partial_upper(y_min=x1,y_max=x2,distribution="Gamma",arg1[(k1+1):(k1+k2)],arg2[(k1+1):(k1+k2)]))*f_order_max(y_max=x1,distribution="Gamma",arg1[1:k1],arg2[1:k1])*
        (cdf_min_given_max_partial_upper(y_min=x2,y_max=x3,distribution="Gamma",arg1[(k1+k2+1):(k1+k2+k3)],arg2[(k1+k2+1):(k1+k2+k3)]))*f_order_max(y_max=x2,distribution="Gamma",arg1[(k1+1):(k1+k2)],arg2[(k1+1):(k1+k2)])*
        (1-cdf_order_r(x3,distribution="Gamma",arg1[(k1+k2+k3+1):(k1+k2+k3+k4)],arg2[(k1+k2+k3+1):(k1+k2+k3+k4)],r=1))*f_order_max(x3,distribution="Gamma",arg1[(k1+k2+1):(k1+k2+k3)],arg2[(k1+k2+1):(k1+k2+k3)])
    }

    integral <- integrate(Vectorize(function(x3) {
      integrate(Vectorize(function(x2) {
        integrate(function(x1) {
          integrand.full(x1,x2,x3)
        },0,x2,rel.tol = 1e-6)$value
      }),0,x3,rel.tol = 1e-6)$value
    }),0,Inf,rel.tol = 1e-6)$value

  }
  return(round(integral,6))
}

#' R function that calculates the minimum of HUMcm under given structure
#'
#' @param num_sub the vector of number of subclasses in each main class
#' @return the minimum of HUMcm
#' @export
humc_min <- function(num_sub) {
  factorials <- factorial(num_sub)
  numerator <- prod(factorials)
  denominator <- factorial(sum(num_sub))
  HUM_min_value <- numerator / denominator
  return(HUM_min_value)

}

#' R function to calculate the standardized HUMcm under given structure
#' @param value the value of HUMcm
#' @param num_sub the vector of number of subclasses in each main class
#' @return The standardized HUMcm
#' @export
humc_standard <- function(value,num_sub){
  HUM_standard_value <- (value-humc_min(num_sub))/(1-humc_min(num_sub))
  return(HUM_standard_value)
}
