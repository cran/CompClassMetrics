#' @import plot3D
#' @importFrom graphics lines
#' @importFrom stats dgamma dnorm integrate pgamma pnorm qgamma qnorm quantile uniroot
NULL
#' R function for plotting the compound ROC surface and chance surface
#'
#' @param k1 number of subclasses in main class-1
#' @param k2 number of subclasses in main class-2
#' @param k3 number of subclasses in main class-3
#' @param distribution the distribution of marker value follows Normal or Gamma
#' @param arg1 if distribution is normal input mean parameters of all subclasses in a vector, if gamma input shape parameters
#' @param arg2 if distribution is gamma input variance parameter, if gamma input rate parameters
#' @return The compound ROC surface and chance surface
#' @importFrom plot3D persp3D
#' @export
rocc_surface <- function(k1,k2,k3,distribution,arg1,arg2){
change_surface <- function(p1,p3){((1-p3^{1/k3}-p1^{1/k1})^{k2})*(1-p3^{1/k3}>p1^{1/k1})}
ROCC <- function(p1,p3){
    if(distribution=="Normal"){
      if(k1==1){
        cut1 <-function(p1){
          qnorm(p1,arg1[1],arg2[1])
        }}
      if(k1>1){
        solve1 <- function(p1,c1){p1-eval(parse(text=paste(sapply(1:k1,function(a) substitute(pnorm(c1,arg1[i],arg2[i]),list(i=a))),collapse="*")))}
        cut1 <-function(p1){
          uniroot(solve1,interval=c(-500,5000),p1=p1)$root
        }}
      Cut1 <- Vectorize(cut1)
      if(k3==1){
        cut2<-function(p3){
          qnorm(1-p3,arg1[k1+k2+1],arg2[k1+k2+1])
        }}
      if(k3>1){
        solve2 <- function(p3,c2){1-p3-eval(parse(text=paste(sapply(1:k3,function(a) substitute(pnorm(c2,arg1[k1+k2+i],arg2[k1+k2+i],lower.tail = F),list(i=a))),collapse="*")))}
        cut2 <-function(p3){
          uniroot(solve2,interval=c(-500,5000),p3=p3)$root
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
          qgamma(p1,shape=arg1[1],rate=arg2[1])
        }}
      if(k1>1){
        solve1 <- function(p1,c1){p1-eval(parse(text=paste(sapply(1:k1,function(a) substitute(pgamma(c1,shape=arg1[i],rate=arg2[i]),list(i=a))),collapse="*")))}
        cut1 <-function(p1){
          uniroot(solve1,interval=c(0,5000),p1=p1)$root
        }}
      Cut1 <- Vectorize(cut1)
      if(k3==1){
        cut2<-function(p3){
          qgamma(1-p3,shape=arg1[k1+k2+1],rate=arg2[k1+k2+1])
        }}
      if(k3>1){
        solve2 <- function(p3,c2){1-p3-eval(parse(text=paste(sapply(1:k3,function(a) substitute(pgamma(c2,shape=arg1[k1+k2+i],rate=arg2[k1+k2+i],lower.tail = F),list(i=a))),collapse="*")))}
        cut2 <-function(p3){
          uniroot(solve2,interval=c(0,5000),p3=p3)$root
        }}
      Cut2 <- Vectorize(cut2)

      if(k2==1){
        integrand.3q1 <- function(p1,p3){(pgamma(Cut2(p3),shape=arg1[k1+1],rate=arg2[k1+1])-pgamma(Cut1(p1),shape=arg1[k1+1],rate=arg2[k1+1]))*(Cut2(p3)>Cut1(p1))}
      }
      if(k2>1){
        integrand.3q1 <- function(p1,p3){eval(parse(text=paste(sapply(1:k2,function(a) substitute((pgamma(Cut2(p3),shape=arg1[k1+i],rate=arg2[k1+i])-pgamma(Cut1(p1),shape=arg1[k1+i],rate=arg2[k1+i])),list(i=a))),collapse="*")))*(Cut2(p3)>Cut1(p1))}
      }
    }

  value <- integrand.3q1(p1,p3)
  return(value)
}

# prepare variables.
p1 <- p3 <- seq(0, 1, length = 30)


TC2 <- outer(p1, p3, change_surface)
TC2.2 <- outer(p1, p3, ROCC)
# plot the 3D surface

persp3D(p1, p3, TC2,theta=290,plot=FALSE,box=T,colkey = F,cex.axis=1.5,cex.lab = 2.2,col='grey',border = "black",xlab='p1',ylab='p3',zlab='p2',scale=T)
persp3D(p1, p3, TC2.2,add=T,colkey = F,box=F,expand = 20, cex.axis=1.5,cex.lab = 2.2,col='pink',border = "black")

}

#' R function for plotting the overall ROC curve and chance curve
#'
#' @param k1 number of subclasses in main class-1
#' @param k2 number of subclasses in main class-2
#' @param distribution the distribution of marker value follows Normal or Gamma
#' @param arg1 if distribution is normal input mean parameters of all subclasses in a vector, if gamma input shape parameters
#' @param arg2 if distribution is gamma input variance parameter, if gamma input rate parameters
#' @return The overall ROC curve and chance curve
#' @export
rocc_curve <- function(k1,k2,distribution,arg1,arg2){

  change_curve <- function(p){(1-(1-p)^(1/k1))^k2}

  ROCC_binary <- function(p){

  if(distribution=="Normal"){
      if(k1==1){
        cut1 <-function(p){
          qnorm(p,arg1[1],arg2[1])
        }}
      if(k1>1){
        solve1 <- function(p,c1){p-eval(parse(text=paste(sapply(1:k1,function(a) substitute(pnorm(c1,arg1[i],arg2[i]),list(i=a))),collapse="*")))}
        cut1 <-function(p){
          uniroot(solve1,interval=c(-500,5000),p=p)$root
        }}

      Cut1 <- Vectorize(cut1)

      part1 <- paste0(
        "(",
        paste(
          sapply(1:k2, function(a) {
            deparse(substitute(pnorm(Cut1(p), arg1[k1 + i], arg2[k1 + i], lower.tail = F), list(i = a)))
          }),
          collapse = " * "
        ),
        ")"
      )

  }

    if(distribution=="Gamma"){
      if(k1==1){
        cut1 <-function(p){
          qgamma(p,arg1[1],arg2[1])
        }}
      if(k1>1){
        solve1 <- function(p,c1){p-eval(parse(text=paste(sapply(1:k1,function(a) substitute(pgamma(c1,arg1[i],arg2[i]),list(i=a))),collapse="*")))}
        cut1 <-function(p){
          uniroot(solve1,interval=c(0.001,5000),p=p)$root
        }}

      Cut1 <- Vectorize(cut1)

      part1 <- paste0(
        "(",
        paste(
          sapply(1:k2, function(a) {
            deparse(substitute(pgamma(Cut1(p), arg1[k1 + i], arg2[k1 + i], lower.tail = F), list(i = a)))
          }),
          collapse = " * "
        ),
        ")"
      )

    }
      final <- function(p){eval(parse(text = part1))}
      value <- final(1-p)
      return(value)
  }

  p <- seq(0, 1, length = 30)

  TC2 <- sapply(p, change_curve)
  TC2.2 <- sapply(p, ROCC_binary)
  # plot the 3D surface

  plot(p, TC2, type = "l", col = "grey", lwd = 2, ylim = c(0, 1),
       xlab = "1-Speo", ylab = "Seno", main = "")

  # Add the second density
  lines(p, TC2.2, col = "pink", lwd = 2)

}
