#' @import plot3D
#' @importFrom graphics lines
#' @importFrom stats dgamma dnorm integrate pgamma pnorm qgamma qnorm quantile uniroot ecdf
NULL
#' R function for plotting the empirical compound ROC surface and chance surface
#'
#' @param dat values in list, each element represents biomarker values for a disease group, ordered in ascending severity
#' @param num_sub a vector of number of subclasses in each subclass
#' @return The empirical compound ROC surface and chance surface
#' @importFrom plot3D persp3D
#' @export


rocc_surface_emp <- function(dat,num_sub){
  k1=num_sub[1]
  k2=num_sub[2]
  k3=num_sub[3]

  chance_surface <- function(p1,p3){((1-p3^{1/k3}-p1^{1/k1})^{k2})*(1-p3^{1/k3}>p1^{1/k1})}

  ROCC <- function(p1,p3){
      if(k1==1){
        cut1 <-function(p1){
          quantile(dat[[1]], probs = p1)
        }}
      if(k1>1){
        solve1 <- function(p1,c1){p1-eval(parse(text=paste(sapply(1:k1,function(a) substitute(ecdf(dat[[i]])(c1),list(i=a))),collapse="*")))}
        cut1 <-function(p1){
          uniroot(solve1,interval=c(-50,500),p1=p1)$root
        }}
      Cut1 <- Vectorize(cut1)
      if(k3==1){
        cut2<-function(p3){
          quantile(dat[[k1+k2+1]],1-p3)
        }}
      if(k3>1){
        solve2 <- function(p3,c2){1-p3-eval(parse(text=paste(sapply(1:k3,function(a) substitute((1-ecdf(dat[[k1+k2+i]])(c2)),list(i=a))),collapse="*")))}
        cut2 <-function(p3){
          uniroot(solve2,interval=c(-50,500),p3=p3)$root
        }}
      Cut2 <- Vectorize(cut2)

      if(k2==1){
        integrand.3q1 <- function(p1,p3){(ecdf(dat[[k1+1]])(Cut2(p3))-ecdf(dat[[k1+1]])(Cut1(p1)))*(Cut2(p3)>Cut1(p1))}
      }
      if(k2>1){
        integrand.3q1 <- function(p1,p3){eval(parse(text=paste(sapply(1:k2,function(a) substitute((ecdf(dat[[k1+i]])(Cut2(p3))-ecdf(dat[[k1+i]])(Cut1(p1))),list(i=a))),collapse="*")))*(Cut2(p3)>Cut1(p1))}
      }

    value <- integrand.3q1(p1,p3)
    return(value)
  }

  # prepare variables.
  p1 <- p3 <- seq(0, 1, length = 30)


  TC2 <- outer(p1, p3, chance_surface)
  TC2.2 <- outer(p1, p3, ROCC)
  # plot the 3D surface

  persp3D(p1, p3, TC2,theta=290,plot=FALSE,box=T,colkey = F,cex.axis=1.5,cex.lab = 2.2,col='grey',border = "black",xlab='p1',ylab='p3',zlab='p2',scale=T)
  persp3D(p1, p3, TC2.2,add=T,colkey = F,box=F,expand = 20, cex.axis=1.5,cex.lab = 2.2,col='pink',border = "black")

}
#' R function for plotting the empirical compound ROC curve and chance curve
#'
#' @param dat values in list, each element represents biomarker values for a disease group, ordered in ascending severity
#' @param num_sub a vector of number of subclasses in each subclass
#' @return The empirical compound ROC curve and chance curve
#' @importFrom plot3D persp3D
#' @export

rocc_curve_emp <- function(dat,num_sub){

  k1=num_sub[1]
  k2=num_sub[2]

  chance_curve <- function(p){(1-(1-p)^(1/k1))^k2}

  ROCC_binary <- function(p){
      if(k1==1){
        cut1 <-function(p){
          quantile(dat[[1]],p)
        }}
      if(k1>1){
        solve1 <- function(p,c1){p-eval(parse(text=paste(sapply(1:k1,function(a) substitute(ecdf(dat[[i]])(c1),list(i=a))),collapse="*")))}

        cut1 <-function(p){
          uniroot(solve1,interval=c(-5000,5000),p=p)$root
        }}

      Cut1 <- Vectorize(cut1)

      part1 <- paste0(
        "(",
        paste(
          sapply(1:k2, function(a) {
            deparse(substitute((1-ecdf(dat[[k1+i]])(Cut1(p))), list(i = a)))
          }),
          collapse = " * "
        ),
        ")"
      )

    final <- function(p){eval(parse(text = part1))}
    value <- final(1-p)
    return(value)
  }

  p <- seq(0, 1, length = 30)

  TC2 <- sapply(p, chance_curve)
  TC2.2 <- sapply(p, ROCC_binary)
  # plot the 3D surface

  plot(p, TC2, type = "l", col = "grey", lwd = 2, ylim = c(0, 1),
       xlab = "1-Speo", ylab = "Seno", main = "")

  # Add the second density
  lines(p, TC2.2, col = "pink", lwd = 2)

}
