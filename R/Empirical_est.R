
#' R function for obtaining all combinations of maximum and minimum of a given dataset
#'
#' @param df Given dataset, in list form
#' @return A list of all combinations of maximum and minimum of df
#' @export
get_max_min_permutations <- function(df) {
  n <- lengths(df)

  seq_list <- lapply(seq_along(n), function(i) seq_len(n[i]))
  result <- do.call(expand.grid, seq_list)

  apply(result,1, function(indices) {
    index <- unlist(indices)
    df0 <- sapply(1:length(index), function(i) df[[i]][index[i]],simplify=T)
    max_val <- max(df0)
    min_val <- min(df0)
    return(c(max = max_val, min = min_val))
  })
}


#' R function that calculates empirical estimates of HUMcm
#'
#' This function provides empirical estimates of HUMcm
#' @param dat test values in list, each element represents biomarker values for a disease group
#' @param num_sub a vector of number of subclasses in each subclass
#' @return The empirical estimate of HUMcm based on given data and num_sub
#' @examples
#' # Create a list of example data
#' Y1 <- c(0.9316, 0.9670, 1.3856, 1.3505, 1.0316, 1.1764, 0.7435, 0.5813, 0.4695, 0.3249)
#' Y2 <- c(1.63950, 1.36535, 1.79859, 0.47961, 1.50978, 1.36525,0.13515, 2.11275, 0.45659)
#' Y3 <- c(1.89856, 1.30920, 2.38615, 2.34785, 2.92493, 2.71615, 2.75243, 0.95060, 0.38964)
#' Y4 <- c(2.580,2.570,2.143,3.079,1.765,3.081,2.175,2.306,2.918,2.507,4.261,3.033,1.836,2.321)
#' Y5 <- c(3.969,3.044,3.318,2.862,3.655,1.523,3.722,4.074,3.662,3.571,5.177,6.321,4.932,4.129)

#' Y.dat <- list(Y1,Y2,Y3,Y4,Y5)

#' num_sub <- c(1,3,1)

#' # calculate HUMcm of Y.dat and num_sub

#' hum.dynamic(Y.dat,num_sub)

#' @export

hum.dynamic <- function(dat, num_sub) {
  # Calculate the maximum and minimum values for each subset
  subsets <- split(dat, rep(1:length(num_sub), times=num_sub))

  max_min_permutations <- lapply(subsets, function(sub) {
    list(
      max = get_max_min_permutations(sub)["max", ],
      min = get_max_min_permutations(sub)["min", ]
    )
  })

  # Initialize lists
  layer_in <- list()
  layer_out <- list()
  m <- c()

  for (i in 1:length(num_sub)) {
    layer_in[[i]] <- max_min_permutations[[i]]$max
    layer_out[[i]] <- max_min_permutations[[i]]$min
    m[i] <- length(layer_in[[i]])

  }

  attri_in <- rep(1,m[1])

  # Iterate over edges dynamically
  for (k in 1:(length(num_sub)-1)) {
    attri_out <- rep(0,m[k+1])
    temp <- c()
    for (i in 1:m[k+1]) {
      for (j in 1:m[k]) {
        temp[j] <- attri_in[j] * (layer_in[[k]][j] < layer_out[[k+1]][i])
      }
      attri_out[i] <- sum(temp)

    }
    if (k < length(num_sub) - 1) {
      attri_in <- attri_out
    }
  }
  # Final computation
  humC0 <- sum(attri_out) / prod(m)

  return(humC0)
}

#' R function that calculates percentile confidence interval given an array of estimates
#'
#' This function provides percentile confidence interval

#' @param x an array of calculated estimates
#' @return The percentile confidence interval of given values
#' @export

CI.func <- function(x){
  cil.1 <- unname(quantile(x,0.025,na.rm = TRUE,type=6))
  ciu.1<- unname(quantile(x,0.975,na.rm = TRUE,type=6))
  est <- mean(x,na.rm=TRUE)
  delta <- ciu.1-cil.1
  n <- sum(!is.na(x))
  return(c(n=n,est=est,cil.1=cil.1,ciu.1=ciu.1,width=delta))
}

#' R function that calculates non-parametric bootstrap percentile confidence interval
#'
#' This function provides non-parametric bootstrap percentile confidence interval of HUMcm

#' @param dat test values in list, each element represents biomarker values for a disease group
#' @param num_sub a vector of number of subclasses in each subclass
#' @param B the number of iteration
#' @return The non-parametric bootstrap percentile confidence interval of HUMcm
#' @examples

#' # Create a list of example data

#' Y1 <- c(0.9316, 0.9670, 1.3856, 1.3505, 1.0316, 1.1764, 0.7435, 0.5813, 0.4695, 0.3249)
#' Y2 <- c(1.63950, 1.36535, 1.79859, 0.47961, 1.50978, 1.36525,0.13515, 2.11275, 0.45659)
#' Y3 <- c(1.89856, 1.30920, 2.38615, 2.34785, 2.92493, 2.71615, 2.75243, 0.95060, 0.38964)
#' Y4 <- c(2.580,2.570,2.143,3.079,1.765,3.081,2.175,2.306,2.918,2.507,4.261,3.033,1.836,2.321)
#' Y5 <- c(3.969,3.044,3.318,2.862,3.655,1.523,3.722,4.074,3.662,3.571,5.177,6.321,4.932,4.129)

#' Y.dat <- list(Y1,Y2,Y3,Y4,Y5)
#' num_sub <- c(1,3,1)

#' # calculate the non-parametric bootstrap percentile confidence interval

#' HUMC_NPCI(Y.dat,num_sub,50)

#' @export

HUMC_NPCI <- function(dat,num_sub,B){

  humC.network <- c()

  repeat{

    databoot <- lapply(dat,function(x){sample(x,size = length(x),replace = TRUE)})
    humC.value <- hum.dynamic(dat=databoot,num_sub=num_sub)
    humC.network <- c(humC.network,humC.value)

    if(length(humC.network)==B) break
  }

  CI <- CI.func(humC.network)

  return(CI)
}
