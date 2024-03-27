




#' Title Normal distribution of base seed water potential
#' @description The distribution of base seed water potential is following the normal distribution. Here, the stress tolerance parameter is characterised by 50 percentile in the distribution of base seed water potential.
#' @param psi Soil water potential
#' @param time Time taken to germinate
#' @param c.germinated Cumulative number of seeds germinated
#' @param N Total number of seeds under each soil water potential
#' @param d Proportion of viability
#' @param theta Hydrotime constant

#' @import dplyr
#' @return
#' \itemize{
#'   \item parameters: mu (stress tolerance parameter)and sigma (uniformity of germination parameter)
#'   \item Result_fitting: Actual cumulative seed germination fraction (Actual_CGfraction) and Fitted cumulative seed germination fraction (Fitted_CGfraction)
#' }

#' @export
#' @usage seed.probit(psi, time, c.germinated, N, d=NULL, theta)
#' @examples
#' psi<-c(rep(0, 19), rep(-0.2, 6), rep(-0.4, 10))
#' time<- c(1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 9, 10, 12, 16, 18, 20, 23, 26,
#' 30, 4.5, 5, 6, 20, 23, 30, 3, 3.5, 4, 4.5, 5, 6, 7, 9, 12, 16)
#' c.germinated<- c(1, 2, 6, 11, 20, 24, 30, 34, 39, 41, 43, 47, 56, 58, 59,
#' 63, 67, 72, 73, 29, 31, 35, 63, 64, 65, 11, 13, 18, 21, 22, 25, 26, 28, 29, 30)
#' d<- c(0.8, 0.8, 0.6)
#' my.probit<-seed.probit(psi= psi, time= time, c.germinated= c.germinated, N=100, d=d, theta= 90)

#' @references
#' \itemize{
#' \item Bradford, K. J. (2002). Applications of Hydrothermal Time to Quantifying and Modeling Seed Germination and Dormancy. Weed Science, 50(2), 248–260. http://www.jstor.org/stable/4046371.

#' \item Bradford, K. J., & Still, D. W. (2004). Applications of Hydrotime Analysis in Seed Testing. Seed Technology, 26(1), 75–85. http://www.jstor.org/stable/23433495

#' }



seed.probit<-function(psi, time, c.germinated, N, d=NULL, theta)
{
  # Actual values to be compared
  actual.cgfraction<-c.germinated/N

  data<-as.data.frame(cbind(psi, time, c.germinated))
  unifirst<- unique(psi)



  for (i in 1:nrow(data)) {
    # Find the index of the unique value in unifirst that matches the value in First_Column
    index <- match(psi[i], unifirst)
    if (!is.null(d)) {
      # Check if the index is found
      if (!is.na(index)) {
        # Apply the dividing factor based on the index


        c.germinated[i] <- c.germinated[i] / (N*d[index])
      }
    }
    else {c.germinated <- actual.cgfraction}

  }



  n.cgfraction<-c.germinated

  n.data<-cbind(psi, time, actual.cgfraction, n.cgfraction)

  probit<-stats::qnorm(n.cgfraction)
  theta<- theta



  psi_b_g<- psi-(theta/(time*24))

  pro.reg<- stats::lm(probit~psi_b_g) # probit regression
  #summary(pro.reg)
  a<-as.numeric(pro.reg$coefficients[1])
  b<-as.numeric(pro.reg$coefficients[2])

  mu<- -(a/b)  # -a/b
  sigma<- 1/b  # 1/b
  x<- (psi-mu-(theta/(time*24)))/sigma

  k<-stats::pnorm(x, 0, 1) # fitted value

  n.data<-cbind(n.data, k)
  n.data <- as.data.frame(n.data)

  for (i in 1:nrow(n.data)) {
    # Find the index of the unique value in unifirst that matches the value in First_Column
    index <- match(n.data$psi[i], unifirst)
    if (!is.null(d)) {
      # Check if the index is found
      if (!is.na(index)) {
        # Apply the dividing factor based on the index
        n.data$k[i] <- n.data$k[i] * d[index]
      }
    }
    else {n.data$k <- n.data$k}
  }


  result.matrix<- matrix(nrow=1, ncol=2)

  colnames(result.matrix)<-c("mu", "sigma")
  result.matrix[1,1]<-mu
  result.matrix[1,2]<-sigma


  # Fitted values to be compared
  fitted.cgfraction<-n.data$k
  result.table<-as.data.frame(cbind(psi=psi, time= time,
                                    Actual_CGfraction=actual.cgfraction,
                                    Fitted_CGfraction=fitted.cgfraction))


  output<-list(parameters=result.matrix, Result_fitting=result.table)
  return(output)

}








#' Title Logistic distribution of base seed water potential
#' @description The distribution of base seed water potential is following the logistic distribution. Here, the stress tolerance parameter is characterised by 50 percentile in the distribution of base seed water potential.
#' @param psi Soil water potential
#' @param time Time taken to germinate
#' @param c.germinated Cumulative number of seeds germinated
#' @param N Total number of seeds under each soil water potential
#' @param d Proportion of viability
#' @param theta Hydrotime constant

#' @import dplyr
#' @return
#' \itemize{
#'   \item parameters: mu (stress tolerance parameter)and sigma (uniformity of germination parameter)
#'   \item Result_fitting: Actual cumulative seed germination fraction (Actual_CGfraction) and Fitted cumulative seed germination fraction (Fitted_CGfraction)
#' }

#' @export
#' @usage seed.logit(psi, time, c.germinated, N, d=NULL, theta)
#' @examples
#' psi<-c(rep(0, 19), rep(-0.2, 6), rep(-0.4, 10))
#' time<- c(1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 9, 10, 12, 16, 18, 20, 23, 26,
#' 30, 4.5, 5, 6, 20, 23, 30, 3, 3.5, 4, 4.5, 5, 6, 7, 9, 12, 16)
#' c.germinated<- c(1, 2, 6, 11, 20, 24, 30, 34, 39, 41, 43, 47, 56, 58, 59,
#' 63, 67, 72, 73, 29, 31, 35, 63, 64, 65, 11, 13, 18, 21, 22, 25, 26, 28, 29, 30)
#' d<- c(0.8, 0.8, 0.6)
#' my.logit<-seed.logit(psi= psi, time= time, c.germinated= c.germinated, N=100, d=d, theta= 90)

#' @references
#' \itemize{
#' \item Bradford, K. J. (2002). Applications of Hydrothermal Time to Quantifying and Modeling Seed Germination and Dormancy. Weed Science, 50(2), 248–260. http://www.jstor.org/stable/4046371.

#' \item Bradford, K. J., & Still, D. W. (2004). Applications of Hydrotime Analysis in Seed Testing. Seed Technology, 26(1), 75–85. http://www.jstor.org/stable/23433495

#' }



seed.logit<-function(psi, time, c.germinated, N, d=NULL, theta)
{
  # Actual values to be compared
  actual.cgfraction<-c.germinated/N

  data<-as.data.frame(cbind(psi, time, c.germinated))
  unifirst<- unique(psi)



  for (i in 1:nrow(data)) {
    # Find the index of the unique value in unifirst that matches the value in First_Column
    index <- match(psi[i], unifirst)
    if (!is.null(d)) {
      # Check if the index is found
      if (!is.na(index)) {
        # Apply the dividing factor based on the index


        c.germinated[i] <- c.germinated[i] / (N*d[index])
      }
    }
    else {c.germinated <- actual.cgfraction}

  }



  n.cgfraction<-c.germinated

  n.data<-cbind(psi, time, actual.cgfraction, n.cgfraction)

  logit<-log(n.cgfraction/(1-n.cgfraction))
  theta<- theta



  psi_b_g<- psi-(theta/(time*24))


  pro.reg<- stats::lm(logit~psi_b_g) # logit regression
  #summary(pro.reg)
  a<-as.numeric(pro.reg$coefficients[1])
  b<-as.numeric(pro.reg$coefficients[2])

  mu<- -(a/b) # -a/b
  sigma<- 1/b  # 1/b
  x<- (psi-mu-(theta/(time*24)))/sigma

  k<-1/(1+exp(-x)) # fitted value

  n.data<-cbind(n.data, k)
  n.data <- as.data.frame(n.data)

  for (i in 1:nrow(n.data)) {
    # Find the index of the unique value in unifirst that matches the value in First_Column
    index <- match(n.data$psi[i], unifirst)
    if (!is.null(d)) {
      # Check if the index is found
      if (!is.na(index)) {
        # Apply the dividing factor based on the index
        n.data$k[i] <- n.data$k[i] * d[index]
      }
    }
    else {n.data$k <- n.data$k}
  }


  result.matrix<- matrix(nrow=1, ncol=2)

  colnames(result.matrix)<-c("mu", "sigma")
  result.matrix[1,1]<-mu
  result.matrix[1,2]<-sigma


  # Fitted values to be compared
  fitted.cgfraction<-n.data$k
  result.table<-as.data.frame(cbind(psi=psi, time= time,
                                    Actual_CGfraction=actual.cgfraction,
                                    Fitted_CGfraction=fitted.cgfraction))


  output<-list(parameters=result.matrix, Result_fitting=result.table)
  return(output)

}









#' Title Extreme value distribution of base seed water potential
#' @description The distribution of base seed water potential is following the extreme value distribution. Here, the stress tolerance parameter is characterised by 63 percentile in the distribution of base seed water potential.
#' @param psi Soil water potential
#' @param time Time taken to germinate
#' @param c.germinated Cumulative number of seeds germinated
#' @param N Total number of seeds under each soil water potential
#' @param d Proportion of viability
#' @param theta Hydrotime constant

#' @import dplyr
#' @return
#' \itemize{
#'   \item parameters: mu (stress tolerance parameter)and sigma (uniformity of germination parameter)
#'   \item Result_fitting: Actual cumulative seed germination fraction (Actual_CGfraction) and Fitted cumulative seed germination fraction (Fitted_CGfraction)
#' }

#' @export
#' @usage seed.extreme(psi, time, c.germinated, N, d=NULL, theta)
#' @examples
#' psi<-c(rep(0, 19), rep(-0.2, 6), rep(-0.4, 10))
#' time<- c(1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 9, 10, 12, 16, 18, 20, 23, 26,
#' 30, 4.5, 5, 6, 20, 23, 30, 3, 3.5, 4, 4.5, 5, 6, 7, 9, 12, 16)
#' c.germinated<- c(1, 2, 6, 11, 20, 24, 30, 34, 39, 41, 43, 47, 56, 58, 59, 63,
#' 67, 72, 73, 29, 31, 35, 63, 64, 65, 11, 13, 18, 21, 22, 25, 26, 28, 29, 30)
#' d<- c(0.8, 0.8, 0.6)
#' my.extreme<-seed.extreme(psi= psi, time= time, c.germinated= c.germinated, N=100, d=d, theta= 90)

#' @references
#' \itemize{
#' \item Bradford, K. J. (2002). Applications of Hydrothermal Time to Quantifying and Modeling Seed Germination and Dormancy. Weed Science, 50(2), 248–260. http://www.jstor.org/stable/4046371.

#' \item Bradford, K. J., & Still, D. W. (2004). Applications of Hydrotime Analysis in Seed Testing. Seed Technology, 26(1), 75–85. http://www.jstor.org/stable/23433495

#' }






seed.extreme<-function(psi, time, c.germinated, N, d=NULL, theta)
{
  # Actual values to be compared
  actual.cgfraction<-c.germinated/N

  data<-as.data.frame(cbind(psi, time, c.germinated))
  unifirst<- unique(psi)



  for (i in 1:nrow(data)) {
    # Find the index of the unique value in unifirst that matches the value in First_Column
    index <- match(psi[i], unifirst)
    if (!is.null(d)) {
      # Check if the index is found
      if (!is.na(index)) {
        # Apply the dividing factor based on the index


        c.germinated[i] <- c.germinated[i] / (N*d[index])
      }
    }
    else {c.germinated <- actual.cgfraction}

  }



  n.cgfraction<-c.germinated

  n.data<-cbind(psi, time, actual.cgfraction, n.cgfraction)

  extreme<-log(-log(1 - n.cgfraction))
  theta<- theta



  psi_b_g<- psi-(theta/(time*24))


  pro.reg<- stats::lm(extreme~psi_b_g) # extreme value regression
  #summary(pro.reg)
  a<-as.numeric(pro.reg$coefficients[1])
  b<-as.numeric(pro.reg$coefficients[2])

  mu<- -(a/b)  # -a/b
  sigma<- 1/b  # 1/b
  x<- (psi-mu-(theta/(time*24)))/sigma

  k<- 1-exp(-exp(x)) # fitted value

  n.data<-cbind(n.data, k)
  n.data <- as.data.frame(n.data)

  for (i in 1:nrow(n.data)) {
    # Find the index of the unique value in unifirst that matches the value in First_Column
    index <- match(n.data$psi[i], unifirst)
    if (!is.null(d)) {
      # Check if the index is found
      if (!is.na(index)) {
        # Apply the dividing factor based on the index
        n.data$k[i] <- n.data$k[i] * d[index]
      }
    }
    else {n.data$k <- n.data$k}
  }


  result.matrix<- matrix(nrow=1, ncol=2)

  colnames(result.matrix)<-c("mu", "sigma")
  result.matrix[1,1]<-mu
  result.matrix[1,2]<-sigma


  # Fitted values to be compared
  fitted.cgfraction<-n.data$k
  result.table<-as.data.frame(cbind(psi=psi, time= time,
                                    Actual_CGfraction=actual.cgfraction,
                                    Fitted_CGfraction=fitted.cgfraction))


  output<-list(parameters=result.matrix, Result_fitting=result.table)
  return(output)

}


