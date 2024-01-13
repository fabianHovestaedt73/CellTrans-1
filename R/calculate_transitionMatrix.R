#' Calculation of the transition matrix.
#'
#' This function derives the transition matrix which contains the transition probabilitites between the distinct cell states.
#' @param M Matrix containing all experimental data: the initial experimental setup matrix and the experimental cell state distribution matrices.
#' @param  t This A vector containing all timepoints at which the cell state distributions have been experimentally measured. ---> realTime
#' @param  used_timepoints A vector with the timepoints utilized to estimate the transition matrix. ---> Tau

#' @keywords transition matrix, Matrix process
#' @export
#' 
calculate_transitionMatrix <- function(M,t,used_timepoints, tau) {
  n=ncol(M) #Rang of Matrix
  transitionMatrix=0
  countQOM=0
  
  for (i in 1:length(t)) {
    # first submatrix in M contains initial matsrix, calculate inverse
    invInitialMatrix=solve(M[1:n,]) #
    
    # After the integration of tau, a different matrix root than the original value must be calculated at a certain point in time, so this value is calculated here. For example, for a value tau=0.5, the 4th matrix root is now calculated instead of the 2nd matrix root, as 2/0.5=4.
    real_time <- t[i]
    k <- real_time/tau
    
    # derive transition matrices beginning with second submatrix in M
    if (t[i] %in% used_timepoints ) {
      
      if (t[i]>1) {
        Ptemp <- expm( (1/(k)) * logm( (invInitialMatrix%*%M[(i*n+1):((i+1)*n), ]) ,method="Eigen"))
      }
      else {    Ptemp=(invInitialMatrix%*%M[(i*n+1):((i+1)*n), ])%^%(1/t[i])}
      
      if (isTrMatrix(Ptemp)==FALSE) {
        assign( paste("P",t[i],sep=""),QOM(Ptemp)  )
        countQOM=countQOM+1
      }
      else {    assign( paste("P",t[i],sep=""),Ptemp)}
      
      transitionMatrix=transitionMatrix+get( paste("P",t[i],sep=""))
    }
  }
  
  transitionMatrix <- transitionMatrix/(length(used_timepoints))
  
  return(transitionMatrix)
}