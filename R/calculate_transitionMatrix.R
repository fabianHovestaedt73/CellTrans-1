#' Calculation of the transition matrix.
#'
#' This function derives the transition matrix which contains the transition probabilitites between the distinct cell states.
#' @param M Matrix containing all experimental data: the initial experimental setup matrix and the experimental cell state distribution matrices.
#' @param  t This A vector containing all timepoints at which the cell state distributions have been experimentally measured. ---> realTime
#' @param  used_timepoints A vector with the timepoints utilized to estimate the transition matrix. ---> Tau

#' @keywords transition matrix, Matrix process
#' @export

  calculate_transitionMatrix <- function(M,t,used_timepoints) {
    
  # step_size <- 0.1  # Kleiner Zeitschritt
  # # Generiere zusätzliche Zeitpunkte mit kleinerem Zeitschritt
  # additional_timepoints <- seq(from = min(t), to = max(t), by = step_size)
  # # Kombiniere vorhandene und zusätzliche Zeitpunkte
  # used_timepoints <- sort(unique(c(t, additional_timepoints)))
  n=ncol(M) #Rang of Matrix
  transitionMatrix=0
  countQOM=0
  
  for (i in 1:length(t)) {
    browser()
    #first submatrix in M contains initial matrix, calculate inverse
    invInitialMatrix=solve(M[1:n,]) #
    #derive transition matrices beginning with second submatrix in M
  if (t[i] %in% used_timepoints ) {
    if (t[i]>1) {
      
    #1. Der Ausdruck "invInitialMatrix %% M[(in+1):((i+1)n), ]" multipliziert die inverse Anfangsmatrix (invInitialMatrix) mit einer Teilmatrix von M. Die Teilmatrix wird aus M ausgewählt, beginnend bei der Position in+1 und endend bei (i+1)*n. Dies ergibt eine Teilmatrix der Größe n x n, die den Übergang von einem Zustand zum nächsten beschreibt.
    #2. Der Ausdruck "logm((invInitialMatrix %% M[(in+1):((i+1)*n), ]), method = "Eigen")" berechnet den Logarithmus der Teilmatrix. Hier wird die Funktion logm() verwendet, um den Logarithmus einer Matrix zu berechnen. Der Parameter "method = "Eigen"" gibt an, dass die Eigenwertmethode zur Berechnung verwendet werden soll.
    #3. Der Ausdruck "(1/t[i]) * logm((invInitialMatrix %% M[(in+1):((i+1)*n), ]), method = "Eigen")" teilt die berechnete logarithmische Matrix durch die Zeit t[i].
    #4. Der Ausdruck "expm((1/t[i]) * logm((invInitialMatrix %% M[(in+1):((i+1)*n), ]), method = "Eigen"))" wendet die Funktion expm() auf die berechnete Matrix an. Die Funktion expm() berechnet die Matrix-Exponentialfunktion.
    #5. Schließlich wird das Ergebnis der berechneten Übergangsmatrix Ptemp zugewiesen.
    #zusammengefasst: berechnet die Zeile Code die Übergangsmatrix Ptemp für den Zeitpunkt t[i], indem sie die inverse Anfangsmatrix mit einer Teilmatrix von M multipliziert, den Logarithmus dieser Teilmatrix berechnet, ihn durch die Zeit t[i] teilt und dann die Matrix-Exponentialfunktion anwendet.
    
    Ptemp=expm ( (1/t[i])*logm( (invInitialMatrix%*%M[(i*n+1):((i+1)*n), ]),method="Eigen"))
    } else {
    Ptemp=(invInitialMatrix%*%M[(i*n+1):((i+1)*n), ])%^%(1/t[i])
    }

    if (isTrMatrix(Ptemp)==FALSE) {
      assign( paste("P",t[i],sep=""),QOM(Ptemp)  )
      countQOM=countQOM+1
    }
    else {    assign( paste("P",t[i],sep=""),Ptemp)}

    transitionMatrix=transitionMatrix+get( paste("P",t[i],sep=""))
  }
  }
  transitionMatrix=transitionMatrix/length(used_timepoints)

  return(transitionMatrix)

}
