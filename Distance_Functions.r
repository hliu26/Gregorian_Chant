
#Distance Metric: MSE, vectorized calculations

mseMatrix <- function(data){
    columnNames <- colnames(data)
    num_loc <- length(columnNames)
    sumSquaresMatrix = matrix(nrow = num_loc - 1, ncol = num_loc - 1)
    sumSquares = 0

    
    for (i in 2:num_loc){
        for (j in 2:num_loc){
            if (i == j){
                sumSquares = 0   
            } else{
                diffVector <- data[i] - data[j]
                diffVector[is.na(diffVector)] <- 1
                diffVector <- diffVector ** 2
                sumSquares <- sum(diffVector) / length(diffVector) 
            }
            sumSquaresMatrix[i-1, j-1] <- sumSquares
        }
    }

#    for (i in 1: num_loc - 1) {
#        for (j in 1: num_loc - 1) {
#            sumSquaresMatrix[i,j] <- (sumSquaresMatrix[i,j] - min(sumSquaresMatrix)) / 
#                                    (max(sumSquaresMatrix) - min(sumSquaresMatrix))
#        }
#    }


    colnames(sumSquaresMatrix) <- columnNames[2:num_loc]
    rownames(sumSquaresMatrix) <- columnNames[2:num_loc]
    sumSquaresMatrix

    #heatmap(sumSquaresMatrix, symm = TRUE)
}



#Distance Metric: Euclidean

euclidMatrix <- function(data){
    columnNames <- colnames(data)
    num_loc <- length(columnNames)
    euclidMatrix = matrix(nrow = num_loc - 1, ncol = num_loc - 1)
    sumSquare = 0


    for (j in 2:num_loc){
        for(k in 2:num_loc){
            for(i in 1:121){
                if (j == k){
                    sumSquare = 0
                } else if (is.na(data[,j][i]) | is.na(data[,k][i])){
                    sumSquare <- sumSquare + 1
                } else{
                    diff <- data[,j][i] - data[,k][i]
                    diffSquare <- diff ^ 2
                    sumSquare <- sumSquare + diffSquare
                }
            }
            sumSquare <- sumSquare ^1/2
            euclidMatrix[j-1, k-1] <- sumSquare
            sumSquare = 0
        }
    }
    

    colnames(euclidMatrix) <- columnNames[2:num_loc]
    rownames(euclidMatrix) <- columnNames[2:num_loc]
    euclidMatrix

    #heatmap(euclidMatrix, symm = TRUE)
}



#Distance Metric: Pearson's Distance

pearsonMatrix <- function(data){
    columnNames <- colnames(data)
    num_loc <- length(columnNames)
    corrMatrix = matrix(nrow = num_loc - 1, ncol = num_loc - 1)
    correlation = 0

    for (j in 2:num_loc){
        for(k in 2:num_loc){
            if (j == k){
                correlation = 1
            } else{
                correlation <- cor(data[j], data[k], use = "pairwise.complete.obs", method = "pearson")
            }
            corrMatrix[j-1, k-1] <- correlation
            correlation = 0 
        }
    }
        



    colnames(corrMatrix) <- columnNames[2:num_loc]
    rownames(corrMatrix) <- columnNames[2:num_loc]
    corrMatrix

    #heatmap(corrMatrix, symm = TRUE)
}



int_to_alpha <- function(data){
    columnNames <- colnames(data)
    num_loc <- length(columnNames)

    convert <- function(x){ 
        if (x == 0){
            x <- "a"
        } else if (x == -1){
            x <- "b"
        } else if (x == -3){
            x <- "d"
        } else if (x== -5){
            x <- "e"
        } else if (x == -7){
            x <- "f"
        } else if (x == -8){
            x <- "g"
        } else if (x == -10){
            x <- "h"
        } else if (x == 2){
            x <- "j"
        } else if (x == -2){
            x <- "c"
        } else if (x == -12){
            x <- "i"
        } else if (x == 5){
            x <- "k"
        } else{
             x <- "|else|"
        }
    }
    
    #remove NAs
    data_list <- c()
    for(i in 1:num_loc){
        data_list[i] <- na.omit(data[i])
    }
    data_list <- data_list[-1]

    #change integer to strings, then collapse
    for(i in 1:length(data_list)){
        data_list[[i]] <- sapply(data_list[[i]], convert)
        data_list[[i]] <- paste(data_list[[i]], collapse="") 
    }
    
    data_list
}

adistMatrix <- function(data){
    columnNames <- colnames(data)
    num_loc <- length(columnNames)
    trans_data <- int_to_alpha(data)
    adistMatrix = matrix(nrow=num_loc - 1, ncol=num_loc - 1)
    adistVector = c()
    adist = 0

    for (i in 1:num_loc - 1){
        for (j in 1:num_loc - 1){
            if (i == j){
                ad = 0
            } else {
                ad <-adist(trans_data[i], trans_data[j])
            }
            adistMatrix[i,j] <- ad
        }
    }

#    for (i in 1:num_loc - 1){
#        for (j in 1:num_loc - 1){
#            adistMatrix[i,j] <- (adistMatrix[i,j] - min(adistMatrix)) /(max(adistMatrix) - min(adistMatrix))
#        }
#    }
    colnames(adistMatrix) <- columnNames[2:num_loc]
    rownames(adistMatrix) <- columnNames[2:num_loc]
    adistMatrix
    

    #heatmap(adistMatrix, symm = TRUE)

}



#Distance Metric: Manhattan Distance

manMatrix <- function(data){
    columnNames <- colnames(data)
    num_loc <- length(columnNames)
    diffMatrix = matrix(nrow = num_loc-1, ncol = num_loc-1)
    diff = 0
    
    for (i in 2:num_loc){
        for (j in 2:num_loc){
            if (i == j){
                diff = 0   
            } else{
                diffVector <- abs(data[i] - data[j])
                diffVector[is.na(diffVector)] <- 1
                diff <- sum(diffVector)
            }
            diffMatrix[i-1, j-1] <- diff
        }
    }

        
#    for (i in 1: num_loc - 1) {
#        for (j in 1: num_loc - 1) {
#            diffMatrix[i,j] <- (diffMatrix[i,j] - min(diffMatrix)) /
#                                (max(diffMatrix) - min(diffMatrix))
#        }
#    }
   
    colnames(diffMatrix) <- columnNames[2:num_loc]
    rownames(diffMatrix) <- columnNames[2:num_loc]
    diffMatrix

    #heatmap(sumSquaresMatrix, symm = TRUE)
}



#Distance Metric: Canberra Distance (weighted Manhattan)

canMatrix <- function(data){
    columnNames <- colnames(data)
    num_loc <- length(columnNames)
    diffMatrix = matrix(nrow = num_loc - 1, ncol = num_loc - 1)
    diff = 0
    
    for (i in 2 : num_loc){
        for (j in 2 : num_loc){
            if (i == j){
                diff = 0   
            } else{
                diffVector <- (abs(data[i] - data[j]) / (abs(data[i]) + (abs(data[j]))))
                diffVector[is.na(diffVector)] <- 1
                diff <- sum(diffVector)
            }
            diffMatrix[i-1, j-1] <- diff
        }
    }


#    for (i in 1 : num_loc - 1){
#        for (j in 1 : num_loc - 1){
#            diffMatrix[i, j] <- (diffMatrix[i, j] - min(diffMatrix)) /
#                                         (max(diffMatrix) - min(diffMatrix))
#        }
#    }
    colnames(diffMatrix) <- columnNames[2:num_loc]
    rownames(diffMatrix) <- columnNames[2:num_loc]
    diffMatrix

    #heatmap(sumSquaresMatrix, symm = TRUE)
}


