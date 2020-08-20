
source("Distance_Functions.r")

dist_funcs = c(mseMatrix, euclidMatrix, pearsonMatrix, adistMatrix, manMatrix, canMatrix)
dist_funcs_name = c("mse", "euclid", "pearson", "adist", "manhattan", "canberra")

for (i in dist_funcs_name) {
    dir.create(paste(i, "_distance", sep=""))
}



for (chants in dir("./Chants")){
    if (!(chants == 'README.md')){ 
        dist = 1
    
        #print(chants)
        curr_chant = read.csv(file=paste("./Chants/",chants, sep="" ), header=TRUE, sep=',')
    
        for (func in dist_funcs){
        
            write.table(func(curr_chant), file=paste(paste(dist_funcs_name[dist], "_distance", sep=""), 
                                                 paste("/", dist_funcs_name[dist], "_", chants, sep=""), sep=""))
            dist = dist + 1
        }
        dist = 1
    }
}

normalize_frame <- function(directory){
    folder = dir(paste("./", directory, "_distance", sep=""), pattern = ".csv")
    largest = -1000
    smallest = 1000
    
    frame_list = list()
    
    for (i in folder){
        currframe = read.table(paste("./", directory, "_distance/", i, sep=""), header=TRUE)
            
        if (largest < max(currframe)){
            largest = max(currframe)
        }
        if (smallest > min(currframe)){
            smallest = min(currframe)
        }
    }
    for (j in folder){
        currframe = read.table(paste("./", directory, "_distance/", j, sep=""), header=TRUE)
        currframe = (currframe - smallest) / (largest - smallest)
        write.table(currframe, file=paste(paste(directory, "_distance", sep=""), paste("/", j, sep=""), sep=""))
    }
    largest = -1000
    smallest = 1000
}

for (i in dist_funcs_name){
    normalize_frame(i)
}


