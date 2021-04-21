## combinig the files
combiningFilesinDirectory <- function(filePath, file_list, n){
  dat1 <- lapply(file_list, readFile, filePath = filePath, colNum = n)
  dat2 <- do.call("cbind", dat1)
  rm.col <- seq(3,ncol(dat2),by=2)
  dat2 <- dat2[,-rm.col]
  rownames(dat2) <- dat2[,1]
  dat2 <- dat2[,-1]
  colnames(dat2) <- file_list
  return(dat2)
  rm(dat1)
  rm(dat2)
}

## read files
readFile <- function(filePath, fileName, colNum){
  fileName <- paste(filePath, fileName, sep = "")
  dat <- read.table(file = fileName, header = FALSE, sep = "\t", quote = "", 
                    na.strings = FALSE, stringsAsFactors = FALSE)
  dat <- dat[,c(1,colNum)]
  return(dat)
  rm(dat)
}
