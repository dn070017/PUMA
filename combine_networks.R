### still need to export the matlab results to know the aggregate network and TF/genes for each prednet, this is now done in calculate_degree.R

setwd("txt")
filenames <- list.files()
# order filenames on nr instead of on string
fileidx <- substr(filenames, 9, nchar(filenames)-4)
fileidx <- as.numeric(fileidx)
filenames <- cbind(filenames, fileidx)
filenames <- as.data.frame(filenames, stringsAsFactors=F)
filenames[,2] <- as.numeric(filenames[,2])
filenames <- filenames[order(filenames[,2]),]

for(i in 1:nrow(filenames)){
	prednet <- read.delim(filenames[i,1], header=F)
	prednet <- unlist(c(prednet))
	if(i==1){
		net <- as.matrix(prednet)
		colnames(net) <- i
	} else {
		net <- cbind(net, prednet)
		colnames(net)[i] <- i
	}
	print(i)
}
save(list=c("net", "filenames"), file="../net.RData")
