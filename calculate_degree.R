# 2018-07-25 (based on script 2017-04-10)
# marieke kuijjer
# this script combines edges with net and then calculates in- and outdegree and saves these in "degree.RData"

# make an edges.RData object
# make sure the ordering of this object is correct
# prednet is of dimensions [regulator, gene], it is then transformed with unlist(c(prednet))
	# reg1	gene1
	# reg2	gene1
	# reg3	gene1, etc
# prior is in:
edges <- read.delim("input/puma_prior.txt", header=F)
# sort these alphabetically as that is being done by panda_run.m as well
edges <- edges[order(edges[,2]),]
test <- unique(edges[,1])
all(test==test[order(test)]) # TRUE
test <- unique(edges[,2])
all(test==test[order(test)]) # TRUE
save(list="edges", file="edges.RData")

library(plyr)
dirname <- getwd()
dirname <- unlist(strsplit(dirname, "/"))
dirname <- dirname[length(dirname)]

load("net.RData")
sampleorderfile <- paste("input/sampleorder.txt", sep="")
samples <- scan(sampleorderfile, what="character")
colnames(net) <- samples
net <- as.data.frame(net)
row.names(net) <- 1:nrow(net)
edges <- edges[,1:2]
edges <- as.data.frame(edges, stringsAsFactors=F)
colnames(edges) <- c("reg","tar")
row.names(edges) <- 1:nrow(edges)
net <- cbind(edges, net)

indegree <- ddply(net[,-1], .(tar), numcolwise(sum))
row.names(indegree) <- indegree[,1]; indegree <- indegree[,-1]
outdegree <- ddply(net[,-2], .(reg), numcolwise(sum))
row.names(outdegree) <- outdegree[,1]; outdegree <- outdegree[,-1]
save(list=c("indegree", "outdegree"), file="degree.RData")
save(list=c("net"), file="net_new.RData")

### now load all degree.RData objects and combine into one
