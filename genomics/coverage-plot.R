x <- read.table('mappings/NODE_20.depth.txt', sep='\t', header=FALSE,  strip.white=TRUE)

# examine the data
head(x)

# calculate average depth
mean(x[,3])

# std dev
sqrt(var(x[,3]))

# plot coverage
plot(x[,2], x[,3], xlab='postion', ylab='coverage')
