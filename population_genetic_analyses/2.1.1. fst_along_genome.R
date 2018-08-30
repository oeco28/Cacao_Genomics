library(ggplot2)
library(RColorBrewer)
lengths <- read.table("path_to/chr_lengths.txt",header=FALSE)
mypalette10 <- c(brewer.pal(10,"Spectral"),colors()[1])

library(RColorBrewer)
group1 <- args[[1]]
group2 <- args[[2]]
jgroup <- paste(group1,"_",group2,sep="")
setwd(paste("/path_to_fsts/",jgroup,sep=""))
#mycolor1 <- c("royalblue4","royalblue1","tomato3","grey50")
mycolor1 <- c(colors()[614],colors()[613],"darkgoldenrod3","grey50")
myquantile_50 <- c()
myquantile_99 <- c()
#pdf(paste("Fst.",jgroup,".Fsts.pdf",sep=""),width=10.491228,height=4.280702)
plot(0,0,xlim=c(0,sum(lengths)),ylim=c(0,1),ylab="pairwise Fst",xlab="genomic position",xaxt='n',col="white",bty='n',main=paste(jgroup," pairwise Fst",sep=""))
data <- data.frame(read.table("chr1.windows.5Kb.windowed.weir.fst",header=TRUE))
midpoint <- c()
for(i in 1:dim(data)[1]){
	midpoint[i] <- ((data$BIN_END[i] - data$BIN_START[i] -1)/2 + data$BIN_START[i])
}
data2 <- cbind(data,midpoint)
myquantile_50 <- c(myquantile_50,quantile(data2$MEAN_FST,0.5))
myquantile_99 <- c(myquantile_99,quantile(data2$MEAN_FST,0.99))
points(data2$midpoint,data2$MEAN_FST, col=mycolor1[1],pch=16)

data <- data.frame(read.table("chr2.windows.5Kb.windowed.weir.fst",header=TRUE))
midpoint <- c()
for(i in 1:dim(data)[1]){
	midpoint[i] <- (((data$BIN_END[i] - data$BIN_START[i] -1)/2 + data$BIN_START[i]) + lengths[1,2])
}
data2 <- cbind(data,midpoint)
points(data2$midpoint,data2$MEAN_FST, col=mycolor1[2],pch=16)
data <- data.frame(read.table("chr3.windows.5Kb.windowed.weir.fst",header=TRUE))
midpoint <- c()
for(i in 1:dim(data)[1]){
	midpoint[i] <- ((data$BIN_END[i] - data$BIN_START[i] -1)/2 + data$BIN_START[i] + sum(lengths[1:2,2]));
}
data2 <- cbind(data,midpoint)
myquantile_50 <- c(myquantile_50,quantile(data2$MEAN_FST,0.5))
myquantile_99 <- c(myquantile_99,quantile(data2$MEAN_FST,0.99))
points(data2$midpoint,data2$MEAN_FST, col=mycolor1[1],pch=16)

data <- data.frame(read.table("chr4.windows.5Kb.windowed.weir.fst",header=TRUE))
midpoint <- c()
for(i in 1:dim(data)[1]){
	midpoint[i] <- ((data$BIN_END[i] - data$BIN_START[i] -1)/2 + data$BIN_START[i] + sum(lengths[1:3,2]));
}
data2 <- cbind(data,midpoint)
myquantile_50 <- c(myquantile_50,quantile(data2$MEAN_FST,0.5))
myquantile_99 <- c(myquantile_99,quantile(data2$MEAN_FST,0.99))
points(data2$midpoint,data2$MEAN_FST, col=mycolor1[2],pch=16)

data <- data.frame(read.table("chr5.windows.5Kb.windowed.weir.fst",header=TRUE))
midpoint <- c()
for(i in 1:dim(data)[1]){
	midpoint[i] <- ((data$BIN_END[i] - data$BIN_START[i] -1)/2 + data$BIN_START[i] + sum(lengths[1:4,2]));
}
data2 <- cbind(data,midpoint)
myquantile_50 <- c(myquantile_50,quantile(data2$MEAN_FST,0.5))
myquantile_99 <- c(myquantile_99,quantile(data2$MEAN_FST,0.99))
points(data2$midpoint,data2$MEAN_FST, col=mycolor1[1],pch=16)

data <- data.frame(read.table("chr6.windows.5Kb.windowed.weir.fst",header=TRUE))
midpoint <- c()
for(i in 1:dim(data)[1]){
	midpoint[i] <- ((data$BIN_END[i] - data$BIN_START[i] -1)/2 + data$BIN_START[i] + sum(lengths[1:5,2]));
}
data2 <- cbind(data,midpoint)
myquantile_50 <- c(myquantile_50,quantile(data2$MEAN_FST,0.5))
myquantile_99 <- c(myquantile_99,quantile(data2$MEAN_FST,0.99))
points(data2$midpoint,data2$MEAN_FST, col=mycolor1[2],pch=16)

data <- data.frame(read.table("chr7.windows.5Kb.windowed.weir.fst",header=TRUE))
midpoint <- c()
for(i in 1:dim(data)[1]){
	midpoint[i] <- ((data$BIN_END[i] - data$BIN_START[i] -1)/2 + data$BIN_START[i] + sum(lengths[1:6,2]));
}
data2 <- cbind(data,midpoint)
myquantile_50 <- c(myquantile_50,quantile(data2$MEAN_FST,0.5))
myquantile_99 <- c(myquantile_99,quantile(data2$MEAN_FST,0.99))
points(data2$midpoint,data2$MEAN_FST, col=mycolor1[1],pch=16)

data <- data.frame(read.table("chr8.windows.5Kb.windowed.weir.fst",header=TRUE))
midpoint <- c()
for(i in 1:dim(data)[1]){
	midpoint[i] <- ((data$BIN_END[i] - data$BIN_START[i] -1)/2 + data$BIN_START[i] + sum(lengths[1:7,2]));
}
data2 <- cbind(data,midpoint)
myquantile_50 <- c(myquantile_50,quantile(data2$MEAN_FST,0.5))
myquantile_99 <- c(myquantile_99,quantile(data2$MEAN_FST,0.99))
points(data2$midpoint,data2$MEAN_FST, col=mycolor1[2],pch=16)

data <- data.frame(read.table("chr9.windows.5Kb.windowed.weir.fst",header=TRUE))
midpoint <- c()
for(i in 1:dim(data)[1]){
	midpoint[i] <- ((data$BIN_END[i] - data$BIN_START[i] -1)/2 + data$BIN_START[i] + sum(lengths[1:8,2]));
}
data2 <- cbind(data,midpoint)
myquantile_50 <- c(myquantile_50,quantile(data2$MEAN_FST,0.5))
myquantile_99 <- c(myquantile_99,quantile(data2$MEAN_FST,0.99))
points(data2$midpoint,data2$MEAN_FST, col=mycolor1[1],pch=16)

data <- data.frame(read.table("chr10.windows.5Kb.windowed.weir.fst",header=TRUE))
midpoint <- c()
for(i in 1:dim(data)[1]){
	midpoint[i] <- ((data$BIN_END[i] - data$BIN_START[i] -1)/2 + data$BIN_START[i] + sum(lengths[1:9,2]));
}
data2 <- cbind(data,midpoint)
myquantile_50 <- c(myquantile_50,quantile(data2$MEAN_FST,0.5))
myquantile_99 <- c(myquantile_99,quantile(data2$MEAN_FST,0.99))
points(data2$midpoint,data2$MEAN_FST, col=mycolor1[2],pch=16)
abline(h=mean(myquantile_50),col=mycolor1[4],lwd=2)
abline(h=max(myquantile_99),col=mycolor1[3],lwd=2)
#axis(1,labels=c("scaffold_1","scaffold_2","scaffold_3","scaffold_4","scaffold_5","scaffold_6","scaffold_7","scaffold_8","scaffold_9","scaffold_10"),at=c(lengths[1,2]/2,(lengths[1,2] + lengths[2,2]/2),(sum(lengths[1:2,2]) + lengths[3,2]/2),(sum(lengths[1:3,2]) + lengths[4,2]/2),(sum(lengths[1:4,2]) + lengths[5,2]/2),(sum(lengths[1:5,2]) + lengths[6,2]/2),(sum(lengths[1:6,2]) + lengths[7,2]/2),(sum(lengths[1:7,2]) + lengths[8,2]/2),(sum(lengths[1:8,2]) + lengths[9,2]/2),(sum(lengths[1:9,2]) + lengths[10,2]/2)),cex.axis=0.7)
axis(1,labels=c("chr 1","chr 2","chr 3","chr 4","chr 5","chr 6","chr 7","chr 8","chr 9","chr 10"),at=c(lengths[1,2]/2,(lengths[1,2] + lengths[2,2]/2),(sum(lengths[1:2,2]) + lengths[3,2]/2),(sum(lengths[1:3,2]) + lengths[4,2]/2),(sum(lengths[1:4,2]) + lengths[5,2]/2),(sum(lengths[1:5,2]) + lengths[6,2]/2),(sum(lengths[1:6,2]) + lengths[7,2]/2),(sum(lengths[1:7,2]) + lengths[8,2]/2),(sum(lengths[1:8,2]) + lengths[9,2]/2),(sum(lengths[1:9,2]) + lengths[10,2]/2)),cex.axis=0.9)
dev.off()
