
library(ggplot2)
setwd("/path_to_xpclr_results/")

lengths <- read.table("/path_to/chr_lengths.txt",header=FALSE)

data1 <- read.table("cacao.xpclr.1.out.xpclr.txt",header=FALSE)
data2 <- read.table("cacao.xpclr.2.out.xpclr.txt",header=FALSE)
data3 <- read.table("cacao.xpclr.3.out.xpclr.txt",header=FALSE)
data4 <- read.table("cacao.xpclr.4.out.xpclr.txt",header=FALSE)
data5 <- read.table("cacao.xpclr.5.out.xpclr.txt",header=FALSE)
data6 <- read.table("cacao.xpclr.6.out.xpclr.txt",header=FALSE)
data7 <- read.table("cacao.xpclr.7.out.xpclr.txt",header=FALSE)
data8 <- read.table("cacao.xpclr.8.out.xpclr.txt",header=FALSE)
data9 <- read.table("cacao.xpclr.9.out.xpclr.txt",header=FALSE)
data10 <- read.table("cacao.xpclr.10.out.xpclr.txt",header=FALSE)
data <- rbind(data1,data2,data3,data4,data5,data6,data7,data8,data9,data10)
colnames(data) <- c("chr","grid#","#ofSNPs_in_window","physical_position","genetic_position","XPCLR_score","max_s")
chr_lab <- c(rep(1,dim(data1)[1]),rep(2,dim(data2)[1]),rep(1,dim(data3)[1]),rep(2,dim(data4)[1]),rep(1,dim(data5)[1]),rep(2,dim(data6)[1]),rep(1,dim(data7)[1]),rep(2,dim(data8)[1]),rep(1,dim(data9)[1]),rep(2,dim(data10)[1]))


norm_XPCLR <- c()
for(i in 1:dim(data)[1]){
		norm_XPCLR[i] <- (data$XPCLR_score[i] - mean(data$XPCLR_score))/(sd(data$XPCLR_score))
}

dataF <- cbind(data,norm_XPCLR)
dim(dataF[which(dataF$norm_XPCLR >= 4.2),])[1]/dim(dataF)[1]
dim(dataF[which(dataF$norm_XPCLR >= 4.2),])[1]

select_criollo <- dataF[which(dataF$norm_XPCLR >= 4),]
select_criollo2 <- cbind(select_criollo$genom_position,select_criollo$genom_position,select_criollo$norm_XPCLR,select_criollo$norm_XPCLR)
for(i in 1:dim(select_criollo)[1]){select_criollo2[i,1] - 1500}
for(i in 1:dim(select_criollo)[1]){select_criollo2[i,2] + 1500}
colnames(select_criollo2) <- c("x1","x2","y1","y2")


annotation <- read.table("annotation",header=FALSE)
genes_SC <- c()
for(i in 1:dim(select_criollo2)[1]){
	temp <- annotation[which(annotation$chr == select_criollo2$chr[i] & annotation$start_position >= (select_criollo2$physical_position[i] - 5000) & annotation$start_position <= (select_criollo2$physical_position[i] + 5000)),]
	temp2 <- cbind(temp,rep(select_criollo2$XPCLR_score[i],dim(temp)[1]),rep(select_criollo2$norm_XPCLR[i],dim(temp)[1]))
	colnames(temp2) <- c(colnames(temp),"XPCLR_score","nomr_XPCLR")
	genes_SC <- rbind(genes_SC,temp2)
}






##########
# I have to work on this manually for now
##########

genom_position <- data1[,4]
#gp <- c()
#for(i in 2:10){
#	temp <- data2[which(data2$chr == i),]
#	for (j in 1:dim(temp)[1])
#	gp[i] <- temp$physical_position[j] + sum(lengths[c(1: i - 1),2])
#	genom_position <- c(genom_position,gp)
#}

gp <- c()
for(i in 1:dim(data2)[1]){gp[i] <- data2$V4[i] + sum(lengths[1,2])}
genom_position <-  c(genom_position,gp)

gp <- c()
for(i in 1:dim(data3)[1]){gp[i] <- data3$V4[i] + sum(lengths[1:2,2])}
genom_position <-  c(genom_position,gp)

gp <- c()
for(i in 1:dim(data4)[1]){gp[i] <- data4$V4[i] + sum(lengths[1:3,2])}
genom_position <-  c(genom_position,gp)

gp <- c()
for(i in 1:dim(data5)[1]){gp[i] <- data5$V4[i] + sum(lengths[1:4,2])}
genom_position <-  c(genom_position,gp)

gp <- c()
for(i in 1:dim(data6)[1]){gp[i] <- data6$V4[i] + sum(lengths[1:5,2])}
genom_position <-  c(genom_position,gp)

gp <- c()
for(i in 1:dim(data7)[1]){gp[i] <- data7$V4[i] + sum(lengths[1:6,2])}
genom_position <-  c(genom_position,gp)

gp <- c()
for(i in 1:dim(data8)[1]){gp[i] <- data8$V4[i] + sum(lengths[1:7,2])}
genom_position <-  c(genom_position,gp)

gp <- c()
for(i in 1:dim(data9)[1]){gp[i] <- data9$V4[i] + sum(lengths[1:8,2])}
genom_position <-  c(genom_position,gp)

gp <- c()
for(i in 1:dim(data10)[1]){gp[i] <- data10$V4[i] + sum(lengths[1:9,2])}
genom_position <-  c(genom_position,gp)

chr_lab <- c(rep(1,dim(data1)[1]),rep(2,dim(data2)[1]),rep(1,dim(data3)[1]),rep(2,dim(data4)[1]),rep(1,dim(data5)[1]),rep(2,dim(data6)[1]),rep(1,dim(data7)[1]),rep(2,dim(data8)[1]),rep(1,dim(data9)[1]),rep(2,dim(data10)[1]))

dataF <- cbind(dataF,genom_position,chr_lab)

select_criollo2 <- dataF[which(dataF$norm_XPCLR >= 3.8),]


#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
mycol <- c("royalblue3","lightblue3")
#add.alpha <- function(col, alpha=1){
#   if(missing(col))
#     stop("Please provide a vector of colours.")
#   apply(sapply(col, col2rgb)/255, 2, 
#                      function(x) 
#                        rgb(x[1], x[2], x[3], alpha=alpha))  
#}

#sc <- scale_fill_gradientn(colours = myPalette(100), limits=c(min(na.omit(data2$Fst)), max(na.omit(data2$Fst))), na.value =add.alpha("white",alpha=0.1))

p <- ggplot(data=dataF,aes(x=genom_position,y=norm_XPCLR, colour=chr_lab)) + geom_point() + ylab("normalized XPCLR") + xlab("") + scale_x_continuous("",breaks = c(lengths[1,2]/2,(lengths[1,2] + lengths[2,2]/2),(sum(lengths[1:2,2]) + lengths[3,2]/2),(sum(lengths[1:3,2]) + lengths[4,2]/2),(sum(lengths[1:4,2]) + lengths[5,2]/2),(sum(lengths[1:5,2]) + lengths[6,2]/2),(sum(lengths[1:6,2]) + lengths[7,2]/2),(sum(lengths[1:7,2]) + lengths[8,2]/2),(sum(lengths[1:8,2]) + lengths[9,2]/2),(sum(lengths[1:9,2]) + lengths[10,2]/2)), labels=c("chr 1","chr 2","chr 3","chr 4","chr 5","chr 6","chr 7","chr 8","chr 9","chr 10")) + guides(colour=FALSE) 
#p + geom_segment(aes(data=select_criollo2, x = x1, y = y1, xend = x2, yend = y2, colour = "segment"))
p + geom_point(aes(x = genom_position, y = norm_XPCLR), color="darkred", data=select_criollo) + geom_hline(yintercept=4.2, color="red")


p <- ggplot(data=dataF,aes(x=genom_position,y=norm_XPCLR, colour=chr_lab)) + geom_point(alpha = 0.3) + ylab("normalized XPCLR") + xlab("") + scale_x_continuous("",breaks = c(lengths[1,2]/2,(lengths[1,2] + lengths[2,2]/2),(sum(lengths[1:2,2]) + lengths[3,2]/2),(sum(lengths[1:3,2]) + lengths[4,2]/2),(sum(lengths[1:4,2]) + lengths[5,2]/2),(sum(lengths[1:5,2]) + lengths[6,2]/2),(sum(lengths[1:6,2]) + lengths[7,2]/2),(sum(lengths[1:7,2]) + lengths[8,2]/2),(sum(lengths[1:8,2]) + lengths[9,2]/2),(sum(lengths[1:9,2]) + lengths[10,2]/2)), labels=c("chr 1","chr 2","chr 3","chr 4","chr 5","chr 6","chr 7","chr 8","chr 9","chr 10")) + guides(colour=FALSE)
p + geom_point(aes(x = genom_position, y = norm_XPCLR), color="tomato3", data=select_criollo2) + geom_hline(yintercept=3.8, color="grey70")




annotation <- read.table("annotation",header=FALSE)
colnames(annotation) <- c("chr","start_position","end_position","geneID")
genes_SC <- c()
for(i in 1:dim(select_criollo2)[1]){
	temp <- annotation[which(annotation$chr == select_criollo2$chr[i] & annotation$start_position >= (select_criollo2$physical_position[i] - 5000) & annotation$start_position <= (select_criollo2$physical_position[i] + 5000)),]
	temp2 <- cbind(temp,rep(select_criollo2$XPCLR_score[i],dim(temp)[1]),rep(select_criollo2$norm_XPCLR[i],dim(temp)[1]),rep(select_criollo2$physical_position[i],dim(temp)[1]))
	colnames(temp2) <- c(colnames(temp),"XPCLR_score","nomr_XPCLR","physical_position")
	genes_SC <- rbind(genes_SC,temp2)
}

genes_SC2 <- c()
for(i in 1:dim(select_criollo2)[1]){
	temp <- annotation[which(annotation$chr == select_criollo2$chr[i] & annotation$end_position >= (select_criollo2$physical_position[i] - 5000) & annotation$end_position <= (select_criollo2$physical_position[i] + 5000)),]
	temp2 <- cbind(temp,rep(select_criollo2$XPCLR_score[i],dim(temp)[1]),rep(select_criollo2$norm_XPCLR[i],dim(temp)[1]),rep(select_criollo2$physical_position[i],dim(temp)[1]))
	colnames(temp2) <- c(colnames(temp),"XPCLR_score","nomr_XPCLR","physical_position")
	genes_SC2 <- rbind(genes_SC2,temp2)
}


# then we can find the genes annotated in regions under selection

library(plyr)
genes_SCT <- join(genes_SC,genes_SC2,by="geneID",type="full")
