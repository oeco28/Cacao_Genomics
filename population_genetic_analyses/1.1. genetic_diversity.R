library("plyr")
library("RColorBrewer")
library("ggplot2")
args<-commandArgs(TRUE)
setwd("/path_to_summaries/")
group = c("Amelonado","Contamana","Criollo","Curaray","Guianna","Iquitos","Maranon","Nacional","Nanay","Purus","Admixed")

data <- c()
for(i in 1:length(group)){
for(chrom in 1:10){
		data_temp <- read.table(paste(group[i],"/chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
		data_temp <- cbind(data_temp,rep(paste(group[i],sep=""),dim(data_temp)[1]))
		colnames(data_temp) <- c("CHROM","BIN_START","BIN_END","N_VARIANTS","PI","Group")
		data <- rbind(data,data_temp)
}
}

mypal <- brewer.pal(10,"Spectral")
mypal2 <- c(mypal[5],mypal[4],mypal[1],mypal[2],mypal[8],mypal[9],mypal[6],mypal[7],mypal[3],mypal[10],"snow4")

pdf(file="summary_piPerGroup3.pdf")
ggplot(data,aes(factor(Group),PI)) + geom_violin(aes(fill = factor(Group))) + scale_fill_manual(values=mypal2) + ylab(expression(pi)) + theme(axis.title.y = element_text(size = rel(1.8)))+ theme(legend.position="none") + xlab("Genetic Group") + theme(axis.title.x = element_text(size = rel(1)), axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

logpi_mean <- aggregate(log(data$PI,base=10),list(data$Group),mean)
logpi_se <- aggregate(log(data$PI,base=10),list(data$Group),sd)
ci_logpi <- logpi_se$x/221157*1.96
logpi_summary <- cbind(logpi_mean$x,logpi_mean$x - ci_logpi,logpi_mean$x + ci_logpi)
rownames(logpi_summary) <- logpi_mean$Group.1
colnames(logpi_summary) <- c(paste("mean log",expression(pi)),"L_CI 95%","H_CI 95%")
write.table(logpi_summary,file="logpi.summary.group.tab",quote=FALSE,row.names=TRUE,col.names=TRUE)

pi_mean <- aggregate(data$PI,list(data$Group),mean)
pi_se <- aggregate(data$PI,list(data$Group),sd)
ci_pi <- pi_se$x/221157*1.96
pi_summary <- cbind(pi_mean$x,pi_mean$x - ci_pi,pi_mean$x + ci_pi)
rownames(pi_summary) <- pi_mean$Group.1
colnames(pi_summary) <- c(paste(expression(pi)),"L_CI 95%","H_CI 95%")
write.table(pi_summary,file="pi.summary.group.tab",quote=FALSE,row.names=TRUE,col.names=TRUE)

log.glm <- glm(log(PI) ~ Group, family=gaussian, data=data)

##################################
##  Diversity plots along chromosomes
##  
##################################

library("plyr")
library("RColorBrewer")
library(ggplot2)
groups = c("Amelonado","Contamana","Criollo","Curaray","Guianna","Iquitos","Maranon","Nacional","Nanay","Purus","Admixed")

lengths <- read.table("/path_to_chr_lengths.txt",header=FALSE)

################
# Amelonado
################

group1 <- group[1]
setwd("/path_to_summaries/Amelonado")
data1 <- c()
midpoint <- c()
chrom <- 1
data1 <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
for(i in 1:dim(data1)[1]){
        midpoint[i] <- ((data1$BIN_END[i] - data1$BIN_START[i] -1)/2 + data1$BIN_START[i]);
    }
data1 <- cbind(data1,midpoint,rep(chrom,dim(data1)[1]),rep(paste(group1,sep=""),dim(data1)[1]))
colnames(data1) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")

for(chrom in 2:10){
    midpoint <- c()
    data_temp <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
    for(i in 1:dim(data_temp)[1]){
        midpoint[i] <- (((data_temp$BIN_END[i] - data_temp$BIN_START[i] -1)/2 + data_temp$BIN_START[i]) + sum(lengths[1:(chrom-1),2]));
    }
  data_temp2 <- cbind(data_temp,midpoint,rep(chrom,dim(data_temp)[1]),rep(paste(group1,sep=""),dim(data_temp)[1]))
  colnames(data_temp2) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")
  data1 <- rbind(data1,data_temp2)
}

#################
# Contamana
#################

group1 <- group[2]
setwd("/path_to_summaries/Contamana")
data2 <- c()
midpoint <- c()
chrom <- 1
data2 <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
for(i in 1:dim(data2)[1]){
        midpoint[i] <- ((data2$BIN_END[i] - data2$BIN_START[i] -1)/2 + data2$BIN_START[i]);
    }
data2 <- cbind(data2,midpoint,rep(chrom,dim(data2)[1]),rep(paste(group1,sep=""),dim(data2)[1]))
colnames(data2) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")

for(chrom in 2:10){
    midpoint <- c()
    data_temp <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
    for(i in 1:dim(data_temp)[1]){
        midpoint[i] <- ((data_temp$BIN_END[i] - data_temp$BIN_START[i] -1)/2 + data_temp$BIN_START[i]) + sum(lengths[1:(chrom-1),2]);
    }
  data_temp <- cbind(data_temp,midpoint,rep(chrom,dim(data_temp)[1]),rep(paste(group1,sep=""),dim(data_temp)[1]))
  colnames(data_temp) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")
  data2 <- rbind(data2,data_temp)
}


#################
# Criollo
#################

group1 <- group[3]
setwd("/path_to_summaries/Criollo")
data3 <- c()
midpoint <- c()
chrom <- 1
data3 <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
for(i in 1:dim(data3)[1]){
        midpoint[i] <- ((data3$BIN_END[i] - data3$BIN_START[i] -1)/2 + data3$BIN_START[i]);
    }
data3 <- cbind(data3,midpoint,rep(chrom,dim(data3)[1]),rep(paste(group1,sep=""),dim(data3)[1]))
colnames(data3) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")

for(chrom in 2:10){
    midpoint <- c()
    data_temp <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
    for(i in 1:dim(data_temp)[1]){
        midpoint[i] <- ((data_temp$BIN_END[i] - data_temp$BIN_START[i] -1)/2 + data_temp$BIN_START[i]) + sum(lengths[1:(chrom-1),2]);
    }
  data_temp <- cbind(data_temp,midpoint,rep(chrom,dim(data_temp)[1]),rep(paste(group1,sep=""),dim(data_temp)[1]))
  colnames(data_temp) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")
  data3 <- rbind(data3,data_temp)
}

#################
# Curaray
#################

group1 <- group[4]
setwd("/path_to_summaries/Curaray")
data4 <- c()
midpoint <- c()
chrom <- 1
data4 <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
for(i in 1:dim(data4)[1]){
        midpoint[i] <- ((data4$BIN_END[i] - data4$BIN_START[i] -1)/2 + data4$BIN_START[i]);
    }
data4 <- cbind(data4,midpoint,rep(chrom,dim(data4)[1]),rep(paste(group1,sep=""),dim(data4)[1]))
colnames(data4) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")

for(chrom in 2:10){
    midpoint <- c()
    data_temp <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
    for(i in 1:dim(data_temp)[1]){
        midpoint[i] <- ((data_temp$BIN_END[i] - data_temp$BIN_START[i] -1)/2 + data_temp$BIN_START[i]) + sum(lengths[1:(chrom-1),2]);
    }
  data_temp <- cbind(data_temp,midpoint,rep(chrom,dim(data_temp)[1]),rep(paste(group1,sep=""),dim(data_temp)[1]))
  colnames(data_temp) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")
  data4 <- rbind(data4,data_temp)
}

#################
# Guianna
#################

group1 <- group[5]
setwd("/path_to_summaries/Guianna")
data5 <- c()
midpoint <- c()
chrom <- 1
data5 <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
for(i in 1:dim(data5)[1]){
        midpoint[i] <- ((data5$BIN_END[i] - data5$BIN_START[i] -1)/2 + data5$BIN_START[i]);
    }
data5 <- cbind(data5,midpoint,rep(chrom,dim(data5)[1]),rep(paste(group1,sep=""),dim(data5)[1]))
colnames(data5) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")

for(chrom in 2:10){
    midpoint <- c()
    data_temp <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
    for(i in 1:dim(data_temp)[1]){
        midpoint[i] <- ((data_temp$BIN_END[i] - data_temp$BIN_START[i] -1)/2 + data_temp$BIN_START[i]) + sum(lengths[1:(chrom-1),2]);
    }
  data_temp <- cbind(data_temp,midpoint,rep(chrom,dim(data_temp)[1]),rep(paste(group1,sep=""),dim(data_temp)[1]))
  colnames(data_temp) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")
  data5 <- rbind(data5,data_temp)
}

#################
# Iquitos
#################

group1 <- group[6]
setwd("/path_to_summaries/Iquitos")
data6 <- c()
midpoint <- c()
chrom <- 1
data6 <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
for(i in 1:dim(data6)[1]){
        midpoint[i] <- ((data6$BIN_END[i] - data6$BIN_START[i] -1)/2 + data6$BIN_START[i]);
    }
data6 <- cbind(data6,midpoint,rep(chrom,dim(data6)[1]),rep(paste(group1,sep=""),dim(data6)[1]))
colnames(data6) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")

for(chrom in 2:10){
    midpoint <- c()
    data_temp <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
    for(i in 1:dim(data_temp)[1]){
        midpoint[i] <- ((data_temp$BIN_END[i] - data_temp$BIN_START[i] -1)/2 + data_temp$BIN_START[i]) + sum(lengths[1:(chrom-1),2]);
    }
  data_temp <- cbind(data_temp,midpoint,rep(chrom,dim(data_temp)[1]),rep(paste(group1,sep=""),dim(data_temp)[1]))
  colnames(data_temp) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")
  data6 <- rbind(data6,data_temp)
}

#################
# Maranon
#################

group1 <- group[7]
setwd("/path_to_summaries/Maranon")
data7 <- c()
midpoint <- c()
chrom <- 1
data7 <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
for(i in 1:dim(data7)[1]){
        midpoint[i] <- ((data7$BIN_END[i] - data7$BIN_START[i] -1)/2 + data7$BIN_START[i]);
    }
data7 <- cbind(data7,midpoint,rep(chrom,dim(data7)[1]),rep(paste(group1,sep=""),dim(data7)[1]))
colnames(data7) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")

for(chrom in 2:10){
    midpoint <- c()
    data_temp <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
    for(i in 1:dim(data_temp)[1]){
        midpoint[i] <- ((data_temp$BIN_END[i] - data_temp$BIN_START[i] -1)/2 + data_temp$BIN_START[i]) + sum(lengths[1:(chrom-1),2]);
    }
  data_temp <- cbind(data_temp,midpoint,rep(chrom,dim(data_temp)[1]),rep(paste(group1,sep=""),dim(data_temp)[1]))
  colnames(data_temp) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")
  data7 <- rbind(data7,data_temp)
}

#################
# Nacional
#################

group1 <- group[8]
setwd("/path_to_summaries/Nacional")
data8 <- c()
midpoint <- c()
chrom <- 1
data8 <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
for(i in 1:dim(data8)[1]){
        midpoint[i] <- ((data8$BIN_END[i] - data8$BIN_START[i] -1)/2 + data8$BIN_START[i]);
    }
data8 <- cbind(data8,midpoint,rep(chrom,dim(data8)[1]),rep(paste(group1,sep=""),dim(data8)[1]))
colnames(data8) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")

for(chrom in 2:10){
    midpoint <- c()
    data_temp <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
    for(i in 1:dim(data_temp)[1]){
        midpoint[i] <- ((data_temp$BIN_END[i] - data_temp$BIN_START[i] -1)/2 + data_temp$BIN_START[i]) + sum(lengths[1:(chrom-1),2]);
    }
  data_temp <- cbind(data_temp,midpoint,rep(chrom,dim(data_temp)[1]),rep(paste(group1,sep=""),dim(data_temp)[1]))
  colnames(data_temp) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")
  data8 <- rbind(data8,data_temp)
}

#################
# Nanay
#################

group1 <- group[9]
setwd("/path_to_summaries/Nanay")
data9 <- c()
midpoint <- c()
chrom <- 1
data9 <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
for(i in 1:dim(data9)[1]){
        midpoint[i] <- ((data9$BIN_END[i] - data9$BIN_START[i] -1)/2 + data9$BIN_START[i]);
    }
data9 <- cbind(data9,midpoint,rep(chrom,dim(data9)[1]),rep(paste(group1,sep=""),dim(data9)[1]))
colnames(data9) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")

for(chrom in 2:10){
    midpoint <- c()
    data_temp <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
    for(i in 1:dim(data_temp)[1]){
        midpoint[i] <- ((data_temp$BIN_END[i] - data_temp$BIN_START[i] -1)/2 + data_temp$BIN_START[i]) + sum(lengths[1:(chrom-1),2]);
    }
  data_temp <- cbind(data_temp,midpoint,rep(chrom,dim(data_temp)[1]),rep(paste(group1,sep=""),dim(data_temp)[1]))
  colnames(data_temp) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")
  data9 <- rbind(data9,data_temp)
}

#################
# Purus
#################

group1 <- group[10]
setwd("/path_to_summaries/Purus")
data10 <- c()
midpoint <- c()
chrom <- 1
data10 <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
for(i in 1:dim(data10)[1]){
        midpoint[i] <- ((data10$BIN_END[i] - data10$BIN_START[i] -1)/2 + data10$BIN_START[i]);
    }
data10 <- cbind(data10,midpoint,rep(chrom,dim(data10)[1]),rep(paste(group1,sep=""),dim(data10)[1]))
colnames(data10) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")

for(chrom in 2:10){
    midpoint <- c()
    data_temp <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
    for(i in 1:dim(data_temp)[1]){
        midpoint[i] <- ((data_temp$BIN_END[i] - data_temp$BIN_START[i] -1)/2 + data_temp$BIN_START[i]) + sum(lengths[1:(chrom-1),2]);
    }
  data_temp <- cbind(data_temp,midpoint,rep(chrom,dim(data_temp)[1]),rep(paste(group1,sep=""),dim(data_temp)[1]))
  colnames(data_temp) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")
  data10 <- rbind(data10,data_temp)
}

#################
# Admixed
#################

group1 <- group[11]
setwd("/path_to_summaries/Admixed")
data11 <- c()
midpoint <- c()
chrom <- 1
data11 <- read.table("chr1.1Kb.windowed.pi",header=TRUE,stringsAsFactors=F)
for(i in 1:dim(data11)[1]){
        midpoint[i] <- ((data11$BIN_END[i] - data11$BIN_START[i] -1)/2 + data11$BIN_START[i]);
    }
data11 <- cbind(data11,midpoint,rep(chrom,dim(data11)[1]),rep(paste(group1,sep=""),dim(data11)[1]))
colnames(data11) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")

for(chrom in 2:10){
    midpoint <- c()
    data_temp <- read.table(paste("chr",chrom,".1Kb.windowed.pi",sep=""),header=TRUE,stringsAsFactors=F)
    for(i in 1:dim(data_temp)[1]){
        midpoint[i] <- ((data_temp$BIN_END[i] - data_temp$BIN_START[i] -1)/2 + data_temp$BIN_START[i]) + sum(lengths[1:(chrom-1),2]);
    }
  data_temp <- cbind(data_temp,midpoint,rep(chrom,dim(data_temp)[1]),rep(paste(group1,sep=""),dim(data_temp)[1]))
  colnames(data_temp) <- c("scaffold","BIN_START","BIN_END","N_VARIANTS","PI","midpoint","chr","group")
  data11 <- rbind(data11,data_temp)
}

dataF <- rbind(data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11)

mypal <- brewer.pal(10,"Spectral")
mypalette_10 <- c(mypal[5],mypal[4],mypal[1],mypal[2],mypal[8],mypal[9],mypal[6],mypal[7],mypal[3],mypal[10],"darkgrey")
c <- ggplot(data = dataF, aes(x = midpoint, y = PI, colour=factor(group))) + geom_point(aes(colour=group), size=0.8) + theme(axis.text.x= element_text(size = rel(1.1), hjust = 1, vjust = 1)) + ylab(expression(pi)) + scale_colour_manual(values=mypalette_10) + guides(colour=FALSE) + scale_x_continuous("",breaks = c(lengths[1,2]/2,(lengths[1,2] + lengths[2,2]/2),(sum(lengths[1:2,2]) + lengths[3,2]/2),(sum(lengths[1:3,2]) + lengths[4,2]/2),(sum(lengths[1:4,2]) + lengths[5,2]/2),(sum(lengths[1:5,2]) + lengths[6,2]/2),(sum(lengths[1:6,2]) + lengths[7,2]/2),(sum(lengths[1:7,2]) + lengths[8,2]/2),(sum(lengths[1:8,2]) + lengths[9,2]/2),(sum(lengths[1:9,2]) + lengths[10,2]/2)), labels=c("chr 1","chr 2","chr 3","chr 4","chr 5","chr 6","chr 7","chr 8","chr 9","chr 10")) + facet_grid(group ~ .) + theme(strip.text.y = element_text(size = 7, colour = "black", angle = 90)) + labs(title = "Genetic diversity along chromosomes per group")
c + stat_smooth(method=lm)

