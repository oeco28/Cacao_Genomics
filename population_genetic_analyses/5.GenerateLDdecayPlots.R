# LD was estimated as D' for all pairs using the --hap-r2 flag in vcftools and filtering out MFA < 0.05


library(ggplot2)
library(RColorBrewer)

####
# Amelonado
####

full_data <- c()
for(chrom in 1:10){
	data <- read.table(paste("Amelonado/chr",chrom,".hap.LD.hap.ld",sep=""),header=TRUE)
	full_data <- rbind(full_data,data)
	}
data2 <- cbind(full_data,full_data$POS2 - full_data$POS1)
colnames(data2) <- c(colnames(full_data),"Dist")

distance <- data2$Dist
LD.data <- data2$R.2
n <- 22
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
#ld.df<-data.frame(distance,fpoints)
ld.df <- cbind(distance,fpoints,rep("Amelonado",length(distance))
colnames(ld.df) <- c("distance","fpoints","Cluster")
ld.df_Amelonado<-ld.df[order(ld.df$distance),]


####
# Contamana
####

full_data <- c()
for(chrom in 1:10){
	data <- read.table(paste("Contamana/chr",chrom,".hap.LD.hap.ld",sep=""),header=TRUE)
	full_data <- rbind(full_data,data)
	}
data2 <- cbind(full_data,full_data$POS2 - full_data$POS1)
colnames(data2) <- c(colnames(full_data),"Dist")

distance <- data2$Dist
LD.data <- data2$R.2
n <- 22
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
#ld.df<-data.frame(distance,fpoints)
ld.df <- cbind(distance,fpoints,rep("Contamana",length(distance))
colnames(ld.df) <- c("distance","fpoints","Cluster")
ld.df_Contamana<-ld.df[order(ld.df$distance),]


####
# Criollo
####

full_data <- c()
for(chrom in 1:10){
	data <- read.table(paste("Criollo/chr",chrom,".hap.LD.hap.ld",sep=""),header=TRUE)
	full_data <- rbind(full_data,data)
	}
data2 <- cbind(full_data,full_data$POS2 - full_data$POS1)
colnames(data2) <- c(colnames(full_data),"Dist")

distance <- data2$Dist
LD.data <- data2$R.2
n <- 22
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
#ld.df<-data.frame(distance,fpoints)
ld.df <- cbind(distance,fpoints,rep("Criollo",length(distance))
colnames(ld.df) <- c("distance","fpoints","Cluster")
ld.df_Criollo<-ld.df[order(ld.df$distance),]

####
# Curaray
####

full_data <- c()
for(chrom in 1:10){
	data <- read.table(paste("Curaray/chr",chrom,".hap.LD.hap.ld",sep=""),header=TRUE)
	full_data <- rbind(full_data,data)
	}
data2 <- cbind(full_data,full_data$POS2 - full_data$POS1)
colnames(data2) <- c(colnames(full_data),"Dist")

distance <- data2$Dist
LD.data <- data2$R.2
n <- 22
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
#ld.df<-data.frame(distance,fpoints)
ld.df <- cbind(distance,fpoints,rep("Curaray",length(distance))
colnames(ld.df) <- c("distance","fpoints","Cluster")
ld.df_Curaray<-ld.df[order(ld.df$distance),]

####
# Guianna
####

full_data <- c()
for(chrom in 1:10){
	data <- read.table(paste("Guianna/chr",chrom,".hap.LD.hap.ld",sep=""),header=TRUE)
	full_data <- rbind(full_data,data)
	}
data2 <- cbind(full_data,full_data$POS2 - full_data$POS1)
colnames(data2) <- c(colnames(full_data),"Dist")

distance <- data2$Dist
LD.data <- data2$R.2
n <- 22
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
#ld.df<-data.frame(distance,fpoints)
ld.df <- cbind(distance,fpoints,rep("Guianna",length(distance))
colnames(ld.df) <- c("distance","fpoints","Cluster")
ld.df_Guianna<-ld.df[order(ld.df$distance),]

####
# Iquitos
####

full_data <- c()
for(chrom in 1:10){
	data <- read.table(paste("Iquitos/chr",chrom,".hap.LD.hap.ld",sep=""),header=TRUE)
	full_data <- rbind(full_data,data)
	}
data2 <- cbind(full_data,full_data$POS2 - full_data$POS1)
colnames(data2) <- c(colnames(full_data),"Dist")

distance <- data2$Dist
LD.data <- data2$R.2
n <- 22
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
#ld.df<-data.frame(distance,fpoints)
ld.df <- cbind(distance,fpoints,rep("Iquitos",length(distance))
colnames(ld.df) <- c("distance","fpoints","Cluster")
ld.df_Iquitos<-ld.df[order(ld.df$distance),]

####
# Maranon
####

full_data <- c()
for(chrom in 1:10){
	data <- read.table(paste("Maranon/chr",chrom,".hap.LD.hap.ld",sep=""),header=TRUE)
	full_data <- rbind(full_data,data)
	}
data2 <- cbind(full_data,full_data$POS2 - full_data$POS1)
colnames(data2) <- c(colnames(full_data),"Dist")

distance <- data2$Dist
LD.data <- data2$R.2
n <- 22
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
#ld.df<-data.frame(distance,fpoints)
ld.df <- cbind(distance,fpoints,rep("Maranon",length(distance))
colnames(ld.df) <- c("distance","fpoints","Cluster")
ld.df_Maranon<-ld.df[order(ld.df$distance),]

####
# Nanay
####

full_data <- c()
for(chrom in 1:10){
	data <- read.table(paste("Nanay/chr",chrom,".hap.LD.hap.ld",sep=""),header=TRUE)
	full_data <- rbind(full_data,data)
	}
data2 <- cbind(full_data,full_data$POS2 - full_data$POS1)
colnames(data2) <- c(colnames(full_data),"Dist")

distance <- data2$Dist
LD.data <- data2$R.2
n <- 22
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
#ld.df<-data.frame(distance,fpoints)
ld.df <- cbind(distance,fpoints,rep("Nanay",length(distance))
colnames(ld.df) <- c("distance","fpoints","Cluster")
ld.df_Nanay<-ld.df[order(ld.df$distance),]

####
# Nacional
####

full_data <- c()
for(chrom in 1:10){
	data <- read.table(paste("Nacional/chr",chrom,".hap.LD.hap.ld",sep=""),header=TRUE)
	full_data <- rbind(full_data,data)
	}
data2 <- cbind(full_data,full_data$POS2 - full_data$POS1)
colnames(data2) <- c(colnames(full_data),"Dist")

distance <- data2$Dist
LD.data <- data2$R.2
n <- 22
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
#ld.df<-data.frame(distance,fpoints)
ld.df <- cbind(distance,fpoints,rep("Nacional",length(distance))
colnames(ld.df) <- c("distance","fpoints","Cluster")
ld.df_Nacional<-ld.df[order(ld.df$distance),]

####
# Purus
####

full_data <- c()
for(chrom in 1:10){
	data <- read.table(paste("Purus/chr",chrom,".hap.LD.hap.ld",sep=""),header=TRUE)
	full_data <- rbind(full_data,data)
	}
data2 <- cbind(full_data,full_data$POS2 - full_data$POS1)
colnames(data2) <- c(colnames(full_data),"Dist")

distance <- data2$Dist
LD.data <- data2$R.2
n <- 22
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
#ld.df<-data.frame(distance,fpoints)
ld.df <- cbind(distance,fpoints,rep("Purus",length(distance))
colnames(ld.df) <- c("distance","fpoints","Cluster")
ld.df_Purus<-ld.df[order(ld.df$distance),]

ld.df_Total <- rbind(ld.df_Amelonado,ld.df_Contamana,ld.df_Criollo,ld.df_Curaray,ld.df_Guianna,ld.df_Iquitos,ld.df_Maranon,ld.df_Nanay,ld.df_Nacional,ld.df_Purus)
ld.df_Total$Cluster <- factor(ld.df_Total$Cluster,levels=c("Amelonado","Contamana","Criollo","Guianna","Iquitos","Maranon","Nacional","Nanay","Purus","Curaray"))
mypalette10 <- brewer.pal(10,"Spectral")

pdf("LD_decay.perCluster.pdf")
ggplot(data = combined, aes(x = distance, y = fpoints)) + geom_line(aes(colour=Cluster)) + scale_fill_manual(values=mypalette) 
dev.off()

pdf("LD_decay.perCluster2.pdf")
plot(ld.df_Amelonado$distance,ld.df_Amelonado$fpoints,lwd=2,col=mypalette10[1])
lines(ld.df_Contamana$distance,ld.df_Contamana$fpoints, lwd=2,col=mypalette10[2])
lines(ld.df_Criollo$distance,ld.df_Criollo$fpoints, lwd=2,col=mypalette10[3])
lines(ld.df_Curaray$distance,ld.df_Curaray$fpoints, lwd=2,col=mypalette10[10])
lines(ld.df_Guianna$distance,ld.df_Guianna$fpoints, lwd=2,col=mypalette10[4])
lines(ld.df_Iquitos$distance,ld.df_Iquitos$fpoints, lwd=2,col=mypalette10[5])
lines(ld.df_Maranon$distance,ld.df_Maranon$fpoints, lwd=2,col=mypalette10[6])
lines(ld.df_Nacional$distance,ld.df_Nacional$fpoints, lwd=2,col=mypalette10[7])
lines(ld.df_Nanay$distance,ld.df_Nanay$fpoints, lwd=2,col=mypalette10[8])
lines(ld.df_Purus$distance,ld.df_Purus$fpoints, lwd=2,col=mypalette10[9])
legend('center','groups',c("Amelonado","Contamana","Criollo","Guianna","Iquitos","Maranon","Nacional","Nanay","Purus","Curaray"), lty = c(1,2,3), col=mypalette10,ncol=5,nrow=2,bty ="n")
dev.off()
