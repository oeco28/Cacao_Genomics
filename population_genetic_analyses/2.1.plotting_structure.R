setwd("/path_to_admixture/")
library("plyr")
library("ggplot2")
library("reshape")
library("RColorBrewer")

mypalette10 <- brewer.pal(10,"Spectral")
mypalette11 <- brewer.pal(11,"Spectral")

mK10 <- read.table("K10_admixture.tab",header=FALSE)
mK11 <- read.table("K11_admixture.tab",header=FALSE)
#combined <- read.table("combined_data.tab",header=FALSE)
combined <- rbind(mK10,mK11)
colnames(combined) <- c("Accessions","ancestry","proportion","K")
combined$Accessions <- factor(combined$Accessions, levels=c("criollo","sp1","sp3","sp9","LCTEEN_141","CUR3_G37_A6","CUR3_G38_A8","SIL-1_G56_A6","CUR3_G39_A10","NA702","NA_286","NA_331","NA_92","POUND_10_B","POUND_7_B","EET_400","NA45","Pound_7","PMF_20","PMF_27","SCA_11","SCA_19","SCA_24.2","SCA_5","T695_SCA-6_A1","NH_40","NH_53","Catongo","Matina","SIAL169","mvP30","SIC806","SIAL70","SIAL84","MATINA_Tica2","REDAMEL_1_31","REDAMEL_1_27","TRD86","PA_120","PA_121","PA_137","PA_56","PA_70","PA289","MO_9","MO_4","PA_169","PA_107","PA_51","PA_150","PA_88","PA_218","Brisas-1","UF273_T2","AM_1_54","UF273_T1","ELP_20_A","GU255_V","GU_114_P","GU_175_P","GU_219_F","GU_222","GU_300_P","GU_307","GU_308A","IMC_14","IMC67","IMC_20","IMC_50","IMC_12","IMC_47","IMC_51","CAB_76_PL3","CAB_77_PL5","RB39PL1","RB_40","RB_47_PL3","CAB_71_PL3","OC_61","ICS_95","ICS_60","SC_1","MXC_67","ICS40","UF_11","UF_668","ICS39","UF_676","SCA6","BR25","KA2","K82","ICS_43","ICS_6","M07","EET_395","M01","ICS_1","MO_109","M05","CC_71","TRD_45","PlayaAlta1","CRU_89","SCA_10","PBC123","SPA_7","Mocorongo","M02","M04","ICA_70","TAP3_G70_A2","CCN51","CCN10","CCAT1119_EET544_A1","UF12","FSC_13","T680_A1","M06","TSH516","B_6_8","LP_3_40","CCAT1858_EET547","AM_2_18","AGU_3339_12","UNAP2_G78_A2","T686_LCT-368_A1","TSH1188","COCA_3370_5","B_6_3","EET_59","EET58_A1","EET95_A2","AMAZ_12","EET_58","EB-19-1S","EET446_A1","CATIE_1000","2126","2462","SNA0707","2699","JA_5_35","2416","CL_27_50","LP_4_48","AMAZ-14_G23_A1","TIP-1_G41_A1","EET400_Arbol-122","T675_A645_A1","EET451_A1","CUR3_G35_A1","LCT_EEN_83_S-8","TAP10_G12_A1","LCT46_A1","TAP6_G12_A1","LCT_EEN_46","T685_EET387_A1","EB2237_A1","MO_99","AMAZ-14_G24_A5","T684_EET233_A1","PA107_A1","T682_A1_D147","AMAZ_15_15","2367","CRUZ_7_14","NA_712","CS_146","CCAT4675_EET575","T678_B60_A1","NA_33_T2","SPEC_194_75","CRU_59","AMAZ-11_G18_A1","AMAZ-11_G21_A10","mvT85","SPEC_54_1","M_8","MAN_15_2","ICS_61","TSA654Zymo","2748","2076","EET397","EET462_A1","IMC_36","PA_13","BE_10","CRU_101","LX_32","LX_43","CCAT4998_EET577","CCAT4688_EET576_A1","SLC_4","JA_5_5","JA_5_36","EET_103","SJ_2_22","EET103_Borde"))

ggplot(data = combined, aes(x = Accessions, y = as.numeric(proportion)), fill=cond) + geom_bar(aes(fill=ancestry),stat="identity") + scale_fill_brewer(palette="Spectral") + theme(axis.text.x= element_text(size = rel(1.1), angle=45, hjust = 1, vjust = 1)) + ylab("Proportion of ancestry") + guides(fill=FALSE) + facet_grid(K ~ .) + theme(axis.text.y= element_text(size = rel(0.7))) + scale_x_discrete(breaks = c("sp1","CUR3_G38_A8","NA_92","SCA_19","SIC806","MO_9","AM_1_54","GU_219_F","IMC_50","RB39PL1","JA_5_35"), labels=c("Criollo","Curaray","Nanay","Contamana","Amelonado","Marañon","Nacional","Guianna","Iquitos","Purús","Admixed"))


####
# MDS analysis
####

library("plyr")
library("ggplot2")
library("reshape")
library("RColorBrewer")

info <- read.table("groups.ids",header=FALSE)
colnames(info) <- c("ID","Cluster")
info <- data.frame(info)
info$Cluster <- factor(info$Cluster,levels=c("Criollo","Curaray","Nanay","Contamana","Amelonado","Marañon","Nacional","Guianna","Iquitos","Purús","Admixed"))
info2 <- info[order(info$Cluster),]
mylevels <- info2$ID

mypal <- brewer.pal(10,"Spectral")

forpc <- read.table("cacao.thin2.012",header=TRUE)

forpc2 <- forpc[match(mylevels,rownames(forpc)),]
forpc2 <- scale(forpc2, center=TRUE,scale=TRUE)
mymds <- cmdscale(forpc2, k=2)

C <- cbind(mymds[,c(1:2)],info2$Cluster)
colnames(C) <-  c("PC1","PC2","Cluster")
mycolours <- rep(NA,200)
mycolours[which(C$Cluster == "Criollo")] <- mypal[1]
mycolours[which(C$Cluster == "Curaray")] <- mypal[2]
mycolours[which(C$Cluster == "Nanay")] <- mypal[3]
mycolours[which(C$Cluster == "Contamana")] <- mypal[4]
mycolours[which(C$Cluster == "Amelonado")] <- mypal[5]
mycolours[which(C$Cluster == "Marañon")] <- mypal[6]
mycolours[which(C$Cluster == "Nacional")] <- mypal[7]
mycolours[which(C$Cluster == "Guianna")] <- mypal[8]
mycolours[which(C$Cluster == "Iquitos")] <- mypal[9]
mycolours[which(C$Cluster == "Purús")] <- mypal[10]
mycolours[which(C$Cluster == "Admixed")] <- "snow4"

mypal2 <- c(brewer.pal(10,"Spectral"),"ivory2")

p <- ggplot(D, aes(x=-PC2, y=-PC1)) + geom_point(aes(colour = factor(Cluster)),size=2.5) + scale_colour_manual(values=mypal2) + labs(colour="Clusters") + xlab("PC2") + ylab("PC1") + theme(axis.text = element_text(size = 11))
p + annotate("text", x=0.11, y=-0.02, label="East") + annotate("text", x=-0.15, y=-0.02, label="West") + annotate("text", x=-0.17, y=-0.05, label="South") + annotate("text", x=-0.17, y=0.3, label="North")


#####################
## calculating centroids and anlysis of genetic diversity decay along longitude West -> East gradient
###################

centroids <- aggregate(cbind(C$PC1,C$PC2)~Cluster,C,mean) # thi is a function to estimate the mean accross clusters.
colnames(centroids) <- c("Cluster","PC1","PC2")

p <- ggplot(D, aes(x=-PC1, y=PC2)) + geom_point(aes(colour = factor(Cluster)),size=2.5) + scale_colour_manual(values=mypal2) + labs(colour="Clusters") + xlab("PC1")
p + geom_point(data=centroids,size=2,colour=mypal2)

pi_cluster <- read.table("/path_to_summaries/pi.summary.group2.tab",header=TRUE)
A <- join(centroids,pi_cluster,by="Cluster",type="inner")
A2 <- A[-1,]
A2 <- A2[-10,]
mypal3 <- mypal2[2:10]
q <- ggplot(A2, aes(x=-PC2, y=pi)) + geom_point(aes(colour = factor(Cluster)),size=3) + scale_colour_manual(values=mypal3) + labs(colour="Clusters") + xlab("PC2") + ylab(expression(pi)) + geom_smooth(method = "lm", se = TRUE,colour="darkgrey") + theme(axis.text = element_text(size = 11),axis.title=element_text(size=13))
q + annotate("text", x=-0.12, y=0.001, label="West") + annotate("text", x=0.10, y=0.001, label="East")


lm(formula = -PC2 ~ pi, data = A2)


######
# Plotting map with coordinates of populations
#####

setwd("/path_to_data/pca")
library(RColorBrewer)
library(maps)

map(database="world",regions=c('brazil','ecuador','bolivia','peru','venezuela','guyana','surinam','colombia','panama'),xlim=c(-90,-30),ylim=c(-18,15),col="grey95",fill=TRUE)
map.axes()
points(q[which(ids2 == "snow4"),1],q[which(ids2 == "snow4"),2],col="snow4",pch=17)
points(q[which(ids2 == mypalette[10]),1],q[which(ids2 == mypalette[10]),2],col=mypalette[10],pch=17)
points(q[which(ids2 == mypalette[2]),1],q[which(ids2 == mypalette[2]),2],col=mypalette[2],pch=17)
points(q[which(ids2 == mypalette[8]),1],q[which(ids2 == mypalette[8]),2],col=mypalette[8],pch=17)
points(q[which(ids2 == mypalette[9]),1],q[which(ids2 == mypalette[9]),2],col=mypalette[9],pch=17)
points(q[which(ids2 == mypalette[6]),1],q[which(ids2 == mypalette[6]),2],col=mypalette[6],pch=17)
points(q[which(ids2 == mypalette[3]),1],q[which(ids2 == mypalette[3]),2],col=mypalette[3],pch=17)
points(q[which(ids2 == mypalette[5]),1],q[which(ids2 == mypalette[5]),2],col=mypalette[5],pch=17)
points(q[which(ids2 == mypalette[1]),1],q[which(ids2 == mypalette[1]),2],col=mypalette[1],pch=17)

poploc <- matrix(NA,8,4)
poploc[1,] <- c(-80.5,10.32,"Criollo","#9E0142")
poploc[2,] <- c(-76.88,-1.955,"Curaray","#D53E4F")
poploc[3,] <- c(-43,-12,"Amelonado","#FEE08B")
poploc[4,] <- c(-55.71,1.29,"Guianna","#ABDDA4")
poploc[5,] <- c(-73,-3.5,"Iquitos","#66C2A5")
poploc[6,] <- c(-75.15,-4.65,"Maranon","#FFFFBF") # #ffefbf #FFFFBF
poploc[7,] <- c(-75.55,-0.73,"Nanay","#F46D43")
poploc[8,] <- c(-68.25,-9.20,"Purus","#3288BD")

# alternatively  

mp <- ggmap(map)
mp + geom_point(data=q2, aes(x=Long, y=Lat, size=2), color=ids3, size=3, alpha=0.8) + annotate("text", x = -67.5, y = 5.32, label = "Criollo",colour ="#9E0142") + annotate("text", x = -76.88, y = -1.955, label = "Curaray", colour = "#D53E4F") + annotate("text", x = -45, y = -17, label = "Amelonado", colour = "yellow") + annotate("text", x = -55.71, y = 1.29, label = "Guianna", colour = "#ABDDA4") + annotate("text", x = -73, y = -3.5, label = "Iquitos", colour = "#66C2A5") + annotate("text", x = -75.15, y = -4.65, label = "Maranon", colour = "#FFFFBF") + annotate("text", x = -75.55, y = -0.73, label = "Nanay", colour = "#F46D43") + annotate("text", x = -68.25, y = -9.20, label = "Purus", colour = "#3288BD")



