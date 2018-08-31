setwd("/working_path/")

####### 
# Plotting PSMC per group
########
library(RColorBrewer)
library(ggplot2)

mypal <- brewer.pal(10,"Spectral")

###############
# load amelonado
################

Catongo <- read.table("Catongo.0.txt",header=FALSE)
MATINA_Tica2 <- read.table("MATINA_Tica2.0.txt",header=FALSE)
Matina <- read.table("Matina.0.txt",header=FALSE)
REDAMEL_1_27  <- read.table("REDAMEL_1_27.0.txt",header=FALSE)
REDAMEL_1_31  <- read.table("REDAMEL_1_31.0.txt",header=FALSE)
SIAL169 <- read.table("SIAL169.0.txt",header=FALSE)
SIAL70 <- read.table("SIAL70.0.txt",header=FALSE)
SIAL84 <- read.table("SIAL84.0.txt",header=FALSE)
SIC806 <- read.table("SIC806.0.txt",header=FALSE)
TRD86 <- read.table("TRD86.0.txt",header=FALSE)
mvP30 <- read.table("mvP30.0.txt",header=FALSE)


Catongo_2 <- cbind(Catongo[,1:2],rep("Catongo",dim(Catongo)[1]),rep("Amelonado",dim(Catongo)[1]))
colnames(Catongo_2) <- c("Time","Ne","accession","Cluster")
MATINA_Tica2_2 <- cbind(MATINA_Tica2[,1:2],rep("MATINA_Tica2",dim(MATINA_Tica2)[1]),rep("Amelonado",dim(MATINA_Tica2)[1]))
colnames(MATINA_Tica2_2) <- c("Time","Ne","accession","Cluster")
Matina_2 <- cbind(Matina[,1:2],rep("Matina",dim(Matina)[1]),rep("Amelonado",dim(Matina)[1]))
colnames(Matina_2) <- c("Time","Ne","accession","Cluster")
REDAMEL_1_27_2 <- cbind(REDAMEL_1_27[,1:2],rep("REDAMEL_1_27",dim(REDAMEL_1_27)[1]),rep("Amelonado",dim(REDAMEL_1_27)[1]))
colnames(REDAMEL_1_27_2) <- c("Time","Ne","accession","Cluster")
REDAMEL_1_31_2 <- cbind(REDAMEL_1_31[,1:2],rep("REDAMEL_1_31",dim(REDAMEL_1_31)[1]),rep("Amelonado",dim(REDAMEL_1_31)[1]))
colnames(REDAMEL_1_31_2) <- c("Time","Ne","accession","Cluster")
SIAL169_2 <- cbind(SIAL169[,1:2],rep("SIAL169",dim(SIAL169)[1]),rep("Amelonado",dim(SIAL169)[1]))
colnames(SIAL169_2) <- c("Time","Ne","accession","Cluster")
SIAL70_2 <- cbind(SIAL70[,1:2],rep("SIAL70",dim(SIAL70)[1]),rep("Amelonado",dim(SIAL70)[1]))
colnames(SIAL70_2) <- c("Time","Ne","accession","Cluster")
SIAL84_2 <- cbind(SIAL84[,1:2],rep("SIAL84",dim(SIAL84)[1]),rep("Amelonado",dim(SIAL84)[1]))
colnames(SIAL84_2) <- c("Time","Ne","accession","Cluster")
SIC806_2 <- cbind(SIC806[,1:2],rep("SIC806",dim(SIC806)[1]),rep("Amelonado",dim(SIC806)[1]))
colnames(SIC806_2) <- c("Time","Ne","accession","Cluster")
TRD86_2 <- cbind(TRD86[,1:2],rep("TRD86",dim(TRD86)[1]),rep("Amelonado",dim(TRD86)[1]))
colnames(TRD86_2) <- c("Time","Ne","accession","Cluster")
mvP30_2 <- cbind(mvP30[,1:2],rep("mvP30",dim(mvP30)[1]),rep("Amelonado",dim(mvP30)[1]))
colnames(mvP30_2) <- c("Time","Ne","accession","Cluster")
amelonado_b <- rbind(Catongo_2,MATINA_Tica2_2,Matina_2,REDAMEL_1_27_2,REDAMEL_1_31_2,SIAL169_2,SIAL70_2,SIAL84_2,SIC806_2,TRD86_2,mvP30_2)


amelonado <- data.frame(matrix(NA,dim(Catongo)[1],3))
colnames(amelonado) <- c("Time","Ne","Cluster")
for(i in 1:dim(Catongo)[1]){amelonado[i,1] <- mean(c(Catongo[i,1],MATINA_Tica2[i,1],Matina[i,1],REDAMEL_1_31[i,1],REDAMEL_1_27[i,1],SIAL169[i,1],SIAL70[i,1],SIAL84[i,1],SIC806[i,1],TRD86[i,1],mvP30[i,1])) + 0.01}
for(i in 1:dim(Catongo)[1]){amelonado[i,2] <- mean(c(Catongo[i,2],MATINA_Tica2[i,2],Matina[i,2],REDAMEL_1_31[i,2],REDAMEL_1_27[i,2],SIAL169[i,2],SIAL70[i,2],SIAL84[i,2],SIC806[i,2],TRD86[i,2],mvP30[i,2]))}
amelonado[,3] <- rep("Amelonado",dim(amelonado)[1])


###############
# load contamana
################

NH_40 <- read.table("NH_40.0.txt",header=FALSE)
NH_53 <- read.table("NH_40.0.txt",header=FALSE)
PMF_20 <- read.table("PMF_20.0.txt",header=FALSE)
PMF_27 <- read.table("PMF_27.0.txt",header=FALSE)
SCA_11 <- read.table("SCA_11.0.txt",header=FALSE)
SCA_19 <- read.table("SCA_19.0.txt",header=FALSE)
SCA_24 <- read.table("SCA_24.0.txt",header=FALSE)
SCA_5 <- read.table("SCA_5.0.txt",header=FALSE)
T695_SCA_6_A1 <- read.table("T695_SCA-6_A1.0.txt",header=FALSE)


NH_40_2 <- cbind(NH_40[,1:2],rep("NH_40",dim(NH_40)[1]),rep("Contamana",dim(NH_40)[1]))
colnames(NH_40_2) <- c("Time","Ne","accession","Cluster")
NH_53_2 <- cbind(NH_53[,1:2],rep("NH_53",dim(NH_53)[1]),rep("Contamana",dim(NH_53)[1]))
colnames(NH_53_2) <- c("Time","Ne","accession","Cluster")
PMF_20_2 <- cbind(PMF_20[,1:2],rep("PMF_20",dim(PMF_20)[1]),rep("Contamana",dim(PMF_20)[1]))
colnames(PMF_20_2) <- c("Time","Ne","accession","Cluster")
PMF_27_2 <- cbind(PMF_27[,1:2],rep("PMF_27",dim(PMF_27)[1]),rep("Contamana",dim(PMF_27)[1]))
colnames(PMF_27_2) <- c("Time","Ne","accession","Cluster")
SCA_11_2 <- cbind(SCA_11[,1:2],rep("SCA_11",dim(SCA_11)[1]),rep("Contamana",dim(SCA_11)[1]))
colnames(SCA_11_2) <- c("Time","Ne","accession","Cluster")
SCA_19_2 <- cbind(SCA_19[,1:2],rep("SCA_19",dim(SCA_19)[1]),rep("Contamana",dim(SCA_19)[1]))
colnames(SCA_19_2) <- c("Time","Ne","accession","Cluster")
SCA_24_2 <- cbind(SCA_24[,1:2],rep("SCA_24",dim(SCA_24)[1]),rep("Contamana",dim(SCA_24)[1]))
colnames(SCA_24_2) <- c("Time","Ne","accession","Cluster")
SCA_5_2 <- cbind(SCA_5[,1:2],rep("SCA_5",dim(SCA_5)[1]),rep("Contamana",dim(SCA_5)[1]))
colnames(SCA_5_2) <- c("Time","Ne","accession","Cluster")
T695_SCA_6_A1_2 <- cbind(T695_SCA_6_A1[,1:2],rep("T695_SCA_6_A1",dim(T695_SCA_6_A1)[1]),rep("Contamana",dim(T695_SCA_6_A1)[1]))
colnames(T695_SCA_6_A1_2) <- c("Time","Ne","accession","Cluster")
contamana_b <- rbind(NH_40_2,NH_53_2,PMF_20_2,PMF_27_2,SCA_11_2,SCA_19_2,SCA_24_2,SCA_5_2,T695_SCA_6_A1_2)


contamana <- data.frame(matrix(NA,dim(NH_40)[1],3))
colnames(contamana) <- c("Time","Ne","Cluster")
for(i in 1:dim(NH_40)[1]){contamana[i,1] <- mean(c(NH_40[i,1],NH_53[i,1],PMF_20[i,1],PMF_27[i,1],SCA_11[i,1],SCA_11[i,1],SCA_19[i,1],SCA_24[i,1],SCA_5[i,1],T695_SCA_6_A1[i,1])) + 0.01}
for(i in 1:dim(NH_40)[1]){contamana[i,2] <- mean(c(NH_40[i,2],NH_53[i,2],PMF_20[i,2],PMF_27[i,2],SCA_11[i,2],SCA_11[i,2],SCA_19[i,2],SCA_24[i,2],SCA_5[i,2],T695_SCA_6_A1[i,2]))}
contamana[,3] <- rep("Contamana",dim(contamana)[1])

###############
# load criollo
################
criollo <- read.table("criollo.0.txt",header=FALSE)
sp1 <- read.table("sp1.0.txt",header=FALSE)
sp3 <- read.table("sp3.0.txt",header=FALSE)
sp9 <- read.table("sp9.0.txt",header=FALSE)

criollo_2 <- cbind(CUR3_G37_A6[,1:2],rep("criollo",dim(criollo)[1]),rep("Criollo",dim(criollo)[1]))
colnames(criollo_2) <- c("Time","Ne","accession","Cluster")
sp1_2 <- cbind(sp1[,1:2],rep("sp1",dim(sp1)[1]),rep("Criollo",dim(sp1)[1]))
colnames(sp1_2) <- c("Time","Ne","accession","Cluster")
sp3_2 <- cbind(sp3[,1:2],rep("sp3",dim(sp3)[1]),rep("Criollo",dim(sp3)[1]))
colnames(sp3_2) <- c("Time","Ne","accession","Cluster")
sp9_2 <- cbind(sp9[,1:2],rep("sp9",dim(sp9)[1]),rep("Criollo",dim(sp9)[1]))
colnames(sp9_2) <- c("Time","Ne","accession","Cluster")
Criollo_b <- rbind(criollo_2,sp1_2,sp3_2,sp9_2)

Criollo <- data.frame(matrix(NA,dim(criollo)[1],3))
colnames(Criollo) <- c("Time","Ne","Cluster")
for(i in 1:dim(criollo)[1]){Criollo[i,1] <- mean(c(criollo[i,1],sp1[i,1],sp3[i,1],sp9[i,1])) + 0.01}
for(i in 1:dim(criollo)[1]){Criollo[i,2] <- mean(c(criollo[i,2],sp1[i,2],sp3[i,2],sp9[i,2]))}
Criollo[,3] <- rep("Criollo",dim(Criollo)[1])

###############
# load curaray
################
CUR3_G37_A6 <- read.table("CUR3_G37_A6.0.txt",header=FALSE)
CUR3_G38_A8 <- read.table("CUR3_G38_A8.0.txt",header=FALSE)
CUR3_G39_A10 <- read.table("CUR3_G39_A10.0.txt",header=FALSE)
LCTEEN_141 <- read.table("LCTEEN_141.0.txt",header=FALSE)
SIL_1_G56_A6 <- read.table("SIL-1_G56_A6.0.txt",header=FALSE)

CUR3_G37_A6_2 <- cbind(CUR3_G37_A6[,1:2],rep("CUR3_G37_A6",dim(CUR3_G37_A6)[1]),rep("Curaray",dim(CUR3_G37_A6)[1]))
colnames(CUR3_G37_A6_2) <- c("Time","Ne","accession","Cluster")
CUR3_G38_A8_2 <- cbind(CUR3_G38_A8[,1:2],rep("CUR3_G38_A8",dim(CUR3_G38_A8)[1]),rep("Curaray",dim(CUR3_G38_A8)[1]))
colnames(CUR3_G38_A8_2) <- c("Time","Ne","accession","Cluster")
CUR3_G39_A10_2 <- cbind(CUR3_G39_A10[,1:2],rep("CUR3_G39_A10",dim(CUR3_G39_A10)[1]),rep("Curaray",dim(CUR3_G39_A10)[1]))
colnames(CUR3_G39_A10_2) <- c("Time","Ne","accession","Cluster")
LCTEEN_141_2 <- cbind(LCTEEN_141[,1:2],rep("LCTEEN_141",dim(LCTEEN_141)[1]),rep("Curaray",dim(LCTEEN_141)[1]))
colnames(LCTEEN_141_2) <- c("Time","Ne","accession","Cluster")
SIL_1_G56_A6_2 <- cbind(SIL_1_G56_A6[,1:2],rep("SIL_1_G56_A6",dim(SIL_1_G56_A6)[1]),rep("Curaray",dim(SIL_1_G56_A6)[1]))
colnames(SIL_1_G56_A6_2) <- c("Time","Ne","accession","Cluster")
curaray_b <- rbind(CUR3_G37_A6_2,CUR3_G38_A8_2,CUR3_G39_A10_2,LCTEEN_141_2,SIL_1_G56_A6_2)

curaray <- data.frame(matrix(NA,dim(CUR3_G37_A6)[1],3))
colnames(curaray) <- c("Time","Ne","Cluster")
for(i in 1:dim(curaray)[1]){curaray[i,1] <- mean(c(CUR3_G37_A6[i,1],CUR3_G38_A8[i,1],CUR3_G39_A10[i,1],LCTEEN_141[i,1],SIL_1_G56_A6[i,1])) + 0.01}
for(i in 1:dim(curaray)[1]){curaray[i,2] <- mean(c(CUR3_G37_A6[i,2],CUR3_G38_A8[i,2],CUR3_G39_A10[i,2],LCTEEN_141[i,2],SIL_1_G56_A6[i,2]))}
curaray[,3] <- rep("Curaray",dim(curaray)[1])


###############
# load guianna
################
ELP_20_A <- read.table("ELP_20_A.0.txt",header=FALSE)
GU255_V <- read.table("GU255_V.2.0.txt",header=FALSE)
GU_114_P <- read.table("GU_114_P.0.txt",header=FALSE)
GU_175_P <- read.table("GU_175_P.2.0.txt",header=FALSE)
GU_219_F <- read.table("GU_219_F.0.txt",header=FALSE)
GU_222 <- read.table("GU_222.0.txt",header=FALSE)
GU_300_P <- read.table("GU_300_P.0.txt",header=FALSE)
GU_307 <- read.table("GU_307.0.txt",header=FALSE)
GU_308A <- read.table("GU_308A.0.txt",header=FALSE)

ELP_20_A_2 <- cbind(ELP_20_A[,1:2],rep("ELP_20_A",dim(ELP_20_A)[1]),rep("Guianna",dim(ELP_20_A)[1]))
colnames(ELP_20_A_2) <- c("Time","Ne","accession","Cluster")
GU255_V_2 <- cbind(GU255_V[,1:2],rep("GU255_V",dim(GU255_V)[1]),rep("Guianna",dim(GU255_V)[1]))
colnames(GU255_V_2) <- c("Time","Ne","accession","Cluster")
GU_114_P_2 <- cbind(GU_114_P[,1:2],rep("GU_114_P",dim(GU_114_P)[1]),rep("Guianna",dim(GU_114_P)[1]))
colnames(GU_114_P_2) <- c("Time","Ne","accession","Cluster")
GU_175_P_2 <- cbind(GU_175_P[,1:2],rep("GU_175_P",dim(GU_175_P)[1]),rep("Guianna",dim(GU_175_P)[1]))
colnames(GU_175_P_2) <- c("Time","Ne","accession","Cluster")
GU_219_F_2 <- cbind(GU_219_F[,1:2],rep("GU_219_F",dim(GU_219_F)[1]),rep("Guianna",dim(GU_219_F)[1]))
colnames(GU_219_F_2) <- c("Time","Ne","accession","Cluster")
GU_222_2 <- cbind(GU_222[,1:2],rep("GU_222",dim(GU_222)[1]),rep("Guianna",dim(GU_222)[1]))
colnames(GU_222_2) <- c("Time","Ne","accession","Cluster")
GU_300_P_2 <- cbind(GU_300_P[,1:2],rep("GU_300_P",dim(GU_300_P)[1]),rep("Guianna",dim(GU_300_P)[1]))
colnames(GU_300_P_2) <- c("Time","Ne","accession","Cluster")
GU_307_2 <- cbind(GU_307[,1:2],rep("GU_307",dim(GU_307)[1]),rep("Guianna",dim(GU_307)[1]))
colnames(GU_307_2) <- c("Time","Ne","accession","Cluster")
GU_308A_2 <- cbind(GU_308A[,1:2],rep("GU_308A",dim(GU_308A)[1]),rep("Guianna",dim(GU_308A)[1]))
colnames(GU_308A_2) <- c("Time","Ne","accession","Cluster")
guianna_b <- rbind(ELP_20_A_2,GU255_V_2,GU_114_P_2,GU_175_P_2,GU_219_F_2,GU_222_2,GU_300_P_2,GU_307_2,GU_308A_2)


guianna <- data.frame(matrix(NA,dim(ELP_20_A)[1],3))
colnames(guianna) <- c("Time","Ne","Cluster")
for(i in 1:dim(guianna)[1]){guianna[i,1] <- mean(c(ELP_20_A[i,1],GU255_V[i,1],GU_114_P[i,1],GU_175_P[i,1],GU_219_F[i,1],GU_222[i,1],GU_300_P[i,1],GU_307[i,1],GU_308A[i,1])) + 0.01}
for(i in 1:dim(guianna)[1]){guianna[i,2] <- mean(c(ELP_20_A[i,2],GU255_V[i,2],GU_114_P[i,2],GU_175_P[i,2],GU_219_F[i,2],GU_222[i,2],GU_300_P[i,2],GU_307[i,2],GU_308A[i,2]))}
guianna[,3] <- rep("Guianna",dim(guianna)[1])


###############
# load iquitos
################
IMC67 <- read.table("IMC67.0.txt",header=FALSE)
IMC_12 <- read.table("IMC_12.0.txt",header=FALSE)
IMC_14 <- read.table("IMC_14.0.txt",header=FALSE)
IMC_20 <- read.table("IMC_20.0.txt",header=FALSE)
IMC_47 <- read.table("IMC_47.0.txt",header=FALSE)
IMC_50 <- read.table("IMC_50.0.txt",header=FALSE)
IMC_51 <- read.table("IMC_51.0.txt",header=FALSE)

IMC67_2 <- cbind(IMC67[,1:2],rep("IMC67",dim(IMC67)[1]),rep("Iquitos",dim(IMC67)[1]))
colnames(IMC67_2) <- c("Time","Ne","accession","Cluster")
IMC_12_2 <- cbind(IMC_12[,1:2],rep("IMC_12",dim(IMC_12)[1]),rep("Iquitos",dim(IMC_12)[1]))
colnames(IMC_12_2) <- c("Time","Ne","accession","Cluster")
IMC_14_2 <- cbind(IMC_14[,1:2],rep("IMC_14",dim(IMC_14)[1]),rep("Iquitos",dim(IMC_14)[1]))
colnames(IMC_14_2) <- c("Time","Ne","accession","Cluster")
IMC_20_2 <- cbind(IMC_20[,1:2],rep("IMC_20",dim(IMC_20)[1]),rep("Iquitos",dim(IMC_20)[1]))
colnames(IMC_20_2) <- c("Time","Ne","accession","Cluster")
IMC_47_2 <- cbind(IMC_47[,1:2],rep("IMC_47",dim(IMC_47)[1]),rep("Iquitos",dim(IMC_47)[1]))
colnames(IMC_47_2) <- c("Time","Ne","accession","Cluster")
IMC_50_2 <- cbind(IMC_50[,1:2],rep("IMC_50",dim(IMC_50)[1]),rep("Iquitos",dim(IMC_50)[1]))
colnames(IMC_50_2) <- c("Time","Ne","accession","Cluster")
IMC_51_2 <- cbind(IMC_51[,1:2],rep("IMC_51",dim(IMC_51)[1]),rep("Iquitos",dim(IMC_51)[1]))
colnames(IMC_51_2) <- c("Time","Ne","accession","Cluster")
iquitos_b <- rbind(IMC67_2,IMC_12_2,IMC_14_2,IMC_20_2,IMC_47_2,IMC_50_2,IMC_51_2)


iquitos <- data.frame(matrix(NA,dim(IMC67)[1],3))
colnames(iquitos) <- c("Time","Ne","Cluster")
for(i in 1:dim(iquitos)[1]){iquitos[i,1] <- mean(c(IMC67[i,1],IMC_12[i,1],IMC_14[i,1],IMC_20[i,1],IMC_47[i,1],IMC_50[i,1],IMC_51[i,1])) + 0.01}
for(i in 1:dim(iquitos)[1]){iquitos[i,2] <- mean(c(IMC67[i,2],IMC_12[i,2],IMC_14[i,2],IMC_20[i,2],IMC_47[i,2],IMC_50[i,2],IMC_51[i,2]))}
iquitos[,3] <- rep("Iquitos",dim(iquitos)[1])


###############
# load maranon
################
MO_4 <- read.table("MO_4.0.txt",header=FALSE)
MO_9 <- read.table("MO_9.0.txt",header=FALSE)
PA289 <- read.table("PA289.0.txt",header=FALSE)
PA_107 <- read.table("PA_107.0.txt",header=FALSE)
PA_120 <- read.table("PA_120.2.0.txt",header=FALSE)
PA_121 <- read.table("PA_121.0.txt",header=FALSE)
PA_137 <- read.table("PA_137.0.txt",header=FALSE)
PA_150 <- read.table("PA_150.0.txt",header=FALSE)
PA_169 <- read.table("PA_169.0.txt",header=FALSE)
PA_218 <- read.table("PA_218.0.txt",header=FALSE)
PA_51 <- read.table("PA_51.0.txt",header=FALSE)
PA_56 <- read.table("PA_56.0.txt",header=FALSE)
PA_70 <- read.table("PA_70.0.txt",header=FALSE)
PA_88 <- read.table("PA_88.0.txt",header=FALSE)

MO_4_2 <- cbind(MO_4[,1:2],rep("MO_4",dim(MO_4)[1]),rep("Maranon",dim(MO_4)[1]))
colnames(MO_4_2) <- c("Time","Ne","accession","Cluster")
MO_9_2 <- cbind(MO_9[,1:2],rep("MO_9",dim(MO_9)[1]),rep("Maranon",dim(MO_9)[1]))
colnames(MO_9_2) <- c("Time","Ne","accession","Cluster")
PA289_2 <- cbind(PA289[,1:2],rep("PA289",dim(PA289)[1]),rep("Maranon",dim(PA289)[1]))
colnames(PA289_2) <- c("Time","Ne","accession","Cluster")
PA_107_2 <- cbind(PA_107[,1:2],rep("PA_107",dim(PA_107)[1]),rep("Maranon",dim(PA_107)[1]))
colnames(PA_107_2) <- c("Time","Ne","accession","Cluster")
PA_120_2 <- cbind(PA_120[,1:2],rep("PA_120",dim(PA_120)[1]),rep("Maranon",dim(PA_120)[1]))
colnames(PA_120_2) <- c("Time","Ne","accession","Cluster")
PA_121_2 <- cbind(PA_121[,1:2],rep("PA_121",dim(PA_121)[1]),rep("Maranon",dim(PA_121)[1]))
colnames(PA_121_2) <- c("Time","Ne","accession","Cluster")
PA_137_2 <- cbind(PA_137[,1:2],rep("PA_137",dim(PA_137)[1]),rep("Maranon",dim(PA_137)[1]))
colnames(PA_137_2) <- c("Time","Ne","accession","Cluster")
PA_150_2 <- cbind(PA_150[,1:2],rep("PA_150",dim(PA_150)[1]),rep("Maranon",dim(PA_150)[1]))
colnames(PA_150_2) <- c("Time","Ne","accession","Cluster")
PA_169_2 <- cbind(PA_169[,1:2],rep("PA_169",dim(PA_169)[1]),rep("Maranon",dim(PA_169)[1]))
colnames(PA_169_2) <- c("Time","Ne","accession","Cluster")
PA_218_2 <- cbind(PA_218[,1:2],rep("PA_218",dim(PA_218)[1]),rep("Maranon",dim(PA_218)[1]))
colnames(PA_218_2) <- c("Time","Ne","accession","Cluster")
PA_51_2 <- cbind(PA_51[,1:2],rep("PA_51",dim(PA_51)[1]),rep("Maranon",dim(PA_51)[1]))
colnames(PA_51_2) <- c("Time","Ne","accession","Cluster")
PA_56_2 <- cbind(PA_56[,1:2],rep("PA_56",dim(PA_56)[1]),rep("Maranon",dim(PA_56)[1]))
colnames(PA_56_2) <- c("Time","Ne","accession","Cluster")
PA_70_2 <- cbind(PA_70[,1:2],rep("PA_70",dim(PA_70)[1]),rep("Maranon",dim(PA_70)[1]))
colnames(PA_70_2) <- c("Time","Ne","accession","Cluster")
PA_88_2 <- cbind(PA_88[,1:2],rep("PA_88",dim(PA_88)[1]),rep("Maranon",dim(PA_88)[1]))
colnames(PA_88_2) <- c("Time","Ne","accession","Cluster")
maranon_b <- rbind(MO_4_2,MO_9_2,PA289_2,PA_107_2,PA_120_2,PA_121_2,PA_137_2,PA_150_2,PA_169_2,PA_218_2,PA_51_2,PA_56_2,PA_70_2,PA_88_2)


maranon <- data.frame(matrix(NA,dim(MO_4)[1],3))
colnames(maranon) <- c("Time","Ne","Cluster")
for(i in 1:dim(maranon)[1]){maranon[i,1] <- mean(c(MO_4[i,1],MO_9[i,1],PA289[i,1],PA_107[i,1],PA_120[i,1],PA_121[i,1],PA_137[i,1],PA_150[i,1],PA_169[i,1],PA_218[i,1],PA_51[i,1],PA_56[i,1],PA_70[i,1],PA_88[i,1])) + 0.01}
for(i in 1:dim(maranon)[1]){maranon[i,2] <- mean(c(MO_4[i,2],MO_9[i,2],PA289[i,2],PA_107[i,2],PA_120[i,2],PA_121[i,2],PA_137[i,2],PA_150[i,2],PA_169[i,2],PA_218[i,2],PA_51[i,2],PA_56[i,2],PA_70[i,2],PA_88[i,2]))}
maranon[,3] <- rep("Maranon",dim(maranon)[1])


###############
# load nacional
################
AM_1_54 <- read.table("AM_1_54.0.txt",header=FALSE)
Brisas_1 <- read.table("Brisas-1.0.txt",header=FALSE)
UF273_T1 <- read.table("UF273_T1.0.txt",header=FALSE)
UF273_T2 <- read.table("UF273_T2.0.txt",header=FALSE)

AM_1_54_2 <- cbind(AM_1_54[,1:2],rep("AM_1_54",dim(AM_1_54)[1]),rep("Nacional",dim(AM_1_54)[1]))
colnames(AM_1_54_2) <- c("Time","Ne","accession","Cluster")
Brisas_1_2 <- cbind(Brisas_1[,1:2],rep("Brisas_1",dim(Brisas_1)[1]),rep("Nacional",dim(Brisas_1)[1]))
colnames(Brisas_1_2) <- c("Time","Ne","accession","Cluster")
UF273_T1_2 <- cbind(UF273_T1[,1:2],rep("UF273_T1",dim(UF273_T1)[1]),rep("Nacional",dim(UF273_T1)[1]))
colnames(UF273_T1_2) <- c("Time","Ne","accession","Cluster")
UF273_T2_2 <- cbind(UF273_T2[,1:2],rep("UF273_T2",dim(UF273_T2)[1]),rep("Nacional",dim(UF273_T2)[1]))
colnames(UF273_T2_2) <- c("Time","Ne","accession","Cluster")
nacional_b <- rbind(AM_1_54_2,Brisas_1_2,UF273_T1_2,UF273_T2_2)

nacional <- data.frame(matrix(NA,dim(AM_1_54)[1],3))
colnames(nacional) <- c("Time","Ne","Cluster")
for(i in 1:dim(nacional)[1]){nacional[i,1] <- mean(c(AM_1_54[i,1],Brisas_1[i,1],UF273_T1[i,1],UF273_T2[i,1])) + 0.01}
for(i in 1:dim(nacional)[1]){nacional[i,2] <- mean(c(AM_1_54[i,2],Brisas_1[i,2],UF273_T1[i,2],UF273_T2[i,2]))}
nacional[,3] <- rep("Nacional",dim(nacional)[1])

###############
# load nanay
################
EET_400 <- read.table("EET_400.0.txt",header=FALSE)
NA45 <- read.table("NA45.0.txt",header=FALSE)
NA702 <- read.table("NA702_.0.txt",header=FALSE)
NA_286 <- read.table("NA_286.0.txt",header=FALSE)
NA_331 <- read.table("NA_331.0.txt",header=FALSE)
NA_92 <- read.table("NA_92.0.txt",header=FALSE)
POUND_10_B <- read.table("POUND_10_B.0.txt",header=FALSE)
POUND_7_B <- read.table("POUND_7_B.0.txt",header=FALSE)
Pound_7 <- read.table("Pound_7.0.txt",header=FALSE)
SPEC_194_75 <- read.table("SPEC_194_75.0.txt",header=FALSE)

EET_400_2 <- cbind(EET_400[,1:2],rep("EET_400",dim(EET_400)[1]),rep("Nanay",dim(EET_400)[1]))
colnames(EET_400_2) <- c("Time","Ne","accession","Cluster")
NA45_2 <- cbind(NA45[,1:2],rep("NA45",dim(NA45)[1]),rep("Nanay",dim(NA45)[1]))
colnames(NA45_2) <- c("Time","Ne","accession","Cluster")
NA702_2 <- cbind(NA702[,1:2],rep("NA702",dim(NA702)[1]),rep("Nanay",dim(NA702)[1]))
colnames(NA702_2) <- c("Time","Ne","accession","Cluster")
NA_286_2 <- cbind(NA_286[,1:2],rep("NA_286",dim(NA_286)[1]),rep("Nanay",dim(NA_286)[1]))
colnames(NA_286_2) <- c("Time","Ne","accession","Cluster")
NA_331_2 <- cbind(NA_331[,1:2],rep("NA_331",dim(NA_331)[1]),rep("Nanay",dim(NA_331)[1]))
colnames(NA_331_2) <- c("Time","Ne","accession","Cluster")
NA_92_2 <- cbind(NA_92[,1:2],rep("NA_92",dim(NA_92)[1]),rep("Nanay",dim(NA_92)[1]))
colnames(NA_92_2) <- c("Time","Ne","accession","Cluster")
POUND_10_B_2 <- cbind(POUND_10_B[,1:2],rep("POUND_10_B",dim(POUND_10_B)[1]),rep("Nanay",dim(POUND_10_B)[1]))
colnames(POUND_10_B_2) <- c("Time","Ne","accession","Cluster")
POUND_7_B_2 <- cbind(POUND_7_B[,1:2],rep("POUND_7_B",dim(POUND_7_B)[1]),rep("Nanay",dim(POUND_7_B)[1]))
colnames(POUND_7_B_2) <- c("Time","Ne","accession","Cluster")
Pound_7_2 <- cbind(Pound_7[,1:2],rep("Pound_7",dim(Pound_7)[1]),rep("Nanay",dim(Pound_7)[1]))
colnames(Pound_7_2) <- c("Time","Ne","accession","Cluster")
SPEC_194_75_2 <- cbind(SPEC_194_75[,1:2],rep("SPEC_194_75",dim(SPEC_194_75)[1]),rep("Nanay",dim(SPEC_194_75)[1]))
colnames(SPEC_194_75_2) <- c("Time","Ne","accession","Cluster")
nanay_b <- rbind(EET_400_2,NA45_2,NA702_2,NA_286_2,NA_331_2,NA_92_2,POUND_10_B_2,Pound_7_2,SPEC_194_75_2)

nanay <- data.frame(matrix(NA,dim(EET_400)[1],3))
colnames(nanay) <- c("Time","Ne","Cluster")
for(i in 1:dim(nanay)[1]){nanay[i,1] <- mean(c(EET_400[i,1],NA45[i,1],NA702[i,1],NA_286[i,1],NA_331[i,1],NA_92[i,1],POUND_10_B[i,1],POUND_7_B[i,1],Pound_7[i,1],SPEC_194_75[i,1])) + 0.01}
for(i in 1:dim(nanay)[1]){nanay[i,2] <- mean(c(EET_400[i,2],NA45[i,2],NA702[i,2],NA_286[i,2],NA_331[i,2],NA_92[i,2],POUND_10_B[i,2],POUND_7_B[i,2],Pound_7[i,2],SPEC_194_75[i,2]))}
nanay[,3] <- rep("Nanay",dim(nanay)[1])


###############
# load purus
################
#CAB_71_PL3 <- read.table("CAB_71_PL3.0.txt",header=FALSE)
CAB_76_PL3 <- read.table("CAB_76_PL3.0.txt",header=FALSE)
CAB_77_PL5 <- read.table("CAB_77_PL5.0.txt",header=FALSE)
RB39PL1 <- read.table("RB39PL1.0.txt",header=FALSE)
RB_40 <- read.table("RB_40.0.txt",header=FALSE)
RB_47_PL3 <- read.table("RB_47_PL3.0.txt",header=FALSE)

CAB_76_PL3_2 <- cbind(CAB_76_PL3[,1:2],rep("CAB_76_PL3",dim(CAB_76_PL3)[1]),rep("Purus",dim(CAB_76_PL3)[1]))
colnames(CAB_76_PL3_2) <- c("Time","Ne","accession","Cluster")
CAB_77_PL5_2 <- cbind(CAB_77_PL5[,1:2],rep("CAB_77_PL5",dim(CAB_77_PL5)[1]),rep("Purus",dim(CAB_77_PL5)[1]))
colnames(CAB_77_PL5_2) <- c("Time","Ne","accession","Cluster")
RB39PL1_2 <- cbind(RB39PL1[,1:2],rep("RB39PL1",dim(RB39PL1)[1]),rep("Purus",dim(RB39PL1)[1]))
colnames(RB39PL1_2) <- c("Time","Ne","accession","Cluster")
RB_40_2 <- cbind(RB_40[,1:2],rep("RB_40",dim(RB_40)[1]),rep("Purus",dim(RB_40)[1]))
colnames(RB_40_2) <- c("Time","Ne","accession","Cluster")
RB_47_PL3_2 <- cbind(RB_47_PL3[,1:2],rep("RB_47_PL3",dim(RB_47_PL3)[1]),rep("Purus",dim(RB_47_PL3)[1]))
colnames(RB_47_PL3_2) <- c("Time","Ne","accession","Cluster")
purus_b <- rbind(CAB_76_PL3_2,CAB_77_PL5_2,RB39PL1_2,RB_40_2,RB_47_PL3_2)

purus <- data.frame(matrix(NA,dim(CAB_76_PL3)[1],3))
colnames(purus) <- c("Time","Ne","Cluster")
for(i in 1:dim(purus)[1]){purus[i,1] <- mean(c(CAB_76_PL3[i,1],CAB_77_PL5[i,1],RB39PL1[i,1],RB_40[i,1],RB_47_PL3[i,1])) + 0.01}
for(i in 1:dim(purus)[1]){purus[i,2] <- mean(c(CAB_76_PL3[i,2],CAB_77_PL5[i,2],RB39PL1[i,2],RB_40[i,2],RB_47_PL3[i,2]))}
purus[,3] <- rep("Purus",dim(purus)[1])








####################
## All combined
####################


B <- rbind(Criollo,curaray,nanay,contamana,amelonado,maranon,nacional,guianna,iquitos,purus)

mypal2 <- c(mypal[5],mypal[4],mypal[1],mypal[2],mypal[8],mypal[9],mypal[6],mypal[7],mypal[3],mypal[10])

ggplot(data = B, aes(x = log(Time,base=10), y = Ne)) + geom_line(aes(colour=Cluster),size=1) + scale_colour_manual(values=mypal2) + ylab(paste("Effective population size Ne*",expression(10^4),sep="")) + guides(fill=FALSE) + ylim(0,40) + scale_x_continuous("Time (years)",breaks = c(2,3,4,5), labels=c(expression(1^2),expression(1^3),expression(1^4),expression(1^5)), limits=c(2,5)) + theme(axis.text.x= element_text(size = rel(1.1), hjust = 1, vjust = 1))

p <- ggplot(data = B, aes(x = log(Time,base=10), y = log(Ne,base=2))) + geom_line(aes(colour=Cluster),size=1) + scale_colour_manual(values=mypal2) + theme(axis.text.x= element_text(size = rel(1.1), angle=45, hjust = 1, vjust = 1)) + ylab("Effective population size") + guides(fill=FALSE)  + xlim(3,6) + xlab("Time (log_10 scale)")


p + scale_x_discrete("",breaks = c(2,3,4,5,6), labels=c(expression(1^2),expression(1^3),expression(1^4),expression(1^5),expression(1^6)))


####################
# Combining the raw for multifacet plot. Not presented in the paper, but another possibility for the plotting
####################

B2 <- rbind(Criollo_b,curaray_b,nanay_b,contamana_b,amelonado_b,maranon_b,nacional_b,guianna_b,iquitos_b,purus_b)


ggplot(data = Matina, aes(x = log(Time,base=10), y = Ne)) + geom_line(size=1.5,colour=mypal2[2]) + ylab(paste("Effective population size Ne*",expression(10^4),sep="")) + guides(fill=FALSE) + ylim(0,0.5) + scale_x_continuous("Time (years)",breaks = c(2,3,4,5), labels=c(expression(1^2),expression(1^3),expression(1^4),expression(1^5)), limits=c(2,5)) + theme(axis.text.x= element_text(size = rel(1.1), hjust = 1, vjust = 1))

pdf("plot_psmc",width=7.096491,height=3.701754)

dev.off()
