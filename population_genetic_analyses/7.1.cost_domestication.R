library(ggplot2)
library(RColorBrewer)
librbary(reshape)
library(plyr)
setwd("/path_to_fitness_analyses/admixture_supervised/")

mypallete <- brewer.pal(10,"Spectral")

data <- data.frame(read.table("cacao2.10.Q",header=FALSE))
ids <- data.frame(read.table("ids2",header=FALSE))
rownames(data) <- ids[,1]
known_list <- c("criollo","sp1","sp3","sp9","LCTEEN_141","CUR3_G37_A6","CUR3_G38_A8","SIL-1_G56_A6","CUR3_G39_A10","NA702","NA_286","NA_331","NA_92","SPEC_194_75","POUND_10_B","POUND_7_B","EET_400","NA45","Pound_7","PMF_20","PMF_27","SCA_11","SCA_19","SCA_24.2","SCA_5","T695_SCA-6_A1","NH_40","NH_53","Catongo","Matina","SIAL169","mvP30","SIC806","SIAL70","SIAL84","MATINA_Tica2","REDAMEL_1_31","REDAMEL_1_27","TRD86","PA_120","PA_121","PA_137","PA_56","PA_70","PA289","MO_9","MO_4","PA_169","PA_107","PA_51","PA_150","PA_88","PA_218","Brisas-1","UF273_T2","AM_1_54","UF273_T1","ELP_20_A","GU255_V","GU_114_P","GU_175_P","GU_219_F","GU_222","GU_300_P","GU_307","GU_308A","IMC_14","IMC67","IMC_20","IMC_50","IMC_12","IMC_47","IMC_51","CAB_76_PL3","CAB_77_PL5","RB39PL1","RB_40","RB_47_PL3","CAB_71_PL3")
data2 <- subset(data,!(rownames(data) %in% known_list))
colnames(data2) <- c("K1","K2","K3","K4","K5","K6","K7","K8","K9","K10")

pheno <- data.frame(read.table("../../TableBLUPs.txt",header=TRUE))
clean_ids <- data.frame(read.table("clean_ids",header=FALSE))
data2.1 <- cbind(rownames(data2),data2)
colnames(data2.1) <- c("Clone","K1","K2","K3","Amelonado_ancestry","K5","K6","K7","K8","K9","Criollo_ancestry")

data_pheno <- join(data2.1,pheno,by="Clone",type="right")

B <- melt(data_pheno)

glm(Yield.kg.ha.year. ~ Criollo_ancestry + F,data=B)
anova(glm(Yield.kg.ha.year. ~ Criollo_ancestry,data=b),glm(Yield.kg.ha.year. ~ Criollo_ancestry + F,data=b),test="Chisqr")

summary(glm(Yield.kg.ha.year. ~ Criollo_ancestry*log(KB + 1),data=B))

glm(Yield.kg.ha.year. ~ Amelonado_ancestry + F,data=B)
anova(glm(Yield.kg.ha.year. ~ Amelonado_ancestry,data=b),glm(Yield.kg.ha.year. ~ Amelonado_ancestry + F,data=b),test="Chisqr")


############
# some additional fiddling with alternative ways to fit the data was attempted and also we ran diagnostics on the fit. IF interested in those, please write to the authors
###########

