data1 <- data1[order(data1$weight),]
data2 <- data2[order(data2$weight),]
data3 <- data3[order(data3$weight),]
data4 <- data4[order(data4$weight),]

data1 <- data1[data1[,3]>0,]
data2 <- data2[data2[,3]>0,]
data3 <- data3[data3[,3]>0,]
data4 <- data4[data4[,3]>0,]

dir.create("all_data_dyngenie")

write.table(data1, file = "all_data_dyngenie/condition_1_dyngenie.txt", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(data2, file = "all_data_dyngenie/condition_2_dyngenie.txt", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(data3, file = "all_data_dyngenie/condition_3_dyngenie.txt", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(data4, file = "all_data_dyngenie/condition_4_dyngenie.txt", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)



setwd("/home/jconnell/projects/niab/gene_regulatory_network/ARACNe-AP")

data1.1 <- read.table("C1_results/network.txt", header = TRUE)
data2.1 <- read.table("C2_results/con2_arac_network.txt", header = TRUE)
data3.1 <- read.table("C3_results/Time_series_1_7_con3_arac.txt", header = TRUE)
data4.1 <- read.table("C4_results/Time_series_1_7_con4_arac.txt", header = TRUE)

data1.1 <- data1.1[,c(1,2,4)]
data2.1 <- data2.1[,c(1,2,4)]
data3.1 <- data3.1[,c(1,2,4)]
data4.1 <- data4.1[,c(1,2,4)]

colnames(data1.1) <- c("regulatory.gene", "target.gene", "weight")
colnames(data2.1) <- c("regulatory.gene", "target.gene", "weight")
colnames(data3.1) <- c("regulatory.gene", "target.gene", "weight")
colnames(data4.1) <- c("regulatory.gene", "target.gene", "weight")

data1.1 <- data1.1[order(data1.1$weight),]
data2.1 <- data2.1[order(data1.1$weight),]
data3.1 <- data3.1[order(data1.1$weight),]
data4.1 <- data4.1[order(data1.1$weight),]

data1.1 <- data1.1[data1.1[,3]>0,]
data2.1 <- data2.1[data2.1[,3]>0,]
data3.1 <- data3.1[data3.1[,3]>0,]
data4.1 <- data4.1[data4.1[,3]>0,]

dir.create("all_data_aracne")

write.table(data1.1, file = "all_data_aracne/condition_1_aracne.txt", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(data2.1, file = "all_data_aracne/condition_2_aracne.txt", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(data3.1, file = "all_data_aracne/condition_3_aracne.txt", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(data4.1, file = "all_data_aracne/condition_4_aracne.txt", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)



