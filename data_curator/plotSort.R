setwd("/projectnb2/bf528/users/dachshund/project_4/project-4-dachsund/data_curator/")

merge_output_3_files <- read.csv("merge_output_3_files.csv", header=FALSE)

m <- merge_output_3_files

#t <- subset(m, nchar(V1) > 33)

plot(ecdf(log(m$V2)), main = "Cumulative Distribution Plot for scRNA-seq Barcodes", xlab = "Cumulative Barcodes", col = "violet", ylab="Frequency")
hist(log(m$V2), main="Barcode Frequencies", xlab = "log(Number of Barcodes)", col = "purple")
summary(m$V2)

filter <- subset(m, V2 > 396)
plot(ecdf(filter$V2), main = "Cumulative Distribution Plot for scRNA-seq Barcodes After Filtering", col = "maroon",ylab="Frequency")
hist(filter$V2, main="Barcode Frequencies", xlab = "Number of Barcodes", col ="red")

summary(filter$V2)

write.csv(filter, "whitelist.csv", row.names = F)

wl <- read.csv("whitelist.csv")

colnames(wl) <- NULL

write.csv(wl[,1], "wl.csv", row.names = F, quote= F)
