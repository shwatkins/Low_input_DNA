library(ggplot2)
library(viridis)
library(ggpubr)

load("betas.Robj")
# loads object norm.beta, cols = samples, rows = cpgs
slides <- colnames(norm.beta)

# separate out meth at 40, 200 and 400ng dilutions
slides_40 <- slides[c(13:16)]
slides_200 <- slides[c(5:8)]
slides_400 <- slides[c(21:24)]

norm.beta <- as.data.frame(norm.beta)

# create separate dataframes for each dilution
meth_40 <- norm.beta[,slides_40]
meth_200 <- norm.beta[,slides_200]
meth_400 <- norm.beta[,slides_400]

# take the mean methylation for each cpg across the 4 pooled samples
# make it an extra column
meth_40$mean_40 <- rowMeans(meth_40)
meth_200$mean_200 <- rowMeans(meth_200)
meth_400$mean_400 <- rowMeans(meth_400)

# Make a df of the means
meth <- merge(meth_40, meth_200, by="row.names")
meth <- merge(meth, meth_400, by.x="Row.names", by.y="row.names")
rownames(meth) <- meth$Row.names
meth_full <- meth
# make df of means only:
meth <- meth[,c(6,11,16)]

library(reshape2)

df_melt <- melt(meth)
# plot density of the mean betas for old pooled for each dilution
jpeg(filename = "pooled_betas_density_plot.jpg", width = 4, height = 3, units = 'in', res=600)
ggplot(df_melt, aes(value, fill = variable, colour = variable)) +
  geom_density(alpha = 0.1) +
  theme_bw()+
  labs(x="DNAm beta value")+
  scale_colour_viridis(discrete = T, name = "Total input\nDNA",labels=c("40ng","200ng","400ng"))+
  scale_fill_viridis(discrete = T, name = "Total input\nDNA",labels=c("40ng","200ng","400ng"))
dev.off()

# Now we look at how the dilutions differ from each other
# To do this we split them in to increments of 5% methylation
n.blocks <- 20
block.size <- ceiling(nrow(meth)/20)
# assign the rows to a block
blocks <- rep(1:n.blocks, each=block.size, len=nrow(meth))
blocks <- split(1:nrow(meth), blocks)
# which row does the block start at?
block.start <- sapply(blocks, head, n=1)
# Size of each block:
block.size <- sapply(blocks, length)

## or do it in 5% increments
j <- seq(from=0,by=0.05, length.out = 20)
k <- seq(to=1,by=0.05, length.out = 20)

meth$block <- NA
meth$block_mean <- NA

# Now we get the mean and SD of cpgs that are in each 5% block of 
# methylation (based on the methylation level of each cpg at 400ng)
mean_sd <- data.frame(matrix(ncol = 7, nrow=0))
names(mean_sd) <- c("block_temp", "mean_40", "sd_40", "mean_200", "sd_200", "mean_400", "sd_400")
for(i in 1:length(blocks)){
  meth_temp <- meth[meth$mean_400>j[i] & meth$mean_400<=k[i],]
  meth_temp
  mean_meth <- mean(meth_temp$mean_400)
  print(mean_meth)
  meth$block_mean[meth$mean_400>j[i] & meth$mean_400<=k[i]] <- mean_meth
  meth$block[meth$mean_400>j[i] & meth$mean_400<=k[i]] <- j[i]
  
  mean_40 <- mean(meth_temp$mean_40)
  sd_40 <- sd(meth_temp$mean_40)
  mean_200 <- mean(meth_temp$mean_200)
  sd_200 <- sd(meth_temp$mean_200)
  mean_400 <- mean(meth_temp$mean_400)
  sd_400 <- sd(meth_temp$mean_400)
  block_temp <- j[i]
  mean_sd_temp <- data.frame(block_temp, mean_40, sd_40, mean_200,
                             sd_200, mean_400, sd_400)
  mean_sd <- rbind(mean_sd,mean_sd_temp)
  
}
head(meth)
head(mean_sd)
write.csv(mean_sd, file = "sd_200_40.csv", quote = F, row.names = F)


# Use Levene's test to test variance differences for each block

library(reshape2)
library(car)

var_diff_40 <- list()
for(i in 1:length(blocks)){
  print(i)
  dat <- meth[meth$block == j[i],]
  print(dim(dat))
  dat <- dat[,c(1,3,4,5)]
  dat_melt <- melt(data = dat, measure.vars=c("mean_40","mean_400"))
  names(dat_melt)[names(dat_melt) == 'variable'] <- 'dilution'
  levene_test_result <- leveneTest(value ~ dilution, data = dat_melt)
  var_diff_40[[i]] <- levene_test_result
}

var_diff_200 <- list()
for(i in 1:length(blocks)){
  print(i)
  dat <- meth[meth$block == j[i],]
  print(dim(dat))
  dat <- dat[,c(2,3,4,5)]
  dat_melt <- melt(data = dat, measure.vars=c("mean_200","mean_400"))
  names(dat_melt)[names(dat_melt) == 'variable'] <- 'dilution'
  levene_test_result <- leveneTest(value ~ dilution, data = dat_melt)
  var_diff_200[[i]] <- levene_test_result
}

names(var_diff_40) <- j
names(var_diff_200) <- j

###

# Create boxplots comparing the methylation of cpgs in the 5% 
# methylation bands with the methylation at 400ng.

## 1. compare 200ng with 400ng 
meth_400_200 <- meth[,c(2:5)]
head(meth_400_200)
meth_400_200_melt <- melt(data = meth_400_200, measure.vars=c("mean_400","mean_200"))
head(meth_400_200_melt)
names(meth_400_200_melt)[names(meth_400_200_melt) == 'variable'] <- 'dilution'
names(meth_400_200_melt)[names(meth_400_200_melt) == 'value'] <- 'mean_meth'

jpeg(filename ="pooled_meth_boxplot_200_400.jpg", width = 5, height = 3, units = "in", res=600)
ggplot(meth_400_200_melt, aes(x=factor(block), y=mean_meth, colour=dilution)) +
  geom_boxplot(outlier.shape = 1,outlier.size = 0.5, outlier.alpha = 0.2)+#(aes(color=dilution))+
  scale_colour_viridis(discrete = T, begin = 0, end=0.75,labels=c("400ng","200ng"),name = "Total input\nDNA")+
  theme(axis.text.x=element_blank())+
  theme_bw()+
  scale_x_discrete(breaks=factor(k)[seq(0,length(k),by=2)])+  
  labs(title="Mean methylation at 200ng and 400ng", x = "Methylation band (400ng)", y = "Mean methylation")
dev.off()

## 1. compare 40ng with 400ng 
meth_400_40 <- meth[,c(1,3,4,5)]
head(meth_400_40)
meth_400_40_melt <- melt(data = meth_400_40, measure.vars=c("mean_400","mean_40"))
head(meth_400_40_melt)
names(meth_400_40_melt)[names(meth_400_40_melt) == 'variable'] <- 'dilution'
names(meth_400_40_melt)[names(meth_400_40_melt) == 'value'] <- 'mean_meth'

jpeg(filename ="pooled_meth_boxplot_40_400.jpg", width = 5, height = 3, units = "in", res=600)
ggplot(meth_400_40_melt, aes(x=factor(block), y=mean_meth, colour=dilution)) +
  geom_boxplot(outlier.shape = 1,outlier.size = 0.5, outlier.alpha = 0.2)+#(aes(color=dilution))+
  scale_colour_viridis(discrete = T, begin = 0, end=0.75,labels=c("400ng","40ng"),name = "Total input\nDNA")+
  theme(axis.text.x=element_blank())+
  theme_bw()+
  scale_x_discrete(breaks=factor(k)[seq(0,length(k),by=2)])+  
  labs(title="Mean methylation at 40ng and 400ng", x = "Methylation band (400ng)", y = "Mean methylation")
dev.off()

####

## Assessing noise for each DNA input level
# we do this separately for each DNA input level

## 1. 40ng 

# take the mean of replicates 1 and 2, and 3 and 4
meth_40$mean_1 <- (meth_40[,1]+meth_40[,2])/2
meth_40$mean_2 <- (meth_40[,3]+meth_40[,4])/2
meth_40$block_mean <- NA
meth_40$block <- NA

# Get mean methylation for replicates 1 and 2, and 3 and 4, for each
# of the 20 blocks
for(i in 1:length(blocks)){
  meth_temp <- meth_40[meth_40$mean_1>j[i] & meth_40$mean_1<=k[i],]
  mean_meth <- mean(meth_temp$mean_1)
  print(mean_meth)
  meth_40$block_mean[meth_40$mean_1>j[i] & meth_40$mean_1<=k[i]] <- mean_meth
  meth_40$block[meth_40$mean_1>j[i] & meth_40$mean_1<=k[i]] <- j[i]
}

noise_40 <- ggplot(meth_40, aes(x=factor(block), y=mean_2)) +
  geom_boxplot(outlier.shape = 1,outlier.size = 0.5, outlier.alpha = 0.2,color="#440154FF", fill="#440154FF", alpha=0.2)+
  theme(axis.text.x=element_blank())+
  theme_bw()+
  scale_x_discrete(breaks=factor(k)[seq(0,length(k),by=2)])+  
  labs(title="Sample noise: 40ng total input DNA", x = "Methylation band", y = "Mean methylation")

## 2. 200ng

# take the mean of replicates 1 and 2, and 3 and 4
meth_200$mean_1 <- (meth_200[,1]+meth_200[,2])/2
meth_200$mean_2 <- (meth_200[,3]+meth_200[,4])/2
head(meth_200)
meth_200$block_mean <- NA
meth_200$block <- NA

# Get mean methylation for replicates 1 and 2, and 3 and 4, for each
# of the 20 blocks
for(i in 1:length(blocks)){
  meth_temp <- meth_200[meth_200$mean_1>j[i] & meth_200$mean_1<=k[i],]
  mean_meth <- mean(meth_temp$mean_1)
  print(mean_meth)
  meth_200$block_mean[meth_200$mean_1>j[i] & meth_200$mean_1<=k[i]] <- mean_meth
  meth_200$block[meth_200$mean_1>j[i] & meth_200$mean_1<=k[i]] <- j[i]
}
head(meth_200)

noise_200 <- ggplot(meth_200, aes(x=factor(block), y=mean_2)) +
  geom_boxplot(outlier.shape = 1,outlier.size = 0.5, outlier.alpha = 0.2,color="#238A8DFF", fill="#238A8DFF", alpha=0.2)+
  theme(axis.text.x=element_blank())+
  theme_bw()+
  scale_x_discrete(breaks=factor(k)[seq(0,length(k),by=2)])+  
  labs(title="Sample noise: 200ng total input DNA", x = "Methylation band", y = "Mean methylation")

### 3. 400ng

# take the mean of replicates 1 and 2, and 3 and 4
meth_400$mean_1 <- (meth_400[,1]+meth_400[,2])/2
meth_400$mean_2 <- (meth_400[,3]+meth_400[,4])/2
head(meth_400)
meth_400$block_mean <- NA
meth_400$block <- NA

# Get mean methylation for replicates 1 and 2, and 3 and 4, for each
# of the 20 blocks
for(i in 1:length(blocks)){
  meth_temp <- meth_400[meth_400$mean_1>j[i] & meth_400$mean_1<=k[i],]
  mean_meth <- mean(meth_temp$mean_1)
  print(mean_meth)
  meth_400$block_mean[meth_400$mean_1>j[i] & meth_400$mean_1<=k[i]] <- mean_meth
  meth_400$block[meth_400$mean_1>j[i] & meth_400$mean_1<=k[i]] <- j[i]
}
head(meth_400)

noise_400 <- ggplot(meth_400, aes(x=factor(block), y=mean_2)) +
  geom_boxplot(outlier.shape = 1,outlier.size = 0.5, outlier.alpha = 0.2,color="#73D055FF", fill="#73D055FF", alpha=0.2)+
  theme(axis.text.x=element_blank())+
  theme_bw()+
  scale_x_discrete(breaks=factor(k)[seq(0,length(k),by=2)])+  
  labs(title="Sample noise: 400ng total input DNA", x = "Methylation band", y = "Mean methylation")

# Now create jpeg with all 3 plots
jpeg(filename = "sample_noise_plots.jpg", width = 8, height = 6, units = 'in', res=600)
ggarrange(noise_40, noise_200, noise_400, 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
dev.off()

#### 

## now compare the noise between dilutions using Levene's test

var_diff_200 <- list()

dat_200 <- meth_200[,c("mean_2","block")]
names(dat_200) <- c("mean_200","block")
dat_400 <- meth_400[,c("mean_2","block")]
names(dat_400) <- c("mean_400","block")

for(i in 1:length(blocks)){
  print(i)
  dat_200_temp <- dat_200[dat_200$block == j[i],]
  colnames(dat_200_temp) <- c("mean","dilution")
  dat_200_temp$dilution <- c("mean_200")
  dat_400_temp <- dat_400[dat_400$block == j[i],]
  colnames(dat_400_temp) <- c("mean","dilution")
  dat_400_temp$dilution <- c("mean_400")
  dat <- rbind(dat_200_temp,dat_400_temp)
  print(dim(dat))
  levene_test_result <- leveneTest(mean ~ dilution, data = dat)
  var_diff_200[[i]] <- levene_test_result
}

names(var_diff_200) <- j

var_diff_40 <- list()

dat_40 <- meth_40[,c("mean_2","block")]
names(dat_40) <- c("mean_40","block")
dat_400 <- meth_400[,c("mean_2","block")]
names(dat_400) <- c("mean_400","block")

for(i in 1:length(blocks)){
  print(i)
  dat_40_temp <- dat_40[dat_40$block == j[i],]
  colnames(dat_40_temp) <- c("mean","dilution")
  dat_40_temp$dilution <- c("mean_40")
  dat_400_temp <- dat_400[dat_400$block == j[i],]
  colnames(dat_400_temp) <- c("mean","dilution")
  dat_400_temp$dilution <- c("mean_400")
  dat <- rbind(dat_40_temp,dat_400_temp)
  print(dim(dat))
  levene_test_result <- leveneTest(mean ~ dilution, data = dat)
  var_diff_40[[i]] <- levene_test_result
}

names(var_diff_40) <- j

summary_40 <- do.call(rbind.data.frame, var_diff_40)
summary_200 <- do.call(rbind.data.frame, var_diff_200)

var_diff_40_200 <- list()

for(i in 1:length(blocks)){
  print(i)
  dat_40_temp <- dat_40[dat_40$block == j[i],]
  colnames(dat_40_temp) <- c("mean","dilution")
  dat_40_temp$dilution <- c("mean_40")
  dat_200_temp <- dat_200[dat_200$block == j[i],]
  colnames(dat_200_temp) <- c("mean","dilution")
  dat_200_temp$dilution <- c("mean_200")
  dat <- rbind(dat_40_temp,dat_200_temp)
  print(dim(dat))
  levene_test_result <- leveneTest(mean ~ dilution, data = dat)
  var_diff_40_200[[i]] <- levene_test_result
}
names(var_diff_40_200) <- j

summary_40_200 <- do.call(rbind.data.frame, var_diff_40_200)

save(summary_40,summary_200,summary_40_200, file="table_noise_variance_test.Rdata")

write.csv(summary_40, file="noise_summary_40.csv",quote = F)
write.csv(summary_200, file="noise_summary_200.csv",quote = F)
write.csv(summary_40_200, file="noise_summary_40_200.csv",quote = F)

