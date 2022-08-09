################################ LOAD LIBRARIES ################################

library(ggplot2)
library(GenomicRanges)

################################ WINDOWS #######################################

# read the windows table 
windows <- read.csv(file = "window.tsv", sep = '\t', header = TRUE)
# remove any rows that are unknown chr 
windows <- windows[!nchar(as.character(windows$chr)) > 5, ]

################################ TP DATA #######################################


# opening the L1PA data 
L1PA_bed <- read.table("Repeatmasker_hg38_4.0.5_L1PA.bed", sep = '\t', header = FALSE)

# extracting the L1PA1-3
sub_L1PA_bed <- L1PA_bed[L1PA_bed$V4 == c("L1HS", "L1PA2", "L1PA3"), ]###

# full length L1PAs only (5kb<)
fl_sub_L1PA_bed <- subset(sub_L1PA_bed, V3 - V2 >= 5000)

# converting to a genomic ranges object 
gr_fl_sub_L1PA_bed <- with(fl_sub_L1PA_bed, GRanges(V1, IRanges(start = V2, end = V3)))

################################# FP DATA ######################################

# opening the L1PA data
L1PA_bed <- read.table("Repeatmasker_hg38_4.0.5_L1PA.bed", sep = '\t', header = FALSE)

# extracting only the full length L1PAs (>5kb)
fl_L1PA_bed <- subset(L1PA_bed, V3 - V2 >= 5000)

# increasing the L1s region by +-50kb
fl_L1PA_bed[, 2] <- fl_L1PA_bed[, 2] - 50000
fl_L1PA_bed[, 3] <- fl_L1PA_bed[, 3] + 50000

# making an empty dataframe to append to 
FP_df <- data.frame(chr = as.character(), start = as.numeric(), end = as.numeric())

# list of all chr 
chr_list = list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

# looping over the list of chr
for (chr in chr_list) {
  
  # subsetting the L1PA df according to chr 
  chr_sub <- subset(fl_L1PA_bed, V1=chr)
  
  # extracting out the chr, start and end column
  chr_sub <- chr_sub[ , c("V1", "V2", "V3")]
  
  # function to shift the start column up by 1
  shift <- function(x, n){
    c(x[-(seq(n))], rep(NA, n))
  }
  
  # shifting up the column
  chr_sub$V2 <- shift(chr_sub$V2, 1)
  
  # removing the final line of the df 
  chr_sub <- head(chr_sub, -1)
  
  # rearrange the columns after the shifts 
  chr_sub <- chr_sub[ , c("V1", "V3", "V2")]
  
  # removing rows that have a start position greater than end (occurs due to 50kb increase)
  chr_sub <- chr_sub[!rowSums(chr_sub["V3"] >= chr_sub["V2"]), ]
  
  # binding the data frames together 
  FP_df <- rbind(FP_df, chr_sub)
  
}

# producing a file for the +-50kb of L1PA
write.table(fl_L1PA_bed, file = "L1PA_50kb.bed", row.names=FALSE, sep="\t", quote = FALSE);

# converting to a genomic ranges object 
gr_FP_df <- with(FP_df, GRanges(V1, IRanges(start = V3, end = V2)))


############

hist(windows$var, breaks = 100)

less_200 <- windows[windows$var <= 60, ]
hist(less_200$var, breaks = 100)

################################################################################
################################ looping it  ###################################
################################################################################

# create a list of threshold values from 0-1 with increments of 0.1
thres_val <- seq(from = 0, to = 70, by = 1)

# produce a df to append to 
df <- data.frame(var = as.numeric(), TPR = as.numeric(), FPR = as.numeric())

# loop through the threshold values 
for (val in thres_val) {
  
  ## TPR calculation 
  # subset the df according to the threshold value 
  thres <- windows[windows$var >= val, ]
  
  # converting to a genomics regions variable 
  gr_thres <- with(thres, GRanges(chr, IRanges(start = start_hg38, end = end_hg38)))
  
  # overlapping the subset L1PA and the threshold df 
  hits <- findOverlaps(gr_fl_sub_L1PA_bed, gr_thres) 
  
  # how many L1 had overlapped a window 
  L1_window_overlap <- length(unique(hits@from)) 
  
  # total number of L1PAs that were subset 
  total_L1 <- nrow(fl_sub_L1PA_bed)  
  
  # calculating the TPR
  TPR <- (L1_window_overlap/total_L1)*100
  
  
  ## FP
  # subsetting the windows df according to variance threshold value 
  thres <- windows[windows$var >= val, ]
  
  # converting the df to a genomic ranges object 
  gr_thres <- with(thres, GRanges(chr, IRanges(start = start_hg38, end = end_hg38)))
  
  # get all the windows that overlap with the L1PA
  hits <- findOverlaps(gr_FP_df, gr_thres)
  
  # index of the windows that do overlap - unique
  non_L1_over <- length(unique(hits@from))
  
  # the total non L1PA regions
  total_non_L1 <- nrow(FP_df)
  
  # calculating the FPR 
  FPR <- (non_L1_over/total_non_L1)*100
  
  
  # adding the variance value, TPR and FDR to the df
  df[nrow(df) + 1, ] = c(val, TPR, FPR)
  
}


### plotting 

plot(df$FPR, df$TPR, xlim = c(0,100), ylim = c(0,100))
text(df$FPR+6, df$TPR, labels=df$var)
title("L1PA1-3")

library(ggplot2)
library(ggrepel)
ggplot(df, aes(x = FPR, y = TPR)) +
  geom_text_repel(col="black", aes(label = var)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "gam", formula = y ~ s(log(x)), color = "#3F93CB") +
  labs(title = "Real data, AUC = 0.715") +
  geom_abline(intercept=0, slope=1, linetype="dashed", color = "black")



### area under the curve 

auc = function(fpr, tpr) {
  fpr.diff = fpr[1:(length(fpr)-1)] - fpr[2:length(fpr)]
  tpr.mean = (tpr[1:(length(tpr)-1)] + tpr[2:length(tpr)])/2
  total.area = sum(fpr.diff*tpr.mean)
  return(total.area/10000)
}

AUC <- auc(df$FPR, df$TPR)
print(AUC)
