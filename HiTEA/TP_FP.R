if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges", force = TRUE)
BiocManager::install("genomation")

#####################################################################

library(GenomicRanges)
library(genomation)

#####################################################################

ref <- read.csv(file = "ref_hg38.txt", sep = '\t', header = FALSE)
non_ref <- read.csv(file = "non_ref_hg38.txt", sep = '\t', header = FALSE)
hitea <- read.table(file = "L1_refined.bed", sep = '\t', header = FALSE)

write.table(ref, file = "ref.bed", row.names=FALSE, sep="\t", quote = FALSE)
write.table(non_ref, file = "non_ref.bed", row.names=FALSE, sep="\t", quote = FALSE)

ref_bed <- with(ref, GRanges(V1, IRanges(start = V2, end = V3)));
non_ref_bed <- with(non_ref, GRanges(V1, IRanges(start = V2, end = V3)));
hitea_bed <- with(hitea, GRanges(V1, IRanges(start = V2, end = V3)));


# reference 

# TPR

hits <- findOverlaps(ref_bed, hitea_bed) 
print(hits)

L1_window_overlap <- length(unique(hits@from)) 
print(L1_window_overlap)
total_L1 <- length(ref_bed@ranges@start)
print(total_L1)
TPR <- (L1_window_overlap/total_L1)*100
print(TPR)

# FDR

total_window <- length(hitea_bed@ranges@start)
print(total_window)

hits <- findOverlaps(ref_bed, hitea_bed) 
print(hits)
window_overlap <- length(unique(hits@to))
print(window_overlap)
window_no_overlap <- total_window - window_overlap 

FDR = (window_no_overlap/total_window)*100
print(FDR)


# non reference 

# TPR

hits <- findOverlaps(non_ref_bed, hitea_bed) 
print(hits)

L1_window_overlap <- length(unique(hits@from)) 
print(L1_window_overlap)
total_L1 <- length(non_ref_bed@ranges@start)
print(total_L1)
TPR <- (L1_window_overlap/total_L1)*100
print(TPR)

# FDR

total_window <- length(hitea_bed@ranges@start)
print(total_window)

hits <- findOverlaps(non_ref_bed, hitea_bed) 
print(hits)
window_overlap <- length(unique(hits@to))
print(window_overlap)
window_no_overlap <- total_window - window_overlap 

FDR = (window_no_overlap/total_window)*100
print(FDR)
