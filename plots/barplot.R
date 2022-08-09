################################################################################
################# with SVA_F included in the Repeatmasker ######################
################################################################################

# subsetting the windows df according to threshold 3 
thres_40 <- windows[windows$var >= 23, ]

# removing any additional columns 
thres_40 <- thres_40[-c(2)]

write.table(thres_40, file = "thres_23.bed", row.names=FALSE, sep="\t", quote = FALSE)

### UNIX code - bedtools interesect 
# bedtools intersect -b input_bedtools/thres_3.bed -a input/Repeatmasker_hg38_4.0.5_L1PA.bed -c > results_bedtools/L1PA_thres3.bed  
# bedtools intersect -b input_bedtools/thres_3.bed -a input/Repeatmasker_hg38_4.0.5_SVA_F.bed -c > results_bedtools/SVA_F_thres3.bed
###

# reading bedtools intersect results 
# L1PA
L1PA_thres3 <- read.table("L1PA_no_blen_thres32.bed", sep = '\t', header = FALSE)
# SVA_F
SVA_F_thres3 <- read.table("SVA_F_no_blen_thres32.bed", sep = '\t', header = FALSE)
# joining the two files together 
L1PA_SVA_F_thres3 <- rbind(L1PA_thres3, SVA_F_thres3)

# subsetting the fl 
# L1PA fl are > 5000
fl_L1PA_thres3 <- subset(L1PA_thres3, V3 - V2 >= 5000)
# SVA_F fl are > 1300
fl_SVA_F_thres3 <- subset(SVA_F_thres3, V3 - V2 >= 1300)
# joining the two files together 
fl_L1PA_SVA_F_thres3 <- rbind(fl_L1PA_thres3, fl_SVA_F_thres3)

# making a list of all the TEs in our file
L1PA_ls <- list("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6", "L1PA7", "L1PA8", "L1PA8A", "L1PA10", "L1PA11", "L1PA12", "L1PA13", "L1PA14", "L1PA15", "L1PA15-16", "L1PA16", "L1PA17", "SVA_F")

# making an empty dataframe to append to 
df_L1PA_thres3 <- data.frame(L1PA = as.character(), L1PA_perc = as.numeric(), all_or_fl = as.character(), overlap = as.numeric(), total_TE = as.numeric())

# looping over the list of L1PAs
for (L1PA in L1PA_ls) {
  
  ## all
  
  # subset according to the L1PA
  sub_L1PA <- subset(L1PA_SVA_F_thres3, V4 == L1PA)
  # total amount of particular L1PA
  total_L1PA <- nrow(sub_L1PA)
  
  # subset where there is an overlap with a window 
  L1PA_reg <- subset(sub_L1PA, V7 != 0)
  # total L1PA that overlapped with a window 
  total_over <- nrow(L1PA_reg)
  
  # calculate the percentage 
  perc_all <- ((total_over/total_L1PA)*100)
  
  # appending the results to the dataframe
  df_L1PA_thres3[nrow(df_L1PA_thres3) + 1, ] = c(L1PA, perc_all, "all", total_over, total_L1PA)
  
  
  ## fl L1PA and SVA_F
  
  # subset according to the L1PA
  sub_L1PA <- subset(fl_L1PA_SVA_F_thres3, V4 == L1PA)
  # total amount of particular L1PA
  total_L1PA <- nrow(sub_L1PA)
  
  # subset where there is an overlap with a window 
  L1PA_reg <- subset(sub_L1PA, V7 != 0)
  # total L1PA that overlapped with a window 
  total_over <- nrow(L1PA_reg)
  
  # calculate the percentage 
  perc_fl <- ((total_over/total_L1PA)*100)
  
  # appending the results to the dataframe  
  df_L1PA_thres3[nrow(df_L1PA_thres3) + 1, ] = c(L1PA, perc_fl, "fl", total_over, total_L1PA)
  
}

# fixing the order of the TEs to be shown on the plot 
df_L1PA_thres3$L1PA = factor(df_L1PA_thres3$L1PA, levels = c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6", "L1PA7", "L1PA8", "L1PA8A", "L1PA10", "L1PA11", "L1PA12", "L1PA13", "L1PA14", "L1PA15", "L1PA15-16", "L1PA16", "L1PA17", "SVA_F") , ordered = TRUE)

library("ggplot2")
# plotting a barplot showing the results for each TE for all TEs and fl TEs
ggplot(df_L1PA_thres3, aes(L1PA, as.numeric(L1PA_perc), fill = all_or_fl)) +
  geom_bar(stat="identity", position = "dodge") +
  # adding labels and titles 
  labs(title = "% of L1PAs that have windows overlapping") +
  xlab("L1PA") +
  ylab("% L1PA regions with windows overlapping") +
  # adding a horizontal line for the FPR found at variance threshold 3
  geom_hline(yintercept = 12.502302, linetype="dashed", color = "black") +
  theme(legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=13)) +
  scale_fill_manual(values = c("#FFA000", "#3F93CB"),
                    name = "All or Full Length",
                    labels = c("All", "Full Length"))
