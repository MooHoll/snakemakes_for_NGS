#------------------------------------------------
# Differential methylation between ages
#------------------------------------------------
setwd("~/Dropbox/Leicester_postdoc/Projects/Snakemake/cov_files")
library(methylKit)
library(readr)


sample.list <- list("old_queen_1.CpG_report.merged_CpG_evidence.cov" ,"old_queen_2.CpG_report.merged_CpG_evidence.cov",
                    "old_queen_3.CpG_report.merged_CpG_evidence.cov","old_queen_4.CpG_report.merged_CpG_evidence.cov",
                    "old_queen_5.CpG_report.merged_CpG_evidence.cov" ,
                    "young_queen_1.CpG_report.merged_CpG_evidence.cov","young_queen_2.CpG_report.merged_CpG_evidence.cov",
                    "young_queen_3.CpG_report.merged_CpG_evidence.cov","young_queen_4.CpG_report.merged_CpG_evidence.cov",
                    "young_queen_5.CpG_report.merged_CpG_evidence.cov","young_queen_6.CpG_report.merged_CpG_evidence.cov")

CPGRaw <- methRead(sample.list, 
                   sample.id = list("Old_Q1", "Old_Q2","Old_Q3","Old_Q4","Old_Q5",
                                    "Young_Q1", "Young_Q2","Young_Q3","Young_Q4","Young_Q5","Young_Q6"),
                   assembly="ant",
                   treatment=c(0,0,0,0,0,1,1,1,1,1,1),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir= path)

# ---

sample.list <- list("old_worker_1.CpG_report.merged_CpG_evidence.cov" ,"old_worker_2.CpG_report.merged_CpG_evidence.cov",
                    "old_worker_3.CpG_report.merged_CpG_evidence.cov","old_worker_4.CpG_report.merged_CpG_evidence.cov",
                    "old_worker_5.CpG_report.merged_CpG_evidence.cov" ,"old_worker_6.CpG_report.merged_CpG_evidence.cov",
                    "old_worker_7.CpG_report.merged_CpG_evidence.cov",
                    "young_worker_1.CpG_report.merged_CpG_evidence.cov","young_worker_2.CpG_report.merged_CpG_evidence.cov",
                    "young_worker_3.CpG_report.merged_CpG_evidence.cov","young_worker_4.CpG_report.merged_CpG_evidence.cov",
                    "young_worker_5.CpG_report.merged_CpG_evidence.cov","young_worker_6.CpG_report.merged_CpG_evidence.cov")

CPGRaw <- methRead(sample.list, 
                   sample.id = list("Old_W1", "Old_W2","Old_W3","Old_W4","Old_W5","Old_W6","Old_W7",
                                    "Young_W1", "Young_W2","Young_W3","Young_W4","Young_W5","Young_W6"),
                   assembly="ant",
                   treatment=c(0,0,0,0,0,0,0,1,1,1,1,1,1),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir= path)


# Filter by coverage
filtered_data <- filterByCoverage(CPGRaw,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)

normalized <- normalizeCoverage(filtered_data)

# Select only CpGs found in all alignments
meth_all_data <- unite(normalized, destrand=F, min.per.group = 3L) 
nrow(meth_all_data) 
# queen: 
# worker: 

## -------------------------------------------------------------------------

# Filter sites using a binomial test so only keep CpGs which are methylated in at least one sample
df_meth_all <- getData(meth_all_data)

a <- df_meth_all[,5:6]
b <- df_meth_all[,8:9]
c <- df_meth_all[,11:12]
d <- df_meth_all[,14:15]
e <- df_meth_all[,17:18]
f <- df_meth_all[,20:21]
g <- df_meth_all[,23:24]
h <- df_meth_all[,26:27]
k <- df_meth_all[,29:30]
l <- df_meth_all[,32:33]
m <- df_meth_all[,35:36]
#n <- df_meth_all[,38:39]
#o <- df_meth_all[,41:42]

# NOTE: p shouold be the average non-conversion rate (proportion of methylated Cs compared to non-meth Cs)
# So if 1000 methylated Cs compared to 200,000 T's then 1000/200,000 = 0.005
# for a paper: 'the success probability is the non-conversion rate'
bt <- function(a, b, p = 0.005) {binom.test(a, b, 0.005, alternative="greater")$p.value}
allrows <- data.frame(CT=numeric(),Ccount=numeric(),pVal=numeric(),FDR=numeric(),row=numeric())

for (df in list(a,b,c,d,e,f,g,h,k,l,m)) {
#for (df in list(a,b,c,d,e,f,g,h,k,l,m,n,o)) {
  colnames(df) <- c("CT", "Ccount")
  df$pVal <- mapply(bt, df$Ccount, df$CT)
  df$FDR <- p.adjust(df$pVal, method = "BH", n = length(df$pVal))
  df$row <- as.numeric(rownames(df))
  dfmeth <- subset(df, df$FDR < 0.05)
  allrows <- rbind(allrows, dfmeth)
}

meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions) 
# queens:
# workers:

subset_methBase <- methylKit::select(meth_all_data, meth_positions)
head(subset_methBase)

# File for Eamonn's clock
perc.meth <- percMethylation(subset_methBase)
write.table(perc.meth, file="percentage_meth_ant_data_queens.txt", sep = '\t',
            quote = F, col.names = T, row.names = F)
#write.table(perc.meth, file="percentage_meth_ant_data_workers.txt", sep = '\t',
#            quote = F, col.names = T, row.names = F)


diff_meth <- calculateDiffMeth(subset_methBase, mc.cores = 2)

diff_meth_5 <- getMethylDiff(diff_meth, difference=10, qvalue=0.01)
nrow(diff_meth_5)

write.csv(diff_meth_5, file="queens_diff_meth_CpGs.csv")
#write.csv(diff_meth_5, file="workers_diff_meth_CpGs.csv")

# queens (+ve = young hyper): 
# workers (+ve = young hyper):  
