## -------------------------------------------------------------------------
## Making Fancy Genome-Methylation Differences Graphs
## -------------------------------------------------------------------------

# Load packages etc.
setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/coverage_files")
library(grid)
library(readr)
library(ggplot2)
library(methylKit)
library(ggdendro)
library(ggfortify)
library(ggrepel)

## -------------------------------------------------------------------------

# Get a methylkit object for all samples
sample.list <- list("old_queen_1.CpG_report.merged_CpG_evidence.cov" ,"old_queen_2.CpG_report.merged_CpG_evidence.cov",
                    "old_queen_3.CpG_report.merged_CpG_evidence.cov","old_queen_4.CpG_report.merged_CpG_evidence.cov",
                    "old_queen_5.CpG_report.merged_CpG_evidence.cov" ,
                    "young_queen_1.CpG_report.merged_CpG_evidence.cov","young_queen_2.CpG_report.merged_CpG_evidence.cov",
                    "young_queen_3.CpG_report.merged_CpG_evidence.cov","young_queen_4.CpG_report.merged_CpG_evidence.cov",
                    "young_queen_5.CpG_report.merged_CpG_evidence.cov","young_queen_6.CpG_report.merged_CpG_evidence.cov",
                    "old_worker_1.CpG_report.merged_CpG_evidence.cov" ,"old_worker_2.CpG_report.merged_CpG_evidence.cov",
                    "old_worker_3.CpG_report.merged_CpG_evidence.cov","old_worker_4.CpG_report.merged_CpG_evidence.cov",
                    "old_worker_5.CpG_report.merged_CpG_evidence.cov" ,"old_worker_6.CpG_report.merged_CpG_evidence.cov",
                    "old_worker_7.CpG_report.merged_CpG_evidence.cov",
                    "young_worker_1.CpG_report.merged_CpG_evidence.cov","young_worker_2.CpG_report.merged_CpG_evidence.cov",
                    "young_worker_3.CpG_report.merged_CpG_evidence.cov","young_worker_4.CpG_report.merged_CpG_evidence.cov",
                    "young_worker_5.CpG_report.merged_CpG_evidence.cov","young_worker_6.CpG_report.merged_CpG_evidence.cov")


CPGRaw <- methRead(sample.list, 
                   sample.id = list("Old_Q1", "Old_Q2","Old_Q3","Old_Q4","Old_Q5",
                                    "Young_Q1", "Young_Q2","Young_Q3","Young_Q4","Young_Q5","Young_Q6",
                                    "Old_W1", "Old_W2","Old_W3","Old_W4","Old_W5","Old_W6","Old_W7",
                                    "Young_W1", "Young_W2","Young_W3","Young_W4","Young_W5","Young_W6"),
                   assembly="ant",
                   treatment=c(0,0,0,0,0,1,1,1,1,1,1,3,3,3,3,3,3,3,4,4,4,4,4,4),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir= path)

filtered_data <- filterByCoverage(CPGRaw,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)
rm(CPGRaw)
normalized <- normalizeCoverage(filtered_data)
rm(filtered_data)
meth_all_data <- unite(normalized, destrand=F) 
rm(normalized)
nrow(meth_all_data) # 1463570

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
n <- df_meth_all[,38:39]
a1 <- df_meth_all[,41:42]
b1 <- df_meth_all[,44:45]
c1 <- df_meth_all[,47:48]
d1 <- df_meth_all[,50:51]
e1 <- df_meth_all[,53:54]
f1 <- df_meth_all[,56:57]
g1 <- df_meth_all[,59:60]
h1 <- df_meth_all[,62:63]
k1 <- df_meth_all[,65:66]
l1 <- df_meth_all[,68:69]
m1 <- df_meth_all[,71:72]
n1 <- df_meth_all[,74:75]

bt <- function(a, b, p = 0.005) {binom.test(a, b, 0.005, alternative="greater")$p.value}
allrows <- data.frame(CT=numeric(),Ccount=numeric(),pVal=numeric(),FDR=numeric(),row=numeric())

for (df in list(a,b,c,d,e,f,g,h,k,l,m,n)) {
  colnames(df) <- c("CT", "Ccount")
  df$pVal <- mapply(bt, df$Ccount, df$CT)
  df$FDR <- p.adjust(df$pVal, method = "BH", n = length(df$pVal))
  df$row <- as.numeric(rownames(df))
  dfmeth <- subset(df, df$FDR < 0.05)
  allrows <- rbind(allrows, dfmeth)
}

meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions) # 6514

subset_methBase <- methylKit::select(meth_all_data, meth_positions)
head(subset_methBase)

## -------------------------------------------------------------------------
# Make a nice PCA

PCA_data <- PCASamples(subset_methBase, screeplot=F, obj.return = T)
PCA_data1 <- as.data.frame(PCA_data$x)
PCA_data1$sample <- row.names(PCA_data1)

PCA_data1$type <- c(rep("Old_Q", 5), rep("Young_Q", 6), rep("Old_W", 7), rep("Young_W", 6))


percentage <- round(PCA_data$sdev / sum(PCA_data$sdev) * 100, 0)
percentage <- paste(colnames(PCA_data), "(", paste( as.character(percentage), "%", ")", sep="") )


ggplot(PCA_data1, aes(PC1, PC2, colour=type))+
  geom_point(size=14)+
  geom_text_repel(aes(label=sample), size=12,show.legend=FALSE, 
                  point.padding = 2, box.padding = 1)+
  theme_bw()+
  xlab(paste0("PC1:",percentage[1],"variance")) +
  ylab(paste0("PC2:",percentage[2],"variance")) +
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=30),
        legend.text=element_text(size=30),
        legend.title=element_blank())+
  scale_colour_manual(values=c("#0066CC","#99CCFF","#663399","#CC99FF"))
