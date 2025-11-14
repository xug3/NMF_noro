library(R.matlab)
library(FAST)
library(openxlsx)
library(tidyverse)
library(lubridate)
library(ComplexHeatmap)
library(circlize)
setwd("D:/RAP UNLV/research projects/waste water/norovirus/data")

all_sample=read.csv("Trimmed_SNP_matrix/matrix_freq_1_7_Days_Wastewater.csv")
all_sample2=read.csv("Trimmed_SNP_matrix/matrix_freq_gt_1_Week_Wastewater.csv")
all_sample3=read.csv("Trimmed_SNP_matrix/matrix_freq_gt_1_Month_Wastewater.csv")
all_sample4=read.csv("Trimmed_SNP_matrix/matrix_freq_gt_6_Months_Wastewater.csv")
all_sample5=read.csv("Trimmed_SNP_matrix/matrix_freq_gt_1_Year_Wastewater.csv")
all_sample=cbind(all_sample,all_sample2[,-c(1:5)],all_sample3[,-c(1:5)],all_sample4[,-c(1:5)],all_sample5[,-c(1:5)])
dim(all_sample)
length(unique(colnames(all_sample)))
all_sample[1:10,1:10]
table(all_sample[,5])
waste_sample=all_sample[all_sample[,5]=="",]
dim(waste_sample)
waste_sample[1:10,1:10]
table(waste_sample[,4])
table(waste_sample[,3])

clinic_sample=all_sample[-which(all_sample[,5]==""),]
dim(clinic_sample)
clinic_sample[1:10,1:10]
table(clinic_sample[,4])
table(clinic_sample[,5])

rownames(all_sample)=all_sample[,1]

# check out data structure
countmatrixt = t(all_sample[,-c(1:5)])
#countmatrixt = t(cross_dsads[[1]])
dim(countmatrixt)
mutationname=rownames(countmatrixt)


date_sample=as.Date(all_sample[,4])


adjacent_matrix=diag(dim(all_sample)[1])

Rank=2:20
Fnorm=rep(0,length(Rank))
for(i in 1:length(Rank)){
  config <- dmain_config
  config$r <- Rank[i]
  config$lambda_1 <- 0.1
  config$lambda_2 <- 0.1
  config$n_iter <- 100000
  set.seed(1)
  res<-dmain(countmatrixt,adjacent_matrix,config)
  reconstr=res$W%*%t(res$H)
  Fnorm[i]=sum((countmatrixt-reconstr)^2)
}
tiff("noro_wastefirstall_rankvsRSS.tiff",width = 8,height = 6,units = "in",res = 300,compression = "lzw",type = "cairo")
plot(Rank,Fnorm,ylab="RSS")
dev.off()

config <- dmain_config
config$r <- 11
config$lambda_1 <- 0.1
config$lambda_2 <- 0.1
config$n_iter <- 100000
set.seed(1)
res<-dmain(countmatrixt,adjacent_matrix,config)
str(res)

#res=readRDS("old2022Q2_morethan2sample_20components.rds")
#saveRDS(res,"allwastewater_11rank.rds")
res=readRDS("allwastewater_11rank.rds")

res$H[1:5,]
loading=res$H
loading[loading<1e-5]=0
for(i in 1:dim(loading)[1]) loading[i,]=loading[i,]/sum(loading[1,])
rownames(loading)=colnames(countmatrixt)
colnames(loading)=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11")
wastewater_loading=loading[which(all_sample[,5]==""),]  
rownames(wastewater_loading)=waste_sample[,4]
order_wastewater_loading=wastewater_loading[order(rownames(wastewater_loading)),]
rownames(order_wastewater_loading)

clinical_loading=loading[-which(all_sample[,5]==""),]  
rownames(clinical_loading)=clinic_sample[,4]
order_clinical_loading=clinical_loading[order(rownames(clinical_loading)),]


trend_table_365_lm_min <- function(M, alpha = 0.05){
  
  dates  <- as.Date(dates_chr)
  t_year <- as.numeric(dates - min(dates)) / 365   
  t_c    <- t_year - mean(t_year, na.rm = TRUE)     
  
  coln <- colnames(M)
  slope <- pval <- rep(NA_real_, length(coln))
  
  for (j in seq_along(coln)) {
    y  <- as.numeric(M[, j])
    ok <- is.finite(y) & is.finite(t_c)
    if (sum(ok) < 3L || var(y[ok], na.rm = TRUE) == 0) next
    
    fit <- lm(y[ok] ~ t_c[ok])
    co  <- summary(fit)$coefficients
    slope[j] <- unname(co["t_c[ok]", "Estimate"])    
    pval[j]  <- unname(co["t_c[ok]", "Pr(>|t|)"])
  }
  
  p_adj <- p.adjust(pval, method = "BH")
  
  trend <- ifelse(is.na(slope), NA_character_,
                  ifelse(slope > 0, "increasing",
                         ifelse(slope < 0, "decreasing", "flat")))
  significance <- ifelse(is.na(p_adj), NA_character_,
                         ifelse(p_adj < alpha, "significant", "not significant"))
  
  out <- data.frame(
    series       = coln,
    slope_per_year = slope,
    p_value      = pval,
    p_adj_BH     = p_adj,
    trend        = trend,
    significance = significance,
    stringsAsFactors = FALSE
  )
  out[order(out$p_adj_BH), ]
}



res_lm <- trend_table_365_lm_min(order_wastewater_loading)
print(res_lm)
write.csv(res_lm,"wastewater_timetrend.csv")

res2 <- trend_table_365_lm_min(order_clinical_loading)
print(res2)
write.csv(res2,"clinical_timetrend.csv")

sourcematrix=res$W
sourcematrix[sourcematrix<1e-5]=0
sourcematrix[1:10,]

dim(sourcematrix)
rownames(sourcematrix)=rownames(countmatrixt)
rownames(sourcematrix)=gsub("[ACGT]", "", rownames(countmatrixt))
colnames(sourcematrix)=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11")
order_source=order(as.numeric(rownames(sourcematrix)))
order_sourcematrix=sourcematrix[order_source,]

reconstr=res$W%*%t(res$H)
sum((countmatrixt-reconstr)^2)



mat <- as.matrix(t(order_sourcematrix))
mat
col_fun = colorRamp2(c(0, max(mat)), c( "white", "red"))

write.csv(order_sourcematrix,"noro_wastefirstall_mutations.csv")

tiff("noro_wastefirstall_mutations.tiff",width = 36,height = 6,units = "in",res = 300,compression = "lzw",type = "cairo")
Heatmap(mat,col=col_fun,
        name = "metric",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8, rot = 90),
        show_row_dend = FALSE,
        show_column_dend = FALSE)

dev.off()

tiff("noro_wastefirstall_wastewater_loading.tiff",width = 18,height = 6,units = "in",res = 300,compression = "lzw",type = "cairo")
col_fun = colorRamp2(c(0, max(order_wastewater_loading)), c( "white", "red"))
Heatmap(t(order_wastewater_loading),col=col_fun,
        name = "metric",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8, rot = 90),
        show_row_dend = FALSE,
        show_column_dend = FALSE)
dev.off()

tiff("noro_wastefirstall_clinical_loading.tiff",width = 18,height = 6,units = "in",res = 300,compression = "lzw",type = "cairo")
col_fun = colorRamp2(c(0, max(order_clinical_loading)), c( "white", "red"))
Heatmap(t(order_clinical_loading),col=col_fun,
        name = "metric",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8, rot = 90),
        show_row_dend = FALSE,
        show_column_dend = FALSE)
dev.off()


component6=sourcematrix[sourcematrix[,6]>0.1,6]
names(component6)

order_clinical_loading[,1]
order_wastewater_loading[,1]

subsample=all_sample[c("Jul24b_NV_GII175","PQ196512.1","PQ215545.1","PQ215546.1"),names(component6)]
subsample=as.matrix(subsample)
rownames(subsample)=all_sample[c("Jul24b_NV_GII175","PQ196512.1","PQ215545.1","PQ215546.1"),4]
tiff("noro_component6_firstdetection.tiff",width = 18,height = 6,units = "in",res = 300,compression = "lzw",type = "cairo")
col_fun = colorRamp2(c(0, max(subsample)), c( "white", "red"))
Heatmap(subsample,col=col_fun,
        name = "metric",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8, rot = 90),
        show_row_dend = FALSE,
        show_column_dend = FALSE)
dev.off()


component1=sourcematrix[sourcematrix[,1]>0.3,1]
names(component1)

order_clinical_loading[order_clinical_loading[,1]>0.01,1]
order_wastewater_loading[order_wastewater_loading[,1]>0.01,1]
firstday=names(c(order_wastewater_loading[order_wastewater_loading[,1]>0.01,1][1],order_clinical_loading[order_clinical_loading[,1]>0.01,1][1]))
firstday=c("2023-02-06","2023-02-08")
subsample=all_sample[,names(component1)]
subsample[c("Jul24b_NV_GII175","PQ196512.1","PQ215545.1","PQ215546.1"),]
subsample=as.matrix(subsample)
rownames(subsample)=all_sample[,4]
subsample[firstday,]
which(rownames(subsample)==firstday[1])
tiff("noro_component6_firstdetection.tiff",width = 18,height = 6,units = "in",res = 300,compression = "lzw",type = "cairo")
col_fun = colorRamp2(c(0, max(subsample)), c( "white", "red"))
Heatmap(subsample[firstday,],col=col_fun,
        name = "metric",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8, rot = 90),
        show_row_dend = FALSE,
        show_column_dend = FALSE)
dev.off()

