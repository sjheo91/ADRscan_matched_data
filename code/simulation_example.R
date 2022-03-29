rm(list=ls())

library(doSNOW)
library(data.table)
source('function.R')

cl <- makeCluster(4) 
registerDoSNOW(cl) # stopCluster(cl)

## data generating
setwd('C:/Users/sjheo/Desktop/교내연구/Drug safety surveillance_2020/ADRscan_matched_data/code')
data <- read.csv('example_data.csv')

AEs <- sort(unique(as.character(data$WHOART2lv)))
SOC_ARRN <- leaf_cut(AEs, lowest_node=T) ## SOC_ARRN level
ARRN <- leaf_cut(AEs) ## ARRN level

node <- c(SOC_ARRN$upper_level,ARRN$upper_level) ## Combine SOC_ARRN, ARRN Node ID
node_ind <- append(SOC_ARRN$cut_ind,ARRN$cut_ind) ## Combine SOC_ARRN, ARRN Node Index

num_AEs <- length(AEs)
num_pairs <- length(unique(data$matched_id))

case <- control <- matrix(0, nrow=num_pairs, ncol=num_AEs)
for(i in 1:num_pairs){
  case[i,which(AEs%in%data[data$DRUG_CHEM=='atorvastatin'&data$matched_id==i,'WHOART2lv'])] <- 1
  control[i,which(AEs%in%data[data$DRUG_CHEM=='rosuvastatin'&data$matched_id==i,'WHOART2lv'])] <- 1
}

mat_info <- extract_mat_inf(case, control, null=T, node_ind[1:num_AEs])

mc_num <- 999
iter_num <- 100
num_pairs <- c(500, 1000, 2000)

for(k in 1:3){
  res <- foreach(i=1:iter_num, .combine='rbind', .multicombine=T) %dopar% {
    simulation(simul_num=i, AEs=AEs, mat_info=mat_info, num_pair=num_pairs[k], num_signal=1, signal=1, mc_num)
  }
}


