# Utils ---------------------------------------------------------------------------------------

my_log <- function(x){
  ifelse(is.infinite(log(x)),0,log(x))
}

calc_RR <- function(signal_vec, signal_AE, node_ind){
  RR <- rep(1,length(node_ind))
  RR[signal_AE] <- signal_vec
  for(i in 1:length(node_ind)){
    RR[i] <- mean(RR[node_ind[[i]]])
  }
  return(RR)
}

comb <- function(x, ...) {  
  mapply(c,x,...,SIMPLIFY=FALSE)
}

# Tree function -----------------------------------------------------------

leaf_cut <- function(leaf, lowest_node=F){
  if(lowest_node==F){
    leaf_list <- strsplit(leaf,"_")
    cut_list <- rep(NA,length(leaf_list))
    for(i in 1:length(leaf_list)){
      cut_list[i] <- paste(leaf_list[[i]][-length(leaf_list[[1]])],collapse="_")
    }  
    upper_level <- unique(cut_list)
  }else{
    cut_list <- leaf
    upper_level <- leaf
  }
  cut_ind <- list()
  for(i in 1:length(upper_level)){
    cut_ind[[i]] <- which(cut_list%in%upper_level[i])
  }
  return(list(upper_level=upper_level, cut_ind=cut_ind))
}

# Matrix Information ------------------------------------------------------

transfrom_matrix <- function(case, control, node){
  if(length(node)==1){
    mat <- matrix(c(sum(case[,node]==1&control[,node]==1),
                    sum(case[,node]!=1&control[,node]==1),
                    sum(case[,node]==1&control[,node]!=1),
                    sum(case[,node]!=1&control[,node]!=1)),
                  nrow=2,ncol=2)
  }else{
    mat <- matrix(c(sum(rowSums(case[,node])!=0&rowSums(control[,node])!=0),
                    sum(rowSums(case[,node])==0&rowSums(control[,node])!=0),
                    sum(rowSums(case[,node])!=0&rowSums(control[,node])==0),
                    sum(rowSums(case[,node])==0&rowSums(control[,node])==0)),
                  nrow=2,ncol=2)
  }
  return(mat)
}

extract_mat_inf <- function(case, control, null, node_ind){
  res <- matrix(nrow=ncol(case), ncol=3)
  for(i in 1:length(node_ind)){
    tab <- transfrom_matrix(case, control, node_ind[[i]])
    if(null==T){
      res[i,] <- c(tab[1,1],(tab[1,2]+tab[2,1])/2,(tab[1,2]+tab[2,1])/2)
    }else{
      res[i,] <- c(tab[1,1],tab[1,2],tab[2,1]) 
    }
  }
  return(res)
}

# Statistics Function -----------------------------------------------------

unconditonal_bernoulli <- function(mat, p=0.5){
  c <- mat[1,2]+mat[1,1]
  n <- mat[2,1]+mat[1,1]
  LLR <- (c*my_log(c/(c+n)) + n*my_log(n/(c+n)) - c*log(p) - n*log(1-p)) * ifelse(c/(c+n)>0.5,1,0)
  return(ifelse(c==0&n==0,0,LLR))
}

Wald <- function(mat){
  b <- mat[1,2]
  c <- mat[2,1]
  LLR <- (my_log(b)-my_log(c))^2/(1/b+1/c)
  LLR <- ifelse(b>c,LLR,0)
  return(LLR)
}

McNemar <- function(mat){
  b <- mat[1,2]
  c <- mat[2,1]
  res <- ((b-c)^2/(b+c))*(b>c)
  return(ifelse(b==0&c==0,0,res))
}

LRT <- function(mat){
  b <- mat[1,2]
  c <- mat[2,1]
  exp_coef <- b/c
  if(c==0){
    res <- -2*((b+c)*log(0.5))
  }else{
    res <- -2*((b+c)*log(0.5)-(b*log(exp_coef/(exp_coef+1))+c*log(1/(exp_coef+1))))
  }
  return(ifelse(b==0&c==0|b<=c,0,res))
}

# Monte Carlo Simulation --------------------------------------------------

calc_LLR <- function(mc_data, node_ind, only_max=T){
  LLR_ub <- LLR_Wald <- LLR_Mc <- LLR_LRT <- rep(NA,length(node_ind))
  for(i in 1:length(node_ind)){
    LLR_ub[i] <- unconditonal_bernoulli(transfrom_matrix(mc_data$case, mc_data$control, node_ind[[i]]))
    LLR_Wald[i] <- Wald(transfrom_matrix(mc_data$case, mc_data$control, node_ind[[i]]))
    LLR_Mc[i] <- McNemar(transfrom_matrix(mc_data$case, mc_data$control, node_ind[[i]]))
    LLR_LRT[i] <- LRT(transfrom_matrix(mc_data$case, mc_data$control, node_ind[[i]]))
  }
  if(only_max==T){
    return(list(LLR_ub=max(LLR_ub), LLR_Wald=max(LLR_Wald), LLR_Mc=max(LLR_Mc), LLR_LRT=max(LLR_LRT)))
  }else{
    return(list(LLR_ub=LLR_ub, LLR_Wald=LLR_Wald, LLR_Mc=LLR_Mc, LLR_LRT=LLR_LRT))
  }
}

mc_gen_pair <- function(mat_info, n_row, n_col){
  type_p <- mat_info[,c(1,2)]/apply(mat_info[,c(1,2)],1,sum)
  
  s_case <- s_control <- matrix(0, nrow=n_row, ncol=n_col)
  case_p <- control_p <- apply(mat_info[,1:2], 1, sum)/sum(mat_info[,1:2])
  
  for(i in 1:n_row){
    case_AE <- sample(1:n_col, 1, prob=case_p)
    s_case[i,case_AE] <- 1
    if(sample(1:2, 1, prob=type_p[case_AE,])==1){
      s_control[i,case_AE] <- 1
    }else{
      control_AE <- sample(c(1:n_col)[-case_AE], 1, prob=control_p[-case_AE]/sum(control_p[-case_AE]))
      s_control[i,control_AE] <- 1
    }
  }
  return(list(case=s_case, control=s_control))
}

mc_gen_unpair <- function(mat_info, n_row, n_col){
  type_p <- mat_info[,c(1,2)]/apply(mat_info[,c(1,2)],1,sum)
  
  s_case <- s_control <- matrix(0, nrow=n_row, ncol=n_col)
  case_p <- control_p <- apply(mat_info[,1:2], 1, sum)/sum(mat_info[,1:2])
  
  for(i in 1:n_row){
    case_AE <- sample(1:n_col, 1, prob=case_p)
    s_case[i,case_AE] <- 1
    control_AE <- sample(1:n_col, 1, prob=control_p)
    s_control[i,control_AE] <- 1  
  }
  return(list(case=s_case, control=s_control))
}

mc_dist <- function(iter, mat_info, node_ind, pair, n_row, n_col){
  
  max_LLR_ub <- max_LLR_Wald <- max_LLR_Mc <- max_LLR_LRT <- rep(0,iter)
  
  for(j in 1:iter){
    if(pair==T){
      mc_data <- mc_gen_pair(mat_info, n_row, n_col)
    }else{
      mc_data <- mc_gen_unpair(mat_info, n_row, n_col)
    }
    max_LLR <- calc_LLR(mc_data, node_ind, only_max=T)
    max_LLR_ub[j] <- max_LLR$LLR_ub
    max_LLR_Wald[j] <- max_LLR$LLR_Wald
    max_LLR_Mc[j] <- max_LLR$LLR_M
    max_LLR_LRT[j] <- max_LLR$LLR_LRT
  }
  return(list(ub_dist=max_LLR_ub, Wald_dist=max_LLR_Wald, Mc_dist=max_LLR_Mc, LRT_dist=max_LLR_LRT))
}


simulation <- function(simul_num, AEs, mat_info, num_pair=1000, num_signal=5, signal=1.5, mc_num=100){
  
  SOC_ARRN <- leaf_cut(AEs, lowest_node=T) ## SOC_ARRN level
  ARRN <- leaf_cut(AEs) ## ARRN level
  
  node <- c(SOC_ARRN$upper_level,ARRN$upper_level) ## Combine SOC_ARRN, ARRN Node ID
  node_ind <- append(SOC_ARRN$cut_ind,ARRN$cut_ind) ## Combine SOC_ARRN, ARRN Node Index
  
  num_AE <- length(AEs)
  
  signal_AE <- sample(1:num_AE, num_signal, replace=F, prob=apply(mat_info,1,sum)/sum(mat_info))
  
  if(signal==-1){
    signal <- runif(num_signal,1,5)
  }
  
  mat_info[signal_AE,2] <- (mat_info[signal_AE,2]*2)*signal/(signal+1)
  mat_info[signal_AE,3] <- (mat_info[signal_AE,3]*2)/(signal+1)
  
  case <- control <- matrix(0, nrow=num_pair, ncol=num_AE)
  
  case_p <- apply(mat_info[,1:2], 1, sum)/sum(mat_info[,1:2])
  control_p <- apply(mat_info[,c(1,3)], 1, sum)/sum(mat_info[,c(1,3)])
  type_p <- mat_info[,c(1,2)]/apply(mat_info[,c(1,2)],1,sum)
  
  for(i in 1:num_pair){
    case_AE <- sample(1:num_AE, 1, prob=case_p)
    case[i,case_AE] <- 1
    if(sample(1:2, 1, prob=type_p[case_AE,])==1){
      control[i,case_AE] <- 1
    }else{
      control_AE <- sample(c(1:num_AE)[-case_AE], 1, prob=control_p[-case_AE]/sum(control_p[-case_AE]))
      control[i,control_AE] <- 1
    }
  }
  
  RR <- calc_RR(signal, signal_AE, node_ind)
  LLR <- calc_LLR(list(case=case,control=control), node_ind, only_max=F)
  
  mc_res <- mc_dist(iter=mc_num, mat_info=mat_info, node_ind=node_ind, pair=T, n_row=num_pair, n_col=num_AE)
  
  cut_ub_pair <- quantile(mc_res$ub_dist,0.95)
  cut_Mc_pair <- quantile(mc_res$Mc_dist,0.95)
  cut_Wald_pair <- quantile(mc_res$Wald_dist,0.95)
  cut_LRT_pair <- quantile(mc_res$LRT_dist,0.95)
  
  detect_ub_pair <- detect_Mc_pair <- detect_Wald_pair <- detect_LRT_pair <-
    rep(NA, length(LLR$LLR_ub))
  
  for(k in 1:length(LLR$LLR_ub)){
    detect_ub_pair[k] <- ifelse(LLR$LLR_ub[k]>=cut_ub_pair,1,0)
    detect_Mc_pair[k] <- ifelse(LLR$LLR_Mc[k]>=cut_Mc_pair,1,0)
    detect_Wald_pair[k] <- ifelse(LLR$LLR_Wald[k]>=cut_Wald_pair,1,0)
    detect_LRT_pair[k] <- ifelse(LLR$LLR_LRT[k]>=cut_LRT_pair,1,0)
  }
  
  return(data.frame(iter=rep(simul_num,length(RR)), node=node, RR=RR, 
                    LLR_ub=LLR$LLR_ub, LLR_Mc=LLR$LLR_Mc, LLR_Wald=LLR$LLR_Wald, LLR_LRT=LLR$LLR_LRT,
                    detect_ub_pair, detect_Mc_pair, detect_Wald_pair, detect_LRT_pair))
}

