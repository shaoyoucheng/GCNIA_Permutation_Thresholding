#Permutation Modelling Benchmark Analysis

#setup
if(T){
  
  library(minet)
  library(netbenchmark)
  library(MASS)
  library(Matrix)
  library(infotheo)
  library(ggplot2)
  
  #n_perms is the number of permuted variables which are added to the dataset. This selection is arbitrary. However,
  #it is recommended that n_perms is set to some value less than the number of original gene variables.
  n_perms <- 250
  
  #Alpha is the statistical significance cutoff. The (1-alpha)th quantile of the null distribution will be used to
  #estimate the hard threshold
  alpha<- 0.01
  
  #Generate grid of experimental parameters to include all datasets in the Netbenchmark package
  param_grid <- expand.grid(Method=c("clr", "aracne", "mrnet"),
                            Dataset=Availabledata[1:5],
                            N_experiments= c(50,250,500),
                            Noise=seq(from=0, to = 40, length.out=10),
                            stringsAsFactors = F)
  
}

#Measure Performance
if(T){
  performance_lists <- sapply(1:nrow(param_grid), function(params_i){
    
    #params
    if(T){
      
      #Assign task parameters
      message(params_i)
      params <- param_grid[params_i,]
      m_i <- params$Method
      dataset_i <- params$Dataset
      n_experiments_i <- params$N_experiments
      noise_i <- params$Noise
      
      #Load benchmark datase
      aux <- grndata::getData(dataset_i)
      datasource <- aux[[1]]
      theta <- aux[[2]]
      
      #Sample dataset with Netbenchmark and add noise.
      X <-  datasource.subsample(datasource, datasets.num = 1, 
                                 experiments = n_experiments_i, local.noise = noise_i, 
                                 global.noise = noise_i)[[1]]
      
    }
    
    #estimate network
    if(T){
      
      #Generate permuted variables
      X_perm <- apply(X[,sample(1:ncol(X),
                                size = n_perms,
                                replace=T)], 2, sample)
      
      #Concatenate the original variables with the permuted variables
      X2 <- cbind(X,X_perm)
      
      #Perform gene network inference on the concatnated dataset
      g <- minet(dataset=X2,
                 method=m_i)
      
      #Isolate the associations between original and permuted variables
      #These associations are assumed to be a sample of the null distribution
      ext_cols <- rep(c(F,T),times=c(ncol(X), n_perms))
      g_ext <- g[!ext_cols, ext_cols]
      
      #Select threshold according to (1-alpha)th quantile of the null distribution
      t_alpha <- quantile(g_ext, 1-alpha)
      
      #Isolate the associations between the original gene expression variables
      g_int <- g[!ext_cols, !ext_cols]
      
    }
    
    #measure performance
    if(T){
      
      #validate performance against ground truth data
      r <- suppressWarnings(validate(inet=g_int,
                                     tnet=theta))
      p_mat <- cbind(pr(r), fscore(r))
      p_mat <- cbind(p_mat, rank(p_mat[,3])/nrow(p_mat))
      
      #find threshold which achieved highest F1 score (oracle)
      opt_row <-  which.max(p_mat[,3])
      
      #find performance of the estimated threshold
      t_alpha_row <- max(which(r[,1]>t_alpha))
      res_mat <- p_mat[c(t_alpha_row, opt_row),]
      dimnames(res_mat)<- list(c("alpha", "opt"),
                               c("precision", "recall", "F1", "Rank"))
      res_list <- c(params, unlist(res_mat))
      res_list <- c(res_list, mean_f1=mean(p_mat[,3]))
      
      return(res_list)
      
    }
    
  })
  
  #configure performance df
  if(T){
    
    performance_df <- as.data.frame(t(performance_lists))
    performance_df[] <- lapply(performance_df, unlist) 
    
    colnames(performance_df)[5:10] <- paste0(rep(c("Precision", "Recall", "F1"), each=2),"_",c("Est", "Opt"))
    #performance_df <- performance_df[,-12]
    
    colnames(performance_df)[11] <- "F1_Rank"
    performance_df$Dataset <- factor(performance_df$Dataset, levels=c(
      "rogers1000","syntren300","syntren1000","gnw1565","gnw2000"
    ))
    performance_df$N_experiments <- factor(paste("N =", performance_df$N_experiments),
                                           levels = paste("N =", unique(performance_df$N_experiments)))
    performance_df$F1_recovery <- (performance_df$F1_Est)/performance_df$F1_Opt
    
    if(F){
      save(performance_df, file="performance_df.RData")
    }
    if(F){
      load("performance_df.RData")
    }

  }


}

#Output Results
if(T){

  p1 <- ggplot()+
    geom_boxplot(data=performance_df,
                 aes(x=Method,
                     y=F1_Rank,
                     colour=Method))+
    facet_grid(Dataset~N_experiments)+
    ylab("F1_Rank")+
    theme(legend.position = "none")
  p1+scale_y_continuous(sec.axis = sec_axis(~.,name = "Dataset",  breaks=NULL))
  
  p2 <- ggplot()+
    geom_boxplot(data=performance_df,
                 aes(x=Method,
                     y=F1_recovery,
                     colour=Method))+
    facet_grid(Dataset~N_experiments, scales = "free_y")+ylab("F1_Recovery")+
    theme(legend.position = "none")
  p2+scale_y_continuous(sec.axis = sec_axis(~.,name = "Dataset",  breaks=NULL)) 
  
  
  results_agg_tbl <- aggregate(performance_df[,5:10], by=performance_df[,3:2], FUN=mean)
  results_agg_tbl <- format(results_agg_tbl, digits=3)
  
  write.csv(results_agg_tbl, "tables/results_agg_tbl.csv")
  save(results_agg_tbl, file="tables/results_agg_tbl.RData")

}




