#Permutation Modelling Implementation

library(minet)


#Parameters

# 1. X is the N x P matrix of gene expression vectors, such that each column is the expression of a single gene
#
# 2. n_perms is the number of permuted variables which are added to the dataset. This selection is arbitrary. However,
# it is recommended that n_perms is set to some value less than the number of original gene variables.
#
# The MINET package provides the "aracne", "clr" and "mrnet" methods

perm_model <- function(X,
                       alpha=0.01,
                       n_perms=250,
                       method="aracne"){
  
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
  
  #Apply the threshold to the association matrix to classify the network
  g_hat <- g_int>t_alpha
  
  return(g_hat)
  
}