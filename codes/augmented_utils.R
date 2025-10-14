library(doParallel)
library(mclust)
library(multicool)
library(transport)
library(randtoolbox)
library(ggplot2)
library(gtools)


################# To compute errors


func_prep_clusts = function(clust, nrow_rep){ #To compute multi-dimensional wasserstein distance, the sample sizes must be equal
  if(length(clust)==0){return(clusts)}
  else{return(clust[rep(1:nrow(clust), nrow_rep%/% nrow(clust)),])}
} 

local_errors = function(clusts, rep, idx, indep = TRUE){ #compute the local errors, i.e. the wasserstein distances between each cluster and the each representative
  error= c()
  for(i in idx){
    if(length(clusts[[i]]) == 0 | length(rep[[i]]) == 0){error = c(error,0)} #if no elements, error is null
    else if(ncol(rep[[i]]) == 1){error= c(error, wasserstein1d(clusts[[i]], rep[[i]],p=2))} #if 1D, juste compute the 1D wasserstein distance 
    else{ #if multi-D
      if(!indep){ #if no independence, compute multi-d wasserstein distance
          clusts[[i]] = func_prep_clusts(clusts[[i]], nrow(rep[[i]])) #enrich cluster size to reach size of rep
          rep[[i]] = rep[[i]][1:nrow(clusts[[i]]),]
          error= c(error, transport::wasserstein(pp(clusts[[i]]), pp(rep[[i]]), p=2)) #compute multi-d wasserstein distance
      }
      else{error= c(error, sqrt(sum(sapply(1:ncol(rep[[i]]), function(k){wasserstein1d(clusts[[i]][,k], rep[[i]][,k],p=2)})^2)))} #If independence between marginals, then aggregate the 1D-wasserstein between the marginals
    }
  }
  return(error)
} 

compute_errors = function(clusts, rep, law, n = 10^6){ #Compute quantization error and global error
    weights = sapply(clusts, nrow)/sum(sapply(clusts, nrow)) #weights of each representative
    mixture = from_law_to_sample(rep = rep, law = law, n = n, quasi_random = FALSE, weights = weights) #Get a large sample for each representative, and size function of weight
    quanti_error = sqrt(weighted.mean(local_errors(clusts, mixture, 1:length(clusts), indep = TRUE)^2, weights)) #compute the quantization error
    mixture = do.call("rbind",mixture) #concatenate the samples of each representative
    samp = do.call("rbind", clusts) #concatenate the clusters
    global_error = local_errors(list(samp), list(mixture), 1, indep = TRUE) #compute global error: wasserstein distance between the sample and the mixture
    return(list(quanti_error = quanti_error, global_error = global_error))
}

wasserstein_between_mixtures = function(mixture1, mixture2, indep = TRUE){ #Compute wasserstein distance between two mixtures
    weights1 = sapply(mixture1$clusts, nrow)/sum(sapply(mixture1$clusts, nrow)) #get weights of first mixture
    mixture1 = do.call("rbind",from_law_to_sample(rep = mixture1$rep, law = mixture1$law, n = 10^5, quasi_random = FALSE, weights = weights1)) #get the large sample representing the first mixture
    weights2 = sapply(mixture2$clusts, nrow)/sum(sapply(mixture2$clusts, nrow)) #same with second mixture
    mixture2 = do.call("rbind",from_law_to_sample(rep = mixture2$rep, law = mixture2$law, n = 10^5, quasi_random = FALSE, weights = weights2))
    return(local_errors(list(mixture1), list(mixture2), 1, indep = indep)) #compute wasserstein distance between these two mixtures
    }

################ For find_c

dist2 = function(x,y){sqrt(sum((x-y)^2))} #euclidean distance

dist_point_sample = function(x, samp){ #Compute the distance between a point and a sample
  dists = apply(samp, 1, function(y){dist2(x,y)})
  return(min(dists))
}

from_unif_to_sample = function(samp, rep, law){ #from a multi-dimensional uniform sample and the parameter of a representative, it generates large sample for this representative
  mat = NULL
  for(ii in 1:length(law)){ #for each marginal
    if(law[ii] == "unif"){mat = cbind(mat,samp[,ii]*(rep[2,ii]-rep[1,ii])+rep[1,ii])}  #if uniform representative, just linear transformation
    if(law[ii] == "normal"){mat = cbind(mat, qnorm(samp[,ii])*rep[2,ii]+rep[1,ii])} #if gaussian, transform with qnorm
    else if(law[ii] == "dirac"){mat = cbind(mat, rep(rep[2,ii], nrow(samp)))} #if dirac, just repeat the same point
  }
  return(mat)
}

from_law_to_sample = function(rep, law = "unif", n, weights = NULL, quasi_random = TRUE){ #from a list of representatives, provides a large sample for each one
  if(quasi_random){sob = matrix(sobol(n, dim = ncol(rep[[1]])), ncol = ncol(rep[[1]]))} # quasi random uniform samples
  else{sob = matrix(runif(n*ncol(rep[[1]])), ncol = ncol(rep[[1]]))} #random uniform samples
  rep_bis = lapply(1:length(rep), function(i){from_unif_to_sample(sob,rep[[i]], law[[i]])}) #for each representative, use from_unif_to_sample
  if(!is.null(weights)){ #if weights are provided, adjust the size of each sample
    thresh = c(1,as.integer(round(cumsum(weights),4)*n))
    rep_bis = lapply(1:length(rep_bis),function(x){as.matrix(rep_bis[[x]][(thresh[x]:thresh[x+1]),])})
  }
  return(rep_bis)
}

find_c = function(samp, rep, weights,law, n_sample = 5000){ #find_c
  rep_large = from_law_to_sample(rep=rep, weights = weights, law=law, n = n_sample, quasi_random = FALSE) #general large sample for the mixture
  nums = lapply(1:length(rep), function(x){c()}) #will store the idx of the points belonging to each cluster
  if(parallel){clust_number = foreach(i=1:nrow(samp),.combine='c')%dopar%{ #parallel computation
    which.min(sapply(rep_large, function(gamm){dist_point_sample(as.numeric(samp[i,]), gamm)}))}} #associate to each point in the sample, its closest point in the mixture, and save the representative it originates from
  else{clust_number = foreach(i=1:nrow(samp),.combine='c')%do%{ #same but no parallel
    which.min(sapply(rep_large, function(gamm){dist_point_sample(as.numeric(samp[i,]), gamm)}))}}
  nums = lapply(1:length(rep), function(j){which(clust_number == j)}) #store the idx of the points belonging to each cluster
  if(0 %in% sapply(nums, length)){ #If one cluster is empty, give him its closest point in the sample
    vec_idx_min = c()
    for(i0 in which(sapply(nums, length)== 0)){ #for each empty cluster
      idx_min = which.min(Vectorize(function(i){dist_point_sample(samp[i,],rep_large[[i0]])})(setdiff(1:nrow(samp), vec_idx_min))) #find its closest point in the sample
      vec_idx_min = c(vec_idx_min, idx_min) #store this point so that it is not used by another empty cluster
      nums[[i0]] = idx_min
      nums[setdiff(1:length(nums),i0)] = lapply(setdiff(1:length(nums),i0), function(x){nums[[x]][nums[[x]] != idx_min]}) #remove this point from its previous cluster, and add it to the empty cluster
    }
  }
  clusts = lapply(nums, function(x){matrix(samp[x,], ncol = ncol(samp))})
  return(list(clusters = clusts,nums = nums))
}

################# For split


f_delete = function(clust, n, n_sample, prop_search = 1){ #f_delete identify the index of the elements to place into the bin
  if(nrow(clust)==1){return(1)}
  else{
    df_res = c() #df_res will store the indexes
    nums = 1:nrow(clust) 
    clust_init = clust
    while(length(df_res) < n){ #will the length of the bin is lower than the objective
      if(prop_search==1){search=1:nrow(clust)} #prop search is the proportion of candidate elements at each iteration
      else{search = sort(sample(1:nrow(clust), size = as.integer(round(prop_search*nrow(clust)))))}  #if prop_search < 1, select the points that will be tested for the transfer
      if(parallel){diff = foreach(j=search, .combine = "c")%dopar%{ #try all points for the transfer to the bin: each point is removed temporarily, and the clustering error is computed
        res_ctor = CtoR(clusts = list(matrix(clust[-j,], nrow = nrow(clust)-1),matrix(clust_init[c(df_res, nums[j]),], nrow = (1+length(df_res)))), only_bornes = FALSE, n_sample = n_sample, return_error = TRUE)
        weighted.mean(res_ctor[[3]]^2,res_ctor[[4]])
        }
      }#for each element we evaluate its adding to the bin
      else{diff = foreach(j=search, .combine = "c")%do%{
        res_ctor = CtoR(clusts = list(matrix(clust[-j,], nrow = nrow(clust)-1),matrix(clust_init[c(df_res, nums[j]),], nrow = (1+length(df_res)))), only_bornes = FALSE, n_sample = n_sample, return_error = TRUE)
        weighted.mean(res_ctor[[3]]^2,res_ctor[[4]])
      }
     }
      best_idx = search[which.min(diff)] #we select the best element
      clust = matrix(clust[-best_idx,], ncol = ncol(clust)) #we remove it from the cluster
      df_res = c(df_res, nums[best_idx]) #we add it to the bin
      nums = nums[(nums %in% nums[best_idx])==FALSE]
    }
    return(df_res) #return the indices of the points that we put in the bin
  }
}

split_clusters = function(clusts, index_throw, pbin, n_sample, prop_search){ 
        l_bin = length(index_throw) #index_throw are the indices of the clusters to split
        num_throw = lapply(index_throw, function(i){f_delete(clusts[[i]], n = pbin*nrow(clusts[[i]]), n_sample = n_sample, prop_search = prop_search)}) #get the index of the points to place into the bins
        clusts_throw = lapply(1:l_bin, function(i){matrix(clusts[[index_throw[[i]]]][num_throw[[i]],], ncol = ncol(clusts[[1]]))}) #get the bins
        clusts_split = lapply(1:l_bin, function(i){matrix(clusts[[index_throw[[i]]]][-num_throw[[i]],], ncol = ncol(clusts[[1]]))}) #remove the elements from the clusters
        if(length(clusts) > l_bin){
          clusts_split = c(clusts_split, lapply(setdiff(1:length(clusts), index_throw), function(i){clusts[[i]]}))
        }
        clusts_split = c(clusts_split, clusts_throw) #get all the clusters
        if(0 %in% sapply(clusts_split, nrow)){
          clusts_split = lapply(which(sapply(clusts_split, nrow)>0), function(idx){clusts_split[[idx]]})
        } #keep only the clusters with elements
        return(clusts_split)
    }

############################### For merge

all_possibilities = function(vec){
  it = 1
  res_list = lapply(1:length(vec), function(i){combinations(sum(vec[i:length(vec)]), vec[i])})
  nb_it = prod(sapply(res_list, nrow))
  res_res = as.list(rep(0, nb_it))
  res_res = lapply(1:nb_it, function(y){as.list(rep(0, length(vec)))})
  l_index = expand.grid(lapply(res_list, function(x){1:nrow(x)}))
  res_res = lapply(1:nrow(l_index),function(x){lapply(1:length(l_index[x,]), function(y){res_list[[y]][l_index[x,y],]})})
  
  return(res_res)
}

arrange_possibilities = function(vec){
  indexes = all_possibilities(vec)
  nb_tot = sum(vec)
  for(i in 1:length(indexes)){
    all_clusts = 1:nb_tot
    for(k in 2:length(indexes[[i]])){
      all_clusts = setdiff(all_clusts, indexes[[i]][[k-1]])
      
      indexes[[i]][[k]] = all_clusts[indexes[[i]][[k]]]
    }
    
  }
  if(sum(duplicated(vec)) > 0){
    df_names = data.frame()
    for(i in 1:length(indexes)){
      ll = indexes[[i]]
      names(ll) = sapply(ll, function(x){paste(x, collapse = "")})
      indexes[[i]] = ll[order(names(ll))]
      df_names = rbind(df_names, as.numeric(sort(names(ll))))
      
    }
    indexes = lapply(which(duplicated(df_names) == FALSE), function(i){indexes[[i]]})
  }
  
  return(indexes)
}

merge_clusts = function(clusts, l, n_sample = 500){ 
  eg <- expand.grid(rep(list(1:length(clusts)), l))
  eg = t(apply(eg, 1, sort))
  eg = unique(eg[which(rowSums(eg)==length(clusts)),]) #eg contains the number of number of clusters that will be merged together. for ex, to go from 4 clusters to 2, it is either 2 merges of two clusters, or 1 cluster unchanges and 3 clusters merged together
  best_error = 10^5
  d = ncol(clusts[[1]])
  it = 0
  rec = list()
  for(row in 1:nrow(eg)){#for each row in eg, i.e. each poss
    indexes = arrange_possibilities(eg[row,]) #We get all the possibilities to merge these clusters. In the case of two merges of two clusters, it can be 1-2 and 3-4, or 1-3 and 2-4, or 1-4 and 2-3
    rec[[row]] = indexes
    for(i in 1:length(indexes)){ #for each merging configuration
      clusts_tilde = lapply(1:length(indexes[[i]]), function(x){do.call("rbind", lapply(indexes[[i]][[x]], function(j){clusts[[j]]}))}) #merge the clusters and store them in clusts_tilde
      res_ctor = CtoR(clusts = clusts_tilde,only_bornes = TRUE, n_sample = n_sample, return_error = TRUE) #Get the representatives and the associated clustering error
      rep = res_ctor[[1]]
      law = res_ctor[[2]]
      error = sqrt(weighted.mean(res_ctor[[3]]^2, res_ctor[[4]])) #clustering error
      
      it=it+1
      if(error < best_error){ #if this merge is the best, store it
        best_mixture = list(clusts = clusts_tilde, rep = rep, law = law, error = error)
        best_error = error
      }
    }
  }
  return(best_mixture) #return the best mixture
}


####################### Full algo


augmented_quanti = function(samp, rep, law = "unif", n_sample = 500, l_bin =2, vec_prop = c(0.4,0.3,0.2,0.1,0.05,0.025),it_lim = 10, threshold = 0.002,prop_search = 1, find_c_bool = TRUE, perturb_bool = TRUE,weights = rep(1/length(rep), length(rep))){
  write.csv(rbind(c(0,0)), "df_follow.csv", row.names = FALSE)
  mixture = list(clusts = NULL, rep = rep, law = law) #initial mixture
  record = list()
  for(pbin in vec_prop){ #for each epoch of perturbation intensity
    record[[paste(pbin)]] = as.list(rep(0,it_lim))
    it = 1
    dist_mixtures = 10^5
    while(it <= it_lim & dist_mixtures > threshold){ #while budget is spent or convergence
      df_follow = as.matrix(read.csv("df_follow.csv"))
      write.csv(rbind(df_follow,cbind(pbin, it)), "df_follow.csv", row.names = FALSE)
        
      if(!(pbin == vec_prop[1] & it == 1)){weights = sapply(mixture$clusts, nrow)/sum(sapply(mixture$clusts, nrow))} #weights for the findclusters step
      res_find_clusts = find_c(samp=samp, rep = mixture$rep, weights = weights, law = mixture$law, n = n_sample*10) #find clusters
      clusts = res_find_clusts[[1]] #get the clusters
      if(pbin == vec_prop[1] & it == 1){ #if initial iteration
         rep_bis = from_law_to_sample(rep, law = law, n_sample) #Associate samples to the representatives
         mixture_old = list(clusts = clusts, rep = rep, law = law) #initial mixture
         error_init = local_errors(clusts = clusts,rep = rep_bis,idx = 1:length(clusts)) #initial local errors
         error = sqrt(weighted.mean(error_init^2, sapply(clusts, nrow))) #initial quantization error
         record = c(list(init = list(c(mixture_old, list(error = error)))), record) #store error
      }
      res_ctor = CtoR(clusts, only_bornes = FALSE, n_sample = n_sample, return_error = TRUE) #ctor to get the local errors of the clusters
      errors = res_ctor[[3]] #local errors
      if(perturb_bool){ #perturb
        index_throw = order(errors, decreasing = T)[1:l_bin] #select the l_bin clusters with the highest local error for the split
        clusts_split = split_clusters(clusts = clusts, index_throw = index_throw, pbin = pbin, n_sample = n_sample, prop_search = prop_search) #split the selected clusters
        merge = merge_clusts(clusts_split, length(rep), n_sample = n_sample) #merge the clusters
        mixture = merge[1:3]
        error = merge[[4]]
      }
      else{
        mixture = list(clusts = clusts, rep = res_ctor[[1]], law = res_ctor[[2]])
        error = sqrt(weighted.mean(errors^2, res_ctor[[4]]))
      }
      record[[paste(pbin)]][[it]] = c(mixture, list(error = error)) #store the clusters before and after the split, after the merge, the representatives, and the quantization error
      dist_mixtures = wasserstein_between_mixtures(mixture_old, mixture) #check if the previous mixture and the new one are close enough
      mixture_old = mixture #update mixture
      it = it + 1 #new iteration
    }
  }
  return(record)
}


############################# For plots and analysis

find_best = function(record){ #from the record of all AQ computations, extract best mixture
  record = lapply(record, function(x){x[sapply(x, length) == 4]}) #get all steps 
  rr = lapply(record, function(x){sapply(x, function(y){y[[4]]})}) #get all quantization errors
  
  idx_1 = which.min(sapply(rr, min)) #get the epoch of pbin with the minimum quantization error
  idx_2 = which.min(rr[[idx_1]]) #get the iteration of this epoch with the minimum quantization error
  
  best = record[[idx_1]][[idx_2]] #get the optimal mixture
  return(best)
}

get_evol_error = function(record){ #from the record of all AQ computations, get all the quantization errors
  record = lapply(record, function(x){x[sapply(x, length) == 4]}) #get all steps
  rr = lapply(record, function(x){sapply(x, function(y){y[[4]]})}) #get all quantization errors
  return(as.numeric(unlist(rr)))
}

create_df_plot = function(list_unifs, list_clusts){ #create dataframe to plot the mixture
  order_clusts = order(sapply(list_clusts, nrow), decreasing = TRUE) #order the clusters using their probabilistic weights
  list_unifs = list_unifs[order_clusts] 
  list_clusts = list_clusts[order_clusts]
  weights = sapply(list_clusts, nrow)/sum(sapply(list_clusts, nrow)) #compute weights
  while(sum(round(weights,2)) != 1){ #just ensure the sum of the weights is equal to one
      delta = 1 - sum(round(weights,2))
      weights = weights + rep(delta/abs(delta)*0.001, l = length(weights)) + c(10^-4, 5*10^-5,0)
      }
  weights = round(weights,2)                         
  res = data.frame()
  for(i in 1:length(list_unifs)){ #for each representative
    for(borne in 1:2){ #for each bound of the uniform
      w = as.character(weights[i]) #get weight
      if(nchar(w) == 3){w = paste(w, "0", sep = "")} #good format for the weight
      res = rbind(res, c(w,  list_unifs[[i]][borne, ])) #concatenate lower and upper bounds of the representatives 
    }   
  }
  names_col = c("Num")
  for(i in 1:ncol(list_unifs[[1]])){
    names_col = c(names_col, paste("X", i, sep = "")) #names of the inputs
  }
  colnames(res) = names_col
  return(res)
}

plot_hybrid_distrib = function(df_plot, names_var,yat = seq(0,1,l=5), ylab = c("0","0.25","0.5","0.75","1"),ysize = 1.4,yname = "Quantiles", cex_legend = 1, inset = c(-0.18,0)){ #plot hybrid mixture
  par(mar=c(5.1, 4.1, 4.1, 7.1)+0.8, xpd=TRUE)
  nbpar2 = (nrow(df_plot)/2)%/%2
  vec = seq(-nbpar2/90,nbpar2/90,l=nrow(df_plot)/2)
  color_palette = rainbow(nrow(df_plot)/2)
  grads = seq(1/2-((ncol(df_plot)-1)%/%2)/8.3, 1/2+((ncol(df_plot)-1)%/%2)/8.3, l = ncol(df_plot)-1)
  ymax = max(1,max(as.numeric(as.matrix(df_plot[,2:ncol(df_plot)]))))
  ymin = min(0,min(as.numeric(as.matrix(df_plot[,2:ncol(df_plot)]))))
  plot(NA, xlim = c(grads[1]-1/10,grads[length(grads)]+1/10), ylim = c(ymin,ymax),yaxt = "n", xaxt = "n",ylab = yname, cex.axis = 1.4,cex.lab=1.8, xlab = "")
  axis(2, at = yat,labels=ylab,cex.axis = ysize)
  for(k in 1:(nrow(df_plot)%/%2)){ #for each representative
    df_sub = df_plot[(2*k-1):(2*k),] #get the bounds of the uniform or dirac
    for(i in 1:length(grads)){
      col = grads[i] 
      if(length(unique(df_sub[,i+1])) == 1){ #if lower bound equal to upper bound: dirac
        points(col + 1.3*vec[k], df_sub[1,i+1], pch = 17, col = color_palette[k], cex=1.4) #plot triangle
      }
      else{
        rect(col + 1.3*vec[k]-+0.005/2, df_sub[1,i+1], col+1.3*vec[k]+0.005/2, df_sub[2,i+1], col = color_palette[k]) #else: plot segment
      }
    }
  }
  legend("right", legend=(df_plot[seq(2, nrow(df_plot), by=2),1]), fill=color_palette, inset = inset,title="Weights",cex=cex_legend) #legend with weights

  axis(1, labels = names_var, at = grads, cex.axis = 1.6)
}

plot_clusts = function(clusts, legend = TRUE, colnames = c("x")){ #plot clusters samples
  clusts = lapply(1:length(clusts), function(x){cbind(clusts[[x]], x)}) #associate cluster with its number
  df_clusts = do.call(rbind, clusts)
  df_clusts[,ncol(df_clusts)] = as.factor(df_clusts[,ncol(df_clusts)]) #cluster number to factor
  colnames(df_clusts) = c(colnames, "Cluster") 
  if(ncol(df_clusts)==2){ #if 1D samples
    plot_return = ggplot(df_clusts) + geom_histogram(aes_string(x = colnames[1], fill = "Cluster", col = "Cluster"), alpha = 0.2)+ theme_bw()}
  else{
    plot_return = ggplot(df_clusts) + geom_point(aes_string(x = colnames[1], y = colnames[2],fill = "Cluster", col = "Cluster"))+ theme_bw()
  }
  if(!legend){plot_return = plot_return+theme(legend.position = "none") + xlab("") + ylab("") }
  return(plot_return)
}


