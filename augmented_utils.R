library(doParallel)
library(sBIC)
library(mclust)
library(multicool)
library(transport)
library(randtoolbox)
library(ggplot2)
library(gtools)

registerDoParallel(cores = 20)


# augmented_quanti = function(samp, rep, law = "unif", n_sample = 500, l_bin =2, vec_prop = c(0.4,0.3,0.2,0.1,0.05,0.025),it_lim = 10, nb_start = 10, threshold = 0.01,prop_search = 1, find_c_bool = TRUE, perturb_bool = TRUE){
#   rep = lapply(order(sapply(rep, function(x){x[1,1]})), function(g){rep[[g]]}) #sort the representatives
#   rep_bis = from_law_to_sample(rep, law = law, n_sample) #Associate samples to the representatives
#   record = list()
#   weights = rep(1/length(rep), length(rep)) #Give equal weights to all the representatives
#   res_find_clusts = find_c(samp=samp, rep = rep, weights = weights,  law = law) #find clusters
#   nums = res_find_clusts[[2]]
#   clusts = lapply(res_find_clusts[[1]], function(x){matrix(x,ncol = ncol(rep[[1]]))})
#   if(0 %in% sapply(clusts, nrow)){ #If one cluster is empty, give him its closest point in the sample
#     for(i0 in which(sapply(clusts, nrow)== 0)){
#       idx_min = which.min(Vectorize(function(i){dist2(samp[i,],rep_bis[[i0]])})(1:nrow(samp)))
#       nums[[i0]] = idx_min
#       
#       nums[setdiff(1:length(clusts),i0)] = lapply(setdiff(1:length(clusts),i0), function(x){nums[[x]][nums[[x]] != idx_min]})
#     }
#     clusts = lapply(nums, function(x){matrix(samp[x,], ncol = ncol(samp))})
#   }
#   error_init = local_errors(clusts = clusts,rep = rep_bis,idx = 1:length(clusts))
#   error_quant = weighted.mean(error_init^2, sapply(clusts, nrow))
#   record[["init"]] = list(list(NULL,NULL, clusts, rep, error_quant))
#   for(pbin in vec_prop){
#     record[[paste(pbin)]] = as.list(rep(0,it_lim))
#     print(pbin)
#     rep_old = rep
#     for(it in 1:it_lim){
#       print(it)
#       if(find_c_bool){
#         res_find_clusts = find_c(samp=samp, rep = rep, weights = weights, law = law) #find clusters
#         clusts = lapply(res_find_clusts[[1]], function(x){matrix(x,ncol = ncol(rep[[1]]))}) #get the clusters
#         nums = res_find_clusts[[2]] #get the indexes of each point of the clusters
#         if(0 %in% sapply(clusts, nrow)){ #If one cluster is empty, give him its closest point in the sample
#           for(i0 in which(sapply(clusts, nrow)== 0)){
#             idx_min = which.min(Vectorize(function(i){dist2(samp[i,],rep_bis[[i0]])})(1:nrow(samp)))
#             nums[[i0]] = idx_min
#             
#             nums[setdiff(1:length(clusts),i0)] = lapply(setdiff(1:length(clusts),i0), function(x){nums[[x]][nums[[x]] != idx_min]})
#           }
#           clusts = lapply(nums, function(x){matrix(samp[x,], ncol = ncol(samp))})
#         }
#       }
#       errors = local_errors(clusts, CtoR(clusts = clusts, only_bornes = FALSE, n_sample = n_sample, law = law), idx = 1:length(clusts)) #get the local errors
#       error = weighted.mean(errors, sapply(clusts, nrow)) #get the quantization error
#       former_clusts = clusts 
#       clusts_split = NULL
#       if(perturb_bool){ #perturb
#         index_throw = order(errors, decreasing = T)[1:l_bin] #select the l_bin clusters with the highest local error for the split
#         num_throw = lapply(index_throw, function(i){f_delete(clusts[[i]], n = pbin*nrow(clusts[[i]]), n_sample = n_sample,law = law, prop_search = prop_search)}) #get the index of the points to place into the bins
#         clusts_throw = lapply(1:l_bin, function(i){matrix(clusts[[index_throw[[i]]]][num_throw[[i]],], ncol = ncol(samp))}) #get the bins
#         clusts_split = lapply(1:l_bin, function(i){matrix(clusts[[index_throw[[i]]]][-num_throw[[i]],], ncol = ncol(samp))}) #remove the elements from the clusters
#         if(length(errors) > l_bin){
#           clusts_split = c(clusts_split, lapply(setdiff(1:length(clusts), index_throw), function(i){clusts[[i]]}))
#         }
#         
#         clusts_split = c(clusts_split, clusts_throw) #get all the clusters
#         if(0 %in% sapply(clusts_split, nrow)){
#           clusts_split = lapply(which(sapply(clusts_split, nrow)>0), function(idx){clusts_split[[idx]]})
#         } #keep only the clusters with elements
#         merge = merge_clusts(clusts_split, length(rep),law = law, n_sample = n_sample) #merge
#         clusts = merge[[2]]
#         rep = merge[[1]]
#         error = merge[[4]]
#       }
#       else{rep = CtoR(clusts = clusts, only_bornes = TRUE, law = law)}
#       
#       record[[paste(pbin)]][[it]] = list(former_clusts, clusts_split, clusts, rep, error) #store the clusters before and after the split, after the merge, the representatives, and the quantization error
#       nums = get_nums(clusts, samp)
#       error_quant = sqrt(error)
#       idx_order = order(sapply(rep, function(x){x[1,1]}))
#       rep = lapply(idx_order, function(g){rep[[g]]}) #sort the representatives
#       weights = sapply(clusts, nrow)/nrow(samp)
#       weights = weights[idx_order] #sort the weights
#       # for(i in 1:length(rep)){
#       #   if(sum(is.na(rep[[i]])) == 4){
#       #     rep[[i]] = matrix(runif(length(rep[[1]])), ncol=2)
#       #   }
#       # }
#       rep_bis = from_law_to_sample(rep, law = law, n_sample)
#       gc()
#       if(sum(Vectorize(function(m){sum(abs(rep[[m]]- rep_old[[m]]))})(1:length(rep))) < threshold){break} #If the change of representatives is lower than the threshold, stop the algorithm
#       rep_old = rep
#       
#       
#       #print(sqrt(error))
#     }
#   }
#   return(record)
# }
dist2 = function(x,y){sqrt(sum((x-y)^2))}

num_func = function(x, rep){
  dists = sapply(rep, function(y){dist2(as.numeric(y), x)})
  return(c(min(dists),which.min(dists)))
}

func_prep_clusts = function(clust, nrow_rep){
  if(length(clust)==0){return(clusts)}
  else{return(clust[rep(1:nrow(clust), nrow_rep%/% nrow(clust)),])}
} 


local_errors = function(clusts, rep, idx){
  error= c()
  for(i in idx){
    if(length(clusts[[i]]) == 0 | length(rep[[i]]) == 0){error = c(error,0)}
    else if(ncol(rep[[i]]) == 1){error= c(error, wasserstein1d(clusts[[i]], rep[[i]],p=2))}
    else{
      clusts[[i]] = func_prep_clusts(clusts[[i]], nrow(rep[[i]]))
      rep[[i]] = rep[[i]][1:nrow(clusts[[i]]),]
      error= c(error, transport::wasserstein(pp(clusts[[i]]), pp(rep[[i]]), p=2))}
  }
  return(error)
} 


vec = seq(10^-10,1-10^-10,l=5000)
vec_erf_inv = Vectorize(function(q){erfinv(2*q-1)})(vec)
mean_vec = mean(vec_erf_inv^2)

try_unif = function(clusts, only_bornes = TRUE, n_sample = NULL){
  l = length(clusts)
  d = ncol(clusts[[1]])
  x_tilde = list()
  for(i in 1:l){
    x_tilde[[i]] = matrix(0, ncol = ncol(clusts[[1]]), nrow = 2)
    for(k in 1:d){
      quantiles_vec = quantile(clusts[[i]][,k], vec)
      borne_sup = mean(quantiles_vec*(6*vec-2))
      borne_inf = mean(quantiles_vec*(-6*vec+4))
      x_tilde[[i]][,k] = c(borne_inf, borne_sup)
    }
  }
  if(only_bornes == FALSE){
    x_tilde = from_law_to_sample(x_tilde, law = as.list(rep("unif",l)), n = n_sample)
  }
  return(x_tilde)
}



try_norm = function(clusts, only_bornes = TRUE, n_sample = NULL){
  x_tilde = list()
  for(i in 1:length(clusts)){
    x_tilde[[i]] = matrix(0, ncol = ncol(clusts[[1]]), nrow = 2)
    for(k in 1:ncol(clusts[[1]])){
      clust = clusts[[i]][,k]
      mu = mean(clust)
      quantiles_vec = quantile(clust, vec)
      sigma = mean(quantiles_vec*vec_erf_inv)/sqrt(2)/mean_vec
      x_tilde[[i]][,k] = c(mu, sigma)
    }
  }
  if(only_bornes == FALSE){
    x_tilde = from_law_to_sample(x_tilde, as.list(rep("normal",length(clusts))), n = n_sample)
  }
  return(x_tilde)
}


regularize_samp = function(samp_clust, rep){
  return(samp_clust[(nrow(samp_clust)-nrow(rep)+2):(nrow(samp_clust)),])
} 


find_c = function(samp, rep, weights,law, n_sample = 10^5){
  rep_large = from_law_to_sample_find_c(rep=rep, weights = weights, law=law, n = n_sample)
  nums = lapply(1:length(rep), function(x){c()})
  clust_number = foreach(i=1:nrow(samp),.combine='c')%dopar%{
    which.min(sapply(rep_large, function(gamm){min(apply(gamm, 1, function(y){dist2(as.numeric(samp[i,]),y)}))}))
  }
  nums = lapply(1:length(rep), function(j){which(clust_number == j)})
  return(list(clusters = lapply(nums, function(x){samp[x,]}),nums = nums))
}

# 
# find_clusters = function(samp, rep, seed = NULL, error_quanti = NULL, nums = NULL){
#   nums_init = nums
#   record = list()
#   error_best = 10^10
#   samp_clust = rep
#   if(!is.null(seed)){set.seed(seed)}
#   idx_mix = sample(1:nrow(samp))
#   samp_mix = matrix(samp[idx_mix,], nrow = nrow(samp))
#   nums = as.list(rep(0,length(rep)))
#   l = length(rep)
#   err_clust_fige = rep(0,l)
#   for(i in 1:nrow(samp_mix)){
#     x = samp_mix[i,]
#     err = c()
#     new_error = c()
#     for(clust in 1:length(rep)){
#       samp_prov = samp_clust
#       if(ncol(samp) > 1 & nrow(rep[[clust]])>1){samp_prov[[clust]] = regularize_samp(samp_prov[[clust]], rep[[clust]])}
#       samp_prov[[clust]] = rbind(samp_prov[[clust]], x)
#       err_clust = err_clust_fige
#       new_error = c(new_error,local_errors(list(samp_prov[[clust]]), list(rep[[clust]]), 1))
#       err_clust[clust] = new_error[clust]
#       err = c(err, sum(err_clust^2))
#     }
#     samp_clust[[which.min(err)]] = rbind(samp_clust[[which.min(err)]], x)
#     nums[[which.min(err)]] = c(nums[[which.min(err)]], idx_mix[i])
#     err_clust_fige[which.min(err)] = new_error[which.min(err)]
#   }
#   nums = lapply(nums, function(x){x[-1]})
#   samp_clust = lapply(nums, function(x){matrix(samp[x,], nrow = length(x))})
#   errors = local_errors(samp_clust, rep, 1:l)
#   error_tot = weighted.mean(errors^2, sapply(nums, length))
#   if(!is.null(error_quanti)){
#     if(sqrt(error_tot)> error_quanti){
#       nums = nums_init
#       error = quanti_error
#     }
#   }
#  
#   return(list(clusters = lapply(nums, function(x){samp[x,]}),nums = nums, error = error_tot))
# }
# 
# find_clusters_multi_start = function(samp, rep, nb_start = 10, seed = NULL, error_quanti = NULL, nums = NULL){
#   if(!is.null(seed)){set.seed(seed,"L'Ecuyer-CMRG")}
#   list_clusters = foreach(test = 1:nb_start) %dopar% {
#     find_clusters(samp, rep)
#   }
#   best = which.min(sapply(list_clusters, function(x){x[["error"]]}))
#   return(list_clusters[[best]])
# }

create_df_plot = function(list_unifs, lists_clusts){
  res = data.frame()
  for(i in 1:length(list_unifs)){
    for(borne in 1:2){
      weight = round(nrow(lists_clusts[[i]])/sum(sapply(lists_clusts, nrow)), 3)
      res = rbind(res, c(paste(i," :",weight), list_unifs[[i]][borne, ]))
    }
    
  }
  names_col = c("Num")
  for(i in 1:ncol(list_unifs[[1]])){
    names_col = c(names_col, paste("X", i, sep = ""))
  }
  colnames(res) = names_col
  return(res)
}
plot_hybrid_distrib = function(df_plot, names_var,yat = seq(0,1,l=5), ylab = c("0","0.25","0.5","0.75","1"),ysize = 1.2,yname = "inputs rank"){
  par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
  nbpar2 = (nrow(df_plot)/2)%/%2
  vec = seq(-nbpar2/90,nbpar2/90,l=nrow(df_plot)/2)
  color_palette = rainbow(nrow(df_plot)/2)
  grads = seq(1/2-((ncol(df_plot)-1)%/%2)/8.3, 1/2+((ncol(df_plot)-1)%/%2)/8.3, l = ncol(df_plot)-1)
  ymax = max(1,max(as.numeric(as.matrix(df_plot[,2:ncol(df_plot)]))))
  ymin = min(0,min(as.numeric(as.matrix(df_plot[,2:ncol(df_plot)]))))
  plot(NA, xlim = c(grads[1]-1/10,grads[length(grads)]+1/10), ylim = c(ymin,ymax),yaxt = "n", xaxt = "n",ylab = yname, cex.axis = 1.4,cex.lab=1.5)
  axis(2, at = yat,labels=ylab,cex.axis = ysize)
  for(k in 1:length(unique(df_plot[,1]))){
    x = unique(df_plot[,1])[k]
    df_sub = df_plot[df_plot[,1]==x,]
    for(i in 1:length(grads)){
      col = grads[i]
      if(length(unique(df_sub[,i+1])) == 1){
        points(col + vec[k], df_sub[1,i+1], pch = 17, col = color_palette[k])
      }
      else{
        rect(col + vec[k], df_sub[1,i+1], col+vec[k]+0.005, df_sub[2,i+1], col = color_palette[k])
      }
    }
  }
  axis(1, labels = names_var, at = grads, cex.axis = 1.3)
  legend("right", legend=sapply(unique(df_plot[,1]), function(x){substr(x,5,nchar(x))}), fill=color_palette, inset = c(-0.18,0),title="Weights",cex=1.25)
}


plot_sample = function(clusts, legend = TRUE, colnames = c("X"), rpz){
  clusts = lapply(1:length(clusts), function(x){cbind(clusts[[x]], x)})
  df_clusts = data.frame()
  for(i in 1:length(clusts)){
    df_clusts = rbind(df_clusts, clusts[[i]])
  }
  df_clusts[1:nrow(clusts[[1]]), 2] = rpz[1]
  df_clusts[(nrow(clusts[[1]])+1):nrow(df_clusts), 2] = rpz[2]
  df_clusts[,ncol(df_clusts)] = as.factor(df_clusts[,ncol(df_clusts)])
  colnames(df_clusts) = c(colnames, "True_representatives")
  fd <<-df_clusts
  plot_return = ggplot(df_clusts) + geom_histogram(aes_string(x = colnames[1], fill = "True_representatives", col = "True_representatives"), alpha = 0.2)+ theme_bw() + labs(fill = "True representatives") + guides(color = FALSE)
  return(plot_return)
}


from_law_to_sample_find_c = function(rep, n, weights, law){
  #set.seed(2)
  thresh = c(1,as.integer(round(cumsum(weights),4)*n))
  rep_bis = list()
  for(x in 1:length(rep)){
    if(law[[x]] == "unif"){
      rep_bis[[x]] = as.matrix(runif(n*ncol(rep[[1]])), ncol = ncol(rep[[1]]))
      rep_bis[[x]] = as.matrix(foreach(k = 1:ncol(rep[[x]]), .combine = "cbind")%do% {rep_bis[[x]][,k]*(rep[[x]][2,k]-rep[[x]][1,k])+rep[[x]][1,k]})
    }
    else if(law[[x]] == "normal"){
      rep_bis[[x]] = as.matrix(rnorm(n*ncol(rep[[1]])), ncol = ncol(rep[[1]]))
      rep_bis[[x]] = as.matrix(foreach(k = 1:ncol(rep[[x]]), .combine = "cbind")%do% {rep_bis[[x]][,k]*rep[[x]][2,k]+rep[[x]][1,k]})
    }
  }
  rep_bis = lapply(1:length(rep_bis),function(x){as.matrix(rep_bis[[x]][(thresh[x]:thresh[x+1]),])})
  return(rep_bis)
}

from_law_to_sample = function(rep, law = "unif", n){
  rep_bis = list()
  for(j in 1:length(rep)){
    if(law[[j]] == "unif"){
      sob = matrix(sobol(n, dim = ncol(rep[[1]])), ncol = ncol(rep[[1]]))
      rep_bis[[j]] = foreach(k = 1:ncol(rep[[j]]), .combine = "cbind")%do% {sob[,k]*(rep[[j]][2,k]-rep[[j]][1,k])+rep[[j]][1,k]}
    }
    else if(law[[j]] == "normal"){
      sob = matrix(sobol(n, normal = TRUE), ncol = ncol(rep[[1]]))
      rep_bis[[j]] = foreach(k = 1:ncol(rep[[j]]), .combine = "cbind")%do% {sob[,k]*rep[[j]][2,k]+rep[[j]][1,k]}
    }
  }
  if(ncol(sob)==1){rep_bis = lapply(rep_bis, function(x){matrix(x, nrow = length(x))})}
  return(rep_bis)
}


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
  eg = unique(eg[which(rowSums(eg)==length(clusts)),])
  best_error = 10^5
  d = ncol(clusts[[1]])
  it = 0
  rec = list()
  for(row in 1:nrow(eg)){
    indexes = arrange_possibilities(eg[row,])
    rec[[row]] = indexes
    for(i in 1:length(indexes)){
      clusts_tilde = lapply(1:length(indexes[[i]]), function(x){do.call("rbind", lapply(indexes[[i]][[x]], function(j){clusts[[j]]}))})
      rep = CtoR(clusts = clusts_tilde)
      rep_bis = from_law_to_sample(rep = rep[[1]], law = rep[[2]], n = n_sample)
      error = local_errors(clusts_tilde, rep_bis, idx = 1:length(rep_bis))
      error = weighted.mean(error^2, sapply(clusts_tilde, nrow))
      it=it+1
      if(error < best_error){
        best_error = error
        best_law = rep
        best_clusts = clusts_tilde
        best_merge = indexes[[i]]
      }
    }
  }
  return(list(best_law, best_clusts, best_merge, best_error))
}

f_delete = function(clust, n, n_sample, prop_search = 1){ #f_delete identify the index of the elements to place into the bin
  if(nrow(clust)==1){return(1)}
  else{
    df_res = c() #df_res will store the indexes
    nums = 1:nrow(clust) 
    clust_init = clust
    while(length(df_res) < n){ #will the length of the bin is lower than the objective
      if(prop_search==1){search=1:nrow(clust)} #prop search is the proportion of candidate elements at each iteration
      else{search = sort(sample(1:nrow(clust), size = as.integer(round(prop_search*nrow(clust)))))}
      diff = foreach(j=search, .combine = "c")%dopar%{
        errors = local_errors(list(matrix(clust[-j,], ncol = ncol(clust_init)), matrix(clust_init[c(df_res, nums[j]),], ncol = ncol(clust_init))), CtoR(clusts = list(matrix(clust[-j,], nrow = nrow(clust)-1),matrix(clust_init[c(df_res, nums[j]),], nrow = (1+length(df_res)))), only_bornes = FALSE, n_sample = n_sample)[[1]], idx = 1:2)
        weighted.mean(errors^2, c(nrow(clust)-1, length(df_res)+1))} #for each element we evaluate its adding to the bin
      best_idx = search[which.min(diff)] #we select the best element
      clust = matrix(clust[-best_idx,], ncol = ncol(clust)) #we remove it from the cluster
      df_res = c(df_res, nums[best_idx]) #we add it to the bin
      nums = nums[(nums %in% nums[best_idx])==FALSE]
    }
    return(df_res)
  }
}

plot_clusts = function(clusts, legend = TRUE, colnames = c("x")){
  clusts = lapply(1:length(clusts), function(x){cbind(clusts[[x]], x)})
  df_clusts = data.frame()
  for(i in 1:length(clusts)){
    df_clusts = rbind(df_clusts, clusts[[i]])
  }
  df_clusts[,ncol(df_clusts)] = as.factor(df_clusts[,ncol(df_clusts)])
  colnames(df_clusts) = c(colnames, "Cluster")
  if(ncol(df_clusts)==2){
    plot_return = ggplot(df_clusts) + geom_histogram(aes_string(x = colnames[1], fill = "Cluster", col = "Cluster"), alpha = 0.2)+ theme_bw()}
  else{
    plot_return = ggplot(df_clusts) + geom_point(aes_string(x = colnames[1], y = colnames[2],fill = "Cluster", col = "Cluster"))+ theme_bw()
  }
  if(!legend){plot_return = plot_return+theme(legend.position = "none") + xlab("") + ylab("") }
  return(plot_return)
}

augmented_quanti = function(samp, rep, law = "unif", n_sample = 500, l_bin =2, vec_prop = c(0.4,0.3,0.2,0.1,0.05,0.025),it_lim = 10, threshold = 0.01,prop_search = 1, find_c_bool = TRUE, perturb_bool = TRUE,weights = rep(1/length(rep), length(rep))){
  #rep = lapply(order(sapply(rep, function(x){x[1,1]})), function(g){rep[[g]]}) #sort the representatives
  rep_bis = from_law_to_sample(rep, law = law, n_sample) #Associate samples to the representatives
  record = list()
  res_find_clusts = find_c(samp=samp, rep = rep, weights = weights,  law = law) #find clusters
  nums = res_find_clusts[[2]]
  clusts = lapply(res_find_clusts[[1]], function(x){matrix(x,ncol = ncol(rep[[1]]))})
  if(0 %in% sapply(clusts, nrow)){ #If one cluster is empty, give him its closest point in the sample
    for(i0 in which(sapply(clusts, nrow)== 0)){
      idx_min = which.min(Vectorize(function(i){dist2(samp[i,],rep_bis[[i0]])})(1:nrow(samp)))
      nums[[i0]] = idx_min
      
      nums[setdiff(1:length(clusts),i0)] = lapply(setdiff(1:length(clusts),i0), function(x){nums[[x]][nums[[x]] != idx_min]})
    }
    clusts = lapply(nums, function(x){matrix(samp[x,], ncol = ncol(samp))})
  }
  error_init = local_errors(clusts = clusts,rep = rep_bis,idx = 1:length(clusts))
  error_quant = weighted.mean(error_init^2, sapply(clusts, nrow))
  record[["init"]] = list(list(NULL,NULL, clusts, list(rep,law), error_quant))
  for(pbin in vec_prop){
    record[[paste(pbin)]] = as.list(rep(0,it_lim))
    rep_old = rep
    for(it in 1:it_lim){
      #print(it)
      #print(list(rep, law))
      if(find_c_bool){
        res_find_clusts = find_c(samp=samp, rep = rep, weights = weights, law = law) #find clusters
        clusts = lapply(res_find_clusts[[1]], function(x){matrix(x,ncol = ncol(rep[[1]]))}) #get the clusters
        nums = res_find_clusts[[2]] #get the indexes of each point of the clusters
        if(0 %in% sapply(clusts, nrow)){ #If one cluster is empty, give him its closest point in the sample
          for(i0 in which(sapply(clusts, nrow)== 0)){
            idx_min = which.min(Vectorize(function(i){dist2(samp[i,],rep_bis[[i0]])})(1:nrow(samp)))
            nums[[i0]] = idx_min
            
            nums[setdiff(1:length(clusts),i0)] = lapply(setdiff(1:length(clusts),i0), function(x){nums[[x]][nums[[x]] != idx_min]})
          }
          clusts = lapply(nums, function(x){matrix(samp[x,], ncol = ncol(samp))})
        }
      }
      errors = local_errors(clusts, CtoR(clusts = clusts, only_bornes = FALSE, n_sample = n_sample)[[1]], idx = 1:length(clusts)) #get the local errors
      error = weighted.mean(errors, sapply(clusts, nrow)) #get the quantization error
      former_clusts = clusts 
      clusts_split = NULL
      if(perturb_bool){ #perturb
        index_throw = order(errors, decreasing = T)[1:l_bin] #select the l_bin clusters with the highest local error for the split
        num_throw = lapply(index_throw, function(i){f_delete(clusts[[i]], n = pbin*nrow(clusts[[i]]), n_sample = n_sample, prop_search = prop_search)}) #get the index of the points to place into the bins
        clusts_throw = lapply(1:l_bin, function(i){matrix(clusts[[index_throw[[i]]]][num_throw[[i]],], ncol = ncol(samp))}) #get the bins
        clusts_split = lapply(1:l_bin, function(i){matrix(clusts[[index_throw[[i]]]][-num_throw[[i]],], ncol = ncol(samp))}) #remove the elements from the clusters
        if(length(errors) > l_bin){
          clusts_split = c(clusts_split, lapply(setdiff(1:length(clusts), index_throw), function(i){clusts[[i]]}))
        }
        
        clusts_split = c(clusts_split, clusts_throw) #get all the clusters
        if(0 %in% sapply(clusts_split, nrow)){
          clusts_split = lapply(which(sapply(clusts_split, nrow)>0), function(idx){clusts_split[[idx]]})
        } #keep only the clusters with elements
        merge = merge_clusts(clusts_split, length(rep), n_sample = n_sample) #merge
        clusts = merge[[2]]
        rep = merge[[1]][[1]]
        law = merge[[1]][[2]]
        error = merge[[4]]
      }
      else{
        rep_and_law = CtoR(clusts = clusts, only_bornes = TRUE, n_sample = n_sample)
        rep = rep_and_law[[1]]
        law = rep_and_law[[2]]
      }

      record[[paste(pbin)]][[it]] = list(former_clusts, clusts_split, clusts, list(rep, law), error) #store the clusters before and after the split, after the merge, the representatives, and the quantization error
      nums = get_nums(clusts, samp)
      error_quant = sqrt(error)
      idx_order = order(sapply(rep, function(x){x[1,1]}))
      rep = lapply(idx_order, function(g){rep[[g]]}) #sort the representatives
      weights = sapply(clusts, nrow)/nrow(samp)
      weights = weights[idx_order]#sort the weights
      law = lapply(idx_order, function(g){law[[g]]})
      # for(i in 1:length(rep)){
      #   if(sum(is.na(rep[[i]])) == 4){
      #     rep[[i]] = matrix(runif(length(rep[[1]])), ncol=2)
      #   }
      # }
      rep_bis = from_law_to_sample(rep, law = law, n_sample)
      gc()
      if(sum(Vectorize(function(m){sum(abs(rep[[m]]- rep_old[[m]]))})(1:length(rep))) < threshold){break} #If the change of representatives is lower than the threshold, stop the algorithm
      rep_old = rep
  

      #print(sqrt(error))
    }
  }
  return(record)
}

find_best = function(record){
  record = lapply(record, function(x){x[sapply(x, length) == 5]})
  rr = lapply(record, function(x){sapply(x, function(y){y[[5]]})})
  
  idx_1 = which.min(sapply(rr, min))
  idx_2 = which.min(rr[[idx_1]])
  
  best = record[[idx_1]][[idx_2]]
  return(best)
}

get_evol_error = function(record){
  record = lapply(record, function(x){x[sapply(x, length) == 5]})
  rr = lapply(record, function(x){sapply(x, function(y){y[[5]]})})
  return(sqrt(unlist(rr)))
}

get_nums = function(clusts, samp){
  nums = lapply(1:length(clusts), function(x){rep(0, length(clusts[[x]]))})
  for(j in 1:length(clusts)){
    for(i in 1:nrow(clusts[[j]])){
      index = which(apply(samp, 1, function(x){sum(abs(x-clusts[[j]][i,]))})==0)[[1]]
      samp[index,] = -10^10
      nums[[j]][i] = index
    }
  }
  return(nums)
}

# 
# get_mixture = function(rep, clusts = NULL, weights = NULL, law, n){
#   if(is.null(weights)){weights = Vectorize(function(j){nrow(clusts[[j]])/sum(sapply(clusts, nrow))})(1:length(clusts))}
#   puiss = log(n)/log(10)
#   res = data.frame(rep(0,n))
#   for(i in 1:ncol(rep[[1]])){
#     res_i = c()
#     for(j in 1:length(rep)){
#       if(law == "unif"){
#         res_i = c(res_i,runif(as.integer(round(round(weights[j],puiss)*n)),rep[[j]][1,i], rep[[j]][2,i]))
#       }
#       if(law == "normal"){
#         res_i = c(res_i,rnorm(as.integer(round(round(weights[j],puiss)*n)),rep[[j]][1,i], rep[[j]][2,i]))
#       }
#     }
#     res = cbind(res, res_i)
#   }
#   return(res[,2:ncol(res)])
# }

get_mixture = function(rep, clusts = NULL, weights = NULL, law, n){
  if(is.null(weights)){weights = Vectorize(function(j){nrow(clusts[[j]])/sum(sapply(clusts, nrow))})(1:length(clusts))}
  puiss = log(n)/log(10)
  res = data.frame(rep(0,n))
  for(i in 1:ncol(rep[[1]])){
    res_i = c()
    for(j in 1:length(rep)){
      if(law[[j]] == "unif"){
        res_i = c(res_i,runif(as.integer(round(round(weights[j],puiss)*n)),rep[[j]][1,i], rep[[j]][2,i]))
      }
      if(law[[j]] == "normal"){
        res_i = c(res_i,rnorm(as.integer(round(round(weights[j],puiss)*n)),rep[[j]][1,i], rep[[j]][2,i]))
      }
    }
    res = cbind(res, res_i)
  }
  return(res[,2:ncol(res)])
}

