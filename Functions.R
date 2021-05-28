get_pairwise_dist<-function(dist_matrix,grouped_list_of_names,treatments){
  
## This function returns a data frame with the average distance of each sample to the other samples in its treatment group, as specified by the parameter grouped_list_of_names.
  
## Arguments:
## dist_matrix = dissimilarity matrix (e.g. as produced by phyloseq::distance)
## grouped_list_of_names = list of vectors, where each vector contains the names of samples all in the same treatment group
## treatments = vector of names for the groups in grouped_list_of_names. Must be the same length as grouped_list_of_names.

  
  # Initialize data frame for pairwise differences
  dispersion_df<-data.frame(matrix(nrow=0,ncol=3))
  
  # Iterate through treatment groups
  for (i in seq(1,length(treatments))){
    
    # Extract sample names belonging to each treatment
    names<-grouped_list_of_names[[i]]
    
    # Calculate average pairwise distance of each sample to other samples in treatment group
    for (sample1 in names){
      BC_distances<-c()
      for (sample2 in names[names!=sample1]){
        BC_distances<-c(BC_distances,dist_matrix[sample1,sample2])
      }
      dispersion_df<-rbind(dispersion_df,c(sample1,treatments[i],mean(BC_distances)))
    }
  }
  
  # After iterating all treatment groups, return data frame
  colnames(dispersion_df)<-c("Sample","Treatment","avg_BC")
  dispersion_df$avg_BC<-as.numeric(dispersion_df$avg_BC)
  return(dispersion_df)
}


## This function returns a vector of total edges (significant co-associations after multiple testing correction), total nodes, sub-categories broken down by positive vs. negative connections and bacterial vs. fungal connections, mean degree of bacterial taxa, mean degree of fungal taxa, and overall clustering coefficient.

## Arguments:
## corr = correlation matrix (e.g. produced by the function rcorr)
edge_metrics<-function(corr){
  
  ##total edges
  g1<-graph.adjacency(corr,weight=T,mode="undirected")
  g1<-simplify(g1)
  num_edges<-gsize(g1)
  num_nodes<-length(V(g1))
  
  #delete positive edges and recount
  corr_neg<-corr
  corr_neg[which(corr_neg>0)]=0
  corr_neg_edges<-gsize(graph.adjacency(corr_neg,weight=T,mode="undirected"))
  
  #delete negative edges and recount
  corr_pos<-corr
  corr_pos[which(corr_pos<0)]=0
  corr_pos_edges<-gsize(graph.adjacency(corr_pos,weight=T,mode="undirected"))
 
  ##bacteria-bacteria edges
  graph_names<-rownames(corr)
  bac_names<-graph_names[substr(graph_names,1,1)=="B"]
  B_only<-corr[bac_names,bac_names]
  gB_only<-graph.adjacency(B_only,weight=T,mode="undirected")
  bac_edges<-gsize(gB_only)
  
  #delete positive edges and recount
  B_neg<-B_only
  B_neg[which(B_only>0)]=0
  bac_neg_edges<-gsize(graph.adjacency(B_neg,weight=T,mode="undirected"))
  
  #delete negative edges and recount
  B_pos<-B_only
  B_pos[which(B_only<0)]=0
  bac_pos_edges<-gsize(graph.adjacency(B_pos,weight=T,mode="undirected"))
  
  ##fungi-fungi edges
  fun_names<-graph_names[substr(graph_names,1,1)=="F"]
  F_only<-corr[fun_names,fun_names]
  gF_only<-graph.adjacency(F_only,weight=T,mode="undirected")
  fun_edges<-gsize(gF_only)

  #delete positive edges and recount
  F_neg<-F_only
  F_neg[which(F_only>0)]=0
  fun_neg_edges<-gsize(graph.adjacency(F_neg,weight=T,mode="undirected"))
  
  #delete negative edges and recount
  F_pos<-F_only
  F_pos[which(F_only<0)]=0
  fun_pos_edges<-gsize(graph.adjacency(F_pos,weight=T,mode="undirected"))
  
  
  mean_bac_degree<-mean(degree(g1)[bac_names])
  mean_fun_degree<-mean(degree(g1)[fun_names])
  
  cc<-transitivity(g1)
  
  output<-c(num_edges,num_nodes,corr_pos_edges,corr_neg_edges,bac_edges,bac_pos_edges,bac_neg_edges,fun_edges,fun_pos_edges,fun_neg_edges,mean_bac_degree,mean_fun_degree,cc)
  return(output)
}

## This function returns a data frame of the degree, degree normalized to network size, and betweenness of each node in a network.
## Arguments:
## corr = correlation matrix (e.g. produced by the function rcorr)
node_metrics<-function(corr){
  g1<-graph.adjacency(corr,weight=T,mode="undirected")
  g1<-simplify(g1)
  
  graph_names<-rownames(corr)
  bac_names<-graph_names[substr(graph_names,1,1)=="B"]
  fun_names<-graph_names[substr(graph_names,1,1)=="F"]
  
  mean_bac_degree<-mean(degree(g1)[bac_names])
  mean_fun_degree<-mean(degree(g1)[fun_names])
  
  E(g1)$weight<-1
  output<-data.frame(degree(g1),degree(g1,normalized=T),betweenness(g1))
  names(output)=c("degree","degree_norm","betweenness")
  return(output)
}
