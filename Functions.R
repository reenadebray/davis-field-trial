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
