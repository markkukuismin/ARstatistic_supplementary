
findCorrelation_length = function(D, thr = NULL){
 
  if(is.null(thr)){
   
    thr = abs(D[upper.tri(D)])
    
    thr = unique(thr)
     
  }
  
  thr = sort(thr)
   
  nmb_of_high_cor_features = rep(0, length(thr))
  
  for(i in 1:length(thr)){
    
    nmb_of_high_cor_features[i] = length(
      findCorrelation(D, cutoff = thr[i])
      )
    
  }
  
  df = data.frame(feature_nmb = nmb_of_high_cor_features,
                  thr = thr)
  
  df = dplyr::distinct(df, feature_nmb, .keep_all = TRUE)
  
  return(df)
  
}

