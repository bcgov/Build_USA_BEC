
# This function is used to calculate accuracy metrics for cv and test estimates

# This currently includes: 
# accuracy, mcc, sensitity, specificity, precision, recall, fmeas, kappa
# aspatial_acc, aspatial_meanacc
# spat_p (same as accuracy), spat_fp (fuzzy), mcc_fp

# for test datasets the following are also generated: 
# spat_pa, spat_fpa,

library(yardstick)

acc_metrix <- function(data){
  acc_bal <- data %>% bal_accuracy(BGC, .pred_class)
  ppv <- data %>% ppv(BGC, .pred_class) # positive predictive value
  precision <- data %>% precision(BGC, .pred_class) # precision 
  recall <- data %>% recall(BGC, .pred_class) # recall
  kap <- data %>% kap(BGC, .pred_class) # kappa
  fmean <- data %>% f_meas(BGC, .pred_class) # f means
  mcc <- data %>%  mcc(BGC, .pred_class)
  sens <- data %>%  sens(BGC, .pred_class)
  spec <- data %>%  spec(BGC, .pred_class)
  acc <- data %>% accuracy(BGC, .pred_class)
  jind <- data %>% j_index(BGC, .pred_class) 
  ###some aspatial metrics
  #data <- test.pred
  aspatial_pred <- data  %>% dplyr::select(.pred_class) %>% group_by(.pred_class) %>% mutate(pred.ratio = n()) %>%ungroup() %>% distinct()
  aspatial_BGC <- data  %>% dplyr::select(BGC) %>% group_by(BGC) %>% mutate(targ.ratio = n())%>%ungroup() %>% distinct()    
  aspatial_sum <- full_join(aspatial_BGC, aspatial_pred, by = c("BGC" = ".pred_class")) %>% mutate_if(is.integer, funs(replace_na(., 0))) %>% rowwise() %>% mutate(Min = min(targ.ratio, pred.ratio))
  .estimate <- colSums(aspatial_sum[,4])/colSums(aspatial_sum[,2])
  .metric= "aspatial_acc"
  .estimator= "aspatial"
  aspatial_acc <- data.frame(.metric, .estimator, .estimate) 
  
  aspatial_sum <- aspatial_sum %>% mutate(unit_pos = Min/targ.ratio)
  mean_acc <- colMeans(aspatial_sum[5])
  .estimate <- mean_acc
  .metric= "aspatial_meanacc"
  .estimator= "aspatial"
  aspatial_meanacc <- data.frame(.metric, .estimator, .estimate) 
  final_metrics <- bind_rows (acc, mcc, aspatial_acc, aspatial_meanacc, acc_bal,ppv, precision, recall, kap, fmean, sens, spec, jind) %>% mutate_if(is.character, as.factor)
  ##______________ add alt-call metrics
  ##______________ add in fuzzy-call metrics     
  ##______________  add in metrics by map unit aspatial
  #unit_metrics
}
