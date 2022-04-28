source('/Users/huongvu/Desktop/PCR_Desktop/Scripts/fitting_rn_curves/find_important_cycles.R')

calc_a_err = function(rn_df, sigmoid_df,param_args){
  # we compare the difference in y between the given cycles in param args
  
  # param_args = matching_points = c(1,2)
  
  rn_range = rn_df[param_args[2],]$normalized - rn_df[param_args[1],]$normalized
  sigmoid_range = sigmoid_df[param_args[2],]$y - sigmoid_df[param_args[1],]$y
  err = (sigmoid_range - rn_range)/rn_range
  # err = sigmoid_range - rn_range
  return(err)
}

calc_r_err = function(rn_df, sigmoid_df, param_args){
  # using convention that err > 0, step back and err < 0, step forward
  # for r, when sigmoid curve has a bigger slope, we want to reduce r
  # we will calculate the error opposite to other param
  
  # param_args = c(inflection_cycle, slope_range)
  
  # for param r which controls the slope, we tune r to have the same slope as 
  # the slope of rn curve between cycle (+- slope range) of inflection cycle
  
  start_cycle = max(param_args[['inflection_cycle']] - param_args[['slope_range']],
                    min(rn_df$cycle_no))
  end_cycle = min(param_args[['inflection_cycle']] + param_args[['slope_range']],
                  max(rn_df$cycle_no))
  rn_m = lm(normalized ~ cycle_no, data = rn_df[start_cycle:end_cycle,])$coefficients[[2]]
  sigmoid_m = lm(y ~ x, sigmoid_df[start_cycle:end_cycle,])$coefficient[[2]]
  err = (rn_m - sigmoid_m)/rn_m
  # err = rn_m - sigmoid_m
  return(err)
}

calc_b_err = function(rn_df, sigmoid_df, param_args){
  # for parameter b which also kind of controls slope of the curve, but we 
  # will control param b by check the error from inflection point to the rest of 
  # of the curve
  
  # param_args = c(inflection_cycle)
  
  # inflection_cycle = param_args[['inflection_cycle']]
  # max_cycle = max(rn_df$cycle_no)
  # diff = rn_df$normalized[inflection_cycle] - sigmoid_df$y[inflection_cycle]
  # sigmoid_df$y = sigmoid_df$y + diff
  # # calculate slopes
  # rn_m = lm(normalized ~ cycle_no,rn_df[(inflection_cycle + 3):max_cycle,])$coefficients[[2]]
  # sigmoid_m = lm(y ~ x, sigmoid_df[(inflection_cycle + 3):max_cycle,])$coefficients[[2]]
  # m_diff = rn_m - sigmoid_m
  # err = mean(abs((rn_df$normalized[inflection_cycle:max_cycle] - 
  #                sigmoid_df$y[inflection_cycle:max_cycle])/rn_df$normalized[inflection_cycle:max_cycle]))
  # return(sign(m_diff)*err)
  
  max_cycle = max(rn_df$cycle_no)
  inflection_cycle = find_switching_cycle_inv(rn_df$normalized,10, max_cycle)
  if(inflection_cycle == max_cycle){
    inflection_cycle = max_cycle - 5
  }
  sigmoid_m = lm(y ~ x, sigmoid_df[c(inflection_cycle,max_cycle),])$coefficients[[2]]
  rn_m = lm(normalized ~ cycle_no, rn_df[c(inflection_cycle,max_cycle),])$coefficients[[2]]
  # in param_optimize function, if err is negative we will step up
  # for param b, if we step up, it means we decrease the slope
  err = (rn_m - sigmoid_m)/rn_m
  # err = rn_m - sigmoid_m
  return(err)
}

err_funcs = c(a = calc_a_err, 
              b = calc_b_err, 
              r = calc_r_err)
#k = calc_k_err)
param_ls = c('a','b','r')#,'k')
# combine_err_func = function(rn_df, sigmoid_df, err_func_list){
#   # err_func_list is a list of funcs of corresponding params
#   # err_func_list$'a' has param_args and param_weights
#   
#   param_weights = c()
#   errs = c()
#   for(p in param_ls){
#     errs[p] = abs(err_funcs[[p]](rn_df, sigmoid_df, err_func_list[[p]]$param_args))
#     param_weights[p] = err_func_list[[p]]$weight
#   }
#   err = weighted.mean(errs, param_weights)
#   return(err)
# }

combine_err_func = function(rn_df, sigmoid_df){
  max_cycle = max(rn_df$cycle_no)
  baseline_end = unique(rn_df$baseline_end)
  baseline_start = unique(rn_df$baseline_start)
  
  inflection_cycle = find_inflection(rn_df)
  diff = find_shift(rn_df, sigmoid_df,inflection_cycle)
  sigmoid_df$y = sigmoid_df$y - diff
  # in case normalized value = 0, then error will be infty
  if(sum(rn_df$normalized == 0) > 0){
    rn_df$normalized = rn_df$normalized + 0.0000001
  }
  # only care about error since cycle_no 30 forward
  err = weighted.mean(abs((rn_df$normalized[30:max_cycle] - sigmoid_df$y[30:max_cycle])/
                   rn_df$normalized[30:max_cycle]),
                   30:max_cycle)
  
  # if(baseline_end > 31){
  #   err = mean(abs((rn_df$normalized[baseline_start:baseline_end] - sigmoid_df$y[baseline_start:baseline_end])/
  #                    rn_df$normalized[baseline_start:baseline_end]))
  # }else{
  #   err = mean(abs((rn_df$normalized[baseline_end:max_cycle] - sigmoid_df$y[baseline_end:max_cycle])/
  #                    rn_df$normalized[baseline_end:max_cycle]))
  # }
  
  return(err)
}

calc_mse = function(truth, prediction){
  tryCatch({
    mse = mean((truth - prediction)^2)
  },warning = function(w){
    cat(w)
  })
  
  return(mse)
}
