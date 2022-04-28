source('/Users/huongvu/Desktop/PCR_Desktop/Scripts/fitting_rn_curves/err_calc.R')
source('/Users/huongvu/Desktop/PCR_Desktop/Scripts/fitting_rn_curves/optimize_utils.R')

err_funcs = c(combine_err_func = combine_err_func,
              calc_mse = calc_mse)
optimize_param = function(target_err_thres = 0.05, 
                          param_thres = 0.01, 
                          combined_err_thres = 0.05,
                          rn_df, param_set, param_step, param_args,
                          target_param){
# 
  # rn_df = train %>% filter(sample_target_idx == 1049)
  # target_err_thres = 0.001
  # param_thres = 0.0001
  # combined_err_thres = 0.001
  # param_args = list(param_args = c(35,40))
  # find_inflection(rn_df)
  # param_set = c(1,0,1,1,0,40)
  # param_step = 1
  # target_param = 'k'

  
  flag = TRUE
  max_cycle = max(rn_df$cycle_no)
  names(param_set) = c('a','b','r','k','s','t')
  first_param_val = param_set[[target_param]]
  # param_args = err_func_inputs[[target_param]]$param_args
  step = param_step
  param_series = c()
  combined_err_series = c()
  target_err_series = c()
  all_combined_errs = c()
  all_params = c()
 
  # step = 2
  i = 0
  while(flag){
    # print(i)
    sigmoid_df = data.frame(x = rn_df$cycle_no,
                            y = general_sigmoid(rn_df$cycle_no,param_set))
    if(is.na(sigmoid_df$y[max_cycle])){
      sigmoid_df$y[max_cycle] = sigmoid_df$y[max_cycle-1] 
    }
    combined_err = combine_err_func(rn_df, sigmoid_df)
    all_combined_errs = c(combined_err, all_combined_errs)
    all_params = c(param_set[[target_param]], all_params)
    if(i == 0){
      combined_err_series = c(combined_err, combined_err_series)
      param_series = c(param_set[[target_param]], param_series)
      catch = TRUE
      while(catch){
        try_check = try(check_both_dirs_1(target_param, param_set, step, combined_err,
                                                      combine_err_func, rn_df))
        if(class(try_check) == 'try-error'){
          step = step*0.9
        }else{
          param_set = check_both_dirs_1(target_param, param_set, step, combined_err,
                                        combine_err_func, rn_df)
          catch = FALSE
        }
      }
    }else if(i > 100){
      warning('When optimizing param', target_param, ', more than 100 loops\n')
      best_err_idx = which.min(combined_err_series)
      combined_err = combined_err_series[best_err_idx]
      param_set[[target_param]] = param_series[best_err_idx]
      flag = FALSE
    }else{
      if(combined_err <= combined_err_thres){
        flag = FALSE
      }else{
        # if the current iteration gives worse result than the previous iteration,
        # we return to the previous iteration including param_set and err
        if(combined_err >= combined_err_series[1]){
          step = step*0.9
          param_set[[target_param]] = param_series[1]
          combined_err = combined_err_series[1]
          catch = TRUE
          while(catch){
            try_check = try(check_both_dirs_1(target_param, param_set, step, combined_err,
                                              combine_err_func, rn_df))
            if(class(try_check) == 'try-error'){
              step = step*0.9
            }else{
              param_set = check_both_dirs_1(target_param, param_set, step, combined_err,
                                            combine_err_func, rn_df)
              catch = FALSE
            }
          }
        }else{
          param_series = c(param_set[[target_param]], param_series)
          combined_err_series = c(combined_err, combined_err_series)
          catch = TRUE
          while(catch){
            try_check = try(check_both_dirs_1(target_param, param_set, step, combined_err,
                                              combine_err_func, rn_df))
            if(class(try_check) == 'try-error'){
              step = step*0.9
            }else{
              param_set = check_both_dirs_1(target_param, param_set, step, combined_err,
                                            combine_err_func, rn_df)
              catch = FALSE
            }
          }
        }
      }
    }
    i = i + 1
  }
  return(c(param_set, combined_err = combined_err))
}


gridSearch = function(rn_df, max_cycle, combined_err_thres, target_err_thres,
                      param_thres, param_step,
                      param_set, param_order){
  # index = 2437
  # combined_err_thres = 0.001
  # target_err_thres = c(a = 0.001, b = 0.001, r = 0.001, k = 0.001)
  # param_thres = c(a = 0.01, b = 0.01, r = 0.01, k = 0.01)
  # param_step = c(a = 0.1, b = 1, r = 1, k = 1)
  # # initialize params
  # param_set = c(1,1,1,1,0,1)

  # print(combined_err_thres)
  # cat('Working on', index, '\n')
  # rn_df = train %>% filter(sample_target_idx == index)
  # max_cycle = max(rn_df$cycle_no)
  names(param_set) = c('a','b','r','k','s','t')
  # create historical param saving variables 
  target_err_series = c()
  combined_err_series = c()
  
  # # first determine inflection_cycle
  t = param_set[['t']]
  # first optimize b
  for(p in param_order){
    cat('Optimizing param', p, ':\n')
    if(p == 'a'){
      param_args = c(max(1,t-5),min(max_cycle,t+5))
    }else if(p == 'r'){
      param_args = c(inflection_cycle = t,
                     slope_range = 1)
    }else{
      param_args = c(inflection_cycle = t)
    }
    optimized_results = optimize_param(combined_err_thres = combined_err_thres, 
                                       param_thres = param_thres[[p]],
                                       target_err_thres = target_err_thres[[p]],
                                       rn_df = rn_df, 
                                       param_set = param_set, 
                                       param_step = param_step[[p]], 
                                       param_args = param_args, 
                                       target_param = p)
    param_set = optimized_results[1:6]
    
    # if error satisfy err_thres, we can stop
    if(optimized_results[['combined_err']] <= combined_err_thres){
      cat('Satisfied error threshold early')
      break
    }
    
  }
  
  sigmoid_df = data.frame(x = rn_df$cycle_no,
                          y = general_sigmoid(x = rn_df$cycle_no,
                                              param_set = param_set))
  shift = find_shift(rn_df = rn_df, 
                     sigmoid_df = sigmoid_df,
                     inflection_cycle = t)
  param_set[['s']] = shift
  sigmoid_df$y = sigmoid_df$y - shift
  err = combine_err_func(rn_df, sigmoid_df)
  
  return(c(param_set, err))
}

fitting_sigmoid = function(index, combined_err_thres, target_err_thres,
                           param_thres, param_step,
                           param_set, param_order){
  cat('Working on sample_target_idx = ', index,'\n')
  rn_df = data %>% filter(sample_target_idx == index)
  if(unique(rn_df$baseline_end) > 31){
    # we want to fit a straight line
    # we will fix b = 0, k = 1 and t = 40
    max_cycle = max(rn_df$cycle_no)
    param_set = c(1,0,1,1,0,max_cycle)
    param_order = c('a','r','a','r','a','r','a','r')
    optimized_result = gridSearch(rn_df, max_cycle, combined_err_thres, 
                                  target_err_thres, param_thres, param_step,
                                  param_set, param_order)
  }else{
    max_cycle = max(rn_df$cycle_no)
    t = find_inflection(rn_df)
    param_set[6] = t
    optimized_result = gridSearch(rn_df, max_cycle, combined_err_thres, 
                                  target_err_thres, param_thres, param_step,
                                  param_set, param_order)
  }
  return(c(index,optimized_result))
}
