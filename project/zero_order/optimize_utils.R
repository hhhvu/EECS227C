normalizing = function(vec){
  std_vec = (vec - min(vec))/(max(vec) - min(vec))
  return(std_vec)
}

extract_sigma = function(x,y){
  df = data.frame(x = x, y = normalizing(y))
  df = df[5:nrow(df),]
  sig = sigma(lm(y ~ x, df))
  return(sig)
}

get_next_param = function(current_val, direction, step, target_param){
  # direction: +1 is forward and -1 is backward
  if(target_param == 'a'){
    next_val = current_val + direction*step
  }else if(target_param == 'b'){
    next_val = current_val + direction*step
    while(next_val < 0){
      step = 0.9*step
      next_val = current_val + direction*step
    }
  }else if(target_param == 'r'){
    next_val = current_val + direction*step
    while(next_val <= 0){
      step = 0.9*step
      next_val = current_val + direction*step
    }
  }else if(target_param == 'k'){
    next_val = current_val + direction*step
    while(next_val <= 0){
      step = 0.9*step
      next_val = current_val + direction*step
    }
  }else{
    next_val = current_val + direction*step
  }
  return(next_val)
}


get_next_param = function(current_val, direction, step, target_param){
  # direction: +1 is forward and -1 is backward
  if(target_param == 'a'){
    next_val = current_val + direction*step
  }else if(target_param == 'b'){
    next_val = current_val + direction*step
    i = 1
    while(next_val < 0){
      next_val = current_val + direction*(step/(2^i))
      i = i +1
    }
  }else if(target_param == 'r'){
    next_val = current_val + direction*step
    i = 1
    while(next_val <= 0){
      next_val = current_val + direction*(step/(2^i))
      i = i +1
    }
  }else if(target_param == 'k'){
    next_val = current_val + direction*step
    i = 1
    while(next_val <= 0){
      next_val = current_val + direction*(step/(2^i))
      i = i +1
    }
  }else{
    next_val = current_val + direction*step
  }
  return(next_val)
}

check_both_dirs = function(param, param_set, step, current_err,
                           err_func, rn_df, param_args){
  # param (character): name of tuning parameter 
  # param_set (array): array of parameters for general_sigmoid function
  # step (float): step size of changing param
  # current_err (float)
  # err_func (function); corresponding function to param
  # rn_df (data.frame)
  # param_args (array): arrays of arguments for err_func
  
  up_param = get_next_param(current_val = param_set[[param]],
                            direction = 1, step = step,
                            target_param = param)
  down_param = get_next_param(current_val = param_set[[param]],
                              direction = -1, step = step,
                              target_param = param)
  up_down_param = c(up_param, down_param)
  param_set[[param]] = up_param
  up_df = data.frame(x = rn_df$cycle_no,
                     y = general_sigmoid(rn_df$cycle_no,param_set))
  param_set[[param]] = down_param
  down_df = data.frame(x = rn_df$cycle_no,
                       y = general_sigmoid(rn_df$cycle_no,param_set))
  up_down_err_gap = current_err - c(err_func(rn_df, up_df, param_args),
                                    err_func(rn_df, down_df, param_args))
  same_sign = sign(current_err) == sign(up_down_err_gap)
  if(sum(same_sign) %in% c(0,2)){
    param_set[[param]] = up_down_param[which.min(abs(up_down_err_gap))]
  }else{param_set[[param]] = up_down_param[same_sign]}
  return(param_set)
}

check_both_dirs_1 = function(param, param_set, step, current_err,
                           err_func, rn_df){
  # param (character): name of tuning parameter 
  # param_set (array): array of parameters for general_sigmoid function
  # step (float): step size of changing param
  # current_err (float)
  # err_func (function); corresponding function to param
  # rn_df (data.frame)
  # param_args (array): arrays of arguments for err_func
  max_cycle = max(rn_df$cycle_no)
  # flag = TRUE
  # while(flag){
  #   up_param = get_next_param(current_val = param_set[[param]],
  #                             direction = 1, step = step,
  #                             target_param = param)
  #   down_param = get_next_param(current_val = param_set[[param]],
  #                               direction = -1, step = step,
  #                               target_param = param)
  #   
  #   up_down_param = c(up_param, down_param)
  #   param_set[[param]] = up_param
  #   up_df = data.frame(x = rn_df$cycle_no,
  #                      y = general_sigmoid(rn_df$cycle_no,param_set))
  #   # impute previous value if the last value is NAN
  #   if(is.na(up_df$y[max_cycle])){
  #     up_df$y[max_cycle] = up_df$y[max_cycle-1]
  #   }
  #   param_set[[param]] = down_param
  #   down_df = data.frame(x = rn_df$cycle_no,
  #                        y = general_sigmoid(rn_df$cycle_no,param_set))
  #   if(is.na(down_df$y[max_cycle])){
  #     down_df$y[max_cycle] = down_df$y[max_cycle-1]
  #   }
  #   if((sum(is.na(up_df$y)) == 0) & (sum(is.na(down_df$y) == 0))){
  #     flag = FALSE
  #   }else{step = step*0.9}
  # }
  up_param = get_next_param(current_val = param_set[[param]],
                            direction = 1, step = step,
                            target_param = param)
  down_param = get_next_param(current_val = param_set[[param]],
                              direction = -1, step = step,
                              target_param = param)
  
  up_down_param = c(up_param, down_param)
  param_set[[param]] = up_param
  up_df = data.frame(x = rn_df$cycle_no,
                     y = general_sigmoid(rn_df$cycle_no,param_set))
  # impute previous value if the last value is NAN
  if(is.na(up_df$y[max_cycle])){
    up_df$y[max_cycle] = up_df$y[max_cycle-1]
  }
  param_set[[param]] = down_param
  down_df = data.frame(x = rn_df$cycle_no,
                       y = general_sigmoid(rn_df$cycle_no,param_set))
  if(is.na(down_df$y[max_cycle])){
    down_df$y[max_cycle] = down_df$y[max_cycle-1]
  }
  # error sign is always positive
  # if error gap is positive, then err is reduced => take the one with largest
  # error gap
  # if err gap is negative, then err is increased => take the one with smallest
  # error gap
  up_down_err_gap = current_err - c(err_func(rn_df, up_df),
                                    err_func(rn_df, down_df))
  pos_sign = sign(up_down_err_gap) == 1
  # if both errors are positive (sum(pos_sign) == 2), take largest err gap
  # if both errors are negative (sum(pos_sign) == 0), take the smallest abs err gap
  # if one error is positive, then take that error
  if(sum(pos_sign) == 0){
    param_set[[param]] = up_down_param[which.min(abs(up_down_err_gap))]
  }else if(sum(pos_sign) == 2){
    param_set[[param]] = up_down_param[which.max(up_down_err_gap)]
  } else{param_set[[param]] = up_down_param[pos_sign]}
  return(param_set)
}


general_sigmoid_1 = function(x,a,b,r,k,t){
  # a is to change the distance between two flat lines of the sigmoid curve (
  #  if a < 0 , then the curve is flipped (left high, right low); 
  #  bigger abs(a) gives bigger distance; changing a just changes the distance 
  #  but still keeps the shape of the curve)
  # b > 0 controls the slopes (b = 0 is a straight line with slope = 4.5/30,
  #  the bigger b value is, the steeper the slope, but the distance between
  #  two horizontal lines decreases as b gets bigger)
  # r > 0 otherwise we will have NA; r also controls the slopes; the larger the
  #  values of r, the flatter the slope
  # k >= 0 otherwise we will have NA (changing k also changes the distance between
  #  the two lines as well as the steepness but more on the distance; increasing
  #  k by one will increase the distance much more than increasing a by 1)
  # s is to shift the curve up and down on y axis (s > 0 is to shift up,
  #  s < 0 is to shift down)
  # t is point of inflection in sigmoid curve
  #  to shift the curve left or right on x axis (t < 0 is to shift right, 
  #  t > 0 is to shift left)
  
  return(a*(x-t)/(r + b*abs(x-t)^k)^(1/k))
}


general_sigmoid = function(x,param_set){
  # a is to change the distance between two flat lines of the sigmoid curve (
  #  if a < 0 , then the curve is flipped (left high, right low); 
  #  bigger abs(a) gives bigger distance; changing a just changes the distance 
  #  but still keeps the shape of the curve)
  # b > 0 controls the slopes (b = 0 is a straight line with slope = 4.5/30,
  #  the bigger b value is, the steeper the slope, but the distance between
  #  two horizontal lines decreases as b gets bigger)
  # r > 0 otherwise we will have NA; r also controls the slopes; the larger the
  #  values of r, the flatter the slope
  #     increase b from 1 to 10 is equivalent to decreasing r from 1 to 0.1
  #     the beginning and the end cycles stay the same in both cases; however,
  #     decrease r from 1 to 0.1 gives a bigger gap of linearly increasing section
  #       biggest gap in y with a complete extension is 3 for r = 0.1 and 1.3 for 
  #       b = 10
  # k >= 0 otherwise we will have NA (changing k also changes the distance between
  #  the two lines as well as the steepness but more on the distance; increasing
  #  k by one will increase the distance much more than increasing a by 1)
  # s is to shift the curve up and down on y axis (s > 0 is to shift up,
  #  s < 0 is to shift down)
  # t is point of inflection in sigmoid curve
  #  to shift the curve left or right on x axis (t < 0 is to shift right, 
  #  t > 0 is to shift left)
  a = param_set[1]
  b = param_set[2]
  r = param_set[3]
  k = param_set[4]
  s = param_set[5]
  t = param_set[6]
  return(s + a*(x-t)/(r + b*abs(x-t)^k)^(1/k))
}


find_shift = function(rn_df, sigmoid_df, inflection_cycle){
  baseline_start = unique(rn_df$baseline_start)
  baseline_end = unique(rn_df$baseline_end)
  max_cycle = max(rn_df$cycle_no)
  period_idx = c(inflection_cycle - baseline_start, max_cycle - inflection_cycle)
  # if baseline_end > 31 then we will consider shift from baseline_start to
  # inflection_cycle; otherwise, baseline_end to max_cycle
  if(baseline_end > 31){
    rn_mean = mean(rn_df$normalized[baseline_start:inflection_cycle])
    sigmoid_mean = mean(sigmoid_df$y[baseline_start:inflection_cycle])
    diff = sigmoid_mean - rn_mean
  }else{
    rn_mean = mean(rn_df$normalized[baseline_end:max_cycle])
    sigmoid_mean = mean(sigmoid_df$y[baseline_end:max_cycle])
    diff = sigmoid_mean - rn_mean
  }
  # if(which.max(period_idx) == 1){
  #   rn_mean = mean(rn_df$normalized[baseline_start:inflection_cycle])
  #   sigmoid_mean = mean(sigmoid_df$y[baseline_start:inflection_cycle])
  #   diff = sigmoid_mean - rn_mean
  # }else{
  #   rn_mean = mean(rn_df$normalized[inflection_cycle:max_cycle])
  #   sigmoid_mean = mean(sigmoid_df$y[inflection_cycle:max_cycle])
  #   diff = sigmoid_mean - rn_mean
  # }
  return(diff)
}



