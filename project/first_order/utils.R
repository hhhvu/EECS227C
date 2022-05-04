denom = function(x,b,r,t,k){
  return(r + b*abs(x - t)**k)
}

loss_fn = function(x,y,a,b,r,k,t){
  return(sum(y - a*(x-t)/denom(x,b,r,t,k)**(1/k))**2)
}

grad_a = function(x,y,a,b,r,k,t){
  return(sum(-2*(x-t)/denom(x,b,r,t,k)**(1/k)*(y - a*(x-t)/denom(x,b,r,t,k)**(1/k))))
}

grad_r = function(x,y,a,b,r,k,t){
  first_chunk = 2*y*a*(x-t)/k*denom(x,b,r,t,k)*(-1/k-1)
  second_chunk = -2*a*(x-t)**2/k*denom(x,b,r,k,t)**(-2/k - 1)
  return(sum(first_chunk + second_chunk))
}

grad_b = function(x,y,a,b,r,k,t){
  first_chunk = 2*a*(x-t)*abs(x-t)**k*denom(a,b,r,t,k)**(-1/k -1)
  second_chunk = y - a*(x-t)*denom(a,b,r,t,k)**(-1/k)
  return(sum(first_chunk*second_chunk)/k)
}

grad_k = function(x,y,a,b,r,k,t){
  first_chunk = a*(x - t)*(b*abs(t-x)**k + r)**(-1/k)
  # cat('first', first_chunk)
  second_chunk = log(b*abs(t-x)**k+r)
  # cat('second', second_chunk)
  third_chunk = (b*k*abs(t-x)**k*log(abs(t-x)))/(b*abs(t-x)**k+r)
  # cat('third', third_chunk)
  return (sum(first_chunk*(second_chunk - third_chunk)/k**2))
  # first_chunk = y - a*(x-t)/denom(x,b,r,t,k)**(1/k)
  # second_chunk = -a*(x-t)*denom(x,b,r,t,k)**(-1/k)*log(denom(x,b,r,t,k))/k**2
  # third_chunk = (denom(x,b,r,t,k) - r)*log(abs(x - t))/(k*denom(x,b,r,t,k))
  # return(sum(2*first_chunk*(second_chunk - third_chunk)))
}

grad = function(x,y,a,b,r,k,t){
  return(c(grad_a(x,y,a,b,r,k,t),
           grad_b(x,y,a,b,r,k,t),
           grad_r(x,y,a,b,r,k,t),
           grad_k(x,y,a,b,r,k,t)))
}

update_x = function(x,y,a,b,r,k,t,gamma,delta){
  x_next = c(a,b,r,k) - gamma*grad(x,y,a,b,r,k,t) + alpha*delta
  return(x_next)
}

calc_rmse = function(truth, prediction){
  # tryCatch({
  #   rmse = sqrt(mean((truth - prediction)^2))
  # },warning = function(w){
  #   cat(unlist(w))
  # })
  rmse = sqrt(mean((truth - prediction)^2))
  return(rmse)
}

find_inflection = function(rn_df){
  max_cycle = max(rn_df$cycle_no)
  x = rn_df$normalized
  diff = c(0,x[-1] - x[-max_cycle])
  inflection_cycle = max(which.max(diff), unique(rn_df$baseline_end))
  # return(inflection_cycle)
  return(which.max(diff))
  # return(round((switch_cycle + max(rn_df$cycle_no))/2))
}


general_sigmoid_frst = function(x,param_set){
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
  # s = param_set[5]
  t = param_set[5]
  return(a*(x-t)/(r + b*abs(x-t)^k)^(1/k))
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

fit_qpcr = function(idx,data){
  # print(idx)
  data = data[data$sample_target_idx == idx,]
  tryCatch({
    fit = try(pcrfit(data,'cycle_no','normalized',model = l5))
  },error=function(e){
    cat(idx,conditionMessage(e),'\n')
  })
  prediction = predict(fit,data.frame(Cycles = data$cycle_no))
  rmse = calc_rmse(data$normalized,prediction$Prediction)
  # fit = pcrfit(data,'cycle_no','rn', model = l5)
  # if(class(fit) == 'try-error'){
  #   data$cycle_no = data$cycle_no + 0.001
  #   fit = pcrfit(data,'cycle_no','normalized', model = l4)
  # }
  # prediction = predict(fit,data.frame(Cycles = data$cycle_no))
  # rmse = calc_rmse(data$normalized,prediction$Prediction)
  return(rmse)
}




