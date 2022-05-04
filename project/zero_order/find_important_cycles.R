
# # find switching cycle
# find_switching_cycle = function(x,baseline_end,max_cycle){
#   # switch_cycle = baseline_end
#   # baseline_mse = mean((x[1:baseline_end] - mean(x[1:baseline_end]))^2)
#   # for (i in (baseline_end + 1):max_cycle){
#   #   mse = mean((x[1:i] - mean(x[1:baseline_end]))^2)
#   #   if (mse > 1.2*baseline_mse){
#   #     switch_cycle = i - 1
#   #     break
#   #   }else{switch_cycle = i}
#   # }
#   
#   # calculate rate of change
#   diff = x[-1] - x[-max_cycle]
#   switch_cycle = baseline_end
#   # if rate of change is doubled then set switch_cycle to be previous cycle  
#   for(i in (baseline_end + 1):max_cycle){
#     if(diff[i] > 2*diff[i-1]){
#       switch_cycle = i - 1
#       break
#     }else{switch_cycle = i}
#   }
#   return(switch_cycle)
# }

find_switching_cycle = function(x,baseline_end,max_cycle){
  switch_cycle = baseline_end
  baseline_mse = mean((x[1:baseline_end] - mean(x[1:baseline_end]))^2)
  for (i in (baseline_end + 1):max_cycle){
    mse = mean((x[1:i] - mean(x[1:baseline_end]))^2)
    if (mse > 1.2*baseline_mse){
      switch_cycle = i - 1
      break
    }else{switch_cycle = i}
  }
  
  # # calculate rate of change
  # diff = x[-1] - x[-max_cycle]
  # switch_cycle = baseline_end
  # # if rate of change is doubled then set switch_cycle to be previous cycle  
  # for(i in (baseline_end + 1):max_cycle){
  #   if(diff[i] > 2*diff[i-1]){
  #     switch_cycle = i - 1
  #     break
  #   }else{switch_cycle = i}
  # }
  return(switch_cycle)
}


find_switching_cycle_inv = function(x,cycle_lim,max_cycle){
  # invert the rn curve
  x = sort(1 - x)
  # calculate rate of change
  diff = x[-1] - x[-max_cycle]
  switch_cycle = 2
  # if rate of change is doubled then set switch_cycle to be previous cycle  
  for(i in 2:cycle_lim){
    if(diff[i] > 2*diff[i-1]){
      switch_cycle = i - 1
      break
    }else{switch_cycle = 0}
  }
  
  return(max_cycle - switch_cycle)
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

