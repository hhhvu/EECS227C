source('/Users/huongvu/Desktop/EECS227C/project/first_order/utils.R')
source('/Users/huongvu/Desktop/EECS227C/project/zero_order/optimize_utils.R')

first_order_fit = function(data,idx,alpha,gamma,err_change_thres,n_rdm,n_rep){
  # first we randomly generate initial points within conditions for each parameter
  # then we optimize each param while holding others constant
  # repeat optimizing params 25 times per param and
  # repeat the whole process with 10 different random initial points
  
  data = data[data$sample_target_idx == idx,]
  x = data$cycle_no
  y = data$normalized
  bestparam_ls = c()
  bestloss_ls = c()
  t = find_inflection(data) + 0.0001
  for(pt in 1:n_rdm){
    a = runif(1, min = -10, max = 10)
    flag = TRUE
    while(flag){
      b = runif(1, min = 0, max = 10)
      if(b > 0){
        flag = FALSE
      }
    }
    flag = TRUE
    while(flag){
      r = runif(1, min = 0, max = 10)
      if(r > 0){
        flag = FALSE
      }
    }
    flag = TRUE
    while(flag){
      k = runif(1, min = 0, max = 3)
      if(k > 0){
        flag = FALSE
      }
    }
    
    delta_a = 0
    delta_b = 0
    delta_r = 0
    delta_k = 0
    a_ls = c()
    b_ls = c()
    r_ls = c()
    k_ls = c()
    best_loss = loss_fn(x,y,a,b,r,k,t)
    rmse = best_loss
    best_params = c(a,b,r,k)
    rmse_ls = c(best_loss)
    patient_cnt = 0
    
    for(p in rep(c('b','r','a','k'),n_rep)){
      if(p == 'b'){
        change = grad_b(x,y,a,b,r,k,t)
        b_next = b - gamma*change + alpha*delta_b
        b_ls = c(b_ls,b)
        if(b_next > 0){
          b_ls = c(b_ls,b)
          rmse = loss_fn(x,y,a,b_next,r,k,t)
          
          if(rmse <= best_loss){
            best_loss = rmse
            best_params = c(a,b,r,k)
            delta_b = b_next - b
            b = b_next
          }
          rmse_ls = c(rmse_ls, rmse)
          if(abs(rmse - rmse_ls[length(rmse_ls) - 1]) <= err_change_thres){
            patient_cnt =+ 1
          }
          if (patient_cnt > 5){
            break
          }
        }else{
          rmse_ls = c(rmse_ls, rmse)
        }
        
      }else if(p == 'r'){
        change = grad_r(x,y,a,b,r,k,t)
        r_next = r - gamma*change + alpha*delta_r
        # print(r_next)
        if(is.infinite(r_next)){
          return(c(x,y,a,b,r,k,t,delta_r))
          break
        }
        r_ls = c(r_ls,r)
        if(r_next > 0){
          rmse = loss_fn(x,y,a,b,r_next,k,t)
          if(rmse <= best_loss){
            best_loss = rmse
            best_params = c(a,b,r,k)
            delta_r = r_next - r
            r = r_next
          }
          rmse_ls = c(rmse_ls, rmse)
          if(abs(rmse - rmse_ls[length(rmse_ls) - 1]) <= err_change_thres){
            patient_cnt =+ 1
          }
          if (patient_cnt > 5){
            break
          }
        }else{
          rmse_ls = c(rmse_ls, rmse)
        }
      }else if(p == 'a'){
        change = grad_a(x,y,a,b,r,k,t)
        a_next = a - gamma*change + alpha*delta_a
        a_ls = c(a_ls,a)
        rmse = loss_fn(x,y,a_next,b,r,k,t)
        if(rmse <= best_loss){
          best_loss = rmse
          best_params = c(a,b,r,k)
          delta_a = a_next - a
          a = a_next
        }
        rmse_ls = c(rmse_ls, rmse)
        if(abs(rmse - rmse_ls[length(rmse_ls) - 1]) <= err_change_thres){
          patient_cnt =+ 1
        }
        if (patient_cnt > 5){
          break
        }
      }else{
        # print(k)
        change = grad_k(x,y,a,b,r,k,t)
        k_next = k - gamma*change + alpha*delta_k
        if(is.nan(k_next)){
          print(c(x,y,a,b,r,k,t,delta_k))
        }
        k_ls = c(k_ls,k)
        if(k_next >= 0){
          rmse = loss_fn(x,y,a,b,r,k_next,t)
          if(rmse <= best_loss){
            best_params = c(a,b,r,k)
            best_loss = rmse
            delta_k = k_next - k
            k = k_next
          }
          rmse_ls = c(rmse_ls, rmse)
          if(abs(rmse - rmse_ls[length(rmse_ls) - 1]) <= err_change_thres){
            patient_cnt =+ 1
          }
          if (patient_cnt > 5){
            break
          }
        }else{
          rmse_ls = c(rmse_ls,rmse)
        }
        
      }
    }
    
    bestloss_ls = c(bestloss_ls, best_loss)
    bestparam_ls = rbind(bestparam_ls, best_params)
    
  }
  
  best_param = bestparam_ls[which.min(bestloss_ls),]
  # prediction = general_sigmoid(x,c(best_param,t))
  # rmse = calc_rmse(y,prediction)
  return(c(best_param,t))
}
