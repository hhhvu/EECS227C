library(qpcR)
library(ggplot2)
library(dplyr)
library(cowplot)
pcr_dir = '/Users/huongvu/Desktop/PCR_desktop'
pcr_data_dir = file.path(pcr_dir, 'Data/New Data')
eecs_dir = '/Users/huongvu/Desktop/EECS227C'
eecs_project_dir = file.path(eecs_dir, 'project')

source(file.path(eecs_project_dir,'first_order','frst_order_method_fit.R'))
source(file.path(eecs_project_dir,'zero_order','param_optimize_1.R'))

full_data = read.csv(file.path(pcr_data_dir, 'full_data.csv'))
full_data = full_data %>%
  group_by(sample_target_idx) %>%
  mutate(normalized = normalizing(rn),
         rsq = extract_sigma(cycle_no, normalized))
idx_45 = unique(full_data[full_data$cycle_no == 45,'sample_target_idx'])
rn_df = full_data %>%
  filter(!(sample_target_idx %in% idx_45)) %>%

gene_idx = unique(rn_df$sample_target_idx)
# remove curves that qPCR can't run using qpcR.R file to detect
singular_vec = c(3,11,38,54,57,66,75,77,78,80,83,107,119,
                 121,125,137,155,168,180,186,189,225,248,266,268,
                 281,302,318,333,336,339,350)
gene_wo_singular = gene_idx[-singular_vec]

### Getting RMSE and running time
##########################################
######## Without singular genes ##########
##########################################
### run zero-order method
running_time = c()
avg_rmse = c()
for(n in c(30,50,75,100,125,150,175,200)){
  n = 1
  ids = gene_wo_singular[1:n]
  start_time = Sys.time()
  out = lapply(ids, 
              FUN = fitting_sigmoid,
              data = rn_df,
              combined_err_thres = 0.001,
              target_err_thres = c(a = 0.001, b = 0.001, r = 0.001,k = 0.001),
              param_thres = c(a = 0.0001, b = 0.0001, r = 0.0001, k = 0.0001),
              param_step = c(a = 0.1, b = 1, r = 1, k = 1),
              param_set = c(1,1,1,1,0,1),
              param_order = c('b','r','a','k',
                              'b','r','a','k',
                              'b','r','a','k',
                              'b','r','a','k'))
  end_time = Sys.time()
  running_time = c(running_time, end_time - start_time)

  result_df = data.frame(do.call(rbind,out))
  avg_rmse = c(avg_rmse,mean(result_df[,5]))
}
# save results
compare_df = data.frame(running_time = running_time,
                        rmse = avg_rmse,
                        method = 'Zero-Order',
                        n = c(30,50,75,100,125,150,175,200))

### run first-order method
running_time = c()
avg_rmse = c()
for(n in c(30,50,75,100,125,150,175,200)){
  ids = gene_wo_singular[1:n]
  start_time = Sys.time()
  out = lapply(ids, 
               FUN = first_order_fit,
               data = rn_df,
               alpha = 0.5,
               gamma = 0.5, 
               err_change_thres = 0.001,
               n_rdm = 10,
               n_rep = 25
               )
  end_time = Sys.time()
  rmse = c()
  for (i in 1:length(out)){
    # print(out[i])
    prediction = general_sigmoid(x,out[[i]])
    rmse = c(rmse,calc_rmse(y,prediction))
  }
  running_time = c(running_time, end_time - start_time)
  avg_rmse = c(avg_rmse, mean(rmse))
}
# save result
compare_df = rbind(compare_df,
                   data.frame(running_time = running_time,
                              rmse = avg_rmse,
                              method = 'First-Order',
                              n = c(30,50,75,100,125,150,175,200)))

### run qPCR result
running_time = c()
avg_rmse = c()
for(n in c(30,50,75,100,125,150,175,200)){
  # cat('n = ',n)
  ids = gene_wo_singular[1:n]
  # print(ids)
  start_time = Sys.time()
  out = lapply(ids, 
               FUN = fit_qpcr,
               data = rn_df)
  end_time = Sys.time()
  print(length(out))
  running_time = c(running_time, end_time - start_time)
  avg_rmse = c(avg_rmse, mean(unlist(out)))
}


compare_df = rbind(compare_df,
                   data.frame(running_time = running_time,
                              rmse = avg_rmse,
                              method = 'qPCR',
                              n = c(30,50,75,100,125,150,175,200)))

#########################################
########## Any gene curves ##############
#########################################

### run zero-order method
running_time = c()
avg_rmse = c()
for(n in c(30,50,75,100,125,150,175,200)){
  ids = gene_idx[1:n]
  start_time = Sys.time()
  out = lapply(ids, 
               FUN = fitting_sigmoid,
               data = rn_df,
               combined_err_thres = 0.001,
               target_err_thres = c(a = 0.001, b = 0.001, r = 0.001,k = 0.001),
               param_thres = c(a = 0.0001, b = 0.0001, r = 0.0001, k = 0.0001),
               param_step = c(a = 0.1, b = 1, r = 1, k = 1),
               param_set = c(1,1,1,1,0,1),
               param_order = c('b','r','a','k',
                               'b','r','a','k',
                               'b','r','a','k',
                               'b','r','a','k'))
  end_time = Sys.time()
  running_time = c(running_time, end_time - start_time)
  
  result_df = data.frame(do.call(rbind,out))
  avg_rmse = c(avg_rmse,mean(result_df[,5]))
}
# save results
all_compare_df = data.frame(running_time = running_time,
                        rmse = avg_rmse,
                        method = 'Zero-Order',
                        n = c(30,50,75,100,125,150,175,200))

### run first-order method
running_time = c()
avg_rmse = c()
for(n in c(30,50,75,100,125,150,175,200)){
  ids = gene_idx[1:n]
  start_time = Sys.time()
  out = lapply(ids, 
               FUN = first_order_fit,
               data = rn_df,
               alpha = 0.5,
               gamma = 0.5, 
               err_change_thres = 0.001,
               n_rdm = 10,
               n_rep = 25
  )
  end_time = Sys.time()
  running_time = c(running_time, end_time - start_time)
  avg_rmse = c(avg_rmse, mean(unlist(out)))
}
# save result
all_compare_df = rbind(all_compare_df,
                   data.frame(running_time = running_time,
                              rmse = avg_rmse,
                              method = 'First-Order',
                              n = c(30,50,75,100,125,150,175,200)))

#### plot results of three methods
a = ggplot(compare_df, aes(as.factor(n),rmse, group = method, color = method)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = 'None') +
  labs(x = 'No. of fitted curves', y ='RMSE')
b = ggplot(compare_df, aes(as.factor(n),running_time, group = method, color = method)) +
  geom_line() +
  theme_bw() +
  labs(x = 'No. of fitted curves', y ='Running Time (s)', color = 'Method')
c = ggplot(all_compare_df, aes(as.factor(n),rmse, group = method, color = method)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = 'None') +
  labs(x = 'No. of fitted curves', y ='RMSE')
d = ggplot(all_compare_df, aes(as.factor(n),running_time, group = method, color = method)) +
  geom_line() +
  theme_bw() +
  labs(x = 'No. of fitted curves', y ='Running Time (s)', color = 'Method')

plot_grid(a,b,c,d,rel_widths = c(0.8,1), nrow = 2,
          labels = c('a','b','c','d'))

### plot excluded curves
singular_5 = gene_idx[singular_vec[1:100]] #43,73, 114, 345
excluded_idx = c(43,73,114,345)
ggplot(full_data %>% filter(sample_target_idx %in% excluded_idx),
       aes(cycle_no, normalized)) +
  geom_line() +
  facet_wrap(sample_target_idx ~ .) +
  theme_bw() +
  theme(legend.position = 'None') +
  labs(x = 'Cycle No.', y = 'Normalized Rn')

#### plot fitted values 
## excluded curves
first_fit = lapply(excluded_idx, 
                   FUN = first_order_fit,
                   data = full_data,
                   alpha = 0.5,
                   gamma = 0.5, 
                   err_change_thres = 0.001,
                   n_rdm = 10,
                   n_rep = 25)
zero_fit = lapply(excluded_idx, 
                  FUN = fitting_sigmoid,
                  data = full_data,
                  combined_err_thres = 0.001,
                  target_err_thres = c(a = 0.001, b = 0.001, r = 0.001,k = 0.001),
                  param_thres = c(a = 0.0001, b = 0.0001, r = 0.0001, k = 0.0001),
                  param_step = c(a = 0.1, b = 1, r = 1, k = 1),
                  param_set = c(1,1,1,1,0,1),
                  param_order = c('b','r','a','k',
                                  'b','r','a','k',
                                  'b','r','a','k',
                                  'b','r','a','k'))

fitted_df = data.frame()
for (i in 1:length(first_fit)){
  data = full_data %>% filter(sample_target_idx == excluded_idx[i]) 
  prediction = general_sigmoid_frst(data$cycle_no, first_fit[[i]])
  fitted_df = rbind(fitted_df,
                       data.frame(cycle_no = data$cycle_no,
                                  pred = prediction,
                                  method = 'First-Order',
                                  sample_target_idx = excluded_idx[i],
                                  truth = data$normalized))
  prediction = general_sigmoid(data$cycle_no, zero_fit[[i]][2:7])
  fitted_df = rbind(fitted_df, data.frame(cycle_no = data$cycle_no, 
                                          pred = prediction,
                                          method = 'Zero-Order',
                                          sample_target_idx = excluded_idx[i],
                                          truth = data$normalized))
}

ggplot(fitted_df,aes(cycle_no, pred, color = method)) +
  geom_line() +
  geom_point(data = fitted_df,
            aes(cycle_no, truth),shape = 1, color = 'black') +
  facet_wrap(sample_target_idx ~.) +
  theme_bw() +
  labs(y = 'Normalized Rn', x = 'Cycle No.', color = 'Method')

## normal sigmoid curves
normal_idx = unique(rn_df$sample_target_idx)[1:20]
normal_idx = c(3456, 3458,3459, 3527)
first_fit = lapply(normal_idx, 
                   FUN = first_order_fit,
                   data = full_data,
                   alpha = 0.5,
                   gamma = 0.5, 
                   err_change_thres = 0.001,
                   n_rdm = 10,
                   n_rep = 25)
zero_fit = lapply(normal_idx, 
                  FUN = fitting_sigmoid,
                  data = full_data,
                  combined_err_thres = 0.001,
                  target_err_thres = c(a = 0.001, b = 0.001, r = 0.001,k = 0.001),
                  param_thres = c(a = 0.0001, b = 0.0001, r = 0.0001, k = 0.0001),
                  param_step = c(a = 0.1, b = 1, r = 1, k = 1),
                  param_set = c(1,1,1,1,0,1),
                  param_order = c('b','r','a','k',
                                  'b','r','a','k',
                                  'b','r','a','k',
                                  'b','r','a','k'))

fitted_df = data.frame()
for (i in 1:length(first_fit)){
  data = full_data %>% filter(sample_target_idx == normal_idx[i]) 
  prediction = general_sigmoid_frst(data$cycle_no, first_fit[[i]])
  fitted_df = rbind(fitted_df,
                    data.frame(cycle_no = data$cycle_no,
                               pred = prediction,
                               method = 'First-Order',
                               sample_target_idx = normal_idx[i],
                               truth = data$normalized))
  prediction = general_sigmoid(data$cycle_no, zero_fit[[i]][2:7])
  fitted_df = rbind(fitted_df, data.frame(cycle_no = data$cycle_no, 
                                          pred = prediction,
                                          method = 'Zero-Order',
                                          sample_target_idx = normal_idx[i],
                                          truth = data$normalized))
  pcrfit = pcrfit(data,'cycle_no', 'normalized', l5)
  prediction = predict(pcrfit, data.frame(Cycles = data$cycle_no))
  fitted_df = rbind(fitted_df, data.frame(cycle_no = data$cycle_no,
                                          pred = prediction$Prediction,
                                          method = 'qPCR',
                                          sample_target_idx = normal_idx[i],
                                          truth = data$normalized))
}

ggplot(fitted_df,aes(cycle_no, pred, color = method)) +
  geom_line() +
  geom_point(data = fitted_df,
            aes(cycle_no, truth), shape = 1, color = 'black') +
  facet_wrap(sample_target_idx ~., scales = 'free') +
  theme_bw() +
  labs(y = 'Normalized Rn', x = 'Cycle No.', color = 'Method')

