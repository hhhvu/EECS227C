library(tidyr)
library(dplyr)
library(pracma)
library(ggplot2)
source('/Users/huongvu/Desktop/PCR_Desktop/Scripts/fitting_rn_curves/param_optimize_1.R')

data = read.csv('/Users/huongvu/Desktop/PCR_Desktop/Data/data_full.csv')
data = data %>% 
  mutate(even_row = mod(row_pos,2),
         even_col = mod(col_pos,2)) %>%
  group_by(sample_target_idx) %>%
  mutate(normalized = normalizing(rn),
         rsq = extract_sigma(cycle_no, normalized))
test = data %>% filter(data_type == 'test')
length(unique(test$sample_target_idx))
fit_out = lapply(unique(test$sample_target_idx), 
              FUN = fitting_sigmoid,
              combined_err_thres = 0.001,
              target_err_thres = c(a = 0.001, b = 0.001, r = 0.001,k = 0.001),
              param_thres = c(a = 0.0001, b = 0.0001, r = 0.0001, k = 0.0001),
              param_step = c(a = 0.1, b = 1, r = 1, k = 1),
              param_set = c(1,1,1,1,0,1),
              param_order = c('b','r','a','k',
                              'b','r','a','k',
                              'b','r','a','k',
                              'b','r','a','k'))

# without k
result_df = data.frame(do.call(rbind,fit_out)) # without k optimize_1
result_df$sample_target_idx = unique(test$sample_target_idx)
write.csv(result_df, 'result_df.csv')

mean(result_df$V7)

#### View some samples
rn_df = data %>% filter(sample_target_idx == 861)
x = 1:nrow(rn_df)
idx = 12
param_set = c(result_df[idx,'a'],result_df[idx,'b'],
              result_df[idx,'r'],result_df[idx,'k'],
              -result_df[idx,'s'],result_df[idx,'t'])
sigmoid_df = data.frame(x = x,
                        y = general_sigmoid(x,param_set))
sigmoid_df$denormalized = sigmoid_df$y*(max(rn_df$rn) - min(rn_df$rn)) + min(rn_df$rn)
ggplot(rn_df, aes(cycle_no, normalized)) + geom_line() +
  geom_line(data = sigmoid_df, aes(x, y), color = 'red')
ggplot(rn_df, aes(cycle_no, rn)) + geom_line() +
  geom_line(data = sigmoid_df, aes(x, denormalized), color = 'red')
### calculate mse on test set
fitting_df = data.frame()
mse = c()
for(i in unique(test$sample_target_idx)){
  rn_df = test %>% filter(sample_target_idx == i)
  x = 1:nrow(rn_df)
  param_set = c(result_df[result_df$sample_target_idx == i,'a'],
                result_df[result_df$sample_target_idx == i,'b'],
                result_df[result_df$sample_target_idx == i,'r'],
                result_df[result_df$sample_target_idx == i,'k'],
                -result_df[result_df$sample_target_idx == i,'s'],
                result_df[result_df$sample_target_idx == i,'t'])
  sigmoid_df = data.frame(sample_target_idx = i,
                          cycle = rn_df$cycle_no,
                          fitted_val = general_sigmoid(x,param_set))
  max_rn = max(rn_df$rn)
  min_rn = min(rn_df$rn)
  sigmoid_df$prediction = sigmoid_df$fitted_val*(max_rn - min_rn) + min_rn
  fitting_df = rbind(fitting_df, sigmoid_df)
  # # denormalizing
  # rn_max = max(rn_df$rn)
  # rn_min = min(rn_df$rn)
  # sigmoid_df$denormalized = sigmoid_df$y*(rn_max - rn_min) + rn_min
  # mse = c(mse, calc_mse(rn_df$rn, sigmoid_df$denormalized))
}
View(fitting_df)
write.csv(fitting_df, 
          '/Users/huongvu/Desktop/PCR_Desktop/Results/fitting_results.csv',
          row.names = FALSE)

result_df$new_mse = mse
mean(result_df$mse)
mean(result_df$new_mse)
names(result_df)[7] = 'fitting_err'
write.csv(result_df, 'fitting_result.csv')

### testing
fitting_sigmoid(index = 2, combined_err_thres = 0.001,
target_err_thres = c(a = 0.001, b = 0.001, r = 0.001,k = 0.001),
param_thres = c(a = 0.0001, b = 0.0001, r = 0.0001, k = 0.0001),
param_step = c(a = 0.1, b = 1, r = 1, k = 1),
param_set = c(1,1,1,1,0,1),
param_order = c('b','r','a','k',
                'b','r','a','k',
                'b','r','a','k',
                'b','r','a','k'))




