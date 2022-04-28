pcr_dir = '/Users/huongvu/Desktop/PCR_desktop'
pcr_data_dir = file.path(pcr_dir, 'Data/New Data')

rn_df = read.csv(file.path(pcr_data_dir, 'full_data.csv'))
gene_idx = unique(rn_df$sample_target_idx)

idx = 1
alpha = 0.05
gamma = 0.05
ITER = 100

y = unique(rn_df[rn_df$sample_target_idx == gene_idx[idx],'rn'])
y = (y - min(y))/(max(y) - min(y))
x = seq(1,length(y) + 1)
