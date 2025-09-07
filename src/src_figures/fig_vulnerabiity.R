## Packages
library(R.matlab)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(viridis)
library(raster)
# dev.off()
##
rm(list = ls())

setwd('/Users/phxtao/Google_Drive/R')

#############################################################################
# Data Path
#############################################################################
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

data_dir_output = '/Users/phxtao/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/phxtao/DATAHUB/ENSEMBLE/INPUT_DATA/'

#############################################################################
# function to increase vertical spacing between legend keys
#############################################################################
# # @clauswilke
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

# register new key drawing function, 
# the effect is global & persistent throughout the R session
GeomBar$draw_key = draw_key_polygon3


#################################################################################
# scenario 3 simulation
#################################################################################
depth_num = 4

cpool_num = 7
year_start = 1900
year_end = 2010
year_series = seq(from = year_start, to = year_end, by = 10)
year_num = length(year_series)
scenario_name = c('Static forcing', 'Input change only', 'Climate change only', 'Input + climate change')

soil_sink_base = array(NA, dim = c(36, 4, 10)) 
#-----------------scenario simulation
cross_valid_num = 1
for (cross_valid_num in 1:10) {
  soc_pool_hist_raw = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_scenario_simu_tot.mat', sep = ''))
  soc_pool_hist_raw = soc_pool_hist_raw$global.soc.sink.tot
  
  soil_sink_middle = c()
  
  iscenario = 1
  for (iscenario in 2:length(scenario_name)) {
    cpool_middle_accu = apply(soc_pool_hist_raw[ , 17:28, iscenario], c(1), sum, na.rm = TRUE)
    cpool_middle_accu = cpool_middle_accu - cpool_middle_accu[1]
    cpool_middle_rate = c(0, diff(cpool_middle_accu))/10
    soil_sink_middle = rbind(soil_sink_middle, 
                             cbind(year_series, cpool_middle_accu, cpool_middle_rate, iscenario))
  }
  
  soil_sink_base[ , , cross_valid_num] = soil_sink_middle
}

soil_sink_summary_base = cbind(soil_sink_base[ , 1, 1],
                               apply(soil_sink_base[ , 2, ], 1, mean, na.rm = TRUE),
                               apply(soil_sink_base[ , 2, ], 1, sd, na.rm = TRUE),
                               apply(soil_sink_base[ , 3, ], 1, mean, na.rm = TRUE),
                               apply(soil_sink_base[ , 3, ], 1, sd, na.rm = TRUE),
                               soil_sink_base[ , 4, 1]                         
)

colnames(soil_sink_summary_base) = c('year', 'sink', 'sink_std', 'rate', 'rate_std', 'scenario')
soil_sink_summary_base = data.frame(soil_sink_summary_base)

sum(soil_sink_summary_base[which(soil_sink_summary_base[ , 'scenario'] == 4 & soil_sink_summary_base[ , 'year'] == 2010), 'sink'])
sum(soil_sink_summary_base[which(soil_sink_summary_base[ , 'scenario'] == 4 & soil_sink_summary_base[ , 'year'] == 2010), 'rate'])

soil_sink_summary_base$scenario[soil_sink_summary_base$scenario == 4] = 0
soil_sink_summary_base$reallocation = 0
#################################################################################
#  response curve simulation
#################################################################################
vert_component = 1
x_axis_label = c('Reallocated aboveground input (%)', 'Enhanced vertical transport (%)')
cpool_num = 7
year_start = 1900
year_end = 2010
year_series = seq(from = year_start, to = year_end, by = 10)
year_num = length(year_series)
control_gradient = seq(from = 0, to = 1, by = 0.1)

response_curve_scenario_raw = array(NA, dim = c(12, 3, 4, 10))

iscenario = 3
for (iscenario in 2:4) {
  cross_valid_num = 1
  for (cross_valid_num in 1:10){
    if (iscenario == 4) {
      soc_hist_raw = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_response_curve_soc_tot.mat', sep = ''))
    } else {
      soc_hist_raw = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_response_curve_soc_tot_scenario_', as.character(iscenario-1), '.mat', sep = ''))
    }
    
    soc_hist_raw = soc_hist_raw$global.soc.sink.tot
    soc_change = apply(soc_hist_raw[ , seq(from = 17, to = 28, by = 1), , vert_component], c(1, 3), sum, na.rm = TRUE)
    
    soc_sink_half_reallocate = soc_change[ , 6] - soc_change[1, 6]
    soc_sink_all_reallocate = soc_change[ , 11] - soc_change[1, 11]
    
    sink_rate_half_reallocate =  c(0, diff(soc_sink_half_reallocate))/10
    sink_rate_all_reallocate =  c(0, diff(soc_sink_all_reallocate))/10
    
    response_curve_scenario_raw[ , iscenario-1, 1, cross_valid_num] = sink_rate_half_reallocate
    response_curve_scenario_raw[ , iscenario-1, 2, cross_valid_num] = sink_rate_all_reallocate
    response_curve_scenario_raw[ , iscenario-1, 3, cross_valid_num] = soc_sink_half_reallocate
    response_curve_scenario_raw[ , iscenario-1, 4, cross_valid_num] = soc_sink_all_reallocate
  } # scenario
  
} # caross valid

response_curve_scenario_raw[, 1, 4, ]
response_curve_scenario_raw[, 2, 4, ]
response_curve_scenario_raw[, 3, 4, ]


soil_sink_summary_response_curve = c()

iscenario = 1
for (iscenario in 1:3) {
  data_middle = 
    rbind(
      cbind(
        year_series,
        apply(response_curve_scenario_raw[ , iscenario, 3, ], 1, mean, na.rm = TRUE), # sink half mean
        apply(response_curve_scenario_raw[ , iscenario, 3, ], 1, sd, na.rm = TRUE), # sink half std
        apply(response_curve_scenario_raw[ , iscenario, 1, ], 1, mean, na.rm = TRUE), # rate half mean
        apply(response_curve_scenario_raw[ , iscenario, 1, ], 1, sd, na.rm = TRUE), # rate half std
        iscenario+1,
        0.5 # fraction of reallocation
      ),
      cbind(
        year_series,
        apply(response_curve_scenario_raw[ , iscenario, 4, ], 1, mean, na.rm = TRUE), # sink all mean
        apply(response_curve_scenario_raw[ , iscenario, 4, ], 1, sd, na.rm = TRUE), # sink all std
        apply(response_curve_scenario_raw[ , iscenario, 2, ], 1, mean, na.rm = TRUE), # rate all mean
        apply(response_curve_scenario_raw[ , iscenario, 2, ], 1, sd, na.rm = TRUE), # rate all std
        iscenario+1,
        1 # fraction of reallocation
      )
    )
  soil_sink_summary_response_curve = rbind(soil_sink_summary_response_curve, data_middle)
}

colnames(soil_sink_summary_response_curve) = c('year', 'sink', 'sink_std', 'rate', 'rate_std', 'scenario', 'reallocation')
soil_sink_summary_response_curve = data.frame(soil_sink_summary_response_curve)
soil_sink_summary_response_curve$scenario[soil_sink_summary_response_curve$scenario == 4] = 0

current_data = rbind(soil_sink_summary_base, soil_sink_summary_response_curve)

if (vert_component == 1) {
  category_labels = c(expression(paste('Natural \n condition', sep = '')), expression(paste('Half input \n reallocation', sep = '')), expression(paste('All input \n reallocation', sep = '')))
}

if (vert_component == 2) {
  category_labels = c(expression(paste('Natural \n condition', sep = '')), expression(paste('50% increased \n vertical transport', sep = '')), expression(paste('100% increased \n vertical transport', sep = '')))
}

p_sink =
  ggplot() +
  geom_bar(data = current_data[current_data$year == 2010, ], aes(x = as.factor(reallocation), y = sink, fill = as.factor(scenario)), position = position_dodge(width = 0.75), width = 0.7, stat = 'identity', size = 0, color = 'white') +
  geom_errorbar(data = current_data[current_data$year == 2010, ], aes(x = as.factor(reallocation), y = sink, ymin = sink - 2*sink_std, ymax = sink + 2*sink_std, group = as.factor(scenario)), position = position_dodge(width = 0.75), stat = 'identity', size = 1, color = 'black', width = 0.2) +
  geom_hline(yintercept = 0, color = 'black', size = 1) +
  scale_fill_manual(name = c(''), labels = c('Input + climate change', 'Input change only', 'Climate change only'), values = c('#0072B2', '#009E73','#D55E00')) +
  scale_x_discrete(labels = category_labels) +
  scale_y_continuous(limits = c(-20, 90), n.breaks = 5) +
  theme_classic() +
  # add title
  labs(title = '', x = '', y = expression(paste('Soil carbon sink (1900 - 2010, Pg C)', sep = ''))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 0), legend.direction = 'vertical')  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) 

p_rate = 
  ggplot() +
  geom_bar(data = current_data[current_data$year == 2010, ], aes(x = as.factor(reallocation), y = rate, fill = as.factor(scenario)), position = position_dodge(width = 0.75), width = 0.7, stat = 'identity', size = 0, color = 'white') +
  geom_errorbar(data = current_data[current_data$year == 2010, ], aes(x = as.factor(reallocation), y = rate, ymin = rate - 2*rate_std, ymax = rate + 2*rate_std, group = as.factor(scenario)), position = position_dodge(width = 0.75), stat = 'identity', size = 1, color = 'black', width = 0.2) +
  geom_hline(yintercept = 0, color = 'black', size = 1) +
  scale_fill_manual(name = c(''), labels = c('Input + climate change', 'Input change only', 'Climate change only'), values = c('#0072B2', '#009E73','#D55E00')) +
  scale_x_discrete(labels = category_labels) +
  scale_y_continuous(limits = c(-0.6, 1.5), n.breaks = 5) +
  theme_classic() +
  # add title
  labs(title = '', x = '', y = expression(paste('Sink rate (2000s, Pg C yr'^'1', ')', sep = ''))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 0), legend.direction = 'vertical')  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.02, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))

jpeg(paste('./Soil_Sink/proda_global_seq_vulnerability_method_', as.character(vert_component), '.jpeg', sep = ''), width = 20, height = 12, units = 'in', res = 300)
plot_grid(
  NULL, p_sink, p_rate, NULL,
  nrow = 1,
  labels = c(' ', 'a', 'b', ' '),
  rel_widths = c(0.02, 1, 1, 0.02),
  label_size = 50,
  label_x = 0., label_y = 1.01,
  label_fontfamily = 'Arial',
  label_fontface = 'bold'
)
dev.off()


#-------------------------rate time series
p_rate_base =
  ggplot() +
  geom_bar(data = current_data[current_data$reallocation == 0, ], aes(x = year, y = rate, color = as.factor(scenario), fill = as.factor(scenario)), alpha = 1, position = position_dodge(width = 8.5), width = 8, stat = 'identity', size = 0) +
  geom_errorbar(data = current_data[current_data$reallocation == 0, ], aes(x = year, y = rate, ymin = rate - 2*rate_std, ymax = rate + 2*rate_std, group = as.factor(scenario)), position = position_dodge(width = 8.5), color = 'black', width = 4, stat = 'identity', size = 1) +
  geom_hline(yintercept = 0, color = 'black', size = 1) +
  scale_fill_manual(name = c(''), labels = c('Input + climate change', 'Input change only', 'Climate change only'), values = c('#0072B2', '#009E73','#D55E00')) +
  scale_color_manual(name = c(''), labels = c('Input + climate change', 'Input change only', 'Climate change only'), values = c('#0072B2', '#009E73','#D55E00')) +
  scale_x_continuous(breaks = c(1910, 1930, 1950, 1970, 1990, 2010)) + 
  scale_y_continuous(limits = c(-0.6, 1.5), n.breaks = 5) +
  # scale_y_continuous(limits = c(-5, 25), n.breaks = 7) +
  theme_classic() +
  # add title
  labs(title = '', x = '', y = expression(paste('Stock change rate (Pg C yr'^'-1', ')', sep = ''))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 0), legend.direction = 'vertical')  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.02, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))


p_rate_half_reallocate =
  ggplot() +
  geom_bar(data = current_data[current_data$reallocation == 0.5, ], aes(x = year, y = rate, color = as.factor(scenario), fill = as.factor(scenario)), alpha = 1, position = position_dodge(width = 8.5), width = 8, stat = 'identity', size = 0) +
  geom_errorbar(data = current_data[current_data$reallocation == 0.5, ], aes(x = year, y = rate, ymin = rate - 2*rate_std, ymax = rate + 2*rate_std, group = as.factor(scenario)), position = position_dodge(width = 8.5), color = 'black', width = 4, stat = 'identity', size = 1) +
  geom_hline(yintercept = 0, color = 'black', size = 1) +
  scale_fill_manual(name = c(''), labels = c('Input + climate change', 'Input change only', 'Climate change only'), values = c('#0072B2', '#009E73','#D55E00')) +
  scale_color_manual(name = c(''), labels = c('Input + climate change', 'Input change only', 'Climate change only'), values = c('#0072B2', '#009E73','#D55E00')) +
  scale_x_continuous(breaks = c(1910, 1930, 1950, 1970, 1990, 2010)) + 
  scale_y_continuous(limits = c(-0.6, 1.5), n.breaks = 5) +
  # scale_y_continuous(limits = c(-5, 25), n.breaks = 7) +
  theme_classic() +
  # add title
  labs(title = '', x = '', y = expression(paste('Stock change rate (Pg C yr'^'-1', ')', sep = ''))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 0), legend.direction = 'vertical')  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.02, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))




p_rate_all_reallocate =
  ggplot() +
  geom_bar(data = current_data[current_data$reallocation == 1, ], aes(x = year, y = rate, color = as.factor(scenario), fill = as.factor(scenario)), alpha = 1, position = position_dodge(width = 8.5), width = 8, stat = 'identity', size = 0) +
  geom_errorbar(data = current_data[current_data$reallocation == 1, ], aes(x = year, y = rate, ymin = rate - 2*rate_std, ymax = rate + 2*rate_std, group = as.factor(scenario)), position = position_dodge(width = 8.5), color = 'black', width = 4, stat = 'identity', size = 1) +
  geom_hline(yintercept = 0, color = 'black', size = 1) +
  scale_fill_manual(name = c(''), labels = c('Input + climate change', 'Input change only', 'Climate change only'), values = c('#0072B2', '#009E73','#D55E00')) +
  scale_color_manual(name = c(''), labels = c('Input + climate change', 'Input change only', 'Climate change only'), values = c('#0072B2', '#009E73','#D55E00')) +
  scale_x_continuous(breaks = c(1910, 1930, 1950, 1970, 1990, 2010)) + 
  scale_y_continuous(limits = c(-0.6, 1.5), n.breaks = 5) +
  # scale_y_continuous(limits = c(-5, 25), n.breaks = 7) +
  theme_classic() +
  # add title
  labs(title = '', x = '', y = expression(paste('Stock change rate (Pg C yr'^'-1', ')', sep = ''))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 0), legend.direction = 'vertical')  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.02, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))


#-------------------------sink time series
p_sink_base =
  ggplot() +
  geom_bar(data = current_data[current_data$reallocation == 0, ], aes(x = year, y = sink, color = as.factor(scenario), fill = as.factor(scenario)), alpha = 1, position = position_dodge(width = 8.5), width = 8, stat = 'identity', size = 0) +
  geom_errorbar(data = current_data[current_data$reallocation == 0, ], aes(x = year, y = sink, ymin = sink - 2*sink_std, ymax = sink + 2*sink_std, group = as.factor(scenario)), position = position_dodge(width = 8.5), color = 'black', width = 4, stat = 'identity', size = 1) +
  geom_hline(yintercept = 0, color = 'black', size = 1) +
  scale_fill_manual(name = c(''), labels = c('Input + climate change', 'Input change only', 'Climate change only'), values = c('#0072B2', '#009E73','#D55E00')) +
  scale_color_manual(name = c(''), labels = c('Input + climate change', 'Input change only', 'Climate change only'), values = c('#0072B2', '#009E73','#D55E00')) +
  scale_x_continuous(breaks = c(1910, 1930, 1950, 1970, 1990, 2010)) + 
  scale_y_continuous(limits = c(-20, 90), n.breaks = 5) +
  # scale_y_continuous(limits = c(-5, 25), n.breaks = 7) +
  theme_classic() +
  # add title
  labs(title = 'Baseline simulation', x = '', y = expression(paste('Cumulative stock change (Pg C)', sep = ''))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 0), legend.direction = 'vertical')  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.02, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))



p_sink_half_reallocate =
  ggplot() +
  geom_bar(data = current_data[current_data$reallocation == 0.5, ], aes(x = year, y = sink, color = as.factor(scenario), fill = as.factor(scenario)), alpha = 1, position = position_dodge(width = 8.5), width = 8, stat = 'identity', size = 0) +
  geom_errorbar(data = current_data[current_data$reallocation == 0.5, ], aes(x = year, y = sink, ymin = sink - 2*sink_std, ymax = sink + 2*sink_std, group = as.factor(scenario)), position = position_dodge(width = 8.5), color = 'black', width = 4, stat = 'identity', size = 1) +
  geom_hline(yintercept = 0, color = 'black', size = 1) +
  scale_fill_manual(name = c(''), labels = c('Input + climate change', 'Input change only', 'Climate change only'), values = c('#0072B2', '#009E73','#D55E00')) +
  scale_color_manual(name = c(''), labels = c('Input + climate change', 'Input change only', 'Climate change only'), values = c('#0072B2', '#009E73','#D55E00')) +
  scale_x_continuous(breaks = c(1910, 1930, 1950, 1970, 1990, 2010)) + 
  scale_y_continuous(limits = c(-20, 90), n.breaks = 5) +
  # scale_y_continuous(limits = c(-5, 25), n.breaks = 7) +
  theme_classic() +
  # add title
  labs(title = 'Half input reallocation', x = '', y = expression(paste('Cumulative stock change (Pg C)', sep = ''))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 0), legend.direction = 'vertical')  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.02, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))

p_sink_all_reallocate =
  ggplot() +
  geom_bar(data = current_data[current_data$reallocation == 1, ], aes(x = year, y = sink, color = as.factor(scenario), fill = as.factor(scenario)), alpha = 1, position = position_dodge(width = 8.5), width = 8, stat = 'identity', size = 0) +
  geom_errorbar(data = current_data[current_data$reallocation == 1, ], aes(x = year, y = sink, ymin = sink - 2*sink_std, ymax = sink + 2*sink_std, group = as.factor(scenario)), position = position_dodge(width = 8.5), color = 'black', width = 4, stat = 'identity', size = 1) +
  geom_hline(yintercept = 0, color = 'black', size = 1) +
  scale_fill_manual(name = c(''), labels = c('Input + climate change', 'Input change only', 'Climate change only'), values = c('#0072B2', '#009E73','#D55E00')) +
  scale_color_manual(name = c(''), labels = c('Input + climate change', 'Input change only', 'Climate change only'), values = c('#0072B2', '#009E73','#D55E00')) +
  scale_x_continuous(breaks = c(1910, 1930, 1950, 1970, 1990, 2010)) + 
  scale_y_continuous(limits = c(-20, 90), n.breaks = 5) +
  # scale_y_continuous(limits = c(-5, 25), n.breaks = 7) +
  theme_classic() +
  # add title
  labs(title = 'All input reallocation', x = '', y = expression(paste('Cumulative stock change (Pg C)', sep = ''))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 0), legend.direction = 'vertical')  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.02, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))


jpeg(paste('./Soil_Sink/proda_global_seq_vulnerability_time_seris_method_', as.character(vert_component), '.jpeg', sep = ''), width = 27, height = 27, units = 'in', res = 300)
plot_grid(
  NULL, p_sink_base, p_rate_base, NULL,
  NULL, p_sink_half_reallocate, p_rate_half_reallocate, NULL,
  NULL, p_sink_all_reallocate, p_rate_all_reallocate, NULL,
  nrow = 3,
  labels = c(' ', 'a', 'b', ' ',
             ' ', 'c', 'd', ' ',
             ' ', 'e', 'f', ' '
  ),
  rel_widths = c(0.02, 1, 1, 0.02),
  label_size = 50,
  label_x = 0., label_y = 1.01,
  label_fontfamily = 'Arial',
  label_fontface = 'bold'
)
dev.off()
